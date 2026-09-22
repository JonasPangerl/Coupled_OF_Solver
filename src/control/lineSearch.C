/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
\*---------------------------------------------------------------------------*/

#include "lineSearch.H"
#include "coupledConstants.H"
#include "coupledDefaults.H"
#include "PstreamReduceOps.H"
#include "DynamicList.H"
#include "ITstream.H"
#include <cmath>

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::lineSearch::lineSearch(const dictionary& coupledDict)
:
    fU_(coupledDefaults::fU),
    fp_(coupledDefaults::fp),
    omegaMin_(coupledDefaults::omegaMin),
    kappa_(coupledDefaults::kappa),
    maxCflCuts_(coupledDefaults::maxCflCuts),
    beta_(coupledDefaults::lineSearchBeta),
    UrefMode_(coupledDefaults::UrefMode),
    UrefExplicit_(0),
    Uref_(0),
    pref_(0),
    refSet_(false),
    UrefSource_(),
    stepMode_(coupledDefaults::UrefStepMode),
    stepCap_(coupledDefaults::UrefStepCap),
    UstepExplicit_(0),
    Ufield0_(0),
    Ustep_(0),
    startupRef_(coupledDefaults::startupReference),
    UrefEff_(0),
    prefEff_(0),
    UstepEff_(0),
    pstepEff_(0),
    excludeDynamic_(false)
{
    // coupled.Uref: a mode word or an explicit velocity scale (D-050)
    if (coupledDict.found("Uref", keyType::LITERAL))
    {
        ITstream& is = coupledDict.lookup("Uref", keyType::LITERAL);
        const token tok(is);
        if (tok.isNumber())
        {
            UrefMode_ = "explicit";
            UrefExplicit_ = tok.number();
            // Negated comparison also rejects non-finite input
            if (!(UrefExplicit_ > 0) || !std::isfinite(UrefExplicit_))
            {
                FatalIOErrorInFunction(coupledDict)
                    << "coupled.Uref must be boundary, field or a finite"
                    << " value > 0, got " << UrefExplicit_
                    << exit(FatalIOError);
            }
        }
        else if
        (
            tok.isWord()
         && (tok.wordToken() == "boundary" || tok.wordToken() == "field")
        )
        {
            UrefMode_ = tok.wordToken();
        }
        else
        {
            FatalIOErrorInFunction(coupledDict)
                << "coupled.Uref must be boundary, field or a value > 0,"
                << " got " << tok << exit(FatalIOError);
        }
    }

    // coupled.UrefStep: a mode word or an explicit step scale (D-057)
    if (coupledDict.found("UrefStep", keyType::LITERAL))
    {
        ITstream& is = coupledDict.lookup("UrefStep", keyType::LITERAL);
        const token tok(is);
        if (tok.isNumber())
        {
            stepMode_ = "explicit";
            UstepExplicit_ = tok.number();
            if (!(UstepExplicit_ > 0) || !std::isfinite(UstepExplicit_))
            {
                FatalIOErrorInFunction(coupledDict)
                    << "coupled.UrefStep must be reference, fieldCapped,"
                    << " field or a finite value > 0, got " << UstepExplicit_
                    << exit(FatalIOError);
            }
        }
        else if
        (
            tok.isWord()
         && (
                tok.wordToken() == "reference"
             || tok.wordToken() == "fieldCapped"
             || tok.wordToken() == "field"
            )
        )
        {
            stepMode_ = tok.wordToken();
        }
        else
        {
            FatalIOErrorInFunction(coupledDict)
                << "coupled.UrefStep must be reference, fieldCapped, field"
                << " or a value > 0, got " << tok << exit(FatalIOError);
        }
    }
    stepCap_ = coupledDict.getOrDefault<scalar>
    (
        "UrefStepCap",
        coupledDefaults::UrefStepCap
    );
    if (!(stepCap_ >= 1) || !std::isfinite(stepCap_))
    {
        FatalIOErrorInFunction(coupledDict)
            << "coupled.UrefStepCap must be a finite value >= 1, got "
            << stepCap_ << exit(FatalIOError);
    }
    startupRef_ = coupledDict.getOrDefault<word>
    (
        "startupReference",
        coupledDefaults::startupReference
    );
    if
    (
        startupRef_ != "none" && startupRef_ != "ramp"
     && startupRef_ != "exclude"
    )
    {
        FatalIOErrorInFunction(coupledDict)
            << "coupled.startupReference must be none, ramp or exclude, got "
            << startupRef_ << exit(FatalIOError);
    }

    const dictionary& d = coupledDict.subOrEmptyDict("lineSearch");
    fU_ = d.getOrDefault<scalar>("fU", coupledDefaults::fU);
    fp_ = d.getOrDefault<scalar>("fp", coupledDefaults::fp);
    omegaMin_ = d.getOrDefault<doubleScalar>("omegaMin", coupledDefaults::omegaMin);
    kappa_ = d.getOrDefault<doubleScalar>("kappa", coupledDefaults::kappa);
    maxCflCuts_ =
        d.getOrDefault<label>("maxCflCuts", coupledDefaults::maxCflCuts);
    beta_ = d.getOrDefault<doubleScalar>("beta", coupledDefaults::lineSearchBeta);

    // Negated comparison also rejects non-finite input
    if (!(beta_ >= 1) || !std::isfinite(beta_))
    {
        FatalIOErrorInFunction(d)
            << "lineSearch.beta (CFL boost on a full step) must be a finite"
            << " value >= 1, got " << beta_ << exit(FatalIOError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::lineSearch::setReference(const volVectorField& U)
{
    const scalar Ufield = gMax(mag(U.primitiveField())());

    scalar Ubnd = 0;
    forAll(U.boundaryField(), patchi)
    {
        const fvPatchVectorField& Up = U.boundaryField()[patchi];
        if (!Up.coupled() && Up.size())
        {
            Ubnd = max(Ubnd, max(mag(Up)()));
        }
    }
    reduce(Ubnd, maxOp<scalar>());

    scalar Umax = 0;
    if (UrefMode_ == "explicit")
    {
        Umax = UrefExplicit_;
        UrefSource_ = "explicit";
    }
    else if (UrefMode_ == "field")
    {
        // The definition before D-050
        Umax = max(Ufield, Ubnd);
        UrefSource_ = "field";
    }
    else if (Ubnd > coupledDefaults::UrefFallbackFactor*Ufield)
    {
        Umax = Ubnd;
        UrefSource_ = "boundary";
    }
    else
    {
        // No driving boundary velocity: the field maximum
        Umax = Ufield;
        UrefSource_ = "fieldFallback";
    }

    // GUARD: a zero reference would make every omega zero
    Uref_ = max(Umax, cfVSmall<scalar>());
    pref_ = 0.5*sqr(Uref_);
    Ufield0_ = max(Ufield, Ubnd);
    refSet_ = true;
    updateStep();
}


void Foam::lineSearch::setFieldScale(const volVectorField& U)
{
    scalar Ubnd = 0;
    forAll(U.boundaryField(), patchi)
    {
        const fvPatchVectorField& Up = U.boundaryField()[patchi];
        if (!Up.coupled() && Up.size())
        {
            Ubnd = max(Ubnd, max(mag(Up)()));
        }
    }
    reduce(Ubnd, maxOp<scalar>());
    Ufield0_ = max(gMax(mag(U.primitiveField())()), Ubnd);
    updateStep();
}


void Foam::lineSearch::setReference
(
    const scalar Uref,
    const scalar pref,
    const scalar Ufield0
)
{
    Uref_ = Uref;
    pref_ = pref;
    // A state without Ufield0 (before D-057): the step scale falls back to
    // Uref for every mode except explicit
    Ufield0_ = (Ufield0 > 0 ? Ufield0 : Uref);
    refSet_ = true;
    UrefSource_ = "restart";
    updateStep();
}


void Foam::lineSearch::updateStep()
{
    if (stepMode_ == "explicit")
    {
        Ustep_ = UstepExplicit_;
    }
    else if (stepMode_ == "field")
    {
        Ustep_ = max(Uref_, Ufield0_);
    }
    else if (stepMode_ == "fieldCapped")
    {
        Ustep_ = max(Uref_, min(Ufield0_, stepCap_*Uref_));
    }
    else
    {
        Ustep_ = Uref_;
    }

    UrefEff_ = Uref_;
    prefEff_ = pref_;
    // reference mode: pstep is pref itself (bitwise D-050 behaviour, also
    // for a restored pref)
    UstepEff_ = Ustep_;
    pstepEff_ = (stepMode_ == "reference" ? pref_ : 0.5*sqr(Ustep_));
    excludeDynamic_ = false;
}


void Foam::lineSearch::setStartup
(
    const scalar betaStartup,
    const bool startupDone
)
{
    updateStep();

    if (startupDone || startupRef_ == "none")
    {
        return;
    }

    if (startupRef_ == "ramp")
    {
        // w = 1 at beta 0 (field scale), 0 at beta 1 (frozen values)
        const scalar w = min(max(1 - betaStartup, scalar(0)), scalar(1));
        UrefEff_ = Uref_ + w*max(Ufield0_ - Uref_, scalar(0));
        prefEff_ = 0.5*sqr(UrefEff_);
        UstepEff_ = max(Ustep_, UrefEff_);
        pstepEff_ = max(pstepEff_, 0.5*sqr(UstepEff_));
    }
    else if (startupRef_ == "exclude" && Ufield0_ > Ustep_)
    {
        // Only a start whose initial field exceeds the step scale (the
        // singular potential-flow peaks) is affected; otherwise a no-op
        UstepEff_ = Ufield0_;
        pstepEff_ = max(pstepEff_, 0.5*sqr(UstepEff_));
        excludeDynamic_ = true;
    }
}


Foam::doubleScalar Foam::lineSearch::omega(const blockScalarUList& dx) const
{
    const label nCells = dx.size()/blockDim;
    doubleScalar om = 1;

    const scalar limU = fU_*UstepEff_;
    const scalar limp = fp_*pstepEff_;

    for (label celli = 0; celli < nCells; ++celli)
    {
        const blockScalar* d = dx.cdata() + celli*blockDim;
        const scalar dU = std::sqrt
        (
            sqr(scalar(d[0])) + sqr(scalar(d[1])) + sqr(scalar(d[2]))
        );
        const scalar dp = std::abs(scalar(d[blockP]));

        // GUARD: line-search denominators (9.2)
        om = min(om, doubleScalar(limU/max(dU, cfVSmall<scalar>())));
        om = min(om, doubleScalar(limp/max(dp, cfVSmall<scalar>())));
    }

    reduce(om, minOp<doubleScalar>());
    return min(om, doubleScalar(1));
}


Foam::labelList Foam::lineSearch::offendingCells
(
    const blockScalarUList& dx
) const
{
    const label nCells = dx.size()/blockDim;
    // GUARD: omegaMin > 0 by construction of the dictionary checks
    const scalar limU = scalar(fU_*UstepEff_/max(omegaMin_, cfVSmall<doubleScalar>()));
    const scalar limp = scalar(fp_*pstepEff_/max(omegaMin_, cfVSmall<doubleScalar>()));

    DynamicList<label> cells;
    for (label celli = 0; celli < nCells; ++celli)
    {
        const blockScalar* d = dx.cdata() + celli*blockDim;
        const scalar dU = std::sqrt
        (
            sqr(scalar(d[0])) + sqr(scalar(d[1])) + sqr(scalar(d[2]))
        );
        const scalar dp = std::abs(scalar(d[blockP]));
        if (dU > limU || dp > limp)
        {
            cells.append(celli);
        }
    }
    return labelList(std::move(cells));
}


void Foam::lineSearch::countViolations
(
    const blockScalarUList& dx,
    const doubleScalar omega,
    label& nU,
    label& np
) const
{
    nU = 0;
    np = 0;
    const label nCells = dx.size()/blockDim;
    const scalar limU = fU_*UstepEff_;
    const scalar limp = fp_*pstepEff_;
    for (label celli = 0; celli < nCells; ++celli)
    {
        const blockScalar* d = dx.cdata() + celli*blockDim;
        const scalar dU = scalar
        (
            omega*std::sqrt
            (
                sqr(scalar(d[0])) + sqr(scalar(d[1])) + sqr(scalar(d[2]))
            )
        );
        const scalar dp = scalar(omega*std::abs(scalar(d[blockP])));
        if (dU > limU)
        {
            ++nU;
        }
        if (dp > limp)
        {
            ++np;
        }
    }
}


void Foam::lineSearch::writeSettings(dictionary& dict) const
{
    dictionary d;
    d.add("fU", fU_);
    d.add("fp", fp_);
    d.add("omegaMin", omegaMin_);
    d.add("kappa", kappa_);
    d.add("maxCflCuts", maxCflCuts_);
    d.add("beta", beta_);
    dict.add("lineSearch", d);
    if (UrefMode_ == "explicit")
    {
        dict.add("Uref", UrefExplicit_);
    }
    else
    {
        dict.add("Uref", UrefMode_);
    }
    if (stepMode_ == "explicit")
    {
        dict.add("UrefStep", UstepExplicit_);
    }
    else
    {
        dict.add("UrefStep", stepMode_);
    }
    dict.add("UrefStepCap", stepCap_);
    dict.add("startupReference", startupRef_);
}


// ************************************************************************* //

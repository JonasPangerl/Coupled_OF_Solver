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
#include "mixedFvPatchFields.H"
#include "directionMixedFvPatchFields.H"
#include <cmath>

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace
{

// Largest velocity magnitude a boundary condition PRESCRIBES on its patch,
// independent of the interior solution; 0 if it prescribes none.
//
// Only this may enter Uref = boundary (D-074). The previous rule took
// max|U_b| over every non-coupled patch, and on a patch that does not
// prescribe U the boundary value is a copy or projection of the interior:
// slip and symmetry take the tangential interior velocity, zeroGradient
// the interior velocity, inletOutlet the interior velocity on outflow. After
// a potential-flow start the interior carries the singular peaks at sharp
// edges, so the F1 half-car got Uref = 321 m/s against a free stream of 50,
// which scaled every safety limit up by 6 (and the pressure ones by 41).
Foam::scalar prescribedMagU(const Foam::fvPatchVectorField& Up)
{
    using namespace Foam;

    if (!Up.fixesValue())
    {
        // zeroGradient, slip, symmetry, calculated, ...
        return 0;
    }

    scalar m = 0;

    if (const auto* mp = isA<mixedFvPatchVectorField>(Up))
    {
        // inletOutlet & co: refValue with weight valueFraction (0 on
        // outflow faces, where the value is the interior's)
        const vectorField& rv = mp->refValue();
        const scalarField& vf = mp->valueFraction();
        forAll(rv, facei)
        {
            m = max(m, vf[facei]*mag(rv[facei]));
        }
    }
    else if (const auto* dp = isA<directionMixedFvPatchVectorField>(Up))
    {
        // pressureInletOutletVelocity & co: the prescribed part of refValue
        // is its projection with the valueFraction tensor
        const vectorField& rv = dp->refValue();
        const symmTensorField& vf = dp->valueFraction();
        forAll(rv, facei)
        {
            m = max(m, mag(vf[facei] & rv[facei]));
        }
    }
    else
    {
        // fixedValue family (fixedValue, noSlip, rotatingWallVelocity,
        // movingWallVelocity, flowRateInletVelocity, ...)
        for (const vector& u : Up)
        {
            m = max(m, mag(u));
        }
    }

    return m;
}

} // End anonymous namespace

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::lineSearch::lineSearch(const dictionary& coupledDict)
:
    fU_(scalar(coupledDefaults::fU)),
    fp_(scalar(coupledDefaults::fp)),
    omegaMin_(coupledDefaults::omegaMin),
    kappa_(coupledDefaults::kappa),
    maxCflCuts_(coupledDefaults::maxCflCuts),
    beta_(coupledDefaults::lineSearchBeta),
    mode_(coupledDefaults::lineSearchMode),
    localFraction_(coupledDefaults::lineSearchLocalFraction),
    UrefMode_(coupledDefaults::UrefMode),
    UrefExplicit_(0),
    Uref_(0),
    pref_(0),
    refSet_(false),
    UrefSource_(),
    stepMode_(coupledDefaults::UrefStepMode),
    stepCap_(coupledDefaults::UrefStepCap),
    UrefFallbackFactor_(scalar(coupledDefaults::UrefFallbackFactor)),
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
    UrefFallbackFactor_ = coupledDict.getOrDefault<scalar>
    (
        "UrefFallbackFactor",
        UrefFallbackFactor_
    );
    if (!(UrefFallbackFactor_ >= 0 && UrefFallbackFactor_ <= 1))
    {
        FatalIOErrorInFunction(coupledDict)
            << "coupled.UrefFallbackFactor must be in [0, 1], got "
            << UrefFallbackFactor_ << exit(FatalIOError);
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
    fU_ = d.getOrDefault<scalar>("fU", scalar(coupledDefaults::fU));
    fp_ = d.getOrDefault<scalar>("fp", scalar(coupledDefaults::fp));
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

    mode_ = d.getOrDefault<word>("mode", word(coupledDefaults::lineSearchMode));
    if (mode_ != "global" && mode_ != "local")
    {
        FatalIOErrorInFunction(d)
            << "lineSearch.mode must be global or local, got " << mode_
            << exit(FatalIOError);
    }
    localFraction_ = d.getOrDefault<doubleScalar>
    (
        "localFraction",
        coupledDefaults::lineSearchLocalFraction
    );
    if (!(localFraction_ >= 0 && localFraction_ <= 1))
    {
        FatalIOErrorInFunction(d)
            << "lineSearch.localFraction must be in [0, 1], got "
            << localFraction_ << exit(FatalIOError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::lineSearch::setReference(const volVectorField& U)
{
    const scalar Ufield = gMax(mag(U.primitiveField())());

    // Only what the boundary conditions prescribe (D-074). A boundary value
    // that merely copies the interior is already in Ufield.
    scalar Ubnd = 0;
    forAll(U.boundaryField(), patchi)
    {
        const fvPatchVectorField& Up = U.boundaryField()[patchi];
        if (!Up.coupled() && Up.size())
        {
            Ubnd = max(Ubnd, prescribedMagU(Up));
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
    else if (Ubnd > UrefFallbackFactor_*Ufield)
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
    pref_ = scalar(0.5*sqr(Uref_));
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
    pstepEff_ = (stepMode_ == "reference" ? pref_ : scalar(0.5*sqr(Ustep_)));
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
        prefEff_ = scalar(0.5*sqr(UrefEff_));
        UstepEff_ = max(Ustep_, UrefEff_);
        pstepEff_ = max(pstepEff_, scalar(0.5*sqr(UstepEff_)));
    }
    else if (startupRef_ == "exclude" && Ufield0_ > Ustep_)
    {
        // Only a start whose initial field exceeds the step scale (the
        // singular potential-flow peaks) is affected; otherwise a no-op
        UstepEff_ = Ufield0_;
        pstepEff_ = max(pstepEff_, scalar(0.5*sqr(UstepEff_)));
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


Foam::doubleScalar Foam::lineSearch::omegaQuantile
(
    const blockScalarUList& dx,
    labelList& violating
) const
{
    // Per-cell step fraction that brings the cell's increment back to the
    // bounds; histogram of log10 over [1e-8, 1), 20 bins per decade
    scalarField omCell;
    const labelList lim(limitedCells(dx, omCell));

    constexpr label nDec = 8, perDec = 20, nBins = nDec*perDec;
    labelList hist(nBins + 1, Zero);     // last bin: below 1e-8
    forAll(omCell, k)
    {
        const scalar o = max(omCell[k], scalar(1e-30));
        label b = label(std::floor(-std::log10(o)*perDec));
        b = min(max(b, label(0)), nBins);
        ++hist[b];
    }
    Pstream::listReduce(hist, sumOp<label>());

    const label nTot = returnReduce(dx.size()/blockDim, sumOp<label>());
    const scalar allowed = localFraction_*scalar(nTot);

    // Bin b holds omega_P in (10^-(b+1)/perDec, 10^-b/perDec]. Going from
    // the smallest omega_P upwards, accumulate the cells that would still
    // violate at omega = the upper edge of the bins above; take the largest
    // omega whose violator count stays within the allowance.
    doubleScalar omega = 1;
    label above = 0;                     // cells with omega_P below the edge
    for (label b = 0; b <= nBins; ++b)
    {
        above += hist[b];
    }
    // above = all limited cells: they violate at omega = 1
    for (label b = 0; b <= nBins && scalar(above) > allowed; ++b)
    {
        // Lower the step to the lower edge of bin b: its cells comply now
        above -= hist[b];
        omega = std::pow(10.0, -doubleScalar(b + 1)/perDec);
    }

    DynamicList<label> v;
    forAll(lim, k)
    {
        if (omCell[k] < omega)
        {
            v.append(lim[k]);
        }
    }
    violating.transfer(v);
    return omega;
}


Foam::labelList Foam::lineSearch::limitedCells
(
    const blockScalarUList& dx,
    scalarField& omegaCell
) const
{
    // The bounds of omega(), applied per cell: omega_P brings the cell's
    // own increment to fU Ustep and fp pstep
    const label nCells = dx.size()/blockDim;
    const scalar limU = fU_*UstepEff_;
    const scalar limp = fp_*pstepEff_;

    DynamicList<label> cells;
    DynamicList<scalar> om;
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
            // GUARD: dU > limU >= 0 or dp > limp >= 0, so the violated
            // denominator is > 0; the other one is guarded
            const scalar o = min
            (
                limU/max(dU, cfVSmall<scalar>()),
                limp/max(dp, cfVSmall<scalar>())
            );
            cells.append(celli);
            om.append(min(o, scalar(1)));
        }
    }
    omegaCell.transfer(om);
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
    d.add("mode", mode_);
    d.add("localFraction", localFraction_);
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
    dict.add("UrefFallbackFactor", UrefFallbackFactor_);
    dict.add("startupReference", startupRef_);
}


// ************************************************************************* //

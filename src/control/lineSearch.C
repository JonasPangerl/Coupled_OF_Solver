/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
\*---------------------------------------------------------------------------*/

#include "lineSearch.H"
#include "coupledDefaults.H"
#include "PstreamReduceOps.H"
#include "DynamicList.H"
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
    Uref_(0),
    pref_(0),
    refSet_(false)
{
    const dictionary& d = coupledDict.subOrEmptyDict("lineSearch");
    fU_ = d.getOrDefault<scalar>("fU", coupledDefaults::fU);
    fp_ = d.getOrDefault<scalar>("fp", coupledDefaults::fp);
    omegaMin_ = d.getOrDefault<scalar>("omegaMin", coupledDefaults::omegaMin);
    kappa_ = d.getOrDefault<scalar>("kappa", coupledDefaults::kappa);
    maxCflCuts_ =
        d.getOrDefault<label>("maxCflCuts", coupledDefaults::maxCflCuts);
    beta_ = d.getOrDefault<scalar>("beta", coupledDefaults::lineSearchBeta);

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
    scalar Umax = gMax(mag(U.primitiveField())());

    forAll(U.boundaryField(), patchi)
    {
        const fvPatchVectorField& Up = U.boundaryField()[patchi];
        if (!Up.coupled() && Up.size())
        {
            Umax = max(Umax, max(mag(Up)()));
        }
    }
    reduce(Umax, maxOp<scalar>());

    // GUARD: a zero reference would make every omega zero
    Uref_ = max(Umax, VSMALL);
    pref_ = 0.5*sqr(Uref_);
    refSet_ = true;
}


void Foam::lineSearch::setReference(const scalar Uref, const scalar pref)
{
    Uref_ = Uref;
    pref_ = pref;
    refSet_ = true;
}


Foam::scalar Foam::lineSearch::omega(const blockScalarUList& dx) const
{
    const label nCells = dx.size()/blockDim;
    scalar om = 1;

    const scalar limU = fU_*Uref_;
    const scalar limp = fp_*pref_;

    for (label celli = 0; celli < nCells; ++celli)
    {
        const blockScalar* d = dx.cdata() + celli*blockDim;
        const scalar dU = std::sqrt
        (
            sqr(scalar(d[0])) + sqr(scalar(d[1])) + sqr(scalar(d[2]))
        );
        const scalar dp = std::abs(scalar(d[blockP]));

        // GUARD: line-search denominators (9.2)
        om = min(om, limU/max(dU, VSMALL));
        om = min(om, limp/max(dp, VSMALL));
    }

    reduce(om, minOp<scalar>());
    return min(om, scalar(1));
}


Foam::labelList Foam::lineSearch::offendingCells
(
    const blockScalarUList& dx
) const
{
    const label nCells = dx.size()/blockDim;
    // GUARD: omegaMin > 0 by construction of the dictionary checks
    const scalar limU = fU_*Uref_/max(omegaMin_, VSMALL);
    const scalar limp = fp_*pref_/max(omegaMin_, VSMALL);

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
}


// ************************************************************************* //

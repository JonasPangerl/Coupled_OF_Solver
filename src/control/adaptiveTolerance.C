/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
\*---------------------------------------------------------------------------*/

#include "adaptiveTolerance.H"
#include "coupledDefaults.H"
#include <cmath>

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::adaptiveTolerance::adaptiveTolerance(const dictionary& d)
:
    enabled_(d.getOrDefault<bool>("adaptiveRelTol", coupledDefaults::adaptiveRelTol)),
    etaMin_(d.getOrDefault<scalar>("etaMin", coupledDefaults::etaMin)),
    etaMax_(d.getOrDefault<scalar>("etaMax", coupledDefaults::etaMax)),
    gamma_(d.getOrDefault<scalar>("gammaEW", coupledDefaults::gammaEW)),
    alpha_(d.getOrDefault<scalar>("alphaEW", coupledDefaults::alphaEW)),
    safeguard_
    (
        d.getOrDefault<scalar>("etaSafeguard", coupledDefaults::etaSafeguard)
    ),
    relTolFixed_(d.getOrDefault<scalar>("relTol", coupledDefaults::relTol)),
    etaPrev_(-1),
    Rprev_(-1)
{
    if (etaMin_ <= 0 || etaMax_ < etaMin_ || etaMax_ >= 1)
    {
        FatalIOErrorInFunction(d)
            << "Require 0 < etaMin <= etaMax < 1 (amendment B2)"
            << exit(FatalIOError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::adaptiveTolerance::eta
(
    const scalar R,
    const bool startupDone
)
{
    if (!enabled_)
    {
        return relTolFixed_;
    }

    // Safeguards 3 and 4
    if (!startupDone || Rprev_ <= 0 || etaPrev_ < 0)
    {
        return etaMax_;
    }

    // GUARD: Rprev > 0 checked above; pow base >= 0
    scalar e = gamma_*std::pow(max(R, scalar(0))/Rprev_, alpha_);

    // Safeguard 1
    const scalar sg = gamma_*std::pow(etaPrev_, alpha_);
    if (sg > safeguard_)
    {
        e = max(e, sg);
    }

    // Safeguard 2
    return min(max(e, etaMin_), etaMax_);
}


void Foam::adaptiveTolerance::accept(const scalar R, const scalar eta)
{
    Rprev_ = R;
    etaPrev_ = eta;
}


void Foam::adaptiveTolerance::writeState(dictionary& dict) const
{
    dict.set("etaPrev", etaPrev_);
    dict.set("RprevEW", Rprev_);
}


void Foam::adaptiveTolerance::readState(const dictionary& dict)
{
    etaPrev_ = dict.getOrDefault<scalar>("etaPrev", -1);
    Rprev_ = dict.getOrDefault<scalar>("RprevEW", -1);
}


void Foam::adaptiveTolerance::writeSettings(dictionary& dict) const
{
    dict.add("adaptiveRelTol", enabled_);
    dict.add("etaMin", etaMin_);
    dict.add("etaMax", etaMax_);
    dict.add("gammaEW", gamma_);
    dict.add("alphaEW", alpha_);
    dict.add("etaSafeguard", safeguard_);
    dict.add("relTol", relTolFixed_);
}


// ************************************************************************* //

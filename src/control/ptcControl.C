/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
\*---------------------------------------------------------------------------*/

#include "ptcControl.H"
#include "coupledDefaults.H"
#include "PstreamReduceOps.H"
#include <cmath>

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ptcControl::ptcControl
(
    const fvMesh& mesh,
    const dictionary& coupledDict
)
:
    mesh_(mesh),
    strategy_(strategy::mRDM),
    strategyName_("mRDM"),
    CFL0_(coupledDefaults::CFL0),
    CFLmin_(coupledDefaults::CFLmin),
    CFLmax_(coupledDefaults::CFLmax),
    gamma_(coupledDefaults::ptcGamma),
    betaMax_(coupledDefaults::betaMax),
    betaExp_(coupledDefaults::betaExp),
    nHold_(coupledDefaults::nHold),
    localLimit_(coupledDefaults::localLimitEnabled),
    fLoc_(coupledDefaults::fLoc),
    CFL_(coupledDefaults::CFL0),
    Rprev_(-1),
    holdRemaining_(0),
    nLocalLimited_(0)
{
    const dictionary& d = coupledDict.subOrEmptyDict("ptc");

    strategyName_ = d.getOrDefault<word>("cflStrategy", "mRDM");
    if (strategyName_ == "mRDM")
    {
        strategy_ = strategy::mRDM;
    }
    else if (strategyName_ == "EXP")
    {
        strategy_ = strategy::EXP;
    }
    else if (strategyName_ == "SER")
    {
        strategy_ = strategy::SER;
    }
    else
    {
        FatalIOErrorInFunction(d)
            << "Unknown cflStrategy " << strategyName_
            << ", valid: mRDM EXP SER" << exit(FatalIOError);
    }

    CFL0_ = d.getOrDefault<scalar>("CFL0", coupledDefaults::CFL0);
    CFLmin_ = d.getOrDefault<scalar>("CFLmin", coupledDefaults::CFLmin);
    CFLmax_ = d.getOrDefault<scalar>("CFLmax", coupledDefaults::CFLmax);
    gamma_ = d.getOrDefault<scalar>("gamma", coupledDefaults::ptcGamma);
    betaMax_ = d.getOrDefault<scalar>("betaMax", coupledDefaults::betaMax);
    betaExp_ = d.getOrDefault<scalar>("betaExp", coupledDefaults::betaExp);
    nHold_ = d.getOrDefault<label>("nHold", coupledDefaults::nHold);

    const dictionary& ll = coupledDict.subOrEmptyDict("localLimit");
    localLimit_ =
        ll.getOrDefault<bool>("enabled", coupledDefaults::localLimitEnabled);
    fLoc_ = ll.getOrDefault<scalar>("fLoc", coupledDefaults::fLoc);

    CFL_ = CFL0_;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField> Foam::ptcControl::rDeltaTV
(
    const surfaceScalarField& phi,
    const volScalarField& nuEff,
    const scalarField& cflFactor
) const
{
    const scalarField& V = mesh_.V();
    auto tr = tmp<scalarField>::New(mesh_.nCells(), Zero);
    scalarField& r = tr.ref();

    // 1/2 sum_f |phi_f|
    scalarField sumPhi(mesh_.nCells(), Zero);
    const labelUList& own = mesh_.owner();
    const labelUList& nei = mesh_.neighbour();
    forAll(own, facei)
    {
        const scalar a = mag(phi[facei]);
        sumPhi[own[facei]] += a;
        sumPhi[nei[facei]] += a;
    }
    forAll(phi.boundaryField(), patchi)
    {
        const labelUList& fc = mesh_.boundary()[patchi].faceCells();
        const scalarField& pp = phi.boundaryField()[patchi];
        forAll(fc, pf)
        {
            sumPhi[fc[pf]] += mag(pp[pf]);
        }
    }

    const scalarField& nu = nuEff.primitiveField();

    forAll(r, celli)
    {
        // GUARD: pow argument V > 0; denominator >= VSMALL (5.4)
        const scalar lambda =
            0.5*sumPhi[celli]
          + nu[celli]*std::cbrt(max(V[celli], VSMALL));
        const scalar cfl = max(CFL_*cflFactor[celli], VSMALL);
        r[celli] = max(lambda, VSMALL)/cfl;
    }

    return tr;
}


Foam::label Foam::ptcControl::applyLocalLimit
(
    scalarField& rDeltaTV,
    const vectorField& rMom,
    const scalar Uref
)
{
    nLocalLimited_ = 0;

    if (!localLimit_)
    {
        return 0;
    }

    const scalar dUmax = fLoc_*Uref;

    forAll(rDeltaTV, celli)
    {
        // dU = |r| dt/V = |r| / (V/dt); GUARD: rDeltaTV >= VSMALL/CFL
        const scalar dU = mag(rMom[celli])/max(rDeltaTV[celli], VSMALL);
        if (dU > dUmax)
        {
            // dt <- dt*dUmax/dU  <=>  V/dt <- V/dt * dU/dUmax
            rDeltaTV[celli] *= dU/max(dUmax, VSMALL);
            ++nLocalLimited_;
        }
    }

    reduce(nLocalLimited_, sumOp<label>());
    return nLocalLimited_;
}


void Foam::ptcControl::update(const scalar R)
{
    const scalar Rold = Rprev_;
    Rprev_ = R;

    if (Rold <= 0)
    {
        return;
    }

    // GUARD: R > 0 before division
    const scalar ratio = Rold/max(R, VSMALL);

    scalar CFLnew = CFL_;

    switch (strategy_)
    {
        case strategy::mRDM:
        {
            if (R <= Rold)
            {
                // GUARD: pow base > 0 (ratio >= 1 here)
                const scalar f = min(betaMax_, std::pow(ratio, gamma_));
                CFLnew = min(CFLmax_, CFL_*max(scalar(1), f));
            }
            break;
        }
        case strategy::EXP:
        {
            CFLnew = min(CFLmax_, CFL_*betaExp_);
            break;
        }
        case strategy::SER:
        {
            CFLnew = min(CFLmax_, max(CFLmin_, CFL_*ratio));
            break;
        }
    }

    if (holdRemaining_ > 0)
    {
        --holdRemaining_;
        // No increase during the hold (SER may still decrease)
        CFLnew = min(CFLnew, CFL_);
    }

    CFL_ = CFLnew;
}


void Foam::ptcControl::decrease(const scalar factor)
{
    CFL_ = max(CFLmin_, factor*CFL_);
    holdRemaining_ = nHold_;
}


void Foam::ptcControl::boost(const scalar factor)
{
    if (holdRemaining_ > 0 || !(factor > 1))
    {
        return;
    }
    CFL_ = min(CFLmax_, factor*CFL_);
}


void Foam::ptcControl::writeState(dictionary& dict) const
{
    dict.set("CFL", CFL_);
    dict.set("Rprev", Rprev_);
    dict.set("nHoldRemaining", holdRemaining_);
}


void Foam::ptcControl::readState(const dictionary& dict)
{
    CFL_ = dict.get<scalar>("CFL");
    Rprev_ = dict.getOrDefault<scalar>("Rprev", -1);
    holdRemaining_ = dict.getOrDefault<label>("nHoldRemaining", 0);
}


void Foam::ptcControl::writeSettings(dictionary& dict) const
{
    dictionary d;
    d.add("cflStrategy", strategyName_);
    d.add("CFL0", CFL0_);
    d.add("CFLmin", CFLmin_);
    d.add("CFLmax", CFLmax_);
    d.add("gamma", gamma_);
    d.add("betaMax", betaMax_);
    d.add("betaExp", betaExp_);
    d.add("nHold", nHold_);
    dict.add("ptc", d);

    dictionary l;
    l.add("enabled", localLimit_);
    l.add("fLoc", fLoc_);
    dict.add("localLimit", l);
}


// ************************************************************************* //

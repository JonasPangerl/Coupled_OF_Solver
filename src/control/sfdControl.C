/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
\*---------------------------------------------------------------------------*/

#include "sfdControl.H"
#include "coupledDefaults.H"
#include "PstreamReduceOps.H"
#include <cmath>

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sfdControl::sfdControl
(
    const fvMesh& mesh,
    const dictionary& coupledDict,
    const label nHold
)
:
    mesh_(mesh),
    enabled_(coupledDefaults::sfdEnabled),
    chi_(coupledDefaults::sfdChi),
    Delta_(coupledDefaults::sfdDelta),
    Lref_(coupledDefaults::sfdLref),
    deactivateBelowR_(coupledDefaults::sfdDeactivateBelowR),
    resetOnFlush_(coupledDefaults::sfdResetOnFlush),
    afterStartup_(coupledDefaults::sfdAfterStartup),
    nHold_(nHold),
    UbarPtr_(nullptr),
    chiStar_(0),
    DeltaStar_(0),
    initialised_(false),
    on_(false),
    quiet_(0),
    nResets_(0)
{
    const dictionary& d = coupledDict.subOrEmptyDict("sfd");
    enabled_ = d.getOrDefault<bool>("enabled", enabled_);
    chi_ = d.getOrDefault<scalar>("chi", chi_);
    Delta_ = d.getOrDefault<scalar>("Delta", Delta_);
    Lref_ = d.getOrDefault<scalar>("Lref", Lref_);
    deactivateBelowR_ =
        d.getOrDefault<scalar>("deactivateBelowR", deactivateBelowR_);
    resetOnFlush_ = d.getOrDefault<bool>("resetOnFlush", resetOnFlush_);
    afterStartup_ = d.getOrDefault<bool>("afterStartup", afterStartup_);

    // Negated comparisons also reject non-finite input
    if
    (
        !(chi_ >= 0) || !std::isfinite(chi_)
     || !(Delta_ > 0) || !std::isfinite(Delta_)
     || !(Lref_ > 0) || !std::isfinite(Lref_)
     || !(deactivateBelowR_ >= 0)
    )
    {
        FatalIOErrorInFunction(d)
            << "sfd: require chi >= 0, Delta > 0, Lref > 0 (finite) and"
            << " deactivateBelowR >= 0; got chi " << chi_ << ", Delta "
            << Delta_ << ", Lref " << Lref_ << ", deactivateBelowR "
            << deactivateBelowR_ << exit(FatalIOError);
    }

    if (enabled_)
    {
        UbarPtr_.reset
        (
            new volVectorField
            (
                IOobject
                (
                    "USFD",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    IOobject::NO_REGISTER
                ),
                mesh_,
                dimensionedVector(dimVelocity, Zero)
            )
        );
        on_ = true;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::volVectorField& Foam::sfdControl::Ubar() const
{
    if (!UbarPtr_)
    {
        FatalErrorInFunction
            << "SFD not enabled" << abort(FatalError);
    }
    return *UbarPtr_;
}


void Foam::sfdControl::setReference(const scalar Uref)
{
    if (!enabled_)
    {
        return;
    }

    // GUARD: Uref, Lref > 0 (Uref >= VSMALL from lineSearch)
    chiStar_ = chi_*Uref/Lref_;
    DeltaStar_ = Delta_*Lref_/max(Uref, VSMALL);

    Info<< "coupledFoam: SFD (7.6) chi* = chi Uref/Lref = " << chiStar_
        << " 1/s, Delta* = Delta Lref/Uref = " << DeltaStar_ << " s"
        << (on_ ? "" : " (switched off, restart)") << nl
        << "coupledFoam: SFD guidance: chi ~ 2x the non-dimensional growth"
        << " rate of the oscillation (per Lref/Uref), Delta ~ 1/(2 f_osc"
        << " Lref/Uref) with f_osc its frequency; if unknown, start with the"
        << " defaults chi " << coupledDefaults::sfdChi << ", Delta "
        << coupledDefaults::sfdDelta << endl;
}


void Foam::sfdControl::begin
(
    const volVectorField& U,
    const bool startupDone
)
{
    if (!enabled_ || initialised_ || !on_ || (afterStartup_ && !startupDone))
    {
        return;
    }
    volVectorField& Ub = *UbarPtr_;
    Ub.primitiveFieldRef() = U.primitiveField();
    Ub.boundaryFieldRef() == U.boundaryField();
    initialised_ = true;
    Info<< "coupledFoam: SFD active from iteration "
        << mesh_.time().timeIndex() << endl;
}


void Foam::sfdControl::reset(const volVectorField& U)
{
    if (!active() || !resetOnFlush_)
    {
        return;
    }
    volVectorField& Ub = *UbarPtr_;
    Ub.primitiveFieldRef() = U.primitiveField();
    Ub.boundaryFieldRef() == U.boundaryField();
    ++nResets_;
}


void Foam::sfdControl::update
(
    const volVectorField& U,
    const scalarField& rDeltaTV
)
{
    if (!active())
    {
        return;
    }

    volVectorField& Ub = *UbarPtr_;
    vectorField& ub = Ub.primitiveFieldRef();
    const vectorField& u = U.primitiveField();
    const scalarField& V = mesh_.V();

    forAll(ub, celli)
    {
        // GUARD: V/dt >= VSMALL (5.4); dt/(Delta* + dt) in [0, 1]
        const scalar dt = V[celli]/max(rDeltaTV[celli], VSMALL);
        const scalar a = dt/(DeltaStar_ + dt);
        ub[celli] += a*(u[celli] - ub[celli]);
    }
    Ub.boundaryFieldRef() == U.boundaryField();
}


bool Foam::sfdControl::checkOff(const scalar R, const label iter)
{
    if (!active())
    {
        return false;
    }

    quiet_ = (R < deactivateBelowR_ ? quiet_ + 1 : 0);

    if (quiet_ >= max(nHold_, label(1)))
    {
        on_ = false;
        Info<< "coupledFoam: SFD-off at iteration " << iter << " (R " << R
            << " < " << deactivateBelowR_ << " for " << quiet_
            << " iterations)" << endl;
        return true;
    }
    return false;
}


Foam::scalar Foam::sfdControl::maxDeviation
(
    const volVectorField& U,
    const scalar Uref
) const
{
    if (!enabled_ || !initialised_ || U.primitiveField().empty())
    {
        return 0;
    }
    const scalar m =
        max(mag(U.primitiveField() - UbarPtr_->primitiveField())());
    // GUARD: Uref >= VSMALL
    return m/max(Uref, VSMALL);
}


void Foam::sfdControl::write() const
{
    if (!enabled_ || !initialised_)
    {
        return;
    }
    volVectorField Uw
    (
        IOobject
        (
            "USFD",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        *UbarPtr_
    );
    Uw.write();
}


void Foam::sfdControl::writeState(dictionary& dict) const
{
    if (!enabled_)
    {
        return;
    }
    dict.set("sfdOn", on_);
    dict.set("sfdQuiet", quiet_);
    dict.set("sfdInitialised", initialised_);
}


void Foam::sfdControl::readState(const dictionary& dict)
{
    if (!enabled_)
    {
        return;
    }
    on_ = dict.getOrDefault<bool>("sfdOn", true);
    quiet_ = dict.getOrDefault<label>("sfdQuiet", 0);

    // Ubar from USFD of the restart time; without it (restart of a run
    // without SFD) Ubar is initialised to U in the first iteration
    IOobject io
    (
        "USFD",
        mesh_.time().timeName(),
        mesh_,
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        IOobject::NO_REGISTER
    );
    if
    (
        dict.getOrDefault<bool>("sfdInitialised", false)
     && io.typeHeaderOk<volVectorField>(true)
    )
    {
        const volVectorField Ur(io, mesh_);
        UbarPtr_->primitiveFieldRef() = Ur.primitiveField();
        UbarPtr_->boundaryFieldRef() == Ur.boundaryField();
        initialised_ = true;
        Info<< "coupledFoam: SFD filter USFD restored from "
            << mesh_.time().timeName() << endl;
    }
}


void Foam::sfdControl::writeSettings(dictionary& dict) const
{
    dictionary s;
    s.add("enabled", enabled_);
    s.add("chi", chi_);
    s.add("Delta", Delta_);
    s.add("Lref", Lref_);
    s.add("deactivateBelowR", deactivateBelowR_);
    s.add("resetOnFlush", resetOnFlush_);
    s.add("afterStartup", afterStartup_);
    s.add("nHold", nHold_);
    dict.add("sfd", s);
}


// ************************************************************************* //

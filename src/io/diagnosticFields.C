/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
\*---------------------------------------------------------------------------*/

#include "diagnosticFields.H"
#include "coupledConstants.H"
#include "coupledDefaults.H"
#include "extrapolatedCalculatedFvPatchFields.H"
#include "fileOperation.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diagnosticFields::diagnosticFields
(
    const fvMesh& mesh,
    const dictionary& coupledDict
)
:
    mesh_(mesh),
    enabled_
    (
        coupledDict.subOrEmptyDict("diagnosticFields").getOrDefault<bool>
        (
            "enabled",
            coupledDefaults::diagnosticFieldsEnabled
        )
    ),
    rDeltaTV_(),
    cflFactor_(),
    beta_(),
    USFD_(nullptr)
{}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::diagnosticFields::makeField
(
    const word& name,
    const word& instance,
    const dimensionSet& dims,
    const scalarField& values
) const
{
    auto tfld = tmp<volScalarField>::New
    (
        IOobject
        (
            name,
            instance,
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        mesh_,
        dimensionedScalar(dims, Zero),
        extrapolatedCalculatedFvPatchScalarField::typeName
    );
    tfld.ref().primitiveFieldRef() = values;
    tfld.ref().correctBoundaryConditions();
    return tfld;
}


void Foam::diagnosticFields::writeInstance(const regIOobject& io)
{
    fileHandler().writeObject
    (
        io,
        IOstreamOption(io.time().writeFormat(), io.time().writeCompression()),
        true
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::diagnosticFields::record
(
    const scalarField& rDeltaTV,
    const scalarField& cflFactor,
    const scalarField& beta
)
{
    if (!enabled_)
    {
        return;
    }
    rDeltaTV_ = rDeltaTV;
    cflFactor_ = cflFactor;
    beta_ = beta;
}


void Foam::diagnosticFields::write
(
    const surfaceScalarField& phi,
    const word& instance
) const
{
    // enabled_ and the recorded state are identical on all ranks (record()
    // is called in the same iteration everywhere): no rank skips the
    // collective boundary evaluation
    if (!enabled_ || rDeltaTV_.size() != mesh_.nCells())
    {
        return;
    }

    const scalarField& V = mesh_.V();

    // dt_P = V_P/(V_P/dt_P); GUARD: rDeltaTV > 0 by construction (5.4)
    scalarField dt(mesh_.nCells());
    forAll(dt, celli)
    {
        dt[celli] = V[celli]/max(rDeltaTV_[celli], cfVSmall<scalar>());
    }

    // 1/2 sum_f |phi_f| (as ptcControl::rDeltaTV)
    scalarField halfSumPhi(mesh_.nCells(), Zero);
    const labelUList& own = mesh_.owner();
    const labelUList& nei = mesh_.neighbour();
    forAll(own, facei)
    {
        const scalar a = mag(phi[facei]);
        halfSumPhi[own[facei]] += a;
        halfSumPhi[nei[facei]] += a;
    }
    forAll(phi.boundaryField(), patchi)
    {
        const labelUList& fc = mesh_.boundary()[patchi].faceCells();
        const scalarField& pp = phi.boundaryField()[patchi];
        forAll(fc, pf)
        {
            halfSumPhi[fc[pf]] += mag(pp[pf]);
        }
    }
    halfSumPhi *= 0.5;

    scalarField cfl(mesh_.nCells());
    forAll(cfl, celli)
    {
        // GUARD: V > 0 for a valid mesh
        cfl[celli] = dt[celli]*halfSumPhi[celli]/max(V[celli], cfVSmall<scalar>());
    }

    writeInstance(makeField("localDt", instance, dimTime, dt)());
    writeInstance(makeField("localCFL", instance, dimless, cfl)());
    writeInstance
    (
        makeField("cflFactorEff", instance, dimless, cflFactor_)()
    );
    writeInstance(makeField("betaEff", instance, dimless, beta_)());

    // C3 hook: SFD filtered velocity while SFD is active
    if (USFD_)
    {
        const volVectorField Ubar
        (
            IOobject
            (
                "USFD",
                instance,
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            *USFD_
        );
        writeInstance(Ubar);
    }
}


void Foam::diagnosticFields::write(const surfaceScalarField& phi) const
{
    write(phi, mesh_.time().timeName());
}


void Foam::diagnosticFields::writeSettings(dictionary& dict) const
{
    dictionary d;
    d.add("enabled", enabled_);
    dict.add("diagnosticFields", d);
}


// ************************************************************************* //

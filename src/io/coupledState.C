/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
\*---------------------------------------------------------------------------*/

#include "coupledState.H"
#include "coupledDefaults.H"
#include "localIOdictionary.H"
#include "fileOperation.H"
#include "OSspecific.H"
#include "surfaceFields.H"
#include <limits>

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coupledState::coupledState
(
    const fvMesh& mesh,
    const dictionary& coupledDict
)
:
    mesh_(mesh),
    writeState_
    (
        coupledDict.getOrDefault<bool>("writeState", coupledDefaults::writeState)
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::coupledState::read(dictionary& state) const
{
    IOobject io
    (
        "coupledState",
        mesh_.time().timeName(),
        mesh_,
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        IOobject::NO_REGISTER
    );

    // localIOdictionary, not IOdictionary: an IOdictionary is a global
    // object, read by the master (from processor0) and broadcast, so in a
    // parallel restart every rank got the state of rank 0 - its dynamic
    // set labels in particular (T0 np4: the lid-corner cells of ranks 2
    // and 3 lost their remediation, R jumped from 5.6e-5 to 1.3e-2)
    if (!io.typeHeaderOk<localIOdictionary>(true))
    {
        return false;
    }

    localIOdictionary d(io);
    state = d;
    return true;
}


void Foam::coupledState::write(const dictionary& state) const
{
    if (!writeState_)
    {
        return;
    }

    const bool uncollated = (fileHandler().type() == "uncollated");

    // Per-processor state: localIOdictionary (see read())
    localIOdictionary d
    (
        IOobject
        (
            uncollated ? "coupledState.tmp" : "coupledState",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        state
    );

    // Always ASCII, with round-trip precision: in a binary-format
    // dictionary an empty labelList "0 ( )" is not read back (the native
    // binary List reader consumes no brackets for zero size: "Entry
    // 'dynamicSet' has 2 excess tokens" on restart), and the scalars of
    // the state (CFL, R1, Rprev, ...) must restart bit-exactly.
    {
        const unsigned oldPrecision = IOstream::defaultPrecision
        (
            std::numeric_limits<doubleScalar>::max_digits10
        );
        d.regIOobject::writeObject
        (
            IOstreamOption(IOstreamOption::ASCII),
            true
        );
        IOstream::defaultPrecision(oldPrecision);
    }

    if (uncollated)
    {
        const fileName src(d.objectPath());
        const fileName dst(src.path()/"coupledState");
        if (!Foam::mv(src, dst))
        {
            FatalErrorInFunction
                << "Could not rename " << src << " to " << dst
                << exit(FatalError);
        }
    }
    else
    {
        Info<< "coupledState: file handler " << fileHandler().type()
            << ", written without tmp+rename (D-017)" << endl;
    }
}


void Foam::coupledState::writeD(const volScalarField& D) const
{
    if (!writeState_)
    {
        return;
    }

    volScalarField Dw
    (
        IOobject
        (
            "coupledD",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        D
    );
    Dw.write();
}


void Foam::coupledState::writeQ(const surfaceScalarField& q) const
{
    if (!writeState_)
    {
        return;
    }

    surfaceScalarField qw
    (
        IOobject
        (
            "coupledQ",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        q
    );
    qw.write();
}


Foam::tmp<Foam::surfaceScalarField> Foam::coupledState::readQ() const
{
    IOobject io
    (
        "coupledQ",
        mesh_.time().timeName(),
        mesh_,
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        IOobject::NO_REGISTER
    );

    if (!io.typeHeaderOk<surfaceScalarField>(true))
    {
        return nullptr;
    }

    return tmp<surfaceScalarField>::New(io, mesh_);
}


Foam::tmp<Foam::volScalarField> Foam::coupledState::readD() const
{
    IOobject io
    (
        "coupledD",
        mesh_.time().timeName(),
        mesh_,
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        IOobject::NO_REGISTER
    );

    if (!io.typeHeaderOk<volScalarField>(true))
    {
        return nullptr;
    }

    return tmp<volScalarField>::New(io, mesh_);
}


// ************************************************************************* //

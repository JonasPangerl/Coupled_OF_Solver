/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
\*---------------------------------------------------------------------------*/

#include "coupledState.H"
#include "coupledDefaults.H"
#include "IOdictionary.H"
#include "fileOperation.H"
#include "OSspecific.H"

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

    if (!io.typeHeaderOk<IOdictionary>(true))
    {
        return false;
    }

    IOdictionary d(io);
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

    IOdictionary d
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
    d.regIOobject::write();

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

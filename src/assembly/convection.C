/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
\*---------------------------------------------------------------------------*/

#include "convection.H"
#include "fvmDiv.H"
#include "fvcDiv.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::convection::convection(const fvMesh& mesh)
:
    mesh_(mesh),
    hoName_("div(phi,U)"),
    udName_("div(phi,U)_upwind")
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::fvVectorMatrix> Foam::convection::implicitPart
(
    const surfaceScalarField& phi,
    const volVectorField& U
) const
{
    return fvm::div(phi, U, udName_);
}


Foam::tmp<Foam::vectorField> Foam::convection::deferredCorrection
(
    const surfaceScalarField& phi,
    const volVectorField& U,
    const scalarField& beta
) const
{
    auto tcorr = tmp<vectorField>::New(mesh_.nCells(), Zero);

    // gMax is already a global (reduced) maximum: no further collective
    const bool anyHO = (gMax(beta) > 0);
    if (!anyHO)
    {
        return tcorr;
    }

    const tmp<volVectorField> tho = fvc::div(phi, U, hoName_);
    const tmp<volVectorField> tud = fvc::div(phi, U, udName_);

    vectorField& corr = tcorr.ref();
    const vectorField& ho = tho().primitiveField();
    const vectorField& ud = tud().primitiveField();

    forAll(corr, celli)
    {
        corr[celli] = beta[celli]*(ho[celli] - ud[celli]);
    }

    return tcorr;
}


// ************************************************************************* //

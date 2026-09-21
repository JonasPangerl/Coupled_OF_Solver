/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
\*---------------------------------------------------------------------------*/

#include "nonOrthCorrection.H"
#include "PstreamReduceOps.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::nonOrthCorrection::nonOrthCorrection
(
    const fvMesh& mesh,
    const scalar limiter,
    const scalar limiterStatic
)
:
    mesh_(mesh),
    limiter_(limiter),
    limiterStatic_(limiterStatic),
    staticFace_(mesh.nInternalFaces(), false),
    staticPatchFace_(mesh.boundary().size()),
    anyStatic_(false),
    orthogonal_(false)
{
    forAll(mesh.boundary(), patchi)
    {
        staticPatchFace_[patchi].resize(mesh.boundary()[patchi].size(), false);
    }

    // Exact test, no threshold: the native correction vectors vanish
    // identically on an orthogonal mesh
    const surfaceVectorField& k = mesh.nonOrthCorrectionVectors();
    scalar kMax = gMax(mag(k.primitiveField())());
    forAll(k.boundaryField(), patchi)
    {
        kMax = max(kMax, gMax(mag(k.boundaryField()[patchi])()));
    }
    reduce(kMax, maxOp<scalar>());
    orthogonal_ = (kMax == 0);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::nonOrthCorrection::setStaticCells(const boolList& isStatic)
{
    const labelUList& own = mesh_.owner();
    const labelUList& nei = mesh_.neighbour();

    bool any = false;

    forAll(staticFace_, facei)
    {
        staticFace_[facei] = isStatic[own[facei]] || isStatic[nei[facei]];
        any = any || staticFace_[facei];
    }

    forAll(mesh_.boundary(), patchi)
    {
        const labelUList& fc = mesh_.boundary()[patchi].faceCells();
        forAll(fc, pf)
        {
            staticPatchFace_[patchi][pf] = isStatic[fc[pf]];
            any = any || staticPatchFace_[patchi][pf];
        }
    }

    anyStatic_ = returnReduceOr(any);
}


// ************************************************************************* //

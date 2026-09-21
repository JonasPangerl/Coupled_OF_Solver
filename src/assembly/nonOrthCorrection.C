/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
\*---------------------------------------------------------------------------*/

#include "nonOrthCorrection.H"
#include "coupledDefaults.H"
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

    // Largest native correction vector (dimensionless, |k| = tan of the
    // non-orthogonality angle) over the internal and all patch faces of
    // this rank, then ONE global reduction. No collective may sit inside the
    // patch loop: the number of patches differs between ranks (processor
    // patches), and a collective called a rank-dependent number of times
    // shifts the collective sequence - the next collective with a different
    // message size then fails with MPI_ERR_TRUNCATE (T0 at 4 ranks).
    const surfaceVectorField& k = mesh.nonOrthCorrectionVectors();
    scalar kMax = 0;
    for (const vector& kf : k.primitiveField())
    {
        kMax = max(kMax, mag(kf));
    }
    for (const fvsPatchVectorField& kp : k.boundaryField())
    {
        for (const vector& kf : kp)
        {
            kMax = max(kMax, mag(kf));
        }
    }
    reduce(kMax, maxOp<scalar>());

    // Exact test kept deliberately (FABLE_REVIEW.md item 1): the tolerance
    // coupledDefaults::orthogonalityTolerance switches the correction off on
    // the T0 cavity (round-off 7.1e-14) and, unexpectedly, changes T0
    // Re 1000 np1 from 58 to 106 outer iterations - to be understood first
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

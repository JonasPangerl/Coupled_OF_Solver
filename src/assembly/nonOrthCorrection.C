/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
\*---------------------------------------------------------------------------*/

#include "nonOrthCorrection.H"
#include "coupledDefaults.H"
#include "PstreamReduceOps.H"
#include "syncTools.H"
#include "linear.H"
#include "snGradScheme.H"

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

void Foam::nonOrthCorrection::setLimitedCells(const boolList& isStatic)
{
    const labelUList& own = mesh_.owner();
    const labelUList& nei = mesh_.neighbour();

    bool any = false;

    forAll(staticFace_, facei)
    {
        staticFace_[facei] = isStatic[own[facei]] || isStatic[nei[facei]];
        any = any || staticFace_[facei];
    }

    // Static flag of the cell on the other side of every boundary face
    // (processor and cyclic patches; own cell elsewhere), so that coupled
    // faces use the same rule as internal faces - static if EITHER cell is
    // static - and both ranks of a processor face compute the same
    // correction. Collective: called on all ranks.
    labelList staticLabel(isStatic.size());
    forAll(isStatic, celli)
    {
        staticLabel[celli] = isStatic[celli] ? 1 : 0;
    }
    const labelList nbrStatic
    (
        syncTools::swapBoundaryCellList(mesh_, staticLabel)
    );

    const label nInternal = mesh_.nInternalFaces();
    forAll(mesh_.boundary(), patchi)
    {
        const fvPatch& patch = mesh_.boundary()[patchi];
        const labelUList& fc = patch.faceCells();
        const label start = patch.start() - nInternal;
        const bool coupled = patch.coupled();
        forAll(fc, pf)
        {
            staticPatchFace_[patchi][pf] =
                isStatic[fc[pf]]
             || (coupled && nbrStatic[start + pf] != 0);
            any = any || staticPatchFace_[patchi][pf];
        }
    }

    anyStatic_ = returnReduceOr(any);
}


Foam::tmp<Foam::surfaceScalarField>
Foam::nonOrthCorrection::correctionFromGrad
(
    const volScalarField& vf,
    const volVectorField& gradVf
) const
{
    if (orthogonal_)
    {
        return correction(vf);
    }

    // Full correction k_f . (grad vf)_f, exactly as native
    // correctedSnGrad::fullGradCorrection but with the given cell gradient
    tmp<surfaceScalarField> tc =
        linear<vector>(mesh_).dotInterpolate
        (
            mesh_.nonOrthCorrectionVectors(),
            gradVf
        );
    surfaceScalarField& c = tc.ref();

    // Uncorrected snGrad with the corrected-scheme deltaCoeffs, as native
    // limitedSnGrad::correction
    const tmp<surfaceScalarField> tsn =
        fv::snGradScheme<scalar>::snGrad
        (
            vf,
            tmp<surfaceScalarField>(mesh_.nonOrthDeltaCoeffs()),
            "SndGrad"
        );
    const surfaceScalarField& sn = tsn();

    // Native limitedSnGrad limiter, lambda per face (static faces use the
    // static limiter)
    auto limit = [](const scalar lambda, const scalar snf, scalar& cf)
    {
        // GUARD: denominator >= SMALL, as native limitedSnGrad
        const scalar l =
            min
            (
                lambda*mag(snf)/((1 - lambda)*mag(cf) + SMALL),
                scalar(1)
            );
        cf *= l;
    };

    scalarField& ci = c.primitiveFieldRef();
    const scalarField& sni = sn.primitiveField();
    forAll(ci, facei)
    {
        const bool st = anyStatic_ && staticFace_[facei];
        limit(st ? limiterStatic_ : limiter_, sni[facei], ci[facei]);
    }

    auto& cbf = c.boundaryFieldRef();
    forAll(cbf, patchi)
    {
        scalarField& cp = cbf[patchi];
        const scalarField& snp = sn.boundaryField()[patchi];
        const boolList& sp = staticPatchFace_[patchi];
        forAll(cp, pf)
        {
            const bool st = anyStatic_ && sp[pf];
            limit(st ? limiterStatic_ : limiter_, snp[pf], cp[pf]);
        }
    }

    return tc;
}


// ************************************************************************* //

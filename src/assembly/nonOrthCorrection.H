/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
-------------------------------------------------------------------------------
Class
    Foam::nonOrthCorrection

Description
    Explicit non-orthogonal correction of face-normal gradients with a
    per-face limiter (spec 5.3a/d nonOrthLimiter, 8.1 static nonOrthLimiter,
    DECISIONS.md D-014).

    The correction is the native limitedSnGrad correction
    ("limited corrected <lambda>"), evaluated for the global limiter and for
    the static-remediation limiter; faces adjacent to a static cell take the
    latter. A face counts as static if EITHER of its cells is static; on
    processor (and cyclic) faces the static flag of the neighbour cell is
    exchanged (syncTools::swapBoundaryCellList), so both ranks of a
    processor face apply the same limiter and compute the same correction.

    correctionFromGrad(p, gradp) evaluates the same limited correction from a
    given cell gradient instead of the grad(p) scheme: the Rhie-Chow terms
    pass the gradient that is cellLimited in static cells (D-018), so the
    non-orthogonal part of q_f uses it too.

    The orthogonal-mesh shortcut is taken only if max |k_f| == 0 exactly
    (all correction vectors bit-zero). A blockMesh cavity is NOT exactly
    orthogonal in double (max |k_f| = 7.1e-14 on T0), so the correction is
    evaluated there, with values at round-off level (FABLE_REVIEW.md item 1).

    Returned per face: corr_f, such that the full face gradient is
        snGrad_f = deltaCoeff_f (phi_N - phi_P) + corr_f.

SourceFiles
    nonOrthCorrection.C

\*---------------------------------------------------------------------------*/

#ifndef coupledFoam_nonOrthCorrection_H
#define coupledFoam_nonOrthCorrection_H

#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class nonOrthCorrection
{
    // Private Data

        const fvMesh& mesh_;

        //- Global limiter
        scalar limiter_;

        //- Limiter on faces adjacent to static-remediation cells
        scalar limiterStatic_;

        //- Internal faces adjacent to a static cell
        boolList staticFace_;

        //- Boundary faces (per patch) whose owner is a static cell, or,
        //  on coupled patches, whose neighbour cell is static
        List<boolList> staticPatchFace_;

        //- Any static face at all (global)
        bool anyStatic_;

        //- All native correction vectors are exactly zero (then the
        //  correction is identically zero and is skipped)
        bool orthogonal_;


    // Private Member Functions

        template<class Type>
        tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> limitedCorrection
        (
            const GeometricField<Type, fvPatchField, volMesh>& vf,
            const scalar lambda
        ) const;


public:

    // Constructors

        nonOrthCorrection
        (
            const fvMesh& mesh,
            const scalar limiter,
            const scalar limiterStatic
        );


    // Member Functions

        //- Set the static cells (from remediation, 8.1). Collective (swaps
        //  the static flag across coupled patches)
        void setStaticCells(const boolList& isStatic);

        //- True if all native correction vectors are exactly zero
        bool orthogonal() const noexcept
        {
            return orthogonal_;
        }

        scalar limiter() const noexcept { return limiter_; }
        scalar limiterStatic() const noexcept { return limiterStatic_; }

        //- Per-face limited correction of snGrad(vf)
        template<class Type>
        tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> correction
        (
            const GeometricField<Type, fvPatchField, volMesh>& vf
        ) const;

        //- Per-face limited correction of snGrad(vf) with the full
        //  correction k_f . (gradVf)_f taken from the given cell gradient
        //  (linear interpolation, coupled boundary values of gradVf must be
        //  up to date); same limiter as correction()
        tmp<surfaceScalarField> correctionFromGrad
        (
            const volScalarField& vf,
            const volVectorField& gradVf
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "nonOrthCorrectionTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

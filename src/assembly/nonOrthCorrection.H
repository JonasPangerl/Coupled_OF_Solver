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
    latter. On processor faces the owner cell decides.

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

        //- Boundary faces (per patch) whose owner is a static cell
        List<boolList> staticPatchFace_;

        //- Any static face at all (global)
        bool anyStatic_;

        //- Mesh is orthogonal (correction identically zero)
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

        //- Set the static cells (from remediation, 8.1)
        void setStaticCells(const boolList& isStatic);

        //- True if the mesh has no non-orthogonality (all corrections 0)
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

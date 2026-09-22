/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
-------------------------------------------------------------------------------
Class
    Foam::coupledAssembler

Description
    Builds the 4x4-block system of one outer iteration (spec 5.3) and the
    right-hand side of the increment solve A dx = b - A x (5.7, D-012).

    Two stages, because the solution-limited local time step (7.3) needs the
    momentum residual before the PTC term, and D = V/abar (5.3d) needs the
    PTC term:

      assembleMomentum()   rows 0-2: native fvVectorMatrix coefficients
                           (implicit upwind convection, uncorrected
                           Laplacian), explicit sources (beta-weighted
                           deferred correction, transpose stress, limited
                           non-orthogonal correction), implicit MRF, pressure
                           gradient column; momentumResidual() available
      assembleContinuity(rDeltaTV)
                           PTC diagonal V/dt, abar, D, D_f, row 3 (divergence,
                           Rhie-Chow Laplacian, explicit RC source), pressure
                           reference, continuity row scaling (5.3f, applied
                           in double before the narrowing to blockScalar),
                           residual and norms

    Diagonal blocks are accumulated in double; every off-diagonal block is
    complete after its own face and narrowed once. The residual b - A x is
    accumulated in double from the same double coefficients (D-012).

    Norms (after row scaling), see DECISIONS.md D-016:
        normFactor = sum_i |(A x)_i| + |b_i| + SMALL  (all rows except the
                                                       replaced row of the
                                                       pressure reference,
                                                       D-021)
        residualL2 = ||b - A x||_2 / normFactor       (the R of 7.1 before
                                                       normalisation by R_1)
        rU, rp     = sum |r_i| / sum(|Ax_i|+|b_i|)   over rows 0-2 / row 3
                     (L1 like native initial residuals, for the CF| line)

SourceFiles
    coupledAssembler.C

\*---------------------------------------------------------------------------*/

#ifndef coupledFoam_coupledAssembler_H
#define coupledFoam_coupledAssembler_H

#include "blockLduMatrix4.H"
#include "boundaryCoupling.H"
#include "convection.H"
#include "nonOrthCorrection.H"
#include "rhieChow.H"
#include "diagnostics.H"
#include "MRFCoupling.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class coupledAssembler
{
public:

    //- Wall-time split of the assembly (diagnostics, TASK 5), seconds,
    //  accumulated until resetTimes()
    struct timings
    {
        doubleScalar momentumOps = 0;   //!< stage 1 without boundary faces
        doubleScalar boundary = 0;      //!< boundary faces, stage 1 and 2
        doubleScalar continuity = 0;    //!< stage 2 without RC, boundary
        doubleScalar rhieChow = 0;      //!< D, D_f and explicit face term
    };


private:

    // Private Data

        //- Timing of the stages on (diagnostics level >= 1 only)
        bool timing_;
        timings times_;

        const fvMesh& mesh_;

        //- Patch classification
        const List<boundaryCoupling::patchKind> kinds_;

        //- Block interface index per patch (-1 if none)
        labelList blockIndex_;

        //- Convection (deferred correction)
        convection conv_;

        //- Non-orthogonal correction with per-face limiter
        nonOrthCorrection noc_;

        //- Rhie-Chow interpolation
        rhieChow rc_;

        //- Implicit MRF (may be null)
        const MRFCoupling* mrfPtr_;

        //- Clamp for the narrowing to blockScalar (9.2)
        const doubleScalar clampValue_;

        //- Pressure reference
        bool needRef_;
        label pRefCell_;
        scalar pRefValue_;

        //- The one-time boundary flux balance check was made (A10)
        bool refFluxChecked_;

        //- The block matrix (float) and the increment right-hand side
        blockLduMatrix4 A_;
        blockScalarList rhs_;

        //- Double work storage
        reduceScalarList Dd_;
        reduceScalarList Ax_;
        reduceScalarList b_;

        //- Momentum residual per cell (rows 0-2, without PTC)
        vectorField rMom_;

        //- Mean momentum diagonal (incl. PTC)
        scalarField abar_;

        //- Fields of the current assembly
        const volVectorField* UPtr_;
        const volScalarField* pPtr_;

        //- Results
        scalar sp_;
        reduceScalar normFactor_;
        reduceScalar residualL2_;
        reduceScalar rU_;
        reduceScalar rp_;
        label nClamped_;


    // Private Member Functions

        //- Index into Dd_ for cell, row, column
        static label di(const label celli, const label r, const label c)
        {
            return celli*blockSize + r*blockDim + c;
        }

        //- Narrow with clamping and counting
        blockScalar store(const doubleScalar v)
        {
            return toBlock(v, clampValue_, nClamped_);
        }

        //- True if the patch carries no flux by construction
        //  (symmetry, symmetryPlane). Wedge is not listed: its native
        //  valueInternal/BoundaryCoeffs linearisation is exact and is what
        //  the flux update (rhieChow::updateFlux) evaluates too.
        static bool zeroFluxPatch(const fvPatch& patch);

        //- One-time Warning if the strong pressure reference is active and
        //  the prescribed boundary fluxes do not sum to ~0 (collective)
        void checkReferenceFluxBalance(const volVectorField& U);


public:

    // Constructors

        //- Construct from mesh, the "coupled" dictionary, p (for the
        //  reference cell) and an optional MRF coupling
        coupledAssembler
        (
            const fvMesh& mesh,
            const dictionary& coupledDict,
            const volScalarField& p,
            const MRFCoupling* mrfPtr
        );


    // Member Functions

        // Set-up

            //- Static remediation cells (8.1): per-face limiter selection
            void setStaticCells(const boolList& isStatic);


        // Assembly

            //- Stage 1: momentum rows
            void assembleMomentum
            (
                const volVectorField& U,
                const volScalarField& p,
                const surfaceScalarField& phi,
                const volScalarField& nuEff,
                const scalarField& beta
            );

            //- Momentum residual per cell (b - A x, rows 0-2, no PTC)
            const vectorField& momentumResidual() const noexcept
            {
                return rMom_;
            }

            //- Stage 2: PTC, D, continuity row, scaling, residual
            //  rDeltaTV = V_P/dt_P per cell (after the local limit)
            void assembleContinuity(const scalarField& rDeltaTV);


        // Diagnostics

            //- Switch the stage timing on (no clock reads when off)
            void setTiming(const bool on) noexcept { timing_ = on; }
            const timings& times() const noexcept { return times_; }
            void resetTimes() noexcept { times_ = timings(); }


        // Access

            const blockLduMatrix4& matrix() const noexcept { return A_; }
            const blockScalarList& rhs() const noexcept { return rhs_; }

            rhieChow& rc() noexcept { return rc_; }
            const rhieChow& rc() const noexcept { return rc_; }
            const nonOrthCorrection& noc() const noexcept { return noc_; }

            const scalarField& abar() const noexcept { return abar_; }

            //- Continuity row scale s_p (5.3f)
            scalar rowScale() const noexcept { return sp_; }

            reduceScalar normFactor() const noexcept { return normFactor_; }

            //- ||b - A x||_2 / normFactor (scaled system)
            reduceScalar residualL2() const noexcept { return residualL2_; }

            reduceScalar rU() const noexcept { return rU_; }
            reduceScalar rp() const noexcept { return rp_; }

            //- Clamped coefficients of the last assembly (this rank)
            label nClamped() const noexcept { return nClamped_; }

            const List<boundaryCoupling::patchKind>& kinds() const noexcept
            {
                return kinds_;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
-------------------------------------------------------------------------------
Class
    Foam::blockGAMG

Description
    Additive-correction block algebraic multigrid (spec 6.3), used as one
    V-cycle per Krylov iteration by the blockGAMG preconditioner.

    Agglomeration: native GAMGAgglomeration::New(mesh, dict) (faceAreaPair
    by default; geometric, so the block coefficients are not needed). Each
    mergeLevels attempt uses its own agglomeration "name" so that a retry
    creates a new cached object instead of returning the previous one.

    Coarse operators: Galerkin by summation, in double, stored as
    blockScalar:
        coarse diag  = sum fine diag of the agglomerated cells
                     + sum fine upper+lower of faces internal to a coarse cell
        coarse upper/lower = sum fine upper/lower of surviving faces
                     (swapped where native faceFlipMap is set)
        coarse interface coeffs = sum over patchFaceRestrictAddressing
    Restriction: sum of residual blocks. Prolongation: piecewise constant.

    Operator complexity C_op = sum_levels nnzBlocks(level)/nnzBlocks(fine)
    (global sums). If C_op > maxOperatorComplexity, mergeLevels is increased
    by 1 and the agglomeration repeated, at most maxCopAttempts times, then
    FatalError.

    V-cycle: nPreSweeps on every level going down (starting from x = 0),
    coarsest level solved with blockBiCGStab + blockDiagonal to
    coarsestTolerance (relative), nPostSweeps going up, nFinestSweeps on the
    finest level.

SourceFiles
    blockGAMG.C

\*---------------------------------------------------------------------------*/

#ifndef coupledFoam_blockGAMG_H
#define coupledFoam_blockGAMG_H

#include "blockLduMatrix4.H"
#include "blockSmoother.H"
#include "blockSolver.H"
#include "GAMGAgglomeration.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class blockGAMG
{
    // Private Data

        //- Finest-level matrix
        const blockLduMatrix4& fine_;

        //- blockGAMG controls
        dictionary dict_;

        label nPreSweeps_;
        label nPostSweeps_;
        label nFinestSweeps_;
        reduceScalar maxCop_;
        label maxCopAttempts_;

        //- Agglomeration in use (cached on the mesh database)
        const GAMGAgglomeration* aggPtr_;

        //- mergeLevels actually used
        label mergeLevelsUsed_;

        //- Operator complexity of the accepted hierarchy
        reduceScalar Cop_;

        //- Coarse matrices, index l = level l+1
        PtrList<blockLduMatrix4> coarse_;

        //- Smoothers, index l = level l (0 = finest)
        PtrList<blockSmoother> smoothers_;

        //- Coarsest-level solver
        autoPtr<blockSolver> coarsestSolver_;

        //- Work vectors per level (x, b, r)
        mutable List<blockScalarList> x_;
        mutable List<blockScalarList> b_;
        mutable List<blockScalarList> r_;

        //- Last coarsest-level iteration count (diagnostics)
        mutable label nCoarsestIters_;


    // Private Member Functions

        //- Choose the agglomeration respecting the C_op cap
        void agglomerate();

        //- Operator complexity of an agglomeration (global)
        reduceScalar operatorComplexity(const GAMGAgglomeration& agg) const;

        //- Build coarse level fineLevel+1 from fineLevel (coefficients only)
        void restrictMatrix(const label fineLevel);

        //- Matrix of level l (0 = finest)
        const blockLduMatrix4& matrixLevel(const label l) const
        {
            return (l == 0 ? fine_ : coarse_[l - 1]);
        }

        //- Restrict r of level l into b of level l+1 (sum)
        void restrictResidual(const label l) const;


public:

    // Constructors

        //- Construct for the finest matrix (must have an lduMesh) and the
        //  blockGAMG sub-dictionary
        blockGAMG(const blockLduMatrix4& fine, const dictionary& dict);

        //- No copy
        blockGAMG(const blockGAMG&) = delete;
        void operator=(const blockGAMG&) = delete;


    // Member Functions

        //- Number of levels including the finest
        label nLevels() const
        {
            return coarse_.size() + 1;
        }

        //- Operator complexity
        reduceScalar operatorComplexity() const noexcept
        {
            return Cop_;
        }

        //- mergeLevels used
        label mergeLevels() const noexcept
        {
            return mergeLevelsUsed_;
        }

        //- Global cell count per level
        labelList globalCellsPerLevel() const;

        //- Guarded pivots summed over all levels at the last update (rank)
        label nSingularDiag() const;

        //- Coarsest-level iterations of the last V-cycle
        label nCoarsestIters() const noexcept
        {
            return nCoarsestIters_;
        }

        //- Recompute coarse operators and smoother factorisations
        void update();

        //- One V-cycle: x = approx A^-1 b (x is overwritten, start from 0)
        void Vcycle(blockScalarUList& x, const blockScalarUList& b) const;

        //- Write hierarchy statistics (levels, cells, C_op) once per run
        void writeStats(Ostream& os) const;

        //- Statistics as a dictionary (for JSON output)
        dictionary statsDict() const;

        //- Effective settings
        void writeSettings(dictionary& dict) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
-------------------------------------------------------------------------------
Class
    Foam::blockGAMG

Description
    Additive-correction block algebraic multigrid (spec 6.3, amendment B3),
    used as the preconditioner of the outer Krylov solver.

    Agglomeration: native GAMGAgglomeration::New(mesh, dict) (faceAreaPair
    by default). Each mergeLevels attempt uses its own agglomeration "name"
    so that a retry creates a new cached object.

    Acceptance of an agglomeration (6.3, 6.3.1), mergeLevels from the
    dictionary value up to maxMergeLevels (4):
      - operator complexity C_op = sum_l nnzBlocks(l)/nnzBlocks(0)
        <= maxOperatorComplexity
      - for cycleType W or K: every measured coarsening ratio
        r_l = nCells(l)/nCells(l+1) >= 3 (global counts)
    otherwise FatalError with the measured values.

    Coarse operators: Galerkin by summation in double, stored as
    blockScalar; restriction: block sum; prolongation: piecewise constant.

    Cycles (6.3.2), level 0 finest, L coarsest:
      V   one recursive coarse correction
      F   F-cycle: coarse correction by an F-cycle followed by a V-cycle on
          the remaining residual
      W   two recursive W-corrections, the second on the updated residual
      K   Notay-Vassilevski K-cycle: on every coarse level the correction is
          one or two GCR steps preconditioned by the recursive cycle; the
          second step only if the first reduced the residual by less than
          kCycleThreshold. Variable preconditioner: requires blockFGMRES.
    Pre-smoothing nPreSweeps, post-smoothing nPostSweeps, on the finest
    level nFinestSweeps for both.

    Coarsest level (6.3.3): blockBiCGStab + blockDiagonal to
    coarsestTolerance / coarsestMaxIter; if the coarsest level lives on one
    rank with at most denseLUMaxCells (64) cells, a dense LU in double of the
    assembled (4n)x(4n) matrix, factorised once per matrix update.

    The cycle type and the sweep counts can be changed between solves by the
    autoTune controller (6.3.5); a promotion to W or K is refused if the
    coarsening ratios do not satisfy the rule of 6.3.1.

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
public:

    //- Cycle types (6.3.2)
    enum class cycleKind { V, F, W, K };

    static cycleKind cycleFromWord(const word& w);
    static word cycleName(const cycleKind c);

private:

    // Private Data

        //- Finest-level matrix
        const blockLduMatrix4& fine_;

        //- blockGAMG controls (effective values)
        dictionary dict_;

        mutable cycleKind cycle_;
        reduceScalar kThreshold_;
        label kMaxSteps_;
        label nPreSweeps_;
        mutable label nPostSweeps_;
        label nFinestSweeps_;
        reduceScalar maxCop_;
        reduceScalar minRatio_;
        label maxMergeLevels_;
        label denseLUMaxCells_;

        //- Agglomeration in use (cached on the mesh database)
        const GAMGAgglomeration* aggPtr_;

        //- mergeLevels actually used
        label mergeLevelsUsed_;

        //- Operator complexity and coarsening ratios of the hierarchy
        reduceScalar Cop_;
        List<reduceScalar> ratios_;

        //- Coarse matrices, index l = level l+1
        PtrList<blockLduMatrix4> coarse_;

        //- Smoothers, index l = level l (0 = finest)
        PtrList<blockSmoother> smoothers_;

        //- Coarsest-level Krylov solver
        autoPtr<blockSolver> coarsestSolver_;

        //- Dense LU of the coarsest level (if selected), row-major, double
        bool useDenseLU_;
        reduceScalarList denseLU_;
        labelList densePivot_;

        //- Work vectors per level
        mutable List<blockScalarList> b_;     // restricted right-hand side
        mutable List<blockScalarList> e_;     // correction (cycle output)
        mutable List<blockScalarList> r_;     // residual scratch
        mutable List<blockScalarList> t_;     // second right-hand side (F, W)
        mutable List<blockScalarList> e2_;    // second correction (F, W)
        mutable List<blockScalarList> z1_;    // K-step workspace
        mutable List<blockScalarList> q1_;
        mutable List<blockScalarList> z2_;
        mutable List<blockScalarList> q2_;
        mutable List<blockScalarList> rk_;

        //- Last coarsest-level iteration count (diagnostics)
        mutable label nCoarsestIters_;


    // Private Member Functions

        //- Choose the agglomeration respecting C_op and the ratio rule
        void agglomerate();

        //- Operator complexity and ratios of an agglomeration (global)
        void measure
        (
            const GAMGAgglomeration& agg,
            reduceScalar& cop,
            List<reduceScalar>& ratios
        ) const;

        //- True if the ratio rule holds
        bool ratiosOk(const List<reduceScalar>& ratios) const;

        //- True if cycle kind needs the ratio rule
        static bool needsRatio(const cycleKind c)
        {
            return c == cycleKind::W || c == cycleKind::K;
        }

        //- Build coarse level fineLevel+1 from fineLevel (coefficients)
        void restrictMatrix(const label fineLevel);

        //- Matrix of level l (0 = finest)
        const blockLduMatrix4& matrixLevel(const label l) const
        {
            return (l == 0 ? fine_ : coarse_[l - 1]);
        }

        //- Coarsest level index
        label L() const
        {
            return coarse_.size();
        }

        //- Allocate per-level work vectors
        void allocateWork();

        //- Factorise the coarsest level densely (if selected)
        void factoriseDense();

        //- Solve the coarsest level
        void solveCoarsest(blockScalarUList& x, const blockScalarUList& b) const;

        //- Restrict r (level l) into bc (level l+1), block sum
        void restrictVector
        (
            const label l,
            const blockScalarUList& r,
            blockScalarUList& bc
        ) const;

        //- x (level l) += P e (level l+1)
        void prolongAdd
        (
            const label l,
            const blockScalarUList& e,
            blockScalarUList& x
        ) const;

        //- One cycle of the given kind at level l: x ~ A_l^-1 b, x from 0
        void cycle
        (
            const label l,
            const blockScalarUList& b,
            blockScalarUList& x,
            const cycleKind kind
        ) const;

        //- K-step at level l (GCR, 1 or 2 steps), output in e
        void kstep
        (
            const label l,
            const blockScalarUList& b,
            blockScalarUList& e
        ) const;


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

        reduceScalar operatorComplexity() const noexcept { return Cop_; }
        label mergeLevels() const noexcept { return mergeLevelsUsed_; }
        const List<reduceScalar>& ratios() const noexcept { return ratios_; }
        bool denseCoarsest() const noexcept { return useDenseLU_; }

        //- Global cell count per level
        labelList globalCellsPerLevel() const;

        //- Guarded pivots summed over all levels at the last update (rank)
        label nSingularDiag() const;

        //- Coarsest-level iterations of the last cycle
        label nCoarsestIters() const noexcept { return nCoarsestIters_; }

        //- Current cycle type and post-sweeps (autoTune)
        cycleKind cycleType() const noexcept { return cycle_; }
        label nPostSweeps() const noexcept { return nPostSweeps_; }

        //- Change the cycle type between solves. Returns false (and changes
        //  nothing) if W/K would violate the coarsening-ratio rule.
        bool setCycleType(const cycleKind c) const;

        //- Change the post-sweeps between solves
        void setNPostSweeps(const label n) const
        {
            nPostSweeps_ = n;
        }

        //- Recompute coarse operators, smoother factorisations, dense LU
        void update();

        //- One cycle: x = approx A^-1 b (x is overwritten)
        void apply(blockScalarUList& x, const blockScalarUList& b) const;

        //- Write hierarchy statistics (levels, cells, ratios, C_op)
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

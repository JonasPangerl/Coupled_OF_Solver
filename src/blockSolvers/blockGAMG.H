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
    Matrix-weighted pair agglomeration (D-039), keyword agglomerationWeights
    geometric | momentum | pressure | combined, rebuilt every
    reagglomerateInterval matrix updates. "agglomerator algebraicPair"
    (amendment C4) is the spec name of agglomerationWeights pressure: face
    weights = p-p block (row 3, col 3) magnitudes max(|U_pp|, |L_pp|), so
    the agglomeration follows the anisotropy of the Rhie-Chow Laplacian;
    before the first matrix exists the hierarchy is faceAreaPair. An
    explicit agglomerationWeights other than pressure together with
    algebraicPair is a FatalIOError.

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

    Processor agglomeration (6.3.4), keyword processorAgglomerator:
      masterCoarsest (default)
          the rule of 6.3.4: a level l >= 1 with N_l global cells lives on
          nRanks ranks if N_l >= procAgglomCellsPerRank*nRanks, on
          max(1, nRanks/procAgglomDivisor) ranks if N_l is below that, and
          on one rank if N_l < procAgglomCellsPerRank (5000, 4).
          Levels 1..L-1 are gathered by the native GAMGProcAgglomeration
          mechanism (runtime type blockGAMGRule, see
          blockGAMGProcAgglomeration.H); the coarsest level L, which the
          native mechanism cannot gather on its own, is gathered here with
          the native lduPrimitiveMesh gathering (gatherCoarsest).
          Level 0 (the decomposed mesh) is never agglomerated.
      none
          no processor agglomeration
      nativeMasterCoarsest
          the native masterCoarsest agglomerator (its own keywords)
      any other word
          passed verbatim to the native selector (manual, procFaces, eager)
          together with the blockGAMG dictionary.
    On a processor-agglomerated level l the ranks of a group restrict their
    level-(l-1) operator locally and send it to the group master, which
    assembles it with the native maps cellOffsets, faceMap, boundaryMap and
    boundaryFaceMap (native GAMGSolver::procAgglomerateMatrix), 16
    coefficients per block; processor faces that become internal go to
    upper (owner side, map >= 0) or lower (neighbour side, map < 0).
    Vectors are gathered/scattered with globalIndex::gather/scatter on the
    group communicator. Ranks without the level skip all work on it and
    below but take part in the gather and the scatter. Every reduction on a
    level uses that level's communicator.

SourceFiles
    blockGAMG.C

\*---------------------------------------------------------------------------*/

#ifndef coupledFoam_blockGAMG_H
#define coupledFoam_blockGAMG_H

#include "blockLduMatrix4.H"
#include "blockSmoother.H"
#include "blockSolver.H"
#include "diagnostics.H"
#include "GAMGAgglomeration.H"
#include "UPtrList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class lduPrimitiveMesh;

class blockGAMG
{
public:

    //- Cycle types (6.3.2)
    enum class cycleKind { V, F, W, K };

    static cycleKind cycleFromWord(const word& w);
    static word cycleName(const cycleKind c);


private:

    // Private Types

        //- Gathering of a processor-agglomerated coarse level from the
        //  local (pre-gather) coarse operators of a rank group (6.3.4)
        struct gatherLevel
        {
            //- Group communicator, rank 0 = group master
            label comm = -1;

            //- True on the group master
            bool master = false;

            //- identity(nProcs(comm)) (procIDs for globalIndex)
            labelList procIDs;

            //- Native maps, master only (index = rank in comm)
            labelList cellOffsets;
            labelList vecOffsets;       // blockDim*cellOffsets
            labelListList faceMap;
            labelListList boundaryMap;
            labelListListList boundaryFaceMap;
        };


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

        //- Coarse-correction scaling: 0 none, 1 finest level, 2 all levels
        //  (minimal-residual factor alpha = <A d, r>/<A d, A d> for the
        //  prolongated correction d, as native GAMG scaleCorrection)
        label scaleMode_;

        //- Agglomeration weights: geometric (native agglomerator) or from
        //  the finest block matrix (momentum | pressure | combined).
        //  "agglomerator algebraicPair" (amendment C4) selects pressure.
        word aggWeights_;
        label reaggInterval_;
        label nUpdates_;
        bool aggFromMatrix_;

        //- Processor agglomeration controls (6.3.4)
        word procAgglomType_;
        label procAgglomCellsPerRank_;
        label procAgglomDivisor_;

        //- Agglomeration in use (cached on the mesh database)
        const GAMGAgglomeration* aggPtr_;

        //- mergeLevels actually used
        label mergeLevelsUsed_;

        //- Operator complexity and coarsening ratios of the hierarchy
        reduceScalar Cop_;
        List<reduceScalar> ratios_;

        //- Coarsest level index on this rank (= agglomeration size)
        label L_;

        //- Global cells and ranks per level (all ranks, cached)
        labelList cellsPerLevel_;
        labelList ranksPerLevel_;

        //- Mesh of each level held by this rank (unset: not held)
        UPtrList<const lduMesh> levelMesh_;

        //- Gathering info, index = coarse level (set where the level is
        //  processor-agglomerated and this rank holds the level below)
        PtrList<gatherLevel> gather_;

        //- Coarsest level gathered by blockGAMG (master only) and its
        //  communicators (-1 if unused)
        autoPtr<lduPrimitiveMesh> coarsestMesh_;
        label coarsestAllComm_;
        label coarsestAgglomComm_;

        //- Coarse matrices, index l = level l+1 (unset: level not held)
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
        mutable List<blockScalarList> gbuf_;  // local pre-gather coarse vector
        mutable List<blockScalarList> sd_;    // scaled correction P e
        mutable List<blockScalarList> sq_;    // A P e

        //- Last coarsest-level iteration count (diagnostics)
        mutable label nCoarsestIters_;

        //- Incremented on every (re)build of the hierarchy (agglomeration,
        //  gathering, coarse levels)
        label hierarchyVersion_;

        //- Deep diagnostics (TASK 5; null: off), set by the preconditioner
        mutable diagnostics* diag_;


    // Private Member Functions

        //- Build agglomeration, gathering, matrices, smoothers, work
        void buildHierarchy();

        //- Release everything built by buildHierarchy (not the agglomeration)
        void clearHierarchy();

        //- Choose the agglomeration respecting C_op and the ratio rule
        void agglomerate();

        //- Operator complexity, ratios and global cells per level of an
        //  agglomeration (collective on the finest communicator)
        void measure
        (
            const GAMGAgglomeration& agg,
            reduceScalar& cop,
            List<reduceScalar>& ratios,
            labelList& cells
        ) const;

        //- Face weights of the finest level from the block matrix
        tmp<scalarField> matrixFaceWeights() const;

        //- True if the coarse correction at level l is scaled
        bool scaleAt(const label l) const
        {
            return scaleMode_ == 2 || (scaleMode_ == 1 && l == 0);
        }

        //- True if the ratio rule holds
        bool ratiosOk(const List<reduceScalar>& ratios) const;

        //- True if cycle kind needs the ratio rule
        static bool needsRatio(const cycleKind c)
        {
            return c == cycleKind::W || c == cycleKind::K;
        }

        //- Rule 6.3.4 active (processorAgglomerator masterCoarsest)
        bool ruleMode() const
        {
            return procAgglomType_ == "masterCoarsest";
        }

        //- Set up gathering info of the natively agglomerated levels
        void setNativeGathering();

        //- Gather the coarsest level if rule 6.3.4 requires fewer ranks
        //  on it than on the level above (native lduPrimitiveMesh gather)
        void gatherCoarsest();

        //- True if this rank holds level l
        bool hasLevel(const label l) const
        {
            return l >= 0 && l <= L_ && levelMesh_.set(l);
        }

        //- Build coarse level fineLevel+1 from fineLevel (coefficients)
        void restrictMatrix(const label fineLevel);

        //- Assemble the gathered coarse matrix of level l on the master
        void gatherMatrix
        (
            const label l,
            const blockScalarList& lDiag,
            const blockScalarList& lUpper,
            const blockScalarList& lLower,
            const labelList& lIfaceIndex,
            const List<blockScalarList>& lIfaceCoeffs
        );

        //- Matrix of level l (0 = finest), level must be held
        const blockLduMatrix4& matrixLevel(const label l) const
        {
            return (l == 0 ? fine_ : coarse_[l - 1]);
        }

        //- Coarsest level index
        label L() const noexcept
        {
            return L_;
        }

        //- Allocate per-level work vectors
        void allocateWork();

        //- Allocate the K-step workspace
        void allocateKWork() const;

        //- Factorise the coarsest level densely (if selected)
        void factoriseDense();

        //- Solve the coarsest level
        void solveCoarsest(blockScalarUList& x, const blockScalarUList& b) const;

        //- Restrict r (level l) into bc (level l+1), block sum; gathers on
        //  a processor-agglomerated level l+1 (bc used on the master only)
        void restrictVector
        (
            const label l,
            const blockScalarUList& r,
            blockScalarUList& bc
        ) const;

        //- x (level l) += P e (level l+1); scatters on a
        //  processor-agglomerated level l+1 (e used on the master only)
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

        //- True if diagnostics at level n are on
        bool diagOn(const label n) const noexcept
        {
            return diag_ && diag_->active(n);
        }

        //- nSweeps single sweeps with the residual norm after each
        //  (diagnostics level 3; same result as one smooth(x, b, nSweeps)
        //  call: the sweeps of both smoothers are independent)
        void smoothLogged
        (
            const label l,
            blockScalarUList& x,
            const blockScalarUList& b,
            const label nSweeps,
            DynamicList<doubleScalar>& norms
        ) const;

        //- r_[l] = b - A_l x (the same float operations as
        //  blockLduMatrix4::residual) and, in the same pass, the local sums
        //  of squares of b and r (diagnostics; no reduction)
        void residualSums
        (
            const label l,
            const blockScalarUList& x,
            const blockScalarUList& b,
            reduceScalar& bb,
            reduceScalar& rr
        ) const;

        //- Residual norm ||b - A_l x|| via the scratch r_[l] (collective
        //  on the level communicator; diagnostics only)
        doubleScalar residualNorm
        (
            const label l,
            const blockScalarUList& x,
            const blockScalarUList& b
        ) const;

        //- Per-level operator statistics after an update (level 3,
        //  rank-local, no reduction)
        void diagOperatorStats() const;

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


    //- Destructor (frees the coarsest-level communicators)
    ~blockGAMG();


    // Member Functions

        //- Number of levels including the finest (global)
        label nLevels() const noexcept
        {
            return cellsPerLevel_.size();
        }

        reduceScalar operatorComplexity() const noexcept { return Cop_; }
        label mergeLevels() const noexcept { return mergeLevelsUsed_; }
        const List<reduceScalar>& ratios() const noexcept { return ratios_; }
        bool denseCoarsest() const noexcept { return useDenseLU_; }

        //- Global cell count per level (cached, not collective)
        const labelList& globalCellsPerLevel() const noexcept
        {
            return cellsPerLevel_;
        }

        //- Number of ranks holding each level (cached, not collective)
        const labelList& ranksPerLevel() const noexcept
        {
            return ranksPerLevel_;
        }

        //- Guarded pivots summed over all levels at the last update (rank)
        label nSingularDiag() const;

        //- Smoother pivot fallbacks (blockILU0 safeguard, D-049) summed
        //  over all levels at the last update (rank)
        label nPivotFallback() const;

        //- Coarsest-level iterations of the last cycle
        label nCoarsestIters() const noexcept { return nCoarsestIters_; }

        //- Hierarchy build counter (changes on every re-agglomeration)
        label hierarchyVersion() const noexcept { return hierarchyVersion_; }

        //- Number of coefficient updates so far
        label nUpdates() const noexcept { return nUpdates_; }

        //- Attach deep diagnostics (null: off)
        void setDiagnostics(diagnostics* d) const noexcept { diag_ = d; }

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

        //- Write hierarchy statistics (levels, cells, ranks, ratios, C_op)
        void writeStats(Ostream& os) const;

        //- Statistics as a dictionary (for JSON output, not collective)
        dictionary statsDict() const;

        //- Effective settings
        void writeSettings(dictionary& dict) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

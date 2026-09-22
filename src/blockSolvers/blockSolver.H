/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
-------------------------------------------------------------------------------
Class
    Foam::blockSolver

Description
    Base class of the block Krylov solvers with runtime selection
    (spec 6.1, keyword "solver": blockBiCGStab | blockGMRES).

    Residual norm (spec 6.1): ||r|| = sqrt(sum_rows r_i^2) over all four
    components, divided by normFactor. normFactor is normally supplied by the
    caller once per outer iteration (sum_i |A x_i| + |b_i| + cfVSmall of the full
    system, computed in double by the assembler). If the caller passes a
    negative value, the solver computes the analogue from its own x and b.

    Convergence: finalResidual < tolerance, or
                 finalResidual < relTol * initialResidual,
    after at least minIter iterations.

    All vectors are blockScalar; all dot products and norms go through
    doubleReduce.

SourceFiles
    blockSolver.C

\*---------------------------------------------------------------------------*/

#ifndef coupledFoam_blockSolver_H
#define coupledFoam_blockSolver_H

#include "blockLduMatrix4.H"
#include "diagnostics.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"
#include "autoPtr.H"
#include "word.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class blockPreconditioner;

/*---------------------------------------------------------------------------*\
                   Class blockSolverPerformance Declaration
\*---------------------------------------------------------------------------*/

struct blockSolverPerformance
{
    word solverName;
    reduceScalar initialResidual = 0;
    reduceScalar finalResidual = 0;
    label nIterations = 0;
    label nRestarts = 0;
    bool converged = false;
    bool breakdown = false;

    //- Preconditioner efficiency on the first application of this solve,
    //  rho = ||r0 - A M^-1 r0|| / ||r0|| (amendment B3 6.3.5); -1 if not
    //  measured
    reduceScalar rho = -1;

    //- Write one-line summary
    void print(Ostream& os) const;
};


/*---------------------------------------------------------------------------*\
                         Class blockSolver Declaration
\*---------------------------------------------------------------------------*/

class blockSolver
{
protected:

    // Protected Data

        //- The matrix
        const blockLduMatrix4& matrix_;

        //- Controls
        dictionary controlDict_;

        //- Absolute tolerance on ||r||/normFactor
        reduceScalar tolerance_;

        //- Relative tolerance to the initial residual of this solve
        //  (mutable: set per outer iteration by Eisenstat-Walker, B2)
        mutable reduceScalar relTol_;

        //- Maximum and minimum number of iterations
        label maxIter_;
        label minIter_;

        //- Maximum number of breakdown restarts (spec 6.2)
        label maxRestarts_;

        //- Preconditioner (may be empty: none)
        mutable autoPtr<blockPreconditioner> preconPtr_;

        //- Matrix version the preconditioner was last updated for
        mutable label preconVersion_;

        //- Deep diagnostics (TASK 5; null: off). Hooks are guarded with
        //  diagActive(level) at every call site.
        mutable diagnostics* diag_;


    // Protected Member Functions

        //- True if converged
        bool converged
        (
            const reduceScalar residual,
            const reduceScalar initialResidual,
            const label nIter
        ) const;

        //- Update the preconditioner if the matrix changed
        void updatePreconditioner() const;

        //- w = M^-1 r (identity if no preconditioner)
        void precondition
        (
            blockScalarUList& w,
            const blockScalarUList& r
        ) const;

        //- True if diagnostics at level l are on
        bool diagActive(const label l) const noexcept
        {
            return diag_ && diag_->active(l);
        }

        //- Level-2 hooks shared by the Krylov solvers
        void diagBeginSolve() const;
        void diagEndSolve(const blockSolverPerformance& perf) const;

        //- Scope of one solve for the level-2 log: begins the solve record
        //  on construction and closes it with the final performance on
        //  destruction (every return path). Declare right after perf.
        class diagSolveScope
        {
            const blockSolver& s_;
            const blockSolverPerformance& perf_;

        public:

            diagSolveScope
            (
                const blockSolver& s,
                const blockSolverPerformance& perf
            )
            :
                s_(s),
                perf_(perf)
            {
                s_.diagBeginSolve();
            }

            ~diagSolveScope()
            {
                s_.diagEndSolve(perf_);
            }

            diagSolveScope(const diagSolveScope&) = delete;
            void operator=(const diagSolveScope&) = delete;
        };

        //- Compute normFactor from x, b (sum |A x| + |b| + cfVSmall)
        reduceScalar normFactor
        (
            const blockScalarUList& x,
            const blockScalarUList& b
        ) const;


public:

    //- Runtime type information
    TypeName("blockSolver");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            blockSolver,
            dictionary,
            (
                const blockLduMatrix4& matrix,
                const dictionary& dict
            ),
            (matrix, dict)
        );


    // Constructors

        //- Construct from matrix and controls
        blockSolver(const blockLduMatrix4& matrix, const dictionary& dict);

        //- No copy
        blockSolver(const blockSolver&) = delete;
        void operator=(const blockSolver&) = delete;


    // Selectors

        //- Select from the "solver" keyword of dict
        static autoPtr<blockSolver> New
        (
            const blockLduMatrix4& matrix,
            const dictionary& dict
        );


    //- Destructor
    virtual ~blockSolver();


    // Member Functions

        const blockLduMatrix4& matrix() const noexcept { return matrix_; }
        const dictionary& controlDict() const noexcept { return controlDict_; }

        //- Relative tolerance of the next solves
        reduceScalar relTol() const noexcept { return relTol_; }
        void setRelTol(const reduceScalar r) const noexcept { relTol_ = r; }

        //- Absolute tolerance
        reduceScalar tolerance() const noexcept { return tolerance_; }

        //- rho = ||r - A M^-1 r|| / ||r|| for a given r (one preconditioner
        //  application and one Amul)
        reduceScalar measureRho(const blockScalarUList& r) const;

        //- Preconditioner (may be null)
        const blockPreconditioner* preconditioner() const
        {
            return preconPtr_.get();
        }

        //- Attach deep diagnostics (null: off); also used by the
        //  preconditioner (blockGAMG level-2/3 hooks)
        void setDiagnostics(diagnostics* d) const noexcept { diag_ = d; }
        diagnostics* diag() const noexcept { return diag_; }

        //- Solve A x = b. x holds the initial guess on input.
        //  normFactor < 0: computed from x and b.
        virtual blockSolverPerformance solve
        (
            blockScalarUList& x,
            const blockScalarUList& b,
            const reduceScalar normFactor = -1
        ) const = 0;

        //- Write the effective controls
        virtual void writeSettings(dictionary& dict) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

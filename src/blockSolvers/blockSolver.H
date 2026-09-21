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
    caller once per outer iteration (sum_i |A x_i| + |b_i| + SMALL of the full
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
        reduceScalar relTol_;

        //- Maximum and minimum number of iterations
        label maxIter_;
        label minIter_;

        //- Maximum number of breakdown restarts (spec 6.2)
        label maxRestarts_;

        //- Preconditioner (may be empty: none)
        mutable autoPtr<blockPreconditioner> preconPtr_;

        //- Matrix version the preconditioner was last updated for
        mutable label preconVersion_;


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

        //- Compute normFactor from x, b (sum |A x| + |b| + SMALL)
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

        //- Preconditioner (may be null)
        const blockPreconditioner* preconditioner() const
        {
            return preconPtr_.get();
        }

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

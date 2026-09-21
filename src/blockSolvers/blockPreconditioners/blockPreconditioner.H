/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
-------------------------------------------------------------------------------
Class
    Foam::blockPreconditioner

Description
    Base class of block preconditioners with runtime selection
    (spec 6.1, keyword "preconditioner": blockDiagonal | blockGAMG).

    update() is called by the owning solver whenever the matrix version has
    changed (new coefficients from the assembler); precondition() applies
    w = M^-1 r and must not change state that depends on r.

SourceFiles
    blockPreconditioner.C

\*---------------------------------------------------------------------------*/

#ifndef coupledFoam_blockPreconditioner_H
#define coupledFoam_blockPreconditioner_H

#include "blockSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class blockPreconditioner
{
protected:

    // Protected Data

        //- The owning solver
        const blockSolver& solver_;


public:

    //- Runtime type information
    TypeName("blockPreconditioner");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            blockPreconditioner,
            dictionary,
            (
                const blockSolver& solver,
                const dictionary& dict
            ),
            (solver, dict)
        );


    // Constructors

        explicit blockPreconditioner(const blockSolver& solver)
        :
            solver_(solver)
        {}

        //- No copy
        blockPreconditioner(const blockPreconditioner&) = delete;
        void operator=(const blockPreconditioner&) = delete;


    // Selectors

        //- Select from the "preconditioner" keyword of the solver dict
        static autoPtr<blockPreconditioner> New
        (
            const blockSolver& solver,
            const dictionary& dict
        );


    //- Destructor
    virtual ~blockPreconditioner() = default;


    // Member Functions

        //- The matrix
        const blockLduMatrix4& matrix() const
        {
            return solver_.matrix();
        }

        //- Recompute after a coefficient change
        virtual void update() = 0;

        //- w = M^-1 r
        virtual void precondition
        (
            blockScalarUList& w,
            const blockScalarUList& r
        ) const = 0;

        //- Number of guarded (near-singular) diagonal blocks at last update
        virtual label nSingularDiag() const
        {
            return 0;
        }

        //- Write the effective settings
        virtual void writeSettings(dictionary& dict) const
        {
            dict.add("type", type());
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

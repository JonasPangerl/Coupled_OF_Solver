/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
-------------------------------------------------------------------------------
Class
    Foam::blockDiagonal

Description
    Block-Jacobi preconditioner: w_P = D_P^-1 r_P with the dense 4x4 inverse
    of each diagonal block (block4Ops::invert, pivot guard of spec 6.3).

    Dictionary (solver dict): pivotGuard (default coupledDefaults::pivotGuard)

SourceFiles
    blockDiagonal.C

\*---------------------------------------------------------------------------*/

#ifndef coupledFoam_blockDiagonal_H
#define coupledFoam_blockDiagonal_H

#include "blockPreconditioner.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class blockDiagonal
:
    public blockPreconditioner
{
    // Private Data

        //- Pivot guard
        const reduceScalar pivotGuard_;

        //- Inverse diagonal blocks, 16 per cell
        blockScalarList rD_;

        //- Guarded pivots at last update (this rank)
        label nSingular_;


public:

    //- Runtime type information
    TypeName("blockDiagonal");


    // Constructors

        blockDiagonal(const blockSolver& solver, const dictionary& dict);


    //- Destructor
    virtual ~blockDiagonal() = default;


    // Member Functions

        virtual void update();

        virtual void precondition
        (
            blockScalarUList& w,
            const blockScalarUList& r
        ) const;

        virtual label nSingularDiag() const
        {
            return nSingular_;
        }

        virtual void writeSettings(dictionary& dict) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

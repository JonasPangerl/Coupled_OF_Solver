/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
-------------------------------------------------------------------------------
Class
    Foam::blockGaussSeidel

Description
    Forward block Gauss-Seidel (spec 6.3):
        x_P <- D_P^-1 (b_P - sum_N A_PN x_N)
    with the precomputed dense 4x4 inverse D_P^-1 (pivot guard, count of
    guarded pivots). Loop structure as native GaussSeidelSmoother: one pass
    over cells in order using ownerStartAddr, the lower contributions of
    already-updated cells are pushed into a work copy of b.

    Processor interfaces are treated Jacobi-style: the interface
    contribution is evaluated once per sweep with the x of the previous
    sweep and subtracted from b (as native).

SourceFiles
    blockGaussSeidel.C

\*---------------------------------------------------------------------------*/

#ifndef coupledFoam_blockGaussSeidel_H
#define coupledFoam_blockGaussSeidel_H

#include "blockSmoother.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class blockGaussSeidel
:
    public blockSmoother
{
    // Private Data

        //- Inverse diagonal blocks, 16 per cell
        blockScalarList rD_;

        //- Work copy of b
        mutable blockScalarList bPrime_;


public:

    //- Runtime type information
    TypeName("blockGaussSeidel");


    // Constructors

        blockGaussSeidel(const blockLduMatrix4& matrix, const dictionary& dict);


    //- Destructor
    virtual ~blockGaussSeidel() = default;


    // Member Functions

        virtual void update();

        virtual void smooth
        (
            blockScalarUList& x,
            const blockScalarUList& b,
            const label nSweeps
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
-------------------------------------------------------------------------------
Class
    Foam::blockILU0

Description
    Block incomplete LU smoother in the LDU-native form (spec 6.3 "LDU-native
    incomplete block factorization with 4x4 block pivots"): the factorisation
    M = (D* + L) D*^-1 (D* + U) keeps the off-diagonal blocks of A and
    modifies only the diagonal blocks,
        D*_u = D_u - sum_{f: u(f)=u} L_f D*_l(f)^-1 U_f,
    the block analogue of native DILUPreconditioner. Each sweep is
        x <- x + M^-1 (b - A x).
    Processor interfaces enter only through the residual (block Jacobi
    across ranks).

    Enabled only if Gate B2 shows >= 15 % wall-clock gain on T2 (spec 6.3),
    otherwise available as an option; the decision is in DECISIONS.md.

SourceFiles
    blockILU0.C

\*---------------------------------------------------------------------------*/

#ifndef coupledFoam_blockILU0_H
#define coupledFoam_blockILU0_H

#include "blockSmoother.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class blockILU0
:
    public blockSmoother
{
    // Private Data

        //- Inverse of the modified diagonal blocks D*^-1, 16 per cell
        blockScalarList rD_;

        //- Work vectors
        mutable blockScalarList r_;
        mutable blockScalarList w_;


public:

    //- Runtime type information
    TypeName("blockILU0");


    // Constructors

        blockILU0(const blockLduMatrix4& matrix, const dictionary& dict);


    //- Destructor
    virtual ~blockILU0() = default;


    // Member Functions

        virtual void update();

        //- w = M^-1 r
        void applyInverse(blockScalarUList& w, const blockScalarUList& r) const;

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

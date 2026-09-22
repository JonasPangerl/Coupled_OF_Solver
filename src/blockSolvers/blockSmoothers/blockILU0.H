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

    Pivot-growth safeguard (D-049): on the saddle-point rows the modified
    pivot D*_c can lose almost all of its dominance (T4b: |D*_c^-1| up to
    220 |D_c^-1| on a coarse level); the sweep then amplifies the residual
    by up to 1e5 and the preconditioned Krylov solve stagnates. A cell with
    max|D*_c^-1| > pivotGrowthLimit * max|D_c^-1| uses the unmodified block
    D_c instead (keyword pivotGrowthLimit in the smoother/GAMG dictionary,
    0 = off; default coupledDefaults::iluPivotGrowthLimit = 5). On T4b the
    healthy factorisations at CFL 500 stay at 3.7-3.9, where nothing changes.

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

        //- Pivot-growth limit (0: safeguard off)
        const doubleScalar pivotGrowthLimit_;

        //- Cells with the pivot fallback at the last update
        label nPivotFallback_;

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

        virtual label nPivotFallback() const noexcept
        {
            return nPivotFallback_;
        }

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

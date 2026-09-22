/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
\*---------------------------------------------------------------------------*/

#include "blockGaussSeidel.H"
#include "block4Ops.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(blockGaussSeidel, 0);
    addToRunTimeSelectionTable(blockSmoother, blockGaussSeidel, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::blockGaussSeidel::blockGaussSeidel
(
    const blockLduMatrix4& matrix,
    const dictionary& dict
)
:
    blockSmoother(matrix, dict),
    rD_(),
    bPrime_()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::blockGaussSeidel::update()
{
    const label nCells = matrix_.nCells();

    rD_.resize_nocopy(blockSize*nCells);
    bPrime_.resize_nocopy(blockDim*nCells);
    nSingular_ = 0;

    for (label celli = 0; celli < nCells; ++celli)
    {
        nSingular_ += block4Ops::invert
        (
            matrix_.diagBlock(celli),
            rD_.data() + celli*blockSize,
            pivotGuard_
        );
    }
}


void Foam::blockGaussSeidel::smooth
(
    blockScalarUList& x,
    const blockScalarUList& b,
    const label nSweeps
) const
{
    const label nCells = matrix_.nCells();
    const lduAddressing& addr = matrix_.lduAddr();

    const label* const __restrict__ uAddr = addr.upperAddr().cdata();
    const label* const __restrict__ ownStart = addr.ownerStartAddr().cdata();

    const blockScalar* const __restrict__ upperPtr = matrix_.upper().cdata();
    const blockScalar* const __restrict__ lowerPtr = matrix_.lower().cdata();
    const blockScalar* const __restrict__ rDPtr = rD_.cdata();

    blockScalar* __restrict__ xPtr = x.data();

    for (label sweep = 0; sweep < nSweeps; ++sweep)
    {
        bPrime_ = b;
        blockScalar* __restrict__ bpPtr = bPrime_.data();

        // Jacobi-style interface contribution with x of the previous sweep
        matrix_.initInterfaces(x);
        matrix_.updateInterfaces(bPrime_, true);

        for (label celli = 0; celli < nCells; ++celli)
        {
            const label fStart = ownStart[celli];
            const label fEnd = ownStart[celli + 1];

            blockScalar t[blockDim];
            for (label k = 0; k < blockDim; ++k)
            {
                t[k] = bpPtr[celli*blockDim + k];
            }

            // Upper neighbours (not yet updated in this sweep)
            for (label facei = fStart; facei < fEnd; ++facei)
            {
                block4Ops::matVecSub
                (
                    upperPtr + facei*blockSize,
                    xPtr + uAddr[facei]*blockDim,
                    t
                );
            }

            if (damped_)
            {
                // x_c = x_c + relax*(x_GS - x_c)
                blockScalar xg[blockDim];
                block4Ops::matVec(rDPtr + celli*blockSize, t, xg);
                blockScalar* xc = xPtr + celli*blockDim;
                for (label k = 0; k < blockDim; ++k)
                {
                    xc[k] += relax_*(xg[k] - xc[k]);
                }
            }
            else
            {
                block4Ops::matVec
                (
                    rDPtr + celli*blockSize,
                    t,
                    xPtr + celli*blockDim
                );
            }

            // Push the lower contribution of the updated x into bPrime
            for (label facei = fStart; facei < fEnd; ++facei)
            {
                block4Ops::matVecSub
                (
                    lowerPtr + facei*blockSize,
                    xPtr + celli*blockDim,
                    bpPtr + uAddr[facei]*blockDim
                );
            }
        }
    }
}


// ************************************************************************* //

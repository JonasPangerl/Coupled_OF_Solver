/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
\*---------------------------------------------------------------------------*/

#include "blockILU0.H"
#include "block4Ops.H"
#include "coupledDefaults.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(blockILU0, 0);
    addToRunTimeSelectionTable(blockSmoother, blockILU0, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::blockILU0::blockILU0
(
    const blockLduMatrix4& matrix,
    const dictionary& dict
)
:
    blockSmoother(matrix, dict),
    pivotGrowthLimit_
    (
        dict.getOrDefault<doubleScalar>
        (
            "pivotGrowthLimit",
            coupledDefaults::iluPivotGrowthLimit
        )
    ),
    nPivotFallback_(0),
    rD_(),
    r_(),
    w_()
{
    // Negated comparison also rejects non-finite input
    if (!(pivotGrowthLimit_ >= 0) || !std::isfinite(pivotGrowthLimit_))
    {
        FatalIOErrorInFunction(dict)
            << "pivotGrowthLimit must be a finite value >= 0 (0: off), got "
            << pivotGrowthLimit_ << exit(FatalIOError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::blockILU0::update()
{
    const label nCells = matrix_.nCells();
    const lduAddressing& addr = matrix_.lduAddr();

    const label* const __restrict__ uAddr = addr.upperAddr().cdata();
    const label* const __restrict__ ownStart = addr.ownerStartAddr().cdata();

    const blockScalar* const __restrict__ upperPtr = matrix_.upper().cdata();
    const blockScalar* const __restrict__ lowerPtr = matrix_.lower().cdata();

    // Modified diagonal D*, starts as D
    blockScalarList Dstar(matrix_.diag());

    rD_.resize_nocopy(blockSize*nCells);
    r_.resize_nocopy(blockDim*nCells);
    w_.resize_nocopy(blockDim*nCells);
    nSingular_ = 0;
    nPivotFallback_ = 0;

    blockScalar LDinv[blockSize];
    blockScalar Dinv[blockSize];
    const bool safeguard = (pivotGrowthLimit_ > 0);
    const blockScalar* const __restrict__ diagPtr = matrix_.diag().cdata();

    // Faces are ordered by owner; when cell c is reached all faces with
    // neighbour c (owners < c) have already been processed, so D*_c is final.
    for (label celli = 0; celli < nCells; ++celli)
    {
        nSingular_ += block4Ops::invert
        (
            Dstar.cdata() + celli*blockSize,
            rD_.data() + celli*blockSize,
            pivotGuard_
        );

        // Pivot-growth safeguard (D-049): fall back to the unmodified
        // block if the modified pivot lost its dominance. The fallback
        // pivot enters the modification of the later cells as well, i.e.
        // the factorisation stays consistent.
        if (safeguard)
        {
            const blockScalar* Dc = diagPtr + celli*blockSize;
            block4Ops::invert(Dc, Dinv, pivotGuard_);
            const reduceScalar nStar =
                block4Ops::maxAbs(rD_.cdata() + celli*blockSize);
            const reduceScalar nOrig = block4Ops::maxAbs(Dinv);
            // Negated comparison also catches non-finite values
            if (!(nStar <= pivotGrowthLimit_*nOrig))
            {
                block4Ops::copy(Dinv, rD_.data() + celli*blockSize);
                ++nPivotFallback_;
            }
        }

        const label fStart = ownStart[celli];
        const label fEnd = ownStart[celli + 1];

        for (label facei = fStart; facei < fEnd; ++facei)
        {
            // D*_u -= L_f D*_c^-1 U_f
            block4Ops::matMul
            (
                lowerPtr + facei*blockSize,
                rD_.cdata() + celli*blockSize,
                LDinv
            );
            block4Ops::matMulSub
            (
                LDinv,
                upperPtr + facei*blockSize,
                Dstar.data() + uAddr[facei]*blockSize
            );
        }
    }
}


void Foam::blockILU0::applyInverse
(
    blockScalarUList& w,
    const blockScalarUList& r
) const
{
    const label nCells = matrix_.nCells();
    const lduAddressing& addr = matrix_.lduAddr();

    const label* const __restrict__ uAddr = addr.upperAddr().cdata();
    const label* const __restrict__ ownStart = addr.ownerStartAddr().cdata();

    const blockScalar* const __restrict__ upperPtr = matrix_.upper().cdata();
    const blockScalar* const __restrict__ lowerPtr = matrix_.lower().cdata();
    const blockScalar* const __restrict__ rDPtr = rD_.cdata();

    // Forward: (D* + L) y = r, y_c = D*_c^-1 (r_c - sum L y_l)
    // t accumulates r_c - sum L y_l, pushed from owner to neighbour
    blockScalarList t(r);
    blockScalar* __restrict__ wPtr = w.data();
    blockScalar* __restrict__ tPtr = t.data();

    for (label celli = 0; celli < nCells; ++celli)
    {
        block4Ops::matVec
        (
            rDPtr + celli*blockSize,
            tPtr + celli*blockDim,
            wPtr + celli*blockDim
        );

        const label fStart = ownStart[celli];
        const label fEnd = ownStart[celli + 1];
        for (label facei = fStart; facei < fEnd; ++facei)
        {
            block4Ops::matVecSub
            (
                lowerPtr + facei*blockSize,
                wPtr + celli*blockDim,
                tPtr + uAddr[facei]*blockDim
            );
        }
    }

    // Backward: w_c = y_c - D*_c^-1 sum_{f owned by c} U_f w_u(f)
    blockScalar s[blockDim];
    blockScalar ds[blockDim];

    for (label celli = nCells - 1; celli >= 0; --celli)
    {
        for (label k = 0; k < blockDim; ++k)
        {
            s[k] = 0;
        }

        const label fStart = ownStart[celli];
        const label fEnd = ownStart[celli + 1];
        for (label facei = fStart; facei < fEnd; ++facei)
        {
            block4Ops::matVecAdd
            (
                upperPtr + facei*blockSize,
                wPtr + uAddr[facei]*blockDim,
                s
            );
        }

        block4Ops::matVec(rDPtr + celli*blockSize, s, ds);

        for (label k = 0; k < blockDim; ++k)
        {
            wPtr[celli*blockDim + k] -= ds[k];
        }
    }
}


void Foam::blockILU0::smooth
(
    blockScalarUList& x,
    const blockScalarUList& b,
    const label nSweeps
) const
{
    const label n = matrix_.nRows();

    for (label sweep = 0; sweep < nSweeps; ++sweep)
    {
        matrix_.residual(r_, x, b);
        applyInverse(w_, r_);

        blockScalar* __restrict__ xPtr = x.data();
        const blockScalar* __restrict__ wPtr = w_.cdata();
        if (damped_)
        {
            const blockScalar a = relax_;
            for (label i = 0; i < n; ++i)
            {
                xPtr[i] += a*wPtr[i];
            }
        }
        else
        {
            for (label i = 0; i < n; ++i)
            {
                xPtr[i] += wPtr[i];
            }
        }
    }
}


// ************************************************************************* //

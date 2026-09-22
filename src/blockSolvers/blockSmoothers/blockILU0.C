/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
\*---------------------------------------------------------------------------*/

#include "blockILU0.H"
#include "block4Ops.H"
#include "coupledDefaults.H"
#include "addToRunTimeSelectionTable.H"
#include "doubleReduce.H"
#include "PstreamReduceOps.H"
#include <cstdlib>

namespace
{
    bool iluDebug()
    {
        static const bool on = (std::getenv("CF_ILU_DEBUG") != nullptr);
        return on;
    }
}

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
{}


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
    // DEBUG (A/B on the exact failing system): CF_PIVOT_START=N enables
    // the safeguard only from the N-th update of this smoother on
    static const char* ps = std::getenv("CF_PIVOT_START");
    static const label pStart = (ps ? label(std::atol(ps)) : 0);
    static label maxCellsDbg = -1;
    static label finestUpd = 0;
    if (nCells >= maxCellsDbg)
    {
        maxCellsDbg = nCells;
        ++finestUpd;
    }
    const bool safeguard = (pivotGrowthLimit_ > 0) && finestUpd >= pStart;
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

        // Pivot-growth safeguard: the modified pivot D*_c can lose its
        // dominance on the saddle-point rows (coarse levels in particular);
        // its inverse then amplifies the residual in every sweep. Fall back
        // to the unmodified block D_c for such cells (their row becomes a
        // block Jacobi/Gauss-Seidel row inside the ILU sweep).
        if (safeguard)
        {
            const blockScalar* Dc = diagPtr + celli*blockSize;
            nSingular_ += block4Ops::invert(Dc, Dinv, pivotGuard_);
            const reduceScalar nStar =
                block4Ops::maxAbs(rD_.cdata() + celli*blockSize);
            const reduceScalar nOrig = block4Ops::maxAbs(Dinv);
            // Negated comparison also catches non-finite values
            if (!(nStar <= pivotGrowthLimit_*nOrig))
            {
                block4Ops::copy(Dinv, rD_.data() + celli*blockSize);
                block4Ops::copy(Dc, Dstar.data() + celli*blockSize);
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

    if (iluDebug())
    {
        const blockScalarList& D = matrix_.diag();
        grow_.resize_nocopy(nCells);
        blockScalar Dinv[blockSize];
        double growMax = 0, condStarMax = 0, condOrigMax = 0;
        label worst = -1, nGrow10 = 0, nGrow1e3 = 0, nCond1e6 = 0;
        for (label celli = 0; celli < nCells; ++celli)
        {
            const blockScalar* Dc = D.cdata() + celli*blockSize;
            block4Ops::invert(Dc, Dinv, pivotGuard_);
            const double nD = block4Ops::maxAbs(Dc);
            const double nDi = block4Ops::maxAbs(Dinv);
            const double nS = block4Ops::maxAbs(Dstar.cdata() + celli*blockSize);
            const double nRD = block4Ops::maxAbs(rD_.cdata() + celli*blockSize);
            const double g = nRD/std::max(nDi, 1e-300);
            const double cs = nS*nRD;
            const double co = nD*nDi;
            grow_[celli] = g;
            if (g > growMax) { growMax = g; worst = celli; }
            if (g > 10) ++nGrow10;
            if (g > 1e3) ++nGrow1e3;
            if (cs > 1e6) ++nCond1e6;
            condStarMax = std::max(condStarMax, cs);
            condOrigMax = std::max(condOrigMax, co);
        }
        const label comm = matrix_.comm();
        const double growLocal = growMax;
        reduce(growMax, maxOp<double>(), UPstream::msgType(), comm);
        reduce(condStarMax, maxOp<double>(), UPstream::msgType(), comm);
        reduce(condOrigMax, maxOp<double>(), UPstream::msgType(), comm);
        reduce(nGrow10, sumOp<label>(), UPstream::msgType(), comm);
        reduce(nGrow1e3, sumOp<label>(), UPstream::msgType(), comm);
        reduce(nCond1e6, sumOp<label>(), UPstream::msgType(), comm);
        label nTot = nCells;
        reduce(nTot, sumOp<label>(), UPstream::msgType(), comm);
        if (UPstream::myProcNo(comm) == 0)
        {
            Pout<< "ILUDBG update nCells " << nTot
                << " growMax " << growMax << " nGrow>10 " << nGrow10
                << " nGrow>1e3 " << nGrow1e3
                << " condStarMax " << condStarMax
                << " nCondStar>1e6 " << nCond1e6
                << " condOrigMax " << condOrigMax << endl;
        }
        if (growLocal == growMax && growMax > 10 && worst >= 0)
        {
            const blockScalar* Dc = D.cdata() + worst*blockSize;
            const blockScalar* Sc = Dstar.cdata() + worst*blockSize;
            Pout<< "ILUDBG worst cell " << worst << " of " << nCells
                << " grow " << growMax << nl << "  D  =";
            for (label i = 0; i < blockSize; ++i) Pout<< ' ' << Dc[i];
            Pout<< nl << "  D* =";
            for (label i = 0; i < blockSize; ++i) Pout<< ' ' << Sc[i];
            Pout<< endl;
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
    const bool dbg = iluDebug();
    double r0 = -1;

    for (label sweep = 0; sweep < nSweeps; ++sweep)
    {
        matrix_.residual(r_, x, b);
        if (dbg && sweep == 0)
        {
            r0 = doubleReduce::norm2(r_, matrix_.comm());
        }
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

    if (dbg && nSweeps > 0)
    {
        const label comm = matrix_.comm();
        matrix_.residual(r_, x, b);
        const double r1 = doubleReduce::norm2(r_, comm);
        // cell with the largest residual after smoothing
        double rmax = 0;
        label cmax = -1;
        const label nCells = matrix_.nCells();
        for (label c = 0; c < nCells; ++c)
        {
            double q = 0;
            for (label k = 0; k < blockDim; ++k)
            {
                q += double(r_[c*blockDim + k])*double(r_[c*blockDim + k]);
            }
            if (q > rmax) { rmax = q; cmax = c; }
        }
        rmax = std::sqrt(rmax);
        double rmaxG = rmax;
        reduce(rmaxG, maxOp<double>(), UPstream::msgType(), comm);
        const double ratio = r1/std::max(r0, 1e-300);
        if (ratio > 2 && rmax == rmaxG && cmax >= 0)
        {
            Pout<< "ILUDBG smooth nCells " << nCells << " sweeps " << nSweeps
                << " |r| " << r0 << " -> " << r1 << " (x" << ratio
                << ") max cell " << cmax << " |r_c| " << rmax
                << " share " << rmax/std::max(r1, 1e-300)
                << " grow_c " << (grow_.size() > cmax ? grow_[cmax] : -1.0)
                << endl;
            const blockScalar* Dc = matrix_.diag().cdata() + cmax*blockSize;
            Pout<< "ILUDBG   D =";
            for (label i = 0; i < blockSize; ++i) Pout<< ' ' << Dc[i];
            Pout<< nl << "ILUDBG   r_c =";
            for (label k = 0; k < blockDim; ++k) Pout<< ' ' << r_[cmax*blockDim + k];
            Pout<< nl << "ILUDBG   x_c =";
            for (label k = 0; k < blockDim; ++k) Pout<< ' ' << x[cmax*blockDim + k];
            Pout<< nl << "ILUDBG   b_c =";
            for (label k = 0; k < blockDim; ++k) Pout<< ' ' << b[cmax*blockDim + k];
            Pout<< nl;
            const lduAddressing& addr = matrix_.lduAddr();
            const labelUList& lA = addr.lowerAddr();
            const labelUList& uA = addr.upperAddr();
            forAll(lA, facei)
            {
                label nb = -1;
                const blockScalar* rowBlk = nullptr;
                if (lA[facei] == cmax)
                {
                    nb = uA[facei];
                    rowBlk = matrix_.upper().cdata() + facei*blockSize;
                }
                else if (uA[facei] == cmax)
                {
                    nb = lA[facei];
                    rowBlk = matrix_.lower().cdata() + facei*blockSize;
                }
                if (nb >= 0)
                {
                    Pout<< "ILUDBG   nb " << nb << " A_c,nb =";
                    for (label i = 0; i < blockSize; ++i) Pout<< ' ' << rowBlk[i];
                    Pout<< " x_nb =";
                    for (label k = 0; k < blockDim; ++k) Pout<< ' ' << x[nb*blockDim + k];
                    Pout<< nl;
                }
            }
            Pout<< "ILUDBG   interfaces " << matrix_.interfaces().size() << endl;
        }
    }
}


// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
-------------------------------------------------------------------------------
Namespace
    Foam::block4Ops

Description
    Inline dense 4x4 operations on row-major blocks (index r*4 + c) stored as
    blockScalar. Arithmetic that feeds a stored result accumulates in double
    where the extra cost is negligible (inversion, block products); the
    matrix-vector kernels used inside Krylov/smoother loops run in float,
    except matVecSubDouble, the kernel of blockLduMatrix4::residualDouble
    (the true residual b - A x of a Krylov solver, D-012).

    invert() uses LU with partial pivoting in double. Pivot guard (spec 6.3):
    a pivot with |pivot| < pivotGuard gets pivotGuard added (with the sign of
    the pivot, or + for an exact zero) and the call reports it, so the caller
    can count nSingularDiag. No NaN or Inf is produced for a singular block.

\*---------------------------------------------------------------------------*/

#ifndef coupledFoam_block4Ops_H
#define coupledFoam_block4Ops_H

#include "blockScalar.H"
#include <cmath>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace block4Ops
{

//- y = A x
inline void matVec
(
    const blockScalar* __restrict__ A,
    const blockScalar* __restrict__ x,
    blockScalar* __restrict__ y
)
{
    for (label r = 0; r < blockDim; ++r)
    {
        const blockScalar* a = A + r*blockDim;
        y[r] = a[0]*x[0] + a[1]*x[1] + a[2]*x[2] + a[3]*x[3];
    }
}


//- y += A x
inline void matVecAdd
(
    const blockScalar* __restrict__ A,
    const blockScalar* __restrict__ x,
    blockScalar* __restrict__ y
)
{
    for (label r = 0; r < blockDim; ++r)
    {
        const blockScalar* a = A + r*blockDim;
        y[r] += a[0]*x[0] + a[1]*x[1] + a[2]*x[2] + a[3]*x[3];
    }
}


//- y -= A x
inline void matVecSub
(
    const blockScalar* __restrict__ A,
    const blockScalar* __restrict__ x,
    blockScalar* __restrict__ y
)
{
    for (label r = 0; r < blockDim; ++r)
    {
        const blockScalar* a = A + r*blockDim;
        y[r] -= a[0]*x[0] + a[1]*x[1] + a[2]*x[2] + a[3]*x[3];
    }
}


//- y -= A x with y in double. The products of two blockScalar values are
//  exact in double and the row sum is accumulated in double, so no
//  cancellation error of the order eps_float*|A||x| enters y. Kernel of
//  blockLduMatrix4::residualDouble (D-012), whose accuracy decides Krylov
//  convergence.
inline void matVecSubDouble
(
    const blockScalar* __restrict__ A,
    const blockScalar* __restrict__ x,
    reduceScalar* __restrict__ y
)
{
    const reduceScalar x0 = toDouble(x[0]);
    const reduceScalar x1 = toDouble(x[1]);
    const reduceScalar x2 = toDouble(x[2]);
    const reduceScalar x3 = toDouble(x[3]);
    for (label r = 0; r < blockDim; ++r)
    {
        const blockScalar* a = A + r*blockDim;
        y[r] -=
            toDouble(a[0])*x0 + toDouble(a[1])*x1
          + toDouble(a[2])*x2 + toDouble(a[3])*x3;
    }
}


//- C = A B (accumulated in double)
inline void matMul
(
    const blockScalar* __restrict__ A,
    const blockScalar* __restrict__ B,
    blockScalar* __restrict__ C
)
{
    for (label r = 0; r < blockDim; ++r)
    {
        for (label c = 0; c < blockDim; ++c)
        {
            reduceScalar s = 0;
            for (label k = 0; k < blockDim; ++k)
            {
                s += toDouble(A[r*blockDim + k])*toDouble(B[k*blockDim + c]);
            }
            C[r*blockDim + c] = narrow(s);
        }
    }
}


//- C -= A B (accumulated in double)
inline void matMulSub
(
    const blockScalar* __restrict__ A,
    const blockScalar* __restrict__ B,
    blockScalar* __restrict__ C
)
{
    for (label r = 0; r < blockDim; ++r)
    {
        for (label c = 0; c < blockDim; ++c)
        {
            reduceScalar s = toDouble(C[r*blockDim + c]);
            for (label k = 0; k < blockDim; ++k)
            {
                s -= toDouble(A[r*blockDim + k])*toDouble(B[k*blockDim + c]);
            }
            C[r*blockDim + c] = narrow(s);
        }
    }
}


//- Copy a block
inline void copy
(
    const blockScalar* __restrict__ A,
    blockScalar* __restrict__ B
)
{
    for (label i = 0; i < blockSize; ++i)
    {
        B[i] = A[i];
    }
}


//- Set a block to zero
inline void zero(blockScalar* A)
{
    for (label i = 0; i < blockSize; ++i)
    {
        A[i] = 0;
    }
}


//- Set a block to the identity
inline void identity(blockScalar* A)
{
    zero(A);
    for (label i = 0; i < blockDim; ++i)
    {
        A[i*blockDim + i] = 1;
    }
}


//- Max-abs norm of a block
inline reduceScalar maxAbs(const blockScalar* A)
{
    reduceScalar m = 0;
    for (label i = 0; i < blockSize; ++i)
    {
        const reduceScalar a = std::abs(toDouble(A[i]));
        if (a > m) m = a;
    }
    return m;
}


//- Invert a 4x4 block by LU with partial pivoting, computed in double.
//  Returns the number of guarded pivots (0 for a regular block).
inline label invert
(
    const blockScalar* __restrict__ A,
    blockScalar* __restrict__ Ainv,
    const reduceScalar pivotGuard
)
{
    reduceScalar lu[blockDim][blockDim];
    label perm[blockDim];

    for (label r = 0; r < blockDim; ++r)
    {
        perm[r] = r;
        for (label c = 0; c < blockDim; ++c)
        {
            lu[r][c] = toDouble(A[r*blockDim + c]);
        }
    }

    label nGuarded = 0;

    for (label k = 0; k < blockDim; ++k)
    {
        // Partial pivoting: largest magnitude in column k at or below k
        label piv = k;
        reduceScalar big = std::abs(lu[k][k]);
        for (label r = k + 1; r < blockDim; ++r)
        {
            const reduceScalar a = std::abs(lu[r][k]);
            if (a > big)
            {
                big = a;
                piv = r;
            }
        }
        if (piv != k)
        {
            for (label c = 0; c < blockDim; ++c)
            {
                const reduceScalar t = lu[k][c];
                lu[k][c] = lu[piv][c];
                lu[piv][c] = t;
            }
            const label t = perm[k];
            perm[k] = perm[piv];
            perm[piv] = t;
        }

        // GUARD: pivot guard (spec 6.3)
        if (std::abs(lu[k][k]) < pivotGuard)
        {
            lu[k][k] += (lu[k][k] < 0 ? -pivotGuard : pivotGuard);
            ++nGuarded;
        }

        const reduceScalar rPiv = 1.0/lu[k][k];  // GUARD: pivot guarded above

        for (label r = k + 1; r < blockDim; ++r)
        {
            lu[r][k] *= rPiv;
            const reduceScalar f = lu[r][k];
            for (label c = k + 1; c < blockDim; ++c)
            {
                lu[r][c] -= f*lu[k][c];
            }
        }
    }

    // Solve LU X = P I column by column
    for (label col = 0; col < blockDim; ++col)
    {
        reduceScalar y[blockDim];

        // Forward substitution with the permuted unit vector
        for (label r = 0; r < blockDim; ++r)
        {
            reduceScalar s = (perm[r] == col ? 1.0 : 0.0);
            for (label c = 0; c < r; ++c)
            {
                s -= lu[r][c]*y[c];
            }
            y[r] = s;
        }

        // Back substitution
        for (label r = blockDim - 1; r >= 0; --r)
        {
            reduceScalar s = y[r];
            for (label c = r + 1; c < blockDim; ++c)
            {
                s -= lu[r][c]*y[c];
            }
            y[r] = s/lu[r][r];  // GUARD: pivot guarded in factorisation
        }

        for (label r = 0; r < blockDim; ++r)
        {
            Ainv[r*blockDim + col] = narrow(y[r]);
        }
    }

    return nGuarded;
}

} // End namespace block4Ops
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

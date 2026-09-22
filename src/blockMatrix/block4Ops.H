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
#include <cstring>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace block4Ops
{

// Row sums of a 4x4 block times a 4-vector (spec 6.5.3).
//
// Every matrix-vector kernel below evaluates row r as
//     ((a_r0 x_0 + a_r1 x_1) + a_r2 x_2) + a_r3 x_3
// (the order of the original scalar code, so results are bit-identical to
// it). With GCC/Clang the 16 products are formed as four 4-wide vector
// multiplies of the row-major blocks (the AoS layout stays, 6.5.4), the
// 4x4 product block is transposed in registers (8 shuffles) and the four
// columns are added in the order above: fixed size, no runtime bounds.
// Measured in Amul on the Test-kernelBandwidth system (Xeon E5-2670, AVX,
// one thread): scalar form 7.6 GB/s; fixed-trip loops with
// `#pragma omp simd` 7.3-8.2 GB/s (GCC 13 does not vectorise the row
// reduction of a row-major block); horizontal adds 10.2-10.5 GB/s (other
// summation order); this transposed form 10.3 GB/s and bit-identical.

#if defined(__GNUC__) || defined(__clang__)
#define COUPLEDFOAM_BLOCK4_VECTOR 1

namespace detail
{

//- blockDim blockScalar values in one vector register
typedef blockScalar vec4
    __attribute__((vector_size(blockDim*sizeof(blockScalar))));

//- Shuffle index vector
typedef int vec4i __attribute__((vector_size(blockDim*sizeof(int))));

inline vec4 load(const blockScalar* p)
{
    vec4 v;
    std::memcpy(&v, p, sizeof(v));
    return v;
}

inline void store(blockScalar* p, const vec4 v)
{
    std::memcpy(p, &v, sizeof(v));
}

//- Two-source shuffle with constant indices (0-3: a, 4-7: b)
#if defined(__clang__)
#define COUPLEDFOAM_SHUFFLE(a, b, i0, i1, i2, i3)                              \
    __builtin_shufflevector(a, b, i0, i1, i2, i3)
#else
#define COUPLEDFOAM_SHUFFLE(a, b, i0, i1, i2, i3)                              \
    __builtin_shuffle(a, b, vec4i{i0, i1, i2, i3})
#endif

//- Row sums of A x, row r = ((p_r0 + p_r1) + p_r2) + p_r3
inline vec4 rowSums
(
    const blockScalar* __restrict__ A,
    const blockScalar* __restrict__ x
)
{
    const vec4 xv = load(x);

    // Products, one row per vector
    const vec4 p0 = load(A)*xv;
    const vec4 p1 = load(A + blockDim)*xv;
    const vec4 p2 = load(A + 2*blockDim)*xv;
    const vec4 p3 = load(A + 3*blockDim)*xv;

    // 4x4 transpose: ck[r] = p_rk
    const vec4 t0 = COUPLEDFOAM_SHUFFLE(p0, p1, 0, 4, 1, 5);
    const vec4 t1 = COUPLEDFOAM_SHUFFLE(p2, p3, 0, 4, 1, 5);
    const vec4 t2 = COUPLEDFOAM_SHUFFLE(p0, p1, 2, 6, 3, 7);
    const vec4 t3 = COUPLEDFOAM_SHUFFLE(p2, p3, 2, 6, 3, 7);
    const vec4 c0 = COUPLEDFOAM_SHUFFLE(t0, t1, 0, 1, 4, 5);
    const vec4 c1 = COUPLEDFOAM_SHUFFLE(t0, t1, 2, 3, 6, 7);
    const vec4 c2 = COUPLEDFOAM_SHUFFLE(t2, t3, 0, 1, 4, 5);
    const vec4 c3 = COUPLEDFOAM_SHUFFLE(t2, t3, 2, 3, 6, 7);

    return ((c0 + c1) + c2) + c3;
}

#undef COUPLEDFOAM_SHUFFLE

} // End namespace detail
#endif


//- y = A x
inline void matVec
(
    const blockScalar* __restrict__ A,
    const blockScalar* __restrict__ x,
    blockScalar* __restrict__ y
)
{
    #ifdef COUPLEDFOAM_BLOCK4_VECTOR
    detail::store(y, detail::rowSums(A, x));
    #else
    for (label r = 0; r < blockDim; ++r)
    {
        const blockScalar* a = A + r*blockDim;
        y[r] = a[0]*x[0] + a[1]*x[1] + a[2]*x[2] + a[3]*x[3];
    }
    #endif
}


//- y += A x
inline void matVecAdd
(
    const blockScalar* __restrict__ A,
    const blockScalar* __restrict__ x,
    blockScalar* __restrict__ y
)
{
    #ifdef COUPLEDFOAM_BLOCK4_VECTOR
    detail::store(y, detail::load(y) + detail::rowSums(A, x));
    #else
    for (label r = 0; r < blockDim; ++r)
    {
        const blockScalar* a = A + r*blockDim;
        y[r] += a[0]*x[0] + a[1]*x[1] + a[2]*x[2] + a[3]*x[3];
    }
    #endif
}


//- y -= A x
inline void matVecSub
(
    const blockScalar* __restrict__ A,
    const blockScalar* __restrict__ x,
    blockScalar* __restrict__ y
)
{
    #ifdef COUPLEDFOAM_BLOCK4_VECTOR
    detail::store(y, detail::load(y) - detail::rowSums(A, x));
    #else
    for (label r = 0; r < blockDim; ++r)
    {
        const blockScalar* a = A + r*blockDim;
        y[r] -= a[0]*x[0] + a[1]*x[1] + a[2]*x[2] + a[3]*x[3];
    }
    #endif
}


//- y -= A x with y in double. The products of two blockScalar values are
//  exact in double and the row sum is accumulated in double (same order
//  as above), so no cancellation error of the order eps_float*|A||x|
//  enters y. Kernel of blockLduMatrix4::residualDouble (D-012), whose
//  accuracy decides Krylov convergence. Fixed-trip loops over the 16
//  products and the 4 rows.
inline void matVecSubDouble
(
    const blockScalar* __restrict__ A,
    const blockScalar* __restrict__ x,
    reduceScalar* __restrict__ y
)
{
    reduceScalar p[blockSize];
    #pragma omp simd
    for (label i = 0; i < blockSize; ++i)
    {
        p[i] = toDouble(A[i])*toDouble(x[i % blockDim]);
    }
    for (label r = 0; r < blockDim; ++r)
    {
        const reduceScalar* q = p + r*blockDim;
        y[r] -= q[0] + q[1] + q[2] + q[3];
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

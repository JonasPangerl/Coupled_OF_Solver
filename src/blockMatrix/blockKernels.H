/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
-------------------------------------------------------------------------------
Namespace
    Foam::blockKernels

Description
    Hand-written vector kernels of the Krylov solvers (spec 6.5.1-6.5.3):
    flat blockScalar arrays, __restrict__ pointers, one pass over the
    operands, `#pragma omp simd` (compiled with -fopenmp-simd, no OpenMP
    runtime). All kernels are local to this rank; the caller reduces the
    returned double partial sums over the communicator (doubleReduce::parSum),
    so several sums share one MPI call.

    Precision (6.0, D-012): vectors are blockScalar, every accumulation is
    reduceScalar (double). The element-wise updates are written exactly as
    the unfused loops they replace, so the iterates are bit-identical; only
    the order in which the double partial sums are added differs (SIMD
    lanes), i.e. the reductions differ in the last bits of a double.

    Fused kernels (6.5.2), bytes moved per entry in brackets:
        axpy_dot              y += a x; <y,y>, <y,z>        [read x,y,z; write y]
        update_residual_norm  r -= a q; ||r||^2              [read r,q; write r]
        mgs_axpy_dot          w -= h v; <w,vNext>            [read w,v,vNext; write w]
            (fused modified Gram-Schmidt step: the axpy with basis vector i
            and the dot product with basis vector i+1 in one pass; the dot
            uses the updated w, so the result is the MGS coefficient h_(i+1))
        dot_sumSqr            <a,b>, <a,a>                   [read a,b]
        axpyMulti(Double)     y += sum_i c_i x_i, chunked so y stays in L1
            and every x_i is read once; per entry the terms are added in the
            order i = 0, 1, ... as in the unfused loop sequence

\*---------------------------------------------------------------------------*/

#ifndef coupledFoam_blockKernels_H
#define coupledFoam_blockKernels_H

#include "blockScalar.H"
#include "coupledDefaults.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace blockKernels
{

//- Two double partial sums of one pass
struct sumPair
{
    reduceScalar first;
    reduceScalar second;
};


//- Local <a,b> in double
inline reduceScalar dot
(
    const label n,
    const blockScalar* __restrict__ a,
    const blockScalar* __restrict__ b
)
{
    reduceScalar s = 0;
    #pragma omp simd reduction(+:s)
    for (label i = 0; i < n; ++i)
    {
        s += toDouble(a[i])*toDouble(b[i]);
    }
    return s;
}


//- Local <a,a> in double
inline reduceScalar sumSqr
(
    const label n,
    const blockScalar* __restrict__ a
)
{
    reduceScalar s = 0;
    #pragma omp simd reduction(+:s)
    for (label i = 0; i < n; ++i)
    {
        const reduceScalar v = toDouble(a[i]);
        s += v*v;
    }
    return s;
}


//- Local (<a,b>, <a,a>) in one pass
inline sumPair dot_sumSqr
(
    const label n,
    const blockScalar* __restrict__ a,
    const blockScalar* __restrict__ b
)
{
    reduceScalar ab = 0;
    reduceScalar aa = 0;
    #pragma omp simd reduction(+:ab, aa)
    for (label i = 0; i < n; ++i)
    {
        const reduceScalar ai = toDouble(a[i]);
        ab += ai*toDouble(b[i]);
        aa += ai*ai;
    }
    return {ab, aa};
}


//- y += a x; returns local (<y,y>, <y,z>) of the updated y
inline sumPair axpy_dot
(
    const label n,
    const blockScalar a,
    const blockScalar* __restrict__ x,
    blockScalar* __restrict__ y,
    const blockScalar* __restrict__ z
)
{
    reduceScalar yy = 0;
    reduceScalar yz = 0;
    #pragma omp simd reduction(+:yy, yz)
    for (label i = 0; i < n; ++i)
    {
        const blockScalar yi = y[i] + a*x[i];
        y[i] = yi;
        const reduceScalar yd = toDouble(yi);
        yy += yd*yd;
        yz += yd*toDouble(z[i]);
    }
    return {yy, yz};
}


//- r -= a q; returns local ||r||^2 of the updated r
inline reduceScalar update_residual_norm
(
    const label n,
    const blockScalar a,
    const blockScalar* __restrict__ q,
    blockScalar* __restrict__ r
)
{
    reduceScalar rr = 0;
    #pragma omp simd reduction(+:rr)
    for (label i = 0; i < n; ++i)
    {
        const blockScalar ri = r[i] - a*q[i];
        r[i] = ri;
        const reduceScalar rd = toDouble(ri);
        rr += rd*rd;
    }
    return rr;
}


//- Fused MGS step: w -= h v; returns local <w,vNext> of the updated w
inline reduceScalar mgs_axpy_dot
(
    const label n,
    const blockScalar h,
    const blockScalar* __restrict__ v,
    blockScalar* __restrict__ w,
    const blockScalar* __restrict__ vNext
)
{
    reduceScalar s = 0;
    #pragma omp simd reduction(+:s)
    for (label i = 0; i < n; ++i)
    {
        const blockScalar wi = w[i] - h*v[i];
        w[i] = wi;
        s += toDouble(wi)*toDouble(vNext[i]);
    }
    return s;
}


//- y = s x
inline void scale
(
    const label n,
    const blockScalar s,
    const blockScalar* __restrict__ x,
    blockScalar* __restrict__ y
)
{
    #pragma omp simd
    for (label i = 0; i < n; ++i)
    {
        y[i] = s*x[i];
    }
}


//- y = x
inline void copy
(
    const label n,
    const blockScalar* __restrict__ x,
    blockScalar* __restrict__ y
)
{
    #pragma omp simd
    for (label i = 0; i < n; ++i)
    {
        y[i] = x[i];
    }
}


//- y = x (widening)
inline void widen
(
    const label n,
    const blockScalar* __restrict__ x,
    reduceScalar* __restrict__ y
)
{
    #pragma omp simd
    for (label i = 0; i < n; ++i)
    {
        y[i] = toDouble(x[i]);
    }
}


//- y = narrow(x)
inline void narrowCopy
(
    const label n,
    const reduceScalar* __restrict__ x,
    blockScalar* __restrict__ y
)
{
    #pragma omp simd
    for (label i = 0; i < n; ++i)
    {
        y[i] = narrow(x[i]);
    }
}


//- y += sum_(i<m) c_i x_i in blockScalar, chunked: y is read and written
//  once, every x_i read once. Per entry the terms are added in the order
//  i = 0 .. m-1 (bit-identical to m separate axpy passes).
inline void axpyMulti
(
    const label n,
    const label m,
    const blockScalar* __restrict__ c,
    const blockScalar* const* __restrict__ x,
    blockScalar* __restrict__ y
)
{
    constexpr label chunk = coupledDefaults::kernelChunk;
    for (label start = 0; start < n; start += chunk)
    {
        const label end = (start + chunk < n ? start + chunk : n);
        for (label i = 0; i < m; ++i)
        {
            const blockScalar ci = c[i];
            const blockScalar* __restrict__ xi = x[i];
            #pragma omp simd
            for (label k = start; k < end; ++k)
            {
                y[k] += ci*xi[k];
            }
        }
    }
}


//- y += sum_(i<m) c_i x_i with a double y and double coefficients (the
//  products are formed in double), chunked as axpyMulti
inline void axpyMultiDouble
(
    const label n,
    const label m,
    const reduceScalar* __restrict__ c,
    const blockScalar* const* __restrict__ x,
    reduceScalar* __restrict__ y
)
{
    constexpr label chunk = coupledDefaults::kernelChunk;
    for (label start = 0; start < n; start += chunk)
    {
        const label end = (start + chunk < n ? start + chunk : n);
        for (label i = 0; i < m; ++i)
        {
            const reduceScalar ci = c[i];
            const blockScalar* __restrict__ xi = x[i];
            #pragma omp simd
            for (label k = start; k < end; ++k)
            {
                y[k] += ci*toDouble(xi[k]);
            }
        }
    }
}


//- STREAM triad a = b + s c (Test-kernelBandwidth reference)
inline void triad
(
    const label n,
    const blockScalar s,
    const blockScalar* __restrict__ b,
    const blockScalar* __restrict__ c,
    blockScalar* __restrict__ a
)
{
    #pragma omp simd
    for (label i = 0; i < n; ++i)
    {
        a[i] = b[i] + s*c[i];
    }
}

} // End namespace blockKernels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

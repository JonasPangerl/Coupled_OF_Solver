/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
\*---------------------------------------------------------------------------*/

#include "doubleReduce.H"
#include "blockKernels.H"
#include "PstreamReduceOps.H"
#include "ops.H"
#include <cmath>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::reduceScalar Foam::doubleReduce::localSum(const blockScalarUList& a)
{
    reduceScalar s = 0;
    const blockScalar* __restrict__ ap = a.cdata();
    const label n = a.size();
    #pragma omp simd reduction(+:s)
    for (label i = 0; i < n; ++i)
    {
        s += toDouble(ap[i]);
    }
    return s;
}


Foam::reduceScalar Foam::doubleReduce::localDot
(
    const blockScalarUList& a,
    const blockScalarUList& b
)
{
    return blockKernels::dot(a.size(), a.cdata(), b.cdata());
}


Foam::reduceScalar Foam::doubleReduce::localSumSqr(const blockScalarUList& a)
{
    return blockKernels::sumSqr(a.size(), a.cdata());
}


Foam::reduceScalar Foam::doubleReduce::localSumMag(const blockScalarUList& a)
{
    reduceScalar s = 0;
    const blockScalar* __restrict__ ap = a.cdata();
    const label n = a.size();
    #pragma omp simd reduction(+:s)
    for (label i = 0; i < n; ++i)
    {
        s += std::abs(toDouble(ap[i]));
    }
    return s;
}


Foam::reduceScalar Foam::doubleReduce::parSum
(
    reduceScalar v,
    const label comm
)
{
    Foam::reduce(v, sumOp<reduceScalar>(), UPstream::msgType(), comm);
    return v;
}


void Foam::doubleReduce::parSum
(
    reduceScalar* v,
    const label n,
    const label comm
)
{
    Foam::reduce(v, int(n), sumOp<reduceScalar>(), UPstream::msgType(), comm);
}


Foam::reduceScalar Foam::doubleReduce::sum
(
    const blockScalarUList& a,
    const label comm
)
{
    return parSum(localSum(a), comm);
}


Foam::reduceScalar Foam::doubleReduce::dot
(
    const blockScalarUList& a,
    const blockScalarUList& b,
    const label comm
)
{
    return parSum(localDot(a, b), comm);
}


Foam::reduceScalar Foam::doubleReduce::sumSqr
(
    const blockScalarUList& a,
    const label comm
)
{
    return parSum(localSumSqr(a), comm);
}


Foam::reduceScalar Foam::doubleReduce::norm2
(
    const blockScalarUList& a,
    const label comm
)
{
    // GUARD: sum of squares is >= 0 by construction; clamp anyway for sqrt
    const reduceScalar s = sumSqr(a, comm);
    return std::sqrt(s > 0 ? s : 0);
}


Foam::reduceScalar Foam::doubleReduce::sumMag
(
    const blockScalarUList& a,
    const label comm
)
{
    return parSum(localSumMag(a), comm);
}


Foam::FixedList<Foam::reduceScalar, 2> Foam::doubleReduce::dot2
(
    const blockScalarUList& a,
    const blockScalarUList& b,
    const blockScalarUList& c,
    const blockScalarUList& d,
    const label comm
)
{
    FixedList<reduceScalar, 2> v;
    v[0] = localDot(a, b);
    v[1] = localDot(c, d);
    parSum(v.data(), 2, comm);
    return v;
}


Foam::FixedList<Foam::reduceScalar, 4> Foam::doubleReduce::componentSumSqr
(
    const blockScalarUList& a,
    const label comm
)
{
    FixedList<reduceScalar, 4> v(Zero);
    const blockScalar* __restrict__ ap = a.cdata();
    const label nCells = a.size()/blockDim;
    reduceScalar s0 = 0, s1 = 0, s2 = 0, s3 = 0;
    #pragma omp simd reduction(+:s0, s1, s2, s3)
    for (label i = 0; i < nCells; ++i)
    {
        const blockScalar* __restrict__ c = ap + i*blockDim;
        const reduceScalar x0 = toDouble(c[0]);
        const reduceScalar x1 = toDouble(c[1]);
        const reduceScalar x2 = toDouble(c[2]);
        const reduceScalar x3 = toDouble(c[3]);
        s0 += x0*x0;
        s1 += x1*x1;
        s2 += x2*x2;
        s3 += x3*x3;
    }
    v[0] = s0;
    v[1] = s1;
    v[2] = s2;
    v[3] = s3;
    parSum(v.data(), blockDim, comm);
    return v;
}


// ************************************************************************* //

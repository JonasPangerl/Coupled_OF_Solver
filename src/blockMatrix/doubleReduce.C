/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
\*---------------------------------------------------------------------------*/

#include "doubleReduce.H"
#include "PstreamReduceOps.H"
#include "ops.H"
#include <cmath>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::reduceScalar Foam::doubleReduce::localSum(const blockScalarUList& a)
{
    reduceScalar s = 0;
    const blockScalar* __restrict__ ap = a.cdata();
    const label n = a.size();
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
    reduceScalar s = 0;
    const blockScalar* __restrict__ ap = a.cdata();
    const blockScalar* __restrict__ bp = b.cdata();
    const label n = a.size();
    for (label i = 0; i < n; ++i)
    {
        s += toDouble(ap[i])*toDouble(bp[i]);
    }
    return s;
}


Foam::reduceScalar Foam::doubleReduce::localSumSqr(const blockScalarUList& a)
{
    reduceScalar s = 0;
    const blockScalar* __restrict__ ap = a.cdata();
    const label n = a.size();
    for (label i = 0; i < n; ++i)
    {
        const reduceScalar v = toDouble(ap[i]);
        s += v*v;
    }
    return s;
}


Foam::reduceScalar Foam::doubleReduce::localSumMag(const blockScalarUList& a)
{
    reduceScalar s = 0;
    const blockScalar* __restrict__ ap = a.cdata();
    const label n = a.size();
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
    for (label i = 0; i < nCells; ++i)
    {
        for (label k = 0; k < blockDim; ++k)
        {
            const reduceScalar x = toDouble(ap[i*blockDim + k]);
            v[k] += x*x;
        }
    }
    parSum(v.data(), blockDim, comm);
    return v;
}


// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
-------------------------------------------------------------------------------
Description
    Templated field reductions of doubleReduce (amendment D4): input float
    or double, accumulator always double, #pragma omp simd reduction with
    the double accumulator, global sum through parSum (MPI_DOUBLE).

\*---------------------------------------------------------------------------*/

#include "PstreamReduceOps.H"
#include <cmath>
#include <type_traits>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace doubleReduce
{

namespace detail
{
    template<class T>
    constexpr void checkType()
    {
        static_assert
        (
            std::is_same<T, float>::value || std::is_same<T, double>::value,
            "doubleReduce: float or double input only"
        );
    }
}


template<class T>
reduceScalar localSum(const UList<T>& a)
{
    detail::checkType<T>();
    reduceScalar s = 0;
    const T* __restrict__ ap = a.cdata();
    const label n = a.size();
    #pragma omp simd reduction(+:s)
    for (label i = 0; i < n; ++i)
    {
        s += static_cast<reduceScalar>(ap[i]);
    }
    return s;
}


template<class T>
reduceScalar localSumMag(const UList<T>& a)
{
    detail::checkType<T>();
    reduceScalar s = 0;
    const T* __restrict__ ap = a.cdata();
    const label n = a.size();
    #pragma omp simd reduction(+:s)
    for (label i = 0; i < n; ++i)
    {
        s += std::abs(static_cast<reduceScalar>(ap[i]));
    }
    return s;
}


template<class T>
reduceScalar localSumSqr(const UList<T>& a)
{
    detail::checkType<T>();
    reduceScalar s = 0;
    const T* __restrict__ ap = a.cdata();
    const label n = a.size();
    #pragma omp simd reduction(+:s)
    for (label i = 0; i < n; ++i)
    {
        const reduceScalar x = static_cast<reduceScalar>(ap[i]);
        s += x*x;
    }
    return s;
}


template<class T>
reduceScalar localDot(const UList<T>& a, const UList<T>& b)
{
    detail::checkType<T>();
    reduceScalar s = 0;
    const T* __restrict__ ap = a.cdata();
    const T* __restrict__ bp = b.cdata();
    const label n = a.size();
    #pragma omp simd reduction(+:s)
    for (label i = 0; i < n; ++i)
    {
        s += static_cast<reduceScalar>(ap[i])*static_cast<reduceScalar>(bp[i]);
    }
    return s;
}


template<class W, class T>
reduceScalar localWeightedSum(const UList<W>& w, const UList<T>& a)
{
    detail::checkType<W>();
    detail::checkType<T>();
    reduceScalar s = 0;
    const W* __restrict__ wp = w.cdata();
    const T* __restrict__ ap = a.cdata();
    const label n = a.size();
    #pragma omp simd reduction(+:s)
    for (label i = 0; i < n; ++i)
    {
        s += static_cast<reduceScalar>(wp[i])*static_cast<reduceScalar>(ap[i]);
    }
    return s;
}


template<class T>
reduceScalar sum(const UList<T>& a, const label comm)
{
    return parSum(localSum(a), comm);
}


template<class T>
reduceScalar sumMag(const UList<T>& a, const label comm)
{
    return parSum(localSumMag(a), comm);
}


template<class T>
reduceScalar sumSqr(const UList<T>& a, const label comm)
{
    return parSum(localSumSqr(a), comm);
}


template<class T>
reduceScalar dot(const UList<T>& a, const UList<T>& b, const label comm)
{
    return parSum(localDot(a, b), comm);
}


template<class T>
reduceScalar norm2(const UList<T>& a, const label comm)
{
    // GUARD: sum of squares >= 0 by construction; clamp anyway for sqrt
    const reduceScalar s = sumSqr(a, comm);
    return std::sqrt(s > 0 ? s : 0);
}


template<class W, class T>
reduceScalar weightedSum
(
    const UList<W>& w,
    const UList<T>& a,
    const label comm
)
{
    return parSum(localWeightedSum(w, a), comm);
}


template<class T>
reduceScalar average(const UList<T>& a, const label comm)
{
    // One reduction for the sum and the count
    reduceScalar v[2] = {localSum(a), static_cast<reduceScalar>(a.size())};
    parSum(v, 2, comm);
    return (v[1] > 0 ? v[0]/v[1] : 0);
}

} // End namespace doubleReduce
} // End namespace Foam

// ************************************************************************* //

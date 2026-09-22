/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
-------------------------------------------------------------------------------
Namespace
    Foam::doubleReduce

Description
    Double-accumulated reductions over blockScalar (float) vectors, reduced
    with MPI_DOUBLE across all ranks of a communicator (spec 6, 9.2).

    Every dot product and norm in the block solvers goes through this
    namespace. Local sums are accumulated in reduceScalar (double); the
    parallel sum reduces the double value.

    The multi-value variants reduce several sums with one MPI call (used by
    the Krylov solvers to save latency).

SourceFiles
    doubleReduce.C

\*---------------------------------------------------------------------------*/

#ifndef coupledFoam_doubleReduce_H
#define coupledFoam_doubleReduce_H

#include "blockScalar.H"
#include "UPstream.H"
#include "FixedList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace doubleReduce
{

//- Local (this rank only) sum in double
reduceScalar localSum(const blockScalarUList& a);

//- Local dot product in double
reduceScalar localDot(const blockScalarUList& a, const blockScalarUList& b);

//- Local sum of squares in double
reduceScalar localSumSqr(const blockScalarUList& a);

//- Local sum of magnitudes in double
reduceScalar localSumMag(const blockScalarUList& a);

//- Sum a double over all ranks of comm
reduceScalar parSum(reduceScalar v, const label comm);

//- Sum n doubles over all ranks of comm in one MPI call
void parSum(reduceScalar* v, const label n, const label comm);

//- Global sum
reduceScalar sum
(
    const blockScalarUList& a,
    const label comm = UPstream::worldComm
);

//- Global dot product
reduceScalar dot
(
    const blockScalarUList& a,
    const blockScalarUList& b,
    const label comm = UPstream::worldComm
);

//- Global sum of squares
reduceScalar sumSqr
(
    const blockScalarUList& a,
    const label comm = UPstream::worldComm
);

//- Global L2 norm, sqrt of the double sum of squares
reduceScalar norm2
(
    const blockScalarUList& a,
    const label comm = UPstream::worldComm
);

//- Global sum of magnitudes
reduceScalar sumMag
(
    const blockScalarUList& a,
    const label comm = UPstream::worldComm
);

//- Two dot products with one reduction: (a.b, c.d)
FixedList<reduceScalar, 2> dot2
(
    const blockScalarUList& a,
    const blockScalarUList& b,
    const blockScalarUList& c,
    const blockScalarUList& d,
    const label comm = UPstream::worldComm
);

//- Per-component sums of squares of a 4-vector field, one reduction.
//  Result index = component (u, v, w, p).
FixedList<reduceScalar, 4> componentSumSqr
(
    const blockScalarUList& a,
    const label comm = UPstream::worldComm
);


// Field data (amendment D4, DECISIONS.md D-064): the same reductions for
// scalar lists of either precision (float in SP, double in DP). Double
// accumulator, SIMD reduction, MPI_DOUBLE. For a blockScalarUList argument
// the non-template overloads above are chosen (identical semantics).

//- Local sum in double
template<class T>
reduceScalar localSum(const UList<T>& a);

//- Local sum of magnitudes in double
template<class T>
reduceScalar localSumMag(const UList<T>& a);

//- Local sum of squares in double
template<class T>
reduceScalar localSumSqr(const UList<T>& a);

//- Local dot product in double
template<class T>
reduceScalar localDot(const UList<T>& a, const UList<T>& b);

//- Local weighted sum sum_i w_i a_i in double (e.g. sum V_P a_P)
template<class W, class T>
reduceScalar localWeightedSum(const UList<W>& w, const UList<T>& a);

//- Global sum
template<class T>
reduceScalar sum(const UList<T>& a, const label comm = UPstream::worldComm);

//- Global sum of magnitudes
template<class T>
reduceScalar sumMag
(
    const UList<T>& a,
    const label comm = UPstream::worldComm
);

//- Global sum of squares
template<class T>
reduceScalar sumSqr
(
    const UList<T>& a,
    const label comm = UPstream::worldComm
);

//- Global dot product
template<class T>
reduceScalar dot
(
    const UList<T>& a,
    const UList<T>& b,
    const label comm = UPstream::worldComm
);

//- Global L2 norm
template<class T>
reduceScalar norm2
(
    const UList<T>& a,
    const label comm = UPstream::worldComm
);

//- Global weighted sum sum_i w_i a_i (volume sum sum V_P a_P)
template<class W, class T>
reduceScalar weightedSum
(
    const UList<W>& w,
    const UList<T>& a,
    const label comm = UPstream::worldComm
);

//- Global average sum_i a_i / N (N the global size; 0 if N == 0), the
//  double-accumulated replacement of gAverage
template<class T>
reduceScalar average
(
    const UList<T>& a,
    const label comm = UPstream::worldComm
);

} // End namespace doubleReduce
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "doubleReduceTemplates.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

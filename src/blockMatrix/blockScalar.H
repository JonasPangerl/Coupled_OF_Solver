/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later (links against OpenFOAM, GPL-3.0)
-------------------------------------------------------------------------------
Description
    Precision types of the block library (see DECISIONS.md D-001).

    blockScalar   storage of block coefficients and of all Krylov/AMG
                  vectors. Always float, independent of WM_PRECISION_OPTION.
    reduceScalar  every accumulation (dot products, norms, Galerkin sums) and
                  every MPI reduction. Always double.

    The spec's "solveScalar" for block data means blockScalar here; the
    OpenFOAM solveScalar is never used for block data, so the library is
    correct under DP, SPDP and SP builds.

    Narrowing from scalar/reduceScalar to blockScalar happens only through
    toBlock(), which clamps |v| > clampValue and counts the clamps
    (spec 9.2, last item).

\*---------------------------------------------------------------------------*/

#ifndef coupledFoam_blockScalar_H
#define coupledFoam_blockScalar_H

#include "scalar.H"
#include "List.H"
#include "label.H"
#include <cmath>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

typedef floatScalar  blockScalar;
typedef doubleScalar reduceScalar;

typedef List<blockScalar>  blockScalarList;
typedef UList<blockScalar> blockScalarUList;
typedef List<reduceScalar> reduceScalarList;

//- Number of unknowns per cell (u, v, w, p)
static constexpr label blockDim = 4;

//- Number of coefficients per 4x4 block
static constexpr label blockSize = 16;

//- Index of the pressure component / continuity row
static constexpr label blockP = 3;


//- Narrow a double value to blockScalar with clamping (GUARD, spec 9.2)
inline blockScalar toBlock
(
    const doubleScalar v,
    const doubleScalar clampValue,
    label& nClamped
)
{
    // GUARD: clamp before narrowing, count every clamp
    if (v > clampValue)
    {
        ++nClamped;
        return static_cast<blockScalar>(clampValue);
    }
    if (v < -clampValue)
    {
        ++nClamped;
        return static_cast<blockScalar>(-clampValue);
    }
    return static_cast<blockScalar>(v);
}


//- As above, and count NaN input in nNonFinite: a NaN passes both clamp
//  comparisons. It is stored unmasked (the solve then fails through the B4
//  path) but reported (D-069 F9); +-inf is clamped and counted in nClamped.
inline blockScalar toBlock
(
    const doubleScalar v,
    const doubleScalar clampValue,
    label& nClamped,
    label& nNonFinite
)
{
    if (std::isnan(v))
    {
        ++nNonFinite;
        return static_cast<blockScalar>(v);
    }
    return toBlock(v, clampValue, nClamped);
}


//- Widen a blockScalar to double (explicit to keep -Wconversion quiet
//  and make every precision crossing grep-able)
inline doubleScalar toDouble(const blockScalar v)
{
    return static_cast<doubleScalar>(v);
}


//- Narrow without clamping, for values already known to be in range
//  (e.g. Krylov scalars that were computed from blockScalar data)
inline blockScalar narrow(const doubleScalar v)
{
    return static_cast<blockScalar>(v);
}

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

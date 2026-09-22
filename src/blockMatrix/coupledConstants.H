/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
-------------------------------------------------------------------------------
Description
    Typed guard constants (amendment D3, DECISIONS.md D-064).

    OpenFOAM's small/very-small/great guard constants depend on the
    precision option of the build (SP: 1e-6 and 1e-37 for the small ones), so a
    stabiliser written with them changes the behaviour between builds. The
    library uses these constants instead, typed by the quantity that is
    guarded:

        cfVSmall<T>()    denominator guard      float 1e-30,  double 1e-300
        cfSmallRel<T>()  relative guard         float 1e-6,   double 1e-14
        cfGreat<T>()     clamp / min-max start  float 1e30,   double 1e300

    T is the type of the guarded quantity: scalar (or blockScalar) for field
    data and block coefficients, doubleScalar (reduceScalar) for the
    quantities that are double in every build (D2.2 list: reductions,
    Krylov/GCR/Anderson small systems, dense LU, residual history, CFL,
    omega, eta, forces).

    The pivot guard (6.3) and the Rhie-Chow pseudo-inverse threshold (C2)
    are dictionary parameters (coupledDefaults) evaluated in double; they
    are not affected (D3).

\*---------------------------------------------------------------------------*/

#ifndef coupledFoam_coupledConstants_H
#define coupledFoam_coupledConstants_H

#include <type_traits>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- Denominator guard of a quantity of type T (1e-30 float, 1e-300 double)
template<class T>
constexpr T cfVSmall()
{
    static_assert
    (
        std::is_same<T, float>::value || std::is_same<T, double>::value,
        "cfVSmall: T must be float or double"
    );
    if constexpr (std::is_same<T, float>::value)
    {
        return 1e-30f;
    }
    else
    {
        return 1e-300;
    }
}


//- Relative guard of a quantity of type T (1e-6 float, 1e-14 double)
template<class T>
constexpr T cfSmallRel()
{
    static_assert
    (
        std::is_same<T, float>::value || std::is_same<T, double>::value,
        "cfSmallRel: T must be float or double"
    );
    if constexpr (std::is_same<T, float>::value)
    {
        return 1e-6f;
    }
    else
    {
        return 1e-14;
    }
}


//- Clamp / start value of a min-max search of type T (1e30 float,
//  1e300 double)
template<class T>
constexpr T cfGreat()
{
    static_assert
    (
        std::is_same<T, float>::value || std::is_same<T, double>::value,
        "cfGreat: T must be float or double"
    );
    if constexpr (std::is_same<T, float>::value)
    {
        return 1e30f;
    }
    else
    {
        return 1e300;
    }
}

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

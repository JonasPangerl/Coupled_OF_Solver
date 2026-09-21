/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
-------------------------------------------------------------------------------
Class
    Foam::lineSearch

Description
    Physicality line search (spec 7.2).

        omega_U = min_P (1, fU Uref / max(|dU_P|, VSMALL))
        omega_p = min_P (1, fp pref / max(|dp_P|, VSMALL))
        omega   = min(omega_U, omega_p), global minimum

    Uref = max(|U_inlet|, gMax|U|) at iteration 1, frozen; pref = Uref^2/2.
    |U_inlet| is taken as the maximum |U| over all physical boundary faces.

    If omega < omegaMin the caller cuts CFL by kappa and repeats the step,
    at most maxCflCuts times; then the step is taken with omegaMin and the
    cells with |dU| > fU Uref/omegaMin or |dp| > fp pref/omegaMin go to the
    dynamic remediation set (offendingCells).

    beta (>= 1, default 1): CFL boost on a full step. After an accepted
    step with omega == 1 and no CFL cut, the caller multiplies CFL by beta
    after the strategy update (ptcControl::boost: capped at CFLmax, not
    during a hold). beta = 1 is the spec default "no additional CFL boost on
    full step; the strategy handles growth" (7.2).

    Dictionary: coupled.lineSearch { fU fp omegaMin kappa maxCflCuts beta }.

SourceFiles
    lineSearch.C

\*---------------------------------------------------------------------------*/

#ifndef coupledFoam_lineSearch_H
#define coupledFoam_lineSearch_H

#include "blockScalar.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class lineSearch
{
    // Private Data

        scalar fU_;
        scalar fp_;
        scalar omegaMin_;
        scalar kappa_;
        label maxCflCuts_;
        scalar beta_;

        //- References (frozen after setReference)
        scalar Uref_;
        scalar pref_;
        bool refSet_;


public:

    // Constructors

        explicit lineSearch(const dictionary& coupledDict);


    // Member Functions

        scalar omegaMin() const noexcept { return omegaMin_; }
        scalar kappa() const noexcept { return kappa_; }
        label maxCflCuts() const noexcept { return maxCflCuts_; }
        scalar beta() const noexcept { return beta_; }
        scalar Uref() const noexcept { return Uref_; }
        scalar pref() const noexcept { return pref_; }
        bool referenceSet() const noexcept { return refSet_; }

        //- Set Uref from U at iteration 1 (and pref = Uref^2/2)
        void setReference(const volVectorField& U);

        //- Restore the references (restart)
        void setReference(const scalar Uref, const scalar pref);

        //- omega for an increment (4 values per cell), global
        scalar omega(const blockScalarUList& dx) const;

        //- Cells violating the omegaMin bound (local labels)
        labelList offendingCells(const blockScalarUList& dx) const;

        //- Effective settings
        void writeSettings(dictionary& dict) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

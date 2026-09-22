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

    Uref (D-050, keyword coupled.Uref, default boundary), frozen after
    iteration 1, pref = Uref^2/2:
        boundary  max |U| over the non-coupled boundary values (inflow,
                  moving walls, free stream); falls back to the field
                  maximum if the boundary maximum is below
                  coupledDefaults::UrefFallbackFactor times the field
                  maximum (no driving boundary velocity)
        field     max(|U_boundary|, gMax|U|) - the definition before D-050;
                  with potentialInit this is the potential-flow peak at
                  sharp edges (T4b: 8.8 U_inf), which effectively disabled
                  the dynamic set (cU Uref) and loosened the line search,
                  the local CFL limit and the sentinel
        <value>   an explicit velocity scale
    UrefSource() records which one was used (boundary, field,
    fieldFallback, explicit, restart).

    If omega < omegaMin the caller cuts CFL by kappa and repeats the step,
    at most maxCflCuts times; then the step is taken with omegaMin and the
    cells with |dU| > fU Uref/omegaMin or |dp| > fp pref/omegaMin go to the
    dynamic remediation set (offendingCells).

    beta (>= 1, default 1): CFL boost on a full step. After an accepted
    step with omega == 1 and no CFL cut, the caller multiplies CFL by beta
    after the strategy update (ptcControl::boost: capped at CFLmax, not
    during a hold). beta = 1 is the spec default "no additional CFL boost on
    full step; the strategy handles growth" (7.2).

    Dictionary: coupled.lineSearch { fU fp omegaMin kappa maxCflCuts beta },
    coupled.Uref.

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

        //- Uref selection (D-050): boundary | field | explicit
        word UrefMode_;
        scalar UrefExplicit_;

        //- References (frozen after setReference)
        scalar Uref_;
        scalar pref_;
        bool refSet_;
        word UrefSource_;


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
        const word& UrefMode() const noexcept { return UrefMode_; }
        const word& UrefSource() const noexcept { return UrefSource_; }

        //- Set Uref from U at iteration 1 according to coupled.Uref
        //  (and pref = Uref^2/2)
        void setReference(const volVectorField& U);

        //- Restore the references (restart)
        void setReference(const scalar Uref, const scalar pref);

        //- omega for an increment (4 values per cell), global
        scalar omega(const blockScalarUList& dx) const;

        //- Cells violating the omegaMin bound (local labels)
        labelList offendingCells(const blockScalarUList& dx) const;

        //- Diagnostics (rank-local, no reduction): cells whose increment
        //  violates the velocity / pressure limit at the given omega
        void countViolations
        (
            const blockScalarUList& dx,
            const scalar omega,
            label& nU,
            label& np
        ) const;

        //- Effective settings
        void writeSettings(dictionary& dict) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

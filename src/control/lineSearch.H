/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
-------------------------------------------------------------------------------
Class
    Foam::lineSearch

Description
    Physicality line search (spec 7.2).

        omega_U = min_P (1, fU Ustep / max(|dU_P|, cfVSmall))
        omega_p = min_P (1, fp pstep / max(|dp_P|, cfVSmall))
        omega   = min(omega_U, omega_p), global minimum

    with the step scale Ustep, pstep = Ustep^2/2 (D-057, below; Ustep = Uref
    and pstep = pref for coupled.UrefStep reference).

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

    Two scales (D-057). Uref/pref above classify the flow physically
    (dynamic set cU/cp/cSpike, increment clipping, sentinel, Anderson
    norm). The STEP limits of this class (omega, offendingCells) and the
    local CFL limit use Ustep/pstep = Ustep^2/2 (keyword coupled.UrefStep):
        reference    Ustep = Uref (D-050 behaviour, bitwise)
        fieldCapped  max(Uref, min(Ufield0, UrefStepCap Uref))
        field        max(Uref, Ufield0)
        <value>      explicit
    Ufield0 = max |U| over the cells and the non-coupled boundaries at
    iteration 1 (after coupled.potentialClip), UrefStepCap =
    coupled.UrefStepCap (default cU = 4).
    coupled.startupReference ramp: Uref and Ustep blend linearly from
    Ufield0 (start-up beta 0) to their frozen values (beta 1); exclude
    (default): if Ufield0 > Ustep, then while beta < 1 the step scale is
    Ufield0 and no cell enters the dynamic set (excludeDynamic; rollback
    marking stays); none: frozen values throughout. Both are no-ops when
    the initial field does not exceed the reference (T0-T3).

    If omega < omegaMin the caller cuts CFL by kappa and repeats the step,
    at most maxCflCuts times (a cut at CFLmin is not repeated, D-057); then
    the step is taken with omegaMin and the cells with
    |dU| > fU Ustep/omegaMin or |dp| > fp pstep/omegaMin go to the dynamic
    remediation set (offendingCells).

    beta (>= 1, default 1): CFL boost on a full step. After an accepted
    step with omega == 1 and no CFL cut, the caller multiplies CFL by beta
    after the strategy update (ptcControl::boost: capped at CFLmax, not
    during a hold). beta = 1 is the spec default "no additional CFL boost on
    full step; the strategy handles growth" (7.2).

    Dictionary: coupled.lineSearch { fU fp omegaMin kappa maxCflCuts beta },
    coupled.Uref, coupled.UrefStep, coupled.UrefStepCap,
    coupled.startupReference.

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
        doubleScalar omegaMin_;
        doubleScalar kappa_;
        label maxCflCuts_;
        doubleScalar beta_;

        //- Uref selection (D-050): boundary | field | explicit
        word UrefMode_;
        scalar UrefExplicit_;

        //- References (frozen after setReference)
        scalar Uref_;
        scalar pref_;
        bool refSet_;
        word UrefSource_;

        //- Step scale selection (D-057): reference | fieldCapped | field
        //  | explicit, cap factor, explicit value
        word stepMode_;
        scalar stepCap_;

        //- Boundary Uref below this fraction of the field maximum: field
        //  maximum (coupled.UrefFallbackFactor, D-050)
        scalar UrefFallbackFactor_;
        scalar UstepExplicit_;

        //- max |U| of the initial field and the non-coupled boundaries
        //  (after potentialClip), frozen with the references
        scalar Ufield0_;

        //- Step scale (frozen), pstep = Ustep^2/2
        scalar Ustep_;

        //- Start-up references (D-057): none | ramp | exclude
        word startupRef_;

        //- Effective values of the current iteration (setStartup)
        scalar UrefEff_;
        scalar prefEff_;
        scalar UstepEff_;
        scalar pstepEff_;
        bool excludeDynamic_;


    // Private Member Functions

        //- Ustep_ from stepMode_, Uref_ and Ufield0_; effective values of
        //  a finished start-up
        void updateStep();


public:

    // Constructors

        explicit lineSearch(const dictionary& coupledDict);


    // Member Functions

        doubleScalar omegaMin() const noexcept { return omegaMin_; }
        doubleScalar kappa() const noexcept { return kappa_; }
        label maxCflCuts() const noexcept { return maxCflCuts_; }
        doubleScalar beta() const noexcept { return beta_; }
        //- Classification reference of the current iteration (dynamic
        //  set, increment clipping, sentinel, Anderson norm); the frozen
        //  Uref except during a startupReference ramp
        scalar Uref() const noexcept { return UrefEff_; }
        scalar pref() const noexcept { return prefEff_; }

        //- Step scale of the current iteration (line search, offending
        //  cells, local CFL limit), D-057
        scalar Ustep() const noexcept { return UstepEff_; }
        scalar pstep() const noexcept { return pstepEff_; }

        //- The frozen values (restart state)
        scalar UrefBase() const noexcept { return Uref_; }
        scalar prefBase() const noexcept { return pref_; }
        scalar Ufield0() const noexcept { return Ufield0_; }
        scalar UstepBase() const noexcept { return Ustep_; }

        bool referenceSet() const noexcept { return refSet_; }
        const word& UrefMode() const noexcept { return UrefMode_; }
        const word& UrefSource() const noexcept { return UrefSource_; }
        const word& stepMode() const noexcept { return stepMode_; }
        const word& startupReference() const noexcept { return startupRef_; }

        //- True while startupReference exclude keeps cells out of the
        //  dynamic set (start-up not finished)
        bool excludeDynamic() const noexcept { return excludeDynamic_; }

        //- Set Uref from U at iteration 1 according to coupled.Uref
        //  (and pref = Uref^2/2), Ufield0 and Ustep
        void setReference(const volVectorField& U);

        //- Recompute Ufield0 (and Ustep) from U, e.g. after potentialClip;
        //  Uref is kept
        void setFieldScale(const volVectorField& U);

        //- Restore the references (restart); Ufield0 <= 0: unknown (a state
        //  written before D-057), Uref is used
        void setReference
        (
            const scalar Uref,
            const scalar pref,
            const scalar Ufield0
        );

        //- Effective references of an iteration from the start-up beta
        //  (startupReference ramp / exclude); none: the frozen values
        void setStartup(const scalar betaStartup, const bool startupDone);

        //- omega for an increment (4 values per cell), global
        doubleScalar omega(const blockScalarUList& dx) const;

        //- Cells violating the omegaMin bound (local labels)
        labelList offendingCells(const blockScalarUList& dx) const;

        //- Diagnostics (rank-local, no reduction): cells whose increment
        //  violates the velocity / pressure limit at the given omega
        void countViolations
        (
            const blockScalarUList& dx,
            const doubleScalar omega,
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

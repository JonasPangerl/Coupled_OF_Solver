/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
-------------------------------------------------------------------------------
Class
    Foam::gamgAutoTune

Description
    Cycle-efficiency controller for blockGAMG (amendment B3, 6.3.5).

    Input per outer iteration: the preconditioner efficiency
        rho = ||r0 - A M^-1 r0|| / ||r0||
    measured on the first preconditioner application of the linear solve
    (rho < 0: not measured, ignored). Every tuneInterval outer iterations
    the exact median rho_med of the window is evaluated:

      - rho_med > tuneRhoHigh (0.7):
            nPostSweeps += 1 (<= nPostSweepsMax); if already at the maximum
            and cycleType V or F: promote V -> F -> W (event
            "promotion-refused" if the coarsening-ratio rule refuses W);
            cycleType K: nothing further, event "kcycle-saturated".
      - rho_med < tuneRhoLow (0.3) and nPostSweeps > minPostSweeps (1):
            nPostSweeps -= 1. Whenever nPostSweeps == 1 (after the
            decrement or already before), rho_med < tuneRhoDemote (0.2)
            and the cycle type was promoted by the controller: demote one
            step (W -> F -> V, never below the user's cycle type).
      - Hysteresis: a change requires the same condition in
        tuneConsecutiveWindows (2) consecutive windows; after a change the
        count restarts.
      - Failure: rho_med > tuneRhoFail (0.9) in tuneConsecutiveWindows
        consecutive windows at
        nPostSweepsMax with cycleType W or K -> Warning, event "failure",
        failed() == true.

    Every change is logged as "GAMG-tune: <what> rho_med=<v>" and stored as
    an event (iteration, what, rho_med). Changes are only applied inside
    record(), i.e. between linear solves.

    Dictionary: solvers.coupled.blockGAMG { autoTune tuneInterval
    nPostSweepsMax }.

SourceFiles
    gamgAutoTune.C

\*---------------------------------------------------------------------------*/

#ifndef coupledFoam_gamgAutoTune_H
#define coupledFoam_gamgAutoTune_H

#include "blockGAMG.H"
#include "dictionary.H"
#include "DynamicList.H"
#include <string>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class gamgAutoTune
{
public:

    //- One controller event
    struct tuneEvent
    {
        label iter;
        std::string what;
        scalar rhoMed;
    };

private:

    //- Window classification
    enum class condition { none, high, low };

    // Private Data

        const blockGAMG& gamg_;

        bool autoTune_;
        label tuneInterval_;
        label nPostSweepsMax_;

        //- Decision thresholds of 6.3.5 (keywords tuneRhoHigh, tuneRhoLow,
        //  tuneRhoDemote, tuneRhoFail, tuneConsecutiveWindows,
        //  nPostSweepsMin; D-066)
        scalar rhoHigh_;
        scalar rhoLow_;
        scalar rhoDemote_;
        scalar rhoFail_;
        label consecutiveWindows_;
        label nPostSweepsMin_;

        //- Cycle type given by the user (lower bound for demotion)
        blockGAMG::cycleKind userCycle_;

        // State

            //- Measured rho values of the current window
            DynamicList<scalar> window_;

            //- Outer iterations recorded in the current window
            label nInWindow_;

            //- Condition of the previous window(s)
            condition prevCondition_;

            //- Consecutive windows with prevCondition_
            label nSameCondition_;

            //- Consecutive windows satisfying the failure condition
            label nFailWindows_;

            //- Promotion steps applied by the controller
            label nPromoted_;

            //- Latches to record repeated saturation events only once
            bool promotionRefused_;
            bool kSaturated_;

            //- Failure detected
            bool failed_;

            //- The last record() changed sweeps or cycle type
            bool changed_;

            //- Event history
            DynamicList<tuneEvent> events_;


    // Private Member Functions

        //- Exact median (copy + sort); requires a non-empty list
        static scalar median(const UList<scalar>& v);

        //- Evaluate a completed window
        void evaluate(const label iter, const scalar rhoMed);

        //- Apply the high / low rule
        void applyHigh(const label iter, const scalar rhoMed);
        void applyLow(const label iter, const scalar rhoMed);

        //- Cycle type one demotion step below the current one if the
        //  controller promoted it (not below the user's type), else the
        //  current type
        blockGAMG::cycleKind demotionTarget() const;

        //- Demote one step (if demotionTarget() differs)
        void demote(const label iter, const scalar rhoMed);

        //- Log and store an event
        void addEvent
        (
            const label iter,
            const std::string& what,
            const scalar rhoMed
        );


public:

    // Constructors

        gamgAutoTune(const blockGAMG& gamg, const dictionary& blockGAMGDict);


    // Member Functions

        bool active() const noexcept { return autoTune_; }
        bool failed() const noexcept { return failed_; }
        label tuneInterval() const noexcept { return tuneInterval_; }
        label nPostSweepsMax() const noexcept { return nPostSweepsMax_; }
        const UList<tuneEvent>& events() const noexcept { return events_; }

        //- Record the rho of outer iteration iter (rho < 0: not measured
        //  or not to be used - rejected/failed solves). Evaluates the
        //  window after every tuneInterval calls. Returns true if the
        //  sweeps or the cycle type were changed.
        bool record(const label iter, const scalar rho);

        //- Restart state (subdictionary gamgAutoTune): sweeps, cycle type,
        //  rho window, hysteresis and failure counters, latches.
        //  readState re-applies sweeps and cycle type to the GAMG.
        void writeState(dictionary& dict) const;
        void readState(const dictionary& dict);

        //- Events as a JSON array [{"iter":..,"what":"..","rhoMed":..}, ...]
        std::string eventsJson() const;

        //- Effective settings
        dictionary settings() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

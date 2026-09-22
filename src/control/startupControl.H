/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
-------------------------------------------------------------------------------
Class
    Foam::startupControl

Description
    Start-up phase of the convection scheme (spec 7.4, D-008, D-048): the
    global blending factor beta (0 = first-order upwind, 1 = the selected
    second-order scheme) as a function of the outer iteration. Remediated
    cells keep their own beta (remediation::beta).

    startupMode
        hybrid  beta = 0 until the ramp start n0, then
                beta_m = min(1, (m - n0)/startupRampLength), i.e. second
                order from n0 + startupRampLength on. n0 is the first
                accepted iteration n with one of the triggers
                  relativeR       R_n/R_1 < startupSwitchR (any n)
                  stagnation      n < startupRampStart and
                                  R_n > startupStagnationFactor R_(n-W)
                  fixed           n = startupRampStart and R_n is not
                                  dropping fast (R_n >= startupFastFactor
                                  R_(n-W))
                  delayedFastDrop startupRampStart < n < startupRampStartMax
                                  and R_n >= startupFastFactor R_(n-W)
                  cap             n = startupRampStartMax
                with W = startupStagnationWindow. While upwind still reduces
                R by more than the factor startupFastFactor over W
                iterations the ramp start is delayed, at most up to
                startupRampStartMax (startupRampStartMax =
                startupRampStart disables the delay).
        upwind  (default) the behaviour before D-048, exactly: beta = 0
                until n >= startupUpwindIters or R_n/R_1 < startupSwitchR,
                then beta = 1 from the next iteration on.
        developedStart (hybrid and upwind): on a fresh start without
                potentialInit and with a non-uniform initial U (a mapped
                or developed initial field), iteration 1 is assembled with
                beta = 1; if max(rU, rp) (the scale-free residuals
                sum|r|/sum(|Ax|+|b|) of the momentum and continuity rows,
                not relative to R_1) is below startupDevelopedTol, the
                start-up is skipped (beta = 1 from iteration 1), otherwise
                iteration 1 is re-assembled with beta = 0. No test case
                starts this way (T0: uniform U; T1-T5: potentialInit), so
                their upwind trajectories are unchanged.
        none    beta = 1 from iteration 1 (mapped or restarted-from-other-
                case initial fields).

    The ramp state (n0, trigger, R history) is part of coupledState, so a
    restart continues the ramp.

    Dictionary coupled: startupMode startupUpwindIters startupSwitchR
    startupRampStart startupRampLength startupRampStartMax
    startupStagnationFactor startupStagnationWindow startupFastFactor
    startupDevelopedTol.

SourceFiles
    startupControl.C

\*---------------------------------------------------------------------------*/

#ifndef coupledFoam_startupControl_H
#define coupledFoam_startupControl_H

#include "dictionary.H"
#include "scalarList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class startupControl
{
public:

    enum class mode { hybrid, upwind, none };

private:

    // Settings

        mode mode_;
        word modeName_;
        label upwindIters_;
        scalar switchR_;
        label rampStart_;
        label rampLength_;
        label rampStartMax_;
        bool stagnationTrigger_;
        scalar stagnationFactor_;
        label stagnationWindow_;
        scalar fastFactor_;
        scalar developedTol_;

    // State

        //- Ramp start iteration n0 (-1: not started)
        label rampStartIter_;

        //- Trigger that set n0 ("" while not started)
        word trigger_;

        //- R of the last accepted iterations (at most window + 1)
        scalarList Rhist_;

        //- Developed-start probe active in the current iteration
        bool probing_;

        //- beta = 1 regardless of the ramp (mode none, developed start,
        //  restart of a finished start-up written before D-048)
        bool full_;


public:

    // Constructors

        explicit startupControl(const dictionary& coupledDict);


    // Member Functions

        //- Global beta of iteration m
        scalar beta(const label m) const;

        //- Start-up finished in iteration m (beta = 1): the former
        //  startupDone flag (Eisenstat-Walker safeguard, phase)
        bool done(const label m) const { return beta(m) >= 1; }

        //- Iteration 1 of a fresh start: should the developed-start probe
        //  run? (hybrid or upwind mode, no potentialInit, non-uniform U).
        //  Sets the probing flag, beta(1) is then 1 until decideDeveloped().
        bool startProbe(const bool candidate);

        bool probing() const noexcept { return probing_; }

        //- Decide the probe from the scale-free residuals of iteration 1;
        //  returns true if the start is developed (keep beta = 1)
        bool decideDeveloped(const scalar rU, const scalar rp);

        //- Accepted iteration n with the normalised residual R: evaluate
        //  the triggers; returns true if the ramp (or jump) starts at n
        bool update(const label n, const scalar R);

        label rampStartIter() const noexcept { return rampStartIter_; }
        label rampEndIter() const;
        const word& trigger() const noexcept { return trigger_; }
        const word& modeName() const noexcept { return modeName_; }

        void writeState(dictionary& dict) const;

        //- Restart; states written before D-048 carry only startupDone
        void readState(const dictionary& dict, const label iter);

        void writeSettings(dictionary& dict) const;
};


} // End namespace Foam

#endif

// ************************************************************************* //

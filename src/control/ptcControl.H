/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
-------------------------------------------------------------------------------
Class
    Foam::ptcControl

Description
    Pseudo-transient continuation (spec 5.4, 7.1, 7.3).

    Local time step (5.4), as V_P/dt_P:
        V/dt = max(1/2 sum_f |phi_f| + nuEff_P V_P^(1/3), VSMALL)
               / (CFL * cflFactor_P)
    cflFactor_P from remediation (static 0.5, dynamic 0.1, else 1).

    Solution-limited local CFL (7.3): with the momentum residual r_P of the
    current assembly, dU_P = |r_P| dt_P / V_P; if dU_P > fLoc*Uref the cell's
    dt is reduced by fLoc*Uref/dU_P.

    Implicit estimate (coupled.localLimit.implicit, D-049): dU_P =
    |r_P| / (V_P/dt_P + a_P) with a_P the mean momentum diagonal of the cell
    (convection and diffusion, without PTC) - the update the implicit step
    makes in the cell with its neighbours frozen (Jacobi estimate). The
    explicit estimate of 7.3 overestimates the update of a stiff cell by
    the factor 1 + a_P dt_P/V_P (hundreds in wall-adjacent cells at CFL
    500), so it limits cells whose implicit update is small.

    Local limiter with memory (coupled.localLimit.memory, D-049): a per-cell
    factor f_P in (0, 1] multiplies dt_P before the check. When the cell is
    limited, f_P is reduced to the factor the check needs. Every line-search
    retrial starts from the factor of the start of the outer iteration and
    the last trial's factor is kept; one outer iteration counts at most one
    limit event per cell. When it is not limited, f_P recovers at the start
    of the next outer iteration by the factor localRecovery per iteration,
    after a hold of localHold iterations since the last limit event. A cell that was limited
    in localStickyAfter iterations (0 = never) keeps its factor for the rest
    of the run: it is never released. memory no is the memoryless limiter of
    7.3 (bit-identical to the implementation before D-049).

    Global CFL strategies (7.1), R_n = ||r_n||/||r_1||:
        mRDM  if R_n <= R_(n-1):
                  CFL <- min(CFLmax, CFL*max(1, min(betaMax, (R_(n-1)/R_n)^gamma)))
              else unchanged
        EXP   CFL <- min(CFLmax, CFL*betaExp)
        SER   CFL <- clamp(CFL*R_(n-1)/R_n, CFLmin, CFLmax)
    After any decrease (line search, sentinel) the strategy may not
    increase CFL for nHold iterations.

    Dictionary: coupled.ptc { cflStrategy CFL0 CFLmin CFLmax gamma betaMax
    betaExp nHold }, coupled.localLimit { enabled fLoc memory localRecovery
    localHold localStickyAfter }.

SourceFiles
    ptcControl.C

\*---------------------------------------------------------------------------*/

#ifndef coupledFoam_ptcControl_H
#define coupledFoam_ptcControl_H

#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class ptcControl
{
public:

    enum class strategy { mRDM, EXP, SER };

private:

    // Private Data

        const fvMesh& mesh_;

        strategy strategy_;
        word strategyName_;
        scalar CFL0_;
        scalar CFLmin_;
        scalar CFLmax_;
        scalar gamma_;
        scalar betaMax_;
        scalar betaExp_;
        label nHold_;

        bool localLimit_;
        scalar fLoc_;
        bool localImplicit_;

        // Local limiter memory (D-049)
        bool localMemory_;
        scalar localRecovery_;
        label localHold_;
        label localStickyAfter_;

        // State

            //- Current global CFL
            scalar CFL_;

            //- R of the previous iteration (-1: none yet)
            scalar Rprev_;

            //- Iterations left in the hold
            label holdRemaining_;

            //- Cells limited by 7.3 in the last call (global)
            label nLocalLimited_;

            //- Cells with a memory factor below 1 after the last call
            //  (global; 0 without memory)
            label nLocalThrottled_;

            //- Cells whose factor is permanently kept (global)
            label nLocalSticky_;

            //- Per-cell dt factor f_P in (0, 1] (empty without memory)
            scalarField localFactor_;

            //- f_P at the start of the outer iteration (base of every
            //  line-search trial)
            scalarField localFactor0_;

            //- Cell limited in a trial of the current outer iteration
            boolList localLimitedNow_;

            //- Per-cell iterations left before the factor may recover
            labelList localHoldRemaining_;

            //- Per-cell number of limit events (iterations limited)
            labelList localCount_;


public:

    // Constructors

        ptcControl(const fvMesh& mesh, const dictionary& coupledDict);


    // Member Functions

        scalar CFL() const noexcept { return CFL_; }
        scalar CFLmin() const noexcept { return CFLmin_; }
        scalar Rprev() const noexcept { return Rprev_; }
        label holdRemaining() const noexcept { return holdRemaining_; }
        label nHold() const noexcept { return nHold_; }
        label nLocalLimited() const noexcept { return nLocalLimited_; }
        label nLocalThrottled() const noexcept { return nLocalThrottled_; }
        label nLocalSticky() const noexcept { return nLocalSticky_; }
        bool localMemory() const noexcept { return localMemory_; }
        bool localImplicit() const noexcept { return localImplicit_; }
        const word& strategyName() const noexcept { return strategyName_; }

        //- Per-cell memory factor (empty field without memory)
        const scalarField& localFactor() const noexcept
        {
            return localFactor_;
        }

        //- Per-cell limit-event counts (empty without memory)
        const labelList& localCount() const noexcept { return localCount_; }

        //- V_P/dt_P per cell (5.4)
        tmp<scalarField> rDeltaTV
        (
            const surfaceScalarField& phi,
            const volScalarField& nuEff,
            const scalarField& cflFactor
        ) const;

        //- Start of an outer iteration: recovery step of the local limiter
        //  memory (no effect without memory). Called once per outer
        //  iteration, before the first assembly.
        void beginIteration();

        //- Apply the solution-limited local CFL (7.3) in place, with the
        //  memory factors if enabled. aMom is the mean momentum diagonal
        //  (implicit estimate only). Returns the global number of limited
        //  cells. If ratio is given (level-3 diagnostics), it receives
        //  dU_P/(fLoc Uref) per cell before the limit (1 = at the limit).
        label applyLocalLimit
        (
            scalarField& rDeltaTV,
            const vectorField& rMom,
            const scalarField& aMom,
            const scalar Uref,
            scalarField* ratio = nullptr
        );

        //- Strategy update after an iteration with normalised residual R
        void update(const scalar R);

        //- Decrease CFL by factor (line search kappa, sentinel), floor CFLmin,
        //  start the hold
        void decrease(const scalar factor);

        //- Increase CFL by factor >= 1 (line-search beta on a full step,
        //  7.2), capped at CFLmax; no effect during a hold
        void boost(const scalar factor);

        //- Restart state
        void writeState(dictionary& dict) const;
        void readState(const dictionary& dict);

        //- Effective settings
        void writeSettings(dictionary& dict) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

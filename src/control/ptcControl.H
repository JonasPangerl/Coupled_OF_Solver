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

    Global CFL strategies (7.1), R_n = ||r_n||/||r_1||:
        mRDM  if R_n <= R_(n-1):
                  CFL <- min(CFLmax, CFL*max(1, min(betaMax, (R_(n-1)/R_n)^gamma)))
              else unchanged
        EXP   CFL <- min(CFLmax, CFL*betaExp)
        SER   CFL <- clamp(CFL*R_(n-1)/R_n, CFLmin, CFLmax)
    After any decrease (line search, sentinel) the strategy may not
    increase CFL for nHold iterations.

    Dictionary: coupled.ptc { cflStrategy CFL0 CFLmin CFLmax gamma betaMax
    betaExp nHold }, coupled.localLimit { enabled fLoc }.

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

        // State

            //- Current global CFL
            scalar CFL_;

            //- R of the previous iteration (-1: none yet)
            scalar Rprev_;

            //- Iterations left in the hold
            label holdRemaining_;

            //- Cells limited by 7.3 in the last call (global)
            label nLocalLimited_;


public:

    // Constructors

        ptcControl(const fvMesh& mesh, const dictionary& coupledDict);


    // Member Functions

        scalar CFL() const noexcept { return CFL_; }
        scalar CFLmin() const noexcept { return CFLmin_; }
        scalar Rprev() const noexcept { return Rprev_; }
        label holdRemaining() const noexcept { return holdRemaining_; }
        label nLocalLimited() const noexcept { return nLocalLimited_; }
        const word& strategyName() const noexcept { return strategyName_; }

        //- V_P/dt_P per cell (5.4)
        tmp<scalarField> rDeltaTV
        (
            const surfaceScalarField& phi,
            const volScalarField& nuEff,
            const scalarField& cflFactor
        ) const;

        //- Apply the solution-limited local CFL (7.3) in place.
        //  Returns the global number of limited cells.
        label applyLocalLimit
        (
            scalarField& rDeltaTV,
            const vectorField& rMom,
            const scalar Uref
        );

        //- Strategy update after an iteration with normalised residual R
        void update(const scalar R);

        //- Decrease CFL by factor (line search kappa, sentinel), floor CFLmin,
        //  start the hold
        void decrease(const scalar factor);

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

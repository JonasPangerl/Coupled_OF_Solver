/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
-------------------------------------------------------------------------------
Class
    Foam::sfdControl

Description
    Selective Frequency Damping (amendment C3, spec 7.6; Akervik et al.,
    Phys. Fluids 18, 2006). Suppresses limit-cycle oscillations of the
    steady iteration by a forcing towards a low-pass filtered velocity:

        chi*   = chi Uref/Lref      [1/s]
        Delta* = Delta Lref/Uref    [s]

    Momentum rows (5.3a, before PTC), set up by coupledAssembler:
        a_P += chi* V_P  on the three momentum diagonals
        b_P += chi* V_P Ubar_P
    i.e. the forcing -chi*(U - Ubar), which vanishes at the steady state.
    The Rhie-Chow abar (5.3d) excludes the SFD term (D-052): D = V/abar would
    otherwise depend on chi and a run converged with SFD still active would
    carry a chi-dependent pressure dissipation.

    Filter, after every accepted outer step, with the local pseudo-time step
    dt_P = V_P/(V_P/dt_P) of the accepted assembly (after the local limit):
        Ubar_P <- Ubar_P + dt_P/(Delta* + dt_P) (U_P - Ubar_P)

    Ubar (volVectorField USFD) is initialised to U at activation (first
    iteration), reset to U whenever the Anderson history is flushed (7.5,
    includes the rollback of 9.3; resetOnFlush no keeps the filter state as
    a variant), written at write times (C6) and read back on restart
    together with the coupledState entries sfdOn, sfdQuiet.

    Deactivation: R_n < deactivateBelowR in nHold consecutive accepted
    iterations sets chi* = 0 for the rest of the run (log "SFD-off"); Ubar
    stays in the state.

    Dictionary: coupled.sfd { enabled chi Delta Lref deactivateBelowR
    resetOnFlush }.

SourceFiles
    sfdControl.C

\*---------------------------------------------------------------------------*/

#ifndef coupledFoam_sfdControl_H
#define coupledFoam_sfdControl_H

#include "fvMesh.H"
#include "volFields.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class sfdControl
{
    // Private Data

        const fvMesh& mesh_;

        // Settings

        bool enabled_;
        scalar chi_;
        scalar Delta_;
        scalar Lref_;
        scalar deactivateBelowR_;
        bool resetOnFlush_;
        label nHold_;

        // State

            //- Filtered velocity Ubar (allocated only if enabled)
            autoPtr<volVectorField> UbarPtr_;

            //- chi*, Delta* (set by setReference)
            scalar chiStar_;
            scalar DeltaStar_;

            //- Ubar initialised
            bool initialised_;

            //- Forcing on (false after SFD-off)
            bool on_;

            //- Consecutive accepted iterations with R < deactivateBelowR
            label quiet_;

            //- Number of resets (Anderson flushes) so far
            label nResets_;


public:

    // Constructors

        //- Construct from the coupled dictionary; nHold is ptc.nHold
        sfdControl
        (
            const fvMesh& mesh,
            const dictionary& coupledDict,
            const label nHold
        );


    // Member Functions

        bool enabled() const noexcept { return enabled_; }

        //- Forcing active (enabled, reference set, not switched off)
        bool active() const noexcept { return enabled_ && on_; }

        //- chi* if active, else 0
        scalar chiStar() const noexcept { return active() ? chiStar_ : 0; }
        scalar DeltaStar() const noexcept { return DeltaStar_; }
        label nResets() const noexcept { return nResets_; }

        //- Filtered velocity (enabled only)
        const volVectorField& Ubar() const;

        //- chi*, Delta* from Uref (7.2); prints the parameter guidance
        void setReference(const scalar Uref);

        //- Initialise Ubar = U unless initialised or restored (start of an
        //  outer iteration)
        void begin(const volVectorField& U);

        //- Anderson flush: Ubar = U (if active and resetOnFlush)
        void reset(const volVectorField& U);

        //- Filter update after an accepted step; rDeltaTV = V/dt per cell
        //  of the accepted assembly
        void update(const volVectorField& U, const scalarField& rDeltaTV);

        //- Deactivation check after an accepted step (true: switched off now)
        bool checkOff(const scalar R, const label iter);

        //- max_P |U_P - Ubar_P| / Uref (rank-local; 0 if not enabled)
        scalar maxDeviation(const volVectorField& U, const scalar Uref) const;

        //- Write USFD at a write time (enabled only)
        void write() const;

        //- Restart state
        void writeState(dictionary& dict) const;
        void readState(const dictionary& dict);

        //- Effective settings (11)
        void writeSettings(dictionary& dict) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //

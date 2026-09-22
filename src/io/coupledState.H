/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
-------------------------------------------------------------------------------
Class
    Foam::coupledState

Description
    Restart state (spec 10): a localIOdictionary "coupledState" in the time
    directory (per processor directory in parallel; a global IOdictionary
    would be read by the master and broadcast, giving every rank the state
    of rank 0) with
        CFL, iter, startupDone, start-up ramp state, R1, Uref, pref
        (frozen at iteration 1, D-050), nHoldRemaining, Rprev,
        dynamicSet, dynamicSetAge, nRollbacks, Eisenstat-Walker state,
        gamgAutoTune (sweeps, cycle type, rho window, hysteresis, latches),
        convergenceMonitor (last forceCoeffsWindow Cd/Cl samples),
        refinementHistory (empty placeholder for Phase 2 AMR)
    plus the volScalarField coupledD (Rhie-Chow D = V/abar of the last
    assembly) and the surfaceScalarField coupledQ (explicit Rhie-Chow face
    term of the last assembly), from which phi is recomputed on restart for
    the consistency check (spec 10; ||phi_read -
    phi_recomputed||_inf/||phi||_inf). The Anderson history is not part of
    the state (D-026).

    The dictionary is always written in ASCII with round-trip precision
    (max_digits10), independent of writeFormat/writePrecision: a binary
    dictionary does not read back an empty labelList, and the scalars must
    restart exactly.

    Atomic write (DECISIONS.md D-017): the state is written to
    "coupledState.tmp" and renamed to "coupledState" (per processor
    directory with the uncollated file handler). It is written after all
    fields of the time step, so its presence marks a complete time
    directory. With other file handlers the rename is skipped and logged.

SourceFiles
    coupledState.C

\*---------------------------------------------------------------------------*/

#ifndef coupledFoam_coupledState_H
#define coupledFoam_coupledState_H

#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class coupledState
{
    // Private Data

        const fvMesh& mesh_;

        //- Write the state at write times
        bool writeState_;


public:

    // Constructors

        coupledState(const fvMesh& mesh, const dictionary& coupledDict);


    // Member Functions

        bool writeState() const noexcept { return writeState_; }

        //- Read the state of the current time; false if absent
        bool read(dictionary& state) const;

        //- Write the state for the current time (atomic, see above)
        void write(const dictionary& state) const;

        //- Write the Rhie-Chow D field for the current time
        void writeD(const volScalarField& D) const;

        //- Read coupledD of the current time if present (null tmp otherwise)
        tmp<volScalarField> readD() const;

        //- Write the explicit Rhie-Chow face term q_f of the last assembly
        //  (the flux of the time directory was built with it, so the
        //  restart consistency check needs it, spec 10)
        void writeQ(const surfaceScalarField& q) const;

        //- Read coupledQ of the current time if present (null tmp otherwise)
        tmp<surfaceScalarField> readQ() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

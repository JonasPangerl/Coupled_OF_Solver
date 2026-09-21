/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
-------------------------------------------------------------------------------
Class
    Foam::coupledState

Description
    Restart state (spec 10): an IOdictionary "coupledState" in the time
    directory with
        CFL, iter, startupDone, R1, Uref, pref, nHoldRemaining, Rprev,
        dynamicSet, dynamicSetAge, nRollbacks,
        refinementHistory (empty placeholder for Phase 2 AMR)
    plus the volScalarField coupledD (Rhie-Chow D = V/abar of the last
    assembly), from which phi is recomputed on restart for the consistency
    check (spec 10; ||phi_read - phi_recomputed||_inf/||phi||_inf).

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
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
-------------------------------------------------------------------------------
Class
    Foam::sentinel

Description
    Finite checks, rollback and last-valid output (spec 9.3).

    store() keeps copies of U, p, phi and the turbulence fields (k, omega,
    nut, when present) of the previous iteration. check() computes in one
    pass gMin/gMax of p, |U|, k, omega and the number of non-finite values.
    Fail condition: nNonFinite > 0, gMax|U| > UFactor*Uref or
    max|p| > pFactor*pref.

    On failure the caller restores the copies (restore()), cuts
    CFL <- max(CFLmin, cflFactor*CFL), marks the offending cells into the
    dynamic remediation set and counts a rollback. After maxRollbacks
    consecutive rollbacks, writeLastValid() writes the restored fields to
    time "<n>_lastValid" together with a volScalarField sentinelFlag of the
    offending cells (the caller adds remediationFlag and the remediation
    cellSets into the same instance), prints the
    offending cell centres (at most maxReport) and the run exits with
    FatalError. A time directory with non-finite values is never written.
    <n> is the outer iteration whose fields are written, i.e. the last
    accepted one (iter - 1 when iteration iter is aborted), see
    lastValidName(). The objects are written through the file handler
    directly (writeInstance): regIOobject::writeObject resets an instance
    that is not the current time name to time().timeName() (v2606), which
    put the fields of iteration n into the directory of iteration n + 1.

    Dictionary: coupled.sentinel { maxRollbacks UFactor pFactor cflFactor
    maxReport }.

SourceFiles
    sentinel.C

\*---------------------------------------------------------------------------*/

#ifndef coupledFoam_sentinel_H
#define coupledFoam_sentinel_H

#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class sentinel
{
public:

    //- Result of one check
    struct checkResult
    {
        bool ok = true;
        label nNonFinite = 0;
        scalar maxMagU = 0;
        scalar minP = 0;
        scalar maxP = 0;
        scalar minK = 0;
        scalar maxK = 0;
        scalar minOmega = 0;
        scalar maxOmega = 0;
        labelList offending;        // local cell labels
    };


private:

    // Private Data

        const fvMesh& mesh_;

        label maxRollbacks_;
        scalar UFactor_;
        scalar pFactor_;
        scalar cflFactor_;
        label maxReport_;

        //- Consecutive and total rollbacks
        label consecutive_;
        label nRollbacks_;

        //- Copies of the previous iteration
        vectorField U0_;
        scalarField p0_;
        scalarField phi0_;
        List<scalarField> phi0Boundary_;
        scalarField k0_;
        scalarField omega0_;
        scalarField nut0_;
        bool stored_;


public:

    // Constructors

        sentinel(const fvMesh& mesh, const dictionary& coupledDict);


    // Member Functions

        scalar cflFactor() const noexcept { return cflFactor_; }
        label nRollbacks() const noexcept { return nRollbacks_; }
        label consecutive() const noexcept { return consecutive_; }
        label maxRollbacks() const noexcept { return maxRollbacks_; }

        //- Keep copies (call before the field update)
        void store
        (
            const volVectorField& U,
            const volScalarField& p,
            const surfaceScalarField& phi,
            const volScalarField* k,
            const volScalarField* omega,
            const volScalarField* nut
        );

        //- Check the fields
        checkResult check
        (
            const volVectorField& U,
            const volScalarField& p,
            const volScalarField* k,
            const volScalarField* omega,
            const scalar Uref,
            const scalar pref
        ) const;

        //- Restore the copies, count the rollback
        void restore
        (
            volVectorField& U,
            volScalarField& p,
            surfaceScalarField& phi,
            volScalarField* k,
            volScalarField* omega,
            volScalarField* nut
        );

        //- A successful iteration resets the consecutive counter
        void accepted() noexcept
        {
            consecutive_ = 0;
        }

        //- True if the consecutive rollback limit is exceeded
        bool exhausted() const noexcept
        {
            return consecutive_ >= maxRollbacks_;
        }

        //- Instance name "<validIter>_lastValid" of the 9.3 abort path;
        //  validIter is the iteration whose fields are written
        static word lastValidName(const label validIter)
        {
            return word(Foam::name(validIter) + "_lastValid");
        }

        //- Write a registered object into its own instance directory,
        //  bypassing the time-name redirection of regIOobject::writeObject
        static bool writeInstance(const regIOobject& io);

        //- Write the last valid fields (those of iteration validIter) to
        //  lastValidName(validIter), sentinelFlag, report (9.3)
        void writeLastValid
        (
            const volVectorField& U,
            const volScalarField& p,
            const surfaceScalarField& phi,
            const volScalarField* k,
            const volScalarField* omega,
            const volScalarField* nut,
            const labelList& offending,
            const label validIter
        ) const;

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

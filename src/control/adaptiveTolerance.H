/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
-------------------------------------------------------------------------------
Class
    Foam::adaptiveTolerance

Description
    Eisenstat-Walker adaptive inner tolerance (amendment B2; Eisenstat and
    Walker, SIAM J. Sci. Comput. 17 (1996), Choice 2).

    Per outer iteration n with the normalised residual R_n (7.1):
        eta_n = gammaEW (R_n / R_(n-1))^alphaEW
    safeguards, in this order:
        1. if gammaEW eta_(n-1)^alphaEW > 0.1:
               eta_n = max(eta_n, gammaEW eta_(n-1)^alphaEW)
        2. eta_n = clamp(eta_n, etaMin, etaMax)
        3. during start-up (7.4): eta_n = etaMax
        4. first iteration: eta_n = etaMax
    The linear solver then uses relTol = eta_n; its absolute tolerance
    stays the floor. With adaptiveRelTol no, relTol of the solver dictionary
    is used unchanged.

    Dictionary: solvers.coupled { adaptiveRelTol etaMin etaMax gammaEW
    alphaEW relTol }.

SourceFiles
    adaptiveTolerance.C

\*---------------------------------------------------------------------------*/

#ifndef coupledFoam_adaptiveTolerance_H
#define coupledFoam_adaptiveTolerance_H

#include "dictionary.H"
#include "scalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class adaptiveTolerance
{
    // Private Data

        bool enabled_;
        doubleScalar etaMin_;
        doubleScalar etaMax_;
        doubleScalar gamma_;
        doubleScalar alpha_;
        doubleScalar safeguard_;
        doubleScalar relTolFixed_;

        //- Previous eta and R (-1: none)
        doubleScalar etaPrev_;
        doubleScalar Rprev_;

        //- Diagnostics of the last eta() call: the raw Eisenstat-Walker
        //  value (-1: not evaluated) and the rule that set the result
        doubleScalar lastRaw_;
        const char* lastClip_;


public:

    // Constructors

        explicit adaptiveTolerance(const dictionary& linearDict);


    // Member Functions

        bool enabled() const noexcept { return enabled_; }

        //- eta for the current outer iteration given R_n (call once per
        //  assembly; a repeated step after a CFL cut passes the same R)
        doubleScalar eta(const doubleScalar R, const bool startupDone);

        //- Raw EW value of the last eta() call (-1: not evaluated)
        doubleScalar lastRaw() const noexcept { return lastRaw_; }

        //- Rule that set the last eta: none | fixed | noHistory |
        //  safeguard1 | etaMin | etaMax
        const char* lastClip() const noexcept { return lastClip_; }

        //- Accept the iteration: R becomes R_(n-1) of the next one
        void accept(const doubleScalar R, const doubleScalar eta);

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

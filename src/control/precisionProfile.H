/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
-------------------------------------------------------------------------------
Class
    Foam::precisionProfile

Description
    Precision profile of the default values (amendment D7, DECISIONS.md
    D-064).

    Keyword coupled.precisionProfile (fvSolution):
        auto   (default) sp if sizeof(scalar) == 4, else dp
        dp     the defaults of the DP/SPDP builds (coupledDefaults)
        sp     the single-precision defaults (coupledDefaults::spProfile)

    The profile only changes DEFAULT values; a keyword that is set
    explicitly in the case always wins:

        keyword                              dp      sp
        coupled.convergence.residualTol      1e-6    1e-5
        solvers.coupled.tolerance            1e-8    1e-6
        solvers.coupled.etaMin               1e-3    1e-2
        coupled.bounds.omegaMin              1e-6    1e-5
        coupled.bounds.kMin                  1e-12   1e-10
        Anderson above anderson.maxCells     no      yes, m <= 4

    select() is called once by coupledFoam before any consumer is
    constructed and logs the profile; consumers read current(). Without a
    select() (test applications) current() is the auto profile.

SourceFiles
    precisionProfile.C

\*---------------------------------------------------------------------------*/

#ifndef coupledFoam_precisionProfile_H
#define coupledFoam_precisionProfile_H

#include "dictionary.H"
#include "word.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class precisionProfile
{
public:

    //- Default values of one profile
    struct values
    {
        word name;
        doubleScalar residualTol;
        doubleScalar tolerance;
        doubleScalar etaMin;
        doubleScalar boundOmegaMin;
        doubleScalar kMin;
        //- Anderson stays enabled above anderson.maxCells ...
        bool andersonAboveMaxCells;
        //- ... with m clamped to at most this value
        label andersonLargeMaxM;
    };


private:

    // Private Static Data

        //- The selected profile (auto until select() is called)
        static values current_;

        //- Keyword value given to select() (auto | dp | sp)
        static word keyword_;


public:

    // Static Member Functions

        //- Profile name of the auto rule: sp if sizeof(scalar) == 4
        static word autoName();

        //- Values of a named profile (dp | sp); FatalError otherwise
        static values named(const word& name);

        //- Resolve coupled.precisionProfile (auto | dp | sp), make it the
        //  current profile and log it once (master)
        static const values& select(const dictionary& coupledDict);

        //- The current profile
        static const values& current() noexcept { return current_; }

        //- Keyword value of the last select() ("auto" without one)
        static const word& keyword() noexcept { return keyword_; }

        //- Profile and its default values as a dictionary (effective
        //  settings)
        static void writeSettings(dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

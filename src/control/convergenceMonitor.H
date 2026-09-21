/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
-------------------------------------------------------------------------------
Class
    Foam::convergenceMonitor

Description
    Solver-internal stop criterion (spec 12.3, DECISIONS.md D-005).

      (i)  R_n < residualTol
      (ii) with a forceCoeffs function object present: over the last
           forceCoeffsWindow iterations
               max(Cd) - min(Cd) <= forceCoeffsTol*|mean(Cd)|
           and the same for Cl.

    mode any (default): stop when (i) holds, or (ii) holds if available.
    mode all: stop when (i) and, if available, (ii) hold.

    Cd and Cl are read from the results of the first forceCoeffs function
    object found in the function-object properties (native
    setResult("Cd"/"Cl")).

    Dictionary: coupled.convergence { mode forceCoeffsWindow forceCoeffsTol
    residualTol }.

SourceFiles
    convergenceMonitor.C

\*---------------------------------------------------------------------------*/

#ifndef coupledFoam_convergenceMonitor_H
#define coupledFoam_convergenceMonitor_H

#include "Time.H"
#include "DynamicList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class convergenceMonitor
{
    // Private Data

        word mode_;
        label window_;
        scalar forceTol_;
        scalar residualTol_;

        DynamicList<scalar> Cd_;
        DynamicList<scalar> Cl_;

        //- Name of the forceCoeffs object used (empty: none found yet)
        word foName_;


    // Private Member Functions

        //- Window criterion on one history
        bool windowConverged(const DynamicList<scalar>& h) const;


public:

    // Constructors

        explicit convergenceMonitor(const dictionary& coupledDict);


    // Member Functions

        scalar residualTol() const noexcept { return residualTol_; }

        //- Append the current Cd, Cl (if a forceCoeffs object exists)
        void record(const Time& runTime);

        //- True if force coefficients are being monitored
        bool haveForces() const noexcept { return !foName_.empty(); }

        //- Last Cd, Cl (0 if none)
        scalar Cd() const { return Cd_.size() ? Cd_.last() : 0; }
        scalar Cl() const { return Cl_.size() ? Cl_.last() : 0; }

        //- Evaluate the criterion for the normalised residual R
        bool converged(const scalar R) const;

        //- Effective settings
        void writeSettings(dictionary& dict) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

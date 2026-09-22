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
           and the same for Cl
           OR, only if forceCoeffsDriftTol > 0 (default 0: off), the
           stationary-mean rule (D-045 b, the rule of D-042): over the last
           forceCoeffsRmsWindow samples, for EVERY monitored coefficient
           (Cd, Cl, and Cm = CmPitch if the object provides it) the means
           of the older and the newer half-window differ by at most
               max(forceCoeffsDriftTol*|window mean|, forceCoeffsDriftAbs)

    mode any (default): stop when (i) holds, or (ii) holds if available.
    mode all: stop when (i) and, if available, (ii) hold.

    Window statistics (displayed on the CF| line, in the diagnostics
    record and in summary.json; D-045 b): per monitored coefficient over
    the last min(forceCoeffsRmsWindow, n) samples the mean, the RMS
    fluctuation (population standard deviation) and the half-window mean
    drift |m_old - m_new|.

    Cd and Cl are read from the results of the first forceCoeffs function
    object found in the function-object properties (native
    setResult("Cd"/"Cl")). One sample per new result: a repeated (Cd, Cl)
    pair (function object not executed this iteration, executeInterval > 1,
    or evaluated on the restored state after a rollback/skipped step) is
    not recorded again.

    Dictionary: coupled.convergence { mode forceCoeffsWindow forceCoeffsTol
    residualTol forceCoeffsRmsWindow forceCoeffsDriftTol forceCoeffsDriftAbs }.

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
public:

    //- Window statistics of one coefficient
    struct coeffStats
    {
        scalar mean = 0;
        scalar rms = 0;         //!< population standard deviation
        scalar drift = 0;       //!< |mean(older half) - mean(newer half)|
        label n = 0;            //!< samples used
    };


private:

    // Private Data

        word mode_;
        label window_;
        scalar forceTol_;
        scalar residualTol_;

        //- Window of the statistics and of the drift rule
        label rmsWindow_;

        //- Drift rule tolerances (driftTol_ 0: rule off)
        scalar driftTol_;
        scalar driftAbs_;

        DynamicList<scalar> Cd_;
        DynamicList<scalar> Cl_;
        DynamicList<scalar> Cm_;

        //- The force object provides CmPitch (Cm is monitored)
        bool haveCm_;

        //- Name of the forceCoeffs object used (empty: none found yet)
        word foName_;

        //- Requested function object (convergence.forceCoeffs; empty:
        //  auto-detect)
        word foRequested_;


    // Private Member Functions

        //- Window criterion on one history
        bool windowConverged(const DynamicList<scalar>& h) const;

        //- Statistics over the last min(rmsWindow_, n) samples of h
        coeffStats windowStats(const DynamicList<scalar>& h) const;

        //- History of monitored coefficient i (0 Cd, 1 Cl, 2 Cm)
        const DynamicList<scalar>& history(const label i) const;


public:

    // Constructors

        explicit convergenceMonitor(const dictionary& coupledDict);


    // Member Functions

        scalar residualTol() const noexcept { return residualTol_; }

        //- Append the current Cd, Cl (if a forceCoeffs object exists and
        //  produced a new result since the last sample, see the .C file)
        void record(const Time& runTime);

        //- Append one sample directly (record() uses it; unit tests).
        //  Without a force object name yet the source is named "samples".
        void addSample(const scalar Cd, const scalar Cl);
        void addSample(const scalar Cd, const scalar Cl, const scalar Cm);

        //- Restart state (subdictionary convergenceMonitor): the last
        //  max(forceCoeffsWindow, forceCoeffsRmsWindow) samples of Cd, Cl
        //  (and Cm)
        void writeState(dictionary& dict) const;
        void readState(const dictionary& dict);

        //- True if force coefficients are being monitored
        bool haveForces() const noexcept { return !foName_.empty(); }

        //- Last Cd, Cl (0 if none)
        scalar Cd() const { return Cd_.size() ? Cd_.last() : 0; }
        scalar Cl() const { return Cl_.size() ? Cl_.last() : 0; }

        //- Number of monitored coefficients (2: Cd, Cl; 3: + Cm)
        label nCoeffs() const noexcept { return haveCm_ ? 3 : 2; }

        //- Name of monitored coefficient i: Cd, Cl, Cm
        static const char* coeffName(const label i);

        //- Last value of coefficient i (0 if none)
        scalar last(const label i) const;

        //- Window statistics of coefficient i
        coeffStats stats(const label i) const;

        //- Stationary-mean (drift) rule holds for every coefficient
        //  (false if the rule is off or a window is not full)
        bool driftConverged() const;

        bool driftEnabled() const noexcept { return driftTol_ > 0; }
        label rmsWindow() const noexcept { return rmsWindow_; }

        //- Diagnostics: max over Cd and Cl of the window spread over its
        //  tolerance, (max - min)/(forceCoeffsTol*|mean|); < 1 means the
        //  window criterion holds. -1 if not available (no forces, window
        //  not full)
        scalar forceCriterionRatio() const;

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

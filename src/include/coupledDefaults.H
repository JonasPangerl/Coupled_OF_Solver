/*---------------------------------------------------------------------------*\
  coupledFoam - block-coupled p-U solver for OpenFOAM
  License: GPL-3.0-or-later
-------------------------------------------------------------------------------
Namespace
    Foam::coupledDefaults

Description
    The single place where default values of runtime dictionary entries are
    defined (spec rule 0.5, Section 11, Section 6.1 and DECISIONS.md D-008).
    Source code reads every parameter with dict.getOrDefault(key, default)
    using a constant from this header; no other numeric parameter literals
    appear in the library.

    Mathematical constants (0.5 in a mean, 2 in a factor of a formula) are
    not parameters and are not listed here.

\*---------------------------------------------------------------------------*/

#ifndef coupledFoam_coupledDefaults_H
#define coupledFoam_coupledDefaults_H

#include "label.H"
#include "scalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace coupledDefaults
{

// Linear solver, fvSolution.solvers.coupled (spec 6.1)
constexpr doubleScalar tolerance = 1e-8;
constexpr doubleScalar relTol = 0.05;
constexpr label maxIter = 200;
constexpr label minIter = 0;
constexpr label gmresRestart = 10;
constexpr label maxRestarts = 3;                 // 6.2 breakdown restarts

// blockGAMG (spec 6.1, 6.3)
constexpr label nCellsInCoarsestLevel = 500;
constexpr label mergeLevels = 1;
constexpr label nPreSweeps = 1;
constexpr label nPostSweeps = 2;
constexpr label nFinestSweeps = 2;
constexpr doubleScalar coarsestTolerance = 1e-3;
constexpr label coarsestMaxIter = 100;           // D-008
constexpr doubleScalar maxOperatorComplexity = 1.5;
constexpr label maxCopAttempts = 3;              // 6.3, D-008
constexpr bool cacheAgglomeration = true;

// Guards (spec 6.3, 9.2; D-008)
constexpr doubleScalar pivotGuard = 1e-30;
constexpr doubleScalar clampValue = 1e30;

// Outer iteration (spec 11)
constexpr label outerMaxIter = 2000;
constexpr bool potentialInit = true;
constexpr label startupUpwindIters = 50;
constexpr doubleScalar startupSwitchR = 1e-2;    // 7.4, D-008
constexpr doubleScalar nonOrthLimiter = 0.5;

// PTC (spec 7.1)
constexpr doubleScalar CFL0 = 5;
constexpr doubleScalar CFLmin = 1;
constexpr doubleScalar CFLmax = 500;
constexpr doubleScalar ptcGamma = 1;
constexpr doubleScalar betaMax = 1.5;
constexpr doubleScalar betaExp = 1.1;
constexpr label nHold = 10;

// Line search (spec 7.2)
constexpr doubleScalar fU = 0.3;
constexpr doubleScalar fp = 0.5;
constexpr doubleScalar omegaMin = 0.1;
constexpr doubleScalar kappa = 0.5;
constexpr label maxCflCuts = 3;
constexpr doubleScalar lineSearchBeta = 1.0;

// Local limit (spec 7.3)
constexpr bool localLimitEnabled = true;
constexpr doubleScalar fLoc = 0.5;

// Remediation, static (spec 8.1)
constexpr bool staticEnabled = true;
constexpr doubleScalar nonOrthThreshold = 70;
constexpr doubleScalar skewThreshold = 4;
constexpr doubleScalar volRatioThreshold = 20;
constexpr doubleScalar aspectThreshold = 1000;
constexpr doubleScalar staticBeta = 0;
constexpr doubleScalar staticNonOrthLimiter = 0.2;
constexpr doubleScalar staticCflFactor = 0.5;

// Remediation, dynamic (spec 8.2)
constexpr bool dynamicEnabled = true;
constexpr doubleScalar cU = 2;
constexpr doubleScalar cp = 5;
constexpr doubleScalar cSpike = 0.5;
constexpr label nLayers = 1;
constexpr label nQuietIters = 20;
constexpr doubleScalar dynamicCflFactor = 0.1;
constexpr bool clipToNeighbourMean = true;
constexpr doubleScalar warnFraction = 0.01;      // D-008
constexpr label warnInterval = 100;              // D-008

// Sentinel (spec 9.3; D-008)
constexpr label maxRollbacks = 5;
constexpr doubleScalar sentinelUFactor = 10;
constexpr doubleScalar sentinelPFactor = 50;
constexpr doubleScalar sentinelCflFactor = 0.25;
constexpr label sentinelMaxReport = 20;

// Bounds (spec 5.8, 9.2)
constexpr doubleScalar kMin = 1e-12;
constexpr doubleScalar boundOmegaMin = 1e-6;
constexpr doubleScalar nutMaxFactor = 1e5;

// Convergence (spec 12.3; D-005)
constexpr label forceCoeffsWindow = 100;
constexpr doubleScalar forceCoeffsTol = 0.002;
constexpr doubleScalar residualTol = 1e-6;

// Misc (spec 9.1, 10)
constexpr bool ftz = true;
constexpr bool writeState = true;

} // End namespace coupledDefaults
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

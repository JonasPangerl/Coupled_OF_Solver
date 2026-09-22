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
#include <cstdint>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace coupledDefaults
{

// Linear solver, fvSolution.solvers.coupled (spec 6.1)
constexpr doubleScalar tolerance = 1e-8;
constexpr doubleScalar relTol = 0.05;
constexpr label maxIter = 200;
constexpr label minIter = 1;                     // B1: never a zero increment
constexpr label restart = 10;                    // B1: FGMRES/GMRES restart
constexpr label restartLarge = 6;                // B7: above restartLargeCells
constexpr label restartLargeCells = 40000000;

// Eisenstat-Walker adaptive inner tolerance (amendment B2)
constexpr bool adaptiveRelTol = true;
constexpr doubleScalar etaMin = 1e-3;
constexpr doubleScalar etaMax = 0.5;
constexpr doubleScalar gammaEW = 0.9;
constexpr doubleScalar alphaEW = 2;
constexpr doubleScalar etaSafeguard = 0.1;      // B2 safeguard 1 threshold
constexpr label maxRestarts = 3;                 // 6.2 breakdown restarts

// blockGAMG (spec 6.1, 6.3)
constexpr label nCellsInCoarsestLevel = 200;     // B3 6.3.1
constexpr label mergeLevels = 2;                 // B3 6.3.1
constexpr label nPreSweeps = 1;
constexpr label nPostSweeps = 2;
constexpr label nFinestSweeps = 2;
constexpr doubleScalar coarsestTolerance = 1e-3;
constexpr label coarsestMaxIter = 50;            // B3 6.3.3
constexpr label denseLUMaxCells = 64;            // B3 6.3.3

// Cycles, ratio rule, controller (amendment B3)
constexpr doubleScalar kCycleThreshold = 0.25;
constexpr label kCycleMaxSteps = 2;
constexpr doubleScalar minCoarseningRatio = 3.0; // 6.3.1, W and K cycles
constexpr label maxMergeLevels = 4;              // 6.3.1
constexpr label procAgglomCellsPerRank = 5000;   // 6.3.4
constexpr label procAgglomDivisor = 4;           // 6.3.4: nRanks/4
constexpr bool autoTune = true;                  // 6.3.5
constexpr label tuneInterval = 50;
constexpr label nPostSweepsMax = 4;
constexpr doubleScalar tuneRhoHigh = 0.7;
constexpr doubleScalar tuneRhoLow = 0.3;
constexpr doubleScalar tuneRhoDemote = 0.2;
constexpr doubleScalar tuneRhoFail = 0.9;
constexpr label tuneConsecutiveWindows = 2;      // 6.3.5 hysteresis/failure
constexpr label minPostSweeps = 1;               // 6.3.5 lower bound
constexpr doubleScalar maxOperatorComplexity = 1.5;
// precond-research: coarse-correction scaling (none | finest | all),
// agglomeration weights (geometric | momentum | pressure | combined),
// re-agglomeration interval of matrix weights in updates (0: once)
constexpr const char* scaleCorrection = "none";
constexpr const char* agglomerationWeights = "geometric";
constexpr label reagglomerateInterval = 0;
// precond-research: blockSimple Schur approximation S = schurScale*C
constexpr doubleScalar schurScale = 2;
constexpr label maxCopAttempts = 3;              // 6.3, D-008
constexpr bool cacheAgglomeration = true;

// Guards (spec 6.3, 9.2; D-008)
constexpr doubleScalar pivotGuard = 1e-30;
constexpr doubleScalar clampValue = 1e30;        // must be <= FLT_MAX

// Strong pressure reference (D-021): the prescribed boundary fluxes of a
// closed domain must balance to this relative tolerance, |sum U_b.S_f| <=
// tol * sum |U_b.S_f|, otherwise a one-time Warning is printed. Same
// relative tolerance as the native adjustPhi continuity check.
constexpr doubleScalar refFluxBalanceTol = 1e-8;

// Outer iteration (spec 11)
constexpr label outerMaxIter = 2000;
constexpr bool potentialInit = true;
constexpr label startupUpwindIters = 50;
constexpr doubleScalar startupSwitchR = 1e-2;    // 7.4, D-008
constexpr doubleScalar nonOrthLimiter = 0.5;

// Under-relaxation of the block smoother update per sweep (OPUS_TASKS
// TASK 1d, D-043): x += smootherRelaxation * (smoothed - x). 1 = the
// undamped smoother, bit-identical to the implementation before D-043.
constexpr doubleScalar smootherRelaxation = 1.0;

// Orthogonal-mesh shortcut of nonOrthCorrection: the mesh counts as
// orthogonal if max |k_f| <= this value, k_f the native (dimensionless)
// non-orthogonal correction vector, |k_f| = tan(non-orthogonality angle).
// Far above the double round-off of an orthogonal mesh (7.1e-14 measured on
// the T0 blockMesh cavity, on which the former exact test kMax == 0 failed)
// and far below any correction the float block matrix could represent
// (eps_float ~1e-7).
constexpr doubleScalar orthogonalityTolerance = 1e-10;

// PTC (spec 7.1)
constexpr doubleScalar CFL0 = 5;
constexpr doubleScalar CFLmin = 1;
constexpr doubleScalar CFLmax = 500;
constexpr doubleScalar ptcGamma = 1;
constexpr doubleScalar betaMax = 1.5;
constexpr doubleScalar betaExp = 1.1;
constexpr label nHold = 10;
constexpr label maxLinFails = 3;                 // B4

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
// Local limiter memory (D-049): per-cell dt factor kept between iterations.
// memory no = the memoryless limiter of 7.3 (bit-identical). localRecovery
// is the per-iteration release factor of a throttled cell (the global PTC
// growth betaMax), localHold the iterations without release after a limit
// event, localStickyAfter the number of limit events after which a cell is
// never released (0 = never sticky).
constexpr bool localMemory = false;
// Implicit (Jacobi) estimate |r_P|/(V_P/dt_P + a_P) of the local limiter
// instead of the explicit |r_P| dt_P/V_P of 7.3 (D-049)
constexpr bool localImplicit = false;
constexpr doubleScalar localRecovery = 1.5;
constexpr label localHold = 0;
constexpr label localStickyAfter = 0;

// Selective Frequency Damping (amendment C3, spec 7.6; D-052)
constexpr bool sfdEnabled = false;
constexpr doubleScalar sfdChi = 0.5;              // chi* = chi Uref/Lref
constexpr doubleScalar sfdDelta = 1.0;            // Delta* = Delta Lref/Uref
constexpr doubleScalar sfdLref = 1.0;
constexpr doubleScalar sfdDeactivateBelowR = 1e-4; // SFD-off after nHold its
// Ubar = U whenever the Anderson history is flushed (C3); no = keep the
// filter state (variant)
constexpr bool sfdResetOnFlush = true;

// Remediation, static (spec 8.1)
constexpr bool staticEnabled = true;
constexpr doubleScalar nonOrthThreshold = 85;     // user, D-047 (was 70)
constexpr doubleScalar skewThreshold = 6;         // user, D-047 (was 4)
constexpr doubleScalar volRatioThreshold = 30;    // user, D-047 (was 20)
constexpr doubleScalar aspectThreshold = 2000;    // user, D-047 (was 1000)
constexpr doubleScalar staticBeta = 0;
constexpr doubleScalar staticNonOrthLimiter = 0.2;
constexpr doubleScalar staticCflFactor = 0.5;

// Remediation, dynamic (spec 8.2)
constexpr bool dynamicEnabled = true;
constexpr doubleScalar cU = 4;    // local U 2-3x U_inf under wings (D-044)
constexpr doubleScalar cp = 15;   // = cU^2-1: Cp down to -15 (D-044)
constexpr doubleScalar cSpike = 0.5;
constexpr label nLayers = 1;
constexpr label nQuietIters = 20;
// A cell that entered the dynamic set stickyAfter times stays in it for the
// rest of the run (D-049); 0 = never (the hysteresis of 8.2 only)
constexpr label dynamicStickyAfter = 3;
// Release ramp (D-049): a released cell returns to the global beta and to
// CFL factor 1 over releaseIters iterations (beta linear, CFL factor
// geometric from dynamic.cflFactor); 0 = instant release of 8.2
constexpr label dynamicReleaseIters = 0;
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
// nut/nu safety cap. This is a divergence guard, NOT a model limiter: it
// must never bind on a physically correct RANS solution. High-Re external
// aerodynamics reaches large ratios routinely - the T3 airFoil2D reference
// (simpleFoam, kOmegaSST, Re_c ~ 9e7) peaks at nut/nu = 3.6e5 in the wake,
// so the previous value of 1e5 clipped the wake eddy viscosity to 28 % of
// its correct value and drove the outer loop into a limit cycle with
// oscillating, sign-flipped force coefficients (D-040).
constexpr doubleScalar nutMaxFactor = 1e8;

// Convergence (spec 12.3; D-005)
constexpr label forceCoeffsWindow = 100;
constexpr doubleScalar forceCoeffsTol = 0.002;
constexpr doubleScalar residualTol = 1e-6;
// Stationary-mean (drift) rule of the force coefficients (D-045 b, the rule
// of D-042): off by default (tolerance 0); absolute floor of the allowed
// half-window mean drift. The statistics window forceCoeffsRmsWindow
// defaults to forceCoeffsWindow.
constexpr doubleScalar forceCoeffsDriftTol = 0;
constexpr doubleScalar forceCoeffsDriftAbs = 0.005;

// Anderson acceleration (amendment B5)
constexpr bool andersonEnabled = false;
constexpr label andersonM = 4;
constexpr doubleScalar andersonBeta = 1.0;
constexpr doubleScalar andersonMaxAlpha = 10;
constexpr label andersonMaxCells = 35000000;     // B7: not enabled above
// Numerical-rank tolerance of the Anderson QR update: sqrt(eps_double) =
// 2^-26, i.e. an orthogonalised column must keep at least half of the
// significant digits of the original column
constexpr doubleScalar andersonRankTol = 1.490116119384765625e-08;

// Zonal factors (amendment B6)
constexpr bool zonalEnabled = false;
constexpr doubleScalar zonalCflFactor = 1;       // entry default: neutral
constexpr doubleScalar zonalBeta = 1;            // entry default: neutral
constexpr label zonalNLayers = 1;                // patchDistance entry default
                                                 // (layer 1 = patch cells)

// Deep diagnostics, coupled.diagnostics (OPUS_TASKS TASK 5, D-045)
constexpr label diagLevel = 0;                   // 0 = off
constexpr label diagMaxLevel = 3;
constexpr bool diagEcho = false;
constexpr std::int64_t diagMaxBytes = 2147483648LL;  // 2 GiB per rank
// Up-leg GAMG level norms at level 2: "first" preconditioner application of
// each linear solve only (they cost two extra residuals per level visit),
// or "all"; level 3 always logs all (D-045)
constexpr const char* diagUpLeg = "first";
// Phase classification (TASK 5.2)
constexpr label diagStallWindow = 50;            // "stalled": no new min R
constexpr doubleScalar diagAsymptoticResidualFactor = 100;  // R <= 100*tol
constexpr doubleScalar diagAsymptoticForceFactor = 10;      // window <= 10*tol
// Level 3: number of locally CFL-limited cells (largest dU/(fLoc Uref))
// whose centres are logged per iteration (D-049)
constexpr label diagTopLimited = 20;

// Misc (spec 9.1, 10)
constexpr bool ftz = true;
constexpr bool writeState = true;

// Hot-loop kernels (spec 6.5; compile-time, not dictionary entries)
// Entries per chunk of blockKernels::axpyMulti(Double): 512 blockScalar/
// reduceScalar entries (2 KiB / 4 KiB) keep the target chunk in L1 while the
// m source vectors stream through
constexpr label kernelChunk = 512;
// Test-kernelBandwidth gates (6.5.6): fraction of the STREAM triad bandwidth
constexpr doubleScalar gateAmulFraction = 0.6;
constexpr doubleScalar gateAxpyDotFraction = 0.8;

} // End namespace coupledDefaults
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

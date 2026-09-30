# coupledFoam keyword reference

Every run-time setting of coupledFoam, with its default, admissible range,
meaning and the decision or spec section behind it. All defaults live in
`src/include/coupledDefaults.H`; the column *constant* names the constant
there. `tests/test_keywords.py` fails if a constant is added to that header
without an entry in this file.

Where the keywords live (all in `system/fvSolution`):

- `solvers.coupled { ... }`: the block linear solver (abbreviated `sc.`).
- `solvers.coupled.blockGAMG { ... }`: the block-GAMG preconditioner, its
  smoothers, the coarsest-level solver and the autoTune controller
  (abbreviated `gamg.`).
- `coupled { ... }`: the outer iteration (abbreviated `c.`).

The solver prints every effective value at start-up under
`coupledFoam: effective settings`. The fvSchemes entries the solver reads
(`div(phi,U)`, `div(phi,U)_upwind`, `grad(U)`, `grad(p)`,
`interpolate(nuEff)`, ...) are the native ones and are not listed here.

Nothing a user may want to change needs a recompile (D-066). The few
compile-time items are listed at the end with the reason.

## Remediation cell categories (D-061, D-066): user guide

Before the first iteration, coupledFoam **pre-selects** cells that get a
special, more careful treatment for the whole run (the *static set*, spec
8.1). Since D-066 these cells fall into four categories. Each category can
be switched off and tuned separately in `coupled.remediation`:

| category | which cells | default treatment |
|---|---|---|
| `meshQuality` | cells with poor mesh quality: a face that is too non-orthogonal or too skewed, a large volume jump to a neighbour, or a very high aspect ratio (spec 8.1) | full: first-order convection (`beta 0`), half the local time step, limited pressure gradient and limited non-orthogonal correction |
| `badMesh` | severe defects: a volume ratio above 50 to a neighbour (C1 volumeJump), open cells (faces that do not close the cell), optionally faces beyond severe non-orthogonality/skewness limits | full, as `meshQuality` |
| `processor` | cells at processor boundaries: cells with both a processor and a cyclicAMI face (C1 procAMI), optionally `nLayers` cell layers next to every processor boundary | mild: full convection scheme (`beta 1`), no limiters, half the local time step |
| `wall` | cells at walls: "starved" wall cells with at most 2 neighbour cells (C1 wallStarved, **off by default**), optionally `nLayers` cell layers next to chosen patches | mild, as `processor` |

Why two kinds of treatment: `meshQuality` and `badMesh` cells actually are
defective, and first-order convection with limiters keeps them stable. The
`processor` and `wall` cells are chosen by *position*, not because they
misbehave. Giving them first-order convection and limiters changes the
solution where it matters most. On the coarse motorBike (T4a), 12
wallStarved cells on the body surface under the full treatment moved the
lift coefficient by -18 % (D-061). Such cells therefore only get a shorter
local time step by default.

Rules:

- A cell in several categories gets the most restrictive setting of each
  parameter: the smallest `beta`, the smallest `cflFactor`, and a limiter if
  any of its categories asks for one.
- Cells that misbehave during the run (the dynamic set, spec 8.2) always
  get the full dynamic treatment (`beta 0`, `dynamic.cflFactor`),
  whatever their category.
- Zonal factors (`zonal`, spec 8.3) multiply on top of everything.

Output:

- Log: one line per category at start-up with its cell count and the
  count of each criterion.
- Cell sets: `remediationStatic` (all categories), `remediationDynamic`,
  `remediationMeshQuality`, `remediationBadMesh`, `remediationProcessor`,
  `remediationWall`.
- The field `remediationFlag`, whose bits are OR-ed:

| bit | meaning |
|---|---|
| 1 | static: in any category |
| 2 | dynamic set |
| 4 | criterion wallStarved (C1) |
| 8 | criterion procAMI (C1) |
| 16 | criterion volumeJump (C1) |
| 32 | category meshQuality |
| 64 | category badMesh |
| 128 | category processor |
| 256 | category wall |
| 512 | criterion closedness (open cell) |
| 1024 | criterion severe non-orthogonality/skewness (badMesh thresholds) |
| 2048 | criterion processor-patch layer |
| 4096 | criterion wall-patch layer |

- Counts: `summary.json` has `staticCells` (union) and `staticCategories`
  (per category). At diagnostics level >= 1, `controls.remediation.nStat`
  and `controls.remediation.nStatCat` are written. `nStat=` on the `CF|`
  line is the union.

Example: put the wall-adjacent cells of the body into the mild wall
category and give them a quarter of the local time step:

    remediation
    {
        wall { patches (motorBikeGroup); nLayers 1; cflFactor 0.25; }
    }

Old layout. The pre-D-066 entries `remediation.static { ... }` still work
and print one deprecation note. They map as follows:

- `enabled` switches all four categories.
- `nonOrthThreshold`, `skewThreshold`, `volRatioThreshold`,
  `aspectThreshold` go to `meshQuality`.
- `beta` and `cflFactor` go to `meshQuality` and `badMesh`.
- `nonOrthLimiter` becomes `limitedNonOrthCoeff`.
- `volumeJump` and `volJumpThreshold` go to `badMesh`.
- `procAMI` goes to `processor`.
- `wallStarved` and `maxWallInternalFaces` go to `wall`, which now means the
  mild treatment.

The same setting given in both the old and the new place with different
values is an error. The pre-release `remediation.mild` dictionary of D-061
is rejected, because it has been replaced by the categories.

## c. - outer iteration (`coupled { }`)

| keyword | default | range | meaning | ref | constant |
|---|---|---|---|---|---|
| `c.maxIter` | 2000 | >= 1 | outer iterations | 11 | `outerMaxIter` |
| `c.potentialInit` | yes | bool | potentialFoam start fields expected | 7.4 | `potentialInit` |
| `c.potentialClip` | 0 | 0 or >= 1 | clip \|U\| of the potential start to f Uref (0 off); U only since D-077 | D-057, D-077 | `potentialClip` |
| `c.potentialPressure` | keep | keep, bernoulli | p of the potential start: as read, or H - \|U\|^2/2 in every cell after the clip (H from the fixed-p boundary) | D-077 | `potentialPressure` |
| `c.nonOrthLimiter` | 0.5 | [0, 1] | limitedSnGrad coefficient of the non-orthogonal correction | 5.3, D-014 | `nonOrthLimiter` |
| `c.orthogonalityTolerance` | 0 | >= 0 | mesh counts as orthogonal if max \|k_f\| <= value (0: exact test) | FABLE item 1, D-066 | `orthogonalityTolerance` |
| `c.ftz` | yes | bool | flush-to-zero/denormals-are-zero | 9.1 | `ftz` |
| `c.writeState` | yes | bool | write the restart state | 10 | `writeState` |
| `c.dumpLinearSystem` | () | list of iterations | dump the linear system (serial, Test-blockSystem) | D-039 | - |
| `c.Uref` | boundary | boundary, field, value | velocity scale of the classification; boundary = largest PRESCRIBED boundary velocity (D-074) | D-050, D-074 | `UrefMode` |
| `c.UrefFallbackFactor` | 1e-6 | [0, 1] | boundary Uref below this fraction of the field maximum: field maximum | D-050, D-066 | `UrefFallbackFactor` |
| `c.UrefStep` | reference | reference, fieldCapped, field, value | velocity scale of the step limits | D-057 | `UrefStepMode` |
| `c.UrefStepCap` | 4 | > 0 | cap of fieldCapped, in Uref | D-057 | `UrefStepCap` |
| `c.startupReference` | exclude | none, ramp, exclude | references during the start-up | D-057 | `startupReference` |
| `c.startupMode` | hybrid | hybrid, upwind, none | start-up of the convection scheme | D-048 | `startupMode` |
| `c.startupUpwindIters` | 50 | >= 0 | upwind iterations (upwind mode) | 7.4 | `startupUpwindIters` |
| `c.startupSwitchR` | 1e-2 | > 0 | R/R1 that ends the upwind phase | 7.4, D-008 | `startupSwitchR` |
| `c.startupRampStart` | 10 | >= 0 | ramp start iteration | D-048 | `startupRampStart` |
| `c.startupRampLength` | 20 | >= 1 | ramp length | D-048 | `startupRampLength` |
| `c.startupRampStartMax` | 10 | >= startupRampStart | latest delayed ramp start | D-048 | `startupRampStartMax` |
| `c.startupStagnationTrigger` | no | bool | start the ramp on stagnation | D-048 | `startupStagnationTrigger` |
| `c.startupStagnationFactor` | 0.9 | (0, 1) | stagnation: R_n > f R_(n-W) | D-048 | `startupStagnationFactor` |
| `c.startupStagnationWindow` | 5 | >= 1 | W of the stagnation test | D-048 | `startupStagnationWindow` |
| `c.startupFastFactor` | 0.7 | (0, 1) | fast drop delays the ramp: R_n < f R_(n-W) | D-048 | `startupFastFactor` |
| `c.startupDevelopedTol` | 1e-2 | > 0 | developed start: keep beta 1 | D-048 | `startupDevelopedTol` |

### c.ptc, c.localLimit

| keyword | default | range | meaning | ref | constant |
|---|---|---|---|---|---|
| `c.ptc.cflStrategy` | mRDM | mRDM, EXP, SER, EFF | global CFL update rule; EFF: hold the linear solve at its optimal share of the wall time (D-083) | 7.1, D-083 | `cflStrategy` |
| `c.ptc.CFL0` | 5 | > 0 | start CFL | 7.1 | `CFL0` |
| `c.ptc.CFLmin` | 1 | > 0 | lower CFL bound | 7.1 | `CFLmin` |
| `c.ptc.CFLmax` | 500 | >= CFLmin | upper CFL bound | 7.1 | `CFLmax` |
| `c.ptc.gamma` | 1 | > 0 | exponent of the CFL rule | 7.1 | `ptcGamma` |
| `c.ptc.betaMax` | 1.5 | > 1 | max CFL growth per iteration | 7.1 | `betaMax` |
| `c.ptc.betaExp` | 1.1 | > 1 | EXP growth factor | 7.1 | `betaExp` |
| `c.ptc.nHold` | 10 | >= 0 | iterations without growth after a cut | 7.1 | `nHold` |
| `c.ptc.continuityFactor` | 1 | >= 0 | continuity pseudo-time term c V/(dtau Uref^2) on the p-p diagonal of the increment system, in the static and dynamic remediation cells only (0 off) | D-075 | `ptcContinuityFactor` |
| `c.ptc.failCeiling` | 0 | [0, 1) | after a failed linear solve the CFL may not grow above this x the failed CFL (0 off) | D-081 | `ptcFailCeiling` |
| `c.ptc.ceilingRelax` | 1.02 | >= 1 | growth of that ceiling per accepted iteration | D-081 | `ptcCeilingRelax` |
| `c.ptc.effTol` | 0.02 | [0, 1) | EFF: relative band around effSolveShare in which the CFL is kept | D-083 | `ptcEffTol` |
| `c.ptc.effSolveShare` | 0.6 | (0, 1) | EFF: target share of the linear solve in the wall time of an iteration; below it the CFL goes up, above it down (step from effExponent) | D-083 | `ptcEffSolveShare` |
| `c.ptc.effWindow` | 4 | >= 1 | EFF: iterations averaged at one CFL before a decision (the first after a change is skipped) | D-083 | `ptcEffWindow` |
| `c.ptc.effExponent` | 1.6 | > 0 | EFF: exponent b of the solve cost ~ CFL^b, sets the step towards the target share (clamped x0.5..x2) | D-083 | `ptcEffExponent` |
| `c.ptc.maxLinFails` | 3 | >= 1 | consecutive linear failures before the abort | B4 | `maxLinFails` |
| `c.ptc.linFailPolicy` | strict | strict, reduction | what counts as a linear failure | B4 | `linFailPolicy` |
| `c.ptc.linAcceptReduction` | 0.9 | (0, 1) | accepted reduction of a capped solve | B4 | `linAcceptReduction` |
| `c.localLimit.enabled` | yes | bool | solution-limited local CFL | 7.3 | `localLimitEnabled` |
| `c.localLimit.fLoc` | 0.5 | > 0 | allowed dU in Ustep | 7.3 | `fLoc` |
| `c.localLimit.implicit` | no | bool | implicit (Jacobi) estimate | D-055 | `localImplicit` |
| `c.localLimit.memory` | no | bool | per-cell dt factor memory | D-055 | `localMemory` |
| `c.localLimit.localRecovery` | 1.5 | > 1 | release factor per iteration | D-055 | `localRecovery` |
| `c.localLimit.localHold` | 0 | >= 0 | iterations without release | D-055 | `localHold` |
| `c.localLimit.localStickyAfter` | 0 | >= 0 | limit events until never released (0 never) | D-055 | `localStickyAfter` |

### c.lineSearch

| keyword | default | range | meaning | ref | constant |
|---|---|---|---|---|---|
| `c.lineSearch.fU` | 0.3 | > 0 | allowed dU per step in Ustep | 7.2 | `fU` |
| `c.lineSearch.fp` | 0.5 | > 0 | allowed dp per step in pstep | 7.2 | `fp` |
| `c.lineSearch.omegaMin` | 0.1 | (0, 1] | smallest step fraction | 7.2 | `omegaMin` |
| `c.lineSearch.kappa` | 0.5 | (0, 1) | CFL cut factor | 7.2 | `kappa` |
| `c.lineSearch.maxCflCuts` | 3 | >= 0 | CFL cuts per iteration | 7.2 | `maxCflCuts` |
| `c.lineSearch.beta` | 1.0 | (0, 1] | line-search blending | D-036 | `lineSearchBeta` |
| `c.lineSearch.mode` | global | global, local | global: one omega, the smallest over all cells; local: omega is the localFraction quantile, the few cells over the bounds go into the dynamic set | D-079 | `lineSearchMode` |
| `c.lineSearch.localFraction` | 0.0005 | [0, 1] | local mode: omega is the quantile at which at most this fraction of the cells violate the bounds; they go into the dynamic set | D-079 | `lineSearchLocalFraction` |

### c.remediation (spec 8, D-066)

| keyword | default | range | meaning | ref | constant |
|---|---|---|---|---|---|
| `c.remediation.meshQuality.enabled` | yes | bool | category switch | 8.1, D-066 | `meshQualityEnabled` |
| `c.remediation.meshQuality.nonOrthThreshold` | 85 | deg | face non-orthogonality | D-047 | `nonOrthThreshold` |
| `c.remediation.meshQuality.skewThreshold` | 6 | > 0 | face skewness | D-047 | `skewThreshold` |
| `c.remediation.meshQuality.volRatioThreshold` | 30 | > 1 | neighbour volume ratio | D-047 | `volRatioThreshold` |
| `c.remediation.meshQuality.aspectThreshold` | 2000 | > 1 | cell aspect ratio | D-047 | `aspectThreshold` |
| `c.remediation.meshQuality.beta` | 0 | [0, 1] | upper bound of the convection blending | 8.1 | `meshQualityBeta` |
| `c.remediation.meshQuality.cflFactor` | 0.5 | (0, 1] | local time step factor | 8.1 | `meshQualityCflFactor` |
| `c.remediation.meshQuality.gradLimiter` | yes | bool | limitedGradScheme for grad p | D-018 | `meshQualityGradLimiter` |
| `c.remediation.meshQuality.nonOrthLimiter` | yes | bool | limitedNonOrthCoeff on the cell faces | 8.1 | `meshQualityNonOrthLimiter` |
| `c.remediation.badMesh.enabled` | yes | bool | category switch | D-066 | `badMeshEnabled` |
| `c.remediation.badMesh.volumeJump` | yes | bool | C1 volumeJump criterion | C1 | `volumeJumpEnabled` |
| `c.remediation.badMesh.volJumpThreshold` | 0.98 | (0, 1) | \|V_P - V_N\|/max(V_P, V_N) | C1 | `volJumpThreshold` |
| `c.remediation.badMesh.closednessThreshold` | 1e-6 | >= 0 (0 off) | open cell: \|sum S_f\|/sum \|S_f\| | D-066 | `closednessThreshold` |
| `c.remediation.badMesh.nonOrthThreshold` | 0 | [0, 180] deg (0 off) | severe face non-orthogonality | D-066 | `severeNonOrthThreshold` |
| `c.remediation.badMesh.skewThreshold` | 0 | >= 0 (0 off) | severe face skewness | D-066 | `severeSkewThreshold` |
| `c.remediation.badMesh.beta` | 0 | [0, 1] | as meshQuality | D-066 | `badMeshBeta` |
| `c.remediation.badMesh.cflFactor` | 0.5 | (0, 1] | as meshQuality | D-066 | `badMeshCflFactor` |
| `c.remediation.badMesh.gradLimiter` | yes | bool | as meshQuality | D-066 | `badMeshGradLimiter` |
| `c.remediation.badMesh.nonOrthLimiter` | yes | bool | as meshQuality | D-066 | `badMeshNonOrthLimiter` |
| `c.remediation.processor.enabled` | yes | bool | category switch | D-066 | `processorEnabled` |
| `c.remediation.processor.procAMI` | yes | bool | C1 procAMI criterion | C1 | `procAMIEnabled` |
| `c.remediation.processor.nLayers` | 0 | >= 0 | cell layers from the processor patches | D-066 | `processorNLayers` |
| `c.remediation.processor.beta` | 1 | [0, 1] | mild | D-061 | `processorBeta` |
| `c.remediation.processor.cflFactor` | 0.5 | (0, 1] | mild | D-061 | `processorCflFactor` |
| `c.remediation.processor.gradLimiter` | no | bool | mild | D-061 | `processorGradLimiter` |
| `c.remediation.processor.nonOrthLimiter` | no | bool | mild | D-061 | `processorNonOrthLimiter` |
| `c.remediation.wall.enabled` | yes | bool | category switch | D-066 | `wallEnabled` |
| `c.remediation.wall.wallStarved` | no | bool | C1 wallStarved criterion (user: off) | C1, D-061 | `wallStarvedEnabled` |
| `c.remediation.wall.maxWallInternalFaces` | 2 | >= 0 | wallStarved: max implicit-neighbour faces | C1 | `maxWallInternalFaces` |
| `c.remediation.wall.patches` | () | wordRes | patches (names, regex, groups) seeding layers | D-066 | - |
| `c.remediation.wall.nLayers` | 1 | >= 1 | cell layers from wall.patches | D-066 | `wallNLayers` |
| `c.remediation.wall.beta` | 1 | [0, 1] | mild | D-061 | `wallBeta` |
| `c.remediation.wall.cflFactor` | 0.5 | (0, 1] | mild | D-061 | `wallCflFactor` |
| `c.remediation.wall.gradLimiter` | no | bool | mild | D-061 | `wallGradLimiter` |
| `c.remediation.wall.nonOrthLimiter` | no | bool | mild | D-061 | `wallNonOrthLimiter` |
| `c.remediation.limitedNonOrthCoeff` | 0.2 | [0, 1] | non-orthogonal limiter of nonOrthLimiter cells (was static.nonOrthLimiter) | 8.1, D-066 | `limitedNonOrthCoeff` |
| `c.remediation.limitedGradScheme` | "cellLimited Gauss linear 1" | gradScheme | grad p of gradLimiter cells | D-018, D-066 | `limitedGradScheme` |
| `c.remediation.dynamic.enabled` | yes | bool | dynamic set | 8.2 | `dynamicEnabled` |
| `c.remediation.dynamic.cU` | 4 | > 0 | \|U\| > cU Uref | D-044 | `cU` |
| `c.remediation.dynamic.cp` | 15 | > 0 | \|p\| > cp pref | D-044 | `cp` |
| `c.remediation.dynamic.cSpike` | 0.5 | > 0 | \|U - Ubar_N\| > cSpike Uref | 8.2 | `cSpike` |
| `c.remediation.dynamic.nLayers` | 1 | >= 0 | neighbour layers added | 8.2 | `nLayers` |
| `c.remediation.dynamic.nQuietIters` | 20 | >= 1 | iterations until release | 8.2 | `nQuietIters` |
| `c.remediation.dynamic.stickyAfter` | 0 | >= 0 (0 never) | entries until never released | D-055 | `dynamicStickyAfter` |
| `c.remediation.dynamic.releaseIters` | 0 | >= 0 | release ramp length | D-055 | `dynamicReleaseIters` |
| `c.remediation.dynamic.cflFactor` | 0.1 | (0, 1] | local time step factor | 8.2 | `dynamicCflFactor` |
| `c.remediation.dynamic.clipToNeighbourMean` | yes | bool | increment clipping | 8.2 | `clipToNeighbourMean` |
| `c.remediation.warnFraction` | 0.01 | (0, 1] | warn above this fraction of cells | D-008 | `warnFraction` |
| `c.remediation.warnInterval` | 100 | >= 0 | warning interval (0 off) | D-008 | `warnInterval` |
| `c.remediation.zonal.enabled` | no | bool | zonal factors | 8.3, B6 | `zonalEnabled` |
| `c.remediation.zonal.zones[i].cellZone` | - | wordRe or built-in | cellZone or `_wallCells`, `_procCells`, `_amiCells`, `_procAMICells`, `_wallStarved`, `_remediationStatic`, `_remediationMeshQuality`, `_remediationBadMesh`, `_remediationProcessor`, `_remediationWall` | C5, D-066 | - |
| `c.remediation.zonal.*[i].cflFactor` | 1 | > 0 | CFL multiplier | B6 | `zonalCflFactor` |
| `c.remediation.zonal.*[i].beta` | 1 | [0, 1] | beta multiplier | B6 | `zonalBeta` |
| `c.remediation.zonal.patchDistance[i].patches` | - | wordRes | seed patches | B6, C5 | - |
| `c.remediation.zonal.patchDistance[i].nLayers` | 1 | >= 1 | cell layers | B6 | `zonalNLayers` |

### c.sentinel, c.bounds, c.guards, c.rhieChow

| keyword | default | range | meaning | ref | constant |
|---|---|---|---|---|---|
| `c.sentinel.maxRollbacks` | 5 | >= 0 | rollbacks before the abort | 9.3 | `maxRollbacks` |
| `c.sentinel.UFactor` | 10 | > 0 | \|U\| > f Uref is non-physical | 9.3 | `sentinelUFactor` |
| `c.sentinel.pFactor` | 50 | > 0 | \|p\| > f pref is non-physical | 9.3 | `sentinelPFactor` |
| `c.sentinel.cflFactor` | 0.25 | (0, 1) | CFL factor after a rollback | 9.3 | `sentinelCflFactor` |
| `c.sentinel.maxReport` | 20 | >= 0 | offending cells printed | 9.3 | `sentinelMaxReport` |
| `c.bounds.kMin` | 1e-12 | > 0 | k lower bound | 5.8 | `kMin` |
| `c.bounds.omegaMin` | 1e-6 | > 0 | omega lower bound | 5.8 | `boundOmegaMin` |
| `c.bounds.nutMaxFactor` | 1e8 | > 0 | nut <= f nu (divergence guard) | D-040 | `nutMaxFactor` |
| `c.guards.clampValue` | 1e30 | (0, FLT_MAX] | clamp of block coefficients | 9.2 | `clampValue` |
| `c.guards.refFluxBalanceTol` | 1e-8 | >= 0 | closed-domain flux balance warning | D-021, D-066 | `refFluxBalanceTol` |
| `c.rhieChow.tensorial` | yes | bool | D = V A^-1 | C2, D-053 | `rhieChowTensorial` |
| `c.rhieChow.detRelTol` | 1e-12 | >= 0 | singular-block test | C2 | `rhieChowDetRelTol` |
| `c.rhieChow.pinvRelTol` | 1e-6 | (0, 1) | pseudo-inverse cut-off | C2 | `rhieChowPinvRelTol` |
| `c.rhieChow.pinvMaxSweeps` | 60 | >= 1 | Jacobi SVD sweeps | C2, D-066 | `rhieChowPinvMaxSweeps` |
| `c.rhieChow.warnInterval` | 100 | >= 1 | pseudo-inverse warning interval | C2, D-066 | `rhieChowWarnInterval` |

### c.sfd, c.anderson, c.convergence

| keyword | default | range | meaning | ref | constant |
|---|---|---|---|---|---|
| `c.sfd.enabled` | no | bool | selective frequency damping | C3, D-052 | `sfdEnabled` |
| `c.sfd.chi` | 0.5 | >= 0 | chi* = chi Uref/Lref | C3 | `sfdChi` |
| `c.sfd.Delta` | 1.0 | > 0 | Delta* = Delta Lref/Uref | C3 | `sfdDelta` |
| `c.sfd.Lref` | 1.0 | > 0 | reference length | C3 | `sfdLref` |
| `c.sfd.deactivateBelowR` | 0 | >= 0 (0 never) | SFD-off threshold | D-058 | `sfdDeactivateBelowR` |
| `c.sfd.nHold` | c.ptc.nHold | >= 1 | iterations below deactivateBelowR | C3, D-066 | `nHold` |
| `c.sfd.resetOnFlush` | yes | bool | Ubar = U on an Anderson flush | C3 | `sfdResetOnFlush` |
| `c.sfd.afterStartup` | yes | bool | activate after the start-up | D-052 | `sfdAfterStartup` |
| `c.sfd.startIter` | 0 | >= 0 | first iteration | D-052 | `sfdStartIter` |
| `c.anderson.enabled` | no | bool | Anderson acceleration | B5 | `andersonEnabled` |
| `c.anderson.m` | 4 | >= 1 | history depth | B5 | `andersonM` |
| `c.anderson.beta` | 1.0 | > 0 | mixing | B5 | `andersonBeta` |
| `c.anderson.maxAlpha` | 10 | > 0 | coefficient bound | B5 | `andersonMaxAlpha` |
| `c.anderson.maxCells` | 35000000 | >= 1 | not enabled above (memory) | B7, D-066 | `andersonMaxCells` |
| `c.anderson.rankTol` | 2^-26 | (0, 1) | numerical-rank tolerance of the QR update | B5, D-066 | `andersonRankTol` |
| `c.convergence.mode` | any | any, all | residual or force criterion | 12.3 | `convergenceMode` |
| `c.convergence.forceCoeffs` | "" | function object name | force coefficients used ("": first with Cd and Cl) | D-066 | - |
| `c.convergence.forceCoeffsWindow` | 100 | >= 1 | force window | D-005 | `forceCoeffsWindow` |
| `c.convergence.forceCoeffsTol` | 0.002 | > 0 | force window tolerance | D-005 | `forceCoeffsTol` |
| `c.convergence.residualTol` | 1e-6 | > 0 | R tolerance | D-005 | `residualTol` |
| `c.convergence.forceCoeffsRmsWindow` | forceCoeffsWindow | >= 1 | statistics window | D-045 | - |
| `c.convergence.forceCoeffsDriftTol` | 0 | >= 0 (0 off) | stationary-mean rule | D-042, D-045 | `forceCoeffsDriftTol` |
| `c.convergence.forceCoeffsDriftAbs` | 0.005 | > 0 | absolute drift floor | D-045 | `forceCoeffsDriftAbs` |

### c.diagnostics, c.diagnosticFields

| keyword | default | range | meaning | ref | constant |
|---|---|---|---|---|---|
| `c.diagnostics.level` | 0 | 0..3 | deep diagnostics level | D-045 | `diagLevel`, bound `diagMaxLevel` |
| `c.diagnostics.echo` | no | bool | echo to the log | D-045 | `diagEcho` |
| `c.diagnostics.maxBytes` | 2 GiB | > 0 | per-rank file cap | D-045 | `diagMaxBytes` |
| `c.diagnostics.upLeg` | first | first, all | level-2 up-leg GAMG norms | D-045 | `diagUpLeg` |
| `c.diagnostics.stallWindow` | 50 | >= 1 | "stalled" phase window | D-045, D-066 | `diagStallWindow` |
| `c.diagnostics.asymptoticResidualFactor` | 100 | > 0 | "asymptotic": R <= f residualTol | D-045, D-066 | `diagAsymptoticResidualFactor` |
| `c.diagnostics.asymptoticForceFactor` | 10 | > 0 | "asymptotic": window <= f tolerance | D-045, D-066 | `diagAsymptoticForceFactor` |
| `c.diagnostics.topLimited` | 20 | >= 0 | level-3 cells logged | D-055, D-066 | `diagTopLimited` |
| `c.diagnosticFields.enabled` | yes | bool | localDt, localCFL, ... at write times | C6 | `diagnosticFieldsEnabled` |

## sc. - block linear solver (`solvers.coupled { }`)

| keyword | default | range | meaning | ref | constant |
|---|---|---|---|---|---|
| `sc.solver` | (required) | blockFGMRES, blockGMRES, blockBiCGStab | Krylov solver | 6.1 | - |
| `sc.preconditioner` | none | none, blockGAMG, blockSmoother, blockDiagonal, blockSimple | preconditioner | 6.1 | `preconditioner` |
| `sc.tolerance` | 1e-8 | >= 0 | absolute tolerance | 6.1 | `tolerance` |
| `sc.relTol` | 0.05 | [0, 1) | relative tolerance (fixed eta) | 6.1 | `relTol` |
| `sc.maxIter` | 200 | >= 1 | Krylov iterations | 6.1 | `maxIter` |
| `sc.minIter` | 1 | >= 0 | minimum iterations | B1 | `minIter` |
| `sc.maxRestarts` | 3 | >= 0 | BiCGStab breakdown restarts | 6.2 | `maxRestarts` |
| `sc.stagnationRatio` | 0 | [0, 1) | FGMRES: stop unconverged when a restart cycle leaves the true residual above this x the previous cycle's (0 off) | D-080 | `stagnationRatio` |
| `sc.restart` (alias `gmresRestart`) | 10 | >= 1 | (F)GMRES restart | B1 | `restart` |
| `sc.restartLarge` | 6 | >= 1 | restart above restartLargeCells unless restart is given | B7, D-066 | `restartLarge` |
| `sc.restartLargeCells` | 40000000 | >= 1 | cell count of the B7 rule | B7, D-066 | `restartLargeCells` |
| `sc.adaptiveRelTol` | yes | bool | Eisenstat-Walker eta | B2 | `adaptiveRelTol` |
| `sc.etaMin` | 1e-3 | (0, etaMax] | eta bounds | B2 | `etaMin` |
| `sc.etaMax` | 0.5 | [etaMin, 1) | eta bounds | B2 | `etaMax` |
| `sc.gammaEW` | 0.9 | (0, 1] | EW gamma | B2 | `gammaEW` |
| `sc.alphaEW` | 2 | (1, 2] | EW alpha | B2 | `alphaEW` |
| `sc.etaSafeguard` | 0.1 | > 0 | EW safeguard 1 threshold | B2 | `etaSafeguard` |
| `sc.pivotGuard` | 1e-30 | > 0 | pivot guard of blockDiagonal/blockSmoother/blockSimple | 6.3, D-008 | `pivotGuard` |
| `sc.smoother` | blockILU0 | smoother | smoother of preconditioner blockSmoother | 6.1 | `smootherPreconSmoother` |
| `sc.nSweeps` | 1 | >= 1 | sweeps of preconditioner blockSmoother | 6.1, D-066 | `smootherPreconSweeps` |
| `sc.smootherRelaxation` | 1.0 | (0, 1] | smoother damping (blockSmoother) | D-043 | `smootherRelaxation` |
| `sc.pivotGrowthLimit` | 20 | >= 0 (0 off) | ILU0 pivot-growth guard (blockSmoother) | D-049 | `iluPivotGrowthLimit` |
| `sc.schurScale` | 2 | > 0 | blockSimple S = s C | precond-research | `schurScale` |
| `sc.simpleMode` | single | single, sequential | blockSimple mode | precond-research | `simpleMode` |

### gamg. (`solvers.coupled.blockGAMG { }`)

| keyword | default | range | meaning | ref | constant |
|---|---|---|---|---|---|
| `gamg.cycleType` | K | V, W, K | cycle | 6.3 | `cycleType` |
| `gamg.agglomerator` | faceAreaPair | faceAreaPair, algebraicPair | agglomeration | 6.3.1, C4 | `agglomerator` |
| `gamg.agglomerationWeights` | geometric | geometric, momentum, pressure, combined | pair weights | D-039 | `agglomerationWeights` |
| `gamg.reagglomerateInterval` | 0 | >= 0 (0 once) | re-agglomeration of matrix weights | D-039 | `reagglomerateInterval` |
| `gamg.nCellsInCoarsestLevel` | 200 | >= 1 | coarsest level size | B3 | `nCellsInCoarsestLevel` |
| `gamg.mergeLevels` | 2 | >= 1 | agglomeration levels merged | B3 | `mergeLevels` |
| `gamg.maxMergeLevels` | 4 | >= mergeLevels | ratio rule retries | 6.3.1 | `maxMergeLevels` |
| `gamg.minCoarseningRatio` | 3.0 | > 1 | ratio rule (W, K) | 6.3.1 | `minCoarseningRatio` |
| `gamg.maxOperatorComplexity` | 1.5 | > 1 | operator complexity limit | 6.3 | `maxOperatorComplexity` |
| `gamg.cacheAgglomeration` | yes | bool | keep the agglomeration | 6.3 | `cacheAgglomeration` |
| `gamg.processorAgglomerator` | masterCoarsest | masterCoarsest, none, ... | processor agglomeration | 6.3.4, D-030 | `processorAgglomerator` |
| `gamg.procAgglomCellsPerRank` | 5000 | >= 1 | trigger | 6.3.4 | `procAgglomCellsPerRank` |
| `gamg.procAgglomDivisor` | 4 | >= 1 | nRanks/divisor | 6.3.4 | `procAgglomDivisor` |
| `gamg.smoother` | blockGaussSeidel | blockGaussSeidel, blockILU0, ... | level smoother | 6.1 | `smoother` |
| `gamg.nPreSweeps` | 1 | >= 0 | pre-smoothing | 6.1 | `nPreSweeps` |
| `gamg.nPostSweeps` | 2 | >= 1 | post-smoothing | 6.1 | `nPostSweeps` |
| `gamg.nFinestSweeps` | 2 | >= 1 | finest-level sweeps | 6.1 | `nFinestSweeps` |
| `gamg.smootherRelaxation` | 1.0 | (0, 1] | smoother damping | D-043 | `smootherRelaxation` |
| `gamg.pivotGrowthLimit` | 20 | >= 0 (0 off) | ILU0 pivot-growth guard | D-049 | `iluPivotGrowthLimit` |
| `gamg.pivotGuard` | 1e-30 | > 0 | pivot guard | 6.3, D-008 | `pivotGuard` |
| `gamg.kCycleThreshold` | 0.25 | (0, 1) | K-cycle second step if r1 > t r0 | 6.3.2 | `kCycleThreshold` |
| `gamg.kCycleMaxSteps` | 2 | 1, 2 | K-cycle GCR steps | 6.3.2, D-066 | `kCycleMaxSteps` |
| `gamg.scaleCorrection` | none | none, finest, all | coarse-correction scaling | precond-research | `scaleCorrection` |
| `gamg.denseLUMaxCells` | 64 | >= 0 | dense LU on the coarsest level up to this size | B3 | `denseLUMaxCells` |
| `gamg.coarsestSolver` | blockBiCGStab | Krylov solver | coarsest-level solver | 6.3.3 | `coarsestSolver` |
| `gamg.coarsestPreconditioner` | blockDiagonal | preconditioner | its preconditioner | 6.3.3, D-066 | `coarsestPreconditioner` |
| `gamg.coarsestTolerance` | 1e-3 | (0, 1) | its relTol | 6.3.3 | `coarsestTolerance` |
| `gamg.coarsestAbsTolerance` | 0 | >= 0 | its absolute tolerance | 6.3.3, D-066 | `coarsestAbsTolerance` |
| `gamg.coarsestMaxIter` | 50 | >= 1 | its maxIter | B3 | `coarsestMaxIter` |
| `gamg.coarsestMinIter` | sc default | >= 0 | its minIter (optional) | D-066 | - |
| `gamg.coarsestMaxRestarts` | sc default | >= 0 | its maxRestarts (optional) | D-066 | - |
| `gamg.autoTune` | yes | bool | cycle controller | 6.3.5 | `autoTune` |
| `gamg.tuneInterval` | 50 | >= 1 | window length | 6.3.5 | `tuneInterval` |
| `gamg.nPostSweepsMax` | 4 | >= nPostSweepsMin | upper sweep bound | 6.3.5 | `nPostSweepsMax` |
| `gamg.nPostSweepsMin` | 1 | >= 1 | lower sweep bound | 6.3.5, D-066 | `minPostSweeps` |
| `gamg.tuneRhoHigh` | 0.7 | (tuneRhoLow, tuneRhoFail] | add a sweep / promote | 6.3.5, D-066 | `tuneRhoHigh` |
| `gamg.tuneRhoLow` | 0.3 | [tuneRhoDemote, tuneRhoHigh) | remove a sweep | 6.3.5, D-066 | `tuneRhoLow` |
| `gamg.tuneRhoDemote` | 0.2 | (0, tuneRhoLow] | demote the cycle | 6.3.5, D-066 | `tuneRhoDemote` |
| `gamg.tuneRhoFail` | 0.9 | [tuneRhoHigh, 1) | failure condition | 6.3.5, D-066 | `tuneRhoFail` |
| `gamg.tuneConsecutiveWindows` | 2 | >= 1 | hysteresis windows | 6.3.5, D-066 | `tuneConsecutiveWindows` |

### c.precisionProfile (amendment D7, D-064)

`coupled.precisionProfile auto | dp | sp` (constant `precisionProfile`,
default `auto` = `sp` when sizeof(scalar) == 4, else `dp`). The profile sets
DEFAULT values only; an explicitly set keyword always wins, so the case
templates, which set these keywords, override it (the SP harness removes
them for SP runs). Values:

| keyword | dp profile (`dpProfile::...`) | sp profile (`spProfile::...`) |
|---|---|---|
| `coupled.convergence.residualTol` | `dpProfile::residualTol` = residualTol (1e-6) | `spProfile::residualTol` = 1e-5 |
| linear `tolerance` (absolute floor) | `dpProfile::tolerance` = tolerance (1e-8) | `spProfile::tolerance` = 1e-6 |
| `etaMin` | `dpProfile::etaMin` = etaMin (1e-3) | `spProfile::etaMin` = 1e-2 |
| `coupled.bounds.omegaMin` | `dpProfile::boundOmegaMin` = boundOmegaMin | `spProfile::boundOmegaMin` = 1e-5 |
| `coupled.bounds.kMin` | `dpProfile::kMin` = kMin (1e-12) | `spProfile::kMin` = 1e-10 |
| Anderson above `anderson.maxCells` | `dpProfile::andersonAboveMaxCells` = no | `spProfile::andersonAboveMaxCells` = yes |
| its history length there | `dpProfile::andersonLargeMaxM` = 4 | `spProfile::andersonLargeMaxM` = 4 |

## Compile-time items (not keywords, by design)

| item | where | reason |
|---|---|---|
| `kernelChunk` (512) | coupledDefaults.H, blockKernels.H | loop chunk of the hot kernels; a constant lets the compiler vectorise and unroll (spec 6.5) |
| `gateAmulFraction` (0.6), `gateAxpyDotFraction` (0.8) | coupledDefaults.H, Test-kernelBandwidth | pass/fail gates of a unit test, not solver settings; changing them is a test-threshold decision (DECISIONS.md) |
| `diagMaxLevel` (3) | coupledDefaults.H | number of implemented diagnostics levels (validation bound of `c.diagnostics.level`) |
| blockDim 4, blockSize 16, blockP 3; blockScalar float, reduceScalar double | blockScalar.H | layout and precision of the 4x4 (u, v, w, p) block system (D-001) |
| alignment 64 | alignedList.H | cache line / SIMD alignment |
| vec4 vector types and shuffle masks | block4Ops.H | compiler vector extensions |
| remediationFlag bit values | staticCriteria.H | output file format (documented above) |
| "Gauss linear uncorrected" momentum Laplacian | coupledAssembler.C | the non-orthogonal part is explicit by design (D-014, `c.nonOrthLimiter`) |
| guard epsilons (SMALL, VSMALL, ...) | everywhere | division-by-zero guards, not parameters (being replaced by typed guard constants, amendment D) |
| coupledFieldCompare large-deviation 0.05, centre tolerance 1e-9, time tolerance 1e-6 | utility | post-processing utility, not the solver |

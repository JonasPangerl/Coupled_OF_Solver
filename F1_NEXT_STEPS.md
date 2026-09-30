# coupledFoam on the F1 half-car: what to build next

Written 2026-09-30 at the end of the day that got coupledFoam running on the
20.6 M-cell F1 half-car at all (D-074..D-083, `F1_RUN_FINDINGS.md`). Goal set
by the user: **CFL 100-120 in the long run, and the fastest turnaround** -
the solver should decide itself which CFL is the most efficient at any
moment. Every number below is measured on that case, 60 ranks, EC2
c7a.16xlarge; the logs are in `~/f1/runs/s5_cpl_diag/log.coupledFoam.*`.

## Where the day ended

| run | step control | reached | limited by |
|---|---|---|---|
| attempt 6 (morning) | mRDM, global line search, continuity term on all cells | CFL 1-10, 300 it., 0 rollbacks | global omega (99 % of the steps < 1), then mRDM never recovered the cuts |
| local, clipping | mRDM + local mode (clipped increments) | CFL 9-23 | clipped dp broke continuity; linear solver fails at CFL 16-23 |
| fixes + mRDM | + failCeiling, stagnation stop | CFL 1.25 -> 1.36 in 26 it. | **mRDM deadlock**: R fell 140x during the post-rollback hold, then plateaued - no growth |
| fixes + EXP | EXP 10 %/iteration | **CFL 60 at iteration 178**, Cd 0.907, Cl -2.098 | the linear solver cost: 57 iterations, 43 s per iteration at CFL 59 |
| EFF, hill climbing | CFL per wall-clock second | CFL 2.7 -> 4.7 in 27 it. | the +-25 % scatter of the wall time at low CFL |
| **EFF, solver share + window** (final) | share of the linear solve -> 0.6, 4-iteration window | **CFL ~40 at iteration 105**, then 29-44 at 3.3-4.4 CFL/s | its own target (solverCost) and 3 rollbacks at CFL 24/32/44 |

Cost of one coupled iteration against the CFL (EXP run, iterations > 30,
without CFL cuts):

| CFL | linear iterations | seconds/iteration | CFL per second |
|---:|---:|---:|---:|
| 7 | 2.6 | 6.4 | 1.12 |
| 13 | 8.9 | 10.6 | 1.21 |
| 23 | 4.5 | 7.8 | 3.01 |
| 43 | 10.8 | 11.8 | **3.67** |
| 55 | 18.4 | 17.0 | 3.25 |
| 59-60 | 28-57 | 24-43 | ~1.4 |

A coupled iteration has a fixed cost of ~5.6 s here (assembly 2.3 s,
turbulence 1.8 s, the rest) and the linear solve on top. Up to CFL ~40 the
CFL is nearly free; above it the linear iterations grow faster than the CFL.
**The optimum of today's solver on this case is CFL ~40-45.** The cost also
depends on the flow state: the start phase needs 9-13 linear iterations at
CFL 13-18, the developed flow 4.5 at CFL 23.

## 1. The linear solver at high CFL - the actual limit for CFL 100+

With the current blockFGMRES + blockGAMG, CFL 100 would need an estimated
100+ linear iterations per step. Options, roughly in order of expected
gain per effort:

- **Measure first**: blockGAMG convergence factor per cycle against the CFL
  (diagnostics level 2 already records restart residuals). Is it the
  coarse levels (agglomeration of the pressure-velocity coupling), the
  smoother, or the Krylov restart?
- **Smoother strength at high CFL**: nFinestSweeps, blockILU0 vs a stronger
  block smoother on the finest level. At high CFL the momentum diagonal
  (V/dtau) shrinks and the saddle-point character of the 4x4 blocks grows.
- **A Schur-complement (SIMPLE-type) block preconditioner**:
  `blockSimplePrecon` exists; it treats exactly the saddle-point structure
  that dominates at high CFL. Compare against blockGAMG at CFL 40, 80, 120.
- **Krylov restart** 30 -> 50 against memory (restart 30 = ~20 GB on this
  mesh; the machine has 123 GB).
- **Inexact Newton**: etaMax 0.5 already; check whether a looser relative
  tolerance at high CFL lowers the cost per unit of residual reduction.

## 2. Step control

- **EFF objective**: pseudo-time per second is right while the flow still
  develops (R plateaus); near convergence residual reduction per second is
  the better measure. Switch or blend the objective by the residual trend.
- **EFF reversal on noise**: solved by the solver-share rule with the
  4-iteration window (D-083): a ratio within one iteration, averaged.
- **Hold after a cut** (nHold 10): too long with the ceiling (D-081) in
  place; make it proportional to how far the cut went.
- **A rollback in the start phase** (iteration 13 of the EFF run: 150 cells
  at |U| 1800 m/s, CFL 10.7). The start-up ramp (upwind -> case scheme,
  iterations 8-28) and the potential start overlap; a lower CFL ceiling
  during the ramp, or the ramp itself made longer, would avoid it.
- **Adaptive CFLmax**: once EFF is trusted, CFLmax 500 is just a safety cap.

### Seen in the hill-climbing EFF runs (afternoon of 2026-09-30, superseded)

- **EFF is too cautious far from the optimum.** After the reversals of the
  start phase its step factor had shrunk to a few per cent, and it needs
  five non-hurting steps to double it again: CFL 2.7 -> 4.7 over
  iterations 23-50, where EXP had reached 9.9 at 50. Better: step size
  from the measured slope of the efficiency (large while it rises
  clearly, fine only near the top), or a Brent/golden-section search on
  log CFL.
- **A ceiling from the start phase brakes the developed flow.** The
  rollback at iteration 13 set the ceiling to 8.6; at 1 % per iteration it
  capped the CFL at 15.6 at iteration 74 while the efficiency was still
  rising (1.4 -> 1.9 CFL/s, logged `capped:failCeiling`). Let the ceiling
  relax faster while the steps succeed and the efficiency rises, or let it
  expire when the start-up ramp ends.
- **The EFF decision lags one iteration** (the cost of an iteration is
  known at its end): `eff` and `linIters` on one CF| line belong to
  different iterations. Put the iteration number of the measurement into
  the reason, e.g. `down:efficiency(solverCost@it9)`.
- **The start phase has a reproducible instability**: with EFF the path to
  iteration 13 is deterministic and ends in the same rollback every time
  (150 cells, |U| 1801 m/s at CFL 10.7, while the ramp from upwind runs,
  iterations 8-28). A lower ceiling during the ramp, or a longer ramp.
- **The underbody flow switches at low CFL.** At CFL 2.7 (iterations 19-26
  of every EFF run) Cl swung -2.6 -> +1.5 -> -5.5 within a few iterations
  with a calm continuity residual and rU jumping 0.017 -> 0.043: the
  coupled solver followed a transient of the floor flow almost time-
  accurately. The EXP run at higher CFL went through without it. The floor
  is also where simpleFoam drifted between 300 and 1200 (-1842 -> -950 N).
  Check with written fields at 18, 20, 24 which region changes.
- **The quantile omega works**: at the iteration where the full-step
  version let 162 091 cells step unbounded (Cl -2.6 -> -7.6), the quantile
  version limited none; the dynamic set stayed at 47 k instead of 215 k.

### Seen in the final EFF run (solver share + window, 2026-09-30 14:52-16:05)

The run log is `~/f1/runs/s5_cpl_diag/log.coupledFoam` (copied to
`OpenFoam_cases/f1_halfcar/runs/s5_ec2_coupledFoam/`). What works: the
controller finds the cost optimum of the EXP table (CFL ~40, 3.7 CFL/s)
by itself, and every change is labelled - `down:...(solverCost)` and
`cut:...` never mixed. What to build next:

- **Most important: the share rule dips when the solve gets dearer for
  another reason.** Iterations ~230-305 (R 0.013 -> 0.002) the linear solve got dearer at the SAME CFL
  (8-21 linear iterations at CFL 29; the tolerance is relative, eta 0.5,
  and the linear residual falls with R, so it is not a fixed tolerance).
  The rule read that as cost of the CFL and stepped down 44 -> 29 -> 22 ->
  16, while the efficiency fell from 3-4 to 1.5 CFL/s: the steps down did
  not make the iterations cheaper. It recovered by itself (16 -> 29 over
  iterations 309-321, 3.1 CFL/s) once the solve got cheaper again - about
  70 iterations at half the efficiency. (The residual fell fastest in
  that dip, 0.007 -> 0.002: see the objective item below.) Fix: **verify every step** - compare
  the efficiency of the window after a change with the one before; a step
  that did not pay is undone and the target share at that state is
  corrected (or b estimated from the windows, see below). The rule finds
  the optimum only while the cost is a function of the CFL alone.
- **A stability limit at the same CFL as the cost optimum.** Sentinel
  rollbacks at CFL 24 (it. 58), 31.6 (it. 115), 43.6 (it. 206), each with
  7-13 k cells at |U| 570-6200 m/s. For CFL 100 the iteration itself must
  get more robust, not only the linear solve cheaper. First step: **log
  where the rollback cells are** (bounding box / patch-nearest of the
  cells over the bound) - the log says how many, not where.
- **Rollback cut too deep, hold too long.** A rollback divides the CFL by
  4 and holds 10 iterations; with the D-081 ceiling in place a factor 2
  and a hold of ~3 would lose less: after the rollback at 115 the
  efficiency sat at 1.3 CFL/s for 10 iterations while the solver share
  was 27 %.
- **The solve cost depends on the flow state, not only on the CFL.**
  Iterations 160-166: two steps down (37 -> 29) and the linear iterations
  stayed at 12-17. The rule assumes cost ~ CFL^b; a down step whose cost
  did not fall should be undone. Measure b from the windows instead of
  assuming 1.6 (a regression of log tSolve on log CFL over the last
  windows), then the target share 1/b follows.
- **The target share 0.6 is conservative in the developed flow.** Cost per
  iteration against the CFL in this run (iterations > 30, no cuts):
  CFL 30-40 10.7 linear it. / 11.3 s, CFL 40-50 10.9 / 11.5 s - flat,
  so the optimum lies above 45; only the rollbacks keep it there.
- **Objective near convergence**: the residual plateaus at 0.01-0.02
  from iteration ~100 while Cl swings +-0.16 (last 50 iterations).
  Pseudo-time per second is the wrong objective then; residual (or force)
  reduction per second, or a switch to a transient/averaging mode for the
  unsteady underbody flow.

## 3. Robustness mechanisms of today, to be settled as defaults

All of today's switches except D-074, D-077 (keep), D-078 are off by
default, so every verified case is bit-identical. Before any of them
becomes a default, run the full battery (T0-T3, restart) with it on:

| switch | F1 setting | default today |
|---|---|---|
| `lineSearch.mode` (quantile, `localFraction` 0.0005) | local | global |
| `ptc.cflStrategy` / `effTol` | EFF / 0.10 | mRDM / 0.02 |
| `ptc.failCeiling` / `ceilingRelax` | 0.8 / 1.03 | 0 (off) |
| `solvers.coupled.stagnationRatio` | 0.9 | 0 (off) |
| `potentialPressure` | bernoulli | keep |

- The **dynamic set grows large with the local mode** (up to 2.4 % of the
  cells in the clipping version, 300 k cells in the full-step version):
  every cell over the bounds joins with a halo layer and stays nQuietIters.
  Those cells run upwind: check the accuracy cost (component forces with
  and without), and consider joining without the halo.

## 4. Case set-up and output (F1 template)

- `forcesPerPatch` with `writeControl timeStep; writeInterval 10;` from the
  start: component histories cost ~0.5 s/iteration then, every iteration
  cost +75 %.
- coupledFoam writes ~30 diagnostic fields per write time (~7 GB on this
  mesh); a production run needs a smaller set.
- The EnSight export and the monitor CSVs (`OpenFoam_cases/tools/
  export_monitors.py`) for coupledFoam runs as for simpleFoam.

## 5. The comparison itself

simpleFoam on the same mesh: Cd 1.000, Cl -2.173 at 1600, Cl stable at
-2.2 +- 0.05 since 1200, Cd still drifting (~5 % over 400 iterations). The
coupledFoam runs of today were at Cd 0.91-0.97, Cl -2.09 to -2.16 after
100-180 iterations. A converged comparison needs both runs to settle; the
component forces (report page 1 of `OpenFoam_cases/f1_halfcar/runs/
s5_ec2_coupledFoam`) show where they differ.

## 6. Housekeeping

- Nothing of today is committed in this repository (the working tree also
  carries uncommitted changes from before 2026-09-30).
- `tests/test_unit.py::test_blockGAMG_cycles` fails as before today
  (V 12, F 8, W 7, K 11); T3-GEKO fails as before (Cd 0.078 vs 0.034).

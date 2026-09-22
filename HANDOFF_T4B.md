# T4b collapse (B4 abort at outer iteration 65) - handoff

State 2026-09-22 12:40. Author: Opus agent (T4b investigation), stopped by
the lead's plan change. Nothing is committed on main; the main install is
untouched.

## 1. Root cause (confidence: high)

**The block ILU0 smoother of blockGAMG breaks down on individual pivots.**
The LDU-native factorisation modifies the diagonal blocks,
D*_c = D_c - sum L D*^-1 U. On the saddle-point block rows the modified
4x4 pivot of a single cell can lose almost all of its dominance. Its
inverse is then up to 223 times larger (max-abs norm) than the inverse of
the unmodified block. That happens mostly on GAMG level 1 (a coarse cell),
sometimes on the finest level. The damped ILU0 sweep (relaxation 0.5)
then amplifies the residual by 1e3 to 1e5 per smoother call, with about
half of the residual norm sitting in that one cell. The V-cycle becomes
garbage, and FGMRES stagnates (0.85 reduction in 400 iterations).

- CFL cuts do not help. The degenerate pivot comes from the continuity
  and pressure couplings, not from V/dt.
- The aggregates stay fixed between re-agglomerations
  (reagglomerateInterval 50). So the same coarse cell degenerates again
  in every retry: 3 failures, then the B4 abort.

Evidence. All numbers are from the instrumented debug build. The
trajectory is bit-identical to the original run: R at iterations 49, 54
and 64 matches the original log in all digits.

- Per-iteration maximum pivot growth max|D*^-1|/max|D^-1| against the rho
  spikes of the original run:

  | iteration | pivot growth | rho |
  |---|---|---|
  | 9 | 14.3 | 19 |
  | 10 | 13.1 | 3027 |
  | 11 | 79 | 945 |
  | 12 | 41 | 36 |
  | 17 | 8.5 | 18 |
  | 27 | 36 | 2.8e5 |
  | 40 | 31 | - |
  | 49 | 38.6 | 973 |
  | 60 | 15.4 | 9.8 |
  | 65 | **223** | abort |

  Healthy iterations stay at 3.7-6.6 (3.7-3.9 at CFL 500).
- Iteration 65:
  - Location: rank 6, GAMG level 1 (18840 local coarse cells), cell
    7565, pivot growth 222.5.
  - Every smoother call amplified the residual by about 900-1000x, once
    by 98776x. The cell held 50-53 % of the residual norm.
  - The unmodified diagonal block of that cell is regular. The
    fallback pivot is therefore fine.
- Level-2 diagnostics (run/t4b_diag1, cf_diag build, level 2):
  - Level-1 post-smoothing grows the residual, e.g. 4.6 -> 807 at iter
    11, 40.9 -> 2.9e5 at iter 27, 3.77 -> 22.7 at iter 60.
  - The finest-level pre-smoother at iterations 9/10 does the same
    (1 -> 6.5).
  - The coarsest BiCGStab is fine: 14-20 iterations, and it keeps its
    best iterate.
- **Decisive A/B on the exact failing system** (run/t4b_reproB):
  - Same trajectory, with the safeguard (pivotGrowthLimit 5) switched
    on only from finest update 65 (env CF_PIVOT_START=65).
  - Iteration 65 converged in 7 FGMRES iterations. Without the
    safeguard it failed 3 x 400 iterations.
  - The run then went on healthy to 75 (3-4 linear its, no cuts).
- Restart side test:
  - A restart from the written iteration-64 fields (run/t4b_rs64)
    re-agglomerates on the iteration-65 matrix.
  - Iteration 65 then passed (6 its): the degenerate coarse cell exists
    only with the old aggregates.
  - Restart from 48 with interval 50 (run/t4b_rs48a) ran 49-75 without
    failure. The trajectory differs from iteration 50 on, so this
    neither confirms nor refutes "stale aggregates" on its own. The
    pivot data does.

## 2. Fix (implemented, T4a-verified, not yet on main)

Branch `t4b-fix`, worktree `/home/jonas/cf_fix`, commit e3d452e ("wip:"),
based on main ba01275 (with the merged diagnostics). It is built into
`/home/jonas/cf_fix_platform` (lib, coupledFoam, test apps, coupledFieldCompare).

1. **blockILU0 pivot-growth safeguard (D-049 candidate)**
   - Where: `src/blockSolvers/blockSmoothers/blockILU0.{C,H}`.
   - Rule: if max|D*_c^-1| > pivotGrowthLimit * max|D_c^-1|, the cell
     uses the unmodified block D_c^-1. The fallback pivot also enters
     the later modifications, so the factorisation stays consistent.
   - Keyword `pivotGrowthLimit` in the blockGAMG dictionary. Default
     `coupledDefaults::iluPivotGrowthLimit = 5`; 0 = off. The value is
     validated (finite, >= 0).
   - Count: `blockSmoother::nPivotFallback()` (virtual, 0 for the other
     smoothers), summed by `blockGAMG::nPivotFallback()`. It is printed
     as `nPivFb=` on the CF| line only when > 0. The summary JSON gets
     `pivotFallbacks`.
   - Cost: one extra 4x4 inversion per cell per update, negligible.
2. **startupControl** (the user's hybrid spec, D-048 candidate)
   - Where: `src/control/startupControl.{C,H}`, wired into coupledFoam.C.
   - `startupMode hybrid | upwind | none`, default hybrid. Keywords:
     `startupRampStart` 10, `startupRampLength` 20,
     `startupRampStartMax` 50, `startupSwitchR` 1e-2,
     `startupStagnationFactor` 0.9, `startupStagnationWindow` 5,
     `startupFastFactor` 0.7, `startupDevelopedTol` 1e-2,
     `startupUpwindIters` 50 (upwind mode only). All are constants in
     coupledDefaults.H.
   - Beta ramp: beta_m = clamp((m - n0)/L, 0, 1).
   - The ramp start n0 is the first accepted iteration with one of these
     triggers:
     - relativeR: R/R1 < 1e-2.
     - stagnation: before 10, R_n > 0.9 R_(n-5).
     - fixed: at 10, if upwind is not dropping R fast.
     - delayedFastDrop: between 10 and 50, once R_n >= 0.7 R_(n-5).
     - cap: at 50.
   - developedStart:
     - Condition: fresh start with potentialInit off and a non-uniform
       U.
     - Iteration 1 is assembled with beta 1. If max(rU, rp) < 1e-2
       (scale-free residuals, not R/R1), beta is 1 from iteration 1.
     - Otherwise iteration 1 is re-assembled with the start-up beta.
   - upwind: the old behaviour exactly (checked: bit-identical
     trajectory in run/t4b_repro). none: beta 1 from iteration 1.
   - The state (n0, trigger, R history) is in coupledState. States
     written before D-048 restart as "restartDone".
   - Logging: one Info line per trigger. The summary JSON gets
     startupMode, startupTrigger, rampStartIter and rampEndIter.
3. **ptc.linFailPolicy strict | reduction** (the lead's maxIter-cap
   request)
   - Default strict, which is the old B4 definition.
   - reduction: a capped solve is accepted if it is finite and
     final <= linAcceptReduction (0.9) * initial.
   - The summary gets linFailPolicy and linAcceptedUnconverged.
   - Built and compiled, **not tested in any run**.

The debug worktree `/home/jonas/cf_t4b` (branch t4b-debug, commit
464028b) has the same features plus temporary instrumentation. Do NOT
merge it. Instrumentation:

- CF_ILU_DEBUG=1: per update, pivot-growth statistics; per smooth call
  with more than 2x residual growth, the worst cell with its block row.
- CF_PIVOT_START=N: safeguard only from the N-th finest update on.

## 3. Verification done

- **T4a np10 pytest with the fix build (cf_fix_platform): PASS.**
  - Hybrid default: the trigger was relativeR at 24 (delayed by the
    fast upwind drop), beta 1 from 44.
  - Result JSON: /home/jonas/cf_fix/results/tests/T4a_np10.json (copy
    /home/jonas/t4b_logs/T4a_np10_fix.json). Run dir
    run/T4a_np10_fix. The baseline run dir is restored as run/T4a_np10.
    Main's results/tests/T4a_np10.json is untouched.

  | metric | fix | baseline (ca33547) |
  |---|---|---|
  | pass | True | False (Cd 2.2 %) |
  | stationary at | 450 its | 450 its |
  | wall / CPU-h to stationary | 422 s / 1.17 | 395 s / 1.10 |
  | solver total wall / CPU-h (800 its) | 753 s / 2.09 | 710 s / 1.97 |
  | Cd / Cl window mean | 0.4032 / 0.0794 | 0.4051 / 0.0771 |
  | vs ref (0.3964 / 0.0768) | Cd 1.71 % (tol 2 %), Cl abs 0.0027 (tol 0.01) | Cd 2.20 %, Cl 0.5 % |
  | field RMS dU / dp | 0.44 % / 0.10 % | 0.46 % / 0.10 % |
  | linear its total / max per solve | 1157 / 5 | 1127 / 8 |
  | cuts / rollbacks | 0 / 0 | 0 / 0 |
  | pivot fallbacks | 1-2 cells in 48 of 800 its | - |
  | speedup vs simpleFoam | 1.86x | 1.99x |

  Some of the ~7 % slower time to stationarity may be machine load:
  another agent's coupledFoam (cf_osc) and a T3 job ran at the same
  time.
- Not yet done:
  - T4b full 800-iteration run with the fix.
  - T0/T1/T2 regression checks.
  - Unit battery on the fix build.
  - (a) plain ramp vs (b) delayed ramp comparison (T1/T2/T4a).
  - The mapped T4a start.
  - W-cycle and maxIter-10 variants.

## 4. Anomaly checks from the brief

- **"nLocLim = 10*iter" is not in the data.** run/T4b_np10/log.coupledFoam
  has nLocLim 2061, 1145, 2470, 4400, 5102 ... 165090 at iteration 29
  (right after the start-up switch) and 190k-260k (12-15 % of the cells)
  at CFL 500. T4a similarly has 80-113k (23-32 %) from iteration 26 on.
  The counter is reset per call and recomputed from a fresh rDeltaTV
  every assembly; no state is carried over. The "938-5993" for T4a are
  only its first iterations. No per-rank one-cell growth exists. The
  lead probably read another column or file.
- Local limit and CFL cuts:
  - A locally limited cell gets dt = dUmax V/|r_P|, which does not
    depend on the CFL. The CFL cuts therefore change only the
    unlimited ~85 % of the cells.
  - The limited cells only increase the diagonal. They are not the
    cause.
- Static set (8251 cells, constant), float range and reagglomeration
  weights: not implicated. Clamped 0, no singular pivots (pivotGuard),
  and the condition of the unmodified blocks is <= 7e6.

## 5. Side findings (not fixed)

1. **Uref comes from the potential-flow field.**
   - lineSearch::setReference takes max|U| over the field at iteration
     1. With potentialInit, that is the singular potential peak at
     sharp edges: T4b Uref 176.5 (8.8 U_inf 20), T4a 85.9 (4.3 U_inf).
   - Consequences:
     - Dynamic-set thresholds cU*Uref = 706 m/s and
       pref = 0.5 Uref^2 = 15571, so the dynamic remediation never
       fires on T4 (nDyn = 0 always). This defeats D-044.
     - Line-search fU*Uref = 53 m/s per step.
     - Local-limit dUmax = 88 m/s.
     - Sentinel limit 10 Uref = 1765 m/s.
   - Proposed: Uref = max|U| over the non-coupled boundary values (the
     free stream), and the field max only as a fallback if that is ~0.
     This changes T1-T5 trajectories, so it needs its own verification.
2. **sentinel::writeLastValid writes into the current time directory, not
   `<iter-1>_lastValid`.** v2606 regIOobject::writeObject resets a
   non-time instance to time().timeName(). The iteration-64 fields of
   the abort therefore land in `processor*/65` next to remediationFlag
   and sentinelFlag, while the message says "written to 64_lastValid".
   A restart would pick up 65 without a coupledState. Fix: write with
   an explicit OFstream to `<case>/<inst>/<name>` (header plus
   operator<<), or rename the directory after writing.
3. T4a with the ca33547 defaults never gets below R ~0.017-0.024 in 800
   iterations; it is "stationary" by the force window only. This is the
   same outer-iteration limit cycle family as FABLE 4/4c.

## 6. Open requests from the lead (status)

- Ramp (a) vs delayed ramp (b) on T1, T2 and T4a, wall and CPU-h.
  - Not run. The runner is ready: scratchpad `t4b_variants.py`. It runs
    serial variants concurrently and records toR1e-5, wall, CPU,
    trigger, pivotFallbacks. It needs CF_RUN_ROOT=/home/jonas/coupledFoam/run
    and the cf_fix_platform env.
  - Defaults stay at (b) delay (rampStartMax 50) until measured. The
    lead's rule: make the delay the default only if it is not slower
    on any case. Otherwise set startupRampStartMax = startupRampStart.
- T0 pytest --ranks 1 with the hybrid default: not run. If it fails,
  pin cases/T0_cavity/system/fvSolution to `startupMode upwind`.
- Mapped T4a start:
  - The case is prepared: run/t4a_mapped. It was created from
    T4a_mesh, potentialInit no, maxIter 300, force stop window 100.
    processor*/0 = fields of run/T4a_np10/800 (the baseline run).
  - Run it with Allrun -restart (no fields() reset). The developed
    trigger should fire at iteration 1: T4a stationary state rU/rp
    0.0033/0.0025 < 1e-2; potential start 0.13/0.033.
- Residual stop criterion for mapped starts (item c, report only):
  - convergence.residualTol is relative to R1, the residual of the
    first iteration. A mapped or near-converged start has a small R1,
    so R/R1 never reaches the tolerance.
  - Proposal: an absolute scale-free criterion on
    Rraw = ||b - Ax||/normFactor, which is already normalised by
    sum|Ax| + |b|. For example, a residualTolAbs next to residualTol.
    Or carry R1 over from the source run in coupledState when mapping.
- W cycle with nPreSweeps 0 / nPostSweeps 2: not tested. The r_l >= 3
  ratio rule should pass on T4b (V hierarchy ratios 8.8/8.7/8.7 with
  mergeLevels 3). With the pivot safeguard the smoother no longer
  amplifies, so the W cycle can be judged on its merits.
- maxIter 10 / minIter 2 + linFailPolicy reduction: code ready (default
  strict), untested.

## 7. What I would do next

1. Merge `t4b-fix` into main:
   - Commit message `fix: blockILU0 pivot-growth safeguard ...` with the
     Co-Authored-By trailer.
   - Keep startupControl and linFailPolicy as separate commits if
     preferred. Their code is in the same WIP commit; split with
     `git checkout e3d452e -- <files>`.
   - Rebuild the main install: wmake lib, then wclean + wmake for the
     apps and tests.
2. Unit battery np1+np4 and T0 --ranks 1. Then T1 and T2 --ranks 1:
   - The safeguard can trigger there too; compare with
     results/tests/*.json.
   - Run the (a)/(b) ramp comparison.
3. T4b pytest (test_T4[b] --heavy, CF_HEAVY_NP=10 CF_MPI_BIND=core).
   Expect no B4 abort; nPivFb shows where the safeguard acts. This also
   runs the 2000-iteration simpleFoam reference continuation, about
   1.6 h.
4. Then the Uref fix (section 5.1) and the lastValid write fix (5.2).

## 8. Run directories (all under /home/jonas/coupledFoam/run)

| dir | build | content |
|---|---|---|
| T4b_np10 | main ca33547 | original failure |
| t4b_diag1 | cf_diag, level 2 | identical failure at 65; diagnostics/diag.rank*.jsonl (iter 65 record missing, abort); written times 16/32/48/64 |
| t4b_rs64 | cf_t4b | restart from 64: iteration 65 OK (re-agglomerated) |
| t4b_rs64b | cf_t4b | prepared (restart 64, pivotGrowthLimit 10), not run |
| t4b_rs48a | cf_t4b | restart 48, interval 50: OK to 75 |
| t4b_rs48b | cf_t4b | prepared (interval 1), not run |
| t4b_repro | cf_t4b, startupMode upwind, CF_ILU_DEBUG | bit-identical reproduction; ILUDBG lines in log.coupledFoam |
| t4b_reproB | cf_t4b, CF_PIVOT_START=65, limit 5 | A/B: 65 converges in 7 its, runs to 75 |
| T4a_np10_fix | cf_fix | T4a pytest with the fix (PASS) |
| t4a_mapped | - | prepared mapped start, not run |

Scratchpad helpers (Windows scratchpad of this session):

- `mkcase.py` (case_from_mesh with sets)
- `t4b_launch.sh` (detached Allrun with a private install)
- `t4b_set.sh` (fvSolution entries)
- `t4b_build.sh` (NJ=jobs; lib + wclean/wmake app [+ tests])
- `t4b_pytest.sh` (detached pytest from a worktree with a private
  install, CF_RUN_ROOT = main run/)
- `t4b_variants.py`, `t4b_cfstat.py`
- `t4b_mapped.sh`
- patch scripts `fix_ilu.py` and `fix_main.py` (the port onto main with
  diagnostics)

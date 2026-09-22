# Work order for Opus - 2026-09-22 morning (written by Fable)

Prescriptive task list for FABLE_REVIEW item 4 and the remaining
pipeline. Follow it literally; where a number or keyword is given, use
exactly that. Fable wrote this document only and touched no code.

User decisions of 2026-09-22 morning (authorised, quote translated):
"For the motorBike and Ahmed body cases we probably have to average and
compare that, and maybe use coarser tolerances, because of the wake
problem. That is okay." -> TASK 2 implements exactly this. No OTHER
threshold or criterion may be loosened.

Read first: FABLE_REVIEW.md (items 2c and 4), DECISIONS.md D-038,
D-039, D-040, HANDOFF.md section 7.

## 0. Live state - check before anything else

- `run/T4b_conservative` (mpirun np10, started 07:38) runs coupledFoam
  on the 1.70 M-cell T4b mesh with the CONSERVATIVE settings:
  cycleType V, agglomerationWeights geometric, nFinestSweeps 2,
  nCellsInCoarsestLevel 200, smoother blockILU0, FGMRES restart 30,
  CFLmax 100, maxIter(outer) 1500, budget endTime 1500.
  At iteration 76 it was healthy: linIters 3, R falling, ~4.5 s/iter
  (tSolve ~2 s), ETA ~09:40.
  - If finished cleanly: record from `log.coupledFoam` the iteration
    count, final R, last tWall; compute mean and std of Cd and Cl over
    the last 500 rows of `postProcessing/forceCoeffs/0/coefficient.dat`;
    keep the case directory (it is the current best T4b evidence).
  - If it aborted: record the abort message verbatim; the failure
    iteration and CFL matter for TASK 1.
  - Do not start any other 10-rank job while it runs (10-core budget is
    the TOTAL for everything).
- No other solver processes should be alive. Kill only PIDs whose
  /proc/PID/cwd is inside a run/ directory you own.

## TASK 1 - Robust linear-solver defaults (FABLE item 4 core)

Problem, precisely: the D-039 settings (K cycle, agglomerationWeights
combined, nFinestSweeps 1, nCellsInCoarsestLevel 20, autoTune no) were
tuned on T1 only. Measured failures elsewhere, all "B4 abort: 3
consecutive linear-solve failures", i.e. FGMRES did not reach
eta*||r0|| within maxIter:
- T2: abort at outer iter 306 (reduction 0.29 vs eta 0.225, 200 its);
  with maxIter 400 it survives to the 1000-iter budget.
- T3-SST: aborts at ~260 and ~434 (reduction 0.51 vs eta 0.5, 200 its).
- T4a (354k snappy): abort at outer iter 51, CFL 500 - reduction only
  0.64 in FOUR HUNDRED linear iterations. Raising maxIter cannot fix
  this; the preconditioner stalls.
- The conservative set (V/geometric/2 sweeps/coarsest 200) needs only
  ~3 linear its per solve on T4b at CFL 20 (see section 0), so the
  regression is caused by one or more of the four D-039 knobs, not by
  the meshes.

### 1a. Dump linear systems from the failing cases (serial only)

The dump feature: dictionary entry `coupled.dumpLinearSystem` = list of
outer iterations, e.g. `(50 200 400 800)`; files appear as
`<case>/linsys/iter<N>`. Serial runs only.

For every foamDictionary call in this task use
`foamDictionary -disableFunctionEntries` (without the flag the
relaxation include handling breaks - hard-learned).

1. T2 dumps: copy `cases/T2_backwardFacingStep2D` to `run/dump_T2`; set
   `coupled.dumpLinearSystem (50 200 400 800)`, `coupled.maxIter 800`,
   controlDict endTime/writeInterval 800; run
   `CF_MPI_BIND=none nice -n 19 ./Allrun -solver coupledFoam -np 1`.
2. T3 dumps: copy `cases/T3_airFoil2D` to `run/dump_T3`; set
   `coupled.dumpLinearSystem (30 100 250 400)`, maxIter 400, endTime
   400, turbulence kOmegaSST (Allrun flag `-turbulence kOmegaSST`). A
   B4 abort after the last dump is acceptable; dumps written before an
   abort survive.
3. T4a dumps: the cached mesh is decomposed for 10 ranks; materialise
   the serial variant first with the existing helper:
   `python3` (venv, from repo root, tests/ on sys.path):
   `from test_T4_motorBike import mesh_dir; mesh_dir("T4_motorBike", "T4a_mesh", ["-mesh", "a"], 1, 10)`
   -> creates `run/T4a_mesh_np1` via reconstructParMesh. Then copy
   `cases/T4_motorBike` to `run/dump_T4a`, copy the mesh in with
   `test_T4_motorBike.copy_mesh(run/T4a_mesh_np1, run/dump_T4a)`, set
   `coupled.dumpLinearSystem (20 40 50)`, maxIter 60, endTime 60,
   numberOfSubdomains 1; run Allrun `-np 1 -mesh a`. (~5 min serial;
   the failing regime was iter ~51 at CFL 500, which is the case
   default CFLmax, so do NOT lower CFLmax here.)
4. T1 dumps already exist: `/home/jonas/cf_precond/run/T1_base/linsys/`
   iter30, iter150, iter300, iter450. Copy that linsys directory to
   `run/dump_T1_linsys/` before deleting worktrees (section 7).

### 1b. Offline variant grid with Test-blockSystem

The harness: `Test-blockSystem -system <linsys file> -variants <dict>`
(also `-only <name>`, `-dominance`, `-diagnose`, `-equilibrate`,
`-scalars "(nf eta cfl)"`, `-deepRelTol <v>`). A variants dictionary
looks like (real example from the D-039 study):

    FoamFile { version 2.0; format ascii; class dictionary; object variants; }
    base   { }
    wF1c20 { blockGAMG { agglomerationWeights combined; nFinestSweeps 1; nCellsInCoarsestLevel 20; } }

Create `run/_variants_robust` with EXACTLY these variants (names
verbatim; every entry inside `blockGAMG { ... }`):

| name | cycleType | agglomerationWeights | nFinestSweeps | nCellsInCoarsestLevel |
|---|---|---|---|---|
| V_geo_f2_c200 (baseline conservative) | V | geometric | 2 | 200 |
| V_comb_f2_c200 | V | combined | 2 | 200 |
| V_comb_f1_c200 | V | combined | 1 | 200 |
| V_comb_f1_c20 | V | combined | 1 | 20 |
| V_geo_f1_c200 | V | geometric | 1 | 200 |
| K_comb_f1_c20 (= D-039) | K | combined | 1 | 20 |
| K_comb_f1_c200 | K | combined | 1 | 200 |
| K_geo_f2_c200 | K | geometric | 2 | 200 |
| K_comb_f2_c200 | K | combined | 2 | 200 |
| GS_control | V | geometric | 2 (smoother blockGaussSeidel) | 200 |

Run every variant against EVERY dump from 1a (T1 iter300+iter450, T2
iter400+iter800, T3 iter250+iter400 - or the last two that exist, T4a
iter40+iter50). Save each harness stdout to
`run/robust_<case>_<iter>.log`. Also run `-dominance` once per case on
the latest dump and keep the output.

### 1c. Selection rule (no judgement calls)

A variant PASSES a dump if the harness's FGMRES solve reaches a 0.1
residual reduction within 30 Krylov iterations AND rhoOpt <= 0.9.
Select the variant that passes ALL dumps and has the smallest total
wall time summed over the dumps. Tie-break: fewer sweeps, then larger
nCellsInCoarsestLevel is preferred (cheaper, more robust setup).

Expected outcome (hypothesis to confirm, not to assume): a V-cycle
combined-weights variant with nFinestSweeps 2 or coarsest 200 passes
everywhere; K_comb_f1_c20 fails T3/T4 dumps.

### 1d. Only if NO variant passes all dumps: damped smoother

Implement under-relaxed smoothing: in `src/blockSolvers/blockSmoothers`
(blockILU0 and blockGaussSeidel), multiply the per-sweep update by a
factor `smootherRelaxation`, read from the blockGAMG dictionary,
default `coupledDefaults::smootherRelaxation = 1.0` (bit-identical to
today at 1.0; new named constant with a comment; no magic numbers).
Add variants `V_comb_f2_c200_r07` and `V_comb_f2_c200_r05`
(smootherRelaxation 0.7 / 0.5) to the grid and re-run 1b/1c. This was
D-039 open idea 5.

### 1e. Apply the winner

1. Set the winning variant's four keywords in the fvSolution of
   cases/T1_pitzDaily, T2, T3, T4_motorBike, T5_ahmed. T0_cavity stays
   UNTOUCHED (validated, uses blockGaussSeidel; see D-039).
2. autoTune stays `no` everywhere (D-039: the unscaled rho misleads
   it). reagglomerateInterval 50 stays. restart 30 stays. Linear
   tolerance and maxIter 400 stay as committed.
3. Exception rule: if the winner makes T1 np1 slower than 55 s wall
   (pytest T1 or a direct Allrun run), keep the D-039 set in
   T1_pitzDaily ONLY, as a documented per-case exception.
4. Write DECISIONS entry D-043 with the full grid table (variant x dump
   -> pass/fail, its, rhoOpt, wall) and the selection.

### 1f. Acceptance battery for TASK 1 (all must be met before push)

Build hygiene first: after any src/ change run wclean+wmake for the
library AND `wclean; wmake` in applications/coupledFoam and every
applications/test/* (a stale app binary produced false results once
already - FABLE item 3).

1. Unit: Test-blockGAMG cycles V,F,W,K x np1,np4 on run/unit_cavity,
   `-tolerance 1e-9 -skipDiagonal` (np4: `mpirun --bind-to none -np 4
   ... -parallel`): all 8 PASS. Test-blockMatrix np1+np4,
   Test-blockFGMRES, Test-block4Ops, Test-doubleReduce, Test-precision:
   PASS. Test-kernelBandwidth: gate stays PASS.
2. pytest tests/test_T0_cavity.py --ranks 1,4: 4/4 PASS. (Iteration
   counts at Re1000 scatter ~2x between runs - that is FABLE item 1
   chaos, NOT a regression; only a FAIL status is a regression.)
3. pytest tests/test_T1_pitzDaily.py --ranks 1: expected to still FAIL
   the 400-iteration budget (~495 needed). Record iterations, wall,
   CPU-h. Do NOT change the budget (user has not approved that).
4. pytest tests/test_T2_backwardFacingStep.py --ranks 1: record. The
   0.5 % xr tolerance stays (T2 is not a wake case).
5. T4a full run (np10, via pytest tests/test_T4_motorBike.py --heavy
   -k a) with TASK 2 in place: must run to budget or stationarity, NO
   B4 abort.
6. T4b: reuse run/T4b_conservative if it completed; else rerun with
   the new defaults. NO B4 abort.

## TASK 2 - Averaging criterion for wake cases (USER-APPROVED)

Scope: T4a, T4b, T5 only (physically oscillating bluff-body wakes).
T0-T3 keep their existing criteria unchanged.

Evidence: the T4b simpleFoam reference (run/ref_T4b_np10, 4000 its,
11876 s, 33.0 CPU-h) oscillates Cd 0.393-0.405, Cl 0.061-0.072
forever; the 12.3(ii) 0.2 % min/max window can never be met by any
steady solver on these cases (FABLE item 2c). T4a reference
(run/ref_T4a_np10, 3000 its, 1516 s, 4.21 CPU-h) likewise. Reference
means, T4b iterations 3500-4000: Cd 0.3993, Cl 0.0671.

Implementation (bench/run_bench.py + tests/test_T4_motorBike.py, the
T5 test inherits via imports; make_report.py + both papers updated):

1. Mark the cases: in run_bench.py CASES, add `"oscillatory": True` to
   T4a, T4b, T5 (absent means False).
2. New function `stationary_mean(hist, quantity)` (place it in
   run_bench.py so tests and bench share it):
   - window W = max(500, n//4), capped at n//2 (n = iterations run).
   - mean m over the last W samples; halves m1 (first W/2), m2 (last
     W/2).
   - stationary iff |m1 - m2| <= max(0.005*|m|, 0.002), evaluated for
     Cd AND Cl separately; both must hold.
   - iters_to_stationary = the smallest N <= n such that the window
     ending at N is stationary (scan in steps of 50); None if never.
3. For oscillatory cases replace, in the test comparison and in the
   benchmark's time-to-convergence:
   - "converged" (both solvers) := stationary by rule 2. Use
     iters_to_stationary for wall/CPU-to-convergence exactly where
     iters_to_conv is used today (potentialFoam offset rule D-033
     unchanged).
   - quantity comparison: window means. PASS iff
     |Cd_mean - Cd_ref_mean| <= max(0.02*|Cd_ref_mean|, 0.002) AND
     |Cl_mean - Cl_ref_mean| <= max(0.02*|Cl_ref_mean|, 0.01).
     (Cd: 2 % relative - the user-approved coarser tolerance, was 1 %.
     Cl: 0.01 absolute floor because Cl_ref ~ 0.067.)
   - record in the JSON: means, stds, W, iters_to_stationary for both
     solvers, and `criterion: "stationaryMean (D-042)"`.
4. The coupledFoam runtime stop criterion (convergence dict) is NOT
   changed by this task; runs stop on budget as configured by the
   tests. Only the evaluation changes.
5. DECISIONS entry D-042, paste verbatim:

   > ## D-042 - Averaged force comparison for wake cases (user decision, 2026-09-22)
   >
   > T4a, T4b and T5 have physically oscillating wakes: the T4b
   > simpleFoam reference oscillates Cd by +-1 % and Cl by +-8 %
   > indefinitely (run/ref_T4b_np10), so criterion 12.3(ii) (0.2 %
   > min/max window) is unsatisfiable for any steady solver on these
   > cases, coupledFoam and simpleFoam alike. The user approved
   > averaging and coarser tolerances for exactly these cases ("bei den
   > motorbike und ahmed body cases muessen wir wahrscheinlich averagen
   > und das vergleichen und vllt auch groebere toleranzen nehmen wegen
   > dem nachlauf problem. das ist okay", 2026-09-22). New evaluation,
   > oscillatory cases only: window W = max(500, n/4) (<= n/2);
   > stationarity = the two half-window means differ by <=
   > max(0.5 % |mean|, 0.002) for Cd and Cl; comparison = window means,
   > Cd within max(2 % relative, 0.002 absolute), Cl within
   > max(2 % relative, 0.01 absolute) of the reference. Time to
   > convergence uses the first stationary window. T0-T3 criteria are
   > unchanged. This is a user-approved criterion change, not an agent
   > loosening.

6. Report/paper: the validation table shows mean +- std for oscillatory
   cases and footnotes the criterion with D-042; the tutorial paper
   gets two sentences in its "how to read" part (its agent's
   placeholder mechanics handle the numbers).

## TASK 3 - T3-SST limit cycle (after TASK 1 lands)

1. Rerun T3-SST with the TASK-1 defaults, serial, maxIter(outer) 3000:
   if R < 1e-6 before 3000, T3 is fixed; run pytest T3 (both models)
   and record.
2. Else run these three probes in parallel (serial each, nice 19,
   budget 1500 outer its each; copies of cases/T3_airFoil2D under
   run/exp_T3_*; remember -disableFunctionEntries):
   a. `coupled/ptc/CFLmax 20`
   b. k/omega relaxation 0.5 (edit the copy's
      system/relaxation.coupledFoam before Allrun: k 0.5; omega 0.5)
   c. `coupled/startup/upwindIters 200` - read the exact keyword names
      from src/include/coupledDefaults.H (startupUpwindIters 50,
      startupSwitchR 1e-2) and applications/coupledFoam createFields /
      coupledFoam.C to find the dictionary path; if there is no
      dictionary override for them, add one (getOrDefault against the
      coupledDefaults value, same pattern as every other control).
   Record for each: R history minimum and final, iterations of
   `bounding k` events (grep -c "bounding k" log), nLocLim trend, Cd/Cl
   window means vs reference (Cd 0.0907, Cl 0.253).
3. If a probe converges: adopt its setting for T3's template ONLY,
   document in D-044 with the probe table. If none: write the findings
   to FABLE_REVIEW item 4 (T3 sub-item) and move on - T3 must not
   block TASK 4.
4. The T3 dumps from 1a plus `-dominance`/`-diagnose` output belong in
   the D-044 (or FABLE) analysis either way.

## TASK 4 - Remaining pipeline (strict order, budget-aware)

1. T5 coarse: `heavy_prep.sh T5coarse` equivalent - mesh once into
   run/T5_mesh_coarse (np10) and simpleFoam reference (5000 its; ~3 h
   at np10). Start it only when no other 10-rank job runs. Then
   coupledFoam T5 coarse run (pytest test_T5 with CF_T5_MESH=coarse)
   under the D-042 criterion. The fine 5-8 M mesh is NOT for today
   unless everything else is green and the user has not said stop.
2. T-scaling on T4b: pytest tests/test_scaling.py --heavy (ranks
   1,2,4,8,10, 300 its each, both solvers). ~2-3 h. Only after T4b
   runs clean.
3. Benchmark: `bench/run_bench.py` configs A-D on T0,T1,T2 (fast trio)
   first, 3 repeats; then E,F,G,H per their B10 scoping (E/H on T2 at
   least); T4b single repeat of A and C if the machine day allows.
4. After every completed stage: `python bench/make_report.py`, then
   `make -C report/paper` (builds BOTH pdfs), copy both to the Windows
   scratchpad and send them to the user with SendUserFile.
5. Update AMENDMENT_B_STATUS.md rows that become "done" (B9 cycle
   study runs on the motorBike mesh once T4a mesh + defaults are
   stable: pytest tests/test_unit.py -k motorBike --heavy, and
   test_procAgglom).

## 5. Reporting duties (every stage)

- Every performance number as wall-clock AND CPU-hours (user's standing
  rule; rank wrapper + summary JSONs provide both).
- Commits per subsection, message style as tonight, attribution line
  per the session reminder. Push only green states.
- Anything not solved cleanly -> FABLE_REVIEW.md "Open" with evidence,
  as tonight.

## 6. Constraints and hard-learned gotchas (do not rediscover these)

- WSL calls: `MSYS_NO_PATHCONV=1 wsl.exe -e bash -lc '<cmd>'`; complex
  commands as script files in the Windows scratchpad, executed via
  `/home/jonas/bin/cfenv sys bash <script>`. Never `wsl.exe -- bash -c`
  with `$vars`.
- Files written through \\wsl.localhost become root-owned: run
  `MSYS_NO_PATHCONV=1 wsl.exe -e /home/jonas/bin/cf-fixown` after every
  UNC write.
- foamDictionary: ALWAYS `-disableFunctionEntries` when setting
  entries. It also quantises scalars to 6 significant digits on write.
- Build: -Wall -Wextra -Werror -Wfloat-conversion; FPE traps are bugs;
  constants in src/include/coupledDefaults.H; never write under
  /usr/lib/openfoam; stale-binary lesson: wclean the apps whenever a
  src/ header changes.
- 10 physical cores TOTAL for everything (user cap, D-031); SMT
  siblings are adjacent (cpu 0,1 = core 0); bind with
  `--bind-to core --map-by core`; `CF_MPI_BIND=none` for small tests.
  Never kill processes not identified by /proc/PID/cwd in your own run
  dirs.
- T0 Re1000 iteration counts are chaotic (2x scatter) - never tune or
  debug against them; use Re100 and the unit tests for regressions.
- Precision model D-001/D-012 (float blocks, double reductions, double
  FGMRES/GMRES iterates) is not negotiable; never loosen any test
  threshold except TASK 2's user-approved change.

## 7. Cleanup (end of day, optional)

- `git worktree remove --force /home/jonas/cf_precond` and
  `/home/jonas/cf_b11`; delete /home/jonas/cf_b11_base,
  /home/jonas/cf_*_platform dirs. FIRST copy
  /home/jonas/cf_precond/run/T1_base/linsys to run/dump_T1_linsys
  (needed by TASK 1) and /home/jonas/cf_precond/run/_variants (variants
  example).
- /home/jonas/cf_b11/applications/test/.rootjunk_tkb is root-owned;
  needs the user (sudo). Harmless, can stay.
- run/exp_re1000 and run/exp_linsolver must STAY
  (bench/exploratory_numbers.py reads them). Other run/exp_*, chaos_*,
  dbg* dirs may be deleted once their DECISIONS entries are written.

## TASK 5 - Deep diagnostics logging (user request, 2026-09-22 ~08:15)

User (translated): "I want the option of a mega-verbose output in the
solver log so every single step can be logged, including all sub-steps
and computations, so that you can analyse everything later and tune
parameters depending on the state the simulation is in - at the start
or end of the run, and depending on current convergence."

Priority: after TASK 1 and TASK 2 are green; implement it while the
long TASK-4 background runs (T5 reference, scaling) occupy the machine.
Exception: the per-patch BC-flip counter (5.4 level 3, item f) may be
pulled forward if the TASK-3 probes fail, since it directly tests the
freestream-chatter hypothesis.

### 5.1 Control and files

- New dictionary block `coupled.diagnostics`:
  `level` 0|1|2|3 (default 0), `echo` yes|no (default no),
  `maxBytes` (default 2 GiB per rank, new constant
  coupledDefaults::diagMaxBytes; when exceeded, write one final
  truncation record and stop logging, never abort the run).
- Output: `<case>/diagnostics/diag.rank<N>.jsonl` - one self-contained
  JSON object per line (JSON Lines), full precision (max_digits10),
  flushed once per outer iteration. The human log (log.coupledFoam and
  its CF| lines) stays EXACTLY as it is; `echo yes` additionally prints
  the level-1 record pretty-printed into the log.
- Implementation: new class `src/io/diagnostics.{H,C}` following the
  jsonWriter/runInfo patterns; modules receive a pointer/reference
  (null object pattern at level 0). HARD RULE: at level 0 no strings
  are built, no fields gathered, no collectives added - guard at the
  call site (`if (diag.active(2)) { ... }`). All new keywords/constants
  in coupledDefaults.H.

### 5.2 Phase classification (logged at every level >= 1)

Field `phase`, computed each outer iteration from existing quantities,
exact rules in this order:
1. "startup"    if beta < 1 (deferred-correction ramp / upwind phase)
2. "stalled"    if R has not decreased in the last 50 iterations
                (compare against min(R) of that window)
3. "asymptotic" if R <= 100 * convergence.residualTol (when
                residualTol > 0) OR the Cd/Cl window criterion is
                within 10x of its tolerance
4. "ramp"       otherwise
This field is the hook for later state-dependent parameter tuning
(gamgAutoTune or a successor will consume it; nothing auto-tunes in
this task - autoTune stays off per D-039).

### 5.3 Record contents, level 1 (one record per outer iteration)

All of: iter, phase, wallTime, and
- residuals: R, rU, rp, R1, normFactor, linear initial/final residual,
  linear its, converged flag, rho; continuity: max and sum |cell mass
  error| from the flux update.
- controls: CFL (value + the local-dt distribution min/median/max and
  count at the local limit), PTC strategy + hold counter + growth
  factor applied this iteration, eta with the raw Eisenstat-Walker
  value and which safeguard clipped it, omega and line-search cuts +
  trial list summary, sentinel checks/rollbacks, remediation
  nStat/nDyn + version counter, Anderson status/m/gamma norm/flush
  reason (if any), beta.
- turbulence: initial/final residuals of the k and omega solves,
  bounding event counts (k, omega), nNutCapped, nClamped.
- timings: tAsm split (stage1 momentum ops, stage2 continuity/scaling,
  boundary coupling, Rhie-Chow), tSolve split (preconditioner setup,
  preconditioner applications, Krylov vector ops), tTurb, tDiag (the
  cost of the diagnostics themselves), tIter, tWall.
- GAMG: nLevels, cellsPerLevel, ranksPerLevel, C_op, and a
  reagglomeration/setup event flag - full hierarchy only when it
  changed, else `hierarchy: "unchanged"`.
- memory: current and peak RSS of this rank (getrusage, cheap).

### 5.4 Additional content, levels 2 and 3

Level 2 (per linear solve):
- a. per-Krylov-iteration preconditioned residual array; at every
  FGMRES restart the true (double) residual as well.
- b. per preconditioner application: per-level residual norm before
  and after smoothing (down and up legs), coarse-level solver its and
  final residual, K-cycle inner GCR steps and acceptance test values.
- c. rhoOpt (the scale-free measure from Test-blockSystem/D-039) for
  the FIRST application of each solve.
Level 3 (per sub-operation; expensive, for short diagnostic runs):
- d. per smoother sweep residual norms; per-level operator stats on
  setup (rows, off-diagonal dominance min/median).
- e. per line-search trial: alpha, number of physicality violations by
  type.
- f. per-patch BC state counters: for every mixed/freestream patch the
  number of faces that switched inflow/outflow since the previous
  iteration (this is the T3 chatter probe).
- g. Anderson internals: QR condition estimate, per-column |gamma|,
  maxAlpha clip events.
Document per level which extra parallel reductions it introduces;
level 1 must add none beyond values already reduced today.

### 5.5 Analysis companion + test

- `bench/diag_tools.py`: `load(case) -> pandas.DataFrame` (level-1
  records; per-rank files merged, master preferred for global fields),
  `linear_history(case, iter)` for level-2 arrays, and three canned
  matplotlib figures: (1) R/CFL/eta/rho vs iteration with phase
  bands, (2) per-level GAMG residual-reduction heatmap over
  iterations, (3) stacked time breakdown. Vector PDF output alongside
  the report figures.
- New quick test `tests/test_diagnostics.py` (NOT heavy): T0 Re100,
  maxIter 20, level 3, np1 and np4: every line of every
  diag.rank*.jsonl parses as JSON; the level-1 required keys above are
  present; tDiag is recorded; with level 0 the diagnostics/ directory
  is not created.

### 5.6 Acceptance

- T0 Re100 np1 wall time: level 0 within +-5 % of the pre-TASK-5
  build (three runs each, compare medians); level 1 overhead <= 2 %;
  level 2 overhead <= 15 % (document the measured numbers).
- Unit battery and pytest T0 unchanged (level 0 default everywhere; no
  case template gets a diagnostics block by default).
- DECISIONS entry D-045: design, field list reference, overhead
  numbers, file sizes per level measured on T0 (level 3) and T1
  (level 2). README section "Diagnostics" with two usage examples
  (turn it on for a case; load and plot in python).

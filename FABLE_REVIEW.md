# Open problems for a follow-up review (Fable)

Started 2026-09-21 ~23:10 at the user's request: everything that causes
problems overnight and is not solved cleanly is listed here with evidence
and a proposed fix, so a stronger model can pick it up in the morning.
Resolved items move to the "Resolved" section with the commit.

## Open

### 1. RESOLVED 2026-09-22 ~04:30 - trajectory chaos, not a code defect

Perturbation experiment (run/chaos_T0_*, committed build, exact
orthogonality test in both runs): T0 Re 1000 np1 with nu perturbed by
1e-5 relative needs 122 (Re 999.99) and 108 (Re 1000.01) outer
iterations instead of 56; unperturbed reruns are bit-identical at 56.
So the outer iteration count at Re 1000 is chaotically sensitive to
round-off-sized perturbations near the start-up switch, and the
orthogonality-flag A/B (58 vs 106) was trajectory scatter with n=1 per
cell, not a causal mechanism. The same explains the np4 decomposition
scatter (182-282).

Consequences:
- The exact test stays (validated baseline); the tolerance would be
  equally valid statistically, but there is no reason to re-validate.
- Iteration counts on T0 Re 1000 (and their wall times) carry a ~2x
  trajectory scatter; the paper must not present single Re-1000
  iteration counts as robust. Re 100 is deterministic and insensitive.
- No further debugging warranted.

### 1-old. Round-off-sized change doubles T0 Re 1000 iterations (np1)

Evidence (debug agent, 2026-09-21, A/B on T0 Re 1000 np1):

| orthogonality test in nonOrthCorrection | FGMRES | outer iterations |
|---|---|---|
| exact `kMax == 0` (committed) | old | 58 |
| exact | mixed precision (committed) | 57 |
| `kMax <= 1e-10` | old | 123 (5 CFL cuts) |
| `kMax <= 1e-10` | mixed precision | 106 (6 CFL cuts) |

The T0 blockMesh cavity has max|k| = 7.1e-14 (round-off), so the only
difference is whether the non-orthogonal correction path runs with
round-off-sized vectors or is skipped. A correction of 1e-14 cannot
change the physics. Most likely `orthogonal_` also switches something
that is not just the k-correction: a different delta-coefficient field,
a different gradient or explicit-source path, or the deferred correction.
It may also be that Re 1000 is extremely sensitive at the start-up switch
(beta 0 -> 1, D-008).

Related: T0 Re 1000 at np4 needs 182-282 outer iterations (np1: 57),
with or without the tolerance. A coupled Newton-like solver should be
nearly independent of the decomposition; only the preconditioner
differs. This suggests the outer convergence at Re 1000 depends on
inner-solve quality (Eisenstat-Walker eta, K-cycle at np4) more than
expected.

To do:
- (a) Diff the code paths guarded by `orthogonal_` / `nonOrthCorrection`
  and find what else the flag controls.
- (b) Log per-iteration R, CFL, eta and cuts for np1 vs np4 at Re 1000
  and find where they diverge (likely around the start-up switch).
- (c) Decide the tolerance: the orthogonality tolerance is conceptually
  right (coupledDefaults::orthogonalityTolerance exists but is unused).
Committed state: exact test (f025546).

### 2. coupledFoam is ~5x SLOWER than simpleFoam on a turbulent case (T1)

STATUS 2026-09-22 (branch precond-research, D-039): fixed for T1 by the
preconditioner study. With agglomerationWeights combined (pair
agglomeration on block-matrix weights), nFinestSweeps 1 (two ILU0 sweeps
amplify: the smoother is a divergent Richardson iteration on this
saddle-point system), nCellsInCoarsestLevel 20, autoTune no and
reagglomerateInterval 50, T1 needs 495 outer iterations in 33 s instead of
518 in 175 s, against 36 s for simpleFoam (np 1, same load; np 2: 31 s
wall, 60 s CPU). Mean linear iterations 25.4 -> 3.3. The text below is the
original problem statement; what remains open is listed in D-039 (the rho
of 6.3.5 is unscaled and misleads autoTune; the startup phase is still the
most expensive part; T2-T5 and the F1 cases have not been re-measured).

This is the most important open problem. It blocks the project's goal
(the fastest solver) on every realistic case, the F1 car included.

Evidence (T1 study, D-038; logs in run/exp_T1_*, summary in
run/exp_T1_results.txt):

| run | outer its | wall |
|---|---|---|
| coupledFoam, converging configuration: ILU0, FGMRES restart 30, K cycle, upwind k/omega, k/omega relax 0.95 | 483 to R < 1e-6 | 364 s (np1) |
| simpleFoam, same schemes | 769 to 1e-8 | 72 s |

Per outer iteration coupledFoam costs 0.35-0.7 s, simpleFoam 0.09 s.
Almost all of that is the linear solve:
- rho (the first preconditioner application) is often > 1.
- A 0.5 residual reduction needs 20-60 FGMRES iterations.
- The block-GAMG K cycle cannot solve the upwind CFL-500 system at all.
- T0 (uniform cavity) shows the opposite (9.6 s vs 173 s), so the
  weakness is specific to graded, high-aspect-ratio meshes and turbulent
  viscosity contrasts.

Suspects / directions:
- (a) Agglomeration: faceAreaPair on the 4x4 system uses geometric
  weights only. The p-p Rhie-Chow block and the strongly anisotropic
  momentum coupling across thin boundary-layer cells probably need
  coefficient-based weights, e.g. the momentum diagonal or
  |a_PN|/max(a_P,a_N).
- (b) The coarse operator: Galerkin summation of the saddle-point blocks
  loses the p-p stabilisation. Check the coarse-level p-row diagonal
  dominance.
- (c) The smoother: point-block ILU0 on the saddle-point system. Try a
  Vanka/SCGS block smoother, or a Schur-complement (SIMPLE-type)
  preconditioner as smoother or preconditioner.
- (d) Literature: Uroic & Jasak 2021 (CPC) report coupled speedups with
  block-selective AMG on OpenFOAM meshes; compare their agglomeration and
  smoother choices.

Also open from the T1 study:
- Whether T1 should stay on upwind k/omega (D-038). With the tutorial's
  limitedLinear k/omega neither solver converges.
- The motorBike, Ahmed and F1 cases may show the same missing steady
  state with second-order turbulence convection. Check each before
  blaming the solver.
- Every case's absolute linear tolerance must stay below
  R1*residualTol. A cheap code guard (a warning, or a tolerance relative
  to R1) is not implemented.

### 2b. Test status after the fixes (final build, 2026-09-22 02:27)

- T0 Re100/Re1000 at np1 and np4: all PASS.

  | run | its | wall |
  |---|---|---|
  | Re100 np1 | 56 | 17 s |
  | Re100 np4 | 56 | 22 s |
  | Re1000 np1 | 56 | 27 s |
  | Re1000 np4 | 214 | 213 s |

  The np4 growth at Re1000 is item 1.
- T1 np1 and np4: FAIL. The spec criterion is R < 1e-6 within 400
  iterations; the final R was 1.5e-5 at 400. The converging
  configuration needs about 483 iterations (item 2). The criterion was
  not loosened.

### 2c. T4b: the simpleFoam reference itself never meets criterion 12.3(ii)

T4b (motorBike, 1.70 M cells, 10 ranks): simpleFoam ran 4000 iterations
in 11876 s wall / 33.0 CPU-h (run/ref_T4b_np10). Cd and Cl keep
oscillating:
- Cd 0.393-0.405 (±1 %) and Cl 0.061-0.072 (±8 %) over iterations
  3500-4000.
- The last 100 iterations span 2.3 % in Cd and 12.5 % in Cl.

Criterion 12.3(ii) needs a 0.2 % window. The flow (a bluff-body wake) is
physically unsteady, so no steady solver can meet it: test_T4 fails on
"referenceConverged" whatever coupledFoam does. The same probably holds
for T5 Ahmed (a 25 degree slant is near the separation switch) and for
the F1 car.

Proposal (not implemented; changing the criterion is the user's call):
- Compare the mean Cd over the last N iterations (e.g. 1000) with its
  statistical uncertainty.
- Record "converged" as a stationary mean (e.g. a drift of the running
  mean below 0.2 % over 1000 iterations) rather than a 0.2 % min/max
  window.

T4b reference, mean over iterations 3500-4000: Cd 0.3993, Cl 0.0671.

### 3. Build hygiene: the installed binary was stale after the C++ fix commit

After 1171d44 the installed coupledFoam binary still referenced the old
rhieChow::updateFlux signature: wmake does not relink an application when
only the library changes. The C++ fix agent's T0 verification may
therefore have run a mixed build. Clean rebuild at 02:10 on 2026-09-22;
T0/T1 are being re-verified.
Recommendation: `./Allwmake` should always rebuild the applications after
the library (check the top-level Allwmake), or run wclean on the apps when
headers in src/ change.

### 4. T3 airFoil2D: SST limit cycle and marginal linear solves (TOP OPEN ITEM)

State 2026-09-22 morning, after the nut-cap fix (D-040) and the D-039
preconditioner settings:
- The nut cap no longer binds (nNutCapped 0) - that bug is fixed.
- T3-SST still limit-cycles: R oscillates 1e-3..1e-2, forces oscillate
  around roughly credible values without settling (old run: Cd swing
  -0.02..+0.45). k bounding fires every iteration; 80-240 cells sit
  permanently in the local CFL limiter.
- Probes with the new build (run/exp_T3_cfl20, exp_T3_cfl500): both
  aborted via B4 at ~434 / ~260 iterations because the linear solve
  reached only a 0.51 reduction vs the eta = 0.5 target in 200
  iterations - the D-039 settings are tuned on T1 and are marginal on
  the stretched T3 C-grid. Linear maxIter is now 400 on T3/T4/T5 as
  mitigation (commit after 83758a3).
- T3-GEKO converges loosely (R 2.3e-4) with Cd 3.9 % off the reference
  (tolerance 0.5 %).
Suspects for the limit cycle, in order: (a) PTC too aggressive for a
case simpleFoam needs 20000 heavily relaxed iterations for (T3-specific
CFLmax, or slower PTC growth after bounding events); (b) segregated
k/omega lag vs the coupled step (relaxation 0.7 -> 0.5, or turbulence
sub-iterations); (c) freestream mixed-BC switching chatter at the far
field; (d) preconditioner on stretched cells - rerun Test-blockSystem
on dumped T3 systems (the harness from D-039 makes this cheap).
Related: T1 needs ~495 its vs the 400-iteration test budget (spec 13);
the outer-loop convergence rate is the remaining lever now that
wall-clock per iteration is fixed.

## Resolved (2026-09-22 morning)

- Item 1: trajectory chaos, evidence in the section above.
- Item 2 wall-clock: D-039 merged - T1 33 s vs simpleFoam 36 s; with
  B11 (D-041) T0 Re100 np1 is at 6.2 s. Open remainder: outer iteration
  count (item 4), scale-free rho so autoTune can return.
- Items 2b: T0 np1+np4 Re100/Re1000 all PASS on the merged build.
- Item 3: clean rebuild done for the merged build; an Allwmake relink
  guard stays a nice-to-have.

### 4a. OPUS_TASKS section 0 result - T4b conservative run (2026-09-22 07:58)

run/T4b_conservative (V / geometric / nFinestSweeps 2 / coarsest 200 /
ILU0 / restart 30 / CFLmax 100, np10, 1.70 M cells) was stopped by the
lead at outer iteration 113 (not aborted): healthy up to ~iter 100 (2-17
linear its, 4-11 s/iter), then at CFL ~31 the solve degraded to 59
linear its and 33 s/iter while R stayed at ~1e-2 (0.0095-0.0177 over
iters 20-113). nLocLim 150-800 cells, rho 0.65-1.37. Extrapolated cost
for the 1500-iteration budget: ~13 h - not viable. Conclusion for
TASK 1: the conservative set also loses robustness as CFL grows on the
snappy mesh, so the TASK-1 grid must include the T4 regime at CFL >= 30
(T4a dumps at iters 20/40/50, CFL 500, cover it). Case dir kept.

### 4b. TASK 1 acceptance (D-043 + addendum + D-044, 2026-09-22 ~10:00)

Unit battery all PASS (Test-blockGAMG V/F/W/K np1+np4, blockMatrix,
blockFGMRES, block4Ops, doubleReduce, precision, kernelBandwidth gates).

| test | result |
|---|---|
| T0 Re100 np1 / np4 | PASS, 57 / 57 its, 5.7 / 5.7 s |
| T0 Re1000 np1 / np4 | PASS, 111 / 257 its, 26 / 71 s (chaotic count, item 1) |
| T1 np1 / np4 | FAIL (400-it budget): final R 2.7e-6 / 4.3e-6 vs 1e-6; dp within 0.057 % of simpleFoam; ~26 s for 400 its |
| T2 np1 | FAIL: no B4 abort any more (1000 its), final R 8.9e-6; xr within 0.68 % vs the 0.5 % tolerance |

Open (not loosened): the T1 400-iteration budget and the T2 0.5 % xr
tolerance are spec thresholds; both runs are physically close but the
outer-iteration convergence rate near R ~1e-5 is the limit. Also: the
test JSONs for T0-T2 carry wallSeconds but no cpuHours (the solver
summary JSON has them) - reporting gap to close in tests/cflib.

### 4c. TASK 3 result - T3-SST limit cycle persists (2026-09-22 ~10:10)

Step 1 (D-043/D-044 defaults, nut cap 1e8 effective, 3000 its, 164 s):
R oscillates 3e-3..8.6e-3 at CFL 500; 1300-2500 of 10720 cells locally
CFL-limited; bounding k in every iteration (3502 events); forces
oscillate around negative values. Step 2 probes (1500 its each,
run/exp_T3_*):

| probe | min R | last R | locally limited (mean last 300) | bounding k | Cd / Cl (mean last 500) |
|---|---|---|---|---|---|
| a CFLmax 20 | 1.5e-4 | 2.6e-3 | 134 | 1750 | 0.110 / 0.055 |
| b k/omega relax 0.5 | 1.3e-4 | 3.3e-3 | 1844 | 1313 | 0.127 / 0.064 |
| c startupUpwindIters 200 | 1.0e-3 | 8.6e-3 | 1807 | 1719 | 0.127 / 0.064 |
| simpleFoam reference | | | | | 0.0907 / 0.253 |

None converges (TASK 3 rule 3: no template change, move on). Key
observation for the next analysis: Cl is ~1/4 of the reference in ALL
probes - a systematic difference, not only a convergence problem
(check: freestream BC linearisation incl. the A3 (1-viC) weighting on
freestreamPressure, the far-field inflow/outflow switch - diagnostics
level 3 item f counts it -, the k bounding every iteration, and whether
the reference itself is converged: simpleFoam needed 20000 relaxed
iterations and the GEKO reference Cl is 0.838 vs SST 0.253).
Note: probe c did not change the switch point materially - the start-up
switch also fires on R < 1e-2 (D-008), so startupUpwindIters alone
cannot extend the upwind phase.

### 5. np4 runs are not bit-reproducible run to run (found by TASK 5, 2026-09-22)

T0 np4 on the main-install binary: three identical runs gave 7 / 9 / 9
linear iterations at outer iteration 1, with different rho. Serial runs
are bit-identical. Suspects: a reduction whose order depends on message
arrival (e.g. a non-deterministic gather/sum in the processor
agglomeration or the coarsest-level gather), or uninitialised data on
interface buffers. D-034 states bit-exact np4 restart round trips; that
was measured within one run pair and may only hold by chance. Relevant
for the chaotic T0-Re1000 np4 counts (item 1). To do: run T0 Re100 np4
three times with -procAgglom none (Test-blockGAMG / processorAgglomerator
none) and with GAMG disabled (blockDiagonal) to bisect.
TASK 5 (branch task5-diagnostics, commits 28fc54b eaedf50 4085f4d)
waits for a quiet window: `/home/jonas/bin/cfenv sys bash
bench/task5_acceptance.sh` (np4 checks + wall-clock overhead gates),
then merge.

### 4d. T3-SST: the reference is steady, only coupledFoam oscillates (2026-09-22 ~12:10)

Read-only D-042 evaluation (window W = n/2): simpleFoam SST reference
(run/ref_T3_kOmegaSST, 20000 its) is perfectly steady - Cd 0.0907,
Cl 0.2527, std ~0 - i.e. strongly separated (Cl vs GEKO 0.84) but NOT
oscillating. coupledFoam SST oscillates massively: task3_T3sst (3000
its) Cd 0.125 +- 0.113, Cl 0.344 +- 0.678; probes a/b/c similar. So the
T3-SST oscillation is a coupledFoam limit cycle, not physics; a
window-mean criterion would hide it (and still fails: Cd 0.125 vs
0.091). A mean criterion for T3-SST was briefly applied on a wrong
premise and reverted before commit. Next steps as in 4c (freestream
BC switching, k bounding every iteration, level-3 diagnostics patch
switch counter now available on main).

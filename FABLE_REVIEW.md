# Open problems for a follow-up review (Fable)

Started 2026-09-21 ~23:10 at the user's request: everything that causes
problems overnight and is not solved cleanly is listed here with evidence
and a proposed fix, so a stronger model can pick it up in the morning.
Resolved items move to the "Resolved" section with the commit.

## Open

### 1. Round-off-sized change doubles T0 Re 1000 iterations (np1)

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

### 3. Build hygiene: the installed binary was stale after the C++ fix commit

After 1171d44 the installed coupledFoam binary still referenced the old
rhieChow::updateFlux signature: wmake does not relink an application when
only the library changes. The C++ fix agent's T0 verification may
therefore have run a mixed build. Clean rebuild at 02:10 on 2026-09-22;
T0/T1 are being re-verified.
Recommendation: `./Allwmake` should always rebuild the applications after
the library (check the top-level Allwmake), or run wclean on the apps when
headers in src/ change.

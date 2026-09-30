# Changelog

All notable changes. Code marked *uncompiled* has only passed
`scripts/syntax-check.sh` (g++ -fsyntax-only) and has not been built with
wmake or run yet.

## [Unreleased]

### F1 robustness - 2026-09-30 (D-074..D-078, F1_RUN_FINDINGS.md)
The first coupledFoam run on the 20.6 M-cell F1 half-car died between
iterations 14 and 18 in every configuration tried. Five solver changes:
- D-075 `ptc.continuityFactor` (default 1): pseudo-time term
  c V/(dtau Uref^2) on the p-p diagonal of the increment system, in the
  static and dynamic remediation cells. Before, a CFL cut made the
  pressure step grow like 1/dtau. Converged solutions unchanged (the term
  is added after the residual). A first version on all cells slowed the
  global pressure modes (T0 Re 100: 68 -> >300 iterations) and was
  withdrawn the same day.
- D-076 a sentinel rollback puts the offending cells into the dynamic
  remediation set immediately (`remediation::activatePending`), one more
  halo layer per consecutive rollback. Before, they were only activated
  after an accepted step, i.e. never in a rollback cascade.
- D-077 `potentialClip` changes U only; new `potentialPressure keep |
  bernoulli` (default keep) sets p = H - |U|^2/2 after the clip. The old
  head-keeping p shift created isolated spikes on the uniform p of
  potentialFoam.
- D-074 `Uref boundary` counts prescribed boundary velocities only (F1:
  50.5 instead of 321.2 from slip/symmetry/inletOutlet copies).
- D-079 `lineSearch.mode local` (default global): full step when only a
  few cells (<= localFraction) violate the step bounds; those cells enter
  the dynamic set, no global CFL cut. `nLocStep` on the CF| line.
- D-080 `solvers.coupled.stagnationRatio` (default 0 = off): FGMRES stops
  unconverged when a restart cycle does not reduce the residual enough.
- D-081 `ptc.failCeiling`, `ptc.ceilingRelax` (default off): after a failed
  linear solve the CFL stays below a ceiling that relaxes slowly. A sentinel rollback sets it too.
- D-079 second version: omega as the localFraction quantile (0.05 %).
- D-082 the potentialClip cells start in the dynamic set.
- D-083 `ptc.cflStrategy EFF`: steers the share of the linear solve in the
  iteration time to `ptc.effSolveShare` (0.6), decided on a window of
  `ptc.effWindow` (4) iterations, step from `ptc.effExponent` (1.6), band
  `ptc.effTol`; `cflWhy=`, `eff=`, `solveShare=`, `ceil=` on the CF| line.
- Test-nonOrthLimiter (unit test of D-078).
- D-078 overflow-free limiter in `nonOrthCorrection::correctionFromGrad`
  (bit-identical where the old form did not overflow; it raised SIGFPE on
  steep pressure gradients).

### Phase 0 / A - 2026-09-21
- Repository skeleton (spec Section 4), README, PLAN, DECISIONS (D-001..D-014),
  spec copy (`SPEC_coupledFoam.html` authoritative, `.md` text extraction).
- `src/blockMatrix`: `blockScalar.H` (float storage / double reductions,
  D-001), `block4Ops.H` (4x4 kernels, LU inversion with pivot guard),
  `doubleReduce` (double-accumulated reductions, MPI_DOUBLE),
  `blockLduInterface` (non-blocking processor exchange of 4-vectors),
  `blockLduMatrix4` (4x4-block LDU storage, Amul, row scaling). *uncompiled*
- `scripts/syntax-check.sh`: syntax-only compile with project warning flags.

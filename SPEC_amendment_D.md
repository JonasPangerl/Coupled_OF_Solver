# Amendment set D — full single precision (SP) as a build option

Received from the user on 2026-09-22 (verbatim below the lead notes).

## Lead notes (deviations and interpretation for this project)

- **Reference precision is DP, not SPDP.** This project builds on the
  read-only system OpenFOAM v2606 **DP** install (D-001: the block linear
  solver is already float internally, reductions in double). The user asked
  for SP "tested against the maximum precision we have", i.e. **SP vs DP**.
  Where the text below says SPDP as the reference, read DP. An SPDP build is
  not made for now (it would be a third OpenFOAM build; it can be added
  later without code changes, D2 makes the source precision-agnostic).
- **Builds:** the system DP-Opt install (unchanged) and a private
  `~/OpenFOAM-v2606-SP` built from the local `openfoam2606-source` package
  with `WM_PRECISION_OPTION=SP`. **Compiler flags identical to the system
  build** (`c++OPT = -O3`, the unmodified linux64Gcc rules), not
  `-march=native`: the system DP build cannot be rebuilt, and SP-vs-DP
  timings are only fair with identical flags (D1 deviation).
- **No SP-Debug build.** FPE trapping (`FOAM_SIGFPE`) works in Opt builds,
  as all project tests already use it (D1 deviation).
- Mesh generation stays in DP (D5.1 with DP in place of SPDP).
- User remarks with the amendment: D6 (double-accumulated forces) is not
  optional — without it SP measures rounding noise instead of convergence on
  the large case, and the comparison with simpleFoam would be unfair. The D5
  origin shift only works if CofR, probes and sampling lines are shifted
  too; that is the typical place where silently wrong moments appear.

---

## Amendment set D — full single precision (SP) as a build option

### D0 — Purpose and scope
Adds `WM_PRECISION_OPTION=SP` (scalar = solveScalar = float) as a supported
build mode alongside SPDP (default) and DP. Goal: minimum wall-clock and
memory. Precision is a build-time choice only; there is no runtime switch.
The same source tree must compile warning-free and pass its gates in DP,
SPDP and SP. Everything not mentioned here is unchanged.

### D1 — Builds (extends Section 3)
Four OpenFOAM v2606 builds side by side:

| Build | WM_PRECISION_OPTION | WM_COMPILE_OPTION | Use |
|---|---|---|---|
| OpenFOAM-v2606-SPDP-Opt | SPDP | Opt | default production + benchmarks |
| OpenFOAM-v2606-SPDP-Debug | SPDP | Debug | FPE test runs |
| OpenFOAM-v2606-SP-Opt | SP | Opt | SP benchmarks |
| OpenFOAM-v2606-SP-Debug | SP | Debug | SP FPE test runs |

All four are compiled with identical optimisation flags. Locate the C++
optimisation rule file of the linux64Gcc rules (`find wmake/rules -name
'c++Opt'`) and set `c++OPT = -O3 -march=native -fno-math-errno` in every
Opt build **before** compiling OpenFOAM. Record the flags in
`results/gates/phase_0.json`. `simpleFoam` from each build is the native
reference for that precision (no cross-precision native comparisons in
speedup numbers).
`Test-precision` extended: in SP, assert sizeof(scalar)==4 and
sizeof(solveScalar)==4; in SPDP 8/4; in DP 8/8.

### D2 — Source-code type rules (mandatory, all of src/)
1. Field data and matrix coefficients use `Foam::scalar` / `Foam::solveScalar`
   exactly as specified elsewhere; never `float` or `double` for them.
2. The following use the literal C++ type `double`, regardless of build
   (these are the complete list; nothing else may use `double`):
   a. all global and local reductions: norms, dot products, sums (D4)
   b. FGMRES Hessenberg matrix, Givens rotations, small least-squares
   c. GCR coefficients in the K-cycle (a1, a2, orthogonalisation)
   d. Anderson QR / Gram system and alpha coefficients
   e. coarsest-level dense LU (6.3.3)
   f. 3×3 inversion / SVD for tensorial Rhie–Chow (C2); result stored as
      scalar after inversion
   g. 4×4 block inversion in block4Ops: computed in double, stored as
      solveScalar
   h. Galerkin coarse-operator summation (6.3)
   i. force and moment integration (D6)
   j. residual-history quantities R_n, R_{n-1}, eta_n, CFL, omega
      (scalars in the outer loop, not fields)
3. When `std::is_same<scalar, solveScalar>::value` is true (SP and DP),
   the assembler writes directly into the matrix arrays with no
   conversion pass and no staging buffer (`if constexpr`). The
   `scalar → solveScalar` clamp counter (9.2 item 10) is then evaluated
   in the same pass.
4. MPI reductions of double quantities use `returnReduce<double>` /
   `MPI_DOUBLE`; never reduce a double through a scalar-typed helper.

### D3 — Guard constants (replaces use of SMALL/VSMALL in src/)
In SP, OpenFOAM's `SMALL` is 1e-6 and `VSMALL` 1e-37; using them as
stabilisers changes behaviour between builds. Therefore:
- Create `src/blockMatrix/coupledConstants.H` with
  `template<class T> constexpr T cfVSmall()` = 1e-30 for float, 1e-300 for
  double; `cfSmallRel<T>()` = 1e-6 for float, 1e-14 for double;
  `cfGreat<T>()` = 1e30 for float, 1e300 for double.
- Every guard in 9.2 uses `cfVSmall<T>()` for denominators and
  `cfGreat<T>()` for clamps, where T is the type of the quantity being
  guarded (double for the D2.2 list).
- Gate check: `grep -nE '\b(SMALL|VSMALL|ROOTVSMALL|GREAT|VGREAT)\b'`
  over src/ must return no matches (applications/test/ excluded).
- Pivot guard (6.3) and pseudo-inverse threshold (C2) stay as written;
  they are evaluated in double per D2.2 f/g.

### D4 — Double-accumulation list (exhaustive)
Implemented once in `doubleReduce` and used everywhere:
`dot(a,b)`, `norm2(a)`, `sumMag(a)`, `sum(a)`, `sumSqr(a)`, weighted
volume sum `Σ V_P a_P`, and max/min are plain (no accumulation issue).
Inputs may be float or double arrays; the accumulator is always double;
loops keep `#pragma omp simd reduction(+:acc)` with a double accumulator.
Users: normFactor (6.1), all Krylov/GCR inner products, rho measurement
(6.3.5), R_n (7.1), line-search reference values, sentinel statistics,
continuity error print, D6 forces. Native OpenFOAM `gSum`/`gSumMag` on
scalar fields must not be used in src/ (grep gate as in D3).

### D5 — Geometry safeguards for SP
1. **Mesh generation and quality in SPDP.** blockMesh/snappyHexMesh/
   cartesianMesh run in the SPDP-Opt build. The SP case is created by
   copying `constant/polyMesh` after writing it with
   `writeFormat ascii; writePrecision 12;` (ASCII avoids binary
   scalar-width mismatches between builds).
2. **Origin shift.** Before any SP run, Allrun executes (SP build)
   `transformPoints -translate '(-cx -cy -cz)'` where (cx,cy,cz) is the
   bounding-box centre printed by `checkMesh` (parse "Overall domain
   bounding box"). The same shift is applied to CofR in D6 and to any
   point-based settings (probes, sampling lines). Record the shift in
   `constant/meshShift` (dictionary with the vector) so results can be
   mapped back.
3. **checkMesh gate.** Run `checkMesh -allGeometry -allTopology` in the SP
   build after the shift. Compare with the SPDP result: the SP run must
   not report any failed check that SPDP does not report, and the count
   of negative-volume / wrong-orientation cells must be 0. Otherwise the
   case is marked `SP-geometry-fail` in results and the SP benchmark for
   that case is skipped (not a solver failure).
4. Wall distance: native `meshWave` in the SP build; no change.
5. The static remediation set (8.1/C1) is recomputed in the SP build;
   its size difference to SPDP is logged and reported (D11).

### D6 — Double-accumulated force coefficients (new, all builds)
Native `forceCoeffs` accumulates in scalar (float in SP) over all wall
faces; with pressure-side cancellation this is too noisy for the 0.2 %
convergence window. Therefore:
- New class `src/io/coupledForces.{H,C}`, evaluated inside coupledFoam
  after every outer iteration, configured by `system/coupledForcesDict`:

      patches (body);          // wordRes
      rhoInf 1.225; magUInf 40; lRef 1.044; Aref 0.056;
      CofR (0 0 0);            // in the unshifted frame; shifted by D5.2
      dragDir (1 0 0); liftDir (0 0 1); pitchAxis (0 1 0);

- Per face: pressure force p_f S_f (kinematic p × rhoInf) plus viscous
  force −(devReff·S_f)·rhoInf with devReff from
  `turbulence->devReff()` boundary field; moment (x_f − CofR) × F_f.
  All products formed in double from the float face values; sums per
  rank in double; `returnReduce<double>`.
- Output: `postProcessing/coupledForces/<startTime>/coeffs.dat` with
  columns `iter Cd Cl Cm Cd_p Cd_v Cl_p Cl_v` (%.10e).
- 12.3(ii) and the benchmark harness (14.4) use coupledForces in all
  builds, and `simpleFoam` runs get the same values through a
  function-object wrapper of the same class (`libs (libcoupledFoam);
  type coupledForcesFO;`) so both solvers are measured identically.
- Unit test `Test-coupledForces`: on T3 airFoil2D (SPDP), coupledForces
  equals native forceCoeffs to 1e-6 relative.

### D7 — Precision profile defaults
`coupled.precisionProfile auto` (default): chooses defaults from
sizeof(scalar). Explicitly set keywords always win.

| Keyword | SPDP / DP | SP |
|---|---|---|
| `convergence.residualTol` | 1e-6 | 1e-5 |
| linear `tolerance` (absolute floor) | 1e-8 | 1e-6 |
| `etaMin` | 1e-3 | 1e-2 |
| `anderson` allowed above 35 M cells | no | yes (m ≤ 4) |
| `bounds.omegaMin` | 1e-6 | 1e-5 |
| `bounds.kMin` | 1e-12 | 1e-10 |

Rationale: in SP the field increment underflows relative to the field
value at ~1e-6 relative; tighter targets only produce stagnation.
The force-window criterion (12.3 ii) is unchanged and remains the
decisive convergence test.

### D8 — SP-specific performance measures
1. Zero-copy assembly (D2.3) — no conversion pass.
2. FTZ/DAZ on (9.1) — mandatory in SP Opt; denormals are far more
   frequent in float and cost up to 100× per operation.
3. Binary I/O (`writeFormat binary`) — half the bytes; restart only
   within the same build (see D5.1 for cross-build transfer).
4. No other SP-specific code paths. Speed comes from halved memory
   traffic in assembly, turbulence and I/O; kernel rules (B11) apply
   unchanged.

### D9 — Memory (45 M cells, SP)

| Item | ~GB |
|---|---|
| Mesh, fields, addressing (float) | 22–25 |
| Block matrix, GAMG hierarchy, Krylov, K-cycle (unchanged) | 55–65 |
| Turbulence, rollback copies, reserve | ~6 |
| **Total** | **~85–95** |
| Anderson m=4 (float: 128 B/cell) | +5.8 |

### D10 — Tests and gates for SP (new Phase F, after Phase E)
All run with the SP builds; SP-Debug with traps for T0, T1, T-fpe.

| Test | Pass criterion |
|---|---|
| Unit tests 6.4 + B9 + D6 | same criteria as SPDP, float tolerances unchanged |
| D3 grep gate | zero matches |
| T0 (Re 100, 1000) | u/v centreline L2 difference to SPDP < 1e-3 |
| T1, T3-SST, T3-GEKO | Δp / Cl / Cd within 0.3 % of SPDP; converged by 12.3(ii) |
| T-fpe | no trap in SP-Debug; recovers and converges |
| T-restart (T1) | as in 13, within the same build |
| T4b, T5 | D5.3 checkMesh gate; Cd within 0.5 % of SPDP; converged by 12.3(ii) |

A failure of D5.3 skips the case (reported); any other failure fails
Phase F. Phase F failing does not invalidate Phases 0–E (SP is an option).

### D11 — Benchmark and report additions
- Benchmark matrix (14.2) gains per case: (F1) simpleFoam SP-Opt,
  (F2) coupledFoam SP-Opt, both 3 runs, same protocol; T4b and T5
  mandatory, T1/T3 optional.
- Metrics added: `precision`, `meshShift`, `checkMeshDiff` (list of
  checks differing from SPDP), `staticSetSizeDiff`, `Cd_rel_to_SPDP`.
- Report new Section 15.10 "Single precision": table of wall-clock,
  time/iteration breakdown (assembly / linear solve / turbulence /
  I/O), peak RSS, and Cd deviation for SPDP vs SP for both solvers;
  convergence-floor plot (R_n SP vs SPDP); explicit verdict per case:
  "SP usable" if Cd within 0.5 % and checkMesh gate passed, else
  "SP not usable" with the reason.

### D12 — DECISIONS.md entries
- Mesh generation always in SPDP; SP is a solve-only mode.
- Expected SP gain is limited to assembly, turbulence and I/O (~15–20 %
  wall-clock); the linear solver is already single precision in SPDP.
- Known limitation: nonlinear convergence floor ~1e-5 normalized residual
  in SP; engineering convergence is judged by the force window only.

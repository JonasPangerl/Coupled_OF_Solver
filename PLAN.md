# coupledFoam - Implementation Plan

This plan implements `SPEC_coupledFoam.md` (v2). The spec's phase/gate structure
(Section 16) is kept unchanged. One phase that the spec does not contain is
added in front of it: **Phase 00, OpenFOAM source builds**. The spec assumes
the builds exist; on this machine they do not.

Deviations and spec defects are recorded in `DECISIONS.md` and referenced here
as `D-nnn`.

## Work order (user instruction, 2026-09-21)

The machine is shared with a long-running 31 M-cell F1 case. Therefore:

1. **Write all code first** (all phases), without `wmake` and without the
   OpenFOAM source build.
2. Allowed while the machine is busy: reading source, writing code, and single
   niced `g++ -fsyntax-only` checks of individual translation units against
   the installed (read-only) v2606 headers.
3. **Compile and run only when the machine is free**, checked immediately
   before every heavy step with

       wsl.exe -- bash -lc "uptime; pgrep -a -f 'snappyHexMesh|simpleFoam|foamRun|potentialFoam|pvbatch|mpirun' | head"

   Busy = load average >= ~16 or any listed process present. When free: every
   build and run uses `nice -n 19` and `WM_NCOMPPROCS=8`.
4. The gates are then passed strictly in order (spec rule 0.1).

## The precision question (D-001) - decision needed from the user

The spec's Phase 0 gate requires `sizeof(solveScalar)==4 && sizeof(scalar)==8`
and says to obtain this with `WM_PRECISION_OPTION=SPDP`. The v2606 source says
otherwise (`src/OpenFOAM/primitives/Scalar/scalar/scalarFwd.H`):

| WM_PRECISION_OPTION | scalar | solveScalar |
|---|---|---|
| DP   | double | double |
| SPDP | **float** | **double** |
| SP   | float  | float  |

No OpenFOAM build gives `scalar=double, solveScalar=float`. An SPDP build is
the opposite of the spec's intent: every field (U, p, phi, k, omega) would be
float and only the native linear solvers would run in double. The native
`simpleFoam` reference would then be a single-precision-field solver.

The spec's *intent* is unambiguous (5.2: "Fields U, p, phi stay scalar
(double)"; block storage float; reductions double). That intent is realised
independently of the OpenFOAM precision option by a library-owned type:

    // src/blockMatrix/blockScalar.H
    typedef floatScalar  blockScalar;     // block matrix storage and Krylov vectors
    typedef doubleScalar reduceScalar;    // all accumulations and MPI reductions

The code is written against `blockScalar`/`reduceScalar`/`scalar` and never
uses `solveScalar` for block data, so it is correct under DP, SPDP and SP.
On a DP build the float/double distinction is visible to the compiler
(`-Wfloat-conversion` is on), which removes the "compiles on DP, breaks on SPDP"
risk: mixing happens only through explicit, grep-able conversion helpers.

**Recommended:** DP builds. The Phase 0 gate is then evaluated as
`sizeof(blockScalar)==4 && sizeof(scalar)==8 && sizeof(reduceScalar)==8`, and
Test-precision additionally prints the literal `sizeof(solveScalar)` value and
records that the spec's literal criterion is unsatisfiable (D-001). The gate
itself is not loosened: it checks the same property the spec wants, on the
type that actually carries it.

## Phase 00 - OpenFOAM source builds (not in the spec)

Target layout (never under `/usr/lib/openfoam/`, never `apt` on `openfoam*`):

    ~/OF/OpenFOAM-v2606-DP-Opt     (only if the user wants no link to the system build, see below)
    ~/OF/OpenFOAM-v2606-DP-Debug   (required: Debug build for FPE-trapping test runs, spec 3)
    ~/OF/venv                      (python venv: pandas, pytest, numpy-stl)

Steps:

1. Ask the user to `sudo apt install bison` (only bison, no `openfoam*` package).
2. `cp -a /usr/lib/openfoam/openfoam2606/{Allwmake,applications,bin,etc,src,wmake,META-INFO}`
   into the target (read-only use of the system tree). `ThirdParty` and
   `tutorials` as needed (tutorials are read from the system tree).
3. Set `WM_PROJECT_DIR`-relative prefs via `etc/prefs.sh` in the copy:
   `WM_COMPILE_OPTION=Debug` (resp. `Opt`), `WM_PRECISION_OPTION=DP`,
   `WM_LABEL_SIZE=32`, `WM_MPLIB=SYSTEMOPENMPI`.
4. In a fresh shell that sources **only** that copy's `etc/bashrc`
   (under `set +u`, AGENTS.md pitfall), run
   `nice -n 19 env WM_NCOMPPROCS=8 ./Allwmake -s -l`.
5. Verify: `foamInstallationTest`, Test-precision, `simpleFoam` pitzDaily
   runs (Phase 0 gate items).

Estimated cost: roughly 2-3 h wall-clock per build at 8 parallel jobs; the
Debug build produces ~5 GB of binaries. Both builds only when the machine is free.

Open choice for the user: the Opt side can either (a) link against the
installed read-only system build `/usr/lib/openfoam/openfoam2606` (no write,
no reconfiguration; saves one build; its `simpleFoam` is the benchmark the user
already uses) or (b) use a private `~/OF/OpenFOAM-v2606-DP-Opt` copy. The
handover rule "Make/options must reference your build" points to (b), which is
the default of this plan. `Make/options` uses only `$(LIB_SRC)`/`$(FOAM_*)`,
so the choice is made purely by which `etc/bashrc` is sourced.

## Phase 0 - toolchain and skeleton (spec 16)

- Repository skeleton (spec 4), README, CHANGELOG, DECISIONS, Allwmake/Allwclean.
- `applications/test/Test-precision`.
- `Allrun-tests.sh` (CI script), `tests/test_env.py`.
- Gate: precision check (D-001), `simpleFoam` pitzDaily runs in both builds,
  `pytest tests/test_env.py` green -> `results/gates/phase_0.json`.

## Phase A - block linear algebra

`src/blockMatrix`: `blockScalar.H`, `block4Ops.H`, `doubleReduce.{H,C}`,
`blockLduMatrix4.{H,C}`, `blockLduInterface.{H,C}` (processor exchange of
4-vectors, native `processorLduInterface` tag/comm, native sign convention for
interface coefficients).
`src/blockSolvers`: `blockSolver` base + runtime selection, `blockBiCGStab`,
`blockDiagonal` preconditioner.
Unit tests: `Test-block4Ops`, `Test-blockMatrix`, `Test-doubleReduce` (1 and 4
ranks). Guard review: `grep -n '/' src/` output reviewed, every division tagged
`// GUARD:`. Gate -> `results/gates/phase_A.json`.

## Phase B1 - assembly and outer loop

`src/assembly`: momentum block from native `fvMatrix` coefficients (upwind
implicit + beta-weighted deferred HO correction, uncorrected Laplacian +
per-face-limited explicit non-orthogonal correction, explicit transpose
stress), pressure-gradient and divergence blocks from face geometry and the
native `valueInternalCoeffs/valueBoundaryCoeffs` of the patch fields (this is
what makes every BC in 5.5 implicit in a uniform way, D-002), Rhie-Chow p-p
block + source and the flux update from the **same** face loop (so that
`sum(phi)` equals the continuity-row residual exactly), row scaling, double-
precision residual `b - A x` computed before the float store.
Outer loop 5.7 with fixed CFL, sentinel 9.3, turbulence 5.8.
Gate: T0 (Re 100, 1000), T1 at fixed CFL 20, Debug build without traps.

## Phase B2 - block-GAMG

Native `GAMGAgglomeration::New` (faceAreaPair) with a per-attempt `name`
entry so that `mergeLevels` retries (C_op cap 6.3) create distinct
agglomerations. Block Galerkin summation in double, `faceFlipMap` swaps
upper/lower blocks, processor interfaces via native coarse `GAMGInterface`
face cells. Block Gauss-Seidel (Jacobi across processor faces), block ILU0,
coarsest-level block-BiCGStab. `blockGMRES`.
Gate: Test-blockGAMG, T1 blockGAMG < 0.5x blockDiagonal wall time, ILU0 decision.

## Phase C - PTC and adaptive control

`src/control/ptcControl` (local dt, mRDM/EXP/SER, hold), `lineSearch`,
local limit 7.3, start-up switch 7.4. Gate: T1, T2 incl. CFL criteria, T-fpe.

## Phase D - remediation, restart, MRF, remaining BCs

`remediation` (static set from native `polyMeshTools`, dynamic set with
hysteresis, `remediationFlag`), `coupledState` (IOdictionary, atomic write),
implicit MRF (D-006), inletOutlet/totalPressure. Gate: T3 SST + GEKO,
T-restart, T4a.

## Phase E - large cases, benchmark, report

T4b, T5 (Ahmed generator + snappyHexMesh), T-scaling, `bench/run_bench.py`,
`rank_wrapper.sh`, `parse_forces.py`, `make_report.py` -> `report/REPORT.md`.
Largest CPU consumers of the whole project; each run is announced with its
expected duration and only started when the machine is free.

## Writing order now (no compiler)

1. `src/blockMatrix`, `src/blockSolvers` (Phase A/B2 code)
2. `src/assembly`, `src/control`, `src/io`, `applications/coupledFoam`
3. `applications/test/*`
4. `cases/T0..T5`, `tests/` (pytest), `bench/`
5. `Make/*`, `Allwmake`, `Allwclean`, `Allrun-tests.sh`

After each block of files: one niced `g++ -fsyntax-only` per translation unit
against the installed DP headers (typos and signatures only), then a commit
`phase-<N>: <what>` with `[uncompiled]` in the body until the first real build.

# coupledFoam

Block-coupled pressure-velocity solver for steady, incompressible RANS in
OpenFOAM (ESI) v2606. The momentum and continuity equations are solved
simultaneously as one 4x4-block system per cell (u, v, w, p), with
single-precision block linear algebra, double-precision reductions,
block-GAMG-preconditioned Krylov solvers and pseudo-transient continuation.

The specification is `SPEC_coupledFoam.md` (v2). Deviations and spec defects
are recorded in `DECISIONS.md`, the work plan in `PLAN.md`.

**Status:** under development - see `CHANGELOG.md` and `results/gates/`.

## What this is (and is not)

coupledFoam is an **add-on**. It never modifies the OpenFOAM installation it
is built against:

- the library `libcoupledFoam.so` goes to `$FOAM_USER_LIBBIN`,
- the solver `coupledFoam` and the test applications go to `$FOAM_USER_APPBIN`,

i.e. the standard per-user locations (`~/OpenFOAM/<user>-v2606/platforms/...`).
Native solvers such as `simpleFoam` are unaffected; `coupledFoam` only runs
when called explicitly.

## Requirements

- OpenFOAM ESI v2606 (any `WM_PRECISION_OPTION`; DP is the reference, see
  `DECISIONS.md` D-001), `WM_LABEL_SIZE=32`
- GCC >= 13, OpenMPI >= 4.1
- Python >= 3.11 with numpy, pandas, matplotlib, pytest, pyyaml, numpy-stl
  (tests, benchmark harness, report)
- For the paper: a LaTeX distribution with `pdflatex`/`latexmk`

## Build

    source /path/to/openfoam2606/etc/bashrc     # exactly one OpenFOAM environment
    ./Allwmake -j 8                             # library, solver, test apps

`./Allwclean` removes all build products. Project code is compiled with
`-Wall -Wextra -Werror`; OpenFOAM headers are included via `-isystem`.

## Run

Every case in `cases/` has `Allrun` and `Allclean`:

    cd cases/T0_cavity && ./Allrun

Solver controls live in `system/fvSolution` (`solvers.coupled` for the linear
solver, top-level `coupled` for the outer iteration); every keyword and its
default is listed in spec Section 11 and printed at start-up under
`coupledFoam: effective settings`. Keywords added after the spec (all with
their defaults in `src/include/coupledDefaults.H`):

- `coupled.startupMode upwind | hybrid | none` (D-048, default `upwind`,
  the original behaviour): `hybrid` ramps the convection blending beta
  linearly from 0 to 1 over `startupRampLength` iterations from a ramp
  start chosen by the residual history (`startupRampStart`,
  `startupRampStartMax`, `startupSwitchR`, `startupStagnationFactor`,
  `startupStagnationWindow`, `startupFastFactor`); `none` starts with beta
  1. A fresh start without `potentialInit` on a non-uniform velocity field
  (a mapped solution) is probed at iteration 1 and keeps beta 1 if it is
  developed (`startupDevelopedTol`). The ramp state is part of the
  restart state.
- `coupled.Uref boundary | field | <value>` (D-050, default `boundary`):
  the velocity scale of the line search, the local CFL limit, the dynamic
  remediation set and the sentinel. `boundary` is the maximum over the
  non-coupled boundary values (inflow, moving walls, free stream);
  `field` is the previous definition (the maximum over the cells as
  well, i.e. the potential-flow peak with `potentialInit`).
- `solvers.coupled.blockGAMG.pivotGrowthLimit` (D-049, default 20, 0 =
  off): the blockILU0 smoother falls back to the unmodified diagonal block
  of a cell whose factorised pivot inverse grew by more than this factor.
  Fallbacks appear as `nPivFb=` on the `CF|` line and as `pivotFallbacks`
  in the summary.
- `coupled.ptc.linFailPolicy strict | reduction` (default `strict`):
  with `reduction` a linear solve that hit `maxIter` counts as a success if
  it is finite and reduced the residual to at most `linAcceptReduction`
  (0.9) times its start value (`linAcceptedUnconverged` in the summary).

An abort (B4 linear-solve failures or the sentinel rollback limit) writes
the last accepted fields, `sentinelFlag`, `remediationFlag` and the
remediation cell sets into `<iter>_lastValid` (per processor directory in
parallel), never into a numbered time directory.

## Tests and benchmarks

    python3 -m venv ~/OF/venv && ~/OF/venv/bin/pip install -r requirements.txt
    ~/OF/venv/bin/pytest tests/                 # writes results/tests/*.json
    ~/OF/venv/bin/python bench/run_bench.py     # writes results/bench/*.json
    ~/OF/venv/bin/python bench/make_report.py   # report/REPORT.md + paper figures

## Diagnostics

`coupled.diagnostics.level` 1-3 writes one JSON object per outer iteration
and rank to `<case>/diagnostics/diag.rank<N>.jsonl` (JSON Lines). The default
is level 0, which is off and adds no work. The `CF|` log lines do not change.

- Level 1: residuals, controls (CFL/PTC, local dt, Eisenstat-Walker,
  line-search trials, sentinel, remediation, Anderson), turbulence, timing
  split, GAMG hierarchy, memory, force window statistics, and the phase
  (`startup | stalled | asymptotic | ramp`).
- Level 2 adds per linear solve: the Krylov history, per-level GAMG norms,
  the K-cycle values and rhoOpt.
- Level 3 adds per-sweep norms, operator dominance, line-search violations,
  inflow/outflow flips on mixed patches, and Anderson internals.

Field list, cost and file sizes: DECISIONS.md D-045.

Turn it on for a case (always use `-disableFunctionEntries`):

    foamDictionary -disableFunctionEntries -entry coupled/diagnostics \
        -set "{ level 2; echo no; }" system/fvSolution
    ./Allrun

Load and plot in Python:

    import sys; sys.path.insert(0, "bench")
    import diag_tools as d
    df = d.load("run/T1_np1")                # level-1 DataFrame, index = iteration
    df[["residuals.R", "controls.CFL", "phase"]].tail()
    d.linear_history("run/T1_np1", 120)      # level-2 solves of iteration 120
    d.make_figures("run/T1_np1")             # 3 figures, PDF + PNG

    ~/OF/venv/bin/python bench/diag_tools.py run/T1_np1 [--report]

With forceCoeffs, the `CF|` line ends with `CdMean= CdRms= ClMean= ClRms=`
(and `CmMean= CmRms=`). These are the window mean and the RMS fluctuation.
An optional stationary-mean stop rule is available:
`coupled.convergence.forceCoeffsDriftTol` (default 0, off).

## Repository layout

    src/            libcoupledFoam (blockMatrix, blockSolvers, assembly, control, io)
    applications/   coupledFoam solver, unit-test applications
    cases/          test cases T0..T5 (Allrun/Allclean), shared case files
    tests/          pytest suite
    bench/          benchmark harness and report generator
    results/        generated JSON (only results/gates/ is versioned)
    report/         generated REPORT.md, LaTeX paper and figures
    scripts/        developer tools (syntax check, environment helpers)

## Report figures

`bench/make_report.py` writes every figure of the two LaTeX documents
(`report/paper/paper.pdf`, `report/paper/paper_tutorial.pdf`) into
`report/paper/figures/`; a figure whose data do not exist yet is shown as a
boxed "results pending" placeholder, so the papers build at every stage.

- **2D flow fields** (`bench/plot_fields2d.py`, matplotlib, vector PDF):
  for T0 (Re 100, 1000), T1, T2, T3 (SST, GEKO) the fields of coupledFoam
  and the simpleFoam reference on the same mesh, same colour scale, plus the
  difference field (`fields_<case>_U.pdf`, `fields_<case>_p.pdf`) and line
  plots (`profiles_<case>.pdf`: velocity profiles, T2 skin friction, T3
  surface C_p). Fields are read read-only with VTK's `vtkOpenFOAMReader`
  from the venv (`/home/jonas/OF/venv`, needs the `vtk` module); nothing is
  written under `run/`.
- **3D renders** (`bench/render_fields.py`, ParaView `pvbatch`, PNG): body
  surface mesh, |U| and C_p slices, surface C_p, and after the heavy runs
  the window-mean fields and the `coupledFieldCompare` delta fields of T4a,
  T4b, T5 (`render_<case>_*.png`). Decomposed cases are read in place.
  Camera, planes and colour ranges are fixed in the script. A case whose
  logs changed in the last 10 minutes is skipped (live run);
  `--reference-only` renders only the simpleFoam side.
- **Speed-up** (`make_report.py:fig_speed`): residual vs wall-clock time,
  wall time and CPU-hours to convergence with speed-up factors, cost per
  outer iteration, iterations vs time trade-off (`speed_*.pdf`), from
  `results/tests/*.json` and the run logs.

Re-render and rebuild:

    ~/OF/venv/bin/python bench/make_report.py     # all figures, calls pvbatch
    CF_NO_RENDER=1 ~/OF/venv/bin/python bench/make_report.py   # without 3D
    nice -n 19 pvbatch --force-offscreen-rendering bench/render_fields.py --cases T4b
    cd report/paper && nice -n 19 make -B

## License

GPL-3.0-or-later (the code links against OpenFOAM, which is GPL-3.0).

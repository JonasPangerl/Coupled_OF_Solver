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
`coupledFoam: effective settings`.

## Tests and benchmarks

    python3 -m venv ~/OF/venv && ~/OF/venv/bin/pip install -r requirements.txt
    ~/OF/venv/bin/pytest tests/                 # writes results/tests/*.json
    ~/OF/venv/bin/python bench/run_bench.py     # writes results/bench/*.json
    ~/OF/venv/bin/python bench/make_report.py   # report/REPORT.md + paper figures

## Repository layout

    src/            libcoupledFoam (blockMatrix, blockSolvers, assembly, control, io)
    applications/   coupledFoam solver, unit-test applications
    cases/          test cases T0..T5 (Allrun/Allclean), shared case files
    tests/          pytest suite
    bench/          benchmark harness and report generator
    results/        generated JSON (only results/gates/ is versioned)
    report/         generated REPORT.md, LaTeX paper and figures
    scripts/        developer tools (syntax check, environment helpers)

## License

GPL-3.0-or-later (the code links against OpenFOAM, which is GPL-3.0).

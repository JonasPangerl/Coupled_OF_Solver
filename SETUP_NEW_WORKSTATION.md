# Setting up coupledFoam on another workstation

Written 2026-09-29 for an agent that continues this project elsewhere.
Everything below is verified on the machine the project was developed on.
Read `HANDOFF.md`, `DECISIONS.md` and `docs/KEYWORDS.md` after the setup.

## 1. What this project is

`coupledFoam` is an **add-on solver** for OpenFOAM ESI v2606: momentum and
continuity are solved as one 4x4 block system per cell (u, v, w, p), with
single-precision block linear algebra, double-precision reductions, a
block-GAMG-preconditioned FGMRES and pseudo-transient continuation. It never
modifies the OpenFOAM installation it is built against: the library goes to
`$FOAM_USER_LIBBIN`, the solver and test apps to `$FOAM_USER_APPBIN`.
Native solvers (`simpleFoam`) stay untouched and serve as the reference.

## 2. State to continue from

- Repository: <https://github.com/JonasPangerl/Coupled_OF_Solver>, branch `main`.
- `main` is the only branch that matters; every development branch is merged
  into it. Two superseded WIP branches (`t4b-debug`, `t4b-fix`) are pushed
  for completeness only - ignore them.
- Tag `campaign-20260923` (= commit cd7ebec) is the frozen commit of the last
  full test campaign; every result file in `results/` carries that commit.
- The last campaign ran: unit battery and gates, T0-T3 (1 and 4 ranks),
  restart, FPE, diagnostics, T4a (motorBike 354 k cells, 10 ranks),
  benchmarks on T1/T2/T3 and T4a, strong scaling on T4a with 1 and 16 ranks.
  T4b (1.7 M cells) and T5 (Ahmed body) were left out by user decision
  (D-063, D-072).
- Both report PDFs build from `results/` (`bench/make_report.py`, then
  `make` in `report/paper`).

## 3. Software the new workstation needs

| Item | Version used here | Note |
|---|---|---|
| OS | Ubuntu 24.04 LTS (under WSL2 on Windows) | any Linux with the packages below |
| OpenFOAM | **ESI OpenFOAM v2606**, package build `2606.0~rc2-1` | from the openfoam.com Debian/Ubuntu repository |
| OpenFOAM build | `linux64GccDPInt32Opt` (DP, `WM_LABEL_SIZE=32`, Opt) | DP is the reference precision (D-001) |
| Compiler | GCC >= 13 | |
| MPI | OpenMPI >= 4.1 (`FOAM_MPI=sys-openmpi`) | |
| Decomposition | `scotch-system` | see the Scotch pitfall in section 8 |
| Python | 3.12 in a venv | numpy 1.26, matplotlib 3.6, pytest 9.1, vtk 9.2, pandas, pyyaml, numpy-stl |
| ParaView | `pvbatch` (system package) | only for the 3D report figures |
| LaTeX | `pdflatex` + `latexmk` | only for the report PDFs |
| Hardware here | 16 physical / 32 logical cores, 94 GB RAM | T4a needs ~10 ranks, T4b ~1.7 M cells |

Install OpenFOAM (Ubuntu):

    curl -s https://dl.openfoam.com/add-debian-repo.sh | sudo bash
    sudo apt-get install openfoam2606-default openfoam2606-dev \
         openfoam2606-source openfoam2606-tools openfoam2606-tutorials

`openfoam2606-source` is only needed for the optional single-precision
build (section 7). The system install lands in
`/usr/lib/openfoam/openfoam2606` and is read-only - never write into it,
never run `apt` on `openfoam*` while work is in progress.

## 4. Setup, step by step

    # 1. repository
    git clone https://github.com/JonasPangerl/Coupled_OF_Solver.git ~/coupledFoam
    cd ~/coupledFoam && git checkout main

    # 2. python venv (the harness, benchmark and report run from it)
    python3 -m venv ~/OF/venv
    ~/OF/venv/bin/pip install numpy pandas matplotlib pytest pyyaml numpy-stl vtk

    # 3. build (exactly one OpenFOAM environment sourced)
    source /usr/lib/openfoam/openfoam2606/etc/bashrc
    ./Allwmake -j 8          # library, solver, test apps, utilities

    # 4. check the binary really is the private install
    ldd $(which coupledFoam) | grep libcoupledFoam
    # -> ~/OpenFOAM/<user>-v2606/platforms/linux64GccDPInt32Opt/lib/libcoupledFoam.so

`./Allwclean` removes all build products. After changing anything in
`src/`, rebuild the library **and** the applications (`Allwmake` does
both); a stale test binary has produced false results before.

## 5. Verify the setup

    cd ~/coupledFoam
    CF_MPI_BIND=none ~/OF/venv/bin/python -m pytest \
        tests/test_env.py tests/test_unit.py tests/test_gates.py \
        tests/test_keywords.py tests/test_harness.py --ranks 1 -q

All of these must pass; they need no long solver runs. Then one real case:

    CF_MPI_BIND=none ~/OF/venv/bin/python -m pytest tests/test_T1_pitzDaily.py --ranks 1 -q

Known failures that are **not** setup problems (see section 9):
`tests/test_unit.py::test_blockGAMG_cycles` (cycle ordering K <= W <= V)
and `tests/test_T3_airFoil.py::test_T3[np1-GEKO]`.

The heavy cases need the `--heavy` flag and build their mesh once:

    CF_FORCE_HEAVY=1 CF_HEAVY_NP=10 ~/OF/venv/bin/python -m pytest \
        "tests/test_T4_motorBike.py::test_T4[a]" --heavy -q

## 6. Two helper scripts (recommended)

`~/bin/cfenv` runs a command in exactly one OpenFOAM environment and
refuses to stack environments - the single most useful guard rail:

    #!/bin/bash
    # cfenv <env> <command...>   env: sys | sp
    set -e
    envName="$1"; shift
    case "$envName" in
        sys) rc=/usr/lib/openfoam/openfoam2606/etc/bashrc ;;
        sp)  rc="$HOME/OpenFOAM-v2606-SP/etc/bashrc" ;;
        *) echo "unknown env '$envName' (sys|sp)" >&2; exit 2 ;;
    esac
    [ -f "$rc" ] || { echo "missing $rc" >&2; exit 2; }
    if [ -n "$WM_PROJECT_DIR" ]; then
        echo "refusing: an OpenFOAM environment is already sourced" >&2; exit 2
    fi
    args=("$@"); set --
    set +eu; source "$rc"; set -e
    export WM_NCOMPPROCS="${WM_NCOMPPROCS:-8}"
    echo "[cfenv] $rc -> $WM_OPTIONS"
    cd "${CFENV_CWD:-$HOME/coupledFoam}"
    exec nice -n 19 "${args[@]}"

`~/bin/cf-fixown` is only needed when files are written into WSL through the
`\\wsl.localhost` share from Windows (they arrive owned by root). It is in
this repository's history; on a native Linux workstation it is unnecessary.

## 7. Optional: the single-precision build (amendment D, D-062/D-064)

Only needed for SP-vs-DP studies. SP is a **solve-only** mode; meshes are
always generated in DP.

    rsync -a --exclude '/platforms' /usr/lib/openfoam/openfoam2606/ ~/OpenFOAM-v2606-SP/
    chmod -R u+w ~/OpenFOAM-v2606-SP
    sed -i 's/^export WM_PRECISION_OPTION=.*/export WM_PRECISION_OPTION=SP/' ~/OpenFOAM-v2606-SP/etc/bashrc
    # GOTCHA: the packaged bashrc hard-codes the system project dir
    sed -i 's|^export WM_PROJECT_DIR=.*|export WM_PROJECT_DIR="'$HOME'/OpenFOAM-v2606-SP"|' ~/OpenFOAM-v2606-SP/etc/bashrc
    source ~/OpenFOAM-v2606-SP/etc/bashrc      # must give linux64GccSPInt32Opt
    (cd $WM_PROJECT_DIR/wmake/src && ./Allmake)
    (cd $WM_PROJECT_DIR/src && ./Allwmake -j 12 -s -q)
    # then the applications you need: simpleFoam potentialFoam decomposePar
    # reconstructPar checkMesh transformPoints foamDictionary blockMesh ...

Full library build takes about 2 h on 12 threads. Afterwards build the
project into the SP install: `CF_REPO=~/coupledFoam cfenv sp bash
bench/build_sp_install.sh`. The DP and SP installs live in different
`$WM_OPTIONS` platform directories and do not collide.

SP status: the solver compiles warning-free, the unit battery passes, and
T0/T1 agree with DP far inside the limits. **Not usable yet:** the
simpleFoam SP references hit the float floor and never satisfy their
residual control, so a fair SP-vs-DP benchmark is still open (D-070).

## 8. Pitfalls that cost time here

- **`mpirun --cpu-set` takes CORE ids (0-15 on a 16-core machine), not
  logical CPU ids.** An out-of-range set (e.g. `0-19`) fails **silently**
  with rc 1 and an empty log. The harness reads `CF_MPI_CPUSET`.
- **Scotch decomposition is not reproducible with threads.** Every case
  `Allrun` sets `SCOTCH_PTHREAD_NUMBER=1`; without it two runs on the same
  mesh differ (this was the cause of "np4 is not bit-reproducible").
- **`foamDictionary` on `fvSolution` must use `-disableFunctionEntries`**,
  otherwise `#sinclude "relaxation"` is silently expanded away. On
  `controlDict`/`fvSchemes` it must **not** be used (they contain `$`
  macros).
- One OpenFOAM environment per shell. Sourcing two (DP then SP) breaks the
  build in confusing ways - that is what `cfenv` guards against.
- Heavy runs: bind ranks (`CF_MPI_CPUSET=0-9`) and let nothing else run
  when a timing matters.

## 9. Known open points (state 2026-09-29)

1. **D-068, the common averaging window** (lead decision, user
   confirmation pending): the stationary-mean criterion of the wake cases
   now uses ONE window for both solvers of a case. This is why T4a shows a
   speed-up of 0.58-0.71 instead of the earlier ~2: the old rule derived
   the window from each run's own iteration budget and let simpleFoam
   count as converged only after 1500 iterations.
2. **T3-GEKO fails its 5 % tolerance.** The case has several steady
   branches; after the pressure-boundary fix the run converges to the
   stalled branch, with a later SFD start to the attached one. Neither is
   within 5 % of simpleFoam. Needs a user decision (test and report
   treatment).
3. **`coupled.remediation.wall.wallStarved` is off by default** (D-061).
   With the mild treatment of the pre-selected cells it is harmless and
   actually matched the reference best on T4a - turning it back on is a
   pending user decision.
4. **SP benchmark** blocked as described in section 7 (D-070).
5. `tests/test_unit.py::test_blockGAMG_cycles` fails (cycle iteration
   ordering) - pre-existing, recorded, not a regression.
6. T4b and T5 have no results in the current campaign by user decision.

## 10. Running a different case (e.g. an external-aerodynamics half-car)

The solver is case-agnostic; nothing is hard-coded for the test cases.

1. **Mesh**: generate as usual (`blockMesh` + `snappyHexMesh`, or cfMesh)
   in the **DP** build, then `decomposePar`. coupledFoam reads an ordinary
   OpenFOAM case.
2. **Copy the solver settings** from the closest template - for an
   external-aerodynamics case that is `cases/T4_motorBike`:
   - `system/fvSolution`: the `solvers.coupled { }` block (FGMRES +
     block-GAMG + blockILU0, V cycle, `agglomerationWeights combined`) and
     the top-level `coupled { }` block (outer iteration, PTC, line search,
     remediation, bounds, convergence).
   - `system/fvSchemes`: the divergence scheme for `div(phi,U)` is the
     high-order target; the solver runs deferred correction against upwind
     internally.
   - `system/coupledForcesDict` (optional but recommended): the
     double-accumulated force coefficients; the same numbers are then used
     by the convergence criterion and by `simpleFoam` through the
     `coupledForcesFO` function object, so both solvers are measured
     identically.
   - `system/relaxation.coupledFoam` if the template's `#sinclude
     "relaxation"` is kept.
3. **Check these keywords for a new geometry** (all documented with
   defaults and ranges in `docs/KEYWORDS.md`):
   - `coupled.Uref` (default `boundary` = the maximum of the non-coupled
     boundary values). For a case driven by an inlet this is right; verify
     the value printed in the effective settings.
   - `coupled.startupMode` (default `hybrid`) and `coupled.potentialInit`
     (default on): with a potential-flow start the peaks can be far above
     the freestream; that is handled (D-057), but watch the first
     iterations.
   - `coupled.remediation.*`: four categories (`meshQuality`, `badMesh`,
     `processor`, `wall`) with individual switches and treatment. On a
     snappyHexMesh mesh the static set should stay around or below 1 % of
     the cells - the start-up log prints the count per category.
   - `coupled.convergence`: for a wake case the force window decides;
     budget 400-800 outer iterations, not thousands.
   - MRF is supported (rotating zones, implicit coupling); the test case is
     `cases/T_mrf_mixerVessel2D`.
4. **Run**: `mpirun -np <N> coupledFoam -parallel` (or the case `Allrun`).
   Every outer iteration prints one `CF|` line: residuals, CFL, line-search
   factor omega, linear iterations, remediation set sizes, timings. An
   abort writes the last accepted state to `<iter>_lastValid` together with
   `sentinelFlag` and `remediationFlag` for ParaView.
5. **Compare against `simpleFoam`** on the same mesh and the same number of
   ranks, and judge both by the same criterion (D-068). Report wall-clock
   **and** CPU-hours.

## 11. Working rules of this project

- Everything written to git is in English; commit per subsection.
- Never loosen a test threshold without an explicit user decision, and
  record every such decision in `DECISIONS.md`.
- Report every performance number as wall-clock **and** CPU-hours.
- Test a change on T4a (354 k cells, ~12 min on 10 ranks) before T4b.
- The system OpenFOAM installation stays untouched.

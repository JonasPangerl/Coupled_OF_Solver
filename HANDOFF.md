# HANDOFF - state of coupledFoam and review list

Written 2026-09-21 (evening) before a session compaction. Read this first,
then `PLAN.md`, `DECISIONS.md` (D-001..D-029), `SPEC_amendment_B.md`,
`AMENDMENT_B_STATUS.md`.

## 1. Where things are

- Repository: WSL `/home/jonas/coupledFoam`, remote
  `https://github.com/JonasPangerl/Coupled_OF_Solver` (branch `main`,
  last pushed commit `d0ace1f`). Commits follow `phase-<N>: <what>`.
- Build: against the read-only system OpenFOAM v2606 DP
  (`/usr/lib/openfoam/openfoam2606`, NEVER modify, NEVER apt on openfoam*).
  Products go to `~/OpenFOAM/jonas-v2606/platforms/linux64GccDPInt32Opt/`.
- Python venv for tests/report: `~/OF/venv` (system-site-packages + pandas,
  pytest, numpy-stl).
- Helper scripts (outside the repo): `~/bin/cfenv <sys|dp-debug|dp-opt> cmd`
  (sources exactly one OpenFOAM bashrc, niced), `~/bin/cf-fixown` (re-owns
  root-owned files written through `\\wsl.localhost`, restores exec bits).
  Scratch scripts: Windows scratchpad `.../scratchpad/build_all.sh`
  (build lib + apps, `build_all.sh 2` = 2 jobs), `pytest_run.sh`,
  `exp_re1000.sh` (variant runs of T0 Re1000).
- Tool gotchas (also in Claude memory `wsl-workflow-gotchas`): use
  `MSYS_NO_PATHCONV=1 wsl.exe -e <abs path> args`; never
  `wsl.exe -- bash -c '...$var...'`; UNC writes are root-owned -> run
  cf-fixown; OpenFOAM `etc/bashrc` sources files passed as positional args
  (cfenv clears them).

## 2. Uncommitted work at handoff time

Nothing any more: the processor agglomeration (6.3.4, D-030) is committed
and verified (serial bit-identical, 4-rank on/off <= 2e-6 for V/F/W and
dense LU), as is the CF_MPI_CPUSET / CF_HEAVY_NP plumbing. New since the
handoff: amendment B11 (hot-loop performance rules, SPEC_amendment_B11.md,
still to implement) and the user directives in D-031.

## 3. Results so far (all 1 rank unless stated)

| Test | Result |
|---|---|
| test_env, Test-precision, block4Ops, doubleReduce, blockMatrix (1+4 ranks), blockGAMG (1+4, tol 1e-9, D-023) | pass |
| T0 Re100 | pass: 58 outer its, 9.6 s; profiles 2.9e-6 / 5.5e-6 vs simpleFoam (1932 its, 173 s) |
| T0 Re1000 | pass: 58 its, 17.9 s; profiles 1.4e-5 |
| T0 np4 (both Re) | **FAIL: MPI_ERR_TRUNCATE** in coupledFoam (see 4.1) |
| T1, T2, T3, T-restart, T-fpe | written, not yet run |
| T4, T5, T-scaling (heavy) | written, not run (user: no heavy tests until the F1 job is stopped and the user says go) |
| Benchmark (A-G), report, paper PDF | harness/generator written; paper needs LaTeX (user is installing via sudo) |

Key finding (D-028): the block Gauss-Seidel smoother is inadequate for the
coupled saddle-point system; ILU0 (block DILU) or FGMRES+K-cycle fixes it.
T0 Re1000 before amendment B: not converged in 300 its (978 s); with
FGMRES+K+Eisenstat-Walker: 58 its / 18 s (GS) or 69 its / 6 s (ILU0).

## 4. Open problems (do these first)

1. **Parallel MPI_ERR_TRUNCATE** (T0 np4, coupledFoam with FGMRES + K-cycle).
   The 4-rank Test-blockGAMG with BiCGStab + V passed earlier. Suspects:
   the build may have contained a half-edited blockGAMG from the
   proc-agglomeration agent; otherwise a collective mismatch in the K-cycle
   or coarsest solve (e.g. useDenseLU_ decided differently on different
   ranks, a rank with 0 coarse cells, K-step dot products on a level
   communicator, message tags of blockLduInterface reused by two
   simultaneous exchanges on different levels). Reproduce with
   `mpirun --bind-to none -np 4 Test-blockGAMG -parallel -cycle K` on
   `run/unit_cavity`, then coupledFoam on `run/T0_Re100_np4`.
   New evidence (D-030 verification): the K cycle stalls just above the
   1e-9 tolerance (final 1.8e-8) on unit_cavity in the *pre-change
   baseline* as well, at 4 ranks with agglomeration off - so there is a
   pre-existing K-cycle/FGMRES defect (or a float-precision stall of the
   K-step GCR) independent of 6.3.4. V/F/W are fine on the same system.
2. **Formal ILU0 decision** (spec 6.3: Gate B2, >= 15 % wall-clock gain on
   T2): run T2 with blockGaussSeidel vs blockILU0 (FGMRES+K), record in
   DECISIONS and set the case default. T0 data already strongly favour ILU0.
3. Run T1, T2, T3 (SST, GEKO), T-restart, T-fpe (pytest), fix failures.
4. Heavy (only after the user's go): mesh T4a, T4b, T5 (`run/*_mesh`
   caches), T4/T5 tests, T-scaling, benchmark A-G (3 repeats; with less time
   at least 1 repeat on the big cases), then `bench/make_report.py` and
   `make -C report/paper`.
5. OpenFOAM DP Debug build in `~/OF/OpenFOAM-v2606-DP-Debug` (needs bison,
   being installed) for the Debug-build FPE tests (spec 9.1/13); several
   hours on 8 cores, only when the machine is free.
6. **F1 readiness** (the user wants to use coupledFoam on the F1 half-car
   cases next): read the F1 case in the case repo
   `C:\Users\Dell T5600\OneDrive\Dokumente\SIMS\OpenFoam_cases`
   (read-only; never touch `f1_halfcar/runs/*`), list its BCs / function
   objects / schemes / MRF or rotating walls / cyclicAMI and check each is
   supported by coupledFoam (D-002 generic BC linearisation; cyclic/AMI are
   explicit per spec 5.5). Prepare a ready-to-use `fvSolution` block for
   coupledFoam on that case and a short how-to.

## 5. Review list for a stronger agent (before the big tests)

Correctness-critical, worth a careful independent read:

1. `src/assembly/coupledAssembler.C`: signs and completeness of the 4x4
   block entries (momentum, grad p, div U, Rhie-Chow p-p), processor
   interface coefficients (native sign convention result -= bouCoeffs),
   boundary linearisation via valueInternal/BoundaryCoeffs, the strong
   pressure reference (D-021), the double residual and normFactor.
   Evidence it is right in serial: T0 matches simpleFoam to 3e-6.
   Not yet verified: processor faces (parallel failed), cyclic/explicit
   patches, MRF, non-orthogonal meshes (T1-T5).
2. `src/assembly/rhieChow.C`: flux update must equal the continuity row
   exactly (same face formula); D_f boundary policy (D-013); cellLimited
   gradient in static cells (D-018).
3. `src/blockSolvers/blockGAMG.C`: K-cycle (Notay-Vassilevski GCR),
   W/F recursion work-vector aliasing, coarse Galerkin sums incl.
   faceFlipMap, dense LU, processor agglomeration (new), parallel
   collectives on level communicators.
4. `src/blockSolvers/blockFGMRES.C`, `blockGMRES.C`, `blockBiCGStab.C`:
   Givens/back-substitution, restart logic, best-iterate return, rho.
5. `applications/coupledFoam/coupledFoam.C`: outer-loop order (5.7),
   B4 failure/abort path, Eisenstat-Walker placement (eta per assembly,
   accept only on accepted steps), Anderson placement and flushes,
   sentinel rollback, restart state completeness (T-restart requires
   +-2 iterations and 1e-5 on dp/Cd).
6. `src/control/anderson.C` (D-026 formula), `gamgAutoTune.C`,
   `remediation.C` zonal part (agent-written, compiled, never run).
7. `bench/run_bench.py`: configuration definitions A-G (D-025), time to
   convergence scaling by tWall/ClockTime, CPU-hours from rank_wrapper.
8. `report/paper/paper.tex`: every claim against DECISIONS.md and the
   code; the paper agent flagged: Re1000 discussion numbers were
   hand-copied from exploratory runs (replace by generated ones), author
   line placeholder, Ceze & Fidkowski 2013 paper number to verify.
9. Spec compliance gaps to re-check: gradient limiter in static cells only
   for p (D-018); collated file handler atomic write (D-017); Anderson
   history not in restart state (D-026).

## 6. User preferences (standing)

- System OpenFOAM stays untouched; coupledFoam is an add-on.
- Machine shared with F1 runs: check load before any build/run; tiny
  tests only while F1 runs; heavy tests only after the user says so.
  Since 2026-09-21 evening (D-031): the F1 job is being stopped; heavy
  runs use at most 10 cores (`CF_HEAVY_NP=10`, now the default) with core
  binding; `CF_MPI_CPUSET` is no longer needed once F1 is gone.
- Every performance result with wall-clock AND CPU-hours; the goal is the
  fastest solver, not the fewest iterations.
- Everything in git in English; final English LaTeX paper with vector
  figures, regenerable from `results/`. The paper is delivered as a
  compiled PDF, with intermediate PDF builds so the user can inspect
  drafts (LaTeX is installed). TWO versions from the same generated
  numbers/figures (D-032): paper.pdf (professional) and
  paper_tutorial.pdf (plain-language, for a reader without CFD
  background - the user wants to learn from it).
- Parallel subagents are welcome to speed up work (give them exclusive
  files; the lead integrates, builds and commits).

## 7a. Work order

OPUS_TASKS.md (written by Fable, 2026-09-22 ~08:00) is the binding,
prescriptive task list for the day: robust preconditioner defaults
(grid + selection rule), the user-approved averaging criterion for
T4/T5 (D-042 text ready to paste), the T3 probes, and the remaining
pipeline. Follow it before this section.

## 7. State at model handover, 2026-09-22 ~07:40 (by Fable)

Read FABLE_REVIEW.md first - it carries the live issue list. Summary:

Done tonight (all committed, main, pushed):
- Bugs fixed: MPI_ERR_TRUNCATE (per-patch collectives), FGMRES/GMRES
  float stall (mixed-precision refinement), parallel restart (per-rank
  state; ASCII header), nut cap 1e5->1e8 (D-040, T3 root cause),
  T4/T5 #include -> #sinclude, B4 linear budgets 200->400 on T2-T5.
- D-039 preconditioner study merged: matrix-weighted agglomeration,
  1 finest ILU0 sweep -> T1 33 s vs simpleFoam 36 s (was ~5x slower).
- B11 kernels merged (D-041): fused SIMD kernels, ~37 % faster linear
  iterations; gate phase_A PASS (Amul 85.6 % of STREAM).
- T0 Re100/1000 x np1/np4 all PASS (Re100 np1: 57 its, 6.2 s).
  T0-Re1000 iteration counts are chaotic (~2x scatter) - do not quote
  single runs (FABLE item 1).
- T4a mesh (354k) + reference (3000 its, 1516 s, 4.2 CPU-h) and
  T4b mesh (1.70 M) + reference (4000 its, 11876 s, 33.0 CPU-h) cached
  in run/. Both references fail criterion 12.3(ii) - physically
  oscillating wakes (FABLE item 2c, user decision needed).
- Both papers rebuilt from generated numbers (report/paper/*.pdf).

Open, in priority order (details in FABLE_REVIEW item 4):
1. D-039 settings are not robust off T1: B4 aborts on T2 (fixed by
   maxIter 400), T3 (both probes), T4a (0.64x reduction in 400 its at
   CFL 500). A conservative T4b attempt (V cycle, ILU0 x2, geometric,
   CFLmax 100) is RUNNING in run/T4b_conservative - check its log.
2. T3-SST limit cycle (nut cap fixed, cycle remains); T3-GEKO 3.9 % Cd.
3. T1 495 its vs 400-it budget; T2 xr 0.56 % vs 0.5 %, R stalls 1e-4.
4. Not run: T5 (mesh pending), T-scaling, benchmark A-H, T-restart/
   T-fpe pytest (slow-case variants), DP-Debug build.
5. Report items: exploratory eta/rho figures fill once E/F/G/H run.

Processes possibly still alive from tonight: run/T4b_conservative
(mpirun np10). Nothing else. All agent worktrees (/home/jonas/cf_precond,
/home/jonas/cf_b11, cf_b11_base*, cf_precond_platform...) are merged and
can be deleted with git worktree remove --force + rm -rf of the
_platform dirs.

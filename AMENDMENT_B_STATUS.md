# Amendment B - implementation status

Source: `SPEC_amendment_B.md` (verbatim). Every item is tracked here; nothing
is dropped silently (spec rule 0.3). Status: `todo`, `code` (written, builds),
`tested` (unit/case test passes), `n/a` (documentation only, done).

| Item | What | Status | Where |
|---|---|---|---|
| B1 | `blockFGMRES` (flexible, V and Z stored, MGS in double, happy-breakdown guard) | code | src/blockSolvers/blockFGMRES |
| B1 | `blockGMRES` standard (right-preconditioned, fixed preconditioner) | code | src/blockSolvers/blockGMRES |
| B1 | keyword `restart` (10; 6 above 40 M cells), `minIter 1` | code | blockSolver, coupledDefaults |
| B1 | FatalError: cycleType K requires blockFGMRES | code | blockGAMGPrecon |
| B2 | Eisenstat-Walker adaptive relTol with safeguards 1-4, `eta=` in CF| | code | control/adaptiveTolerance, coupledFoam |
| B3 6.3.1 | defaults mergeLevels 2, nCellsInCoarsestLevel 200; ratio rule r_l >= 3 for W/K (mergeLevels up to 4, FatalError) | code | blockGAMG |
| B3 6.3.2 | cycles V, F, W, K (Notay-Vassilevski GCR K-step, cached workspace) | code | blockGAMG |
| B3 6.3.3 | coarsest: BiCGStab 1e-3 / 50 its; dense LU (double) if <= 64 cells on one rank | code | blockGAMG, blockDenseLU |
| B3 6.3.4 | processor agglomeration (masterCoarsest, rule 5000 cells/rank), block coefficient gather | done (D-030) | blockGAMG, blockGAMGProcAgglomeration.H |
| B3 6.3.5 | autoTune controller (rho on first preconditioner application, windows, hysteresis, events) | code | blockGAMG, coupledFoam |
| B4 | linear-solve failure trigger, maxLinFails abort via 9.3 path | code (D-029) | coupledFoam |
| B5 | Anderson acceleration (Type II, safeguards, flush rules) | code (D-026) | control/anderson |
| B6 | zonal factors (cellZones, patch-distance layers by BFS) | code (D-027) | control/remediation |
| B7 | memory table | code: paper Table tab:memory (Anderson row corrected to 320 B/cell = 14.4 GB, D-026, value generated); bench_memory figure draws the B7 budget scaled per cell against measured peak RSS (sum over ranks). Measured comparison pending benchmark runs | report/paper/paper.tex, bench/make_report.py |
| B8 | Section 11 keywords and defaults | code (coupledDefaults.H, cases T0-T3) | coupledDefaults.H, case fvSolution |
| B9 | Test-blockGAMG cycle comparison (cavity + motorBike 350k) | code | applications/test |
| B9 | Test-blockFGMRES (variable preconditioner) | code | applications/test |
| B9 | Test-procAgglom (16 ranks, motorBike) | code (heavy test) | applications/test |
| B10 | benchmark configs E (fixed V: cycleType V + autoTune no), F (adaptiveRelTol no), G (anderson yes), plus H (fixed K, autoTune no, reference for E); default case scoping E/H: T2,T4b,T5, F/G: T1,T3-SST; acceptance evaluation (F: C <= 5 % slower in wall AND CPU-h, Cd/dp to 1e-4; E/G/H deltas) into summary.json, b10_acceptance.csv, table b10_acceptance + macros | code (synthetic check), runs pending | bench/run_bench.py, bench/make_report.py |
| B10 | report 15.7: eta, rho + tune events, per-level table (cells, ranks, ratios, cycle/nPostSweeps start/end, coarsest solver), cycle bars (unit test + benchmark C/H/E), Anderson (iterations, wall, CPU-h) | code; figures appear once test/benchmark JSONs exist (T0 rows of the level table already generated) | bench/make_report.py, tests/cflib/logs.py |

## Amendment B11 (6.5 hot-loop performance rules, SPEC_amendment_B11.md)

| item | what | status | where |
|---|---|---|---|
| B11 6.5.1 | no Field algebra / tmp<> in solver kernels (flat blockScalar arrays) | done (audit: no Field/tmp<> was present; the remaining List operations are flat fills/copies) | blockMatrix, blockSolvers |
| B11 6.5.2 | fused kernels: axpy_dot, update_residual_norm, fused MGS; double reductions | done | blockMatrix/blockKernels.H, blockBiCGStab, blockGMRES, blockFGMRES |
| B11 6.5.3 | omp simd (-fopenmp-simd), __restrict__, 64-byte alignedList, fixed-trip 4x4 matvec | done (Krylov work vectors alignedList; matrix arrays keep List, v2606 aligns them to 256 B, checked at runtime by the gate) | blockMatrix/alignedList.H, block4Ops.H, Make/options |
| B11 6.5.4 | block AoS layout stays (no SoA) | holds already (5.2) | - |
| B11 6.5.5 | assembly writes flat arrays in one owner/neighbour pass; fvm:: only to obtain coefficients | verified, see FABLE_REVIEW / final report of the B11 session | coupledAssembler |
| B11 6.5.6 | Test-kernelBandwidth: 5 M cells, Amul >= 60 % STREAM, axpy_dot >= 80 % -> results/gates/phase_A.json, report 15.5 | done, gate PASS (Amul 85 %, axpy_dot 135 % of triad) | applications/test/Test-kernelBandwidth |
| B11 6.5.7 | DECISIONS entry: expression templates out of scope | text proposed in the B11 session report (not yet in DECISIONS.md) | DECISIONS.md |

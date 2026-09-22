# Precision fields of the result records (amendment D11, D-065)

Producer: `tests/cflib/precision.py` (tests, `precision.annotate`) and
`bench/run_bench.py` (`precision_fields`). Consumer: the report (Section
15.10 "Single precision"). Reference precision is **DP** (D-062): where the
amendment says `Cd_rel_to_SPDP`, the key is `Cd_rel_to_DP`.

## Where the records are

| kind | DP record | SP record |
|---|---|---|
| test | `results/tests/<name>.json` | `results/tests/<name>_sp.json` |
| benchmark | `results/bench/<case>_<cfg>_<run>.json`, cfg A..H, E-* | `results/bench/<case>_F1_<run>.json` (simpleFoam SP, settings of B), `<case>_F2_<run>.json` (coupledFoam SP, settings of C) |
| run directory | `run/<name>` | `run/<name>_sp` (holds `spHarness.json`, `constant/meshShift`, `log.checkMesh.dp`, `log.checkMesh.sp`) |

F1/F2 scope (D-063): T1, T3-SST, T4a; one repeat by default. The DP
counterpart of F1 is B, of F2 is C (same case). A precision-harness
self-test (SP procedure with the DP build, `CF_SP_ALLOW_DP_BUILD=1`) writes
`precision = "dp-shifted"`; it is not an SP result.

## Fields (every record written after D-065)

| key | type | DP record | SP record |
|---|---|---|---|
| `precision` | string | `"dp"` | `"sp"` (`"dp-shifted"`: harness self-test) |
| `buildPrecision` | string | `"DP"` | `"SP"` (WM_PRECISION_OPTION of the build) |
| `meshShift` | [sx, sy, sz] or null | null | origin shift applied to the DP mesh: points_SP = points_DP + shift (= minus the bounding-box centre of the DP checkMesh) |
| `checkMeshDiff` | list or null | null | failed-check kinds that differ between DP and SP `checkMesh -allGeometry -allTopology`: `[{"check": "<message without numbers>", "only": "sp" or "dp"}]`; `[]` = identical |
| `nCheckMeshDiff` | int or null (bench) | null | `len(checkMeshDiff)` |
| `spGeometry` | string or null | null | `"ok"` or `"SP-geometry-fail"` (D5.3; the run is then skipped: bench record `skipped: true`, no timings) |
| `staticSetSize` | int or null | cells in the static remediation set (first `nStat` of the coupledFoam log) | the same, computed by the SP build (D5.5) |
| `staticSetSizeDP` | int or null | null | `staticSetSize` of the DP counterpart |
| `staticSetSizeDiff` | int or null | null | SP - DP |
| `<m>_rel_to_DP` | float or null | null (tests: only the metrics of the test) | (SP - DP)/abs(DP), signed; tests: `Cd`, `Cl` (T3), `dp` (T1), `xr` (T2); bench: `Cd`, `Cl`, `dp` (at convergence) and `Cd_final`, `Cl_final`, `dp_final` (end of the fixed budget) against the median of the current DP counterpart records |
| `<m>_DP` | float or null | absent | the DP value used |
| `Cd_rel_to_DP` | float or null | null | always present (null for the dp-monitored T1/T2) |
| `spVsDpTol` | float | absent/null | 0.003 (D10: dp/Cl/Cd within 0.3 % of DP) |
| `spVsDpPass` | bool or null | null | all monitored `<m>_rel_to_DP` within `spVsDpTol` (T0: centreline L2 < 1e-3, keys `l2rel_u_vs_DP`, `l2rel_v_vs_DP`); null if the DP counterpart is missing |
| `speedupWall_DP_over_SP`, `speedupCpu_DP_over_SP` | float or null (bench) | absent | DP counterpart median wall / CPU-h to convergence over this run's |
| `dpCounterpart` | object | absent | tests: `{"record", "case"}`; bench: `{"config", "runs"}` |
| `spHarness` | object | absent | `bboxCentreDP`, `shiftedEntries` (file, key, old, new), `unclassifiedVectors` (three-component entries the harness did not classify - should be empty), `reasons` (gate), `wallSecondsDPMesh`, `wallSecondsPrepare`; tests also `checkMeshFailedDP`, `checkMeshFailedSP` (raw *** lines) |
| `forceSource` | string | `"coupledForces"` or `"forceCoeffs"` (which file the forces came from, D6) | same |

Summary (`results/bench/summary.csv`/`.json`, per case and configuration):
`precision`, `Cd_rel_to_DP`, `Cl_rel_to_DP`, `dp_rel_to_DP`,
`staticSetSizeDiff`, `spGeometry`, `nCheckMeshDiff` (medians over the
repeats), `speedup_wall_DP_over_SP`, `speedup_cpu_DP_over_SP` (DP
counterpart median over SP median, F1 and F2 rows).

## Force coefficients (D6)

`postProcessing/coupledForces/<startTime>/coeffs.dat`: header lines start
with `#`, the column line is `# iter Cd Cl Cm Cd_p Cd_v Cl_p Cl_v`, values
`%.10e` (iter `%.10g`). Written by coupledFoam (every outer iteration) and
by the `coupledForcesFO` function object (simpleFoam, every time step).
`cflib.post.force_coeffs` prefers it over the native `forceCoeffs` file and
returns the same keys (`Time` = iter, `CmPitch` = Cm).

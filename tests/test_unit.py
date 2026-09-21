"""Phase A / B2 unit tests (spec 6.4): the Test-* applications on 1 and 4
ranks, with cross-rank comparison.

  Test-block4Ops     ||A A^-1 - I||_inf < 1e-5, singular block guarded
  Test-doubleReduce  double sum exact to 1e-9 (D-004), float sum wrong
  Test-blockMatrix   block Amul vs native Amul < 1e-6 ||x||; 1 vs 4 ranks 1e-6
  Test-blockGAMG     converged to 1e-8 in <= 20 iterations; 1 vs 4 ranks 1e-5

Amendment B9:
  test_blockGAMG_cycles            cycleType V, F, W, K on the cavity mesh,
                                   mergeLevels 2: all converge, iterations
                                   K <= W <= V, ratios r_l >= 3 for W and K,
                                   K on 1 vs 4 ranks identical to 1e-5
  test_blockFGMRES                 variable preconditioner: FGMRES converges,
                                   GMRES outcome recorded
  test_blockGAMG_cycles_motorBike  (heavy) the cycle study on the motorBike
                                   mesh (run/T4a_mesh)
  test_procAgglom                  (heavy) CF_HEAVY_NP ranks, motorBike mesh:
                                   processorAgglomerator on vs off to 1e-5,
                                   ranks per level follow rule 6.3.4
"""

from __future__ import annotations

import json
import os
import re
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pytest

from cflib import case as cfcase
from cflib import env as cfenv
from cflib import results

TOL_CROSS_MATRIX = 1e-6
TOL_CROSS_GAMG = 1e-5
# Solve tolerance of Test-blockGAMG: spec 1e-8, tightened to 1e-9 (D-023)
GAMG_TOLERANCE = 1e-9
RANKS = (1, 4)


def _app(cmd: list[str], cwd: Path, nprocs: int, json_out: Path,
         extra: list[str] | None = None) -> dict:
    full = list(cmd) + ["-json", str(json_out)] + list(extra or [])
    if nprocs > 1:
        full = cfenv.mpirun_prefix(nprocs) + full + ["-parallel"]
    log = json_out.with_suffix(".log")
    rc = cfenv.run(full, cwd=cwd, log=log)
    data = json.loads(json_out.read_text()) if json_out.exists() else {}
    data["rc"] = rc
    return data


@pytest.fixture(scope="module")
def cavity_mesh(foam):
    """T0 cavity mesh (128x128), decomposed for 4 ranks."""
    case = cfcase.prepare("T0_cavity", "unit_cavity")
    assert cfcase.allrun(case, ["-mesh-only"]) == 0
    cfcase.set_entry(case, "system/decomposeParDict", "numberOfSubdomains", 4)
    assert cfenv.run(["decomposePar", "-force"], cwd=case,
                     log=case / "log.decomposePar") == 0
    return case


@pytest.mark.unit
def test_block4Ops(foam, tmp_path):
    d = _app(["Test-block4Ops"], tmp_path, 1, tmp_path / "b4.json")
    results.write("tests", "Test-block4Ops", d)
    assert d["rc"] == 0 and d.get("pass"), d


@pytest.mark.unit
def test_doubleReduce(foam, tmp_path):
    out = {}
    for n in RANKS:
        out[f"np{n}"] = _app(["Test-doubleReduce"], tmp_path, n,
                             tmp_path / f"dr{n}.json")
    passed = all(v["rc"] == 0 and v.get("pass") for v in out.values())
    results.write("tests", "Test-doubleReduce", {"runs": out, "pass": passed})
    assert passed, out


@pytest.mark.unit
def test_blockMatrix(foam, cavity_mesh):
    out = {}
    for n in RANKS:
        out[f"np{n}"] = _app(["Test-blockMatrix"], cavity_mesh, n,
                             cavity_mesh / f"bm{n}.json")
    a = np.array(out["np1"]["sumMagAx"])
    b = np.array(out["np4"]["sumMagAx"])
    cross = float(np.max(np.abs(a - b) / np.abs(a)))
    passed = (all(v["rc"] == 0 and v.get("pass") for v in out.values())
              and cross < TOL_CROSS_MATRIX)
    results.write("tests", "Test-blockMatrix",
                  {"runs": out, "crossRankRelDiff": cross,
                   "crossRankTol": TOL_CROSS_MATRIX, "pass": passed})
    assert passed, (cross, out)


def _sorted_solution(path: Path) -> np.ndarray:
    """Solution rows sorted by cell centre (decomposition independent)."""
    d = np.loadtxt(path)
    order = np.lexsort((d[:, 2], d[:, 1], d[:, 0]))
    return d[order]


@pytest.mark.unit
def test_blockGAMG(foam, cavity_mesh):
    out = {}
    for n in RANKS:
        out[f"np{n}"] = _app(["Test-blockGAMG"], cavity_mesh, n,
                             cavity_mesh / f"bg{n}.json",
                             ["-dumpSolution", str(cavity_mesh / f"bg{n}.xyz"),
                              "-tolerance", str(GAMG_TOLERANCE)])
    # "Identical solution 1 vs 4 ranks to 1e-5" (6.4): relative L2 difference
    # of the whole solution vector, cells matched by their centres
    s1 = _sorted_solution(cavity_mesh / "bg1.xyz")
    s4 = _sorted_solution(cavity_mesh / "bg4.xyz")
    assert np.allclose(s1[:, :3], s4[:, :3]), "cell centres do not match"
    x1, x4 = s1[:, 3:], s4[:, 3:]
    cross = float(np.linalg.norm(x1 - x4) / np.linalg.norm(x1))
    max_abs = float(np.max(np.abs(x1 - x4)))
    passed = (all(v["rc"] == 0 and v.get("pass") for v in out.values())
              and cross < TOL_CROSS_GAMG)
    results.write("tests", "Test-blockGAMG",
                  {"runs": out, "crossRankRelDiff": cross,
                   "crossRankMaxAbsDiff": max_abs,
                   "crossRankMetric": "||x1-x4||_2/||x1||_2, cells matched by centre",
                   "crossRankTol": TOL_CROSS_GAMG, "pass": passed})
    assert passed, (cross, out)


# ---------------------------------------------------------------------------
# Amendment B9: cycle types, flexible GMRES, processor agglomeration
# ---------------------------------------------------------------------------

CYCLES = ("V", "F", "W", "K")
MERGE_LEVELS = 2
MIN_RATIO = 3.0                 # coarsening-ratio rule 6.3.1 (W and K)
PROC_AGGLOM_NPROCS = int(os.environ.get("CF_HEAVY_NP", "10"))
PROC_AGGLOM_CELLS_PER_RANK = 5000   # rule 6.3.4
PROC_AGGLOM_SINGLE_RANK = 5000      # rule 6.3.4
RANKS_LINE = "blockGAMG: ranks per level"
MOTORBIKE_MESH_CASE = "T4a_mesh"    # cases/T4_motorBike: Allrun -mesh a -mesh-only


def _compare_solutions(f1: Path, f2: Path) -> dict:
    """Relative L2 and max abs difference of two -dumpSolution files, cells
    matched by their centres (decomposition independent)."""
    s1 = _sorted_solution(f1)
    s2 = _sorted_solution(f2)
    if s1.shape != s2.shape or not np.allclose(s1[:, :3], s2[:, :3]):
        return {"centresMatch": False, "relDiff": None, "maxAbsDiff": None}
    x1, x2 = s1[:, 3:], s2[:, 3:]
    return {"centresMatch": True,
            "relDiff": float(np.linalg.norm(x1 - x2) / np.linalg.norm(x1)),
            "maxAbsDiff": float(np.max(np.abs(x1 - x2)))}


def _cycle_study(case: Path, tag: str) -> dict:
    """Test-blockGAMG with every cycle type (serial) and cycle K on 1 vs 4
    ranks. Returns the result record including "pass" and "failures"."""
    runs = {}
    failures = []
    common = ["-tolerance", str(GAMG_TOLERANCE),
              "-mergeLevels", str(MERGE_LEVELS)]

    for c in CYCLES:
        extra = common + ["-cycle", c, "-skipDiagonal"]
        if c == "K":
            extra += ["-dumpSolution", str(case / f"{tag}_K_np1.xyz")]
        d = _app(["Test-blockGAMG"], case, 1, case / f"{tag}_{c}_np1.json",
                 extra)
        d["processWallSeconds"] = cfenv.last_timing.wall
        d["processCpuSeconds"] = cfenv.last_timing.cpu
        runs[c] = d
        if not d.get("converged"):
            failures.append(f"cycle {c} did not converge (rc {d['rc']})")
        if d.get("cycleType") not in (None, c):
            failures.append(f"cycle {c}: hierarchy reports cycleType "
                            f"{d.get('cycleType')}")
        if c in ("W", "K"):
            ratios = d.get("ratios") or []
            if not ratios:
                failures.append(f"cycle {c}: no coarsening ratios recorded")
            low = [r for r in ratios if r < MIN_RATIO]
            if low:
                failures.append(f"cycle {c}: ratios {ratios} below {MIN_RATIO}")

    iters = {c: runs[c].get("nIterations") for c in CYCLES}
    if all(isinstance(iters[c], int) for c in ("V", "W", "K")):
        if not (iters["K"] <= iters["W"] <= iters["V"]):
            failures.append(f"iterations not K <= W <= V: {iters}")
    else:
        failures.append(f"missing iteration counts: {iters}")

    # Cycle K: 1 vs 4 ranks (6.4: identical solution to 1e-5)
    d4 = _app(["Test-blockGAMG"], case, 4, case / f"{tag}_K_np4.json",
              common + ["-cycle", "K", "-skipDiagonal",
                        "-dumpSolution", str(case / f"{tag}_K_np4.xyz")])
    d4["processWallSeconds"] = cfenv.last_timing.wall
    if not d4.get("converged"):
        failures.append(f"cycle K on 4 ranks did not converge (rc {d4['rc']})")
    f1, f4 = case / f"{tag}_K_np1.xyz", case / f"{tag}_K_np4.xyz"
    if f1.exists() and f4.exists():
        cross = _compare_solutions(f1, f4)
    else:
        cross = {"centresMatch": False, "relDiff": None, "maxAbsDiff": None}
    if not cross["centresMatch"]:
        failures.append("K 1 vs 4 ranks: solution dumps missing or cell "
                        "centres do not match")
    elif cross["relDiff"] >= TOL_CROSS_GAMG:
        failures.append(f"K 1 vs 4 ranks: rel. difference {cross['relDiff']:.3e}"
                        f" >= {TOL_CROSS_GAMG}")

    return {
        "mesh": str(case),
        "tolerance": GAMG_TOLERANCE,
        "mergeLevels": MERGE_LEVELS,
        "iterations": iters,
        "wallSeconds": {c: runs[c].get("wallSeconds") for c in CYCLES},
        "processWallSeconds": {c: runs[c].get("processWallSeconds")
                               for c in CYCLES},
        "processCpuSeconds": {c: runs[c].get("processCpuSeconds")
                              for c in CYCLES},
        "ratios": {c: runs[c].get("ratios") for c in CYCLES},
        "mergeLevelsUsed": {c: runs[c].get("mergeLevels") for c in CYCLES},
        "converged": {c: runs[c].get("converged") for c in CYCLES},
        "runs": runs,
        "K_np4": d4,
        "crossRankK": cross,
        "crossRankMetric": "||x1-x4||_2/||x1||_2, cells matched by centre",
        "crossRankTol": TOL_CROSS_GAMG,
        "failures": failures,
        "pass": not failures,
    }


@pytest.mark.unit
def test_blockGAMG_cycles(foam, cavity_mesh):
    """B9: cycleType V, F, W, K on the 128x128 cavity mesh (D-023 tolerance)."""
    rec = _cycle_study(cavity_mesh, "cyc")
    results.write("tests", "Test-blockGAMG_cycles", rec)
    assert rec["pass"], rec["failures"]


@pytest.mark.unit
def test_blockFGMRES(foam, cavity_mesh):
    """B9: FGMRES converges with a variable preconditioner; standard GMRES
    on the same set-up is only recorded (allowed to stall)."""
    d = _app(["Test-blockFGMRES"], cavity_mesh, 1,
             cavity_mesh / "fgmres.json")
    passed = (d["rc"] == 0 and bool(d.get("pass"))
              and bool(d.get("fgmresConverged")))
    results.write("tests", "Test-blockFGMRES",
                  {"run": d,
                   "fgmres": d.get("fgmres"),
                   "gmres": d.get("gmres"),
                   "gmresFixed": d.get("gmresFixed"),
                   "gmresOutcome": d.get("gmresOutcome"),
                   "pass": passed})
    assert passed, d


# --- heavy: motorBike mesh (350 k cells) ----------------------------------

def _motorbike_mesh(nprocs: int) -> Path:
    """Copy of the motorBike mesh (run/T4a_mesh, produced by
    cases/T4_motorBike 'Allrun -mesh a -mesh-only'), decomposed for nprocs.
    The source case is never modified. Skips if the mesh does not exist."""
    src = cfenv.RUN_ROOT / MOTORBIKE_MESH_CASE
    poly = src / "constant" / "polyMesh"
    owner = [p for p in (poly / "owner", poly / "owner.gz") if p.exists()]
    if not owner:
        pytest.skip(f"motorBike mesh not available: {poly}/owner missing "
                    "(produce it with cases/T4_motorBike: "
                    f"Allrun -mesh a -mesh-only into run/{MOTORBIKE_MESH_CASE})")
    cfenv.foam_env()
    dst = cfenv.RUN_ROOT / f"unit_motorBike_np{nprocs}"
    stamp = dst / ".meshSource"
    src_stamp = str(owner[0].stat().st_mtime)
    fresh = (dst.is_dir() and stamp.exists()
             and stamp.read_text() == src_stamp
             and (dst / f"processor{nprocs - 1}").is_dir()
             and not (dst / f"processor{nprocs}").exists())
    if fresh:
        return dst
    if dst.exists():
        shutil.rmtree(dst)
    dst.mkdir(parents=True)
    ignore = shutil.ignore_patterns("triSurface", "extendedFeatureEdgeMesh",
                                    "*.eMesh", "log.*")
    shutil.copytree(src / "system", dst / "system", symlinks=True)
    shutil.copytree(src / "constant", dst / "constant", symlinks=True,
                    ignore=ignore)
    cfcase.set_entry(dst, "system/decomposeParDict", "numberOfSubdomains",
                     nprocs)
    cfcase.set_entry(dst, "system/decomposeParDict", "method", "scotch")
    rc = cfenv.run(["decomposePar", "-force"], cwd=dst,
                   log=dst / "log.decomposePar")
    assert rc == 0, f"decomposePar failed, see {dst / 'log.decomposePar'}"
    stamp.write_text(src_stamp)
    return dst


@pytest.mark.heavy
def test_blockGAMG_cycles_motorBike(foam):
    """B9: the cycle study on the motorBike-tutorial mesh (350 k cells)."""
    case = _motorbike_mesh(4)
    rec = _cycle_study(case, "cyc")
    results.write("tests", "Test-blockGAMG_cycles_motorBike", rec)
    assert rec["pass"], rec["failures"]


def _parse_ranks_per_level(log: Path) -> list[int] | None:
    """First line starting with 'blockGAMG: ranks per level' -> integers."""
    if not log.exists():
        return None
    for line in log.read_text(errors="replace").splitlines():
        s = line.strip()
        if s.startswith(RANKS_LINE):
            return [int(v) for v in re.findall(r"\d+", s[len(RANKS_LINE):])]
    return None


def _expected_ranks(cells_per_level: list[int], nranks: int) -> list[int]:
    """Rule 6.3.4 for levels >= 1; level 0 is never agglomerated (native
    restriction, D-030). nCells(l) < 5000 -> 1 rank; nCells(l) <
    5000*nRanks -> max(1, nRanks/4) ranks; otherwise all ranks."""
    out = []
    for i, n in enumerate(cells_per_level):
        if i == 0:
            out.append(nranks)
        elif n < PROC_AGGLOM_SINGLE_RANK:
            out.append(1)
        elif n < PROC_AGGLOM_CELLS_PER_RANK * nranks:
            out.append(max(1, nranks // 4))
        else:
            out.append(nranks)
    return out


@pytest.mark.heavy
def test_procAgglom(foam):
    """B9 Test-procAgglom: CF_HEAVY_NP ranks, motorBike mesh; processor
    agglomeration on vs off identical to 1e-5; ranks per level follow
    rule 6.3.4 (level 0 excepted, D-030)."""
    n = PROC_AGGLOM_NPROCS
    case = _motorbike_mesh(n)
    common = ["-tolerance", str(GAMG_TOLERANCE),
              "-mergeLevels", str(MERGE_LEVELS),
              "-cycle", "K", "-skipDiagonal"]
    runs = {}
    for mode in ("on", "off"):
        d = _app(["Test-blockGAMG"], case, n, case / f"pa_{mode}.json",
                 common + ["-procAgglom", mode,
                           "-dumpSolution", str(case / f"pa_{mode}.xyz")])
        d["processWallSeconds"] = cfenv.last_timing.wall
        runs[mode] = d

    failures = []
    for mode, d in runs.items():
        if not d.get("converged"):
            failures.append(f"procAgglom {mode} did not converge (rc {d['rc']})")

    fon, foff = case / "pa_on.xyz", case / "pa_off.xyz"
    if fon.exists() and foff.exists():
        cross = _compare_solutions(foff, fon)
    else:
        cross = {"centresMatch": False, "relDiff": None, "maxAbsDiff": None}
    if not cross["centresMatch"]:
        failures.append("solution dumps missing or cell centres do not match")
    elif cross["relDiff"] >= TOL_CROSS_GAMG:
        failures.append(f"on vs off rel. difference {cross['relDiff']:.3e} "
                        f">= {TOL_CROSS_GAMG}")

    ranks = _parse_ranks_per_level(case / "pa_on.log")
    cells = runs["on"].get("gamgCellsPerLevel") or []
    expected = _expected_ranks(cells, n) if cells else None
    if ranks is None:
        failures.append(f"no '{RANKS_LINE}' line in {case / 'pa_on.log'}")
    elif not cells:
        failures.append("no gamgCellsPerLevel in the JSON of the 'on' run")
    elif ranks != expected:
        failures.append(f"ranks per level {ranks} != rule 6.3.4 {expected} "
                        f"(cells per level {cells})")

    rec = {"nProcs": n, "mesh": str(case), "runs": runs,
           "ranksPerLevel": ranks, "ranksPerLevelExpected": expected,
           "cellsPerLevel": cells,
           "rule": "nCells<5000 -> 1; nCells<5000*nRanks -> max(1,nRanks/4); "
                   "else nRanks",
           "solutionDiff": cross,
           "solutionMetric": "||x_off-x_on||_2/||x_off||_2, cells matched by centre",
           "tol": TOL_CROSS_GAMG,
           "failures": failures, "pass": not failures}
    results.write("tests", "Test-procAgglom", rec)
    assert rec["pass"], failures

"""Phase A / B2 unit tests (spec 6.4): the Test-* applications on 1 and 4
ranks, with cross-rank comparison.

  Test-block4Ops     ||A A^-1 - I||_inf < 1e-5, singular block guarded
  Test-doubleReduce  double sum exact to 1e-9 (D-004), float sum wrong
  Test-blockMatrix   block Amul vs native Amul < 1e-6 ||x||; 1 vs 4 ranks 1e-6
  Test-blockGAMG     converged to 1e-8 in <= 20 iterations; 1 vs 4 ranks 1e-5
"""

from __future__ import annotations

import json
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

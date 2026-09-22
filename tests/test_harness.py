"""Unit tests of the test/benchmark harness itself (no solver runs, no
OpenFOAM environment needed): D-068 fixes of the harness review.

Run: pytest tests/test_harness.py
"""

from __future__ import annotations

import math
import sys

import pytest

from cflib import env as cfenv

sys.path.insert(0, str(cfenv.REPO / "bench"))
import run_bench  # noqa: E402


def _wake(n: int, settle: int = 150, amp: float = 0.004, period: int = 37,
          mean: float = 0.4, cl: float = 0.07) -> dict[str, list[float]]:
    """Synthetic wake history: a start-up transient that decays by
    `settle`, then a stationary oscillation."""
    cd, cl_ = [], []
    for i in range(1, n + 1):
        tr = math.exp(-i / (settle / 5.0))
        cd.append(mean + 0.2 * tr + amp * math.sin(2 * math.pi * i / period))
        cl_.append(cl - 0.1 * tr + 3 * amp * math.cos(2 * math.pi * i / period))
    return {"Cd": cd, "Cl": cl_}


# --------------------------------------------------------------------------- #
# C1: one averaging window per wake case (D-068)
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("case", ["T4a", "T4b", "T5"])
def test_case_window_rule(case):
    spec = run_bench.CASES[case]
    assert spec.get("oscillatory")
    rule = max(run_bench.STAT_WINDOW_MIN,
               spec["iters"]["coupledFoam"] // run_bench.STAT_WINDOW_DIV)
    assert spec["statWindow"] == rule == 400
    for solver in ("simpleFoam", "coupledFoam"):
        n = spec["iters"][solver]
        assert run_bench.stat_window(n, case) == 400
    # capped at the iterations run
    assert run_bench.stat_window(250, case) == 250


def test_per_run_window_unchanged_without_case():
    assert run_bench.stat_window(3000) == 1500
    assert run_bench.stat_window(800) == 400
    assert run_bench.stat_window(400) == 300
    assert run_bench.stat_window(200) == 200


@pytest.mark.parametrize("name, key", [
    ("T4a", "T4a"), ("T4a_np10", "T4a"), ("ref_T4a_np10", "T4a"),
    ("bench_T4b_C_1", "T4b"), ("T4b_H_2", "T4b"),
    ("T5_coarse_np10", "T5"), ("ref_T5_coarse_np10", "T5"), ("T5_np10", "T5"),
    ("T3_kOmegaSST_np1", "T3-SST"), ("ref_T3_GEKO", "T3-GEKO"),
    ("bench_T3-SST_E-noSFD_1", "T3-SST"), ("T1_np4", "T1"), ("ref_T2", "T2"),
    ("T0_Re100_np1", None), ("T-fpe_linux64GccDPInt32Opt", None), (None, None)])
def test_case_of_run(name, key):
    assert run_bench.case_of_run(name) == key


def test_run_name_windows_equal_case_window():
    assert run_bench.stat_window(3000, "ref_T4a_np10") == 400
    assert run_bench.stat_window(800, "T4a_np10") == 400
    # no common window: per-run rule
    assert run_bench.stat_window(3000, "T3_kOmegaSST_np1") == 1500


def test_field_average_start_common_window():
    assert run_bench.field_average_start(800, "T4a") == 401
    assert run_bench.field_average_start(3000, "T4a") == 2601
    assert run_bench.field_average_start(3000) == 1501


def test_common_window_independent_of_budget():
    """The same flow history evaluated with two budgets gives the same
    convergence point under the common window; under the per-run window the
    longer budget can only converge after n/2 (the review C1 artefact)."""
    short, long_ = _wake(800), _wake(3000)
    n_short = run_bench.iters_to_stationary(short, case="T4a")
    n_long = run_bench.iters_to_stationary(long_, case="T4a")
    assert n_short is not None and n_short == n_long
    assert n_short <= 550
    n_long_run = run_bench.iters_to_stationary(long_)
    assert n_long_run is not None and n_long_run >= 1500


def test_stationary_eval_records_both_windows():
    h = _wake(3000)
    st = run_bench.stationary_eval(h, case="T4a")
    assert st["W"] == 400 and st["W_perRun"] == 1500
    assert st["windowRule"].startswith("common")
    assert st["stationary"] and st["stationary_perRun"]
    assert st["iters_to_stationary"] < st["iters_to_stationary_perRun"]
    assert abs(st["Cd_mean"] - 0.4) < 1e-3
    assert abs(st["Cd_mean_perRun"] - 0.4) < 1e-3
    # window statistics over exactly the last W samples
    tail = h["Cd"][-400:]
    assert abs(st["Cd_mean"] - sum(tail) / 400) < 1e-12
    st0 = run_bench.stationary_eval(h)
    assert st0["W"] == 1500 and st0["windowRule"].startswith("per run")


def test_config_hash_depends_on_case_window(monkeypatch):
    h0 = run_bench.config_hash("T4a", "C")
    monkeypatch.setitem(run_bench.CASES["T4a"], "statWindow", 500)
    assert run_bench.config_hash("T4a", "C") != h0

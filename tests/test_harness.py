"""Unit tests of the test/benchmark harness itself (no solver runs, no
OpenFOAM environment needed): D-068 fixes of the harness review.

Run: pytest tests/test_harness.py
"""

from __future__ import annotations

import json
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


# --------------------------------------------------------------------------- #
# C2: every benchmark configuration is a distinct run (D-068)
# --------------------------------------------------------------------------- #

def _in_scope(case):
    return [c for c in run_bench.CONFIGS if run_bench.in_scope(case, c)]


@pytest.mark.parametrize("case", sorted(run_bench.CASES))
def test_configurations_distinct(case):
    """All configurations in the scope of a case differ in the settings the
    solver actually uses (template + sets), and in their config hash."""
    cfgs = _in_scope(case)
    eff = {c: run_bench.effective_settings(case, c) for c in cfgs}
    for i, a in enumerate(cfgs):
        for b in cfgs[i + 1:]:
            assert eff[a] != eff[b], f"{case}: {a} and {b} are the same run"
    hashes = {run_bench.config_hash(case, c) for c in cfgs}
    assert len(hashes) == len(cfgs)


def test_configuration_semantics():
    """H is a fixed K-cycle, H-tune the K-cycle with the controller, C the
    template (V-cycle, autoTune off, D-043) on every case."""
    for case in run_bench.CASES:
        c = run_bench.effective_settings(case, "C")["solvers.coupled"]["blockGAMG"]
        h = run_bench.effective_settings(case, "H")["solvers.coupled"]["blockGAMG"]
        ht = run_bench.effective_settings(case, "H-tune")["solvers.coupled"]["blockGAMG"]
        assert (c["cycleType"], c["autoTune"]) == ("V", False), case
        assert (h["cycleType"], h["autoTune"]) == ("K", False), case
        assert (ht["cycleType"], ht["autoTune"]) == ("K", True), case


def test_t1_native_configurations_identical():
    """T1: the pitzDaily tutorial is SIMPLEC with p unrelaxed and U, k,
    omega 0.9, i.e. configuration B: B is out of the T1 scope."""
    assert (run_bench.effective_settings("T1", "A")
            == run_bench.effective_settings("T1", "B"))
    assert not run_bench.in_scope("T1", "B")


def test_t3_sfd_variant_is_not_a_noop():
    for case in ("T3-SST", "T3-GEKO"):
        c = run_bench.effective_settings(case, "C")["coupled"]["sfd"]["enabled"]
        v = run_bench.effective_settings(case, "E-noSFD")["coupled"]["sfd"]["enabled"]
        assert c is True and v is False


def test_scope_d063():
    assert set(_in_scope("T4a")) == {"A", "B", "C", "H"}
    assert _in_scope("T4b") == ["C"]
    assert _in_scope("T5") == []
    assert run_bench.MAX_REPEATS["T4a"] == run_bench.MAX_REPEATS["T4b"] == 1


def test_alias_e_is_c():
    cfgs, notes = run_bench.resolve_configs(["C", "E", "H"])
    assert cfgs == ["C", "H"] and notes
    assert "E" not in run_bench.CONFIGS


# --------------------------------------------------------------------------- #
# M3: failed runs are detected, recorded as failed and rerun
# --------------------------------------------------------------------------- #

from cflib import case as cfcase  # noqa: E402
from cflib import results as cfresults  # noqa: E402

TIME_REPORT = """\tCommand being timed: "x"
\tUser time (seconds): {user}
\tSystem time (seconds): 0.50
\tElapsed (wall clock) time (h:mm:ss or m:ss): 0:{wall:05.2f}
\tMaximum resident set size (kbytes): 102400
"""


def _timing(tmp_path, solver, walls, incomplete=()):
    t = tmp_path / "timing"
    t.mkdir(exist_ok=True)
    for i, w in enumerate(walls):
        f = t / f"{solver}.rank{i}.time"
        if i in incomplete:
            f.write_text('\tCommand being timed: "x"\n')     # killed rank
        else:
            f.write_text(TIME_REPORT.format(user=w, wall=w))
    return tmp_path


def test_rank_times_all_incomplete_no_keyerror(tmp_path):
    case = _timing(tmp_path, "coupledFoam", [10, 11], incomplete=(0, 1))
    rt = run_bench.rank_times(case, "coupledFoam")
    assert rt["complete"] is False and "wallSeconds" not in rt
    assert len(rt["incompleteReports"]) == 2
    assert run_bench.to_convergence(rt, 0.5) == {}


def test_rank_times_partial(tmp_path):
    case = _timing(tmp_path, "coupledFoam", [10, 12, 11], incomplete=(2,))
    _timing(tmp_path, "potentialFoam", [2, 2])
    rt = run_bench.rank_times(case, "coupledFoam")
    assert rt["complete"] is False
    assert rt["solverWallSeconds"] == pytest.approx(12)
    assert rt["wallSeconds"] == pytest.approx(14)
    (tmp_path / "b").mkdir()
    full = run_bench.rank_times(_timing(tmp_path / "b", "coupledFoam", [5]),
                                "coupledFoam")
    assert full["complete"] is True and full["incompleteReports"] == []


def test_run_failure(tmp_path):
    assert "Allrun rc 1" in cfcase.run_failure(tmp_path, "coupledFoam", 1)
    assert "no log.coupledFoam" in cfcase.run_failure(tmp_path, "coupledFoam", 1)
    (tmp_path / "log.coupledFoam").write_text("CF| iter=1\n")
    r = cfcase.run_failure(tmp_path, "coupledFoam", 0)
    assert r == ["log.coupledFoam has no normal end"]
    (tmp_path / "log.coupledFoam").write_text("CF| iter=1\nEnd\n")
    assert cfcase.run_failure(tmp_path, "coupledFoam", 0) == []
    assert cfcase.run_failure(tmp_path, "coupledFoam", 0, 500, 800) == \
        ["500 of 800 iterations"]
    (tmp_path / "log.coupledFoam").write_text("--> FOAM FATAL ERROR\nEnd\n")
    assert cfcase.run_failure(tmp_path, "coupledFoam", 0)
    (tmp_path / "log.Allrun").write_text("mpirun: invalid cpu-set\n")
    assert "invalid cpu-set" in cfcase.log_tail(tmp_path, "coupledFoam")


@pytest.fixture
def bench_results(tmp_path, monkeypatch):
    monkeypatch.setattr(cfresults, "RESULTS", tmp_path)
    (tmp_path / "bench").mkdir()
    return tmp_path / "bench"


def _write_bench(d, case, cfg, run, **kw):
    rec = {"case": case, "config": cfg, "run": run,
           "configHash": run_bench.config_hash(case, cfg)}
    rec.update(kw)
    (d / f"{case}_{cfg}_{run}.json").write_text(json.dumps(rec))
    return rec


def test_failed_records_excluded_and_listed(bench_results):
    _write_bench(bench_results, "T1", "C", 1, failed=False,
                 wall_to_conv_s=10.0, cpu_to_conv_h=0.01, iters_to_conv=300)
    _write_bench(bench_results, "T1", "A", 1, failed=True, rc=1,
                 failure=["Allrun rc 1"])
    _write_bench(bench_results, "T1", "H", 1, rc=1, error="monitor: x")  # old style
    recs, stale = run_bench.load_current(["T1"])
    assert [r["config"] for r in recs] == ["C"] and not stale
    failed = {r["config"] for r in run_bench.load_failed(["T1"])}
    assert failed == {"A", "H"}
    out = run_bench.summarise(["T1"], ["A", "C", "H"], repeats=1)
    s = json.loads(out.with_name("summary.json").read_text())
    assert {f["tag"] for f in s["failed"]} == {"T1_A_1", "T1_H_1"}
    assert set(s["missing"]) == {"T1_A_1", "T1_H_1"}
    assert [t["config"] for t in s["table"]] == ["C"]


def test_failed_record_is_rerun(bench_results, monkeypatch):
    class _Stop(Exception):
        pass

    def _prepare(*a, **k):
        raise _Stop

    class _State:
        def as_dict(self):
            return {}

    monkeypatch.setattr(run_bench.cfcase, "prepare", _prepare)
    monkeypatch.setattr(run_bench.cfenv, "machine_state", lambda: _State())
    ok = _write_bench(bench_results, "T1", "C", 1, failed=False,
                      wall_to_conv_s=1.0)
    assert run_bench.run_one("T1", "C", 1, 1, force=False) == ok
    _write_bench(bench_results, "T1", "C", 2, failed=True,
                 failure=["Allrun rc 1"])
    with pytest.raises(_Stop):
        run_bench.run_one("T1", "C", 2, 1, force=False)
    # the failed record was moved aside, not kept as a result
    assert not (bench_results / "T1_C_2.json").exists()


# --------------------------------------------------------------------------- #
# M1: the speed-up times both solvers to the same criterion
# --------------------------------------------------------------------------- #

def _make_report():
    pytest.importorskip("matplotlib")
    import make_report  # noqa: PLC0415
    return make_report


def test_speed_criterion_per_case():
    mr = _make_report()
    assert mr.speed_criterion("T0_Re100_np1", {"Rtarget": 1e-8}) == ("residual", 1e-8)
    assert mr.speed_criterion("T1_np1", {"Rtarget": 1e-5}) == ("residual", 1e-5)
    assert mr.speed_criterion("T2_np1", {}) == ("residual", 1e-5)
    assert mr.speed_criterion("T3_GEKO_np1", {})[0] == "forceWindow"
    assert mr.speed_criterion("T4a_np10", {})[0] == "stationary"
    assert mr.speed_criterion("T5_coarse_np10", {})[0] == "stationary"


def test_simplefoam_residual_iteration_all_fields():
    """simpleFoam converges when EVERY initial residual is below R (not at
    its residualControl stop)."""
    import numpy as np  # noqa: PLC0415
    mr = _make_report()
    n = 100
    sft = {"t": np.arange(1, n + 1, dtype=float), "n": n,
           "res": {"p": np.logspace(-2, -9, n), "Ux": np.logspace(-3, -9, n),
                   "k": np.logspace(-1, -6, n)}}
    it = mr._sf_residual_iteration(sft, 1e-5)
    k = sft["res"]["k"]
    assert it == int(np.nonzero(k < 1e-5)[0][0]) + 1
    assert mr._sf_residual_iteration(sft, 1e-12) is None
    assert mr._frac(sft, it, n) == pytest.approx(it / n)


# --------------------------------------------------------------------------- #
# M2: fairness flags of the speed-up against a cached reference
# --------------------------------------------------------------------------- #

def test_reference_timing_flags(tmp_path):
    f = run_bench.reference_timing_flags(tmp_path, {})
    assert f["referenceNoPotentialStart"] and f["referenceTimingConditionsUnknown"]
    (tmp_path / "log.potentialFoam").write_text("End\n")
    f = run_bench.reference_timing_flags(tmp_path, {"machineBefore": {"loadavg": 1}})
    assert not f["referenceNoPotentialStart"]
    assert not f["referenceTimingConditionsUnknown"]


# --------------------------------------------------------------------------- #
# M6: T3 requires convergence and compares window means
# --------------------------------------------------------------------------- #

def _force_case(tmp_path, cd, cl):
    d = tmp_path / "postProcessing" / "forceCoeffs" / "0"
    d.mkdir(parents=True)
    rows = ["# Time Cd Cl"] + [f"{i + 1} {a!r} {b!r}"
                                for i, (a, b) in enumerate(zip(cd, cl))]
    (d / "coefficient.dat").write_text("\n".join(rows) + "\n")
    return tmp_path


def test_t3_limit_cycle_does_not_pass_on_last_sample(tmp_path):
    import test_T3_airFoil as t3  # noqa: PLC0415
    n = 600
    # limit cycle around 0.12 +- 0.1; the last sample is exactly 0.10
    cd = [0.12 + 0.1 * math.sin(2 * math.pi * i / 40) for i in range(n - 1)]
    cd.append(0.10)
    cl = [0.5] * n
    r = t3.converged_coeffs(_force_case(tmp_path, cd, cl), solver_stop=False)
    assert r["converged"] is False and r["itersToConv"] is None
    assert r["CdLast"] == 0.10 and abs(r["Cd"] - 0.10) > 1e-3


def test_t3_converged_run(tmp_path):
    import test_T3_airFoil as t3  # noqa: PLC0415
    n = 500
    cd = [0.09 + 0.05 * math.exp(-i / 30) for i in range(n)]
    cl = [0.25 - 0.1 * math.exp(-i / 30) for i in range(n)]
    r = t3.converged_coeffs(_force_case(tmp_path, cd, cl), solver_stop=False)
    assert r["converged"] and r["finalWindowOk"]
    assert r["itersToConv"] == run_bench.iters_to_conv({"Cd": cd, "Cl": cl})
    assert r["Cd"] == pytest.approx(sum(cd[-100:]) / 100)
    # the solver's own stop counts as converged even without the window
    r2 = t3.converged_coeffs(_force_case(tmp_path / "s", cd[:120], cl[:120]),
                             solver_stop=True)
    assert r2["converged"] and r2["itersToConv"] == 120


# --------------------------------------------------------------------------- #
# M7: a continued reference is read up to its original budget only
# --------------------------------------------------------------------------- #

def test_reference_continuation_excluded(tmp_path):
    import user_convergence as ucv  # noqa: PLC0415
    n, n0 = 450, 300
    case = _force_case(tmp_path / "ref_T4a_np10", [0.4] * n, [0.07] * n)
    (case / "reference.json").write_text(json.dumps(
        {"continuation": {"startTime": float(n0), "ok": True}}))
    assert run_bench.reference_t_max(case) == n0
    assert len(run_bench.force_history(case)["Cd"]) == n0
    assert len(run_bench.force_history(case, original_only=False)["Cd"]) == n
    assert len(ucv.force_hist(case)["Cd"]) == n0
    ev = ucv.evaluate(case, "simpleFoam", 400, {})
    assert "ignored" in ev and "original reference budget" in ev["ignored"]
    # a run without continuation is not cut
    other = _force_case(tmp_path / "T4a_np10", [0.4] * n, [0.07] * n)
    assert len(run_bench.force_history(other)["Cd"]) == n

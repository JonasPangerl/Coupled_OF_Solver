"""Deep diagnostics logging (OPUS_TASKS TASK 5.5, DECISIONS D-045). Quick.

T0 Re100, maxIter 20, coupled.diagnostics.level 3, np1 and np4:
- every line of every diagnostics/diag.rank*.jsonl parses as JSON
- one file per rank, a header record, one "iter" record per iteration
- all level-1 keys of TASK 5.3 are present, tDiag is recorded (>= 0)
- level-2 ("linear") and level-3 ("gamgSetup", "bcFlips") content present
- global fields identical on every rank; bench/diag_tools.load merges the
  rank files
- the diagnostics only read the solver state: the CF| log lines (without
  the timing fields) equal those of a level-0 run (asserted for np1; T0
  np4 is not reproducible run-to-run even at level 0, see D-045)
With level 0 (default, no diagnostics block): no diagnostics/ directory.
"""

from __future__ import annotations

import json
import subprocess
import sys

import pytest

from cflib import case as cfcase
from cflib import env as cfenv
from cflib import logs, results

sys.path.insert(0, str(cfenv.REPO / "bench"))
import diag_tools  # noqa: E402

MAX_ITERS = 20
LEVEL = 3
PHASES = {"startup", "stalled", "asymptotic", "ramp"}

# TASK 5.3, level 1 (dotted paths into the record)
REQUIRED = [
    "iter", "phase", "wallTime",
    "residuals.R", "residuals.rU", "residuals.rp", "residuals.R1",
    "residuals.normFactor", "residuals.linInitial", "residuals.linFinal",
    "residuals.linIts", "residuals.linConverged", "residuals.rho",
    "residuals.massErrMax", "residuals.massErrSum",
    "controls.CFL", "controls.dt.min", "controls.dt.median",
    "controls.dt.max", "controls.nLocLim", "controls.strategy",
    "controls.hold", "controls.growth", "controls.eta", "controls.etaRaw",
    "controls.etaClip", "controls.omega", "controls.cuts", "controls.trials",
    "controls.sentinel.checks", "controls.sentinel.nRollbacks",
    "controls.remediation.nStat", "controls.remediation.nDyn",
    "controls.remediation.version", "controls.anderson.enabled",
    "controls.beta",
    "turbulence.k", "turbulence.omega", "turbulence.nBoundK",
    "turbulence.nBoundOmega", "turbulence.nNutCapped", "turbulence.nClamped",
    "timings.tAsm", "timings.tMomentumOps", "timings.tContinuity",
    "timings.tBoundary", "timings.tRhieChow", "timings.tSolve",
    "timings.tPrecSetup", "timings.tPrecApply", "timings.tKrylov",
    "timings.tTurb", "timings.tDiag", "timings.tIter", "timings.tWall",
    "gamg.nLevels", "gamg.Cop", "gamg.hierarchy", "gamg.reagglomerated",
    "memory.rssKB", "memory.peakRssKB",
]
# identical on all ranks (reduced by the solver)
GLOBAL = ["residuals.R", "residuals.linIts", "controls.CFL", "controls.eta",
          "controls.omega", "phase", "gamg.nLevels"]
# CF| fields that are timings (differ between runs)
CF_TIMING = {"tAsm", "tSolve", "tTurb", "tIter", "tWall"}


def _get(rec: dict, path: str):
    cur = rec
    for k in path.split("."):
        if not isinstance(cur, dict) or k not in cur:
            raise KeyError(path)
        cur = cur[k]
    return cur


def _set(case, entry: str, value) -> None:
    """foamDictionary -disableFunctionEntries -set on system/fvSolution"""
    subprocess.run(
        ["foamDictionary", "-disableFunctionEntries", "-entry", entry,
         "-set", str(value), "system/fvSolution"],
        cwd=case, check=True, stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )


def _set_diag(case, level: int) -> None:
    _set(case, "coupled/diagnostics", f"{{ level {level}; }}")


def _run(name: str, nprocs: int, level: int | None):
    case = cfcase.prepare("T0_cavity", name)
    _set(case, "coupled/maxIter", MAX_ITERS)
    if level is not None:
        _set_diag(case, level)
    rc = cfcase.allrun(case, ["-solver", "coupledFoam", "-Re", "100",
                              "-np", str(nprocs)], fpe=True,
                       extra_env={"CF_MPI_BIND": "none"})
    return case, rc


def _cf_rows(case):
    rows = logs.parse_cf(case / "log.coupledFoam")
    return [{k: v for k, v in r.items() if k not in CF_TIMING} for r in rows]


@pytest.mark.case
def test_diagnostics_level0(foam):
    """Default (no diagnostics block): nothing is written."""
    case, rc = _run("diag_T0_l0_np1", 1, None)
    rec = {"rc": rc, "diagnosticsDir": (case / "diagnostics").exists()}
    rec["pass"] = rc == 0 and not rec["diagnosticsDir"]
    results.write("tests", "diagnostics_level0", rec)
    assert rc == 0
    assert not (case / "diagnostics").exists(), \
        "level 0 must not create diagnostics/"


@pytest.mark.case
def test_diagnostics(foam, nprocs):
    name = f"diag_T0_l{LEVEL}_np{nprocs}"
    case, rc = _run(name, nprocs, LEVEL)
    assert rc == 0, "coupledFoam failed"
    assert not logs.fpe_trapped(case / "log.coupledFoam"), "FPE trap"

    files = sorted((case / "diagnostics").glob("diag.rank*.jsonl"))
    assert len(files) == nprocs, [f.name for f in files]

    n_lines = 0
    per_rank = {}
    missing = set()
    for f in files:
        lines = f.read_text().splitlines()
        parsed = []
        for i, line in enumerate(lines):
            try:
                parsed.append(json.loads(line))
            except json.JSONDecodeError as err:
                pytest.fail(f"{f.name}:{i + 1} is not JSON: {err}")
        n_lines += len(parsed)
        assert parsed[0]["type"] == "header"
        assert parsed[0]["level"] == LEVEL
        assert not any(r.get("type") == "truncated" for r in parsed)
        its = [r for r in parsed if r.get("type") == "iter"]
        assert [r["iter"] for r in its] == list(range(1, MAX_ITERS + 1))
        for r in its:
            for k in REQUIRED:
                try:
                    _get(r, k)
                except KeyError:
                    missing.add(k)
            assert r["phase"] in PHASES, r["phase"]
            assert r["timings"]["tDiag"] >= 0
            assert isinstance(r["linear"], list) and r["linear"], "level 2"
            s = r["linear"][0]
            assert s["krylov"] and s["precon"], "level 2 arrays"
            assert s["precon"][0]["levels"], "per-level norms"
            assert "gamgSetup" in r and "bcFlips" in r, "level 3"
        # full hierarchy in the first record, "unchanged" afterwards
        assert isinstance(its[0]["gamg"]["hierarchy"], dict)
        assert "cellsPerLevel" in its[0]["gamg"]["hierarchy"]
        per_rank[f.name] = its
    assert not missing, f"missing level-1 keys: {sorted(missing)}"

    # global fields identical on every rank
    ranks = list(per_rank.values())
    for its in ranks[1:]:
        for a, b in zip(ranks[0], its):
            for k in GLOBAL:
                assert _get(a, k) == _get(b, k), (k, a["iter"])

    # analysis companion
    df = diag_tools.load(case)
    assert len(df) == MAX_ITERS and int(df["nRanks"].iloc[0]) == nprocs
    hist = diag_tools.linear_history(case, MAX_ITERS)
    assert hist and hist[0]["krylov"].size == hist[0]["its"]
    red = diag_tools.level_reduction(case)
    assert not red.empty

    # diagnostics only read the state: CF| lines equal to a level-0 run.
    # Asserted serially only: T0 np4 is not reproducible run-to-run even at
    # level 0 with the pre-TASK-5 binary (iteration-1 linIters 7/9/9 in
    # three identical runs, 2026-09-22, D-045), so the comparison is only
    # recorded there.
    ref, rc0 = _run(f"diag_T0_l0cmp_np{nprocs}", nprocs, None)
    assert rc0 == 0
    same = _cf_rows(ref) == _cf_rows(case)

    tdiag = [r["timings"]["tDiag"] for r in ranks[0]]
    titer = [r["timings"]["tIter"] for r in ranks[0]]
    sizes = {f.name: f.stat().st_size for f in files}
    rec = {
        "nProcs": nprocs, "level": LEVEL, "iterations": MAX_ITERS,
        "lines": n_lines, "bytesPerRank": sizes,
        "tDiagSum": sum(tdiag), "tIterSum": sum(titer),
        "phases": df["phase"].value_counts().to_dict(),
        "cfLogIdenticalToLevel0": same,
    }
    rec["pass"] = bool(same) or nprocs > 1
    results.write("tests", f"diagnostics_np{nprocs}", rec)
    if nprocs == 1:
        assert same, "diagnostics changed the solver trajectory (CF| lines)"


# * * * * * * * * * * Force-coefficient window statistics * * * * * * * * * //
# (D-045 b: window mean/RMS on the CF| line, in the diagnostics record and in
# summary.json; optional drift stop rule, off by default)

FORCE_ITERS = 30
RESTART_ITERS = 5


@pytest.mark.unit
def test_convergence_monitor_unit(foam, tmp_path):
    """Synthetic sequences: RMS 0 / drift 0 exactly for a constant, RMS of
    a sinusoid, drift rule on/off, ramp, restart round trip."""
    out = tmp_path / "cm.json"
    rc = cfenv.run(["Test-convergenceMonitor", "-json", str(out)],
                   cwd=tmp_path, log=tmp_path / "cm.log")
    data = json.loads(out.read_text()) if out.exists() else {}
    data["rc"] = rc
    results.write("tests", "Test-convergenceMonitor", data)
    assert rc == 0 and data.get("pass"), (tmp_path / "cm.log").read_text()


@pytest.mark.case
def test_force_stats(foam):
    """T3 airFoil2D (kOmegaSST, forceCoeffs), np1, 30 iterations and a
    restart to 35: CdMean/CdRms/ClMean/ClRms on every CF| line (appended
    after Cl), the forces block in the diagnostics record, forceStats and
    forceHistory in summary.json; the window continues across the
    restart."""
    case = cfcase.prepare("T3_airFoil2D", "diag_T3_forces")
    _set(case, "coupled/maxIter", FORCE_ITERS)
    _set_diag(case, 1)
    args = ["-solver", "coupledFoam", "-turbulence", "kOmegaSST", "-np", "1"]
    rc = cfcase.allrun(case, args, fpe=True, extra_env={"CF_MPI_BIND": "none"})
    assert rc == 0, "coupledFoam failed"

    rows = logs.parse_cf(case / "log.coupledFoam")
    assert len(rows) == FORCE_ITERS
    # The forceCoeffs result of iteration k is first available to the
    # monitor in iteration k+1: rows with Cd carry the statistics
    rows = [r for r in rows if "Cd" in r]
    assert len(rows) >= FORCE_ITERS - 1
    keys = ["CdMean", "CdRms", "ClMean", "ClRms"]
    for r in rows:
        for k in keys:
            assert isinstance(r.get(k), (int, float)), (k, r)
        assert r["CdRms"] >= 0 and r["ClRms"] >= 0
    # appended at the end: the existing fields keep their order
    line = [ln for ln in (case / "log.coupledFoam").read_text().splitlines()
            if ln.startswith("CF|")][-1]
    toks = [t.split("=")[0] for t in line.split()[1:]]
    assert toks.index("Cl") < toks.index("CdMean") and toks[-1].endswith("Rms")

    summ = logs.coupled_summary(case)
    fs, fh = summ.get("forceStats"), summ.get("forceHistory")
    assert fs and fh, "summary.json without forceStats/forceHistory"
    assert fs["driftEnabled"] is False
    assert len(fh["iter"]) == len(rows) == len(fh["CdMean"]) == len(fh["ClRms"])
    assert abs(fs["Cd"]["mean"] - rows[-1]["CdMean"]) <= 1e-5 * abs(fs["Cd"]["mean"])

    recs = [r for r in diag_tools.records(case)[0] if r["forces"]]
    n_samples = [r["forces"]["Cd"]["n"] for r in recs]
    assert n_samples == list(range(1, len(rows) + 1)), n_samples

    # restart: the statistics window continues (restored samples + new)
    _set(case, "coupled/maxIter", FORCE_ITERS + RESTART_ITERS)
    rc2 = cfcase.allrun(case, args + ["-restart"], fpe=True,
                        extra_env={"CF_MPI_BIND": "none"})
    assert rc2 == 0, "restart failed"
    recs = diag_tools.records(case)[0]
    first = [r for r in recs if r["iter"] == FORCE_ITERS + 1][0]
    n_restart = first["forces"]["Cd"]["n"]
    # all samples of the first part were restored, plus the new one
    ok = n_restart == n_samples[-1] + 1
    results.write("tests", "diagnostics_force_stats", {
        "rows": len(rows), "finalForceStats": fs,
        "nBeforeRestart": n_samples[-1], "nAfterRestart": n_restart,
        "pass": ok,
    })
    assert ok, (n_samples[-1], n_restart)

"""User-judged convergence iterations of the force cases (D-060).

The automatic criteria (12.3(ii) for T3, the D-042 stationary mean for the
wake cases T4/T5) decide the pass/fail of the tests. For the numbers the
report quotes - iterations, wall time and CPU-hours to convergence, speed-ups
and the force coefficients - the user may instead name an earlier
iteration from which they consider a run converged, after looking at the
load histories. The runs keep their stop criteria (user, 2026-09-22), so a
user iteration can only move the convergence point earlier: an iteration
beyond the last one of the run is ignored and listed as such.

File: report/user_convergence.json (tracked; edited by hand)

    {
      "T4a": {"coupledFoam": 450, "simpleFoam": 1500, "C": 500, "A": 1600},
      "T3-SST": {"coupledFoam": 650, "simpleFoam": 4000},
      ...
    }

Case keys: T3-SST, T3-GEKO, T4a, T4b, T5. Within a case:
  "coupledFoam" / "simpleFoam"  the test run and its simpleFoam reference
  "A" ... "H", "H-tune", "E-..." one benchmark configuration (all repeats)
A missing or null entry keeps the automatic value. Keys starting with "_"
are comments.

A cached reference that was continued with averaging (ref_T4a, ref_T4b)
is read up to the end of its original budget only (M7): a user iteration
beyond it is ignored and flagged like one beyond the end of a run.

With a user iteration N of a run with n iterations:
  iterations to convergence  = N
  wall / CPU-h to convergence = pre-processing + solver total x fraction of
                                the solver wall time spent up to N
                                (run_bench.to_convergence)
  Cd, Cl                      = mean (and std) over iterations N..n
The automatic values are kept in rec["auto"]; rec["convergenceSource"] is
"user" or "auto".
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
FILE = REPO / "report" / "user_convergence.json"
RUN = REPO / "run"
CASES = ("T3-SST", "T3-GEKO", "T4a", "T4b", "T5")


def load(path: Path = FILE) -> dict:
    if not path.exists():
        return {}
    d = json.loads(path.read_text())
    return {k: {kk: vv for kk, vv in v.items() if not kk.startswith("_")}
            for k, v in d.items() if not k.startswith("_") and isinstance(v, dict)}


def case_key(test_name: str) -> str | None:
    """Test record name -> case key (T3_kOmegaSST_np1 -> T3-SST, ...)."""
    if test_name.startswith("T3_kOmegaSST"):
        return "T3-SST"
    if test_name.startswith("T3_GEKO"):
        return "T3-GEKO"
    for k in ("T4a", "T4b"):
        if test_name.startswith(k + "_"):
            return k
    if test_name.startswith("T5_"):
        return "T5"
    return None


def lookup(table: dict, case: str | None, key: str) -> int | None:
    v = (table.get(case) or {}).get(key) if case else None
    return int(v) if v is not None else None


def force_hist(case: Path) -> dict[str, np.ndarray]:
    """Cd / Cl history of a run. A continued simpleFoam reference (ref_*
    with reference.json "continuation") is cut at the end of its ORIGINAL
    budget (M7): the continuation samples are not part of the reference run,
    as in the test's own evaluation (t_max)."""
    import run_bench  # noqa: PLC0415
    h = run_bench.force_history(case)
    if not h:
        return {}
    return {"Cd": np.asarray(h["Cd"], float), "Cl": np.asarray(h["Cl"], float)}


def _rt(rec: dict) -> dict:
    """rank_times-like totals of a record (T3 records carry only the
    solver totals)."""
    rt = rec.get("ranks") or rec.get("total")
    if rt and "solverWallSeconds" in rt:
        return rt
    if rec.get("wallSecondsSolver") is not None:
        return {"solverWallSeconds": rec["wallSecondsSolver"],
                "solverCpuHours": rec.get("cpuHoursSolver") or 0.0}
    # the CPU part is flagged unavailable by evaluate() when missing
    return {}


def guard_ok(case: Path, rec: dict) -> bool:
    """The run directory may be read for record `rec` (M9): it passes the
    active staleness guard of the report (provenance at the target commit
    and build; run/ref_* exempt) and is still the run the record was made
    from. Always true without an active guard."""
    import run_bench  # noqa: PLC0415, F401  (puts tests/ on sys.path)
    from cflib import provenance  # noqa: PLC0415
    return provenance.active().check_record_run(Path(case).name, rec,
                                                Path(case))


def evaluate(case: Path, solver: str, it: int, rec: dict) -> dict:
    """Numbers of a run converged at iteration `it` (see module doc). A run
    directory rejected by the staleness guard (guard_ok) is not read: the
    result is flagged "ignored"."""
    import run_bench  # noqa: PLC0415
    out: dict = {"iters": it}
    if not guard_ok(case, rec):
        out["ignored"] = (f"run directory {Path(case).name} rejected by the "
                          "staleness guard (listed in the appendix)")
        return out
    h = force_hist(case)
    if h:
        n = min(len(h["Cd"]), len(h["Cl"]))
        out["iterationsRun"] = n
        if it > n:
            cont = run_bench.reference_t_max(case) is not None
            out["ignored"] = (f"iteration {it} > {n} iterations run"
                              + (" (the original reference budget; the "
                                 "continuation is not part of the "
                                 "reference, D-060)" if cont else ""))
            return out
        k = min(max(it, 1), n) - 1
        for q in ("Cd", "Cl"):
            w = h[q][k:n]
            out[q] = float(np.mean(w))
            out[q + "_std"] = float(np.std(w))
        out["samplesAfter"] = n - k
    frac = run_bench.progress_fraction(case, solver, it)
    out["progressFraction"] = frac
    rt = _rt(rec)
    tc = run_bench.to_convergence(rt, frac)
    no_cpu = rt.get("solverCpuHours") is None or (
        "ranks" not in rec and "total" not in rec
        and rec.get("cpuHoursSolver") is None)
    out["wall_to_conv_s"] = tc.get("wall_to_conv_s")
    out["cpu_to_conv_h"] = None if no_cpu else tc.get("cpu_to_conv_h")
    out["solver_wall_to_conv_s"] = tc.get("solver_wall_to_conv_s")
    out["solver_cpu_to_conv_h"] = None if no_cpu else tc.get("solver_cpu_to_conv_h")
    return out


# --------------------------------------------------------------------------- #
# benchmark records
# --------------------------------------------------------------------------- #

BENCH_KEYS = ("iters_to_conv", "wall_to_conv_s", "cpu_to_conv_h",
              "time_per_iter_s", "cpu_per_iter_s")


def apply_bench(rec: dict, table: dict) -> dict:
    case, cfg = rec.get("case"), rec.get("config")
    it = lookup(table, case, cfg)
    rec.setdefault("convergenceSource", "auto")
    if it is None:
        return rec
    d = RUN / f"bench_{case}_{cfg}_{rec.get('run', 1)}"
    if not d.is_dir():
        return rec
    ev = evaluate(d, rec.get("solver", "coupledFoam"), it, rec)
    if ev.get("ignored"):
        rec["userIgnored"] = ev["ignored"]
        return rec
    rec["auto"] ={k: rec.get(k) for k in BENCH_KEYS + ("Cd_mean", "Cl_mean")}
    rec["iters_to_conv"] = it
    rec["wall_to_conv_s"] = ev["wall_to_conv_s"]
    rec["cpu_to_conv_h"] = ev["cpu_to_conv_h"]
    if ev.get("solver_wall_to_conv_s"):
        rec["time_per_iter_s"] = ev["solver_wall_to_conv_s"] / it
    if ev.get("solver_cpu_to_conv_h"):
        rec["cpu_per_iter_s"] = ev["solver_cpu_to_conv_h"] * 3600.0 / it
    if "Cd" in ev:
        rec["Cd_mean"], rec["Cl_mean"] = ev["Cd"], ev["Cl"]
    rec["user"] = ev
    rec["convergenceSource"] = "user"
    return rec


# --------------------------------------------------------------------------- #
# test records (T3, T4a, T4b, T5 and their simpleFoam references)
# --------------------------------------------------------------------------- #

TEST_KEYS = ("itersToConv", "wallToConv_s", "cpuHoursToConv", "Cd", "Cl",
             "CdRelDiff", "ClRelDiff", "speedupWall", "speedupCpu")


def ref_dir(test_name: str) -> Path:
    if test_name.startswith("T3_"):
        return RUN / ("ref_" + test_name.rsplit("_np", 1)[0])
    return RUN / ("ref_" + test_name)


def auto_iteration(rec: dict, case: Path) -> int | None:
    """Automatic convergence iteration of a record (T3 records do not store
    one: criterion 12.3(ii) on the history)."""
    if rec.get("itersToConv") is not None:
        return rec["itersToConv"]
    if rec.get("convergedAt") is not None:
        return rec["convergedAt"]
    import run_bench  # noqa: PLC0415
    if not Path(case).is_dir() or not guard_ok(case, rec):
        return rec.get("iterations")
    h = force_hist(case)
    it = (run_bench.iters_to_conv({k: v.tolist() for k, v in h.items()})
          if h else None)
    # T3: the solver stops on its own criterion (residual / force window),
    # so that stop IS the automatic convergence point
    return it if it is not None else rec.get("iterations")


def apply_test(name: str, rec: dict, table: dict) -> dict:
    case = case_key(name)
    rec.setdefault("convergenceSource", "auto")
    if case is None:
        return rec
    cdir, rdir = RUN / name, ref_dir(name)
    ref = rec.get("reference") or {}
    rec.setdefault("autoIteration", auto_iteration(rec, cdir))
    ref.setdefault("autoIteration", auto_iteration(ref, rdir))
    itc = lookup(table, case, "coupledFoam")
    its = lookup(table, case, "simpleFoam")
    if itc is None and its is None:
        return rec
    rec["auto"] = {k: rec.get(k) for k in TEST_KEYS}
    rec["auto"]["iteration"] = rec["autoIteration"]
    if itc is not None and cdir.is_dir():
        ev = evaluate(cdir, "coupledFoam", itc, rec)
        if ev.get("ignored"):
            rec["userIgnored"] = ev["ignored"]
            itc = None
    if itc is not None and cdir.is_dir():
        rec["user"] = ev
        rec["itersToConv"] = itc
        rec["wallToConv_s"] = ev["wall_to_conv_s"]
        rec["cpuHoursToConv"] = ev["cpu_to_conv_h"]
        if "Cd" in ev:
            rec["Cd"], rec["Cl"] = ev["Cd"], ev["Cl"]
            # wake-case tables read the *_mean keys: keep both consistent
            # (review M9: the reference side was updated, this side not)
            if rec.get("Cd_mean") is not None:
                rec["Cd_mean"], rec["Cl_mean"] = ev["Cd"], ev["Cl"]
                rec["Cd_std"], rec["Cl_std"] = ev["Cd_std"], ev["Cl_std"]
    if its is not None and rdir.is_dir():
        ev = evaluate(rdir, "simpleFoam", its, ref)
        if ev.get("ignored"):
            ref["userIgnored"] = ev["ignored"]
            its = None
    if its is not None and rdir.is_dir():
        ref["auto"] = {k: ref.get(k) for k in TEST_KEYS}
        ref["user"] = ev
        ref["itersToConv"] = its
        ref["wallToConv_s"] = ev["wall_to_conv_s"]
        ref["cpuHoursToConv"] = ev["cpu_to_conv_h"]
        if "Cd" in ev:
            ref["Cd"], ref["Cl"] = ev["Cd"], ev["Cl"]
            ref["Cd_mean"], ref["Cl_mean"] = ev["Cd"], ev["Cl"]
            ref["Cd_std"], ref["Cl_std"] = ev["Cd_std"], ev["Cl_std"]
            rec["CdRef"], rec["ClRef"] = ev["Cd"], ev["Cl"]
    if rec.get("Cd") is not None and rec.get("CdRef"):
        rec["CdRelDiff"] = abs(rec["Cd"] - rec["CdRef"]) / abs(rec["CdRef"])
        rec["ClRelDiff"] = abs(rec["Cl"] - rec["ClRef"]) / max(
            abs(rec["ClRef"]), 1e-12)
    if rec.get("wallToConv_s") and ref.get("wallToConv_s"):
        rec["speedupWall"] = ref["wallToConv_s"] / rec["wallToConv_s"]
    if rec.get("cpuHoursToConv") and ref.get("cpuHoursToConv"):
        rec["speedupCpu"] = ref["cpuHoursToConv"] / rec["cpuHoursToConv"]
    rec["reference"] = ref
    rec["convergenceSource"] = ("user" if itc is not None or its is not None
                                else "auto")
    return rec


def template(tests: dict, bench: list[dict]) -> dict:
    """Skeleton of the user file with the automatic iterations as hints."""
    out: dict = {"_doc": "Iteration from which a run counts as converged "
                 "(null = automatic criterion). See bench/user_convergence.py."}
    for name, rec in sorted(tests.items()):
        case = case_key(name)
        if case is None or not name.endswith(("_np1", "_np10")):
            continue
        e = out.setdefault(case, {})
        e["coupledFoam"] = None
        e["simpleFoam"] = None
        e[f"_auto coupledFoam ({name})"] = rec.get("autoIteration")
        e["_auto simpleFoam"] = (rec.get("reference") or {}).get("autoIteration")
    for r in bench:
        if r.get("case") in CASES:
            e = out.setdefault(r["case"], {})
            e.setdefault(r["config"], None)
            e.setdefault(f"_auto {r['config']}", r.get("iters_to_conv")
                         if r.get("convergenceSource") != "user"
                         else (r.get("auto") or {}).get("iters_to_conv"))
    return out

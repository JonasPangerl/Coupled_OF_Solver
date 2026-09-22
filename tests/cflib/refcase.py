"""Shared helpers for case tests: cached native references, coupled runs,
and the per-run record (iterations, wall and CPU time, memory)."""

from __future__ import annotations

import json
import re
from pathlib import Path

from . import case as cfcase
from . import env as cfenv
from . import logs

# ExecutionTime (CPU time of the solver process) and ClockTime (wall) of an
# OpenFOAM solver log; ClockTime is printed in whole seconds
_EXEC_CLOCK = re.compile(r"ExecutionTime = ([\d.eE+-]+) s +ClockTime = ([\d.eE+-]+) s")


def _last_exec_clock(log: Path) -> tuple[float, float] | None:
    last = None
    if not log.exists():
        return None
    with open(log, errors="replace") as fh:
        for line in fh:
            if "ExecutionTime" in line:
                m = _EXEC_CLOCK.search(line)
                if m:
                    last = (float(m.group(1)), float(m.group(2)))
    return last


def coupled_timing(summ: dict, timing: dict | None = None) -> dict:
    """wallSeconds, cpuSeconds, cpuHours of a coupledFoam run: the solver's
    own summary.json (solver loop, CPU summed over ranks); if it is missing
    (aborted run) the Allrun rusage (includes the pre-processing)."""
    wall, cpus, cpuh = (summ.get("wallSeconds"), summ.get("cpuSeconds"),
                        summ.get("cpuHours"))
    if cpuh is None and cpus is not None:
        cpuh = cpus / 3600.0
    if cpus is None and cpuh is not None:
        cpus = cpuh * 3600.0
    src = "coupledFoam summary.json"
    if wall is None or cpuh is None:
        t = timing or {}
        wall = wall if wall is not None else t.get("wallSeconds")
        cpus = t.get("cpuSeconds") if cpuh is None else cpus
        cpuh = t.get("cpuHours") if cpuh is None else cpuh
        src = "Allrun rusage (solver summary missing)"
    return {"wallSeconds": wall, "cpuSeconds": cpus, "cpuHours": cpuh,
            "timingSource": src}


def native_timing(log: Path, nprocs: int = 1,
                  timing: dict | None = None) -> dict:
    """wallSeconds, cpuSeconds, cpuHours of a native OpenFOAM solver run
    from its log: wall = final ClockTime, CPU = final ExecutionTime (CPU
    time of the master process) times the rank count. Without timing lines
    in the log: the Allrun rusage."""
    ec = _last_exec_clock(log)
    if ec is not None:
        cpus = ec[0] * max(1, int(nprocs))
        return {"wallSeconds": ec[1], "cpuSeconds": cpus,
                "cpuHours": cpus / 3600.0,
                "timingSource": "solver log (ClockTime, ExecutionTime x ranks)"}
    t = timing or {}
    return {"wallSeconds": t.get("wallSeconds"),
            "cpuSeconds": t.get("cpuSeconds"), "cpuHours": t.get("cpuHours"),
            "timingSource": "Allrun rusage (no timing in the solver log)"}


def reference(template: str, name: str, args: list[str],
              sets: dict | None = None) -> tuple[Path, dict]:
    """Native simpleFoam reference, serial, cached in run/<name>.

    Returns the case path and a record with iterations, wall and CPU time.
    The record is stored next to the case (reference.json) so a cached
    reference keeps its timing."""
    case = cfcase.RUN_ROOT / name
    meta = case / "reference.json"
    if cfcase.solver_ok(case, "simpleFoam") and meta.exists():
        rec = json.loads(meta.read_text())
        # records written before the timing fields: complete them from the
        # log (in memory; reference.json stays as it is)
        if rec.get("cpuHours") is None or rec.get("wallSeconds") is None:
            rec.update(native_timing(case / "log.simpleFoam", 1,
                                     rec.get("timingAllrun")))
        return case, rec

    case = cfcase.prepare(template, name, sets)
    rc = cfcase.allrun(case, ["-solver", "simpleFoam"] + list(args), fpe=False)
    timing = cfenv.last_timing.as_dict()
    native = logs.parse_native(case / "log.simpleFoam")
    rec = {
        "solver": "simpleFoam", "rc": rc,
        "iterations": native["iterations"],
        "convergedAt": native["convergedAt"],
        "wallSecondsSolver": native["wall"],
        "timingAllrun": timing,
        "args": args,
    }
    rec.update(native_timing(case / "log.simpleFoam", 1, timing))
    meta.write_text(json.dumps(rec, indent=2))
    assert rc == 0 and cfcase.solver_ok(case, "simpleFoam"), \
        f"simpleFoam reference {name} failed"
    return case, rec


def coupled(template: str, name: str, args: list[str],
            sets: dict | None = None, fpe: bool = True) -> tuple[Path, dict]:
    """coupledFoam run with the usual record."""
    case = cfcase.prepare(template, name, sets)
    rc = cfcase.allrun(case, ["-solver", "coupledFoam"] + list(args), fpe=fpe)
    timing = cfenv.last_timing.as_dict()
    log = case / "log.coupledFoam"
    rows = logs.parse_cf(log) if log.exists() else []
    summ = logs.coupled_summary(case)
    rec = {
        "solver": "coupledFoam", "rc": rc,
        "iterations": len(rows),
        "finalR": rows[-1]["R"] if rows else None,
        "fpeTrap": logs.fpe_trapped(log),
        "nClampedMax": max((r.get("nClamped", 0) for r in rows), default=None),
        "rollbacks": summ.get("rollbacks"),
        "wallSecondsSolver": summ.get("wallSeconds"),
        "cpuHoursSolver": summ.get("cpuHours"),
        "peakRSS_MB_sum": summ.get("peakRSS_MB_sum"),
        "peakRSS_MB_maxRank": summ.get("peakRSS_MB_maxRank"),
        "timingAllrun": timing,
        **coupled_timing(summ, timing),
        "history": {k: [r.get(k) for r in rows]
                    for k in ("R", "CFL", "omega", "cuts", "linIters",
                              "tIter", "tWall", "nDyn", "nLocLim", "eta",
                              "rho", "aa")},
        "gamgTuneEvents": summ.get("gamgTuneEvents"),
        "gamg": {k: summ.get(k) for k in ("gamgLevels", "gamgMergeLevels",
                                          "gamgCop", "gamgCellsPerLevel")},
        "args": args,
    }
    # ranks per level, coarsening ratios, coarsest solver and the final
    # cycle type / nPostSweeps after autoTune are only in the log (D-030)
    if log.exists():
        for k, v in logs.gamg_log_stats(log).items():
            if rec["gamg"].get(k) is None:
                rec["gamg"][k] = v
    return case, rec

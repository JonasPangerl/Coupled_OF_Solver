"""Shared helpers for case tests: cached native references, coupled runs,
and the per-run record (iterations, wall and CPU time, memory)."""

from __future__ import annotations

import json
from pathlib import Path

from . import case as cfcase
from . import env as cfenv
from . import logs


def reference(template: str, name: str, args: list[str],
              sets: dict | None = None) -> tuple[Path, dict]:
    """Native simpleFoam reference, serial, cached in run/<name>.

    Returns the case path and a record with iterations, wall and CPU time.
    The record is stored next to the case (reference.json) so a cached
    reference keeps its timing."""
    case = cfcase.RUN_ROOT / name
    meta = case / "reference.json"
    if cfcase.solver_ok(case, "simpleFoam") and meta.exists():
        return case, json.loads(meta.read_text())

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
        # C2: pseudo-inverse fallbacks of the tensorial Rhie-Chow D
        "nPseudoInverse": summ.get("nPseudoInverse", 0),
        "wallSecondsSolver": summ.get("wallSeconds"),
        "cpuHoursSolver": summ.get("cpuHours"),
        "peakRSS_MB_sum": summ.get("peakRSS_MB_sum"),
        "peakRSS_MB_maxRank": summ.get("peakRSS_MB_maxRank"),
        "timingAllrun": timing,
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

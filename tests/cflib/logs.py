"""Log parsers: coupledFoam CF| lines (spec 12.1) and native simpleFoam."""

from __future__ import annotations

import json
import re
from pathlib import Path


def _num(v: str):
    try:
        if re.fullmatch(r"[-+]?\d+", v):
            return int(v)
        return float(v)
    except ValueError:
        return v


def parse_cf(log: Path) -> list[dict]:
    """All CF| lines of a coupledFoam log as dicts of numbers."""
    rows = []
    with open(log, errors="replace") as fh:
        for line in fh:
            if not line.startswith("CF|"):
                continue
            row = {}
            for tok in line.split()[1:]:
                if "=" in tok:
                    k, v = tok.split("=", 1)
                    row[k] = _num(v)
            rows.append(row)
    return rows


_FPE_TRAP = re.compile(r"sigFpe::sigHandler|Floating point exception \(core dumped\)")


def fpe_trapped(log: Path) -> bool:
    """True if a floating-point exception was trapped. The start-up line
    "trapFpe: Floating point exception trapping enabled" is not a trap."""
    if not log.exists():
        return False
    return bool(_FPE_TRAP.search(log.read_text(errors="replace")))


def coupled_summary(case: Path) -> dict:
    f = case / "postProcessing" / "coupledFoam" / "summary.json"
    return json.loads(f.read_text()) if f.exists() else {}


_TIME = re.compile(r"^Time = (\S+)")
_RES = re.compile(
    r"^(\w+): +Solving for (\w+), Initial residual = ([^,]+), "
    r"Final residual = ([^,]+), No Iterations (\d+)"
)
_EXEC = re.compile(r"ExecutionTime = ([\d.eE+-]+) s +ClockTime = ([\d.eE+-]+) s")


def parse_native(log: Path) -> dict:
    """Iterations, per-field initial residual history (first solve of each
    field per iteration), timing and the SIMPLE convergence flag."""
    iters = []
    res: dict[str, list[float]] = {}
    clock = []
    converged_at = None
    seen: set[str] = set()
    with open(log, errors="replace") as fh:
        for line in fh:
            m = _TIME.match(line)
            if m:
                iters.append(_num(m.group(1)))
                seen = set()
                continue
            m = _RES.match(line)
            if m:
                field = m.group(2)
                if field not in seen:
                    res.setdefault(field, []).append(float(m.group(3)))
                    seen.add(field)
                continue
            m = _EXEC.search(line)
            if m:
                clock.append(float(m.group(2)))
                continue
            if "SIMPLE solution converged in" in line:
                converged_at = _num(line.split()[-2])
    return {
        "iterations": len(iters),
        "times": iters,
        "residuals": res,
        "clock": clock,
        "wall": clock[-1] if clock else None,
        "convergedAt": converged_at,
    }


def native_iterations_to(res: dict, tol: float, fields=("p", "Ux", "Uy")) -> int | None:
    """First iteration at which all listed initial residuals are < tol."""
    hist = [res["residuals"].get(f) for f in fields if f in res["residuals"]]
    if not hist:
        return None
    n = min(len(h) for h in hist)
    for i in range(n):
        if all(h[i] < tol for h in hist):
            return i + 1
    return None


def coupled_iterations_to(rows: list[dict], tol: float) -> int | None:
    for r in rows:
        if r.get("R", 1.0) < tol:
            return int(r["iter"])
    return None


_GAMG_HEAD = re.compile(r"^blockGAMG: levels (\d+), mergeLevels (\d+), "
                        r"C_op ([\d.eE+-]+), cycle (\w+)")
_GAMG_LIST = {
    "gamgCellsPerLevel": re.compile(r"^blockGAMG: cells per level \(([^)]*)\)"),
    "gamgRanksPerLevel": re.compile(r"^blockGAMG: ranks per level \(([^)]*)\)"),
    "gamgRatios": re.compile(r"^blockGAMG: coarsening ratios \(([^)]*)\)"),
}
_GAMG_COARSEST = re.compile(r"^blockGAMG: coarsest level \d+ \(.*\): (.+)$")
_NPOST = re.compile(r"^\s*nPostSweeps\s+(\d+);")
_TUNE = re.compile(r"^GAMG-tune: (cycleType|nPostSweeps) (\w+)->(\w+)")


def gamg_log_stats(log: Path) -> dict:
    """Block-GAMG hierarchy and autoTune outcome from a solver log.

    The first hierarchy report of the run is used (levels, mergeLevels,
    C_op, cycle, cells/ranks per level, coarsening ratios, coarsest-level
    solver). The initial nPostSweeps comes from the printed effective
    settings; the final cycle type and nPostSweeps replay the
    'GAMG-tune: <what> a->b' events of the autoTune controller. Keys that
    the log does not contain are absent (older builds print no ranks line).
    """
    out: dict = {}
    if not log.exists():
        return out
    ntune = 0
    with open(log, errors="replace") as fh:
        for line in fh:
            if "gamgLevels" not in out:
                m = _GAMG_HEAD.match(line)
                if m:
                    out.update({"gamgLevels": int(m.group(1)),
                                "gamgMergeLevels": int(m.group(2)),
                                "gamgCop": float(m.group(3)),
                                "gamgCycleInitial": m.group(4)})
                    continue
            for key, rx in _GAMG_LIST.items():
                if key not in out:
                    m = rx.match(line)
                    if m:
                        out[key] = [_num(v) for v in m.group(1).split()]
            if "gamgCoarsestSolver" not in out:
                m = _GAMG_COARSEST.match(line)
                if m:
                    out["gamgCoarsestSolver"] = m.group(1).strip()
            if "gamgNPostSweepsInitial" not in out:
                m = _NPOST.match(line)
                if m:
                    out["gamgNPostSweepsInitial"] = int(m.group(1))
            m = _TUNE.match(line)
            if m:
                ntune += 1
                if m.group(1) == "cycleType":
                    out["gamgCycleFinal"] = m.group(3)
                else:
                    out["gamgNPostSweepsFinal"] = int(m.group(3))
    if out:
        out["gamgTuneChanges"] = ntune
        out.setdefault("gamgCycleFinal", out.get("gamgCycleInitial"))
        out.setdefault("gamgNPostSweepsFinal", out.get("gamgNPostSweepsInitial"))
    return out

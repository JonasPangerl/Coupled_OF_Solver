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

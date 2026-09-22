"""JSON result files (spec rule 0.4): every measured number goes to
results/<kind>/<name>.json together with the provenance of the run.

Every record carries gitCommit (short hash, "-dirty" if tracked files
outside results/ and report/ were modified) - the key of the report's
staleness guard (TASK 6) - and, for solver runs, wallSeconds AND cpuHours
of coupledFoam and of the reference (user rule: every performance number
as wall-clock time and CPU-hours; completed from the solver-specific
fields by complete_timing() when a test did not set them)."""

from __future__ import annotations

import json
import math
import os
import time
from pathlib import Path

from .env import RESULTS, machine_state
from . import provenance


def git_commit() -> str:
    head = provenance.git_head()
    if not head:
        return "unknown"
    return head[:7] + ("-dirty" if provenance.git_dirty_files() else "")


def _clean(v):
    """JSON cannot hold NaN/Inf: map to None; numpy -> python."""
    try:
        import numpy as np
        if isinstance(v, np.generic):
            v = v.item()
        if isinstance(v, np.ndarray):
            v = v.tolist()
    except ImportError:
        pass
    if isinstance(v, float) and not math.isfinite(v):
        return None
    if isinstance(v, dict):
        return {k: _clean(x) for k, x in v.items()}
    if isinstance(v, (list, tuple)):
        return [_clean(x) for x in v]
    return v


def _scalar(v) -> bool:
    return v is None or isinstance(v, (int, float))


def _fill(d: dict) -> None:
    """wallSeconds / cpuSeconds / cpuHours of one run record from its
    solver-specific fields (…Solver, ranks, timingAllrun), in place."""
    if not isinstance(d, dict):
        return
    if not all(_scalar(d.get(k)) for k in ("wallSeconds", "cpuHours")):
        return      # per-rank-count dicts (T-scaling summary): not a run
    t = d.get("timingAllrun") if isinstance(d.get("timingAllrun"), dict) else {}
    r = d.get("ranks") if isinstance(d.get("ranks"), dict) else {}
    if d.get("wallSeconds") is None:
        for v in (d.get("wallSecondsSolver"), r.get("wallSeconds"),
                  t.get("wallSeconds")):
            if v is not None:
                d["wallSeconds"] = v
                break
    if d.get("cpuHours") is None:
        for v in (d.get("cpuHoursSolver"), r.get("cpuHours"),
                  t.get("cpuHours")):
            if v is not None:
                d["cpuHours"] = v
                break
    if d.get("cpuSeconds") is None and d.get("cpuHours") is not None:
        d["cpuSeconds"] = d["cpuHours"] * 3600.0


def complete_timing(record: dict) -> dict:
    """Fill wallSeconds and cpuHours of the coupledFoam run (top level) and
    of the reference (record["reference"]) where the test left them out."""
    if any(k in record for k in ("wallSeconds", "wallSecondsSolver",
                                 "timingAllrun", "cpuHoursSolver")):
        _fill(record)
    if isinstance(record.get("reference"), dict):
        _fill(record["reference"])
    # T-scaling: runs[solver][ranks] are run records
    runs = record.get("runs")
    if isinstance(runs, dict):
        for per_solver in runs.values():
            if isinstance(per_solver, dict):
                for r in per_solver.values():
                    _fill(r)
    return record


def _tag(kind: str, name: str) -> str:
    """Test records of an SP run are <name>_sp (cflib.precision, D11); the
    benchmark has its own SP configurations (F1, F2) and is not tagged."""
    from . import precision  # noqa: PLC0415
    return precision.tag(name) if kind == "tests" else name


def write(kind: str, name: str, data: dict) -> Path:
    name = _tag(kind, name)
    out = RESULTS / kind / f"{name}.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    head = provenance.git_head()
    record = {
        "name": name,
        "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "gitCommit": git_commit(),
        "gitCommitFull": head,
        "WM_OPTIONS": os.environ.get("WM_OPTIONS"),
        "machine": machine_state().as_dict(),
    }
    record.update(data)
    complete_timing(record)
    out.write_text(json.dumps(_clean(record), indent=2) + "\n")
    return out


def read(kind: str, name: str) -> dict | None:
    """Record of the current precision (tests: <name>_sp in SP)."""
    return read_exact(kind, _tag(kind, name))


def read_exact(kind: str, name: str) -> dict | None:
    """Record <name> as given (no precision tag)."""
    f = RESULTS / kind / f"{name}.json"
    return json.loads(f.read_text()) if f.exists() else None

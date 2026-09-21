"""JSON result files (spec rule 0.4): every measured number goes to
results/<kind>/<name>.json together with the provenance of the run."""

from __future__ import annotations

import json
import math
import os
import subprocess
import time
from pathlib import Path

from .env import REPO, RESULTS, machine_state


def git_commit() -> str:
    try:
        out = subprocess.run(["git", "rev-parse", "--short", "HEAD"],
                             cwd=REPO, capture_output=True, text=True,
                             check=True).stdout.strip()
        dirty = subprocess.run(["git", "status", "--porcelain", "--untracked-files=no"],
                               cwd=REPO, capture_output=True, text=True,
                               check=True).stdout.strip()
        return out + ("-dirty" if dirty else "")
    except (FileNotFoundError, subprocess.CalledProcessError):
        return "unknown"


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


def write(kind: str, name: str, data: dict) -> Path:
    out = RESULTS / kind / f"{name}.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    record = {
        "name": name,
        "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "gitCommit": git_commit(),
        "WM_OPTIONS": os.environ.get("WM_OPTIONS"),
        "machine": machine_state().as_dict(),
    }
    record.update(data)
    out.write_text(json.dumps(_clean(record), indent=2) + "\n")
    return out


def read(kind: str, name: str) -> dict | None:
    f = RESULTS / kind / f"{name}.json"
    return json.loads(f.read_text()) if f.exists() else None

"""Case preparation and execution.

A run is a copy of a template under cases/ into RUN_ROOT/<name>, modified
with foamDictionary, and executed through the template's Allrun. Templates
are never modified.
"""

from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

from .env import REPO, RUN_ROOT, foam_env, run

CASES = REPO / "cases"


def prepare(template: str, name: str, sets: dict | None = None,
            reuse: bool = False) -> Path:
    """Copy cases/<template> to RUN_ROOT/<name> and apply dictionary sets.

    sets: {"system/fvSolution": {"coupled.maxIter": "300", ...}, ...}
    reuse: keep an existing run directory (e.g. a cached reference run).
    """
    foam_env()
    src = CASES / template
    dst = RUN_ROOT / name
    if dst.exists() and reuse:
        return dst
    if dst.exists():
        shutil.rmtree(dst)
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copytree(src, dst, symlinks=True)
    for fname, entries in (sets or {}).items():
        for entry, value in entries.items():
            set_entry(dst, fname, entry, value)
    return dst


def set_entry(case: Path, fname: str, entry: str, value) -> None:
    """foamDictionary -entry <entry> -set <value> <file>"""
    subprocess.run(
        ["foamDictionary", "-entry", entry, "-set", str(value), fname],
        cwd=case, check=True, stdout=subprocess.DEVNULL,
    )


def get_entry(case: Path, fname: str, entry: str) -> str:
    out = subprocess.run(
        ["foamDictionary", "-entry", entry, "-value", fname],
        cwd=case, check=True, capture_output=True, text=True,
    ).stdout
    return out.strip()


def allrun(case: Path, args: list[str] | None = None, fpe: bool = True,
           timeout: float | None = None, extra_env: dict | None = None) -> int:
    """Run the case's Allrun. fpe: FOAM_SIGFPE/FOAM_SETNAN (spec 9.1)."""
    env = {
        "FOAM_SIGFPE": "true" if fpe else "false",
        "FOAM_SETNAN": "true" if fpe else "false",
    }
    if extra_env:
        env.update(extra_env)
    return run(["./Allrun"] + list(args or []), cwd=case,
               log=case / "log.Allrun", env=env, timeout=timeout)


def time_dirs(case: Path, processor: bool = False) -> list[float]:
    base = case / "processor0" if processor else case
    times = []
    if not base.is_dir():
        return times
    for d in base.iterdir():
        if d.is_dir():
            try:
                times.append(float(d.name))
            except ValueError:
                pass
    return sorted(times)


def latest_time(case: Path) -> str | None:
    """Latest time directory name (serial or decomposed)."""
    for proc in (False, True):
        t = time_dirs(case, processor=proc)
        if t and t[-1] > 0:
            name = f"{t[-1]:g}"
            base = case / "processor0" if proc else case
            # Keep the directory's own spelling
            for d in base.iterdir():
                try:
                    if float(d.name) == t[-1]:
                        return d.name
                except ValueError:
                    pass
            return name
    return None


def postprocess(case: Path, args: list[str], nprocs: int = 1,
                log_name: str = "log.postProcess") -> int:
    """postProcess (serial) or with -parallel on the decomposed case."""
    from .env import mpirun_prefix
    cmd = ["postProcess"] + args
    if nprocs > 1:
        cmd = mpirun_prefix(nprocs) + cmd + ["-parallel"]
    return run(cmd, cwd=case, log=case / log_name)


def solver_ok(case: Path, solver: str) -> bool:
    log = case / f"log.{solver}"
    if not log.exists():
        return False
    text = log.read_text(errors="replace")
    return "FOAM FATAL" not in text and "\nEnd" in text

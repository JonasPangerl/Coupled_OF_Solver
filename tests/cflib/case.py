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

from . import precision
from .env import REPO, RUN_ROOT, foam_env, run

CASES = REPO / "cases"


def run_dir(name: str) -> Path:
    """RUN_ROOT/<name> of the current precision (<name>_sp in SP, D5)."""
    return RUN_ROOT / precision.tag(name)


def prepare(template: str, name: str, sets: dict | None = None,
            reuse: bool = False) -> Path:
    """Copy cases/<template> to RUN_ROOT/<name> and apply dictionary sets.

    sets: {"system/fvSolution": {"coupled.maxIter": "300", ...}, ...}
    reuse: keep an existing run directory (e.g. a cached reference run).
    SP (CF_PRECISION=sp): the run directory is <name>_sp; the mesh is made
    in DP and shifted by the first allrun() (cflib.precision).
    """
    foam_env()
    precision.check_environment()
    src = CASES / template
    dst = run_dir(name)
    if dst.exists() and reuse:
        return dst
    if dst.exists():
        shutil.rmtree(dst)
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copytree(src, dst, symlinks=True)
    for fname, entries in (sets or {}).items():
        for entry, value in entries.items():
            set_entry(dst, fname, entry, value)
    # the explicit sets of this run (the SP harness keeps them, D-065)
    import json  # noqa: PLC0415
    (dst / "harnessSets.json").write_text(
        json.dumps(sets or {}, indent=2, default=str) + "\n")
    return dst


def set_entry(case: Path, fname: str, entry: str, value) -> None:
    """foamDictionary -entry <entry> -set <value> <file>

    fvSolution: with -disableFunctionEntries, otherwise foamDictionary
    rewrites the file with directives expanded and silently drops the
    `#sinclude "relaxation"` of the case templates (k/omega unrelaxed;
    found 2026-09-22 on the T4 runs). Other dictionaries keep the plain
    call: controlDict/fvSchemes use $-macros ($inletP, $turbulence) that
    must be expanded, and the flag would write them back quoted.

    Missing sub-dictionaries of a dotted entry are created (foamDictionary
    -set refuses a path whose parent does not exist): e.g.
    coupled.rhieChow.tensorial on a template without coupled.rhieChow sets
    coupled.rhieChow to "{ tensorial <value>; }". Existing sibling entries
    of an existing parent are never replaced.
    """
    base = ["foamDictionary"]
    if Path(fname).name == "fvSolution":
        base.append("-disableFunctionEntries")

    def _run(ent: str, val: str) -> int:
        return subprocess.run(base + ["-entry", ent, "-set", val, fname],
                              cwd=case, stdout=subprocess.DEVNULL,
                              stderr=subprocess.DEVNULL).returncode

    def _exists(ent: str) -> bool:
        return subprocess.run(base + ["-entry", ent, fname], cwd=case,
                              stdout=subprocess.DEVNULL,
                              stderr=subprocess.DEVNULL).returncode == 0

    if _run(entry, str(value)) == 0:
        return
    parts = entry.split(".")
    # deepest existing parent: parts[:k] (k = 0: top level)
    k = len(parts) - 1
    while k > 0 and not _exists(".".join(parts[:k])):
        k -= 1
    if k == len(parts) - 1:
        # the parent exists: the failure was not a missing sub-dictionary
        subprocess.run(base + ["-entry", entry, "-set", str(value), fname],
                       cwd=case, check=True, stdout=subprocess.DEVNULL)
        return
    nested = f"{parts[-1]} {value};"
    for p in reversed(parts[k + 1:-1]):
        nested = f"{p} {{ {nested} }}"
    subprocess.run(base + ["-entry", ".".join(parts[:k + 1]),
                           "-set", f"{{ {nested} }}", fname],
                   cwd=case, check=True, stdout=subprocess.DEVNULL)


def get_entry(case: Path, fname: str, entry: str) -> str:
    out = subprocess.run(
        ["foamDictionary", "-entry", entry, "-value", fname],
        cwd=case, check=True, capture_output=True, text=True,
    ).stdout
    return out.strip()


def allrun(case: Path, args: list[str] | None = None, fpe: bool = True,
           timeout: float | None = None, extra_env: dict | None = None) -> int:
    """Run the case's Allrun. fpe: FOAM_SIGFPE/FOAM_SETNAN (spec 9.1).

    Writes <case>/provenance.json first (git commit and dirty flag, the
    coupledFoam binary and libcoupledFoam.so actually used, host, date,
    CF_* environment; cflib.provenance), which the report's staleness
    guard reads (TASK 6).

    SP (CF_PRECISION=sp, amendment D5): the first call on a case generates
    the mesh in DP, shifts it to the origin, shifts the point settings and
    runs the checkMesh gate (cflib.precision.prepare_case; raises
    SPGeometryFail), then runs Allrun with -keep-mesh. A -mesh-only call is
    completed by that preparation (Allrun is not run in SP)."""
    from . import provenance  # noqa: PLC0415
    if precision.is_sp():
        new_args = precision.prepare_case(case, list(args or []))
        if new_args is None:
            return 0
        args = new_args
    provenance.write(case, list(args or []))
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


def log_tail(case: Path, solver: str | None = None, n: int = 25) -> str:
    """The last n lines of log.Allrun and of log.<solver> (for failure
    messages: a failed mpirun, e.g. an invalid --cpu-set, leaves rc 1 and
    nothing else)."""
    out = []
    names = ["log.Allrun"] + ([f"log.{solver}"] if solver else [])
    for name in names:
        f = Path(case) / name
        if not f.exists():
            out.append(f"--- {name}: missing")
            continue
        lines = f.read_text(errors="replace").splitlines()
        out.append(f"--- {name} (last {min(n, len(lines))} of {len(lines)} "
                   "lines)")
        out.extend(lines[-n:])
    return "\n".join(out)


def run_failure(case: Path, solver: str, rc: int | None,
                iterations: int | None = None,
                budget: int | None = None) -> list[str]:
    """Reasons why a solver run failed; empty if it ended normally.

    Checks: Allrun return code, the solver log exists and ends normally
    ("End", no FOAM FATAL), and - if both are given - the run reached its
    iteration budget (a harness run with the solver's own stop disabled
    must run all of it)."""
    reasons = []
    if rc != 0:
        reasons.append(f"Allrun rc {rc}")
    log = Path(case) / f"log.{solver}"
    if not log.exists():
        reasons.append(f"no log.{solver}")
    else:
        text = log.read_text(errors="replace")
        if "FOAM FATAL" in text:
            reasons.append(f"FOAM FATAL in log.{solver}")
        if "\nEnd" not in text:
            reasons.append(f"log.{solver} has no normal end")
    if iterations is not None and budget is not None and iterations < budget:
        reasons.append(f"{iterations} of {budget} iterations")
    return reasons

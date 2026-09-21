"""OpenFOAM environment, machine load and process execution.

Nothing here sources an OpenFOAM environment: the caller must run pytest or
the benchmark harness from a shell that has sourced exactly one
etc/bashrc (see README). This keeps the choice of build explicit.
"""

from __future__ import annotations

import os
import resource
import shlex
import subprocess
import time
from dataclasses import dataclass, field
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
RUN_ROOT = Path(os.environ.get("CF_RUN_ROOT", REPO / "run"))
RESULTS = REPO / "results"

# Processes of other jobs that make the machine "busy" (shared workstation)
BUSY_PATTERN = r"snappyHexMesh|simpleFoam|foamRun|potentialFoam|pvbatch|mpirun"


def foam_env() -> dict:
    """Return the OpenFOAM environment variables or raise with a clear hint."""
    need = ("WM_PROJECT_DIR", "WM_OPTIONS", "FOAM_USER_APPBIN")
    missing = [k for k in need if not os.environ.get(k)]
    if missing:
        raise RuntimeError(
            "No OpenFOAM environment: source exactly one OpenFOAM v2606 "
            "etc/bashrc before running (missing %s)" % ", ".join(missing)
        )
    return {k: os.environ[k] for k in os.environ if k.startswith(("WM_", "FOAM_"))}


def load_average() -> float:
    try:
        return os.getloadavg()[0]
    except OSError:
        return 0.0


def other_jobs() -> list[str]:
    """Command lines of other solver/meshing/render processes (not ours)."""
    try:
        out = subprocess.run(
            ["pgrep", "-a", "-f", BUSY_PATTERN],
            capture_output=True, text=True, check=False,
        ).stdout
    except FileNotFoundError:
        return []
    mine = str(RUN_ROOT)
    lines = []
    for line in out.splitlines():
        if "pgrep" in line or mine in line:
            continue
        lines.append(line)
    return lines


def ncores_physical() -> int:
    """Physical core count (SMT excluded), fallback to logical/2."""
    try:
        out = subprocess.run(["lscpu", "-p=CORE,SOCKET"], capture_output=True,
                             text=True, check=True).stdout
        cores = {l for l in out.splitlines() if not l.startswith("#")}
        return max(1, len(cores))
    except (FileNotFoundError, subprocess.CalledProcessError):
        return max(1, (os.cpu_count() or 2) // 2)


@dataclass
class MachineState:
    load: float
    jobs: list[str]
    physical_cores: int

    @property
    def busy(self) -> bool:
        return bool(self.jobs) or self.load >= 0.8 * self.physical_cores

    def as_dict(self) -> dict:
        return {"loadavg": self.load, "otherJobs": self.jobs,
                "physicalCores": self.physical_cores, "busy": self.busy}


def machine_state() -> MachineState:
    return MachineState(load_average(), other_jobs(), ncores_physical())


def mpirun_prefix(nprocs: int, busy: bool | None = None) -> list[str]:
    """mpirun with core binding (spec 3); unbound when the machine is shared.

    CF_MPI_BIND overrides (e.g. "none" or "core").
    """
    if busy is None:
        busy = machine_state().busy
    bind = os.environ.get("CF_MPI_BIND", "none" if busy else "core")
    cmd = ["mpirun", "-np", str(nprocs), "--bind-to", bind]
    if bind == "core":
        cmd += ["--map-by", "core"]
    return cmd


@dataclass
class Timing:
    """Wall-clock and CPU time of one command, summed over all processes it
    started (MPI ranks are children of mpirun and are included once they
    are reaped). CPU time = user + system."""
    wall: float = 0.0
    cpu: float = 0.0

    @property
    def cpu_hours(self) -> float:
        return self.cpu / 3600.0

    def as_dict(self) -> dict:
        return {"wallSeconds": self.wall, "cpuSeconds": self.cpu,
                "cpuHours": self.cpu_hours}


# Timing of the last run() call (per process; tests run sequentially)
last_timing = Timing()


def _children_cpu() -> float:
    r = resource.getrusage(resource.RUSAGE_CHILDREN)
    return r.ru_utime + r.ru_stime


def run(cmd, cwd: Path, log: Path | None = None, env: dict | None = None,
        timeout: float | None = None, nice: bool = True) -> int:
    """Run a command (list or string), optionally logging stdout+stderr.
    Wall and CPU time are stored in env.last_timing."""
    global last_timing
    if isinstance(cmd, str):
        cmd = shlex.split(cmd)
    if nice:
        cmd = ["nice", "-n", os.environ.get("CF_NICE", "19")] + list(cmd)
    full_env = dict(os.environ)
    if env:
        full_env.update(env)
    cpu0 = _children_cpu()
    t0 = time.perf_counter()
    if log is None:
        proc = subprocess.run(cmd, cwd=cwd, env=full_env, timeout=timeout,
                              check=False)
    else:
        with open(log, "w") as fh:
            proc = subprocess.run(cmd, cwd=cwd, env=full_env, stdout=fh,
                                  stderr=subprocess.STDOUT, timeout=timeout,
                                  check=False)
    last_timing = Timing(time.perf_counter() - t0, _children_cpu() - cpu0)
    return proc.returncode

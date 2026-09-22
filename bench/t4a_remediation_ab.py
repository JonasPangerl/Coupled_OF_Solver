"""T4a A/B of the remediation categories (D-061, D-066).

Runs T4a (coarse motorBike, 800 iterations, 10 ranks) with the wall
category variants and compares the D-042 window means with the simpleFoam
reference (Cd 0.3964, Cl 0.0768; tolerance max(2 %, 0.01) abs on Cl):

    default  template defaults (wall.wallStarved no)
    wsMild   wall.wallStarved yes with the mild wall treatment (D-061)
    wsFull   wall.wallStarved yes with the full treatment (beta 0, both
             limiters): the pre-D-061 behaviour (main: Cl 0.0630)

Adapted from the lead's t4a_cl_isolate.py. Uses the harness
(tests/test_T4_motorBike.py run_solver) of THIS worktree and the solver in
PATH; the cached mesh run/T4a_mesh may be a read-only link to another
worktree's cache. No CPU set (CF_MPI_BIND=none): timing is not the point.
Before each run it waits (bounded) while another 10-rank coupledFoam of a
different case directory is running.

Usage: python bench/t4a_remediation_ab.py [variant ...]
"""
from __future__ import annotations

import json
import os
import sys
import time
from pathlib import Path

R = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(R / "tests"), str(R / "bench")]
os.environ.pop("CF_MPI_CPUSET", None)
os.environ.setdefault("CF_MPI_BIND", "none")
import test_T4_motorBike as T4  # noqa: E402

REF = {"Cd": 0.39640, "Cl": 0.07675}
FULL_WALL = {"coupled.remediation.wall.wallStarved": "yes",
             "coupled.remediation.wall.beta": 0,
             "coupled.remediation.wall.cflFactor": 0.5,
             "coupled.remediation.wall.gradLimiter": "yes",
             "coupled.remediation.wall.nonOrthLimiter": "yes"}
VARIANTS = {
    "default": {},
    "wsMild": {"coupled.remediation.wall.wallStarved": "yes"},
    "wsFull": FULL_WALL,
}
OUT = R / "run" / "T4a_rem_ab.json"


def foreign_big_runs(min_ranks: int = 8) -> dict[str, int]:
    """Case directories (not ours) with >= min_ranks coupledFoam processes."""
    counts: dict[str, int] = {}
    mine = str(R / "run")
    for pid in os.listdir("/proc"):
        if not pid.isdigit():
            continue
        try:
            with open(f"/proc/{pid}/comm") as fh:
                if fh.read().strip() != "coupledFoam":
                    continue
            cwd = os.readlink(f"/proc/{pid}/cwd")
        except OSError:
            continue
        if cwd.startswith(mine):
            continue
        case = cwd.split("/processor")[0]
        counts[case] = counts.get(case, 0) + 1
    return {c: n for c, n in counts.items() if n >= min_ranks}


def wait_free(max_minutes: int = 240) -> bool:
    for i in range(max_minutes):
        busy = foreign_big_runs()
        if not busy:
            return True
        if i % 10 == 0:
            print(f"waiting: other 10-rank coupledFoam running {busy}",
                  flush=True)
        time.sleep(60)
    return False


def summary_categories(case: Path) -> dict:
    f = case / "postProcessing" / "coupledFoam" / "summary.json"
    try:
        s = json.loads(f.read_text())
    except (OSError, ValueError):
        return {}
    return {"staticCategories": s.get("staticCategories"),
            "staticCells": s.get("staticCells"),
            "dynamicCells": s.get("dynamicCells")}


def main() -> None:
    only = sys.argv[1:] or list(VARIANTS)
    budget = T4.BUDGET["a"]["coupledFoam"]
    mesh = T4.mesh_dir(T4.TEMPLATE, "T4a_mesh", ["-mesh", "a"], T4.NP)
    out = json.loads(OUT.read_text()) if OUT.exists() else {}
    for name in only:
        if not wait_free():
            print(f"{name}: skipped, machine still busy", flush=True)
            continue
        extra = VARIANTS[name]
        sets = T4.budget_sets("coupledFoam", budget,
                              {"system/fvSolution": extra} if extra else None)
        case, rec = T4.run_solver(T4.TEMPLATE, mesh, f"T4a_rem_{name}",
                                  "coupledFoam", ["-mesh", "a"], sets,
                                  oscillatory=True)
        st = {k: rec.get(k) for k in (
            "Cd_mean", "Cl_mean", "Cd_std", "Cl_std", "stationary",
            "iters_to_stationary", "staticCells", "staticFraction",
            "rollbacks", "wallSeconds", "cpuHours")}
        st.update(summary_categories(case))
        st["CdRelDiff"] = (abs(st["Cd_mean"] - REF["Cd"]) / REF["Cd"]
                           if st["Cd_mean"] is not None else None)
        st["ClAbsDiff"] = (abs(st["Cl_mean"] - REF["Cl"])
                           if st["Cl_mean"] is not None else None)
        st["rc"] = rec["rc"]
        st["sets"] = extra
        out[name] = st
        OUT.write_text(json.dumps(out, indent=2))
        print(f"{name:10s} Cd {st['Cd_mean']} Cl {st['Cl_mean']} "
              f"dCd {st['CdRelDiff']} dCl(abs) {st['ClAbsDiff']} "
              f"static {st['staticCells']} {st.get('staticCategories')} "
              f"rollbacks {st['rollbacks']} wall {st['wallSeconds']} "
              f"cpuh {st['cpuHours']} rc {st['rc']}", flush=True)
    print("written", OUT)


if __name__ == "__main__":
    main()

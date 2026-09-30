"""T-scaling - strong scaling on T4b (spec 13). Heavy.

CF_SCALING_MESH=a runs it on T4a instead, CF_SCALING_ITERS overrides the
iteration count (final re-run budget of 12 h, user 2026-09-22: T4a, 150
iterations, ranks 1,2,4,8,12,16 on physical cores).

Ranks 1, 2, 4, 8, 16 and NP = CF_HEAVY_NP (default 48: the user caps this
machine at 48 of its 64 cores, D-073; CF_SCALING_RANKS overrides),
300 iterations each,
coupledFoam and simpleFoam (SIMPLEC), same mesh and decomposition method
(scotch). Time per iteration is measured inside the solver loop from the
logs (coupledFoam: CF| tWall, or the sum of tIter; simpleFoam: ClockTime),
excluding the first SKIP iterations (start-up, upwind start of coupledFoam);
wall time, CPU-hours and peak RSS of every run come from the per-rank
/usr/bin/time -v reports.

Parallel efficiency E(N) = t(N0) * N0 / (t(N) * N) with N0 the smallest rank
count (1 in the spec run).
Pass: E_coupledFoam(N) >= 0.8 * E_simpleFoam(N), N the largest rank count.

The T4b mesh is built once with NP ranks into run/T4b_mesh and reused; the
other rank counts are derived from it (run/T4b_mesh_np<N>, see
test_T4_motorBike.mesh_dir).
"""

from __future__ import annotations

import os

import numpy as np
import pytest

from cflib import logs, results

from test_T4_motorBike import NP, TEMPLATE, mesh_dir, run_solver

N_ITER = int(os.environ.get("CF_SCALING_ITERS", "300"))
SKIP = 50
SOLVERS = ("coupledFoam", "simpleFoam")
MIN_RELATIVE_EFFICIENCY = 0.8


def ranks() -> list[int]:
    r = os.environ.get("CF_SCALING_RANKS")
    if r:
        return sorted({int(x) for x in r.split(",") if x.strip()})
    return sorted({n for n in (1, 2, 4, 8, 16) if n < NP} | {NP})


def fixed_iterations(solver: str) -> dict:
    """Exactly N_ITER iterations, every stop criterion disabled."""
    sets = {"system/controlDict": {"endTime": N_ITER,
                                   "writeInterval": N_ITER}}
    if solver == "coupledFoam":
        sets["system/fvSolution"] = {
            "coupled.maxIter": N_ITER,
            "coupled.convergence.residualTol": 0,
            "coupled.convergence.forceCoeffsWindow": 10 * N_ITER,
        }
    else:
        sets["system/fvSolution"] = {"SIMPLE.residualControl": "{}"}
    return sets


def time_per_iteration(case, solver: str) -> dict:
    """Seconds per iteration over iterations SKIP+1..end, from the log."""
    if solver == "coupledFoam":
        rows = logs.parse_cf(case / "log.coupledFoam")
        n = len(rows)
        if n <= SKIP:
            return {"iterations": n}
        if all("tWall" in r for r in rows):
            t = rows[-1]["tWall"] - rows[SKIP - 1]["tWall"]
            src = "tWall"
        else:
            t = sum(r.get("tIter", 0.0) for r in rows[SKIP:])
            src = "tIter"
        per = t / (n - SKIP)
        med = float(np.median([r.get("tIter", np.nan) for r in rows[SKIP:]]))
        return {"iterations": n, "timePerIter_s": per, "source": src,
                "tIterMedian_s": med}
    nat = logs.parse_native(case / f"log.{solver}")
    clk = nat["clock"]
    n = len(clk)
    if n <= SKIP:
        return {"iterations": n}
    per = (clk[-1] - clk[SKIP - 1]) / (n - SKIP)
    return {"iterations": n, "timePerIter_s": per, "source": "ClockTime"}


@pytest.mark.heavy
def test_scaling(foam):
    variant = os.environ.get("CF_SCALING_MESH", "b")
    mesh_args = ["-mesh", variant]
    rs = ranks()
    runs: dict[str, dict[int, dict]] = {s: {} for s in SOLVERS}
    for n in rs:
        mesh = mesh_dir(TEMPLATE, f"T4{variant}_mesh", mesh_args, n, NP)
        for solver in SOLVERS:
            name = f"scaling_T4{variant}_{solver}_np{n}"
            case, rec = run_solver(TEMPLATE, mesh, name, solver, mesh_args,
                                   fixed_iterations(solver), fpe=False)
            rec.pop("history", None)
            rec.update(time_per_iteration(case, solver))
            runs[solver][n] = rec
            assert rec["rc"] == 0 and not rec.get("failed"), \
                (f"{name} failed ({'; '.join(rec.get('failure') or [])}):\n"
                 f"{rec.get('logTail')}")

    n0 = rs[0]
    summary: dict = {"mesh": f"T4{variant}", "ranks": rs, "baseRanks": n0, "iterations": N_ITER,
                     "skip": SKIP, "meshCells": None}
    for solver in SOLVERS:
        t0 = runs[solver][n0].get("timePerIter_s")
        eff, speedup = {}, {}
        for n in rs:
            t = runs[solver][n].get("timePerIter_s")
            if t and t0:
                speedup[n] = t0 / t
                eff[n] = t0 * n0 / (t * n)
        summary[solver] = {
            "timePerIter_s": {n: runs[solver][n].get("timePerIter_s")
                              for n in rs},
            "speedup": speedup,
            "efficiency": eff,
            "wallSeconds": {n: runs[solver][n].get("wallSeconds") for n in rs},
            "cpuHours": {n: runs[solver][n].get("cpuHours") for n in rs},
            "peakRSS_GB_sum": {n: runs[solver][n].get("peakRSS_GB_sum")
                               for n in rs},
        }
        summary["meshCells"] = runs[solver][n0].get("meshCells")

    nmax = rs[-1]
    ec = summary["coupledFoam"]["efficiency"].get(nmax)
    es = summary["simpleFoam"]["efficiency"].get(nmax)
    rel = ec / es if ec is not None and es else None
    summary.update({
        "efficiencyRanks": nmax,
        "efficiencyCoupled": ec, "efficiencySimple": es,
        "relativeEfficiency": rel,
        "minRelativeEfficiency": MIN_RELATIVE_EFFICIENCY,
        "specRanks": rs == [1, 2, 4, 8, 16],
        "pass": rel is not None and rel >= MIN_RELATIVE_EFFICIENCY,
        "runs": {s: {str(n): r for n, r in runs[s].items()} for s in SOLVERS},
    })
    # JSON keys must be strings
    for solver in SOLVERS:
        for k in ("timePerIter_s", "speedup", "efficiency", "wallSeconds",
                  "cpuHours", "peakRSS_GB_sum"):
            summary[solver][k] = {str(n): v
                                  for n, v in summary[solver][k].items()}
    results.write("tests", f"T_scaling_T4{variant}", summary)

    assert ec is not None and es is not None, "efficiency not measured"
    assert rel >= MIN_RELATIVE_EFFICIENCY, (
        f"coupledFoam efficiency {ec:.3f} at {nmax} ranks is "
        f"{rel:.2f} of simpleFoam's {es:.3f} (< {MIN_RELATIVE_EFFICIENCY})")

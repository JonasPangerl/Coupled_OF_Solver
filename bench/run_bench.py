#!/usr/bin/env python3
"""Benchmark harness (spec 14): native simpleFoam vs coupledFoam.

For each case and configuration the solver runs with its own stop criteria
disabled (fixed iteration budget). The harness then finds the convergence
iteration with the identical criterion for all solvers (spec 12.3 ii): over
the last `window` iterations max - min <= tol*|mean| for every monitored
quantity (Cd and Cl, or the pressure drop for T1/T2).

Configurations (DECISIONS.md D-025):
    A  simpleFoam, the case's native settings (tutorial solver settings and
       relaxation factors)
    B  simpleFoam SIMPLEC, consistent yes, relaxation p 1.0 / U 0.9 / k,omega 0.9
    C  coupledFoam defaults
    D  coupledFoam with preconditioner blockDiagonal (isolates the AMG gain)

Every rank runs under bench/rank_wrapper.sh (/usr/bin/time -v), so for every
run the harness records wall time, CPU time summed over all ranks (user +
sys, reported as CPU-hours), per-rank peak RSS (sum and max). Wall and CPU
time "to convergence" are the totals scaled by the solver's own time
progression at the convergence iteration (simpleFoam ExecutionTime/ClockTime,
coupledFoam tWall).

Usage (from the repository root, one OpenFOAM environment sourced):
    bench/run_bench.py --cases T1,T2 --configs A,B,C,D --repeats 3 --np 1
    bench/run_bench.py --list
Results: results/bench/<case>_<config>_<run>.json, results/bench/summary.csv

Machine sharing: the harness refuses to start while other solver jobs run
(timings would be meaningless) unless --allow-busy is given; the machine
state before each run is stored in the JSON.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import statistics
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "tests"))

from cflib import case as cfcase          # noqa: E402
from cflib import env as cfenv            # noqa: E402
from cflib import logs, post, results     # noqa: E402

WRAPPER = REPO / "bench" / "rank_wrapper.sh"

# Case table: template, extra Allrun args, monitor, iteration budgets, ranks
CASES = {
    "T1": {"template": "T1_pitzDaily", "args": [], "monitor": "dp",
           "iters": {"simpleFoam": 4000, "coupledFoam": 1000}, "np": 1},
    "T2": {"template": "T2_backwardFacingStep2D", "args": [], "monitor": "dp",
           "iters": {"simpleFoam": 6000, "coupledFoam": 1500}, "np": 1},
    "T3-SST": {"template": "T3_airFoil2D", "args": ["-turbulence", "kOmegaSST"],
               "monitor": "forces",
               "iters": {"simpleFoam": 6000, "coupledFoam": 2000}, "np": 1},
    "T3-GEKO": {"template": "T3_airFoil2D", "args": ["-turbulence", "GEKO"],
                "monitor": "forces",
                "iters": {"simpleFoam": 6000, "coupledFoam": 2000}, "np": 1},
    "T4a": {"template": "T4_motorBike", "args": ["-mesh", "a"],
            "monitor": "forces",
            "iters": {"simpleFoam": 3000, "coupledFoam": 1500}, "np": 16},
    "T4b": {"template": "T4_motorBike", "args": ["-mesh", "b"],
            "monitor": "forces",
            "iters": {"simpleFoam": 4000, "coupledFoam": 2000}, "np": 16},
    "T5": {"template": "T5_ahmed", "args": [], "monitor": "forces",
           "iters": {"simpleFoam": 5000, "coupledFoam": 2500}, "np": 16},
}

WINDOW = 100
TOL = 0.002

RELAX_B = """// Benchmark configuration B: SIMPLEC, p 1.0 / U 0.9 / k,omega 0.9
relaxationFactors
{
    fields
    {
        p               1;
    }
    equations
    {
        U               0.9;
        ".*"            0.9;
    }
}
"""


def config_sets(cfg: str, spec: dict) -> tuple[str, dict, dict]:
    """(solver, foamDictionary sets, files to write) of a configuration."""
    solver = "simpleFoam" if cfg in ("A", "B") else "coupledFoam"
    n = spec["iters"][solver]
    sets = {"system/controlDict": {"endTime": n, "writeInterval": n}}
    files = {}
    if solver == "simpleFoam":
        sets["system/fvSolution"] = {"SIMPLE.residualControl": "{}"}
        if cfg == "B":
            sets["system/fvSolution"]["SIMPLE.consistent"] = "yes"
            files["system/relaxation.simpleFoam"] = RELAX_B
    else:
        fv = {
            "coupled.maxIter": n,
            # internal stop disabled: the harness decides (spec 14)
            "coupled.convergence.residualTol": 0,
            "coupled.convergence.forceCoeffsWindow": 10 * n,
        }
        if cfg == "D":
            fv["solvers.coupled.preconditioner"] = "blockDiagonal"
        sets["system/fvSolution"] = fv
    return solver, sets, files


def monitor_history(case: Path, kind: str) -> dict[str, list[float]]:
    if kind == "dp":
        pin = post.surface_value(case, "inletP")
        pout = post.surface_value(case, "outletP")
        key = [k for k in pin if k.startswith("areaAverage")][0]
        n = min(len(pin[key]), len(pout[key]))
        return {"dp": (pin[key][:n] - pout[key][:n]).tolist()}
    fc = post.force_coeffs(case)
    return {"Cd": fc["Cd"].tolist(), "Cl": fc["Cl"].tolist()}


def iters_to_conv(hist: dict[str, list[float]]) -> int | None:
    """First iteration at which every monitored quantity satisfies the
    window criterion (spec 12.3 ii) at the same time."""
    n = min(len(h) for h in hist.values())
    for i in range(WINDOW, n + 1):
        if all(_window_ok(h[:i]) for h in hist.values()):
            return i
    return None


def _window_ok(h: list[float]) -> bool:
    if len(h) < WINDOW:
        return False
    w = h[-WINDOW:]
    mean = sum(w) / WINDOW
    return max(w) - min(w) <= TOL * abs(mean)


def rank_times(case: Path, solver: str) -> dict:
    """Parse the /usr/bin/time -v reports of all ranks."""
    reps = sorted((case / "timing").glob(f"{solver}.rank*.time"))
    wall, cpu, rss = [], [], []
    for r in reps:
        t = r.read_text()
        u = float(re.search(r"User time \(seconds\): ([\d.]+)", t).group(1))
        s = float(re.search(r"System time \(seconds\): ([\d.]+)", t).group(1))
        m = int(re.search(r"Maximum resident set size \(kbytes\): (\d+)", t).group(1))
        e = re.search(r"Elapsed \(wall clock\) time .*: ([\d:.]+)", t).group(1)
        parts = [float(p) for p in e.split(":")]
        w = sum(p * 60 ** i for i, p in enumerate(reversed(parts)))
        wall.append(w)
        cpu.append(u + s)
        rss.append(m)
    if not reps:
        return {}
    return {
        "ranks": len(reps),
        "wallSeconds": max(wall),
        "cpuSeconds": sum(cpu),
        "cpuHours": sum(cpu) / 3600.0,
        "peakRSS_GB_sum": sum(rss) / 1024**2,
        "peakRSS_GB_max_rank": max(rss) / 1024**2,
    }


def progress_fraction(case: Path, solver: str, it: int) -> float | None:
    """Fraction of the solver's run time spent up to iteration `it`."""
    log = case / f"log.{solver}"
    if solver == "coupledFoam":
        rows = logs.parse_cf(log)
        if not rows or "tWall" not in rows[-1]:
            return None
        return rows[it - 1]["tWall"] / rows[-1]["tWall"]
    nat = logs.parse_native(log)
    clk = nat["clock"]
    if not clk or clk[-1] <= 0:
        return None
    return clk[min(it, len(clk)) - 1] / clk[-1]


def coupled_breakdown(case: Path, it: int) -> dict:
    rows = logs.parse_cf(case / "log.coupledFoam")[:it]
    return {
        "t_assembly": sum(r.get("tAsm", 0) for r in rows),
        "t_linsolve": sum(r.get("tSolve", 0) for r in rows),
        "t_turb": sum(r.get("tTurb", 0) for r in rows),
        "linIters_total": sum(r.get("linIters", 0) for r in rows),
    }


def run_one(name: str, cfg: str, run: int, nprocs: int, force: bool) -> dict:
    spec = CASES[name]
    tag = f"{name}_{cfg}_{run}"
    if not force and results.read("bench", tag):
        print(f"skip {tag} (exists)")
        return results.read("bench", tag)

    solver, sets, files = config_sets(cfg, spec)
    state = cfenv.machine_state().as_dict()
    case = cfcase.prepare(spec["template"], f"bench_{tag}", sets)
    for rel, text in files.items():
        (case / rel).write_text(text)

    rc = cfcase.allrun(
        case, ["-solver", solver, "-np", str(nprocs)] + spec["args"],
        fpe=False,
        extra_env={"CF_RANK_WRAPPER": str(WRAPPER),
                   "CF_TIMING_DIR": str(case / "timing")},
    )
    rec = {"case": name, "config": cfg, "run": run, "solver": solver,
           "nProcs": nprocs, "rc": rc, "machineBefore": state}
    rt = rank_times(case, solver)
    rec["total"] = rt
    try:
        hist = monitor_history(case, spec["monitor"])
    except (FileNotFoundError, IndexError) as err:
        rec["error"] = f"monitor: {err}"
        results.write("bench", tag, rec)
        return rec

    it = iters_to_conv(hist)
    rec["iters_to_conv"] = it
    rec["iterations_run"] = min(len(h) for h in hist.values())
    if it is not None and rt:
        frac = progress_fraction(case, solver, it)
        rec["progressFraction"] = frac
        if frac is not None:
            rec["wall_to_conv_s"] = rt["wallSeconds"] * frac
            rec["cpu_to_conv_h"] = rt["cpuHours"] * frac
            rec["time_per_iter_s"] = rec["wall_to_conv_s"] / it
        rec.update({k: v[it - 1] for k, v in hist.items()})
        if solver == "coupledFoam":
            rec.update(coupled_breakdown(case, it))
    rec["peakRSS_GB_sum"] = rt.get("peakRSS_GB_sum")
    rec["peakRSS_GB_max_rank"] = rt.get("peakRSS_GB_max_rank")
    results.write("bench", tag, rec)
    print(f"{tag}: iters {it}, wall {rec.get('wall_to_conv_s')}, "
          f"CPU-h {rec.get('cpu_to_conv_h')}")
    return rec


def summarise(names: list[str], cfgs: list[str]) -> Path:
    rows = []
    for f in sorted((results.RESULTS / "bench").glob("*_*_*.json")):
        d = json.loads(f.read_text())
        if "case" in d:
            rows.append(d)
    out = results.RESULTS / "bench" / "summary.csv"
    fields = ["case", "config", "n", "iters_median", "wall_median",
              "wall_min", "wall_max", "cpuh_median", "rss_sum_GB", "Cd", "Cl",
              "dp", "speedup_wall_B_over_C", "speedup_cpu_B_over_C"]
    groups: dict = {}
    for d in rows:
        groups.setdefault((d["case"], d["config"]), []).append(d)

    def med(key, g):
        v = [x[key] for x in g if x.get(key) is not None]
        return statistics.median(v) if v else None

    table = []
    for (c, cfg), g in sorted(groups.items()):
        walls = [x["wall_to_conv_s"] for x in g if x.get("wall_to_conv_s")]
        table.append({
            "case": c, "config": cfg, "n": len(g),
            "iters_median": med("iters_to_conv", g),
            "wall_median": med("wall_to_conv_s", g),
            "wall_min": min(walls) if walls else None,
            "wall_max": max(walls) if walls else None,
            "cpuh_median": med("cpu_to_conv_h", g),
            "rss_sum_GB": med("peakRSS_GB_sum", g),
            "Cd": med("Cd", g), "Cl": med("Cl", g), "dp": med("dp", g),
        })
    by = {(t["case"], t["config"]): t for t in table}
    for t in table:
        b, cc = by.get((t["case"], "B")), by.get((t["case"], "C"))
        if b and cc and b["wall_median"] and cc["wall_median"]:
            t["speedup_wall_B_over_C"] = b["wall_median"] / cc["wall_median"]
        if b and cc and b["cpuh_median"] and cc["cpuh_median"]:
            t["speedup_cpu_B_over_C"] = b["cpuh_median"] / cc["cpuh_median"]
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        for t in table:
            w.writerow({k: t.get(k) for k in fields})
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cases", default="T1,T2,T3-SST,T3-GEKO")
    ap.add_argument("--configs", default="A,B,C,D")
    ap.add_argument("--repeats", type=int, default=3)
    ap.add_argument("--np", type=int, default=0,
                    help="ranks (0: the case default)")
    ap.add_argument("--force", action="store_true", help="rerun existing")
    ap.add_argument("--allow-busy", action="store_true")
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--summary-only", action="store_true")
    a = ap.parse_args()

    if a.list:
        for k, v in CASES.items():
            print(f"{k:8s} {v['template']:26s} np={v['np']} monitor={v['monitor']}")
        return 0
    names = a.cases.split(",")
    cfgs = a.configs.split(",")
    if not a.summary_only:
        cfenv.foam_env()
        for name in names:
            nprocs = a.np or CASES[name]["np"]
            for run in range(1, a.repeats + 1):
                for cfg in cfgs:
                    st = cfenv.machine_state()
                    if st.busy and not a.allow_busy:
                        print("machine busy (other jobs: %d, load %.1f): "
                              "refusing to benchmark; --allow-busy overrides"
                              % (len(st.jobs), st.load))
                        return 3
                    run_one(name, cfg, run, nprocs, a.force)
    print("summary:", summarise(names, cfgs))
    return 0


if __name__ == "__main__":
    sys.exit(main())

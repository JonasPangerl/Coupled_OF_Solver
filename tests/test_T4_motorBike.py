"""T4 - motorBike, kOmegaSST, snappyHexMesh (spec 13). Heavy,
CF_HEAVY_NP ranks (default 10, D-031).

Two meshes: (a) tutorial refinement (~350 k cells), (b) surface level
(6 6), features 7, refinementBox 5 (target 1-2 M cells, D-031; the
actual count is recorded).

Pass (per mesh), averaged force criterion for the oscillating wake (D-042
and its addendum, user decisions; the old 12.3(ii) 0.2 % min/max window is
unsatisfiable for any steady solver here, the simpleFoam reference
included): both solvers have a stationary window mean of Cd and Cl at the
end of their run (bench/run_bench.py:stationary_mean - ONE window W per
case for both solvers, run_bench.CASES[case]["statWindow"] = 400 on T4a
and T4b (D-068; the per-run window max(300, n/2) is recorded as the
*_perRun sensitivity value), half-window means differ by
<= max(1 % |mean|, 0.005)); the window means agree with simpleFoam
(SIMPLEC) on the same mesh: Cd within max(2 %, 0.002 absolute), Cl within
max(2 %, 0.01 absolute); the mean fields over the same window agree:
volume RMS of |dUMean|/U_inf and of |dpMean|/p_ref at most 0.02 each
(proposed thresholds; applications/utilities/coupledFieldCompare writes
the delta fields and fieldCompare.json; a cached reference without mean
fields is continued once with averaging, continue_reference_with_average);
static remediation set <= 1 % of the cells; nRollbacks == 0. Time to convergence
uses the first stationary window (run_bench.iters_to_stationary). Peak RSS,
wall time and CPU-hours of both solvers are recorded (every rank runs under
bench/rank_wrapper.sh). The coupledFoam runtime stop (convergence dict) is
unchanged: runs stop on their budget or their own 12.3(ii) stop.

Meshing is expensive, so every mesh is built once into run/<mesh name>
(Allrun -mesh-only) and reused while it exists: solver runs are fresh copies
of the template (current dictionaries) with the mesh copied in
(case_from_mesh). Delete run/T4a_mesh or run/T4b_mesh to re-mesh. The
helpers in this module are shared with test_T5_ahmed.py and
test_scaling.py.
"""

from __future__ import annotations

import json
import os
import shutil
import sys
from pathlib import Path

import numpy as np
import pytest

from cflib import case as cfcase
from cflib import env as cfenv
from cflib import logs, post, results

sys.path.insert(0, str(cfenv.REPO / "bench"))
import run_bench  # noqa: E402  (harness criterion, rank timing)

TEMPLATE = "T4_motorBike"
# Ranks of the heavy cases (the user caps this machine at 10 cores, D-031)
NP = int(os.environ.get("CF_HEAVY_NP", "10"))
WINDOW = run_bench.WINDOW           # 100
TOL = run_bench.TOL                 # 0.002
TOL_CD = 0.01
MAX_STATIC_FRACTION = 0.015   # user, D-047 (was 0.01)
WRAPPER = cfenv.REPO / "bench" / "rank_wrapper.sh"

# Iteration budgets: those of the benchmark harness (bench/run_bench.py)
BUDGET = {v: run_bench.CASES[f"T4{v}"]["iters"] for v in ("a", "b")}

# Files that make up a finished mesh (besides the polyMesh directories)
MESH_FILES = ("mesh.info", "mesh_cells.txt", "log.checkMesh",
              "log.snappyHexMesh", "meshing.json")


# ---------------------------------------------------------------------------
# Mesh cache
# ---------------------------------------------------------------------------

def mesh_info(d: Path) -> dict:
    """mesh.info written by the Allrun: variant, np, reconstructed, cells."""
    f = d / "mesh.info"
    out: dict[str, str] = {}
    if f.exists():
        for line in f.read_text().splitlines():
            parts = line.split(None, 1)
            if len(parts) == 2:
                out[parts[0]] = parts[1].strip()
    return out


def mesh_cells(d: Path) -> int | None:
    try:
        return int(mesh_info(d).get("cells", ""))
    except ValueError:
        return None


def _variant(mesh_args: list[str]) -> str | None:
    if "-mesh" in mesh_args:
        return mesh_args[mesh_args.index("-mesh") + 1]
    return None


def _mesh_ready(d: Path, variant: str | None, nprocs: int) -> bool:
    info = mesh_info(d)
    if not info or info.get("np") != str(nprocs) or not info.get("cells"):
        return False
    if variant is not None and info.get("variant") != variant:
        return False
    if nprocs > 1:
        return all((d / f"processor{i}" / "constant" / "polyMesh" / "owner")
                   .exists() for i in range(nprocs))
    return (d / "constant" / "polyMesh" / "owner").exists()


def copy_mesh(src: Path, dst: Path) -> None:
    """Copy the mesh of a meshed case (serial and decomposed polyMesh, mesh
    bookkeeping) into another case directory."""
    for f in MESH_FILES:
        if (src / f).exists():
            shutil.copy2(src / f, dst / f)
    subs = [Path("constant") / "polyMesh"]
    subs += [Path(p.name) / "constant" / "polyMesh"
             for p in sorted(src.glob("processor*"))]
    for sub in subs:
        s = src / sub
        if not s.is_dir():
            continue
        d = dst / sub
        if d.exists():
            shutil.rmtree(d)
        d.parent.mkdir(parents=True, exist_ok=True)
        shutil.copytree(s, d)


def mesh_dir(template: str, mesh_name: str, mesh_args: list[str],
             nprocs: int = NP, base_np: int = NP,
             extra_env: dict | None = None) -> Path:
    """Meshed case run/<mesh_name> for `nprocs` ranks, built once and reused.

    The mesh is generated with `base_np` ranks into run/<mesh_name>
    (Allrun <mesh_args> -np base_np -mesh-only). Other rank counts get
    run/<mesh_name>_np<N>: the serial one is reconstructed from the base
    (reconstructParMesh -constant), every other one is decomposed from the
    serial one; all by the Allrun (-mesh-only on a copied mesh)."""
    variant = _variant(mesh_args)
    name = mesh_name if nprocs == base_np else f"{mesh_name}_np{nprocs}"
    d = cfcase.RUN_ROOT / name
    if _mesh_ready(d, variant, nprocs):
        return d

    if nprocs == base_np:
        d = cfcase.prepare(template, name)
        args = list(mesh_args) + ["-np", str(nprocs), "-mesh-only", "-remesh"]
    else:
        src_np = base_np if nprocs == 1 else 1
        src = mesh_dir(template, mesh_name, mesh_args, src_np, base_np,
                       extra_env)
        d = cfcase.prepare(template, name)
        copy_mesh(src, d)
        (d / "meshing.json").unlink(missing_ok=True)
        args = list(mesh_args) + ["-np", str(nprocs), "-mesh-only"]

    rc = cfcase.allrun(d, args, fpe=False, extra_env=extra_env)
    timing = cfenv.last_timing.as_dict()
    (d / "meshing.json").write_text(json.dumps(
        {"args": args, "rc": rc, "timing": timing, "mesh": mesh_info(d)},
        indent=2) + "\n")
    assert rc == 0 and _mesh_ready(d, variant, nprocs), \
        f"meshing {name} failed, see {d}/log.Allrun"
    return d


def case_from_mesh(template: str, mesh: Path, name: str,
                   sets: dict | None = None) -> Path:
    """Fresh copy of the template (the current dictionaries, never those of
    the mesh directory) with the mesh of `mesh` copied in."""
    case = cfcase.prepare(template, name, sets)
    copy_mesh(mesh, case)
    np_mesh = mesh_info(mesh).get("np", "1")
    cfcase.set_entry(case, "system/decomposeParDict", "numberOfSubdomains",
                     np_mesh)
    return case


# ---------------------------------------------------------------------------
# Solver runs
# ---------------------------------------------------------------------------

def budget_sets(solver: str, n: int, extra: dict | None = None) -> dict:
    """Iteration budget. coupledFoam stops on its own criterion 12.3(ii)
    (window 100, tol 0.002 on Cd and Cl) - the residual criterion (i) is
    disabled so the stop is by (ii); simpleFoam keeps residualControl 1e-8
    (spec 13) and normally runs the full budget."""
    sets: dict = {"system/controlDict": {"endTime": n, "writeInterval": n}}
    if solver == "coupledFoam":
        sets["system/fvSolution"] = {
            "coupled.maxIter": n,
            "coupled.convergence.residualTol": 0,
            "coupled.convergence.forceCoeffsWindow": WINDOW,
            "coupled.convergence.forceCoeffsTol": TOL,
        }
    for f, entries in (extra or {}).items():
        sets.setdefault(f, {}).update(entries)
    return sets


def force_history(case: Path, t_min: float | None = None,
                  t_max: float | None = None) -> dict[str, list[float]]:
    """Cd / Cl history (restarts merged); optionally only the samples with
    t_min < Time <= t_max (a reference continuation is evaluated over its
    own iterations, the original reference over its original ones)."""
    try:
        fc = post.force_coeffs(case)
    except (FileNotFoundError, KeyError):
        return {}
    keep = np.ones(len(fc["Cd"]), dtype=bool)
    if (t_min is not None or t_max is not None) and "Time" in fc:
        t = fc["Time"]
        if t_min is not None:
            keep &= t > t_min
        if t_max is not None:
            keep &= t <= t_max
    return {"Cd": fc["Cd"][keep].tolist(), "Cl": fc["Cl"][keep].tolist()}


def evaluate_history(rec: dict, case: Path, solver: str,
                     hist: dict[str, list[float]], ranks: dict,
                     oscillatory: bool = False,
                     case_name: str | None = None) -> int | None:
    """Convergence, coefficients and time to convergence of a force history
    into `rec`; returns the convergence iteration.

    oscillatory (D-042): convergence := stationary final window
    (run_bench.stationary_eval with the common window of `case_name`,
    D-068); the convergence iteration is iters_to_stationary, Cd / Cl are
    the final-window means (with std, W, half-window drift). The per-run
    window is evaluated for information (iters_to_stationary_perRun,
    wallToConv_s_perRun, cpuHoursToConv_perRun). Otherwise criterion
    12.3(ii) (run_bench.iters_to_conv, window means at convergence)."""
    rec["iterationsRun"] = min((len(h) for h in hist.values()), default=0)
    if hist:
        rec["CdFinal"] = hist["Cd"][-1]
        rec["ClFinal"] = hist["Cl"][-1]
    if oscillatory:
        st = run_bench.stationary_eval(hist, case=case_name) \
            if rec["iterationsRun"] else \
            {"criterion": run_bench.STAT_CRITERION, "stationary": False,
             "iters_to_stationary": None}
        rec.update(st)
        it = st["iters_to_stationary"]
        rec["itersToConv"] = it
        # the 12.3(ii) result for information only (not a pass condition)
        rec["itersToConvWindow"] = (run_bench.iters_to_conv(hist)
                                    if rec["iterationsRun"] else None)
        if st.get("Cd_mean") is not None:
            rec["Cd"] = st["Cd_mean"]
            rec["Cl"] = st["Cl_mean"]
    else:
        rec["criterion"] = f"window {WINDOW}, tol {TOL} (12.3 ii)"
        it = run_bench.iters_to_conv(hist) if rec["iterationsRun"] else None
        rec["itersToConv"] = it
        if it is not None:
            rec["Cd"] = float(np.mean(hist["Cd"][it - WINDOW:it]))
            rec["Cl"] = float(np.mean(hist["Cl"][it - WINDOW:it]))
            # Cd alone with the post.py implementation (cross-check)
            rec["itersToConvCdOnly"] = post.window_converged(
                np.asarray(hist["Cd"]), WINDOW, TOL)
    if it is not None:
        frac = run_bench.progress_fraction(case, solver, it)
        rec["progressFraction"] = frac
        if frac is not None and ranks:
            # potentialFoam is a fixed offset, the solver part is scaled
            tc = run_bench.to_convergence(ranks, frac)
            rec["wallToConv_s"] = tc.get("wall_to_conv_s")
            rec["cpuHoursToConv"] = tc.get("cpu_to_conv_h")
    it_run = rec.get("iters_to_stationary_perRun") if oscillatory else None
    if it_run is not None and ranks:
        # D-068 sensitivity: time to the per-run window's stationary point
        tc = run_bench.to_convergence(
            ranks, run_bench.progress_fraction(case, solver, it_run))
        rec["wallToConv_s_perRun"] = tc.get("wall_to_conv_s")
        rec["cpuHoursToConv_perRun"] = tc.get("cpu_to_conv_h")
    return it


def run_solver(template: str, mesh: Path, name: str, solver: str,
               mesh_args: list[str], sets: dict | None = None,
               fpe: bool | None = None,
               extra_env: dict | None = None,
               oscillatory: bool = False,
               case_name: str | None = None) -> tuple[Path, dict]:
    """Run one solver on a copy of a cached mesh; returns the case and a
    record with convergence (harness criterion; D-042 stationary mean if
    `oscillatory`), coefficients, wall time, CPU-hours and peak RSS
    (per-rank /usr/bin/time -v).

    oscillatory: the fieldAverage function object of the template averages
    UMean/pMean over the force window, timeStart = n - W + 1 with n the
    endTime of `sets` and W the common window of `case_name` (D-068;
    run_bench.field_average_start, D-042 addendum). case_name defaults to
    the run name (run_bench.case_of_run)."""
    info = mesh_info(mesh)
    nprocs = int(info.get("np", "1"))
    case_name = case_name or run_bench.case_of_run(name)
    case = case_from_mesh(template, mesh, name, sets)
    avg_start = None
    n_budget = (sets or {}).get("system/controlDict", {}).get("endTime")
    if oscillatory and n_budget is not None:
        avg_start = run_bench.field_average_start(int(n_budget), case_name)
        run_bench.set_field_average_start(case, avg_start)
    env = {"CF_RANK_WRAPPER": str(WRAPPER),
           "CF_TIMING_DIR": str(case / "timing")}
    env.update(extra_env or {})
    if fpe is None:
        fpe = solver == "coupledFoam"
    args = ["-solver", solver, "-np", str(nprocs)] + list(mesh_args)
    rc = cfcase.allrun(case, args, fpe=fpe, extra_env=env)

    ranks = run_bench.rank_times(case, solver)
    rec: dict = {
        "solver": solver, "rc": rc, "nProcs": nprocs, "args": args,
        "meshCells": mesh_cells(mesh), "meshVariant": info.get("variant"),
        "timingAllrun": cfenv.last_timing.as_dict(),
        "ranks": ranks,
        "wallSeconds": ranks.get("wallSeconds"),
        "cpuHours": ranks.get("cpuHours"),
        "peakRSS_GB_sum": ranks.get("peakRSS_GB_sum"),
        "peakRSS_GB_max_rank": ranks.get("peakRSS_GB_max_rank"),
        "fieldAverageTimeStart": avg_start,
    }

    hist = force_history(case)
    it = evaluate_history(rec, case, solver, hist, ranks, oscillatory,
                          case_name)

    if solver == "coupledFoam":
        log = case / "log.coupledFoam"
        rows = logs.parse_cf(log) if log.exists() else []
        summ = logs.coupled_summary(case)
        ncells = summ.get("nCells") or rec["meshCells"]
        stat = summ.get("staticCells")
        rec.update({
            "fpeTrap": logs.fpe_trapped(log),
            "iterations": len(rows),
            "converged": summ.get("converged"),
            "rollbacks": summ.get("rollbacks"),
            "staticCells": stat,
            "staticFraction": (stat / ncells) if stat is not None and ncells
            else None,
            "peakRSS_MB_sum": summ.get("peakRSS_MB_sum"),
            "peakRSS_MB_maxRank": summ.get("peakRSS_MB_maxRank"),
            "wallSecondsSolver": summ.get("wallSeconds"),
            "cpuHoursSolver": summ.get("cpuHours"),
            "gamgTuneFailed": summ.get("gamgTuneFailed"),
            "history": {k: [r.get(k) for r in rows]
                        for k in ("R", "CFL", "omega", "cuts", "linIters",
                                  "tIter", "nStat", "nDyn", "nRollback")},
        })
        if it is not None:
            rec.update(run_bench.coupled_breakdown(case, it))
    else:
        log = case / f"log.{solver}"
        nat = logs.parse_native(log) if log.exists() else {}
        rec.update({
            "iterations": nat.get("iterations"),
            "convergedAt": nat.get("convergedAt"),
            "wallSecondsSolver": nat.get("wall"),
            "ok": cfcase.solver_ok(case, solver),
        })
    return case, rec


def reevaluate_reference(case: Path, rec: dict,
                         oscillatory: bool = False,
                         case_name: str | None = None) -> dict:
    """Re-evaluate a cached reference record with the current criterion
    (read-only: reference.json is not rewritten; the common window of
    `case_name`, D-068, default the case named by the directory). The force
    history and the rank timing reports are re-read from the case; records
    written before the pre-processing split lack solverWallSeconds, so the
    rank timing is re-parsed when needed."""
    case_name = case_name or run_bench.case_of_run(Path(case).name)
    rec = dict(rec)
    ranks = rec.get("ranks") or {}
    if ranks and "solverWallSeconds" not in ranks:
        ranks = run_bench.rank_times(case, "simpleFoam") or ranks
        rec["ranks"] = ranks
    for k in ("wallToConv_s", "cpuHoursToConv", "progressFraction", "Cd",
              "Cl", "itersToConvCdOnly", "wallToConv_s_perRun",
              "cpuHoursToConv_perRun"):
        rec.pop(k, None)
    # a continuation (continue_reference_with_average) extends the force
    # history; the reference itself is its original run only
    cont = rec.get("continuation") or {}
    hist = force_history(case, t_max=cont.get("startTime"))
    evaluate_history(rec, case, "simpleFoam", hist, ranks, oscillatory,
                     case_name)
    return rec


def reference(template: str, mesh: Path, name: str, mesh_args: list[str],
              sets: dict | None = None,
              oscillatory: bool = False,
              case_name: str | None = None) -> tuple[Path, dict]:
    """simpleFoam (SIMPLEC) reference on the same mesh, cached in
    run/<name> together with its record (reference.json). A cached record
    is re-evaluated with the current criterion (D-042 for `oscillatory`,
    common window of `case_name`, D-068)."""
    case = cfcase.RUN_ROOT / name
    meta = case / "reference.json"
    if cfcase.solver_ok(case, "simpleFoam") and meta.exists():
        rec = json.loads(meta.read_text())
        if rec.get("meshCells") == mesh_cells(mesh):
            return case, reevaluate_reference(case, rec, oscillatory,
                                              case_name)
    case, rec = run_solver(template, mesh, name, "simpleFoam", mesh_args,
                           sets, fpe=False, oscillatory=oscillatory,
                           case_name=case_name)
    meta.write_text(json.dumps(results._clean(rec), indent=2) + "\n")
    assert rec["rc"] == 0 and cfcase.solver_ok(case, "simpleFoam"), \
        f"simpleFoam reference {name} failed"
    return case, rec


# ---------------------------------------------------------------------------
# Mean fields of the reference (D-042 addendum)
# ---------------------------------------------------------------------------

CONTINUATION_TIMING = "timing_continuation"
CONTINUATION_CASES = ("T4a", "T4b", "T5")


def default_n_extra(case_name: str) -> int:
    """Continuation length of a cached reference: the per-run window of its
    budget, max(300, n/2) (T4a 1500, T4b 2000, T5 1000). Deliberately NOT
    the common force window of D-068 (400): the existing continuations of
    ref_T4a / ref_T4b have these lengths and are reused (a different
    length would be refused as chaining), and a longer reference average
    is the better estimate of the reference mean field."""
    return run_bench.per_run_window(
        run_bench.CASES[case_name]["iters"]["simpleFoam"])


def _case_name_of(ref_dir: Path) -> str | None:
    """run/ref_T4a_np10 -> T4a, run/ref_T5_coarse_np10 -> T5."""
    for name in CONTINUATION_CASES:
        if ref_dir.name.startswith(f"ref_{name}_"):
            return name
    return None


def continue_reference_with_average(case: Path, n_extra: int | None = None
                                    ) -> dict:
    """Continue a cached simpleFoam reference from its final time for
    n_extra iterations with the fieldAverage function object active from
    the first continued iteration, and evaluate the forces over the
    continuation only (D-042 addendum; the T4a/T4b references were run
    before the averaging existed and are not recomputed from scratch).

    n_extra default: default_n_extra of the case named by the directory
    (ref_T4a_* 1500, ref_T4b_* 2000, ref_T5* 2500). The run uses the
    Allrun -restart path (log.<solver>.restart) with endTime = t0 + n_extra,
    writeInterval n_extra, fieldAverage timeStart t0 + 1 and
    restartOnRestart true. Its rank timing goes to <case>/timing_continuation
    (the reference's own timing is untouched) and is recorded separately
    (continuationWallSeconds, continuationCpuHours). The record is stored in
    reference.json under "continuation" and returned. A completed
    continuation of the same length is reused; any other earlier
    continuation is an error (no silent chaining)."""
    case = Path(case)
    meta = case / "reference.json"
    rec = json.loads(meta.read_text()) if meta.exists() else {}
    if n_extra is None:
        name = _case_name_of(case)
        if name is None:
            raise ValueError(f"{case}: n_extra needed (unknown case)")
        n_extra = default_n_extra(name)
    n_extra = int(n_extra)
    t0_name = cfcase.latest_time(case)
    if t0_name is None:
        raise FileNotFoundError(f"{case}: no time directory to continue from")
    t0 = float(t0_name)
    old = rec.get("continuation")
    if old:
        done = (old.get("ok") and old.get("endTimeName")
                in run_bench.mean_field_times(case))
        if done and old.get("nExtra") == n_extra:
            return old
        if done or float(old.get("startTime", -1)) != t0:
            raise RuntimeError(
                f"{case}: earlier continuation {old.get('startTime')} -> "
                f"{old.get('endTime')} (ok {old.get('ok')}) exists; "
                "refusing to chain another one")

    solver = rec.get("solver", "simpleFoam")
    nprocs = int(rec.get("nProcs") or 1)
    args = list(rec.get("args") or ["-solver", solver, "-np", str(nprocs)])
    t_end = t0 + n_extra
    end_name = f"{t_end:g}"
    for entry, value in (("startFrom", "latestTime"), ("endTime", end_name),
                         ("writeInterval", str(n_extra))):
        run_bench.foam_dictionary(case, "system/controlDict",
                                  ["-entry", entry, "-set", value])
    run_bench.ensure_field_average(case, t0 + 1, restart_on_restart=True)
    # output of an interrupted earlier attempt from the same start time
    for d in (case / "postProcessing").glob(f"*/{t0_name}"):
        shutil.rmtree(d)
    tdir = case / CONTINUATION_TIMING
    if tdir.exists():
        shutil.rmtree(tdir)
    rc = cfcase.allrun(case, args + ["-restart"], fpe=False,
                       extra_env={"CF_RANK_WRAPPER": str(WRAPPER),
                                  "CF_TIMING_DIR": str(tdir)})
    timing_allrun = cfenv.last_timing.as_dict()
    log = case / f"log.{solver}.restart"
    text = log.read_text(errors="replace") if log.exists() else ""
    ok = rc == 0 and "FOAM FATAL" not in text and "\nEnd" in text
    rt = run_bench.rank_times(case, solver, CONTINUATION_TIMING)

    hist = force_history(case, t_min=t0, t_max=t_end)
    n = min((len(h) for h in hist.values()), default=0)
    cont: dict = {
        "method": "Allrun -restart, fieldAverage from the first continued "
                  "iteration (D-042 addendum)",
        "startTime": t0, "startTimeName": t0_name, "endTime": t_end,
        "endTimeName": end_name, "nExtra": n_extra,
        "fieldAverageTimeStart": t0 + 1, "rc": rc, "ok": ok,
        "log": str(log), "timingDir": str(tdir),
        "timingAllrun": timing_allrun, "ranks": rt,
        "continuationWallSeconds": rt.get("wallSeconds"),
        "continuationCpuHours": rt.get("cpuHours"),
        "iterationsRun": n,
        "meanFieldTimes": run_bench.mean_field_times(case),
    }
    # forces over the continuation window only (the whole continuation)
    for q in run_bench.STAT_QUANTITIES:
        if n < 2 or q not in hist:
            continue
        s = run_bench.stationary_mean(hist, q, n, w=n)
        cont.update({f"{q}_mean": s["mean"], f"{q}_std": s["std"],
                     f"{q}_drift": s["drift"], f"{q}_driftTol": s["driftTol"],
                     f"{q}_stationary": s["stationary"]})
    cont["stationary"] = bool(n >= 2 and all(
        cont.get(f"{q}_stationary") for q in run_bench.STAT_QUANTITIES))
    rec["continuation"] = cont
    rec["continuationWallSeconds"] = cont["continuationWallSeconds"]
    rec["continuationCpuHours"] = cont["continuationCpuHours"]
    meta.write_text(json.dumps(results._clean(rec), indent=2) + "\n")
    return cont


def reference_mean_time(ref_case: Path, ref: dict, name: str) -> str | None:
    """Time of the reference's mean fields. A reference run with the
    averaging has them at its final time; a cached reference without them
    is continued once (continue_reference_with_average, n_extra =
    default_n_extra(name)), which updates ref["continuation"]."""
    times = run_bench.mean_field_times(ref_case)
    cont = ref.get("continuation")
    if cont and cont.get("ok") and cont.get("endTimeName") in times:
        return cont["endTimeName"]
    latest = cfcase.latest_time(ref_case)
    if latest in times and not cont:
        return latest
    cont = continue_reference_with_average(ref_case, default_n_extra(name))
    ref["continuation"] = cont
    ref["continuationWallSeconds"] = cont.get("continuationWallSeconds")
    ref["continuationCpuHours"] = cont.get("continuationCpuHours")
    return cont["endTimeName"] if cont.get("ok") else None


def mean_field_comparison(case: Path, rec: dict, ref_case: Path, ref: dict,
                          name: str) -> dict:
    """coupledFieldCompare of the coupledFoam mean fields (final time)
    against the reference mean fields (D-042 addendum); both cases share
    the cached mesh and its decomposition."""
    own = run_bench.mean_field_times(case)
    ref_time = reference_mean_time(ref_case, ref, name)
    if not own or ref_time is None:
        out: dict = {"rc": None,
                     "error": "no UMean/pMean in "
                     + ("the coupledFoam run" if not own else "the reference")}
        out.update(run_bench.field_checks(out))
        return out
    out = run_bench.field_delta_compare(case, ref_case,
                                        int(rec.get("nProcs") or 1),
                                        time=own[-1], ref_time=ref_time)
    out["ownAverageTimeStart"] = rec.get("fieldAverageTimeStart")
    cont = ref.get("continuation") or {}
    out["referenceAverageTimeStart"] = (cont.get("fieldAverageTimeStart")
                                        if cont.get("endTimeName") == ref_time
                                        else ref.get("fieldAverageTimeStart"))
    return out


STAT_KEYS = ("W", "iterations_run", "Cd_mean", "Cd_std", "Cl_mean", "Cl_std",
             "Cd_drift", "Cd_driftTol", "Cl_drift", "Cl_driftTol",
             "stationary", "iters_to_stationary", "windowRule",
             "W_perRun", "iters_to_stationary_perRun", "stationary_perRun")


def compare(rec: dict, ref: dict, tol_cd: float = TOL_CD,
            oscillatory: bool = False, field: dict | None = None) -> dict:
    """Solver-to-solver comparison and the common pass conditions.
    oscillatory: the averaged force criterion of D-042 and its addendum
    (stationary window means, Cd within max(2 %, 0.002), Cl within
    max(2 %, 0.01)) plus the mean-field delta checks of `field`
    (mean_field_comparison: volume RMS |dUMean|/U_inf and |dpMean|/p_ref
    <= 0.02 each; missing = fail)."""
    if oscillatory:
        return _compare_mean(rec, ref, field)
    out: dict = {"tolCd": tol_cd, "window": WINDOW, "tolWindow": TOL}
    if rec.get("Cd") is not None and ref.get("Cd") is not None:
        out["CdRef"] = ref["Cd"]
        out["ClRef"] = ref["Cl"]
        out["CdRelDiff"] = abs(rec["Cd"] - ref["Cd"]) / abs(ref["Cd"])
        out["ClRelDiff"] = abs(rec["Cl"] - ref["Cl"]) / max(abs(ref["Cl"]),
                                                            1e-12)
    if rec.get("wallToConv_s") and ref.get("wallToConv_s"):
        out["speedupWall"] = ref["wallToConv_s"] / rec["wallToConv_s"]
    if rec.get("cpuHoursToConv") and ref.get("cpuHoursToConv"):
        out["speedupCpu"] = ref["cpuHoursToConv"] / rec["cpuHoursToConv"]
    out["checks"] = {
        "rc": rec["rc"] == 0 and not rec.get("fpeTrap"),
        "converged": rec.get("itersToConv") is not None,
        "referenceConverged": ref.get("itersToConv") is not None,
        "Cd": out.get("CdRelDiff") is not None and out["CdRelDiff"] <= tol_cd,
        "staticSet": rec.get("staticFraction") is not None
        and rec["staticFraction"] <= MAX_STATIC_FRACTION,
        "rollbacks": rec.get("rollbacks") == 0,
        "peakRSS": rec.get("peakRSS_MB_sum") is not None
        or rec.get("peakRSS_GB_sum") is not None,
    }
    out["pass"] = all(out["checks"].values())
    return out


def _speedups(out: dict, rec: dict, ref: dict) -> None:
    """Speed-ups simpleFoam / coupledFoam to the common-window stationary
    point (D-068) and, for information, to the per-run-window points."""
    for key, sfx in (("", ""), ("_perRun", "_perRun")):
        w, c = f"wallToConv_s{key}", f"cpuHoursToConv{key}"
        if rec.get(w) and ref.get(w):
            out[f"speedupWall{sfx}"] = ref[w] / rec[w]
        if rec.get(c) and ref.get(c):
            out[f"speedupCpu{sfx}"] = ref[c] / rec[c]


def _compare_mean(rec: dict, ref: dict, field: dict | None = None) -> dict:
    """D-042 comparison: both runs stationary, window means within the
    user-approved tolerances (run_bench.mean_comparison); mean-field delta
    checks from `field` (D-042 addendum, proposed thresholds)."""
    mc = run_bench.mean_comparison(rec, ref)
    out: dict = {
        "criterion": run_bench.STAT_CRITERION,
        "tolCd": run_bench.OSC_TOL["Cd"], "tolCl": run_bench.OSC_TOL["Cl"],
        "meanComparison": mc,
        "stationaryMean": {
            rec.get("solver", "coupledFoam"):
            {k: rec.get(k) for k in STAT_KEYS},
            ref.get("solver", "simpleFoam"):
            {k: ref.get(k) for k in STAT_KEYS},
        },
    }
    if rec.get("Cd_mean") is not None and ref.get("Cd_mean") is not None:
        out.update({
            "CdRef": ref["Cd_mean"], "ClRef": ref["Cl_mean"],
            "CdRefStd": ref.get("Cd_std"), "ClRefStd": ref.get("Cl_std"),
            "CdRelDiff": mc.get("Cd_relDiff"), "ClRelDiff": mc.get("Cl_relDiff"),
            "CdAbsDiff": mc.get("Cd_absDiff"), "ClAbsDiff": mc.get("Cl_absDiff"),
            "CdTol": mc.get("Cd_tol"), "ClTol": mc.get("Cl_tol"),
        })
    _speedups(out, rec, ref)
    field = field or {}
    fchk = field.get("fieldChecks") or {}
    out.update({
        "fieldCompare": field or None,
        "fieldTol": dict(run_bench.FIELD_TOL),
        "fieldRmsU": field.get("volRmsMagUDeltaRel"),
        "fieldRmsP": field.get("volRmsPDeltaRel"),
        "fieldMaxU": field.get("volMaxMagUDeltaRel"),
        "fieldCellFracU": field.get("cellFractionMagUDeltaAbove"),
    })
    out["checks"] = {
        "fieldU": field.get("rc") == 0
        and bool(fchk.get("volRmsMagUDeltaRel")),
        "fieldP": field.get("rc") == 0 and bool(fchk.get("volRmsPDeltaRel")),
        "rc": rec["rc"] == 0 and not rec.get("fpeTrap"),
        "converged": bool(rec.get("stationary"))
        and rec.get("iters_to_stationary") is not None,
        "referenceConverged": bool(ref.get("stationary"))
        and ref.get("iters_to_stationary") is not None,
        "Cd": bool(mc.get("Cd_pass")),
        "Cl": bool(mc.get("Cl_pass")),
        "staticSet": rec.get("staticFraction") is not None
        and rec["staticFraction"] <= MAX_STATIC_FRACTION,
        "rollbacks": rec.get("rollbacks") == 0,
        "peakRSS": rec.get("peakRSS_MB_sum") is not None
        or rec.get("peakRSS_GB_sum") is not None,
    }
    out["pass"] = all(out["checks"].values())
    return out


def assert_checks(cmp: dict, rec: dict, ref: dict) -> None:
    c = cmp["checks"]
    crit = cmp.get("criterion", "12.3(ii)")
    assert c["rc"], f"coupledFoam failed (rc {rec['rc']}, fpe {rec.get('fpeTrap')})"
    assert c["referenceConverged"], \
        f"simpleFoam reference did not meet {crit} within its budget"
    assert c["converged"], f"coupledFoam did not meet {crit}"
    assert c["Cd"], (rec.get("Cd"), ref.get("Cd"), cmp.get("CdRelDiff"))
    if "Cl" in c:
        assert c["Cl"], (rec.get("Cl"), ref.get("Cl"), cmp.get("ClAbsDiff"))
    if "fieldU" in c:
        fc = cmp.get("fieldCompare") or {}
        assert c["fieldU"], ("mean-field RMS |dU|/U_inf", cmp.get("fieldRmsU"),
                             cmp.get("fieldTol"), fc.get("error"),
                             fc.get("log"))
        assert c["fieldP"], ("mean-field RMS |dp|/p_ref", cmp.get("fieldRmsP"),
                             cmp.get("fieldTol"), fc.get("error"),
                             fc.get("log"))
    assert c["staticSet"], (rec.get("staticCells"), rec.get("staticFraction"))
    assert c["rollbacks"], rec.get("rollbacks")
    assert c["peakRSS"], "no peak RSS recorded"


# ---------------------------------------------------------------------------
# Test
# ---------------------------------------------------------------------------

@pytest.mark.heavy
@pytest.mark.parametrize("variant", ["a", "b"])
def test_T4(foam, variant):
    mesh_args = ["-mesh", variant]
    mesh = mesh_dir(TEMPLATE, f"T4{variant}_mesh", mesh_args, NP)
    budget = BUDGET[variant]
    # wake case: averaged force criterion (D-042)
    osc = run_bench.is_oscillatory(f"T4{variant}")

    ref_case, ref = reference(
        TEMPLATE, mesh, f"ref_T4{variant}_np{NP}", mesh_args,
        budget_sets("simpleFoam", budget["simpleFoam"]), oscillatory=osc,
        case_name=f"T4{variant}")

    name = f"T4{variant}_np{NP}"
    case, rec = run_solver(
        TEMPLATE, mesh, name, "coupledFoam", mesh_args,
        budget_sets("coupledFoam", budget["coupledFoam"]), oscillatory=osc,
        case_name=f"T4{variant}")

    # mean-field delta comparison (D-042 addendum); a cached reference
    # without mean fields is continued once with averaging
    field = (mean_field_comparison(case, rec, ref_case, ref, f"T4{variant}")
             if osc else None)
    cmp = compare(rec, ref, oscillatory=osc, field=field)
    meshing = cfcase.RUN_ROOT / f"T4{variant}_mesh" / "meshing.json"
    rec.update(cmp)
    rec.update({
        "case": f"T4{variant}", "meshDir": str(mesh),
        "meshing": json.loads(meshing.read_text()) if meshing.exists() else None,
        "reference": ref,
    })
    results.write("tests", name, rec)
    assert_checks(cmp, rec, ref)

"""T4 - motorBike, kOmegaSST, snappyHexMesh (spec 13). Heavy, 16 ranks.

Two meshes: (a) tutorial refinement (~350 k cells), (b) refinement levels +1
on the motorBike surface and the refinement box (target 3-5 M cells; the
actual count is recorded).

Pass (per mesh): converges by criterion 12.3(ii) with window 100 and tol
0.002 on Cd and Cl, evaluated with the benchmark harness function
bench/run_bench.py:iters_to_conv; Cd within 1 % of simpleFoam (SIMPLEC) on
the same mesh (window means at convergence); static remediation set
<= 1 % of the cells; nRollbacks == 0. Peak RSS, wall time and CPU-hours of
both solvers are recorded (every rank runs under bench/rank_wrapper.sh).

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
# Ranks of the heavy cases; CF_HEAVY_NP (e.g. 6) next to another job
NP = int(os.environ.get("CF_HEAVY_NP", "16"))
WINDOW = run_bench.WINDOW           # 100
TOL = run_bench.TOL                 # 0.002
TOL_CD = 0.01
MAX_STATIC_FRACTION = 0.01
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


def force_history(case: Path) -> dict[str, list[float]]:
    try:
        fc = post.force_coeffs(case)
    except (FileNotFoundError, KeyError):
        return {}
    return {"Cd": fc["Cd"].tolist(), "Cl": fc["Cl"].tolist()}


def run_solver(template: str, mesh: Path, name: str, solver: str,
               mesh_args: list[str], sets: dict | None = None,
               fpe: bool | None = None,
               extra_env: dict | None = None) -> tuple[Path, dict]:
    """Run one solver on a copy of a cached mesh; returns the case and a
    record with convergence (harness criterion), coefficients, wall time,
    CPU-hours and peak RSS (per-rank /usr/bin/time -v)."""
    info = mesh_info(mesh)
    nprocs = int(info.get("np", "1"))
    case = case_from_mesh(template, mesh, name, sets)
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
    }

    hist = force_history(case)
    rec["iterationsRun"] = min((len(h) for h in hist.values()), default=0)
    it = run_bench.iters_to_conv(hist) if rec["iterationsRun"] else None
    rec["itersToConv"] = it
    if hist:
        rec["CdFinal"] = hist["Cd"][-1]
        rec["ClFinal"] = hist["Cl"][-1]
    if it is not None:
        rec["Cd"] = float(np.mean(hist["Cd"][it - WINDOW:it]))
        rec["Cl"] = float(np.mean(hist["Cl"][it - WINDOW:it]))
        # Cd alone with the post.py implementation (cross-check)
        rec["itersToConvCdOnly"] = post.window_converged(
            np.asarray(hist["Cd"]), WINDOW, TOL)
        frac = run_bench.progress_fraction(case, solver, it)
        rec["progressFraction"] = frac
        if frac is not None and ranks:
            rec["wallToConv_s"] = ranks["wallSeconds"] * frac
            rec["cpuHoursToConv"] = ranks["cpuHours"] * frac

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


def reference(template: str, mesh: Path, name: str, mesh_args: list[str],
              sets: dict | None = None) -> tuple[Path, dict]:
    """simpleFoam (SIMPLEC) reference on the same mesh, cached in
    run/<name> together with its record (reference.json)."""
    case = cfcase.RUN_ROOT / name
    meta = case / "reference.json"
    if cfcase.solver_ok(case, "simpleFoam") and meta.exists():
        rec = json.loads(meta.read_text())
        if rec.get("meshCells") == mesh_cells(mesh):
            return case, rec
    case, rec = run_solver(template, mesh, name, "simpleFoam", mesh_args,
                           sets, fpe=False)
    meta.write_text(json.dumps(results._clean(rec), indent=2) + "\n")
    assert rec["rc"] == 0 and cfcase.solver_ok(case, "simpleFoam"), \
        f"simpleFoam reference {name} failed"
    return case, rec


def compare(rec: dict, ref: dict, tol_cd: float = TOL_CD) -> dict:
    """Solver-to-solver comparison and the common pass conditions."""
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


def assert_checks(cmp: dict, rec: dict, ref: dict) -> None:
    c = cmp["checks"]
    assert c["rc"], f"coupledFoam failed (rc {rec['rc']}, fpe {rec.get('fpeTrap')})"
    assert c["referenceConverged"], \
        "simpleFoam reference did not meet 12.3(ii) within its budget"
    assert c["converged"], "coupledFoam did not meet 12.3(ii)"
    assert c["Cd"], (rec.get("Cd"), ref.get("Cd"), cmp.get("CdRelDiff"))
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

    ref_case, ref = reference(
        TEMPLATE, mesh, f"ref_T4{variant}_np{NP}", mesh_args,
        budget_sets("simpleFoam", budget["simpleFoam"]))

    name = f"T4{variant}_np{NP}"
    case, rec = run_solver(
        TEMPLATE, mesh, name, "coupledFoam", mesh_args,
        budget_sets("coupledFoam", budget["coupledFoam"]))

    cmp = compare(rec, ref)
    meshing = cfcase.RUN_ROOT / f"T4{variant}_mesh" / "meshing.json"
    rec.update(cmp)
    rec.update({
        "case": f"T4{variant}", "meshDir": str(mesh),
        "meshing": json.loads(meshing.read_text()) if meshing.exists() else None,
        "reference": ref,
    })
    results.write("tests", name, rec)
    assert_checks(cmp, rec, ref)

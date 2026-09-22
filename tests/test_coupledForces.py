"""D6 unit test Test-coupledForces (amendment D6, DECISIONS.md D-065).

On T3 airFoil2D (kOmegaSST):

1. simpleFoam runs N_ITERS iterations with the native forceCoeffs and the
   coupledForcesFO function object (libcoupledFoam) side by side: the two
   histories must agree at every iteration to TOL relative to the scale of
   the history (max |a - b| / max |b|) and at the final iteration to TOL
   relative (the function-object path in a native solver).
2. Test-coupledForces on the final state, serial and on NP_PAR ranks:
   coupledForces equals native forceCoeffs to TOL for Cd, Cl and Cm (native
   CmPitch), evaluated on the same fields in one process.

TOL = 1e-6 (spec D6). In SP (CF_PRECISION=sp) the same test runs on the
shifted SP case; the native object then sums in float.
"""

from __future__ import annotations

import json

import numpy as np
import pytest

from cflib import case as cfcase
from cflib import env as cfenv
from cflib import post, precision, results

TOL = 1e-6
N_ITERS = 200
NP_PAR = 4


def _app(case, nprocs: int, tag: str) -> dict:
    out = case / f"coupledForces_{tag}.json"
    cmd = ["Test-coupledForces", "-latestTime", "-json", str(out)]
    if nprocs > 1:
        cmd = cfenv.mpirun_prefix(nprocs) + cmd + ["-parallel"]
    rc = cfenv.run(cmd, cwd=case, log=case / f"log.Test-coupledForces.{tag}")
    d = json.loads(out.read_text()) if out.exists() else {}
    d["rc"] = rc
    return d


@pytest.mark.unit
def test_coupledForces(foam):
    case = cfcase.prepare("T3_airFoil2D", "unit_T3_forces", {
        "system/controlDict": {"endTime": N_ITERS, "writeInterval": N_ITERS},
        "system/fvSolution": {"SIMPLE.residualControl": "{}"},
    })
    rc = cfcase.allrun(case, ["-solver", "simpleFoam",
                              "-turbulence", "kOmegaSST"], fpe=True)
    assert rc == 0 and cfcase.solver_ok(case, "simpleFoam"), \
        f"simpleFoam failed, see {case}/log.simpleFoam"

    # 1. function-object histories, row by row
    cf = post.coupled_force_coeffs(case)
    nat = post.force_coeffs_native(case)
    common, ic, inat = np.intersect1d(cf["Time"], nat["Time"],
                                      return_indices=True)
    hist = {"nCommon": int(common.size), "nCoupled": int(cf["Time"].size),
            "nNative": int(nat["Time"].size)}
    failures = []
    if common.size < N_ITERS - 1:
        failures.append(f"only {common.size} common iterations")
    for k, kn in (("Cd", "Cd"), ("Cl", "Cl"), ("Cm", "CmPitch")):
        a, b = cf[k][ic], nat[kn][inat]
        scale = float(np.max(np.abs(b))) if b.size else 0.0
        hist[f"{k}_maxDiffOverScale"] = (float(np.max(np.abs(a - b))) / scale
                                         if scale > 0 else None)
        hist[f"{k}_finalRelDiff"] = (float(abs(a[-1] - b[-1]) / abs(b[-1]))
                                     if b.size and b[-1] != 0 else None)
        for key in (f"{k}_maxDiffOverScale", f"{k}_finalRelDiff"):
            if hist[key] is None or not hist[key] < TOL:
                failures.append(f"history {key} = {hist[key]}")

    # 2. same fields, one process: serial and parallel
    runs = {"np1": _app(case, 1, "np1")}
    cfcase.set_entry(case, "system/decomposeParDict", "numberOfSubdomains",
                     NP_PAR)
    rcd = cfenv.run(["decomposePar", "-force", "-latestTime"], cwd=case,
                    log=case / "log.decomposePar")
    runs[f"np{NP_PAR}"] = (_app(case, NP_PAR, f"np{NP_PAR}") if rcd == 0
                           else {"rc": rcd, "error": "decomposePar failed"})
    for tag, d in runs.items():
        if d.get("rc") != 0 or not d.get("pass"):
            failures.append(f"Test-coupledForces {tag}: rc {d.get('rc')}, "
                            f"maxRelDiff {d.get('maxRelDiff')}")
    # serial vs parallel coupledForces (double sums): recorded
    if all(runs[t].get("Cd") is not None for t in runs):
        a, b = runs["np1"], runs[f"np{NP_PAR}"]
        runs["crossRankRelDiff"] = max(
            abs(a[k] - b[k]) / max(abs(a[k]), 1e-300) for k in ("Cd", "Cl", "Cm"))

    rec = {"case": str(case), "iterations": N_ITERS, "tol": TOL,
           "history": hist, "runs": runs, "failures": failures,
           "precision": precision.label(),
           "meshShift": precision.read_mesh_shift(case),
           "pass": not failures}
    results.write("tests", "Test-coupledForces", rec)
    assert not failures, failures

"""T2 - backwardFacingStep2D, k-omega SST: reattachment (spec 13).

Pass: reattachment length within 2 % of simpleFoam; at least 2x fewer outer
iterations than simpleFoam to R < 1e-5 (D-024: for simpleFoam the first
iteration at which all initial residuals p, Ux, Uy, k, omega are < 1e-5).
Wall and CPU time of both solvers are recorded (the user's primary metric).

Reattachment: wall shear stress on lowerWall; x_r is the downstream-most
point where tau_x changes to the sign it has at the outlet, linearly
interpolated; the corner eddy at the step foot is excluded that way.
"""

from __future__ import annotations

import numpy as np
import pytest

from cflib import logs, post, precision, refcase, results

STEP_HEIGHT = 0.0127
TOL_REATTACH = 0.02
ITER_FACTOR = 2.0
R_TARGET = 1e-5


def reattachment(case) -> float:
    base = case / "postProcessing" / "lowerWallShear"
    tdir = max((d for d in base.iterdir()), key=lambda d: float(d.name))
    f = next(tdir.rglob("*wallShearStress*"))
    d = np.loadtxt(f, comments="#")
    # x in the DP frame (an SP run is shifted to the origin, D5.2): the
    # step is at x = 0
    x, tx = post.unshift(case, d[:, 0], "x"), d[:, 3]
    order = np.argsort(x)
    x, tx = x[order], tx[order]
    down = x > 0
    x, tx = x[down], tx[down]
    attached = np.sign(tx[-1])
    xr = None
    for i in range(len(x) - 1):
        if np.sign(tx[i]) != attached and np.sign(tx[i + 1]) == attached:
            # linear zero crossing
            xr = x[i] - tx[i] * (x[i + 1] - x[i]) / (tx[i + 1] - tx[i])
    if xr is None:
        raise AssertionError("no reattachment point found")
    return float(xr)


@pytest.mark.case
def test_T2(foam, nprocs):
    ref, ref_rec = refcase.reference("T2_backwardFacingStep2D", "ref_T2", [])
    xr_ref = reattachment(ref)
    native = logs.parse_native(ref / "log.simpleFoam")
    it_native = logs.native_iterations_to(
        native, R_TARGET, fields=("p", "Ux", "Uy", "k", "omega"))
    # simpleFoam may never reach R_TARGET within its budget: then its
    # iteration count is a lower bound of the iterations it needs, and the
    # ratio computed from it a lower bound of the true ratio (D-046)
    native_reached = it_native is not None
    if not native_reached:
        it_native = native.get("iterations")

    name = f"T2_np{nprocs}"
    case, rec = refcase.coupled("T2_backwardFacingStep2D", name,
                                ["-np", str(nprocs)])
    xr = reattachment(case)
    R = rec["history"]["R"]
    it_coupled = next((i + 1 for i, r in enumerate(R) if r < R_TARGET), None)

    rec.update({
        "nProcs": nprocs,
        "xr": xr, "xrRef": xr_ref,
        "xr_over_h": xr / STEP_HEIGHT, "xrRef_over_h": xr_ref / STEP_HEIGHT,
        "xrRelDiff": abs(xr - xr_ref) / xr_ref, "tolXr": TOL_REATTACH,
        "iterationsToR_coupled": it_coupled,
        "iterationsToR_native": it_native,
        "nativeReachedR": native_reached,
        "iterationRatioIsLowerBound": not native_reached,
        "iterationRatio": (it_native / it_coupled
                           if it_native and it_coupled else None),
        "requiredRatio": ITER_FACTOR,
        "reference": ref_rec,
    })
    precision.annotate(rec, case, {"xr": xr})
    passed = (
        rec["rc"] == 0 and not rec["fpeTrap"]
        and rec["xrRelDiff"] < TOL_REATTACH
        and rec["iterationRatio"] is not None
        and rec["iterationRatio"] >= ITER_FACTOR
        and rec["nPseudoInverse"] == 0
    )
    rec["pass"] = passed
    results.write("tests", name, rec)

    assert rec["rc"] == 0 and not rec["fpeTrap"], \
        (rec.get("failure"), rec.get("logTail"))
    assert rec["nPseudoInverse"] == 0, rec["nPseudoInverse"]  # C2
    assert rec["xrRelDiff"] < TOL_REATTACH, (xr, xr_ref)
    assert rec["iterationRatio"] is not None, (it_coupled, it_native)
    assert rec["iterationRatio"] >= ITER_FACTOR, rec["iterationRatio"]

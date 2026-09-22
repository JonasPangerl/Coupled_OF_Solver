"""T-fpe - torture test on T1 (spec 13).

Init U = 0 everywhere, potentialInit no, CFL0 200, no start-up ramp:
coupled.startupMode none (full second-order convection and the full
pseudo-time step from iteration 1). startupUpwindIters 0 alone had no
effect: it is read only in startupMode upwind, while the default is the
hybrid beta ramp (D-048; review m12, D-068).
Pass (Debug build, traps on): no trap; the run recovers (rollbacks >= 1
allowed) and converges to R < 1e-5. Pass (Opt build): identical converged
dp to 1e-4 relative to the Debug result.

Each build writes results/tests/T-fpe_<WM_OPTIONS>.json; the Opt/Debug
comparison is made by whichever run comes second.
"""

from __future__ import annotations

import os

import pytest

from cflib import post, refcase, results

R_TARGET = 1e-5
TOL_DP_BUILDS = 1e-4
MAX_ITERS = 1000


def pressure_drop(case) -> float:
    pin = post.surface_value(case, "inletP")
    pout = post.surface_value(case, "outletP")
    key = [k for k in pin if k.startswith("areaAverage")][0]
    return float(pin[key][-1] - pout[key][-1])


@pytest.mark.case
def test_fpe(foam):
    options = os.environ.get("WM_OPTIONS", "unknown")
    debug = options.endswith("Debug")
    sets = {
        "system/fvSolution": {
            "coupled.potentialInit": "no",
            "coupled.ptc.CFL0": 200,
            # no start-up ramp at all (the hybrid default would ramp)
            "coupled.startupMode": "none",
            "coupled.startupUpwindIters": 0,
            "coupled.maxIter": MAX_ITERS,
            "coupled.convergence.residualTol": R_TARGET,
        },
        # U = 0 everywhere including the inlet start value of the interior
        "0.orig/U": {"internalField": "uniform (0 0 0)"},
    }
    case, rec = refcase.coupled("T1_pitzDaily", f"T-fpe_{options}", [],
                                sets=sets, fpe=True)
    R = rec["history"]["R"]
    it_conv = next((i + 1 for i, r in enumerate(R) if r < R_TARGET), None)
    dp = pressure_drop(case) if rec["rc"] == 0 else None

    rec.update({
        "build": options, "debugBuild": debug,
        "iterationsToR": it_conv, "Rtarget": R_TARGET, "dp": dp,
    })

    other_name = None
    for f in (results.RESULTS / "tests").glob("T-fpe_*.json"):
        if f.stem != f"T-fpe_{options}":
            other_name = f.stem
    if other_name and dp is not None:
        other = results.read("tests", other_name)
        if other and other.get("dp"):
            rec["otherBuild"] = other_name
            rec["dpRelDiffBuilds"] = abs(dp - other["dp"]) / abs(other["dp"])

    passed = (
        rec["rc"] == 0 and not rec["fpeTrap"]
        and it_conv is not None
        and rec.get("dpRelDiffBuilds", 0) < TOL_DP_BUILDS
    )
    rec["pass"] = passed
    results.write("tests", f"T-fpe_{options}", rec)

    assert not rec["fpeTrap"], "floating-point exception trapped"
    assert rec["rc"] == 0
    assert it_conv is not None, f"R < {R_TARGET} not reached"
    assert rec.get("dpRelDiffBuilds", 0) < TOL_DP_BUILDS, rec.get("dpRelDiffBuilds")

"""T1 - pitzDaily, k-omega SST: PTC and turbulence coupling (spec 13).

Pass: R < 1e-6 in <= 400 iterations; pressure drop inlet -> outlet within
1 % of simpleFoam (converged to 1e-8); CFL >= 100 reached, and no CFL cuts
after iteration 100; nRollbacks == 0. 1 vs 4 ranks: 1e-4 relative on dp.
Timing (wall, CPU-h) of both solvers is recorded.
"""

from __future__ import annotations

import pytest

from cflib import post, precision, refcase, results

MAX_ITERS = 400
R_TARGET = 1e-5      # user-approved relaxation 1e-6 -> 1e-5 (D-046)
TOL_DP = 0.01
CFL_MIN_REACHED = 100
NO_CUTS_AFTER = 100
TOL_CROSS = 1e-4


def pressure_drop(case) -> float:
    pin = post.surface_value(case, "inletP")
    pout = post.surface_value(case, "outletP")
    key = [k for k in pin if k.startswith("areaAverage")][0]
    return float(pin[key][-1] - pout[key][-1])


@pytest.mark.case
def test_T1(foam, nprocs):
    ref, ref_rec = refcase.reference("T1_pitzDaily", "ref_T1", [])
    dp_ref = pressure_drop(ref)

    name = f"T1_np{nprocs}"
    case, rec = refcase.coupled("T1_pitzDaily", name, ["-np", str(nprocs)])
    dp = pressure_drop(case)

    R = rec["history"]["R"]
    cfl = rec["history"]["CFL"]
    cuts = rec["history"]["cuts"]
    it_conv = next((i + 1 for i, r in enumerate(R) if r < R_TARGET), None)
    cfl_max = max(cfl) if cfl else 0
    cuts_late = sum(c for i, c in enumerate(cuts) if i + 1 > NO_CUTS_AFTER)

    rec.update({
        "nProcs": nprocs, "iterationsToR": it_conv, "Rtarget": R_TARGET,
        "maxIters": MAX_ITERS, "dp": dp, "dpRef": dp_ref,
        "dpRelDiff": abs(dp - dp_ref) / abs(dp_ref), "tolDp": TOL_DP,
        "CFLmaxReached": cfl_max, "cutsAfter100": cuts_late,
        "reference": ref_rec,
    })
    precision.annotate(rec, case, {"dp": dp})
    if nprocs > 1:
        serial = results.read("tests", "T1_np1")
        if serial:
            rec["crossRankDp"] = abs(dp - serial["dp"]) / abs(serial["dp"])

    passed = (
        rec["rc"] == 0 and not rec["fpeTrap"]
        and it_conv is not None and it_conv <= MAX_ITERS
        and rec["dpRelDiff"] < TOL_DP
        and cfl_max >= CFL_MIN_REACHED and cuts_late == 0
        and rec["rollbacks"] == 0
        and rec["nPseudoInverse"] == 0
        and rec.get("crossRankDp", 0) < TOL_CROSS
    )
    rec["pass"] = passed
    results.write("tests", name, rec)

    assert rec["rc"] == 0 and not rec["fpeTrap"]
    assert it_conv is not None and it_conv <= MAX_ITERS, \
        f"R < {R_TARGET} not reached in {MAX_ITERS} (final {rec['finalR']})"
    assert rec["dpRelDiff"] < TOL_DP, (dp, dp_ref)
    assert cfl_max >= CFL_MIN_REACHED, cfl_max
    assert cuts_late == 0, cuts_late
    assert rec["rollbacks"] == 0
    assert rec["nPseudoInverse"] == 0, rec["nPseudoInverse"]  # C2
    assert rec.get("crossRankDp", 0) < TOL_CROSS

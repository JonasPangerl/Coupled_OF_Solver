"""T-restart - exact resume (spec 13), on T1 and T3-SST.

Run 150 iterations and write; restart from 150 to convergence; compare with
the uninterrupted run: iteration count to convergence identical +-2; final
dp (T1) or Cd/Cl (T3) identical to 1e-5 relative; phi consistency on restart
(spec 10) < 1e-6.
"""

from __future__ import annotations

import pytest

from cflib import case as cfcase
from cflib import env as cfenv
from cflib import logs, post, refcase, results

SPLIT_AT = 150
TOL_ITERS = 2
TOL_VALUE = 1e-5
TOL_PHI = 1e-6

CASES = {
    "T1": ("T1_pitzDaily", []),
    "T3-SST": ("T3_airFoil2D", ["-turbulence", "kOmegaSST"]),
}


def quantities(tag: str, case) -> dict:
    if tag == "T1":
        pin = post.surface_value(case, "inletP")
        pout = post.surface_value(case, "outletP")
        key = [k for k in pin if k.startswith("areaAverage")][0]
        return {"dp": float(pin[key][-1] - pout[key][-1])}
    fc = post.force_coeffs(case)
    return {"Cd": float(fc["Cd"][-1]), "Cl": float(fc["Cl"][-1])}


@pytest.mark.case
@pytest.mark.parametrize("tag", list(CASES))
def test_restart(foam, tag):
    template, args = CASES[tag]

    full, full_rec = refcase.coupled(template, f"restart_{tag}_full", args)
    q_full = quantities(tag, full)

    # First part: stop at SPLIT_AT (maxIter), write state
    part = cfcase.prepare(template, f"restart_{tag}_split", {
        "system/fvSolution": {"coupled.maxIter": SPLIT_AT},
    })
    rc1 = cfcase.allrun(part, ["-solver", "coupledFoam"] + args)
    t1 = cfenv.last_timing.as_dict()
    # Second part: original maxIter, continue from the latest time
    maxiter = cfcase.get_entry(cfcase.CASES / template, "system/fvSolution",
                               "coupled.maxIter")
    cfcase.set_entry(part, "system/fvSolution", "coupled.maxIter", maxiter)
    rc2 = cfcase.allrun(part, ["-solver", "coupledFoam", "-restart"] + args)
    t2 = cfenv.last_timing.as_dict()

    rows2 = logs.parse_cf(part / "log.coupledFoam.restart")
    text2 = (part / "log.coupledFoam.restart").read_text(errors="replace")
    summ2 = logs.coupled_summary(part)
    q_split = quantities(tag, part)

    it_full = full_rec["iterations"]
    it_split = int(rows2[-1]["iter"]) if rows2 else None
    rel = {k: abs(q_split[k] - q_full[k]) / abs(q_full[k]) for k in q_full}
    phi_c = summ2.get("phiConsistency")

    rec = {
        "tag": tag, "rc": [full_rec["rc"], rc1, rc2],
        "iterationsFull": it_full, "iterationsSplit": it_split,
        "quantitiesFull": q_full, "quantitiesSplit": q_split,
        "relDiff": rel, "tolValue": TOL_VALUE,
        "phiConsistency": phi_c, "tolPhi": TOL_PHI,
        "restartedFromState": "restart from" in text2,
        "timingFull": full_rec["timingAllrun"], "timingPart1": t1,
        "timingPart2": t2,
    }
    passed = (
        all(r == 0 for r in rec["rc"])
        and it_split is not None and abs(it_split - it_full) <= TOL_ITERS
        and all(v < TOL_VALUE for v in rel.values())
        and phi_c is not None and 0 <= phi_c < TOL_PHI
        and rec["restartedFromState"]
    )
    rec["pass"] = passed
    results.write("tests", f"T-restart_{tag}", rec)

    assert all(r == 0 for r in rec["rc"]), rec["rc"]
    assert rec["restartedFromState"], "coupledState was not read"
    assert it_split is not None and abs(it_split - it_full) <= TOL_ITERS, \
        (it_split, it_full)
    assert all(v < TOL_VALUE for v in rel.values()), rel
    assert phi_c is not None and 0 <= phi_c < TOL_PHI, phi_c

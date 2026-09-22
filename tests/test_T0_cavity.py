"""T0 - lid-driven cavity 128x128, laminar, Re 100 and 1000 (spec 13).

Pass: R < 1e-8 within 300 iterations; centreline u(y), v(x) (129 points)
vs. the simpleFoam SIMPLEC reference (converged to 1e-8): L2 relative
difference < 1e-3; no FPE trap (FOAM_SIGFPE on); nClamped == 0; 1 vs 4
ranks: 1e-4 relative on the profiles.
"""

from __future__ import annotations

import pytest

from cflib import case as cfcase
from cflib import env as cfenv
from cflib import logs, post, precision, refcase, results

MAX_ITERS = 300
R_TARGET = 1e-8
TOL_PROFILE = 1e-3
TOL_CROSS = 1e-4
REYNOLDS = (100, 1000)


def _reference(re: int):
    """simpleFoam reference, serial, cached in run/ref_T0_Re<re>."""
    case = cfcase.prepare("T0_cavity", f"ref_T0_Re{re}", reuse=True)
    if not cfcase.solver_ok(case, "simpleFoam"):
        case = cfcase.prepare("T0_cavity", f"ref_T0_Re{re}")
        assert cfcase.allrun(case, ["-solver", "simpleFoam", "-Re", str(re)],
                             fpe=False) == 0
    assert cfcase.solver_ok(case, "simpleFoam")
    return case


def _profiles(case):
    """(coordinate, u) on the vertical and (coordinate, v) on the horizontal
    centreline, as written (unsorted, with duplicates in parallel)."""
    v = post.read_xy(post.sets_file(case, "centreLines", "vertical"))
    h = post.read_xy(post.sets_file(case, "centreLines", "horizontal"))
    # columns: coordinate, p, Ux, Uy, Uz; coordinates in the DP frame (an
    # SP run is shifted to the origin, D5.2)
    return ((post.unshift(case, v[:, 0], "y"), v[:, 2]),
            (post.unshift(case, h[:, 0], "x"), h[:, 3]))


def _on(profile, reference):
    """Profile and reference at their common sample points: sorted,
    de-duplicated; points a parallel sets run dropped are left out (see
    post.match_profiles). Returns (values, reference values, number of
    reference points missing in the profile)."""
    vals, ref, n_missing = post.match_profiles(profile, reference)
    assert vals.size > 0, "no common sample points"
    return vals, ref, n_missing


@pytest.mark.case
@pytest.mark.parametrize("re", REYNOLDS, ids=[f"Re{r}" for r in REYNOLDS])
def test_T0(foam, re, nprocs):
    ref = _reference(re)
    ref_u, ref_v = _profiles(ref)
    ref_log = logs.parse_native(ref / "log.simpleFoam")

    name = f"T0_Re{re}_np{nprocs}"
    case = cfcase.prepare("T0_cavity", name, {
        "system/fvSolution": {"coupled.maxIter": MAX_ITERS},
    })
    rc = cfcase.allrun(case, ["-solver", "coupledFoam", "-Re", str(re),
                              "-np", str(nprocs)], fpe=True)
    timing = cfenv.last_timing.as_dict()
    text = (case / "log.coupledFoam").read_text(errors="replace")
    rows = logs.parse_cf(case / "log.coupledFoam")
    summ = logs.coupled_summary(case)

    fpe_trap = logs.fpe_trapped(case / "log.coupledFoam")
    it_conv = logs.coupled_iterations_to(rows, R_TARGET)
    n_clamped = max((r.get("nClamped", 0) for r in rows), default=-1)

    prof_u, prof_v = _profiles(case)
    u, ru, fill_u = _on(prof_u, ref_u)
    v, rv, fill_v = _on(prof_v, ref_v)
    du = post.l2rel(u, ru)
    dv = post.l2rel(v, rv)

    rec = {
        "Re": re, "nProcs": nprocs, "rc": rc,
        "iterationsToR": it_conv, "maxIters": MAX_ITERS, "Rtarget": R_TARGET,
        "finalR": rows[-1]["R"] if rows else None,
        "iterations": len(rows),
        "l2rel_u": du, "l2rel_v": dv, "tolProfile": TOL_PROFILE,
        "profilePointsMissing": {"u": fill_u, "v": fill_v},
        "fpeTrap": fpe_trap, "fpeEnabled": "trapFpe" in text,
        "nClampedMax": n_clamped,
        "rollbacks": summ.get("rollbacks"),
        "nPseudoInverse": summ.get("nPseudoInverse", 0),
        # wall-clock AND CPU-hours of both solvers (solver loop only)
        **refcase.coupled_timing(summ, timing),
        "timingAllrun": timing,
        "peakRSS_MB_sum": summ.get("peakRSS_MB_sum"),
        "reference": {"solver": "simpleFoam", "iterations": ref_log["iterations"],
                      "convergedAt": ref_log["convergedAt"],
                      **refcase.native_timing(ref / "log.simpleFoam", 1)},
        "history": {k: [r.get(k) for r in rows]
                    for k in ("R", "CFL", "omega", "cuts", "linIters", "tIter",
                              "tWall", "eta", "rho")},
        "gamg": {k: summ.get(k) for k in ("gamgLevels", "gamgMergeLevels",
                                          "gamgCop", "gamgCellsPerLevel")},
        "cpuHoursSolver": summ.get("cpuHours"),
    }

    if nprocs > 1:
        serial = results.read("tests", f"T0_Re{re}_np1")
        if serial is not None:
            su, sv = _profiles(cfcase.run_dir(f"T0_Re{re}_np1"))
            pu, psu, fu = _on(prof_u, su)
            pv, psv, fv = _on(prof_v, sv)
            rec["crossRank_u"] = post.l2rel(pu, psu)
            rec["crossRank_v"] = post.l2rel(pv, psv)
            rec["crossRankPointsMissing"] = {"u": fu, "v": fv}

    # D10 (SP only): u/v centreline L2 difference to the DP run of the same
    # name < 1e-3. Sample points matched in the DP frame to 1e-5 of the
    # line length (the SP coordinates carry float rounding).
    precision.annotate(rec, case)
    dp_case = cfcase.RUN_ROOT / precision.dp_name(case.name)
    if precision.is_sp() and (dp_case / "postProcessing").is_dir():
        du_ref, dv_ref = _profiles(dp_case)
        a, b, _ = post.match_profiles(prof_u, du_ref, rel_tol=1e-5)
        c, d, _ = post.match_profiles(prof_v, dv_ref, rel_tol=1e-5)
        rec["l2rel_u_vs_DP"] = post.l2rel(a, b) if b.size else None
        rec["l2rel_v_vs_DP"] = post.l2rel(c, d) if d.size else None
        rec["spVsDpTol"] = TOL_PROFILE
        rec["spVsDpPass"] = (None if None in (rec["l2rel_u_vs_DP"],
                                              rec["l2rel_v_vs_DP"])
                             else max(rec["l2rel_u_vs_DP"],
                                      rec["l2rel_v_vs_DP"]) < TOL_PROFILE)

    passed = (
        rc == 0 and it_conv is not None and it_conv <= MAX_ITERS
        and du < TOL_PROFILE and dv < TOL_PROFILE
        and not fpe_trap and n_clamped == 0
        and rec["nPseudoInverse"] == 0
        and rec.get("crossRank_u", 0) < TOL_CROSS
        and rec.get("crossRank_v", 0) < TOL_CROSS
    )
    rec["pass"] = passed
    results.write("tests", name, rec)

    assert rc == 0, "coupledFoam failed"
    assert not fpe_trap, "FPE trap"
    assert n_clamped == 0, "clamped coefficients"
    assert rec["nPseudoInverse"] == 0, \
        f"Rhie-Chow pseudo-inverse used {rec['nPseudoInverse']} times (C2)"
    assert it_conv is not None and it_conv <= MAX_ITERS, \
        f"R < {R_TARGET} not reached in {MAX_ITERS} iterations (final {rec['finalR']})"
    assert du < TOL_PROFILE and dv < TOL_PROFILE, (du, dv)
    assert rec.get("crossRank_u", 0) < TOL_CROSS, rec.get("crossRank_u")
    assert rec.get("crossRank_v", 0) < TOL_CROSS, rec.get("crossRank_v")

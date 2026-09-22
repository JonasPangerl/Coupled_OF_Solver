"""T3 - airFoil2D, kOmegaSST and GEKO (spec 13).

Pass (review M6, D-068): both solvers converged, and the converged
coefficients agree.
- coupledFoam converged: the spec 12.3(ii) force window (the last 100
  iterations, max - min <= 0.2 % |mean| for Cd AND Cl) holds at the end of
  the run, or the solver stopped on its own criterion (summary.json
  "converged": R below its tolerance or its own force window).
- simpleFoam reference converged: its residualControl stop, or the 12.3(ii)
  window at the end of its run.
- Cd and Cl = the means over the final 100-iteration window of each run
  (not the last sample: a limit cycle whose last sample happens to lie
  within the tolerance must not pass); within 5 % of simpleFoam
  (user-approved: 0.5 % spec -> 2 % D-046 -> 5 % D-058).
- GEKO with the same solver settings as SST (no per-model tuning); nDyn at
  the end == 0; no Rhie-Chow pseudo-inverse fallback (C2).
itersToConv (first 12.3(ii) window, or the solver's stop) of both solvers
is recorded for the speed-up (make_report._speed_record). Wall and CPU
time of both solvers are recorded.
"""

from __future__ import annotations

import sys

import pytest

from cflib import env as cfenv
from cflib import logs, post, precision, refcase, results

sys.path.insert(0, str(cfenv.REPO / "bench"))
import run_bench  # noqa: E402  (12.3(ii) window, force history)

TOL_COEFF = 0.05     # user-approved: 0.5 % -> 2 % (D-046) -> 5 % (D-058)
TOL_IDENTITY = 1e-6  # D6: coupledForces vs native forceCoeffs, same run
MODELS = ("kOmegaSST", "GEKO")
WINDOW = run_bench.WINDOW           # 100
TOL = run_bench.TOL                 # 0.002
CRITERION = (f"12.3(ii): window {WINDOW}, tol {TOL} on Cd and Cl at the end "
             "of the run, or the solver's own stop; coefficients = final-"
             "window means")


def converged_coeffs(case, solver_stop: bool) -> dict:
    """Convergence and converged coefficients of a run (see module doc).
    solver_stop: the solver ended on its own convergence criterion."""
    h = run_bench.force_history(case, original_only=False)
    n = min((len(v) for v in h.values()), default=0)
    out: dict = {"criterion": CRITERION, "iterationsRun": n,
                 "solverStop": bool(solver_stop)}
    if n == 0:
        out.update({"converged": False, "itersToConv": None})
        return out
    first = run_bench.iters_to_conv(h)
    final_ok = n >= WINDOW and all(run_bench._window_ok(h[q])
                                   for q in ("Cd", "Cl"))
    converged = final_ok or bool(solver_stop)
    w = min(WINDOW, n)
    out.update({
        "finalWindowOk": final_ok,
        "itersToConvFirstWindow": first,
        "converged": converged,
        # the convergence point: the first 12.3(ii) window, else the stop
        "itersToConv": (first if first is not None else n) if converged
        else None,
        "Cd": sum(h["Cd"][-w:]) / w, "Cl": sum(h["Cl"][-w:]) / w,
        "CdLast": h["Cd"][-1], "ClLast": h["Cl"][-1],
        "CdRange": (max(h["Cd"][-w:]) - min(h["Cd"][-w:])),
        "ClRange": (max(h["Cl"][-w:]) - min(h["Cl"][-w:])),
    })
    return out


@pytest.mark.case
@pytest.mark.parametrize("model", MODELS)
def test_T3(foam, model, nprocs):
    ref, ref_rec = refcase.reference(
        "T3_airFoil2D", f"ref_T3_{model}", ["-turbulence", model])
    ref_rec = dict(ref_rec)
    ref_rec.update(converged_coeffs(
        ref, ref_rec.get("convergedAt") is not None))

    name = f"T3_{model}_np{nprocs}"
    case, rec = refcase.coupled(
        "T3_airFoil2D", name, ["-turbulence", model, "-np", str(nprocs)])
    summ = logs.coupled_summary(case)
    conv = converged_coeffs(case, summ.get("converged") is True)
    ndyn = rec["history"]["nDyn"]
    cd, cl = conv.get("Cd"), conv.get("Cl")
    cd_ref, cl_ref = ref_rec.get("Cd"), ref_rec.get("Cl")

    rec.update(conv)
    # D6: the forces come from coupledForces (double accumulation) when the
    # run wrote them; in DP the native forceCoeffs of the same run must
    # agree (same fields, same iteration) to TOL_IDENTITY at the final
    # iteration
    ident = post.force_identity(case)
    ident_ref = post.force_identity(ref)
    ref_rec.update(forceSource=post.force_source(ref),
                   forceIdentity=ident_ref)
    rec.update({
        "model": model, "nProcs": nprocs,
        "CdRef": cd_ref, "ClRef": cl_ref,
        "CdRelDiff": (abs(cd - cd_ref) / abs(cd_ref)
                      if None not in (cd, cd_ref) else None),
        "ClRelDiff": (abs(cl - cl_ref) / abs(cl_ref)
                      if None not in (cl, cl_ref) else None),
        "tol": TOL_COEFF,
        "nDynFinal": ndyn[-1] if ndyn else None,
        "forceSource": post.force_source(case),
        "forceIdentity": ident,
        "tolIdentity": TOL_IDENTITY,
        "reference": ref_rec,
    })
    precision.annotate(rec, case, {"Cd": cd, "Cl": cl})
    # In SP the native object accumulates in float: the difference is the
    # effect D6 removes, recorded but not a criterion there
    identity_ok = precision.is_sp() or all(
        d.get(f"{k}_finalRelDiff", 0.0) < TOL_IDENTITY
        for d in (ident, ident_ref) for k in ("Cd", "Cl", "Cm"))
    rec["forceIdentityPass"] = None if precision.is_sp() else identity_ok
    checks = {
        "rc": rec["rc"] == 0 and not rec["fpeTrap"],
        "converged": bool(rec["converged"]),
        "referenceConverged": bool(ref_rec["converged"]),
        "Cd": rec["CdRelDiff"] is not None and rec["CdRelDiff"] < TOL_COEFF,
        "Cl": rec["ClRelDiff"] is not None and rec["ClRelDiff"] < TOL_COEFF,
        "nDynFinal": rec["nDynFinal"] == 0,
        "nPseudoInverse": rec["nPseudoInverse"] == 0,
        "forceIdentity": identity_ok,
    }
    rec["checks"] = checks
    rec["pass"] = all(checks.values())
    results.write("tests", name, rec)

    assert checks["rc"], (rec.get("failure"), rec.get("logTail"))
    assert rec["nPseudoInverse"] == 0, rec["nPseudoInverse"]  # C2
    assert checks["referenceConverged"], \
        f"simpleFoam reference ref_T3_{model} did not converge ({CRITERION})"
    assert checks["converged"], \
        (f"coupledFoam did not converge ({CRITERION}): final Cd range "
         f"{conv.get('CdRange')}, Cl range {conv.get('ClRange')}")
    assert checks["Cd"], (cd, cd_ref)
    assert checks["Cl"], (cl, cl_ref)
    assert checks["nDynFinal"], rec["nDynFinal"]
    assert checks["forceIdentity"], (ident, ident_ref)

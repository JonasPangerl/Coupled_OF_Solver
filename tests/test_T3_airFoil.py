"""T3 - airFoil2D, kOmegaSST and GEKO (spec 13).

Pass: Cl, Cd within 0.5 % of simpleFoam for each model; GEKO with the same
solver settings as SST (no per-model tuning); nDyn at convergence == 0.
Wall and CPU time of both solvers are recorded.
"""

from __future__ import annotations

import pytest

from cflib import post, precision, refcase, results

TOL_COEFF = 0.05     # user-approved: 0.5 % -> 2 % (D-046) -> 5 % (D-058)
TOL_IDENTITY = 1e-6  # D6: coupledForces vs native forceCoeffs, same run
MODELS = ("kOmegaSST", "GEKO")


def final_coeffs(case) -> tuple[float, float]:
    fc = post.force_coeffs(case)
    return float(fc["Cd"][-1]), float(fc["Cl"][-1])


@pytest.mark.case
@pytest.mark.parametrize("model", MODELS)
def test_T3(foam, model, nprocs):
    ref, ref_rec = refcase.reference(
        "T3_airFoil2D", f"ref_T3_{model}", ["-turbulence", model])
    cd_ref, cl_ref = final_coeffs(ref)

    name = f"T3_{model}_np{nprocs}"
    case, rec = refcase.coupled(
        "T3_airFoil2D", name, ["-turbulence", model, "-np", str(nprocs)])
    cd, cl = final_coeffs(case)
    ndyn = rec["history"]["nDyn"]

    # D6: the forces come from coupledForces (double accumulation) when the
    # run wrote them; in DP the native forceCoeffs of the same run must
    # agree (same fields, same iteration) to TOL_IDENTITY at the final
    # iteration
    ident = post.force_identity(case)
    ident_ref = post.force_identity(ref)
    ref_rec = dict(ref_rec, forceSource=post.force_source(ref),
                   forceIdentity=ident_ref)
    rec.update({
        "model": model, "nProcs": nprocs,
        "Cd": cd, "Cl": cl, "CdRef": cd_ref, "ClRef": cl_ref,
        "CdRelDiff": abs(cd - cd_ref) / abs(cd_ref),
        "ClRelDiff": abs(cl - cl_ref) / abs(cl_ref),
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
    passed = (
        rec["rc"] == 0 and not rec["fpeTrap"]
        and rec["CdRelDiff"] < TOL_COEFF and rec["ClRelDiff"] < TOL_COEFF
        and rec["nDynFinal"] == 0
        and rec["nPseudoInverse"] == 0
        and identity_ok
    )
    rec["pass"] = passed
    results.write("tests", name, rec)

    assert rec["rc"] == 0 and not rec["fpeTrap"]
    assert rec["nPseudoInverse"] == 0, rec["nPseudoInverse"]  # C2
    assert rec["CdRelDiff"] < TOL_COEFF, (cd, cd_ref)
    assert rec["ClRelDiff"] < TOL_COEFF, (cl, cl_ref)
    assert rec["nDynFinal"] == 0, rec["nDynFinal"]
    assert identity_ok, (ident, ident_ref)

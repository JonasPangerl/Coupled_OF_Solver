"""T3 - airFoil2D, kOmegaSST and GEKO (spec 13).

Pass: Cl, Cd within 0.5 % of simpleFoam for each model; GEKO with the same
solver settings as SST (no per-model tuning); nDyn at convergence == 0.
Wall and CPU time of both solvers are recorded.
"""

from __future__ import annotations

import pytest

from cflib import post, refcase, results

TOL_COEFF = 0.02     # user-approved relaxation 0.5 % -> 2 % (D-046)
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

    rec.update({
        "model": model, "nProcs": nprocs,
        "Cd": cd, "Cl": cl, "CdRef": cd_ref, "ClRef": cl_ref,
        "CdRelDiff": abs(cd - cd_ref) / abs(cd_ref),
        "ClRelDiff": abs(cl - cl_ref) / abs(cl_ref),
        "tol": TOL_COEFF,
        "nDynFinal": ndyn[-1] if ndyn else None,
        "reference": ref_rec,
    })
    passed = (
        rec["rc"] == 0 and not rec["fpeTrap"]
        and rec["CdRelDiff"] < TOL_COEFF and rec["ClRelDiff"] < TOL_COEFF
        and rec["nDynFinal"] == 0
    )
    rec["pass"] = passed
    results.write("tests", name, rec)

    assert rec["rc"] == 0 and not rec["fpeTrap"]
    assert rec["CdRelDiff"] < TOL_COEFF, (cd, cd_ref)
    assert rec["ClRelDiff"] < TOL_COEFF, (cl, cl_ref)
    assert rec["nDynFinal"] == 0, rec["nDynFinal"]

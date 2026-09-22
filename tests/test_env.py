"""Phase 0 gate: toolchain and environment (spec 3, 16).

- Python >= 3.11 with the required packages
- OpenFOAM v2606, label size 32
- Test-precision: per build sizeof(scalar)/sizeof(solveScalar) (DP 8/8,
  SP 4/4, SPDP 4/8), blockScalar 4 and reduceScalar 8 always (D-001,
  amendment D1, D-064); the binary's build matches WM_PRECISION_OPTION
- coupledFoam and the test applications are built
- native simpleFoam runs the pitzDaily tutorial (short run)
"""

from __future__ import annotations

import importlib
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

from cflib import env as cfenv
from cflib import results

APPS = ["coupledFoam", "Test-precision", "Test-block4Ops", "Test-doubleReduce",
        "Test-blockMatrix", "Test-blockGAMG"]

PACKAGES = ["numpy", "pandas", "matplotlib", "yaml", "stl", "pytest"]

# Short tutorial run: enough iterations to show the solver works
TUTORIAL_ITERS = 20


def test_env(foam, tmp_path):
    rec: dict = {"checks": {}}
    checks = rec["checks"]

    checks["python"] = sys.version.split()[0]
    checks["pythonOk"] = sys.version_info >= (3, 11)

    missing = []
    for pkg in PACKAGES:
        try:
            importlib.import_module(pkg)
        except ImportError:
            missing.append(pkg)
    checks["missingPackages"] = missing

    checks["WM_PROJECT_VERSION"] = os.environ.get("WM_PROJECT_VERSION")
    checks["WM_LABEL_SIZE"] = os.environ.get("WM_LABEL_SIZE")
    checks["WM_PRECISION_OPTION"] = os.environ.get("WM_PRECISION_OPTION")
    checks["versionOk"] = str(checks["WM_PROJECT_VERSION"]).startswith("v2606") \
        or str(checks["WM_PROJECT_VERSION"]) == "2606"
    checks["labelOk"] = checks["WM_LABEL_SIZE"] == "32"

    appbin = Path(os.environ["FOAM_USER_APPBIN"])
    checks["appsMissing"] = [a for a in APPS if not (appbin / a).exists()]

    # Precision (D-001)
    pj = tmp_path / "precision.json"
    rc = subprocess.run(["Test-precision", "-json", str(pj)],
                        capture_output=True, text=True).returncode
    checks["precisionRc"] = rc
    if pj.exists():
        import json
        checks["precision"] = json.loads(pj.read_text())
    # The binary was built for the sourced precision option (D1)
    checks["precisionMatchesEnv"] = (
        checks.get("precision", {}).get("precision")
        == os.environ.get("WM_PRECISION_OPTION"))

    # simpleFoam pitzDaily tutorial (short)
    tut = Path(os.environ["FOAM_TUTORIALS"]) / "incompressible/simpleFoam/pitzDaily"
    case = tmp_path / "pitzDaily"
    shutil.copytree(tut, case)
    subprocess.run(["foamDictionary", "-entry", "endTime", "-set",
                    str(TUTORIAL_ITERS), "system/controlDict"], cwd=case,
                   check=True, stdout=subprocess.DEVNULL)
    subprocess.run(["foamDictionary", "-entry", "functions", "-set", "{}",
                    "system/controlDict"], cwd=case, check=True,
                   stdout=subprocess.DEVNULL)
    rc_mesh = cfenv.run(["blockMesh"], cwd=case, log=case / "log.blockMesh")
    rc_sf = cfenv.run(["simpleFoam"], cwd=case, log=case / "log.simpleFoam")
    log = (case / "log.simpleFoam").read_text(errors="replace")
    checks["pitzDailyRc"] = [rc_mesh, rc_sf]
    checks["pitzDailyOk"] = rc_mesh == 0 and rc_sf == 0 and "\nEnd" in log

    passed = (
        checks["pythonOk"] and not missing and checks["versionOk"]
        and checks["labelOk"] and not checks["appsMissing"]
        and rc == 0 and checks["precisionMatchesEnv"]
        and checks["pitzDailyOk"]
    )
    rec["pass"] = passed
    results.write("tests", "test_env", rec)

    assert checks["pythonOk"], checks["python"]
    assert not missing, f"missing Python packages: {missing}"
    assert checks["versionOk"], checks["WM_PROJECT_VERSION"]
    assert checks["labelOk"]
    assert not checks["appsMissing"], f"not built: {checks['appsMissing']}"
    assert rc == 0, "Test-precision gate failed (D-001, D1)"
    assert checks["precisionMatchesEnv"], (
        "Test-precision build does not match WM_PRECISION_OPTION")
    assert checks["pitzDailyOk"], "simpleFoam pitzDaily tutorial failed"

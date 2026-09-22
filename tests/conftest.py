"""pytest configuration for coupledFoam (spec 13).

Run from the repository root in a shell with exactly one OpenFOAM v2606
environment sourced:

    pytest tests/                      # everything that is not heavy
    pytest tests/ --heavy              # also T4, T5, T-scaling
    pytest tests/ -k T0                # one case
    pytest tests/ --ranks 1            # serial variants only
    pytest tests/ --precision sp       # SP runs (amendment D5; needs the SP
                                       # environment, see cflib/precision.py)

Every test writes results/tests/<name>.json (measured quantities and
pass/fail). Case runs go to run/ (git-ignored); simpleFoam references are
cached there and reused (delete run/ref_* to recompute).

Shared machine: if other solver/meshing jobs run or the load is high, MPI
ranks are not bound to cores (CF_MPI_BIND=none) and heavy tests are skipped
unless CF_FORCE_HEAVY=1. With CF_MPI_CPUSET=10-15 (a core set disjoint from
the other job, which must itself run with --cpu-set) ranks are pinned there
and heavy tests run; use --ranks / case np not larger than the set.

Heavy cases run with CF_HEAVY_NP ranks (default 10: the user caps this
machine at 10 cores, D-031); on a free machine ranks are bound to cores.
"""

from __future__ import annotations

import os

import pytest

from cflib import env as cfenv
from cflib import precision


def pytest_addoption(parser):
    parser.addoption("--heavy", action="store_true", default=False,
                     help="run heavy tests (T4, T5, T-scaling)")
    parser.addoption("--ranks", default="1,4",
                     help="comma-separated rank counts for the case tests")
    parser.addoption("--precision", default=None, choices=("dp", "sp"),
                     help="dp (default) or sp: sets CF_PRECISION; sp needs "
                          "the SP OpenFOAM environment (amendment D5, "
                          "cflib/precision.py)")


def pytest_configure(config):
    if config.getoption("--precision"):
        os.environ["CF_PRECISION"] = config.getoption("--precision")
    config.addinivalue_line("markers", "heavy: long-running, many cores")
    config.addinivalue_line("markers", "unit: unit-test applications")
    config.addinivalue_line("markers", "case: solver case tests")


def pytest_collection_modifyitems(config, items):
    state = cfenv.machine_state()
    if state.busy and "CF_MPI_BIND" not in os.environ:
        os.environ["CF_MPI_BIND"] = "none"
    run_heavy = config.getoption("--heavy")
    # A dedicated core set (CF_MPI_CPUSET, disjoint from the other job)
    # makes heavy runs acceptable on a shared machine
    force = (os.environ.get("CF_FORCE_HEAVY") == "1"
             or bool(os.environ.get("CF_MPI_CPUSET")))
    for item in items:
        if "heavy" in item.keywords:
            if not run_heavy:
                item.add_marker(pytest.mark.skip(reason="heavy: use --heavy"))
            elif state.busy and not force:
                item.add_marker(pytest.mark.skip(
                    reason="machine busy (other jobs: %d, load %.1f)"
                    % (len(state.jobs), state.load)))


def pytest_generate_tests(metafunc):
    if "nprocs" in metafunc.fixturenames:
        ranks = [int(r) for r in metafunc.config.getoption("--ranks").split(",")]
        metafunc.parametrize("nprocs", ranks, ids=[f"np{r}" for r in ranks])


@pytest.fixture(scope="session")
def foam():
    try:
        env = cfenv.foam_env()
        precision.check_environment()
        return env
    except (RuntimeError, ValueError) as err:
        pytest.exit(str(err), returncode=2)


@pytest.fixture(scope="session")
def machine():
    return cfenv.machine_state()

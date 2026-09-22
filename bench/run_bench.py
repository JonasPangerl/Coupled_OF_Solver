#!/usr/bin/env python3
"""Benchmark harness (spec 14): native simpleFoam vs coupledFoam.

For each case and configuration the solver runs with its own stop criteria
disabled (fixed iteration budget). The harness then finds the convergence
iteration with the identical criterion for all solvers (spec 12.3 ii): over
the last `window` iterations max - min <= tol*|mean| for every monitored
quantity (Cd and Cl, or the pressure drop for T1/T2).

Exception, the wake cases marked "oscillatory" (T4a, T4b, T5; D-042 and its
addendum): their forces oscillate physically for ever, so 12.3(ii) is
unsatisfiable. There "converged" means a stationary window mean
(stationary_mean: window W = max(1000, n/2), capped at n, the two
half-window means of Cd and of Cl differ by at most max(1 % |mean|, 0.005));
the time to convergence uses the first stationary window
(iters_to_stationary, scanned in steps of 50) in place of iters_to_conv,
and Cd / Cl are the final-window means. The mean fields (fieldAverage
function object of the T4/T5 controlDict) are averaged over the same
window: set_field_average_start sets its timeStart to n - W + 1.

Configurations (DECISIONS.md D-025, amendment B10):
    A  simpleFoam, the tutorial's solver settings and relaxation factors.
       Per-case overrides where the test case differs from the tutorial:
       T3 (airFoil2D tutorial): plain SIMPLE (consistent no), p 0.3,
       U / turbulence 0.7. T5 has no tutorial: A = the motorBike tutorial
       settings (SIMPLEC, U 0.9, k/omega 0.7), as in T4.
    B  simpleFoam SIMPLEC, consistent yes, relaxation p 1.0 / U 0.9 / k,omega 0.9
    C  coupledFoam defaults (K-cycle, autoTune on, adaptive relTol)
    D  coupledFoam with preconditioner blockDiagonal (isolates the AMG gain)
    E  coupledFoam with a fixed V-cycle: cycleType V, autoTune no (autoTune
       would promote V -> F -> W). B10 scope: T2, T4b, T5.
    F  coupledFoam with adaptiveRelTol no (isolates Eisenstat-Walker).
       B10 scope: T1, T3-SST.
    G  coupledFoam with anderson enabled. B10 scope: T1, T3-SST.
    H  coupledFoam with autoTune no (fixed K-cycle, fixed nPostSweeps): the
       clean reference for E, which differs from C also by the controller.
       Same scope as E.
    E variants (amendment C7): configuration E plus one change each,
    compared with C (and E) in the B10 table. Scope in CONFIG_SCOPE:
       E-rcScalar   coupled.rhieChow.tensorial no (scalar D, C2). T2, T4b, T5
       E-nonOrth60  remediation.static.nonOrthThreshold 60 (C1). T4a, T4b, T5
       E-nonOrth65  remediation.static.nonOrthThreshold 65 (C1 value; the
                    default stays 85, D-047/D-051). T4a, T4b, T5
       E-algPair    blockGAMG agglomerator algebraicPair (C4; the templates'
                    agglomerationWeights combined is replaced by pressure,
                    which algebraicPair means). T4a, T4b, T5
       E-sfd        coupled.sfd.enabled yes (C3, T3 only). T3-SST, T3-GEKO
       E-eta07      etaMax 0.7 with minIter 2 ("solve loosely, iterate
                    often"). T2, T4b, T5
    T4a is in the scope of the heavy-case variants because every change is
    tried on T4a before T4b (user rule, 2026-09-22).
    F1   simpleFoam in single precision (amendment D11): the settings of B
         (SIMPLEC), SP build. Compared with B.
    F2   coupledFoam in single precision: the settings of C, SP build.
         Compared with C.
         Scope T1, T3-SST, T4a (D-063: SP vs DP on a few cases only), one
         repeat by default. F1/F2 run from a shell with the SP OpenFOAM
         environment and the SP private install (they are refused in a DP
         shell, and the DP configurations in an SP shell); the case is
         prepared by tests/cflib/precision.py (mesh in DP, origin shift,
         checkMesh gate; SP-geometry-fail -> recorded and skipped). Records
         carry the D11 fields (bench/SCHEMA_precision.md).
Configurations with a scope run only on the cases of their scope unless
--no-scope is given.

Timing. Every rank of every application started through the case's runApp
hook runs under bench/rank_wrapper.sh (/usr/bin/time -v). For coupledFoam
the Allrun script runs potentialFoam first (coupled.potentialInit yes); its
wall and CPU time belong to the cost of the coupledFoam solution and are
added to the totals (wallSeconds, cpuHours) and, as a fixed offset, to the
time to convergence. Both parts are recorded separately. Any other
pre-processing application listed in PRE_APPS is treated the same way,
whichever solver it precedes.

Time to convergence = pre-processing + solver total x progress fraction,
where the progress fraction is the solver's own elapsed wall time at the
convergence iteration divided by its total (simpleFoam ClockTime,
coupledFoam tWall). The same wall-clock fraction scales the solver's CPU
time (assumption: constant CPU/wall ratio during the run).

Usage (from the repository root, one OpenFOAM environment sourced):
    bench/run_bench.py --cases T1,T2 --configs A,B,C,D --repeats 3 --np 1
    bench/run_bench.py --list
    bench/run_bench.py --summary-only
Results: results/bench/<case>_<config>_<run>.json, results/bench/summary.csv,
results/bench/summary.json (table and the B10 acceptance evaluation).
Results whose configuration hash differs from the current definition are
stale: they are excluded from the summary (and listed), and a new run moves
them to results/bench/stale/.

Machine sharing: the harness refuses to start while other solver jobs run
(timings would be meaningless) unless --allow-busy is given; the machine
state before each run is stored in the JSON.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import shutil
import statistics
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "tests"))

from cflib import case as cfcase          # noqa: E402
from cflib import env as cfenv            # noqa: E402
from cflib import logs, post, precision, results   # noqa: E402

WRAPPER = REPO / "bench" / "rank_wrapper.sh"

# Ranks of the heavy benchmark cases (user cap: 10 cores, D-031)
HEAVY_NP = int(os.environ.get("CF_HEAVY_NP", "10"))

# Case table: template, extra Allrun args, monitor, iteration budgets, ranks
CASES = {
    "T1": {"template": "T1_pitzDaily", "args": [], "monitor": "dp",
           "iters": {"simpleFoam": 4000, "coupledFoam": 1000}, "np": 1},
    "T2": {"template": "T2_backwardFacingStep2D", "args": [], "monitor": "dp",
           "iters": {"simpleFoam": 6000, "coupledFoam": 1500}, "np": 1},
    "T3-SST": {"template": "T3_airFoil2D", "args": ["-turbulence", "kOmegaSST"],
               "monitor": "forces",
               "iters": {"simpleFoam": 6000, "coupledFoam": 2000}, "np": 1},
    "T3-GEKO": {"template": "T3_airFoil2D", "args": ["-turbulence", "GEKO"],
                "monitor": "forces",
                "iters": {"simpleFoam": 6000, "coupledFoam": 2000}, "np": 1},
    "T4a": {"template": "T4_motorBike", "args": ["-mesh", "a"],
            "monitor": "forces",
            "iters": {"simpleFoam": 3000, "coupledFoam": 800}, "np": HEAVY_NP,
            "oscillatory": True},
    "T4b": {"template": "T4_motorBike", "args": ["-mesh", "b"],
            "monitor": "forces",
            "iters": {"simpleFoam": 4000, "coupledFoam": 800}, "np": HEAVY_NP,
            "oscillatory": True},
    # CF_T5_MESH=coarse: development mesh as in tests/test_T5_ahmed.py (D-059)
    "T5": {"template": "T5_ahmed",
           "args": ["-mesh", os.environ.get("CF_T5_MESH", "fine")],
           "monitor": "forces",
           "iters": {"simpleFoam": 2000, "coupledFoam": 800}, "np": HEAVY_NP,  # user: 2000 for now (D-042 add. 3)
           "oscillatory": True},
}

# Amendment C7: variants of configuration E, one change each (foamDictionary
# sets on top of E; missing sub-dictionaries are created by
# cflib.case.set_entry)
E_VARIANTS = {
    "E-rcScalar": {"coupled.rhieChow.tensorial": "no"},
    "E-nonOrth60": {"coupled.remediation.static.nonOrthThreshold": 60},
    "E-nonOrth65": {"coupled.remediation.static.nonOrthThreshold": 65},
    "E-algPair": {"solvers.coupled.blockGAMG.agglomerator": "algebraicPair",
                  "solvers.coupled.blockGAMG.agglomerationWeights":
                  "pressure"},
    "E-sfd": {"coupled.sfd.enabled": "yes"},
    "E-eta07": {"solvers.coupled.etaMax": 0.7,
                "solvers.coupled.minIter": 2},
}

CONFIGS = (("A", "B", "C", "D", "E", "F", "G", "H") + tuple(E_VARIANTS)
           + ("F1", "F2"))
NATIVE_CONFIGS = frozenset({"A", "B", "F1"})

# Amendment D11: single-precision configurations and their DP counterpart
# (same settings, DP build); Cd_rel_to_DP etc. are relative to it
SP_BASE = {"F1": "B", "F2": "C"}
SP_CONFIGS = frozenset(SP_BASE)
# Repeats when --repeats is not given (D-063: SP on a few cases, 1 repeat)
DEFAULT_REPEATS = 3
CONFIG_REPEATS = {"F1": 1, "F2": 1}

# B10 case scoping (default when E, F, G, H are requested)
CONFIG_SCOPE = {
    "E": frozenset({"T2", "T4b", "T5"}),
    "F": frozenset({"T1", "T3-SST"}),
    "G": frozenset({"T1", "T3-SST"}),
    "H": frozenset({"T2", "T4b", "T5"}),
    # C7 E variants (C1/C4 name T4b/T5; C3 names T3)
    "E-rcScalar": frozenset({"T2", "T4b", "T5"}),
    "E-nonOrth60": frozenset({"T4a", "T4b", "T5"}),
    "E-nonOrth65": frozenset({"T4a", "T4b", "T5"}),
    "E-algPair": frozenset({"T4a", "T4b", "T5"}),
    "E-sfd": frozenset({"T3-SST", "T3-GEKO"}),
    "E-eta07": frozenset({"T2", "T4b", "T5"}),
    # D11 single precision, scoped per D-063
    "F1": frozenset({"T1", "T3-SST", "T4a"}),
    "F2": frozenset({"T1", "T3-SST", "T4a"}),
}

# Per-case overrides of a configuration (foamDictionary sets).
# A must be the tutorial settings: the airFoil2D tutorial is plain SIMPLE
# with p 0.3 (relaxation.simpleFoam of T3 has the tutorial factors), while
# the T3 case fvSolution has consistent yes for the spec-13 reference.
CASE_OVERRIDES = {
    ("T3-SST", "A"): {"system/fvSolution": {"SIMPLE.consistent": "no"}},
    ("T3-GEKO", "A"): {"system/fvSolution": {"SIMPLE.consistent": "no"}},
}

# Cases whose tutorial configuration (A) runs potentialFoam before the
# native solver (motorBike tutorial Allrun); its cost is counted like
# coupledFoam's (PRE_APPS)
NATIVE_POTENTIAL_CASES = ("T4a", "T4b", "T5")

# Applications that the Allrun scripts run under the rank wrapper before the
# solver (coupledFoam: potentialFoam with coupled.potentialInit yes)
PRE_APPS = ("potentialFoam",)

WINDOW = 100
TOL = 0.002

# Oscillatory wake cases ("oscillatory": True in CASES; D-042 and its
# addendum of 2026-09-22): criterion 12.3(ii) is unsatisfiable there,
# convergence := stationary window mean
STAT_CRITERION = "stationaryMean (D-042 addendum)"
STAT_WINDOW_MIN = 300       # W = max(300, n // STAT_WINDOW_DIV), capped at n (D-042 add. 2)
STAT_WINDOW_DIV = 2
STAT_REL = 0.01             # half-window means differ <= max(1 % |m|,
STAT_ABS = 0.005            #                                  0.005)
STAT_STEP = 50              # scan step of iters_to_stationary
STAT_QUANTITIES = ("Cd", "Cl")
# comparison of the window means with the reference: (relative, absolute)
OSC_TOL = {"Cd": (0.02, 0.002), "Cl": (0.02, 0.01)}

# Mean-field delta comparison of the oscillatory cases (D-042 addendum,
# proposed pass thresholds): fieldAverage function object in the case
# controlDict (entry functions.<FIELD_AVERAGE_FO>), averaged over the same
# window W as the forces; applications/utilities/coupledFieldCompare writes
# the delta fields and <case>/fieldCompare.json. Pass: volume-weighted RMS
# of |dUMean|/U_inf and of |dpMean|/p_ref (p_ref = 0.5 U_inf^2) at most
# FIELD_TOL each.
FIELD_AVERAGE_FO = "fieldAverage"
FIELD_TOL = {"volRmsMagUDeltaRel": 0.02, "volRmsPDeltaRel": 0.02}
FIELD_COMPARE_APP = "coupledFieldCompare"
FIELD_COMPARE_JSON = "fieldCompare.json"

# B10 acceptance (F): adaptive (C) not more than 5 % slower than fixed
# relTol (F), and identical Cd / dp to 1e-4 (relative)
B10_MAX_SLOWDOWN = 0.05
B10_MONITOR_TOL = 1e-4

# Bump when the harness changes how a run is set up or evaluated
HARNESS_VERSION = 3

RELAX_B = """// Benchmark configuration B: SIMPLEC, p 1.0 / U 0.9 / k,omega 0.9
relaxationFactors
{
    fields
    {
        p               1;
    }
    equations
    {
        U               0.9;
        ".*"            0.9;
    }
}
"""


def solver_of(cfg: str) -> str:
    if cfg not in CONFIGS:
        raise ValueError(f"unknown configuration {cfg!r} "
                         f"(known: {','.join(CONFIGS)})")
    return "simpleFoam" if cfg in NATIVE_CONFIGS else "coupledFoam"


def in_scope(name: str, cfg: str) -> bool:
    scope = CONFIG_SCOPE.get(cfg)
    return scope is None or name in scope


def config_sets(cfg: str, spec: dict, name: str | None = None
                ) -> tuple[str, dict, dict]:
    """(solver, foamDictionary sets, files to write) of a configuration,
    including the per-case overrides of CASE_OVERRIDES for case `name`.
    F1/F2 (D11) have the settings of their DP counterpart SP_BASE."""
    solver = solver_of(cfg)
    cfg = SP_BASE.get(cfg, cfg)
    n = spec["iters"][solver]
    sets = {"system/controlDict": {"endTime": n, "writeInterval": n}}
    files = {}
    if solver == "simpleFoam":
        sets["system/fvSolution"] = {"SIMPLE.residualControl": "{}"}
        if cfg == "B":
            sets["system/fvSolution"]["SIMPLE.consistent"] = "yes"
            files["system/relaxation.simpleFoam"] = RELAX_B
    else:
        fv = {
            "coupled.maxIter": n,
            # internal stop disabled: the harness decides (spec 14)
            "coupled.convergence.residualTol": 0,
            "coupled.convergence.forceCoeffsWindow": 10 * n,
        }
        if cfg == "D":
            fv["solvers.coupled.preconditioner"] = "blockDiagonal"
            # the K-cycle needs FGMRES; blockDiagonal is a fixed operator
        if cfg == "E" or cfg in E_VARIANTS:
            # fixed V-cycle: autoTune would promote V -> F -> W (6.3.5)
            fv["solvers.coupled.blockGAMG.cycleType"] = "V"
            fv["solvers.coupled.blockGAMG.autoTune"] = "no"
        fv.update(E_VARIANTS.get(cfg, {}))
        if cfg == "F":
            fv["solvers.coupled.adaptiveRelTol"] = "no"
        if cfg == "G":
            fv["coupled.anderson.enabled"] = "yes"
        if cfg == "H":
            # fixed K-cycle reference for E
            fv["solvers.coupled.blockGAMG.autoTune"] = "no"
        sets["system/fvSolution"] = fv
    for fname, entries in CASE_OVERRIDES.get((name, cfg), {}).items():
        sets.setdefault(fname, {}).update(entries)
    return solver, sets, files


def config_hash(name: str, cfg: str) -> str:
    """Hash of everything that defines a run of (case, configuration)."""
    spec = CASES[name]
    solver, sets, files = config_sets(cfg, spec, name)
    d = {"v": HARNESS_VERSION, "case": name, "cfg": cfg,
         "template": spec["template"], "args": spec["args"],
         "monitor": spec["monitor"], "solver": solver,
         "sets": sets, "files": files,
         "window": WINDOW, "tol": TOL}
    if cfg in SP_CONFIGS:
        # only the SP configurations: the DP hashes stay unchanged
        d["precision"] = "sp"
        d["spHarness"] = {"tolSpVsDp": precision.TOL_SP_VS_DP}
    if spec.get("oscillatory"):
        # D-042 evaluation: only the wake cases' records become stale
        d["criterion"] = {"name": STAT_CRITERION, "wmin": STAT_WINDOW_MIN,
                          "wdiv": STAT_WINDOW_DIV,
                          "rel": STAT_REL, "abs": STAT_ABS,
                          "step": STAT_STEP}
        d["fieldAverageStart"] = field_average_start(spec["iters"][solver])
    blob = json.dumps(d, sort_keys=True, default=str)
    return hashlib.sha256(blob.encode()).hexdigest()[:12]


def is_current(rec: dict) -> bool:
    """True if a result record was produced with the current definition of
    its (case, configuration)."""
    try:
        return rec.get("configHash") == config_hash(rec["case"], rec["config"])
    except (KeyError, ValueError):
        return False


def monitor_history(case: Path, kind: str) -> dict[str, list[float]]:
    if kind == "dp":
        pin = post.surface_value(case, "inletP")
        pout = post.surface_value(case, "outletP")
        key = [k for k in pin if k.startswith("areaAverage")][0]
        n = min(len(pin[key]), len(pout[key]))
        return {"dp": (pin[key][:n] - pout[key][:n]).tolist()}
    fc = post.force_coeffs(case)
    return {"Cd": fc["Cd"].tolist(), "Cl": fc["Cl"].tolist()}


def iters_to_conv(hist: dict[str, list[float]]) -> int | None:
    """First iteration at which every monitored quantity satisfies the
    window criterion (spec 12.3 ii) at the same time."""
    n = min(len(h) for h in hist.values())
    for i in range(WINDOW, n + 1):
        if all(_window_ok(h[:i]) for h in hist.values()):
            return i
    return None


def _window_ok(h: list[float]) -> bool:
    if len(h) < WINDOW:
        return False
    w = h[-WINDOW:]
    mean = sum(w) / WINDOW
    return max(w) - min(w) <= TOL * abs(mean)


# --------------------------------------------------------------------------- #
# oscillatory wake cases (D-042): stationary window mean instead of 12.3(ii)
# --------------------------------------------------------------------------- #

def is_oscillatory(name: str) -> bool:
    """True for the cases with a physically oscillating wake (D-042)."""
    return bool(CASES.get(name, {}).get("oscillatory", False))


def stat_window(n: int) -> int:
    """Averaging window of a run of n iterations: max(300, n//2), capped at
    n (D-042 addenda; the minimum was 1000, which only ever bound the short
    coupledFoam runs - the simpleFoam references have n//2 >= 1500)."""
    return max(0, min(max(STAT_WINDOW_MIN, n // STAT_WINDOW_DIV), n))


def field_average_start(n: int) -> int:
    """timeStart of the fieldAverage function object for a run of n
    iterations (deltaT 1, time = iteration): n - W + 1, so the mean fields
    cover exactly the force window, the last W iterations (the timeControl
    of a function object is active from time >= timeStart - 0.5 deltaT)."""
    return n - stat_window(n) + 1


def foam_dictionary(case: Path, fname: str, args: list[str],
                    capture: bool = False) -> str:
    """foamDictionary -disableFunctionEntries <args> <fname> in `case`.
    ALWAYS with -disableFunctionEntries: without it foamDictionary expands
    and drops #include/#sinclude directives when it rewrites a file."""
    out = subprocess.run(
        ["foamDictionary", "-disableFunctionEntries"] + list(args) + [fname],
        cwd=case, check=True, capture_output=True, text=True)
    return out.stdout if capture else ""


def _function_objects(case: Path) -> list[str]:
    return foam_dictionary(case, "system/controlDict",
                           ["-entry", "functions", "-keywords"],
                           capture=True).split()


def set_field_average_start(case: Path, time_start: int | float) -> None:
    """Set functions.<FIELD_AVERAGE_FO>.timeStart in <case>/system/controlDict
    (the entry must exist: T4/T5 templates, D-042 addendum)."""
    keys = _function_objects(case)
    if FIELD_AVERAGE_FO not in keys:
        raise KeyError(f"{case}/system/controlDict has no functions."
                       f"{FIELD_AVERAGE_FO} (functions: {keys})")
    foam_dictionary(case, "system/controlDict",
                    ["-entry", f"functions.{FIELD_AVERAGE_FO}.timeStart",
                     "-set", f"{time_start:g}"])


# The fieldAverage entry of the T4/T5 templates (for cases prepared before
# it existed: the cached references, continue_reference_with_average)
FIELD_AVERAGE_DICT = (
    "{ type fieldAverage; libs (fieldFunctionObjects); timeStart %s; "
    "writeControl writeTime; restartOnRestart %s; restartOnOutput false; "
    "periodicRestart false; log false; "
    "fields ( U { mean on; prime2Mean on; base iteration; } "
    "p { mean on; prime2Mean off; base iteration; } ); }")


def ensure_field_average(case: Path, time_start: int | float,
                         restart_on_restart: bool = False) -> None:
    """Add the fieldAverage function object to <case>/system/controlDict if
    it is missing, and set its timeStart (and restartOnRestart: a
    continuation must not resume an earlier average)."""
    if FIELD_AVERAGE_FO not in _function_objects(case):
        foam_dictionary(case, "system/controlDict",
                        ["-entry", f"functions.{FIELD_AVERAGE_FO}", "-set",
                         FIELD_AVERAGE_DICT % (f"{time_start:g}", "false")])
    set_field_average_start(case, time_start)
    foam_dictionary(case, "system/controlDict",
                    ["-entry",
                     f"functions.{FIELD_AVERAGE_FO}.restartOnRestart",
                     "-set", "true" if restart_on_restart else "false"])


def mean_field_times(case: Path) -> list[str]:
    """Time directories (serial or processor0) that hold UMean and pMean."""
    out = []
    for base in (case, case / "processor0"):
        if not base.is_dir():
            continue
        for d in base.iterdir():
            try:
                float(d.name)
            except ValueError:
                continue
            if (d / "UMean").exists() and (d / "pMean").exists():
                out.append(d.name)
    return sorted(set(out), key=float)


def field_delta_compare(case: Path, ref_case: Path, nprocs: int = 1,
                        time: str | None = None,
                        ref_time: str | None = None,
                        log_name: str = "log.coupledFieldCompare") -> dict:
    """Mean-field delta comparison (D-042 addendum): run
    applications/utilities/coupledFieldCompare in `case` against `ref_case`
    (serial, or -parallel with nprocs ranks on identically decomposed
    cases) and return <case>/fieldCompare.json plus rc, log and the checks
    of field_checks. U_inf is the magUInf of the case's forceCoeffs
    function object, p_ref = 0.5 U_inf^2 (kinematic pressure)."""
    uinf = float(foam_dictionary(
        case, "system/controlDict",
        ["-entry", "functions.forceCoeffs.magUInf", "-value"],
        capture=True).split()[0])
    pref = 0.5 * uinf ** 2
    cmd = [FIELD_COMPARE_APP, "-reference", str(Path(ref_case).resolve()),
           "-Uinf", repr(uinf), "-pref", repr(pref)]
    if time is not None:
        cmd += ["-time", str(time)]
    if ref_time is not None:
        cmd += ["-referenceTime", str(ref_time)]
    if nprocs > 1:
        cmd = cfenv.mpirun_prefix(nprocs) + cmd + ["-parallel"]
    out_json = case / FIELD_COMPARE_JSON
    out_json.unlink(missing_ok=True)
    rc = cfenv.run(cmd, cwd=case, log=case / log_name)
    rec: dict = {"rc": rc, "log": str(case / log_name), "command": cmd,
                 "Uinf": uinf, "pref": pref, "nProcs": nprocs}
    if rc == 0 and out_json.exists():
        rec.update(json.loads(out_json.read_text()))
    rec.update(field_checks(rec))
    return rec


def field_checks(metrics: dict) -> dict:
    """Proposed pass criteria of the mean-field comparison (D-042
    addendum): volume RMS |dUMean|/U_inf <= 0.02 and volume RMS
    |dpMean|/p_ref <= 0.02."""
    checks = {}
    for key, tol in FIELD_TOL.items():
        v = metrics.get(key)
        checks[key] = v is not None and v <= tol
    return {"fieldTol": dict(FIELD_TOL), "fieldChecks": checks,
            "fieldPass": metrics.get("rc") == 0 and all(checks.values())}


def stationary_mean(hist: dict[str, list[float]], quantity: str,
                    n: int | None = None, w: int | None = None) -> dict:
    """Window statistics of one monitored quantity (D-042).

    n is the number of iterations run (default: the shortest history in
    `hist`, i.e. all quantities are evaluated over the same iterations);
    the window is the last W samples up to n, W = stat_window(n) unless
    given (iters_to_stationary slides the run's window W over the history).
    Returns W, mean, std (population), the half-window means mean1 (first
    W/2) and mean2 (last W/2), the drift |mean1 - mean2|, its tolerance
    max(STAT_REL |mean|, STAT_ABS) = max(1 % |mean|, 0.005) and `stationary`
    (drift <= tolerance). For an
    odd W the middle sample belongs to neither half; a window that does
    not fit (W > n) or has fewer than 2 samples is never stationary."""
    if n is None:
        n = min(len(h) for h in hist.values())
    h = hist[quantity][:n]
    n = len(h)
    if w is None:
        w = stat_window(n)
    half = w // 2
    out: dict = {"n": n, "W": w, "mean": None, "std": None, "mean1": None,
                 "mean2": None, "drift": None, "driftTol": None,
                 "stationary": False}
    if half < 1 or w > n:
        return out
    win = h[n - w:]
    m = sum(win) / w
    m1 = sum(win[:half]) / half
    m2 = sum(win[w - half:]) / half
    tol = max(STAT_REL * abs(m), STAT_ABS)
    out.update({
        "mean": m,
        "std": (sum((x - m) ** 2 for x in win) / w) ** 0.5,
        "mean1": m1, "mean2": m2, "drift": abs(m1 - m2), "driftTol": tol,
        "stationary": abs(m1 - m2) <= tol,
    })
    return out


def is_stationary(hist: dict[str, list[float]], n: int | None = None,
                  quantities: tuple[str, ...] = STAT_QUANTITIES,
                  w: int | None = None) -> bool:
    """Both Cd and Cl stationary in the window (W = w, default
    stat_window(n)) ending at iteration n."""
    return all(stationary_mean(hist, q, n, w)["stationary"]
               for q in quantities)


def iters_to_stationary(hist: dict[str, list[float]],
                        quantities: tuple[str, ...] = STAT_QUANTITIES
                        ) -> int | None:
    """Smallest N <= n (iterations run) such that the window ending at N is
    stationary for all quantities; None if never (D-042).

    "The window" is the run's window W = stat_window(n), slid along the
    history: N is scanned in steps of STAT_STEP (50) from the first
    multiple of 50 >= W, and finally N = n itself, so a stationary final
    window always yields a value. (Re-deriving W from N instead would
    shrink the window early in the run: under the original D-042 rule, on
    the T4a reference, a 75-iteration window at N = 150 was stationary.)"""
    if not hist:
        return None
    n = min(len(h) for h in hist.values())
    w = stat_window(n)
    if w < 2:
        return None
    first = -(-w // STAT_STEP) * STAT_STEP
    cand = list(range(first, n + 1, STAT_STEP))
    if not cand or cand[-1] != n:
        cand.append(n)
    for N in cand:
        if is_stationary(hist, N, quantities, w):
            return N
    return None


def stationary_eval(hist: dict[str, list[float]],
                    quantities: tuple[str, ...] = STAT_QUANTITIES) -> dict:
    """Record of the D-042 evaluation of a force history: W, means, stds,
    half-window drifts per quantity, `stationary` (final window, all
    quantities) and iters_to_stationary."""
    if not hist or not all(q in hist for q in quantities):
        return {"criterion": STAT_CRITERION, "stationary": False,
                "iters_to_stationary": None}
    n = min(len(h) for h in hist.values())
    out: dict = {"criterion": STAT_CRITERION, "iterations_run": n,
                 "W": stat_window(n)}
    for q in quantities:
        s = stationary_mean(hist, q, n)
        out.update({f"{q}_mean": s["mean"], f"{q}_std": s["std"],
                    f"{q}_mean1": s["mean1"], f"{q}_mean2": s["mean2"],
                    f"{q}_drift": s["drift"], f"{q}_driftTol": s["driftTol"],
                    f"{q}_stationary": s["stationary"]})
    out["stationary"] = all(out[f"{q}_stationary"] for q in quantities)
    out["iters_to_stationary"] = iters_to_stationary(hist, quantities)
    return out


def mean_comparison(rec: dict, ref: dict) -> dict:
    """D-042 comparison of window means (keys Cd_mean, Cl_mean of
    stationary_eval): PASS iff |Cd - Cd_ref| <= max(2 % |Cd_ref|, 0.002)
    and |Cl - Cl_ref| <= max(2 % |Cl_ref|, 0.01)."""
    out: dict = {"criterion": STAT_CRITERION}
    ok = True
    for q, (rel, ab) in OSC_TOL.items():
        a, b = rec.get(f"{q}_mean"), ref.get(f"{q}_mean")
        if a is None or b is None:
            out[f"{q}_pass"] = False
            ok = False
            continue
        tol = max(rel * abs(b), ab)
        out.update({f"{q}_absDiff": abs(a - b),
                    f"{q}_relDiff": abs(a - b) / max(abs(b), 1e-12),
                    f"{q}_tol": tol, f"{q}_pass": abs(a - b) <= tol})
        ok = ok and out[f"{q}_pass"]
    out["pass"] = ok
    return out


def _parse_time_reports(reps: list[Path]) -> dict:
    """wall (max over ranks), CPU (user + sys summed) and peak RSS of a set
    of /usr/bin/time -v reports of one application."""
    wall, cpu, rss, incomplete = [], [], [], []
    for r in reps:
        t = r.read_text(errors="replace")
        mu = re.search(r"User time \(seconds\): ([\d.]+)", t)
        ms = re.search(r"System time \(seconds\): ([\d.]+)", t)
        mm = re.search(r"Maximum resident set size \(kbytes\): (\d+)", t)
        me = re.search(r"Elapsed \(wall clock\) time .*: ([\d:.]+)", t)
        if not (mu and ms and mm and me):
            # a rank killed by MPI_Abort (e.g. a B4 abort) may leave an
            # incomplete report: skip it, record it, never crash the caller
            incomplete.append(r.name)
            continue
        parts = [float(p) for p in me.group(1).split(":")]
        w = sum(p * 60 ** i for i, p in enumerate(reversed(parts)))
        wall.append(w)
        cpu.append(float(mu.group(1)) + float(ms.group(1)))
        rss.append(int(mm.group(1)))
    if not wall:
        return {"ranks": len(reps), "incompleteReports": incomplete}
    return {
        "incompleteReports": incomplete,
        "ranks": len(reps),
        "wallSeconds": max(wall),
        "cpuSeconds": sum(cpu),
        "cpuHours": sum(cpu) / 3600.0,
        "peakRSS_GB_sum": sum(rss) / 1024**2,
        "peakRSS_GB_max_rank": max(rss) / 1024**2,
    }


def rank_times(case: Path, solver: str, timing_dir: str = "timing") -> dict:
    """Wall time, CPU time and peak RSS from the /usr/bin/time -v reports in
    <case>/<timing_dir> (default timing; a reference continuation keeps its
    own reports in timing_continuation).

    The totals (wallSeconds, cpuSeconds, cpuHours) include the
    pre-processing applications of PRE_APPS (potentialFoam), which run
    sequentially before the solver in the same case: wall times add, CPU
    times add. The parts are recorded separately:
    solverWallSeconds / solverCpuHours and <app>WallSeconds / <app>CpuHours
    (e.g. potentialFoamWallSeconds). Peak RSS is the larger of the solver's
    and the pre-processing's (they never run at the same time); the
    solver's own values are kept as solverPeakRSS_GB_*.
    Returns {} if the solver left no report (run failed before start)."""
    tdir = case / timing_dir
    reps = sorted(tdir.glob(f"{solver}.rank*.time"))
    if not reps:
        return {}
    sol = _parse_time_reports(reps)
    out = {
        "ranks": sol["ranks"],
        "solverWallSeconds": sol["wallSeconds"],
        "solverCpuSeconds": sol["cpuSeconds"],
        "solverCpuHours": sol["cpuHours"],
        "solverPeakRSS_GB_sum": sol["peakRSS_GB_sum"],
        "solverPeakRSS_GB_max_rank": sol["peakRSS_GB_max_rank"],
        "preApps": [],
        "preWallSeconds": 0.0,
        "preCpuSeconds": 0.0,
    }
    rss_sum = sol["peakRSS_GB_sum"]
    rss_max = sol["peakRSS_GB_max_rank"]
    for app in PRE_APPS:
        if app == solver:
            continue
        preps = sorted(tdir.glob(f"{app}.rank*.time"))
        if not preps:
            continue
        p = _parse_time_reports(preps)
        out["preApps"].append(app)
        out[f"{app}Ranks"] = p["ranks"]
        out[f"{app}WallSeconds"] = p["wallSeconds"]
        out[f"{app}CpuSeconds"] = p["cpuSeconds"]
        out[f"{app}CpuHours"] = p["cpuHours"]
        out[f"{app}PeakRSS_GB_sum"] = p["peakRSS_GB_sum"]
        out["preWallSeconds"] += p["wallSeconds"]
        out["preCpuSeconds"] += p["cpuSeconds"]
        rss_sum = max(rss_sum, p["peakRSS_GB_sum"])
        rss_max = max(rss_max, p["peakRSS_GB_max_rank"])
    out["preCpuHours"] = out["preCpuSeconds"] / 3600.0
    out["wallSeconds"] = sol["wallSeconds"] + out["preWallSeconds"]
    out["cpuSeconds"] = sol["cpuSeconds"] + out["preCpuSeconds"]
    out["cpuHours"] = out["cpuSeconds"] / 3600.0
    out["peakRSS_GB_sum"] = rss_sum
    out["peakRSS_GB_max_rank"] = rss_max
    return out


def progress_fraction(case: Path, solver: str, it: int) -> float | None:
    """Fraction of the solver's own run time spent up to iteration `it`
    (wall clock). Pre-processing is not part of it (see to_convergence).

    CPU time to convergence is scaled with this wall-clock fraction too,
    i.e. a constant CPU/wall ratio over the run is assumed."""
    if it is None or it < 1:
        return None
    log = case / f"log.{solver}"
    if solver == "coupledFoam":
        rows = [r for r in logs.parse_cf(log) if "tWall" in r]
        if not rows or rows[-1]["tWall"] <= 0:
            return None
        return rows[min(it, len(rows)) - 1]["tWall"] / rows[-1]["tWall"]
    nat = logs.parse_native(log)
    clk = nat["clock"]
    if not clk or clk[-1] <= 0:
        return None
    return clk[min(it, len(clk)) - 1] / clk[-1]


def to_convergence(rt: dict, frac: float | None) -> dict:
    """Wall seconds and CPU-hours to convergence: pre-processing (fixed
    offset, not scaled) + solver totals x progress fraction."""
    if not rt or frac is None:
        return {}
    return {
        "wall_to_conv_s": rt.get("preWallSeconds", 0.0)
        + rt["solverWallSeconds"] * frac,
        "cpu_to_conv_h": rt.get("preCpuHours", 0.0)
        + rt["solverCpuHours"] * frac,
        "solver_wall_to_conv_s": rt["solverWallSeconds"] * frac,
        "solver_cpu_to_conv_h": rt["solverCpuHours"] * frac,
    }


def coupled_breakdown(case: Path, it: int) -> dict:
    rows = logs.parse_cf(case / "log.coupledFoam")[:it]
    return {
        "t_assembly": sum(r.get("tAsm", 0) for r in rows),
        "t_linsolve": sum(r.get("tSolve", 0) for r in rows),
        "t_turb": sum(r.get("tTurb", 0) for r in rows),
        "linIters_total": sum(r.get("linIters", 0) for r in rows),
        "eta_median": sorted(r.get("eta", 0) for r in rows)[len(rows) // 2]
        if rows else None,
    }


def _stale_move(tag: str) -> None:
    src = results.RESULTS / "bench" / f"{tag}.json"
    if not src.exists():
        return
    old = json.loads(src.read_text())
    dst = results.RESULTS / "bench" / "stale" / \
        f"{tag}.{old.get('configHash', 'nohash')}.json"
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.move(str(src), str(dst))
    print(f"stale result {tag} moved to {dst}")


# --------------------------------------------------------------------------- #
# amendment D11: single-precision record fields (bench/SCHEMA_precision.md)
# --------------------------------------------------------------------------- #

# Quantities compared with the DP counterpart: value at convergence and at
# the end of the fixed budget
D11_METRICS = ("Cd", "Cl", "dp", "Cd_final", "Cl_final", "dp_final")


def dp_counterpart(name: str, cfg: str) -> list[dict]:
    """Current records of the DP configuration of an SP configuration
    (F1 -> B, F2 -> C) on case `name`."""
    base = SP_BASE.get(cfg)
    if base is None:
        return []
    recs, _ = load_current([name], [base])
    return [r for r in recs if not r.get("skipped")]


def precision_fields(name: str, cfg: str, case: Path, rec: dict) -> dict:
    """D11 fields of a benchmark record, in place. DP configurations get
    precision "dp" and None for the SP-only fields."""
    sp = cfg in SP_CONFIGS
    rec["precision"] = precision.label() if sp else "dp"
    rec["buildPrecision"] = precision.build_precision()
    if not sp:
        rec.update({"meshShift": None, "checkMeshDiff": None,
                    "nCheckMeshDiff": None,
                    "spGeometry": None, "staticSetSizeDP": None,
                    "staticSetSizeDiff": None, "Cd_rel_to_DP": None})
        return rec
    info = precision.read_info(case)
    rec["meshShift"] = info.get("shift")
    rec["checkMeshDiff"] = info.get("checkMeshDiff")
    rec["nCheckMeshDiff"] = (len(rec["checkMeshDiff"])
                             if rec["checkMeshDiff"] is not None else None)
    rec["spGeometry"] = info.get("status")
    rec["spHarness"] = {k: info.get(k) for k in (
        "bboxCentreDP", "shiftedEntries", "unclassifiedVectors", "reasons",
        "gateOverridden",
        "wallSecondsDPMesh", "wallSecondsPrepare")}
    dps = dp_counterpart(name, cfg)
    rec["dpCounterpart"] = {"config": SP_BASE[cfg],
                            "runs": [r.get("run") for r in dps]}
    sizes = [r["staticSetSize"] for r in dps
             if r.get("staticSetSize") is not None]
    rec["staticSetSizeDP"] = sizes[0] if sizes else None
    s = rec.get("staticSetSize")
    rec["staticSetSizeDiff"] = (s - rec["staticSetSizeDP"]
                                if None not in (s, rec["staticSetSizeDP"])
                                else None)
    for m in D11_METRICS:
        v, ref = rec.get(m), _med(dps, m)
        rec[f"{m}_rel_to_DP"] = ((v - ref) / abs(ref)
                                 if v is not None and ref not in (None, 0)
                                 else None)
        rec[f"{m}_DP"] = ref
    mon = ("dp",) if CASES[name]["monitor"] == "dp" else ("Cd", "Cl")
    rels = [rec.get(f"{m}_rel_to_DP") for m in mon]
    rec["spVsDpTol"] = precision.TOL_SP_VS_DP
    rec["spVsDpPass"] = (None if None in rels else
                         all(abs(r) <= precision.TOL_SP_VS_DP for r in rels))
    for k, kk in (("wall_to_conv_s", "Wall"), ("cpu_to_conv_h", "Cpu")):
        ref, v = _med(dps, k), rec.get(k)
        rec[f"speedup{kk}_DP_over_SP"] = (ref / v if ref and v else None)
    return rec


def run_one(name: str, cfg: str, run: int, nprocs: int, force: bool) -> dict:
    spec = CASES[name]
    tag = f"{name}_{cfg}_{run}"
    chash = config_hash(name, cfg)
    old = results.read("bench", tag)
    if old and not force:
        if old.get("configHash") == chash:
            print(f"skip {tag} (exists)")
            return old
        print(f"{tag}: existing result has configHash "
              f"{old.get('configHash')}, current {chash}: rerun")
    if old:
        _stale_move(tag)

    solver, sets, files = config_sets(cfg, spec, name)
    state = cfenv.machine_state().as_dict()
    case = cfcase.prepare(spec["template"], f"bench_{tag}", sets)
    for rel, text in files.items():
        (case / rel).write_text(text)
    if is_oscillatory(name):
        # mean fields over the force window (D-042 addendum)
        set_field_average_start(case, field_average_start(spec["iters"][solver]))

    try:
        rc = cfcase.allrun(
            case, ["-solver", solver, "-np", str(nprocs)] + spec["args"],
            fpe=False,
            extra_env={"CF_RANK_WRAPPER": str(WRAPPER),
                       "CF_TIMING_DIR": str(case / "timing"),
                       # motorBike tutorial: potentialFoam before simpleFoam
                       "CF_NATIVE_POTENTIAL":
                       "yes" if (cfg == "A" and name in NATIVE_POTENTIAL_CASES)
                       else "no"},
        )
    except precision.SPGeometryFail as err:
        # D5.3: not a solver failure; the SP run of this case is skipped
        rec = {"case": name, "config": cfg, "run": run, "solver": solver,
               "nProcs": nprocs, "rc": None, "machineBefore": state,
               "configHash": chash, "harnessVersion": HARNESS_VERSION,
               "sets": sets, "skipped": True, "error": str(err)}
        precision_fields(name, cfg, case, rec)
        results.write("bench", tag, rec)
        print(f"{tag}: {precision.SP_GEOMETRY_FAIL}, skipped ({err})")
        return rec
    rec = {"case": name, "config": cfg, "run": run, "solver": solver,
           "nProcs": nprocs, "rc": rc, "machineBefore": state,
           "configHash": chash, "harnessVersion": HARNESS_VERSION,
           "sets": sets}
    rt = rank_times(case, solver)
    rec["total"] = rt
    # totals of the whole solution (pre-processing included) and the parts
    for k in ("wallSeconds", "cpuHours", "solverWallSeconds",
              "solverCpuHours", "preWallSeconds", "preCpuHours"):
        rec[k] = rt.get(k)
    for app in rt.get("preApps", []):
        rec[f"{app}WallSeconds"] = rt.get(f"{app}WallSeconds")
        rec[f"{app}CpuHours"] = rt.get(f"{app}CpuHours")
    rec["peakRSS_GB_sum"] = rt.get("peakRSS_GB_sum")
    rec["peakRSS_GB_max_rank"] = rt.get("peakRSS_GB_max_rank")
    if solver == "coupledFoam":
        summ = logs.coupled_summary(case)
        rec["nCells"] = summ.get("nCells")
        rec["gamg"] = {k: summ.get(k) for k in (
            "gamgLevels", "gamgMergeLevels", "gamgCop", "gamgCellsPerLevel")}
        for k, v in logs.gamg_log_stats(case / "log.coupledFoam").items():
            if rec["gamg"].get(k) is None:
                rec["gamg"][k] = v
        # final cycle type and nPostSweeps after autoTune (E must stay V)
        rec["cycleTypeFinal"] = rec["gamg"].get("gamgCycleFinal")
        rec["nPostSweepsFinal"] = rec["gamg"].get("gamgNPostSweepsFinal")
        rec["gamgTuneEvents"] = summ.get("gamgTuneEvents")
        # static remediation set of this build (D5.5, D11 staticSetSizeDiff)
        rows0 = logs.parse_cf(case / "log.coupledFoam")[:1]
        rec["staticSetSize"] = rows0[0].get("nStat") if rows0 else None
        if "andersonApplied" in summ:
            rec["anderson"] = {k: summ.get(k) for k in (
                "andersonApplied", "andersonSkipped", "andersonRejected",
                "andersonFlushed")}
    try:
        hist = monitor_history(case, spec["monitor"])
    except (FileNotFoundError, IndexError) as err:
        rec["error"] = f"monitor: {err}"
        precision_fields(name, cfg, case, rec)
        results.write("bench", tag, rec)
        return rec

    osc = is_oscillatory(name)
    if osc:
        # D-042: convergence := stationary window mean of Cd and Cl; the
        # first stationary window gives the time to convergence
        st = stationary_eval(hist)
        rec.update(st)
        it = st["iters_to_stationary"]
    else:
        rec["criterion"] = f"window {WINDOW}, tol {TOL} (12.3 ii)"
        it = iters_to_conv(hist)
    rec["iters_to_conv"] = it
    rec["iterations_run"] = min(len(h) for h in hist.values())
    # final values at the end of the fixed budget (B10 identity check)
    rec.update({f"{k}_final": v[-1] for k, v in hist.items() if v})
    if it is not None and rt:
        frac = progress_fraction(case, solver, it)
        rec["progressFraction"] = frac
        tc = to_convergence(rt, frac)
        rec.update(tc)
        if tc:
            # per iteration: the solver part only (pre-processing is a
            # one-off cost)
            rec["time_per_iter_s"] = tc["solver_wall_to_conv_s"] / it
            rec["cpu_per_iter_s"] = tc["solver_cpu_to_conv_h"] * 3600.0 / it
        if osc:
            # compared quantities: the final-window means (D-042)
            rec.update({k: rec.get(f"{k}_mean") for k in hist})
        else:
            rec.update({k: v[it - 1] for k, v in hist.items()})
            rec.update({f"{k}_windowMean": sum(v[it - WINDOW:it]) / WINDOW
                        for k, v in hist.items()})
        if solver == "coupledFoam":
            rec.update(coupled_breakdown(case, it))
    rec["forceSource"] = (post.force_source(case)
                          if spec["monitor"] == "forces" else None)
    precision_fields(name, cfg, case, rec)
    results.write("bench", tag, rec)
    print(f"{tag}: iters {it}, wall {rec.get('wall_to_conv_s')}, "
          f"CPU-h {rec.get('cpu_to_conv_h')}")
    return rec


# --------------------------------------------------------------------------- #
# summary and B10 acceptance
# --------------------------------------------------------------------------- #

def _med(g: list[dict], key: str):
    v = [x[key] for x in g if x.get(key) is not None]
    return statistics.median(v) if v else None


def _rng(g: list[dict], key: str, f):
    v = [x[key] for x in g if x.get(key) is not None]
    return f(v) if v else None


def summary_rows(recs: list[dict]) -> list[dict]:
    """Per (case, configuration): medians, min and max of the repeats."""
    groups: dict = {}
    for d in recs:
        groups.setdefault((d["case"], d["config"]), []).append(d)
    table = []
    for (c, cfg), g in sorted(groups.items()):
        table.append({
            "case": c, "config": cfg, "n": len(g),
            "iters_median": _med(g, "iters_to_conv"),
            "wall_median": _med(g, "wall_to_conv_s"),
            "wall_min": _rng(g, "wall_to_conv_s", min),
            "wall_max": _rng(g, "wall_to_conv_s", max),
            "cpuh_median": _med(g, "cpu_to_conv_h"),
            "cpuh_min": _rng(g, "cpu_to_conv_h", min),
            "cpuh_max": _rng(g, "cpu_to_conv_h", max),
            "pre_wall_median": _med(g, "preWallSeconds"),
            "pre_cpuh_median": _med(g, "preCpuHours"),
            "time_per_iter_median": _med(g, "time_per_iter_s"),
            "cpu_per_iter_median": _med(g, "cpu_per_iter_s"),
            "rss_sum_GB": _med(g, "peakRSS_GB_sum"),
            "nCells": _med(g, "nCells"),
            "Cd": _med(g, "Cd"), "Cl": _med(g, "Cl"), "dp": _med(g, "dp"),
            "Cd_final": _med(g, "Cd_final"), "Cl_final": _med(g, "Cl_final"),
            "dp_final": _med(g, "dp_final"),
            "Cd_std": _med(g, "Cd_std"), "Cl_std": _med(g, "Cl_std"),
            "W": _med(g, "W"),
            "criterion": ",".join(sorted({str(x["criterion"]) for x in g
                                          if x.get("criterion")})) or None,
            "cycleTypeFinal": ",".join(sorted({str(x.get("cycleTypeFinal"))
                                               for x in g
                                               if x.get("cycleTypeFinal")}))
            or None,
            "nPostSweepsFinal": _med(g, "nPostSweepsFinal"),
            # D11 (single precision)
            "precision": ",".join(sorted({str(x.get("precision", "dp"))
                                          for x in g})),
            "Cd_rel_to_DP": _med(g, "Cd_rel_to_DP"),
            "Cl_rel_to_DP": _med(g, "Cl_rel_to_DP"),
            "dp_rel_to_DP": _med(g, "dp_rel_to_DP"),
            "staticSetSizeDiff": _med(g, "staticSetSizeDiff"),
            "spGeometry": ",".join(sorted({str(x["spGeometry"]) for x in g
                                           if x.get("spGeometry")})) or None,
            "nCheckMeshDiff": _med(g, "nCheckMeshDiff"),
        })
    by = {(t["case"], t["config"]): t for t in table}
    for t in table:
        b, cc = by.get((t["case"], "B")), by.get((t["case"], "C"))
        if b and cc and b["wall_median"] and cc["wall_median"]:
            t["speedup_wall_B_over_C"] = b["wall_median"] / cc["wall_median"]
        if b and cc and b["cpuh_median"] and cc["cpuh_median"]:
            t["speedup_cpu_B_over_C"] = b["cpuh_median"] / cc["cpuh_median"]
        # D11: DP counterpart over SP (F1 vs B, F2 vs C)
        base = by.get((t["case"], SP_BASE.get(t["config"], "")))
        if base:
            for k, kk in (("wall_median", "wall"), ("cpuh_median", "cpu")):
                if base[k] and t[k]:
                    t[f"speedup_{kk}_DP_over_SP"] = base[k] / t[k]
    return table


def _rel(a, b):
    """a/b - 1, None if not computable."""
    if a is None or not b:
        return None
    return a / b - 1.0


def b10_evaluate(table: list[dict]) -> list[dict]:
    """Amendment B10 comparisons against configuration C, wall and CPU-hours.

    C vs E (V-cycle), C vs H (fixed K; H vs E is the clean cycle
    comparison), C vs G (Anderson): deltas only, B10 gives no threshold.
    C vs F (fixed relTol): pass if adaptive C is at most 5 % slower than F
    in wall time AND in CPU-hours, and the monitored quantity (Cd for force
    cases, dp for T1/T2) at the end of the fixed budget agrees to 1e-4
    (relative). Deltas are X/C - 1 (positive: X slower than C)."""
    by = {(t["case"], t["config"]): t for t in table}
    out = []
    for (case, cfg), t in sorted(by.items()):
        if cfg not in ("E", "F", "G", "H") and cfg not in E_VARIANTS:
            continue
        ref_cfg = "C"
        c = by.get((case, ref_cfg))
        row = {"case": case, "config": cfg, "reference": ref_cfg,
               "wall_X": t["wall_median"], "cpuh_X": t["cpuh_median"],
               "iters_X": t["iters_median"]}
        if c is None:
            row.update({"status": "no reference C", "pass": None})
            out.append(row)
            continue
        row.update({
            "wall_C": c["wall_median"], "cpuh_C": c["cpuh_median"],
            "iters_C": c["iters_median"],
            "dWall_X_vs_C": _rel(t["wall_median"], c["wall_median"]),
            "dCpu_X_vs_C": _rel(t["cpuh_median"], c["cpuh_median"]),
        })
        mon = "dp" if CASES.get(case, {}).get("monitor") == "dp" else "Cd"
        xv, cv = t.get(f"{mon}_final"), c.get(f"{mon}_final")
        row["monitor"] = mon
        row["monitorRelDiff"] = (abs(cv - xv) / abs(xv)
                                 if xv not in (None, 0) and cv is not None
                                 else None)
        if cfg == "F":
            # C slower than F by (C/F - 1)
            sw = _rel(c["wall_median"], t["wall_median"])
            sc = _rel(c["cpuh_median"], t["cpuh_median"])
            row["slowdownWall_C_vs_F"] = sw
            row["slowdownCpu_C_vs_F"] = sc
            ok = (sw is not None and sc is not None
                  and row["monitorRelDiff"] is not None
                  and sw <= B10_MAX_SLOWDOWN and sc <= B10_MAX_SLOWDOWN
                  and row["monitorRelDiff"] <= B10_MONITOR_TOL)
            row["criterion"] = (f"C at most {B10_MAX_SLOWDOWN:.0%} slower than "
                                f"F (wall and CPU-h); {mon} identical to "
                                f"{B10_MONITOR_TOL:g}")
            row["pass"] = ok if None not in (sw, sc, row["monitorRelDiff"]) \
                else None
            row["status"] = ("pass" if row["pass"] else
                             "FAIL" if row["pass"] is False else "incomplete")
        else:
            row["criterion"] = "none (B10: report the delta)"
            row["pass"] = None
            row["status"] = "reported"
        if cfg == "E":
            h = by.get((case, "H"))
            if h:
                row["dWall_E_vs_H"] = _rel(t["wall_median"], h["wall_median"])
                row["dCpu_E_vs_H"] = _rel(t["cpuh_median"], h["cpuh_median"])
            row["cycleTypeFinal"] = t.get("cycleTypeFinal")
        if cfg in E_VARIANTS:
            # C7: the variant against plain E isolates its one change
            e = by.get((case, "E"))
            if e:
                row["dWall_X_vs_E"] = _rel(t["wall_median"], e["wall_median"])
                row["dCpu_X_vs_E"] = _rel(t["cpuh_median"], e["cpuh_median"])
                row["iters_E"] = e["iters_median"]
            row["cycleTypeFinal"] = t.get("cycleTypeFinal")
        out.append(row)
    return out


def load_current(names: list[str] | None = None,
                 cfgs: list[str] | None = None) -> tuple[list[dict], list[str]]:
    """Current benchmark records (optionally filtered) and the stale tags."""
    recs, stale = [], []
    d = results.RESULTS / "bench"
    for f in sorted(d.glob("*_*_*.json")) if d.is_dir() else []:
        try:
            rec = json.loads(f.read_text())
        except json.JSONDecodeError:
            continue
        if "case" not in rec or "config" not in rec:
            continue
        if names and rec["case"] not in names:
            continue
        if cfgs and rec["config"] not in cfgs:
            continue
        if not is_current(rec):
            stale.append(f.stem)
            continue
        recs.append(rec)
    return recs, stale


SUMMARY_FIELDS = [
    "case", "config", "n", "iters_median", "wall_median", "wall_min",
    "wall_max", "cpuh_median", "cpuh_min", "cpuh_max", "pre_wall_median",
    "pre_cpuh_median", "time_per_iter_median", "cpu_per_iter_median",
    "rss_sum_GB", "nCells", "Cd", "Cl", "dp", "Cd_final", "Cl_final",
    "dp_final", "Cd_std", "Cl_std", "W", "criterion", "cycleTypeFinal",
    "nPostSweepsFinal",
    "speedup_wall_B_over_C", "speedup_cpu_B_over_C",
    "precision", "Cd_rel_to_DP", "Cl_rel_to_DP", "dp_rel_to_DP",
    "staticSetSizeDiff", "spGeometry", "nCheckMeshDiff",
    "speedup_wall_DP_over_SP", "speedup_cpu_DP_over_SP"]
B10_FIELDS = [
    "case", "config", "reference", "status", "pass", "criterion",
    "wall_X", "wall_C", "dWall_X_vs_C", "cpuh_X", "cpuh_C", "dCpu_X_vs_C",
    "slowdownWall_C_vs_F", "slowdownCpu_C_vs_F", "monitor", "monitorRelDiff",
    "dWall_E_vs_H", "dCpu_E_vs_H", "dWall_X_vs_E", "dCpu_X_vs_E",
    "cycleTypeFinal", "iters_X", "iters_C", "iters_E"]


def summarise(names: list[str], cfgs: list[str]) -> Path:
    recs, stale = load_current(names, cfgs)
    for s in stale:
        print(f"summary: stale result {s} excluded (configHash mismatch)")
    table = summary_rows(recs)
    b10 = b10_evaluate(table)
    out = results.RESULTS / "bench" / "summary.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=SUMMARY_FIELDS)
        w.writeheader()
        for t in table:
            w.writerow({k: t.get(k) for k in SUMMARY_FIELDS})
    with open(out.with_name("b10_acceptance.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=B10_FIELDS)
        w.writeheader()
        for t in b10:
            w.writerow({k: t.get(k) for k in B10_FIELDS})
    (out.with_name("summary.json")).write_text(json.dumps(
        {"cases": names, "configs": cfgs, "table": table, "b10": b10,
         "stale": stale}, indent=2, default=str) + "\n")
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cases", default="T1,T2,T3-SST,T3-GEKO")
    ap.add_argument("--configs", default=",".join(CONFIGS),
                    help="E, F, G, H and the C7 E variants run only on "
                         "their scope (see --no-scope)")
    ap.add_argument("--no-scope", action="store_true",
                    help="run E/F/G/H/E-* on every requested case")
    ap.add_argument("--repeats", type=int, default=None,
                    help=f"repeats per configuration (default "
                         f"{DEFAULT_REPEATS}; F1, F2: 1)")
    ap.add_argument("--np", type=int, default=0,
                    help="ranks (0: the case default)")
    ap.add_argument("--force", action="store_true", help="rerun existing")
    ap.add_argument("--allow-busy", action="store_true")
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--summary-only", action="store_true")
    a = ap.parse_args()

    if a.list:
        for k, v in CASES.items():
            cf = [c for c in CONFIGS if in_scope(k, c)]
            print(f"{k:8s} {v['template']:26s} np={v['np']} "
                  f"monitor={v['monitor']} configs={','.join(cf)}")
        return 0
    names = [n for n in a.cases.split(",") if n]
    cfgs = [c for c in a.configs.split(",") if c]
    for n in names:
        if n not in CASES:
            ap.error(f"unknown case {n!r} (known: {','.join(CASES)})")
    for c in cfgs:
        if c not in CONFIGS:
            ap.error(f"unknown configuration {c!r} (known: {','.join(CONFIGS)})")
    if not a.summary_only:
        cfenv.foam_env()
        # D11: precision is a build choice - SP configurations from an SP
        # shell only, DP configurations from a DP shell only
        sp_cfgs = [c for c in cfgs if c in SP_CONFIGS]
        if sp_cfgs and len(sp_cfgs) != len(cfgs):
            if "--configs" in " ".join(sys.argv):
                ap.error("run the SP configurations (F1, F2) separately, "
                         "from a shell with the SP environment")
            # default list: keep what the sourced build can run
            sp_env = precision.build_precision() == "SP"
            cfgs = sp_cfgs if sp_env else [c for c in cfgs
                                           if c not in SP_CONFIGS]
            sp_cfgs = sp_cfgs if sp_env else []
        if sp_cfgs:
            os.environ["CF_PRECISION"] = "sp"
        else:
            os.environ["CF_PRECISION"] = "dp"
        try:
            precision.check_environment()
        except RuntimeError as err:
            ap.error(str(err))
        rep_max = a.repeats or max([CONFIG_REPEATS.get(c, DEFAULT_REPEATS)
                                    for c in cfgs] or [DEFAULT_REPEATS])
        for name in names:
            nprocs = a.np or CASES[name]["np"]
            for run in range(1, rep_max + 1):
                for cfg in cfgs:
                    if not a.no_scope and not in_scope(name, cfg):
                        continue
                    if run > (a.repeats
                              or CONFIG_REPEATS.get(cfg, DEFAULT_REPEATS)):
                        continue
                    st = cfenv.machine_state()
                    if st.busy and not a.allow_busy:
                        print("machine busy (other jobs: %d, load %.1f): "
                              "refusing to benchmark; --allow-busy overrides"
                              % (len(st.jobs), st.load))
                        return 3
                    run_one(name, cfg, run, nprocs, a.force)
    print("summary:", summarise(names, cfgs))
    return 0


if __name__ == "__main__":
    sys.exit(main())

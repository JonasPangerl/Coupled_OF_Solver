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
(stationary_mean: the two half-window means of Cd and of Cl differ by at
most max(1 % |mean|, 0.005)); the time to convergence uses the first
stationary window (iters_to_stationary, scanned in steps of 50) in place
of iters_to_conv, and Cd / Cl are the final-window means. The window W is
ONE value per case for both solvers (D-068): CASES[<case>]["statWindow"]
= max(300, coupledFoam budget // 2) = 400 on T4a, T4b and T5 (capped at
the iterations run). The per-run rule of D-042 addendum 2,
W = max(300, n/2) from each run's own budget n, made the simpleFoam window
(1500 / 2000) and therefore its earliest possible convergence point a
function of its budget; it is kept only as an informational sensitivity
value (keys *_perRun in the records). The mean fields (fieldAverage
function object of the T4/T5 controlDict) are averaged over the same
window: set_field_average_start sets its timeStart to n - W + 1.

Configurations (DECISIONS.md D-025, amendment B10, redefined in D-068;
every configuration is a distinct run on every case of its scope, which
tests/test_harness.py checks against the templates):
    A  simpleFoam, the tutorial's solver settings and relaxation factors.
       Per-case overrides where the test case differs from the tutorial:
       T3 (airFoil2D tutorial): plain SIMPLE (consistent no), p 0.3,
       U / turbulence 0.7. T5 has no tutorial: A = the motorBike tutorial
       settings (SIMPLEC, U 0.9, k/omega 0.7), as in T4.
    B  simpleFoam SIMPLEC, consistent yes, relaxation p 1.0 / U 0.9 / k,omega 0.9.
       Not on T1: the pitzDaily tutorial IS this setting (SIMPLEC, p
       unrelaxed, U and turbulence 0.9), so B == A there.
    C  coupledFoam, the case template as it is: since D-043 a V-cycle with
       autoTune off, adaptive relTol (Eisenstat-Walker).
    D  C with preconditioner blockDiagonal (isolates the AMG gain)
    F  C with adaptiveRelTol no (isolates Eisenstat-Walker). B10: T1, T3-SST.
    G  C with anderson enabled. B10: T1, T3-SST.
    H  C with a fixed K-cycle: cycleType K, autoTune no. H vs C is the
       cycle comparison of B10 (K vs V, controller off in both).
    H-tune  C with the K-cycle and the autoTune controller (cycleType K,
       autoTune yes): the pre-D-043 default; H-tune vs H isolates the
       controller.
    E  no longer a configuration: the old E (fixed V-cycle, autoTune no) is
       identical to C since D-043. "E" is accepted on the command line as an
       alias of C (CONFIG_ALIASES).
    E variants (amendment C7): C plus one change each (the name is kept
    from the time E was the V-cycle variant), compared with C in the B10
    table:
       E-rcScalar   coupled.rhieChow.tensorial no (scalar D, C2)
       E-nonOrth60  remediation.meshQuality.nonOrthThreshold 60 (C1)
       E-nonOrth65  remediation.meshQuality.nonOrthThreshold 65 (C1 value; the
                    default stays 85, D-047/D-051)
       E-algPair    blockGAMG agglomerator algebraicPair (C4; the templates'
                    agglomerationWeights combined is replaced by pressure,
                    which algebraicPair means)
       E-noSFD      coupled.sfd.enabled no: the T3 template enables SFD
                    (D-058), so the variant is SFD OFF (was E-sfd = "SFD on",
                    a no-op on T3)
       E-eta07      etaMax 0.7 with minIter 2 ("solve loosely, iterate
                    often")
Scope (CONFIG_SCOPE, D-063): most configurations on the light cases T1,
T2, T3-SST, T3-GEKO; T4a only A, B, C, H; T4b only C; T5 none. The heavy
cases run one repeat (MAX_REPEATS, D-059). E-nonOrth60/65 change only the
static remediation of badly non-orthogonal cells, which exist on the
snappyHexMesh meshes only; with the heavy cases restricted by D-063 their
scope is empty (run them with --no-scope). --no-scope runs every
requested configuration on every requested case with --repeats repeats.

Single precision (amendment D11, D-065):
    SPn  simpleFoam in single precision: the settings of B (SIMPLEC),
         SP build. Compared with B (T1: with A, as B == A there).
    SPc  coupledFoam in single precision: the settings of C, SP
         build. Compared with C.
         Scope T1, T3-SST, T4a (D-063: SP vs DP on a few cases
         only), one repeat. SPn/SPc run from a shell with the SP
         OpenFOAM environment and the SP private install (refused
         in a DP shell, and the DP configurations in an SP shell);
         the case is prepared by tests/cflib/precision.py (mesh in
         DP, origin shift, checkMesh gate, SP solver settings;
         SP-geometry-fail -> recorded and skipped). Records carry
         the D11 fields (bench/SCHEMA_precision.md).

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

# Ranks of the heavy benchmark cases (user cap: 48 of the 64 cores, D-073)
HEAVY_NP = int(os.environ.get("CF_HEAVY_NP", "48"))

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
    # statWindow: the common D-042 averaging window of both solvers of the
    # case (D-068) = max(STAT_WINDOW_MIN, coupledFoam budget // 2)
    "T4a": {"template": "T4_motorBike", "args": ["-mesh", "a"],
            "monitor": "forces",
            "iters": {"simpleFoam": 3000, "coupledFoam": 800}, "np": HEAVY_NP,
            "oscillatory": True, "statWindow": 400},
    "T4b": {"template": "T4_motorBike", "args": ["-mesh", "b"],
            "monitor": "forces",
            "iters": {"simpleFoam": 4000, "coupledFoam": 800}, "np": HEAVY_NP,
            "oscillatory": True, "statWindow": 400},
    # CF_T5_MESH=coarse: development mesh as in tests/test_T5_ahmed.py (D-059)
    "T5": {"template": "T5_ahmed",
           "args": ["-mesh", os.environ.get("CF_T5_MESH", "fine")],
           "monitor": "forces",
           "iters": {"simpleFoam": 2000, "coupledFoam": 800}, "np": HEAVY_NP,  # user: 2000 for now (D-042 add. 3)
           "oscillatory": True, "statWindow": 400},
}

GAMG = "solvers.coupled.blockGAMG"

# coupledFoam configurations (D-068): fvSolution sets on top of the case
# template (missing sub-dictionaries are created by cflib.case.set_entry).
# C is the template itself: since D-043 every template runs a V-cycle with
# autoTune off.
COUPLED_CONFIGS = {
    "C": {},
    "D": {"solvers.coupled.preconditioner": "blockDiagonal"},
    "F": {"solvers.coupled.adaptiveRelTol": "no"},
    "G": {"coupled.anderson.enabled": "yes"},
    # fixed K-cycle: the B10 cycle comparison against C (V-cycle)
    "H": {f"{GAMG}.cycleType": "K", f"{GAMG}.autoTune": "no"},
    # K-cycle with the autoTune controller: the pre-D-043 default
    "H-tune": {f"{GAMG}.cycleType": "K", f"{GAMG}.autoTune": "yes"},
}

# Amendment C7: variants of the defaults C, one change each (named E-* since
# they were variants of the then-distinct V-cycle configuration E)
E_VARIANTS = {
    "E-rcScalar": {"coupled.rhieChow.tensorial": "no"},
    "E-nonOrth60": {"coupled.remediation.meshQuality.nonOrthThreshold": 60},
    "E-nonOrth65": {"coupled.remediation.meshQuality.nonOrthThreshold": 65},
    "E-algPair": {f"{GAMG}.agglomerator": "algebraicPair",
                  f"{GAMG}.agglomerationWeights": "pressure"},
    # the T3 template enables SFD (D-058): the variant switches it OFF
    "E-noSFD": {"coupled.sfd.enabled": "no"},
    "E-eta07": {"solvers.coupled.etaMax": 0.7,
                "solvers.coupled.minIter": 2},
}

CONFIGS = (("A", "B") + tuple(COUPLED_CONFIGS) + tuple(E_VARIANTS)
           + ("SPn", "SPc"))
NATIVE_CONFIGS = frozenset({"A", "B", "SPn"})
# Names accepted on the command line for a configuration of another name:
# the old E (fixed V-cycle, autoTune no) is the template default C since
# D-043
CONFIG_ALIASES = {"E": "C"}

# Amendment D11: single-precision configurations and their DP counterpart
# (same settings, DP build); Cd_rel_to_DP etc. are relative to it. SPn
# falls back to A where B is out of scope (T1: B == A).
SP_BASE = {"SPn": "B", "SPc": "C"}
SP_BASE_FALLBACK = {"SPn": "A"}
SP_CONFIGS = frozenset(SP_BASE)
# Repeats of the SP configurations (D-063: a few cases, 1 repeat)
CONFIG_REPEATS = {"SPn": 1, "SPc": 1}

LIGHT_CASES = frozenset({"T1", "T2", "T3-SST", "T3-GEKO"})

# Case scope of every configuration (D-063: most configurations on the light
# cases; T4a only A, B, C, H; T4b only C; T5 none). --no-scope overrides.
CONFIG_SCOPE = {
    "A": LIGHT_CASES | {"T4a"},
    # T1: B == A (the pitzDaily tutorial is SIMPLEC with p unrelaxed and
    # U, k, omega 0.9), so B would repeat A
    "B": frozenset({"T2", "T3-SST", "T3-GEKO", "T4a"}),
    "C": LIGHT_CASES | {"T4a", "T4b"},
    "D": LIGHT_CASES,
    "F": frozenset({"T1", "T3-SST"}),          # B10
    "G": frozenset({"T1", "T3-SST"}),          # B10
    "H": LIGHT_CASES | {"T4a"},
    "H-tune": frozenset({"T1", "T2"}),
    # D-063: the detailed variants on a few cases only (T2 was in the C7
    # scope of all of them; the heavy cases are restricted)
    "E-rcScalar": frozenset({"T2"}),
    # static remediation thresholds: only the snappyHexMesh meshes have such
    # cells, and D-063 restricts the heavy cases
    "E-nonOrth60": frozenset(),
    "E-nonOrth65": frozenset(),
    "E-algPair": frozenset({"T2"}),
    "E-noSFD": frozenset({"T3-SST", "T3-GEKO"}),   # SFD is on only in T3
    "E-eta07": frozenset({"T2"}),
    # D11 single precision, scoped per D-063
    "SPn": frozenset({"T1", "T3-SST", "T4a"}),
    "SPc": frozenset({"T1", "T3-SST", "T4a"}),
}

# Repeats of the heavy cases (D-059: one); the light cases use --repeats
MAX_REPEATS = {"T4a": 1, "T4b": 1, "T5": 1}

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
# Common window per case (D-068): CASES[case]["statWindow"] =
# max(STAT_WINDOW_MIN, coupledFoam budget // STAT_WINDOW_DIV), capped at the
# iterations run. The per-run rule max(300, n // 2) (D-042 add. 2) remains
# for cases without statWindow and as the *_perRun sensitivity value.
STAT_WINDOW_MIN = 300
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
# (4: D-068 - common wake window, configurations redefined, failed runs)
HARNESS_VERSION = 4

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
                         f"(known: {','.join(CONFIGS)}; aliases: "
                         f"{', '.join(f'{a}={b}' for a, b in CONFIG_ALIASES.items())})")
    return "simpleFoam" if cfg in NATIVE_CONFIGS else "coupledFoam"


def in_scope(name: str, cfg: str) -> bool:
    scope = CONFIG_SCOPE.get(cfg)
    return scope is None or name in scope


def resolve_configs(cfgs: list[str]) -> tuple[list[str], list[str]]:
    """Configuration names with aliases replaced (E -> C), duplicates
    removed in order; returns (names, notes about the aliases)."""
    out, notes = [], []
    for c in cfgs:
        if c in CONFIG_ALIASES:
            notes.append(f"configuration {c} is an alias of "
                         f"{CONFIG_ALIASES[c]} (D-068: identical since D-043)")
            c = CONFIG_ALIASES[c]
        if c not in out:
            out.append(c)
    return out, notes


def config_sets(cfg: str, spec: dict, name: str | None = None
                ) -> tuple[str, dict, dict]:
    """(solver, foamDictionary sets, files to write) of a configuration,
    including the per-case overrides of CASE_OVERRIDES for case `name`.
    SPn/SPc (D11) have the settings of their DP counterpart SP_BASE."""
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
        fv.update(COUPLED_CONFIGS.get(cfg, {}))
        fv.update(E_VARIANTS.get(cfg, {}))
        sets["system/fvSolution"] = fv
    for fname, entries in CASE_OVERRIDES.get((name, cfg), {}).items():
        sets.setdefault(fname, {}).update(entries)
    return solver, sets, files


def t5_variant() -> str:
    """T5 mesh variant of this process: CF_T5_MESH (default fine)."""
    return os.environ.get("CF_T5_MESH", "fine")


def t5_run_name(nprocs: int = HEAVY_NP, variant: str | None = None) -> str:
    """Name of the T5 test run (and, with the prefix ref_, of its
    reference): T5_np10 on the fine mesh, T5_<variant>_np10 otherwise (as
    tests/test_T5_ahmed.py names them)."""
    v = variant or t5_variant()
    return f"T5{'' if v == 'fine' else '_' + v}_np{nprocs}"


def case_args(name: str, variant: str | None = None) -> list[str]:
    """Allrun arguments of a case; for T5 the mesh variant `variant`
    (default: CF_T5_MESH of this process)."""
    if name == "T5":
        return ["-mesh", variant or t5_variant()]
    return list(CASES[name]["args"])


def mesh_variant(name: str, variant: str | None = None) -> str | None:
    """-mesh value of a case's Allrun arguments (T4a: a, T5: fine/coarse)."""
    args = case_args(name, variant)
    return args[args.index("-mesh") + 1] if "-mesh" in args else None


def config_hash(name: str, cfg: str, variant: str | None = None) -> str:
    """Hash of everything that defines a run of (case, configuration).
    variant: the T5 mesh variant of the run (default: CF_T5_MESH of this
    process); records store it (meshVariant) so the hash does not depend on
    the environment of the process that reads them."""
    spec = CASES[name]
    solver, sets, files = config_sets(cfg, spec, name)
    d = {"v": HARNESS_VERSION, "case": name, "cfg": cfg,
         "template": spec["template"], "args": case_args(name, variant),
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
                          "caseWindow": case_window(name),
                          "rel": STAT_REL, "abs": STAT_ABS,
                          "step": STAT_STEP}
        d["fieldAverageStart"] = field_average_start(spec["iters"][solver],
                                                     name)
    blob = json.dumps(d, sort_keys=True, default=str)
    return hashlib.sha256(blob.encode()).hexdigest()[:12]


def effective_settings(name: str, cfg: str) -> dict:
    """The normalised solver settings a run of (case, configuration)
    actually uses: the template dictionaries with the configuration's sets
    applied (tests/cflib/foamdict, no OpenFOAM needed). coupledFoam: the
    solvers.coupled and coupled dictionaries of fvSolution; simpleFoam:
    SIMPLE.consistent, the relaxation factors (p unrelaxed = 1) and the
    potentialFoam start. Two configurations with equal effective settings
    on a case are the same run (review C2)."""
    from cflib import foamdict  # noqa: PLC0415
    spec = CASES[name]
    tdir = cfcase.CASES / spec["template"] / "system"
    solver, sets, files = config_sets(cfg, spec, name)
    fvs = foamdict.read(tdir / "fvSolution")
    for entry, value in sets.get("system/fvSolution", {}).items():
        foamdict.set_dotted(fvs, entry, str(value))
    if solver == "coupledFoam":
        eff = {"solvers.coupled": foamdict.get(fvs, "solvers.coupled", {}),
               "coupled": foamdict.get(fvs, "coupled", {})}
    else:
        text = (files.get("system/relaxation.simpleFoam")
                or (tdir / "relaxation.simpleFoam").read_text())
        relax = foamdict.parse(text).get("relaxationFactors", {})
        fields = dict(relax.get("fields") or {})
        fields.setdefault("p", "1")          # no factor: p is not relaxed
        eff = {"consistent": foamdict.get(fvs, "SIMPLE.consistent", "no"),
               "fields": fields, "equations": dict(relax.get("equations")
                                                   or {}),
               "potentialStart": str(cfg == "A"
                                     and name in NATIVE_POTENTIAL_CASES)}
    return {"solver": solver, **foamdict.normalise(eff)}


def is_failed(rec: dict) -> bool:
    """True for the record of a failed run (M3): flagged failed, or an
    older record with a non-zero rc or an "error" and no time to
    convergence."""
    if rec.get("failed"):
        return True
    if "failed" in rec:
        return False
    return bool(rec.get("rc")) or ("error" in rec
                                   and rec.get("wall_to_conv_s") is None)


def is_current(rec: dict) -> bool:
    """True if a result record was produced with the current definition of
    its (case, configuration)."""
    try:
        return rec.get("configHash") == config_hash(
            rec["case"], rec["config"],
            rec.get("meshVariant") if rec["case"] == "T5" else None)
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
    window criterion (spec 12.3 ii) at the same time. None for an empty
    history."""
    if not hist:
        return None
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


def case_of_run(name: str | None) -> str | None:
    """CASES key of a case key, test/run name or benchmark tag:
    T4a, T4a_np10, ref_T4a_np10, bench_T4a_C_1, T4a_C_1 -> T4a;
    T5_coarse_np10, ref_T5_coarse_np10 -> T5; T3_kOmegaSST_np1 -> T3-SST,
    T3_GEKO_np1 -> T3-GEKO; T1_np4 -> T1. None if no case matches (T0)."""
    if not name:
        return None
    s = str(name)
    if s in CASES:
        return s
    for pre in ("ref_", "bench_"):
        if s.startswith(pre):
            s = s[len(pre):]
    if s in CASES:
        return s
    for pre, key in (("T3_kOmegaSST", "T3-SST"), ("T3_GEKO", "T3-GEKO")):
        if s.startswith(pre):
            return key
    for key in sorted(CASES, key=len, reverse=True):
        if s == key or s.startswith(key + "_"):
            return key
    return None


def case_window(case: str | None) -> int | None:
    """Common D-042 averaging window of a case (D-068), None for cases
    without one. `case` may be a CASES key or any run name (case_of_run)."""
    key = case_of_run(case)
    return CASES[key].get("statWindow") if key else None


def per_run_window(n: int) -> int:
    """The per-run window of D-042 addendum 2: max(300, n // 2) capped at n
    (informational since D-068)."""
    return max(0, min(max(STAT_WINDOW_MIN, n // STAT_WINDOW_DIV), n))


def stat_window(n: int, case: str | None = None) -> int:
    """Averaging window of a run of n iterations.

    With `case` (a CASES key or a run name such as T4a_np10 or
    ref_T4a_np10) that has a common window (D-068): that window, capped at
    n - the SAME window for both solvers of the case. Otherwise the per-run
    rule max(300, n // 2), capped at n (D-042 addendum 2). Callers that
    evaluate a wake case must pass the case; without it the result is the
    per-run sensitivity value."""
    w = case_window(case)
    if w is not None:
        return max(0, min(int(w), n))
    return per_run_window(n)


def field_average_start(n: int, case: str | None = None) -> int:
    """timeStart of the fieldAverage function object for a run of n
    iterations (deltaT 1, time = iteration): n - W + 1 with
    W = stat_window(n, case), so the mean fields cover exactly the force
    window, the last W iterations (the timeControl of a function object is
    active from time >= timeStart - 0.5 deltaT)."""
    return n - stat_window(n, case) + 1


def foam_dictionary(case: Path, fname: str, args: list[str],
                    capture: bool = False) -> str:
    """foamDictionary <args> <fname> in `case` (CLAUDE.md rule, as
    cflib.case.set_entry): fvSolution WITH -disableFunctionEntries (without
    it foamDictionary expands and drops the `#sinclude "relaxation"` when it
    rewrites the file); every other dictionary (controlDict, fvSchemes)
    WITHOUT it, because the flag writes their $-macros ($inletP,
    $turbulence) back quoted and broken."""
    flag = (["-disableFunctionEntries"] if Path(fname).name == "fvSolution"
            else [])
    out = subprocess.run(
        ["foamDictionary"] + flag + list(args) + [fname],
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
                        quantities: tuple[str, ...] = STAT_QUANTITIES,
                        case: str | None = None, w: int | None = None
                        ) -> int | None:
    """Smallest N <= n (iterations run) such that the window ending at N is
    stationary for all quantities; None if never (D-042).

    "The window" is W = w if given, else stat_window(n, case) - the case's
    common window (D-068) when `case` is given, the per-run window
    otherwise - slid along the history: N is scanned in steps of STAT_STEP
    (50) from the first multiple of 50 >= W, and finally N = n itself, so a
    stationary final window always yields a value. (Re-deriving W from N
    instead would shrink the window early in the run: under the original
    D-042 rule, on the T4a reference, a 75-iteration window at N = 150 was
    stationary.)"""
    if not hist:
        return None
    n = min(len(h) for h in hist.values())
    if w is None:
        w = stat_window(n, case)
    w = min(int(w), n)
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
                    quantities: tuple[str, ...] = STAT_QUANTITIES,
                    case: str | None = None) -> dict:
    """Record of the D-042 evaluation of a force history: W, means, stds,
    half-window drifts per quantity, `stationary` (final window, all
    quantities) and iters_to_stationary.

    W = stat_window(n, case): pass the case (CASES key or run name) so that
    both solvers of a wake case use its common window (D-068). The per-run
    window max(300, n // 2) is evaluated as well, for information only:
    W_perRun, iters_to_stationary_perRun, stationary_perRun and
    <q>_mean_perRun (the report may show both)."""
    if not hist or not all(q in hist for q in quantities):
        return {"criterion": STAT_CRITERION, "stationary": False,
                "iters_to_stationary": None,
                "iters_to_stationary_perRun": None}
    n = min(len(h) for h in hist.values())
    w = stat_window(n, case)
    w_run = per_run_window(n)
    out: dict = {"criterion": STAT_CRITERION, "iterations_run": n, "W": w,
                 "windowRule": ("common per case (D-068)"
                                if case_window(case) is not None
                                else "per run (D-042 addendum 2)"),
                 "W_perRun": w_run}
    for q in quantities:
        s = stationary_mean(hist, q, n, w)
        out.update({f"{q}_mean": s["mean"], f"{q}_std": s["std"],
                    f"{q}_mean1": s["mean1"], f"{q}_mean2": s["mean2"],
                    f"{q}_drift": s["drift"], f"{q}_driftTol": s["driftTol"],
                    f"{q}_stationary": s["stationary"]})
    out["stationary"] = all(out[f"{q}_stationary"] for q in quantities)
    out["iters_to_stationary"] = iters_to_stationary(hist, quantities, w=w)
    # sensitivity: the per-run window of D-042 addendum 2
    st_run = [stationary_mean(hist, q, n, w_run) for q in quantities]
    for q, s in zip(quantities, st_run):
        out[f"{q}_mean_perRun"] = s["mean"]
    out["stationary_perRun"] = all(s["stationary"] for s in st_run)
    out["iters_to_stationary_perRun"] = iters_to_stationary(hist, quantities,
                                                            w=w_run)
    return out


def reference_t_max(case: Path) -> float | None:
    """Last iteration of the ORIGINAL run of a cached simpleFoam reference
    that was continued with averaging (reference.json
    continuation.startTime; tests/test_T4_motorBike.py
    continue_reference_with_average), None for every other directory. Every
    evaluation of the reference must stop there (M7, D-060): the
    continuation samples are not part of the reference run."""
    meta = Path(case) / "reference.json"
    try:
        rec = json.loads(meta.read_text())
    except (OSError, ValueError):
        return None
    t = (rec.get("continuation") or {}).get("startTime")
    return float(t) if t is not None else None


def force_history(case: Path, t_min: float | None = None,
                  t_max: float | None = None,
                  original_only: bool = True) -> dict[str, list[float]]:
    """Cd / Cl history of a run (restarts merged), only samples with
    t_min < Time <= t_max. original_only (default): t_max defaults to
    reference_t_max(case), i.e. a continued reference is cut at the end of
    its original budget. {} if the case has no force coefficients."""
    try:
        fc = post.force_coeffs(case)
    except (FileNotFoundError, KeyError, OSError, ValueError, IndexError):
        return {}
    if "Cd" not in fc or "Cl" not in fc:
        return {}
    if original_only and t_max is None:
        t_max = reference_t_max(case)
    cd, cl = fc["Cd"], fc["Cl"]
    if "Time" in fc and (t_min is not None or t_max is not None):
        t = fc["Time"]
        keep = [(t_min is None or x > t_min) and (t_max is None or x <= t_max)
                for x in t]
        cd = [v for v, k in zip(cd, keep) if k]
        cl = [v for v, k in zip(cl, keep) if k]
    return {"Cd": [float(v) for v in cd], "Cl": [float(v) for v in cl]}


def stationary_point(case_dir: Path, solver: str, case: str | None) -> dict:
    """D-042 convergence point of a wake-case run directory, read-only:
    iterations, wall seconds and CPU-hours to the first stationary window
    under the case's common window (D-068) and, for information, under the
    per-run window (keys *_perRun). A continued reference is evaluated over
    its original budget only. {} if the run has no force history."""
    hist = force_history(case_dir)
    if not hist or not hist["Cd"]:
        return {}
    st = stationary_eval(hist, case=case)
    rt = rank_times(Path(case_dir), solver)
    out = {"n": st.get("iterations_run"), "W": st.get("W"),
           "W_perRun": st.get("W_perRun"),
           "stationary": st.get("stationary"),
           "Cd_mean": st.get("Cd_mean"), "Cl_mean": st.get("Cl_mean"),
           "wallTotal": rt.get("wallSeconds"), "cpuhTotal": rt.get("cpuHours")}
    for sfx, it in (("", st.get("iters_to_stationary")),
                    ("_perRun", st.get("iters_to_stationary_perRun"))):
        tc = (to_convergence(rt, progress_fraction(Path(case_dir), solver, it))
              if it else {})
        out[f"iters{sfx}"] = it
        out[f"wall{sfx}"] = tc.get("wall_to_conv_s")
        out[f"cpuh{sfx}"] = tc.get("cpu_to_conv_h")
    return out


def reference_timing_flags(ref_dir: Path | None, ref: dict) -> dict:
    """Fairness flags of a speed-up against a cached simpleFoam reference
    (review M2; the report must state them):
    referenceNoPotentialStart   the reference ran without the tutorial's
                                potentialFoam start (the harness never set
                                CF_NATIVE_POTENTIAL for the test references)
    referenceTimingConditionsUnknown  the reference record carries no
                                machine state (load, other jobs) of its run
    referenceSingleConfig       one native configuration (the tutorial
                                relaxation with SIMPLEC), not the best of the
                                benchmark configurations A/B"""
    ranks = ref.get("ranks") or {}
    pot = bool(ref_dir is not None
               and (Path(ref_dir) / "log.potentialFoam").exists()) \
        or "potentialFoam" in (ranks.get("preApps") or [])
    return {
        "referenceNoPotentialStart": not pot,
        "referenceTimingConditionsUnknown":
            not (ref.get("machineBefore") or ref.get("machine")),
        "referenceSingleConfig": True,
    }


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
    Returns {} if the solver left no report (run failed before start).
    Incomplete reports (a rank killed by MPI_Abort) never raise: the
    complete ones are used, their names are listed in incompleteReports and
    `complete` is False; if no report of the solver is complete, only
    {"ranks", "incompleteReports", "complete": False} is returned."""
    tdir = case / timing_dir
    reps = sorted(tdir.glob(f"{solver}.rank*.time"))
    if not reps:
        return {}
    sol = _parse_time_reports(reps)
    if "wallSeconds" not in sol:
        return {"ranks": sol["ranks"],
                "incompleteReports": sol["incompleteReports"],
                "complete": False}
    incomplete = list(sol["incompleteReports"])
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
        incomplete += p["incompleteReports"]
        if "wallSeconds" not in p:
            out[f"{app}Incomplete"] = True
            continue
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
    out["incompleteReports"] = incomplete
    out["complete"] = not incomplete
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
    if not rt or frac is None or rt.get("solverWallSeconds") is None:
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
    (SPn -> B, SPc -> C) on case `name`."""
    for base in (SP_BASE.get(cfg), SP_BASE_FALLBACK.get(cfg)):
        if base is None:
            continue
        recs, _ = load_current([name], [base])
        recs = [r for r in recs if not r.get("skipped") and not r.get("failed")]
        if recs:
            return recs
    return []


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
        "gateOverridden", "spSolverSettings",
        "wallSecondsDPMesh", "wallSecondsPrepare")}
    dps = dp_counterpart(name, cfg)
    rec["dpCounterpart"] = {"config": (dps[0].get("config") if dps
                                       else SP_BASE[cfg]),
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
        if old.get("configHash") == chash and not is_failed(old):
            print(f"skip {tag} (exists)")
            return old
        if old.get("configHash") == chash:
            print(f"{tag}: existing result is a FAILED run "
                  f"({'; '.join(old.get('failure') or ['?'])}): rerun")
        else:
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
        set_field_average_start(
            case, field_average_start(spec["iters"][solver], name))

    try:
        rc = cfcase.allrun(
            case, ["-solver", solver, "-np", str(nprocs)] + case_args(name),
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
           "meshVariant": mesh_variant(name),
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
    # M3: a failed run is recorded as failed, never as a timing (and is
    # rerun by the next invocation): Allrun rc, normal end of the solver
    # log, the whole budget run (the solver's own stop is disabled), all
    # timing reports complete
    hist: dict = {}
    failure: list[str] = []
    try:
        hist = monitor_history(case, spec["monitor"])
    except (FileNotFoundError, IndexError, KeyError, ValueError) as err:
        failure.append(f"monitor: {err}")
    n_run = min((len(h) for h in hist.values()), default=0)
    failure = cfcase.run_failure(case, solver, rc, n_run,
                                 spec["iters"][solver]) + failure
    if not rt:
        failure.append("no timing reports")
    elif not rt.get("complete", True):
        failure.append("incomplete timing reports: "
                       + ", ".join(rt.get("incompleteReports") or []))
    rec["iterations_run"] = n_run
    if failure:
        rec["failed"] = True
        rec["failure"] = failure
        rec["logTail"] = cfcase.log_tail(case, solver)
        precision_fields(name, cfg, case, rec)
        results.write("bench", tag, rec)
        print(f"{tag}: FAILED ({'; '.join(failure)})")
        return rec
    rec["failed"] = False

    osc = is_oscillatory(name)
    if osc:
        # D-042: convergence := stationary window mean of Cd and Cl; the
        # first stationary window gives the time to convergence. W is the
        # case's common window (D-068); the per-run window is recorded as
        # a sensitivity value (*_perRun)
        st = stationary_eval(hist, case=name)
        rec.update(st)
        it = st["iters_to_stationary"]
        it_run = st.get("iters_to_stationary_perRun")
        if it_run is not None and rt:
            tc_run = to_convergence(rt, progress_fraction(case, solver, it_run))
            rec["wall_to_conv_s_perRun"] = tc_run.get("wall_to_conv_s")
            rec["cpu_to_conv_h_perRun"] = tc_run.get("cpu_to_conv_h")
    else:
        rec["criterion"] = f"window {WINDOW}, tol {TOL} (12.3 ii)"
        it = iters_to_conv(hist)
    rec["iters_to_conv"] = it
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
            # D-068 sensitivity: the per-run window of D-042 addendum 2
            "W_perRun": _med(g, "W_perRun"),
            "iters_perRun_median": _med(g, "iters_to_stationary_perRun"),
            "wall_perRun_median": _med(g, "wall_to_conv_s_perRun"),
            "cpuh_perRun_median": _med(g, "cpu_to_conv_h_perRun"),
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
        if b and cc and b["wall_perRun_median"] and cc["wall_perRun_median"]:
            t["speedup_wall_B_over_C_perRun"] = (b["wall_perRun_median"]
                                                 / cc["wall_perRun_median"])
        # D11: DP counterpart over SP (SPn vs B, SPc vs C)
        base = (by.get((t["case"], SP_BASE.get(t["config"], "")))
                or by.get((t["case"], SP_BASE_FALLBACK.get(t["config"], ""))))
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
    """Amendment B10 comparisons against configuration C, wall and CPU-hours
    (configurations of D-068).

    H vs C: fixed K-cycle vs the default V-cycle (controller off in both):
    the cycle comparison. H-tune vs C, and vs H (dWall_X_vs_H): the autoTune
    controller. D (blockDiagonal), G (Anderson) and the E-* variants (one
    change of the defaults each) vs C: deltas only, B10 gives no threshold.
    F vs C (fixed relTol): pass if adaptive C is at most 5 % slower than F
    in wall time AND in CPU-hours, and the monitored quantity (Cd for force
    cases, dp for T1/T2) at the end of the fixed budget agrees to 1e-4
    (relative). Deltas are X/C - 1 (positive: X slower than C)."""
    by = {(t["case"], t["config"]): t for t in table}
    out = []
    for (case, cfg), t in sorted(by.items()):
        if cfg in NATIVE_CONFIGS or cfg == "C" or cfg not in CONFIGS:
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
        if cfg == "H-tune":
            # the controller alone: H-tune vs the fixed K-cycle H
            h = by.get((case, "H"))
            if h:
                row["dWall_X_vs_H"] = _rel(t["wall_median"], h["wall_median"])
                row["dCpu_X_vs_H"] = _rel(t["cpuh_median"], h["cpuh_median"])
                row["iters_H"] = h["iters_median"]
        row["cycleTypeFinal"] = t.get("cycleTypeFinal")
        out.append(row)
    return out


def _bench_records(names: list[str] | None = None,
                   cfgs: list[str] | None = None) -> list[tuple[str, dict]]:
    out = []
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
        out.append((f.stem, rec))
    return out


def load_failed(names: list[str] | None = None,
                cfgs: list[str] | None = None) -> list[dict]:
    """Current records of FAILED runs (M3): excluded from load_current and
    the summary, to be listed as failed/missing in the report and rerun."""
    return [rec for _, rec in _bench_records(names, cfgs)
            if is_current(rec) and is_failed(rec)]


def expected_runs(names: list[str], cfgs: list[str], repeats: int,
                  no_scope: bool = False) -> list[tuple[str, str, int]]:
    """(case, configuration, run) that an invocation with these arguments
    produces (scope and MAX_REPEATS as in main)."""
    out = []
    for name in names:
        rep = repeats if no_scope else min(repeats,
                                           MAX_REPEATS.get(name, repeats))
        for run in range(1, rep + 1):
            for cfg in cfgs:
                if no_scope or in_scope(name, cfg):
                    out.append((name, cfg, run))
    return out


def load_current(names: list[str] | None = None,
                 cfgs: list[str] | None = None) -> tuple[list[dict], list[str]]:
    """Current benchmark records of successful runs (optionally filtered)
    and the stale tags. Records of failed runs are left out (load_failed)."""
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
        if is_failed(rec):
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
    "W_perRun", "iters_perRun_median", "wall_perRun_median",
    "cpuh_perRun_median", "speedup_wall_B_over_C_perRun",
    "precision", "Cd_rel_to_DP", "Cl_rel_to_DP", "dp_rel_to_DP",
    "staticSetSizeDiff", "spGeometry", "nCheckMeshDiff",
    "speedup_wall_DP_over_SP", "speedup_cpu_DP_over_SP"]
B10_FIELDS = [
    "case", "config", "reference", "status", "pass", "criterion",
    "wall_X", "wall_C", "dWall_X_vs_C", "cpuh_X", "cpuh_C", "dCpu_X_vs_C",
    "slowdownWall_C_vs_F", "slowdownCpu_C_vs_F", "monitor", "monitorRelDiff",
    "dWall_X_vs_H", "dCpu_X_vs_H",
    "cycleTypeFinal", "iters_X", "iters_C", "iters_H"]


def summarise(names: list[str], cfgs: list[str], repeats: int | None = None,
              no_scope: bool = False) -> Path:
    """summary.csv, b10_acceptance.csv, summary.json. Failed runs (M3) are
    excluded from the table and listed under "failed"; with `repeats`, the
    expected (case, config, run) without a successful current record are
    listed under "missing"."""
    recs, stale = load_current(names, cfgs)
    for s in stale:
        print(f"summary: stale result {s} excluded (configHash mismatch)")
    failed = [{"tag": f"{r['case']}_{r['config']}_{r.get('run')}",
               "failure": r.get("failure") or
               ([f"rc {r.get('rc')}"] + ([r["error"]] if r.get("error") else []))}
              for r in load_failed(names, cfgs)]
    for f in failed:
        print(f"summary: FAILED run {f['tag']} excluded "
              f"({'; '.join(map(str, f['failure']))})")
    missing = []
    if repeats:
        have = {(r["case"], r["config"], r.get("run")) for r in recs}
        missing = [f"{c}_{k}_{r}" for c, k, r in
                   expected_runs(names, cfgs, repeats, no_scope)
                   if (c, k, r) not in have]
        for m in missing:
            print(f"summary: MISSING run {m} (no successful current record)")
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
         "stale": stale, "failed": failed, "missing": missing},
        indent=2, default=str) + "\n")
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cases", default="T1,T2,T3-SST,T3-GEKO")
    ap.add_argument("--configs", default=",".join(CONFIGS),
                    help="every configuration runs only on the cases of "
                         "its scope (CONFIG_SCOPE, D-063; see --no-scope); "
                         "E is an alias of C")
    ap.add_argument("--no-scope", action="store_true",
                    help="run every requested configuration on every "
                         "requested case, with --repeats repeats also on "
                         "the heavy cases")
    ap.add_argument("--repeats", type=int, default=3,
                    help="repeats (heavy cases: at most MAX_REPEATS, "
                         "D-059, unless --no-scope)")
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
                  f"monitor={v['monitor']} "
                  f"repeats<={MAX_REPEATS.get(k, a.repeats)} "
                  f"configs={','.join(cf) or '-'}")
        return 0
    names = [n for n in a.cases.split(",") if n]
    cfgs, alias_notes = resolve_configs([c for c in a.configs.split(",") if c])
    for s in alias_notes:
        print(s)
    for n in names:
        if n not in CASES:
            ap.error(f"unknown case {n!r} (known: {','.join(CASES)})")
    for c in cfgs:
        if c not in CONFIGS:
            ap.error(f"unknown configuration {c!r} (known: {','.join(CONFIGS)}"
                     f"; aliases {CONFIG_ALIASES})")
    n_failed = 0
    if not a.summary_only:
        cfenv.foam_env()
        # D11: precision is a build choice - SP configurations from an SP
        # shell only, DP configurations from a DP shell only
        sp_cfgs = [c for c in cfgs if c in SP_CONFIGS]
        if sp_cfgs and len(sp_cfgs) != len(cfgs):
            if "--configs" in " ".join(sys.argv):
                ap.error("run the SP configurations (SPn, SPc) separately, "
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
        for name in names:
            nprocs = a.np or CASES[name]["np"]
            repeats = (a.repeats if a.no_scope
                       else min(a.repeats, MAX_REPEATS.get(name, a.repeats)))
            for run in range(1, repeats + 1):
                for cfg in cfgs:
                    if not a.no_scope and not in_scope(name, cfg):
                        continue
                    if (not a.no_scope and cfg in CONFIG_REPEATS
                            and run > CONFIG_REPEATS[cfg]):
                        continue
                    st = cfenv.machine_state()
                    if st.busy and not a.allow_busy:
                        print("machine busy (other jobs: %d, load %.1f): "
                              "refusing to benchmark; --allow-busy overrides"
                              % (len(st.jobs), st.load))
                        return 3
                    if run_one(name, cfg, run, nprocs, a.force).get("failed"):
                        n_failed += 1
    print("summary:", summarise(names, cfgs, a.repeats, a.no_scope))
    if n_failed:
        print(f"{n_failed} benchmark run(s) FAILED (see summary.json "
              "\"failed\"); they are rerun by the next invocation")
        return 4
    return 0


if __name__ == "__main__":
    sys.exit(main())

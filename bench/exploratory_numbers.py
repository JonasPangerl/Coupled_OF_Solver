#!/usr/bin/env python3
"""Numbers quoted in the paper's discussion that do not come from the test
protocol: exploratory cavity runs, the deferred-correction symbol analysis,
and measured values recorded in DECISIONS.md.

Writes results/exploratory/<name>.json. Every value carries its source:
a log file under run/ (read-only), this script's own computation, or a
DECISIONS.md entry where the log no longer exists. bench/make_report.py
turns the values into numbers.tex macros (\\cfExp...), so the paper carries
no hand-copied numbers.

    results/exploratory/re1000.json      run/exp_re1000/<variant>/log.coupledFoam
    results/exploratory/linsolver.json   run/exp_linsolver/<variant>/log.coupledFoam
    results/exploratory/symbol.json      Fourier symbol of the deferred-
                                         correction iteration (paper eq. dcG)
    results/exploratory/decisions.json   D-022 and D-023 measured values

Usage: bench/exploratory_numbers.py        (no OpenFOAM environment needed)
A run directory that is missing is reported and its values are omitted
(never invented).
"""

from __future__ import annotations

import json
import re
import sys
import time
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "tests"))
from cflib import logs  # noqa: E402

RUN = REPO / "run"
OUT = REPO / "results" / "exploratory"

R_TARGET = 1e-8


# --------------------------------------------------------------------------- #
# exploratory runs
# --------------------------------------------------------------------------- #

def _settings(log: Path) -> dict:
    """Linear-solver and PTC settings from the 'effective settings' dump."""
    text = log.read_text(errors="replace")
    out = {}
    i = text.find("linearSolver")
    if i >= 0:
        block = text[i:i + 3000]
        for key in ("solver", "relTol", "maxIter", "smoother", "cycleType"):
            m = re.search(rf"^\s+{key}\s+(\S+);", block, re.M)
            if m:
                out[key] = logs._num(m.group(1))
    for key in ("CFLmax", "startupUpwindIters"):
        m = re.search(rf"^\s+{key}\s+(\S+);", text, re.M)
        if m:
            out[key] = logs._num(m.group(1))
    m = re.search(r"start-up phase done at iteration (\d+)", text)
    out["startupDoneAt"] = int(m.group(1)) if m else None
    return out


def run_stats(case: Path) -> dict | None:
    log = case / "log.coupledFoam"
    if not log.exists():
        return None
    rows = logs.parse_cf(log)
    if not rows:
        return None
    s = _settings(log)
    R = [r["R"] for r in rows]
    cfl = [r["CFL"] for r in rows]
    cflmax = s.get("CFLmax")
    first_max = next((r["iter"] for r in rows
                      if cflmax is not None and r["CFL"] >= cflmax), None)
    below = next((r["iter"] for r in rows if r["R"] < R_TARGET), None)
    sd = s.get("startupDoneAt")
    post = [r["CFL"] for r in rows if sd is not None and r["iter"] > sd]
    maxit = s.get("maxIter")
    lin = [r.get("linIters", 0) for r in rows]
    return {
        "source": str(log.relative_to(REPO)),
        "settings": s,
        "iterations": len(rows),
        "firstIterBelowTarget": below,
        "target": R_TARGET,
        "finalR": R[-1],
        "minR": min(R),
        "cflCutsTotal": int(sum(r.get("cuts", 0) for r in rows)),
        "firstIterAtCFLmax": first_max,
        "finalCFL": cfl[-1],
        "cflMinAfterStartup": min(post) if post else None,
        "cflMaxAfterStartup": max(post) if post else None,
        "cflP05AfterStartup": float(np.percentile(post, 5)) if post else None,
        "cflP95AfterStartup": float(np.percentile(post, 95)) if post else None,
        "linItersTotal": int(sum(lin)),
        "linSolvesAtMaxIter": int(sum(1 for n in lin
                                      if maxit is not None and n >= maxit)),
        "tSolveTotal_s": float(sum(r.get("tSolve", 0.0) for r in rows)),
        "tIterTotal_s": float(sum(r.get("tIter", 0.0) for r in rows)),
    }


def exploratory(dirname: str) -> dict:
    base = RUN / dirname
    out = {"source": f"run/{dirname}", "runs": {}, "missing": []}
    if not base.is_dir():
        out["missing"].append(str(base))
        return out
    for d in sorted(p for p in base.iterdir() if p.is_dir()):
        st = run_stats(d)
        if st is None:
            out["missing"].append(str(d.relative_to(REPO)))
        else:
            out["runs"][d.name] = st
    return out


# --------------------------------------------------------------------------- #
# symbol analysis of the deferred-correction iteration (paper eq. dcG)
# --------------------------------------------------------------------------- #

def dc_symbol_rho(pe: float, cfl: float, theta_deg: float = 45.0,
                  n_cells: int = 128, m: int = 2049) -> dict:
    """max |G(xi, eta)| of G = (A_UD + T)^-1 (T - C) for constant convection
    at angle theta on a uniform 2D grid (h = 1, |U| = 1, nu = 1/Pe_h).

    A_UD: first-order upwind convection + central diffusion (5-point).
    A_HO: linearUpwind = upwind cell value + (h/2) x central gradient of the
    upwind cell (Fromm in 1D); C = A_HO - A_UD.
    T = lambda / CFL with lambda = 0.5 sum_f |phi_f| + nu V^(1/3)
      = |u| + |v| + nu (paper eq. dt with h = 1).
    Modes with wavelength longer than the domain (|k| < 2 pi / n_cells) are
    excluded; the mode set is otherwise continuous (m x m samples of
    [-pi, pi]^2)."""
    u = np.cos(np.radians(theta_deg))
    v = np.sin(np.radians(theta_deg))
    nu = 1.0 / pe
    k = np.linspace(-np.pi, np.pi, m)
    X, Y = np.meshgrid(k, k)

    def ud(a, x):
        return a * (1 - np.exp(-1j * x))

    def corr(a, x):
        return a * (np.exp(1j * x) - np.exp(-1j * x) - 1 + np.exp(-2j * x)) / 4

    diff = nu * (4 - 2 * np.cos(X) - 2 * np.cos(Y))
    a_ud = ud(u, X) + ud(v, Y) + diff
    c = corr(u, X) + corr(v, Y)
    t = 0.0 if np.isinf(cfl) else (abs(u) + abs(v) + nu) / cfl
    mask = np.hypot(X, Y) >= 2 * np.pi / n_cells
    with np.errstate(divide="ignore", invalid="ignore"):
        g = np.abs((t - c) / (a_ud + t))
    g = np.where(mask, g, 0.0)
    i = np.unravel_index(np.argmax(g), g.shape)
    return {"rho": float(g[i]), "argmaxKnorm": float(np.hypot(X[i], Y[i])),
            "kMin": 2 * np.pi / n_cells}


def symbol() -> dict:
    re_cav, n = 1000.0, 128
    pe1 = re_cav / n                  # Re 1000 cavity, lid velocity, 128^2
    pe2 = 10 * pe1
    cases = {
        "peLow": pe1, "peHigh": pe2,
        "rhoInfPeLow": dc_symbol_rho(pe1, np.inf, n_cells=n),
        "rhoInfPeHigh": dc_symbol_rho(pe2, np.inf, n_cells=n),
        "rhoCflFiveHundredPeLow": dc_symbol_rho(pe1, 500.0, n_cells=n),
    }
    # 1D check: symbol of A_UD^-1 C for pure convection is (i/2) sin(xi)
    xi = np.linspace(1e-3, np.pi, 1001)
    oned = ((np.exp(1j * xi) - np.exp(-1j * xi) - 1 + np.exp(-2j * xi)) / 4
            / (1 - np.exp(-1j * xi)))
    cases["oneDimMaxAbs"] = float(np.max(np.abs(oned)))
    cases["oneDimMaxErrToHalfSin"] = float(np.max(np.abs(oned - 0.5j * np.sin(xi))))
    cases["source"] = ("bench/exploratory_numbers.py:dc_symbol_rho "
                       "(theta 45 deg, uniform grid, modes |k| >= 2 pi/128)")
    return cases


# --------------------------------------------------------------------------- #
# values measured earlier and recorded only in DECISIONS.md
# --------------------------------------------------------------------------- #

def decisions() -> dict:
    return {
        "dTwoTwoFrozenR": {
            "value": 7.0e-8,
            "source": "DECISIONS.md D-022 (observed R frozen at 7.0e-8 on T0; "
                      "log not kept)"},
        "dTwoThreeRelDiffTolEightMinus": {
            "value": 1.32e-5,
            "source": "DECISIONS.md D-023 table (tolerance 1e-8, np1 vs np4)"},
        "dTwoThreeRelDiffTolNineMinus": {
            "value": 9.3e-7,
            "source": "DECISIONS.md D-023 table (tolerance 1e-9)"},
        "dTwoThreeRelDiffFloor": {
            "value": 6.0e-7,
            "source": "DECISIONS.md D-023 table (tolerance 1e-10 and 1e-11)"},
    }


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    stamp = time.strftime("%Y-%m-%dT%H:%M:%S")
    items = {
        "re1000": exploratory("exp_re1000"),
        "linsolver": exploratory("exp_linsolver"),
        "symbol": symbol(),
        "decisions": decisions(),
    }
    for name, d in items.items():
        d["generatedBy"] = "bench/exploratory_numbers.py"
        d["timestamp"] = stamp
        (OUT / f"{name}.json").write_text(json.dumps(d, indent=2) + "\n")
        miss = d.get("missing") or []
        print(f"{name}: written" + (f", missing {miss}" if miss else ""))
    return 0


if __name__ == "__main__":
    sys.exit(main())

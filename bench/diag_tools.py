#!/usr/bin/env python3
"""Analysis companion of the coupledFoam deep diagnostics (TASK 5, D-045).

The solver writes <case>/diagnostics/diag.rank<N>.jsonl when
coupled.diagnostics.level >= 1 (one JSON object per line; line 1 is a
"header" record, then one "iter" record per outer iteration, possibly a
final "truncated" record).

API
    records(case)              {rank: [iteration records]} (raw dicts)
    load(case)                 pandas.DataFrame of the level-1 records,
                               index = outer iteration; per-rank files
                               merged: global fields from the master (rank
                               0), rank-local fields (header "localKeys")
                               aggregated over the ranks (see LOCAL_AGG)
    linear_history(case, it)   level-2 linear-solve records of iteration
                               `it` (list of dicts; arrays as numpy)
    level_reduction(case)      DataFrame iteration x GAMG level: geometric
                               mean over the preconditioner applications of
                               the per-level reduction ||r_post||/||b||
    fig_convergence, fig_gamg_heatmap, fig_time_breakdown
                               the three canned figures (matplotlib)
    make_figures(case, out)    all three as vector PDF + PNG

Command line
    bench/diag_tools.py <case> [--out DIR]    figures into DIR (default
                                              <case>/diagnostics/figures)
    bench/diag_tools.py <case> --report       figures next to the report
                                              figures: report/figures/*.png
                                              and report/paper/figures/*.pdf
                                              (prefix diag_<case name>_)
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]

PHASE_COLOR = {
    "startup": "#fdd0a2",
    "ramp": "#c6dbef",
    "asymptotic": "#c7e9c0",
    "stalled": "#fcbba1",
}

# Aggregation of rank-local level-1 fields over the ranks
LOCAL_AGG = {
    "residuals.massErrMax": "max",
    "residuals.massErrSum": "sum",
    "controls.dt.min": "min",
    "controls.dt.max": "max",
    # median of the rank medians (not the global median: no gather)
    "controls.dt.median": "median",
    "turbulence.nBoundK": "sum",
    "turbulence.nBoundOmega": "sum",
    "memory.rssKB": "sum",
    "memory.peakRssKB": "sum",
}
# timings.*: the slowest rank (critical path)
TIMING_AGG = "max"

# Nested level-1 blocks that are not flattened into columns
SKIP_KEYS = {"linear", "gamgSetup", "bcFlips", "andersonInternals"}


# * * * * * * * * * * * * * * * * * Reading * * * * * * * * * * * * * * * * //

def _diag_dir(case) -> Path:
    case = Path(case)
    d = case / "diagnostics" if (case / "diagnostics").is_dir() else case
    if not any(d.glob("diag.rank*.jsonl")):
        raise FileNotFoundError(f"no diag.rank*.jsonl under {d}")
    return d


def _rank_of(path: Path) -> int:
    return int(path.stem.replace("diag.rank", ""))


def read_rank_file(path) -> tuple[dict, list[dict], dict | None]:
    """(header, iteration records, truncation record or None) of one file.
    A restarted run appends a second header; later records of the same
    iteration replace earlier ones."""
    header: dict = {}
    recs: dict[int, dict] = {}
    trunc = None
    with open(path) as f:
        for line in f:
            if not line.strip():
                continue
            r = json.loads(line)
            t = r.get("type")
            if t == "header":
                header = header or r
            elif t == "iter":
                recs[int(r["iter"])] = r
            elif t == "truncated":
                trunc = r
    return header, [recs[k] for k in sorted(recs)], trunc


def records(case) -> dict[int, list[dict]]:
    """{rank: [iteration records]} of all rank files of a case."""
    out = {}
    for p in sorted(_diag_dir(case).glob("diag.rank*.jsonl"), key=_rank_of):
        out[_rank_of(p)] = read_rank_file(p)[1]
    return out


def headers(case) -> dict[int, dict]:
    out = {}
    for p in sorted(_diag_dir(case).glob("diag.rank*.jsonl"), key=_rank_of):
        out[_rank_of(p)] = read_rank_file(p)[0]
    return out


def _flatten(d: dict, prefix: str = "", out: dict | None = None) -> dict:
    out = {} if out is None else out
    for k, v in d.items():
        if not prefix and k in SKIP_KEYS:
            continue
        key = f"{prefix}{k}"
        if isinstance(v, dict):
            _flatten(v, key + ".", out)
        else:
            out[key] = v
    return out


def _bc_flips(rec: dict) -> float:
    bc = rec.get("bcFlips")
    if not bc:
        return math.nan
    vals = [p.get("nFlips") for p in bc.values()]
    vals = [v for v in vals if v is not None]
    return float(sum(vals)) if vals else math.nan


def load(case) -> pd.DataFrame:
    """Level-1 records as a DataFrame, index = outer iteration.

    Global fields (identical on all ranks, already reduced by the solver)
    come from the master file; rank-local fields are aggregated with
    LOCAL_AGG, timings with the maximum over the ranks. Lists (trials,
    cellsPerLevel, ...) stay as Python objects in their columns. Level >= 3:
    column bcFlips.total = inflow/outflow switches summed over the patches
    and ranks.
    """
    per_rank = records(case)
    if not per_rank:
        raise ValueError("no iteration records")
    frames = {}
    for rank, recs in per_rank.items():
        rows = []
        for r in recs:
            flat = _flatten(r)
            flat["bcFlips.total"] = _bc_flips(r)
            rows.append(flat)
        frames[rank] = pd.DataFrame(rows).set_index("iter")
    master = frames[min(frames)].copy()
    if len(frames) > 1:
        cols = set(LOCAL_AGG) | {c for c in master.columns
                                 if c.startswith("timings.")}
        cols.add("bcFlips.total")
        for c in cols:
            if c not in master.columns:
                continue
            how = LOCAL_AGG.get(c, TIMING_AGG if c.startswith("timings.")
                                else "sum")
            stack = pd.concat([f[c] for f in frames.values() if c in f],
                              axis=1)
            master[c] = getattr(stack, how)(axis=1)
    master["nRanks"] = len(frames)
    return master


def linear_history(case, it: int, rank: int = 0) -> list[dict]:
    """Level-2 linear solves of outer iteration `it` (one entry per solve
    attempt; CFL cuts give several). Arrays converted to numpy."""
    recs = records(case).get(rank, [])
    for r in recs:
        if r["iter"] == it:
            if "linear" not in r:
                raise ValueError("no level-2 data (diagnostics level < 2)")
            out = []
            for s in r["linear"]:
                s = dict(s)
                s["krylov"] = np.asarray(s["krylov"], dtype=float)
                s["trueResidualAtRestart"] = np.asarray(
                    s["trueResidualAtRestart"], dtype=float)
                out.append(s)
            return out
    raise KeyError(f"iteration {it} not in the diagnostics of rank {rank}")


def level_reduction(case, rank: int = 0) -> pd.DataFrame:
    """Iteration x level: geometric mean over all visits of the level in
    all preconditioner applications of the iteration of
    ||r after post-smoothing|| / ||b|| (the reduction the cycle achieved on
    that level). Needs level >= 2; with upLeg first (default) only the
    first application of each solve carries the up leg."""
    rows = {}
    for r in records(case).get(rank, []):
        if "linear" not in r:
            continue
        acc: dict[int, list[float]] = {}
        for s in r["linear"]:
            for app in s["precon"]:
                for v in app["levels"]:
                    if not v.get("post"):
                        continue    # up leg not logged (upLeg first)
                    b, a = v["pre"][0], v["post"][1]
                    if b and a is not None and b > 0 and a > 0:
                        acc.setdefault(v["l"], []).append(math.log10(a / b))
        if acc:
            rows[r["iter"]] = {l: 10 ** float(np.mean(x))
                               for l, x in acc.items()}
    return pd.DataFrame.from_dict(rows, orient="index").sort_index(axis=1)


# * * * * * * * * * * * * * * * * * Figures * * * * * * * * * * * * * * * * //

def _plt():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    return plt


def _phase_bands(ax, df):
    if "phase" not in df:
        return
    it = df.index.to_numpy()
    ph = df["phase"].to_numpy()
    start = 0
    for i in range(1, len(ph) + 1):
        if i == len(ph) or ph[i] != ph[start]:
            ax.axvspan(it[start] - 0.5, it[i - 1] + 0.5,
                       color=PHASE_COLOR.get(ph[start], "#eeeeee"),
                       lw=0, zorder=0)
            start = i


def fig_convergence(df: pd.DataFrame, title: str = ""):
    """(1) R, CFL, eta, rho vs iteration with phase bands."""
    plt = _plt()
    fig, axes = plt.subplots(4, 1, figsize=(7, 8), sharex=True)
    series = [("residuals.R", "R", True), ("controls.CFL", "CFL", True),
              ("controls.eta", "eta", True), ("residuals.rho", "rho", True)]
    for ax, (col, lab, logy) in zip(axes, series):
        _phase_bands(ax, df)
        if col in df:
            y = pd.to_numeric(df[col], errors="coerce")
            if logy:
                y = y.where(y > 0)
            ax.plot(df.index, y, lw=1.0, color="#1f5fbf")
            if logy and y.notna().any():
                ax.set_yscale("log")
        ax.set_ylabel(lab)
        ax.grid(alpha=0.3, which="both", lw=0.4)
    axes[-1].set_xlabel("outer iteration")
    handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in PHASE_COLOR.values()]
    axes[0].legend(handles, list(PHASE_COLOR), ncol=4, fontsize=7,
                   loc="upper right")
    if title:
        axes[0].set_title(title)
    fig.tight_layout()
    return fig


def fig_gamg_heatmap(red: pd.DataFrame, title: str = ""):
    """(2) per-level GAMG residual reduction over the iterations."""
    plt = _plt()
    fig, ax = plt.subplots(figsize=(7, 3.2))
    if red.empty:
        ax.text(0.5, 0.5, "no level-2 data", ha="center", va="center")
        ax.set_axis_off()
        return fig
    data = np.log10(red.to_numpy(dtype=float).T)
    im = ax.imshow(data, aspect="auto", origin="lower", cmap="viridis_r",
                   interpolation="nearest",
                   extent=[red.index.min() - 0.5, red.index.max() + 0.5,
                           -0.5, data.shape[0] - 0.5])
    ax.set_yticks(range(data.shape[0]))
    ax.set_yticklabels([str(c) for c in red.columns])
    ax.set_xlabel("outer iteration")
    ax.set_ylabel("GAMG level")
    cb = fig.colorbar(im, ax=ax)
    cb.set_label("log10 ||r_post|| / ||b||")
    if title:
        ax.set_title(title)
    fig.tight_layout()
    return fig


TIME_PARTS = [
    ("timings.tMomentumOps", "momentum ops"),
    ("timings.tBoundary", "boundary coupling"),
    ("timings.tContinuity", "continuity/scaling"),
    ("timings.tRhieChow", "Rhie-Chow"),
    ("asmOther", "assembly other (dt, local limit)"),
    ("timings.tPrecSetup", "precond. setup"),
    ("timings.tPrecApply", "precond. apply"),
    ("timings.tKrylov", "Krylov vector ops"),
    ("timings.tFlux", "flux update"),
    ("timings.tTurb", "turbulence"),
    ("timings.tDiag", "diagnostics"),
    ("iterOther", "other"),
]


def time_breakdown(df: pd.DataFrame) -> pd.DataFrame:
    """Time parts per iteration [s]; asmOther and iterOther close the sums
    to tAsm and tIter."""
    t = pd.DataFrame(index=df.index)
    for c, _ in TIME_PARTS:
        if c in df:
            t[c] = pd.to_numeric(df[c], errors="coerce").fillna(0.0)
    asm = [c for c in ("timings.tMomentumOps", "timings.tBoundary",
                       "timings.tContinuity", "timings.tRhieChow") if c in t]
    if "timings.tAsm" in df:
        t["asmOther"] = (pd.to_numeric(df["timings.tAsm"]) - t[asm].sum(axis=1)
                         ).clip(lower=0)
    if "timings.tIter" in df:
        known = [c for c, _ in TIME_PARTS if c in t and c != "iterOther"]
        t["iterOther"] = (pd.to_numeric(df["timings.tIter"])
                          - t[known].sum(axis=1)).clip(lower=0)
    return t[[c for c, _ in TIME_PARTS if c in t]]


def fig_time_breakdown(df: pd.DataFrame, title: str = ""):
    """(3) stacked time breakdown per iteration."""
    plt = _plt()
    t = time_breakdown(df)
    fig, ax = plt.subplots(figsize=(7, 3.6))
    labels = dict(TIME_PARTS)
    cmap = plt.get_cmap("tab20")
    ax.stackplot(t.index, [t[c].to_numpy() for c in t.columns],
                 labels=[labels[c] for c in t.columns],
                 colors=[cmap(i) for i in range(len(t.columns))], lw=0)
    ax.set_xlabel("outer iteration")
    ax.set_ylabel("time per iteration [s]")
    ax.legend(fontsize=6, ncol=2, loc="upper right")
    ax.grid(alpha=0.3, lw=0.4)
    if title:
        ax.set_title(title)
    fig.tight_layout()
    return fig


def make_figures(case, out=None, prefix: str = "diag_") -> list[Path]:
    """Write the three figures as vector PDF and PNG; returns the paths.
    out: a directory, or the tuple (png_dir, pdf_dir)."""
    plt = _plt()
    case = Path(case)
    if out is None:
        out = case / "diagnostics" / "figures"
    png_dir, pdf_dir = (out if isinstance(out, tuple) else (out, out))
    png_dir, pdf_dir = Path(png_dir), Path(pdf_dir)
    png_dir.mkdir(parents=True, exist_ok=True)
    pdf_dir.mkdir(parents=True, exist_ok=True)
    df = load(case)
    name = case.name
    figs = {
        "convergence": fig_convergence(df, name),
        "gamg_levels": fig_gamg_heatmap(level_reduction(case), name),
        "time_breakdown": fig_time_breakdown(df, name),
    }
    written = []
    for key, fig in figs.items():
        for d, ext in ((pdf_dir, "pdf"), (png_dir, "png")):
            p = d / f"{prefix}{key}.{ext}"
            fig.savefig(p, dpi=150)
            written.append(p)
        plt.close(fig)
    return written


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("case")
    ap.add_argument("--out", default=None, help="output directory")
    ap.add_argument("--report", action="store_true",
                    help="write next to the report figures")
    a = ap.parse_args()
    case = Path(a.case)
    if a.report:
        out = (REPO / "report" / "figures", REPO / "report" / "paper" / "figures")
        paths = make_figures(case, out, prefix=f"diag_{case.name}_")
    else:
        paths = make_figures(case, a.out)
    df = load(case)
    print(f"{case}: {len(df)} iterations, {int(df['nRanks'].iloc[0])} rank(s),"
          f" phases {df['phase'].value_counts().to_dict()}")
    for p in paths:
        print(p)


if __name__ == "__main__":
    main()

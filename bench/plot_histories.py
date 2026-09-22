#!/usr/bin/env python3
"""Per-case iteration histories for the report appendix (matplotlib).

For every test run (T0 Re100/Re1000 np1+np4, T1 np1+np4, T2 np1+np4, T3
SST/GEKO, T4a, T4b, T5) and its simpleFoam reference this writes vector
PDFs into report/paper/figures/:

    hist_<run>_loads.pdf      monitored quantities per iteration (left) and
                              over wall-clock time (right), both solvers:
                              Cd, Cl, CmPitch (force cases), inlet-outlet
                              pressure difference (T1, T2); running mean
                              and +-RMS band over the D-042 window W,
                              shaded final averaging window, vertical line
                              at iters_to_stationary
    hist_<run>_residuals.pdf  coupledFoam R, rU, rp and the simpleFoam
                              initial residuals of every solved field,
                              per iteration (left) and over wall time
                              (right), same y range for both solvers
    hist_<run>_wall.pdf       cumulative wall time per iteration with a
                              linear fit, time per iteration (raw and
                              rolling median), coupledFoam split
                              tAsm/tSolve/tTurb/other per iteration;
                              totals as wall-clock AND CPU-hours
    hist_<run>_linear.pdf     linear solver and controls per iteration:
                              coupledFoam linIters, linRes, eta, rho, GAMG
                              cycle type and nPostSweeps (autoTune events),
                              CFL, omega, CFL cuts and rollbacks, nLocLim,
                              nStat/nDyn, nClamped, turbulence bounding
                              events, every further numeric CF| field that
                              varies; simpleFoam linear iterations per field

and report/paper/figures/histories.tex, the appendix text (one subsection
per run, \\cffigure entries and generated key numbers), \\input by both
papers.

Sources (read-only): run/<run>/log.coupledFoam (all key=value fields of the
CF| lines are parsed generically, so new fields such as CdMean/CdRms appear
automatically), postProcessing/coupledFoam/summary.json,
postProcessing/forceCoeffs*/<t>/coefficient.dat, postProcessing/inletP and
outletP (surfaceFieldValue), run/ref_*/log.simpleFoam and reference.json.
A directory whose log.* files changed in the last BUSY_MINUTES is a live
run and is skipped (its part stays "pending"). If bench/diag_tools.py
exists (TASK 5 diagnostics) and the run has diagnostics/diag.rank*.jsonl,
its phase classification is drawn as background bands.

Staleness guard (TASK 6): a coupledFoam run directory is read only if its
provenance.json names the target commit (tests/cflib/provenance.py); a
stale run is treated like a missing one (pending). Called from
bench/make_report.py the guard of the report run is used; standalone the
target is --commit (default HEAD), --allow-stale disables the check.
simpleFoam reference directories are exempt (system OpenFOAM, cached).

Usage: ~/OF/venv/bin/python bench/plot_histories.py [--commit SHA]
           [--allow-stale] [run ...]
"""

from __future__ import annotations

import json
import re
import sys
import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

REPO = Path(__file__).resolve().parents[1]
RUN = REPO / "run"
FIG = REPO / "report" / "paper" / "figures"
sys.path.insert(0, str(REPO / "tests"))
sys.path.insert(0, str(REPO / "bench"))
from cflib import logs, provenance  # noqa: E402
import run_bench  # noqa: E402

BUSY_MINUTES = 10.0
# One colour pair everywhere (Okabe-Ito, colour-blind safe): simpleFoam
# vermillion, coupledFoam blue
C_NATIVE = "#D55E00"
C_COUPLED = "#0072B2"
# per-field series (residuals, linear iterations), Okabe-Ito order
SERIES = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9",
          "#000000"]
NMAX = 2500      # plotted points per curve (dense histories are binned)
SPLIT = {"tAsm": "#9ecae1", "tSolve": "#1f5fbf", "tTurb": "#fdae6b",
         "other": "#d9d9d9"}

# run name: coupledFoam dir, simpleFoam reference dir, title, monitor
# T5 on the coarse development mesh (CF_T5_MESH=coarse, D-059) unless only
# the fine-mesh run exists
T5_RUN = ("T5_np10" if (RUN / "T5_np10").is_dir()
          and not (RUN / "T5_coarse_np10").is_dir() else "T5_coarse_np10")
RUNS = {
    "T0_Re100_np1": ("T0_Re100_np1", "ref_T0_Re100", "T0 cavity, Re 100, 1 rank", None),
    "T0_Re100_np4": ("T0_Re100_np4", "ref_T0_Re100", "T0 cavity, Re 100, 4 ranks", None),
    "T0_Re1000_np1": ("T0_Re1000_np1", "ref_T0_Re1000", "T0 cavity, Re 1000, 1 rank", None),
    "T0_Re1000_np4": ("T0_Re1000_np4", "ref_T0_Re1000", "T0 cavity, Re 1000, 4 ranks", None),
    "T1_np1": ("T1_np1", "ref_T1", "T1 pitzDaily, 1 rank", "dp"),
    "T1_np4": ("T1_np4", "ref_T1", "T1 pitzDaily, 4 ranks", "dp"),
    "T2_np1": ("T2_np1", "ref_T2", "T2 backward-facing step, 1 rank", "dp"),
    "T2_np4": ("T2_np4", "ref_T2", "T2 backward-facing step, 4 ranks", "dp"),
    "T3_kOmegaSST_np1": ("T3_kOmegaSST_np1", "ref_T3_kOmegaSST", "T3 airFoil2D, SST", "forces"),
    "T3_GEKO_np1": ("T3_GEKO_np1", "ref_T3_GEKO", "T3 airFoil2D, GEKO", "forces"),
    "T4a_np10": ("T4a_np10", "ref_T4a_np10", "T4a motorBike, 354k cells, 10 ranks", "forces"),
    "T4b_np10": ("T4b_np10", "ref_T4b_np10", "T4b motorBike, 1.70M cells, 10 ranks", "forces"),
    T5_RUN: (T5_RUN, "ref_" + T5_RUN, "T5 Ahmed body, 10 ranks"
             + (" (coarse mesh)" if "coarse" in T5_RUN else ""), "forces"),
}
OSCILLATORY = ("T4", "T5")      # D-042 wake cases

plt.rcParams.update({
    "figure.dpi": 150, "savefig.dpi": 300, "font.size": 7,
    "axes.grid": True, "grid.alpha": 0.3, "legend.frameon": False,
    "pdf.fonttype": 42, "ps.fonttype": 42, "lines.linewidth": 0.8,
})


# --------------------------------------------------------------------------- #
# reading
# --------------------------------------------------------------------------- #

def busy(case: Path) -> bool:
    now = time.time()
    return any(now - f.stat().st_mtime < BUSY_MINUTES * 60
               for f in case.glob("log.*"))


_TUNE = re.compile(r"^GAMG-tune: (\S+)(?: (\w+)->(\w+))?")
_CYCLE0 = re.compile(r"^blockGAMG: levels \d+, mergeLevels \d+, C_op \S+, cycle (\w+)")
_NPOST = re.compile(r"^\s*nPostSweeps\s+(\d+);")
_RES = re.compile(r"^\w+: +Solving for (\w+), Initial residual = ([^,]+), "
                  r"Final residual = ([^,]+), No Iterations (\d+)")
_EXEC = re.compile(r"ExecutionTime = ([\d.eE+-]+) s +ClockTime = ([\d.eE+-]+) s")


def read_cf(case: Path) -> dict | None:
    """All CF| fields as arrays, GAMG-tune events, bounding counts, wall axis,
    summary."""
    log = case / "log.coupledFoam"
    if not log.exists():
        return None
    rows, events, bounds = [], [], []
    cyc0, npost0 = None, None
    nb = 0
    with open(log, errors="replace") as fh:
        for line in fh:
            if line.startswith("CF|"):
                r = {}
                for tok in line.split()[1:]:
                    if "=" in tok:
                        k, v = tok.split("=", 1)
                        try:
                            r[k] = float(v)
                        except ValueError:
                            pass
                if "iter" in r:
                    rows.append(r)
                    bounds.append(nb)
                    nb = 0
                continue
            if line.startswith("bounding "):
                nb += 1
                continue
            m = _TUNE.match(line)
            if m:
                it = (rows[-1]["iter"] + 1) if rows else 1
                events.append((it, m.group(1), m.group(2), m.group(3)))
                continue
            if cyc0 is None:
                m = _CYCLE0.match(line)
                if m:
                    cyc0 = m.group(1)
                    continue
            if npost0 is None:
                m = _NPOST.match(line)
                if m:
                    npost0 = int(m.group(1))
    if not rows:
        return None
    keys = []
    for r in rows:
        for k in r:
            if k not in keys:
                keys.append(k)
    d = {k: np.array([r.get(k, np.nan) for r in rows]) for k in keys}
    d["nBound"] = np.array(bounds, dtype=float) / 2.0   # min and max line
    summ = logs.coupled_summary(case)
    t = np.nancumsum(d["tIter"]) if "tIter" in d else np.arange(len(rows), dtype=float)
    if summ.get("wallSeconds") and t[-1] > 0:
        t = t * (summ["wallSeconds"] / t[-1])
    d["_t"] = t
    d["_summary"] = summ
    d["_events"] = events
    d["_cycle0"] = cyc0
    d["_npost0"] = npost0
    return d


def read_sf(case: Path) -> dict | None:
    log = case / "log.simpleFoam"
    if not log.exists():
        return None
    res: dict[str, list] = {}
    fin: dict[str, list] = {}
    nit: dict[str, list] = {}
    ex, clk, bnd = [], [], []
    cur_r, cur_f, cur_n, seen = {}, {}, {}, set()
    nb = 0
    started = False

    def flush():
        n = len(ex)
        for f in set(res) | set(cur_r):
            for store, cur in ((res, cur_r), (fin, cur_f), (nit, cur_n)):
                lst = store.setdefault(f, [np.nan] * (n - 1))
                lst.append(cur.get(f, np.nan))

    with open(log, errors="replace") as fh:
        for line in fh:
            if line.startswith("Time = "):
                started = True
                cur_r, cur_f, cur_n, seen = {}, {}, {}, set()
                nb = 0
                continue
            if not started:
                continue
            m = _RES.match(line)
            if m:
                f = m.group(1)
                if f not in seen:
                    cur_r[f] = float(m.group(2))
                    seen.add(f)
                cur_f[f] = float(m.group(3))
                cur_n[f] = cur_n.get(f, 0) + int(m.group(4))
                continue
            if line.startswith("bounding "):
                nb += 1
                continue
            m = _EXEC.search(line)
            if m:
                ex.append(float(m.group(1)))
                clk.append(float(m.group(2)))
                bnd.append(nb / 2.0)
                flush()
    if not ex:
        return None
    t = np.array(ex)
    if clk[-1] > 0 and t[-1] > 0:
        t = t * (clk[-1] / t[-1])
    n = len(t)
    ref = {}
    rj = case / "reference.json"
    if rj.exists():
        try:
            ref = json.loads(rj.read_text())
        except json.JSONDecodeError:
            ref = {}
    return {"n": n, "t": t, "res": {k: np.array(v[:n]) for k, v in res.items()},
            "fin": {k: np.array(v[:n]) for k, v in fin.items()},
            "nit": {k: np.array(v[:n]) for k, v in nit.items()},
            "nBound": np.array(bnd), "ref": ref}


def _dat(path: Path) -> np.ndarray | None:
    try:
        a = np.loadtxt(path, comments="#", usecols=None, dtype=str)
    except (ValueError, OSError):
        return None
    return a


def forces(case: Path) -> dict | None:
    """Cd, Cl, CmPitch per iteration, restarts concatenated (later start
    directories override earlier iterations)."""
    pp = case / "postProcessing"
    fos = sorted(pp.glob("forceCoeffs*")) if pp.is_dir() else []
    for fo in fos:
        out: dict[float, dict] = {}
        starts = sorted((d for d in fo.iterdir() if d.is_dir()),
                        key=lambda d: _float(d.name))
        for s in starts:
            f = s / "coefficient.dat"
            if not f.exists():
                continue
            hdr = None
            with open(f, errors="replace") as fh:
                for line in fh:
                    if line.startswith("# Time"):
                        hdr = line[1:].split()
                        continue
                    if line.startswith("#") or not line.strip() or hdr is None:
                        continue
                    v = line.split()
                    try:
                        out[float(v[0])] = {h: float(x) for h, x in zip(hdr, v)}
                    except ValueError:
                        continue
        if out:
            it = np.array(sorted(out))
            keys = [k for k in ("Cd", "Cl", "CmPitch") if k in out[it[0]]]
            return {"iter": it, **{k: np.array([out[i][k] for i in it])
                                   for k in keys}}
    return None


def _float(s: str) -> float:
    try:
        return float(s)
    except ValueError:
        return -1.0


def pressure_drop(case: Path) -> dict | None:
    """areaAverage(p) inlet - outlet per iteration (surfaceFieldValue)."""
    vals = {}
    for side in ("inletP", "outletP"):
        d = case / "postProcessing" / side
        if not d.is_dir():
            return None
        h: dict[float, float] = {}
        for s in sorted((x for x in d.iterdir() if x.is_dir()),
                        key=lambda x: _float(x.name)):
            f = s / "surfaceFieldValue.dat"
            if not f.exists():
                continue
            for line in open(f, errors="replace"):
                if line.startswith("#") or not line.strip():
                    continue
                v = line.split()
                try:
                    h[float(v[0])] = float(v[1])
                except (ValueError, IndexError):
                    continue
        vals[side] = h
    it = np.array(sorted(set(vals["inletP"]) & set(vals["outletP"])))
    if not len(it):
        return None
    return {"iter": it, "dp": np.array([vals["inletP"][i] - vals["outletP"][i]
                                        for i in it])}


def phases(case: Path):
    """[(start, end, phase)] from the TASK 5 diagnostics, if available."""
    if not (REPO / "bench" / "diag_tools.py").exists():
        return [], {}
    if not any((case / "diagnostics").glob("diag.rank*.jsonl")):
        return [], {}
    try:
        import diag_tools  # noqa: PLC0415
        df = diag_tools.load(case)
        col = getattr(diag_tools, "PHASE_COLOR", {})
    except Exception:  # noqa: BLE001 - optional input, never fatal
        return [], {}
    if "phase" not in df:
        return [], col
    it = np.asarray(df.index, dtype=float)
    ph = df["phase"].astype(str).to_numpy()
    bands, s = [], 0
    for i in range(1, len(ph) + 1):
        if i == len(ph) or ph[i] != ph[s]:
            bands.append((it[s] - 0.5, it[i - 1] + 0.5, ph[s]))
            s = i
    return bands, col


# --------------------------------------------------------------------------- #
# plotting helpers
# --------------------------------------------------------------------------- #

def roll(y: np.ndarray, w: int, fn=np.nanmedian) -> np.ndarray:
    """Trailing rolling statistic (window w, shorter at the start)."""
    y = np.asarray(y, dtype=float)
    out = np.full_like(y, np.nan)
    if w <= 1:
        return y.copy()
    if fn is np.nanmedian and len(y) > w:
        from numpy.lib.stride_tricks import sliding_window_view  # noqa: PLC0415
        pad = np.concatenate([np.full(w - 1, np.nan), y])
        step = max(1, len(y) // 4000)       # median every `step` iterations
        idx = np.arange(0, len(y), step)
        win = sliding_window_view(pad, w)[idx]
        with np.errstate(all="ignore"):
            import warnings  # noqa: PLC0415
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                med = np.nanmedian(win, axis=1)
        return np.interp(np.arange(len(y)), idx, med)
    if fn is np.nanmean:
        c = np.nancumsum(np.nan_to_num(y))
        cnt = np.cumsum(np.isfinite(y))
        i = np.arange(len(y))
        j = np.maximum(0, i - w + 1)
        cs = np.concatenate([[0.0], c])
        cn = np.concatenate([[0], cnt])
        k = cn[i + 1] - cn[j]
        with np.errstate(all="ignore"):
            out = np.where(k > 0, (cs[i + 1] - cs[j]) / np.maximum(k, 1), np.nan)
        return out
    for i in range(len(y)):
        out[i] = fn(y[max(0, i - w + 1):i + 1])
    return out


def roll_std(y: np.ndarray, w: int) -> np.ndarray:
    m = roll(y, w, np.nanmean)
    m2 = roll(np.asarray(y) ** 2, w, np.nanmean)
    return np.sqrt(np.clip(m2 - m ** 2, 0, None))


def bands(ax, ph):
    b, col = ph
    for x0, x1, name in b:
        ax.axvspan(x0, x1, color=col.get(name, "#eeeeee"), alpha=0.5, lw=0,
                   zorder=0)


def events(ax, ev):
    for it, what, a, b in ev:
        ax.axvline(it, color="k", lw=0.4, ls=":", zorder=1)


def small_legend(ax, **kw):
    h, _ = ax.get_legend_handles_labels()
    if h:
        ax.legend(fontsize=5.5, ncol=kw.pop("ncol", 3), **kw)


def iter_axis(ax, lengths: list[int]) -> None:
    """Log iteration axis when the run lengths differ by more than 3x (the
    shorter coupledFoam history would otherwise be squeezed into the left
    edge)."""
    lengths = [n for n in lengths if n]
    if len(lengths) > 1 and max(lengths) > 3 * min(lengths):
        ax.set_xscale("log")
        ax.set_xlim(0.9, max(lengths) * 1.1)


def _bins(n: int, nmax: int = NMAX) -> list[slice]:
    edges = np.linspace(0, n, min(n, nmax) + 1).astype(int)
    return [slice(a, b) for a, b in zip(edges[:-1], edges[1:]) if b > a]


def decimate(x, y, nmax: int = NMAX):
    """Smooth curve with at most nmax points (bin means); exact below."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    if len(x) <= nmax:
        return x, y
    with np.errstate(all="ignore"):
        bx = np.array([np.nanmean(x[b]) for b in _bins(len(x), nmax)])
        by = np.array([np.nanmean(y[b]) if np.any(np.isfinite(y[b])) else np.nan
                       for b in _bins(len(x), nmax)])
    return bx, by


NHEAD = 400      # leading iterations always drawn exactly (log-x axes)


def raw(ax, x, y, color, label=None, log=False, lw=0.7, alpha=0.5):
    """Raw per-iteration data drawn light but visible: a line when short;
    for long histories the first NHEAD iterations exactly, then the per-bin
    min/max envelope plus the bin-mean line (vector, at most NMAX bins) -
    never rasterized, so nothing can disappear in the PDF."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    # simpleFoam is usually the smoother history: draw it above the noisier
    # coupledFoam raw data so it is never hidden
    z = 2.6 if color == C_NATIVE else 2.0
    if len(x) <= NMAX:
        ax.plot(x, y, color=color, lw=lw, alpha=alpha, label=label, zorder=z)
        return
    ax.plot(x[:NHEAD], y[:NHEAD], color=color, lw=lw, alpha=alpha, zorder=z)
    x, y = x[NHEAD - 1:], y[NHEAD - 1:]
    bx, lo, hi = [], [], []
    for b in _bins(len(x), NMAX):
        yy = y[b][np.isfinite(y[b])]
        if log:
            yy = yy[yy > 0]
        if not yy.size:
            continue
        bx.append(np.nanmean(x[b]))
        lo.append(yy.min())
        hi.append(yy.max())
    ax.fill_between(bx, lo, hi, color=color, alpha=alpha * 0.6, lw=0,
                    label=label, zorder=z)
    bm = decimate(x, np.where(y > 0, y, np.nan) if log else y)[1]
    ax.plot(decimate(x, x)[0], bm, color=color, lw=lw, alpha=min(1.0, alpha + 0.2),
            zorder=z)


def main(ax, x, y, color, label=None, lw=1.4, ls="-", zorder=4):
    """Main curve (mean, median, cumulative, fit): strong, on top."""
    bx, by = decimate(x, y)
    ax.plot(bx, by, color=color, lw=lw, ls=ls, label=label, zorder=zorder)


def save(fig, stem: str, outdir: Path) -> str:
    fig.savefig(outdir / f"{stem}.pdf", bbox_inches="tight")
    plt.close(fig)
    return stem


# --------------------------------------------------------------------------- #
# figures
# --------------------------------------------------------------------------- #

def load_series(cf, sf, cfd, sfd, monitor, info):
    """[(solver, colour, iterations, wall time at each, {quantity: values})]"""
    series = []
    for solver, d, case, col, taxis in (
            ("simpleFoam", sf, sfd, C_NATIVE, sf["t"] if sf else None),
            ("coupledFoam", cf, cfd, C_COUPLED, cf["_t"] if cf else None)):
        if d is None:
            continue
        h = forces(case) if monitor == "forces" else pressure_drop(case)
        if h is None:
            continue
        it = h["iter"]
        qs = {k: v for k, v in h.items() if k != "iter"}
        # wall time at each monitored iteration (monitor iteration i is the
        # i-th solver iteration of the run); iterations beyond the logged
        # ones (a continuation whose log is elsewhere) have no wall time
        idx = np.clip(it.astype(int) - 1, 0, len(taxis) - 1)
        tt = np.where(it <= len(taxis), taxis[idx], np.nan)
        if np.any(it > len(taxis)) and info is not None:
            info["notes"].append(
                f"{solver}: the monitor has {int(it.max())} iterations, the "
                f"solver log {len(taxis)}; the later iterations have no "
                "wall-clock time (right column)")
        series.append((solver, col, it, tt, qs))
    return series


def _user_iteration(run: str, solver: str) -> int | None:
    """User-judged convergence iteration of a test run (D-060)."""
    try:
        import user_convergence as ucv  # noqa: PLC0415
        return ucv.lookup(ucv.load(), ucv.case_key(run), solver)
    except Exception:  # noqa: BLE001 - the plot must not fail on the file
        return None


def fig_loads(run, title, cf, sf, cfd, sfd, monitor, ph, outdir, info):
    series = load_series(cf, sf, cfd, sfd, monitor, info)
    extra = ({k: cf[k] for k in cf if re.match(r"^(Cd|Cl|Cm)\w*(Mean|Rms)$", k)}
             if cf is not None else {})
    if not series:
        info["notes"].append("no monitored quantity history")
        return None
    quants = []
    for s_ in series:
        for q in s_[4]:
            if q not in quants:
                quants.append(q)
    labels = {"Cd": "$C_d$", "Cl": "$C_l$", "CmPitch": "$C_{m}$ (pitch)",
              "dp": r"$\Delta p$ inlet$-$outlet [m$^2$/s$^2$]"}
    fig, axs = plt.subplots(len(quants), 2, figsize=(6.5, 1.9 * len(quants) + 0.5),
                            squeeze=False)
    for r, q in enumerate(quants):
        for solver, col, it, tt, qs in series:
            if q not in qs:
                continue
            y = qs[q]
            n = len(y)
            W = run_bench.stat_window(n)
            m = roll(y, W, np.nanmean)
            sd = roll_std(y, W)
            m[:max(W - 1, 0)] = np.nan   # only full windows (D-042 W)
            sd[:max(W - 1, 0)] = np.nan
            for c, x in ((0, it), (1, tt)):
                ax = axs[r, c]
                ok = np.isfinite(x)
                raw(ax, x[ok], y[ok], col, label=f"{solver} per iteration" if r == 0 and c == 0 else None)
                bx, bm = decimate(x[ok], m[ok])
                _, bs = decimate(x[ok], sd[ok])
                ax.fill_between(bx, bm - bs, bm + bs, color=col, alpha=0.18, lw=0,
                                zorder=3)
                ax.plot(bx, bm, color=col, lw=1.6, ls="--", zorder=5,
                        label=f"{solver} running mean $\\pm$RMS (W={W})"
                        if r == 0 and c == 0 else None)
                if W >= 2 and n >= W and np.isfinite(x[n - W]):
                    ax.axvspan(x[n - W], np.nanmax(x), color=col, alpha=0.08,
                               lw=0, zorder=1)
            if q in ("Cd", "Cl") and monitor == "forces":
                hist = {k: list(v) for k, v in qs.items() if k in ("Cd", "Cl")}
                try:
                    its = run_bench.iters_to_stationary(hist)
                    ev = run_bench.stationary_eval(hist)
                except Exception:  # noqa: BLE001
                    its, ev = None, {}
                if its:
                    axs[r, 0].axvline(its, color=col, lw=1.2, ls="-.", zorder=6)
                    k = min(int(np.searchsorted(it, its)), len(tt) - 1)
                    if np.isfinite(tt[k]):
                        axs[r, 1].axvline(tt[k], color=col, lw=1.2, ls="-.",
                                          zorder=6)
                # user-judged convergence iteration (D-060): solid line
                uit = _user_iteration(run, solver)
                if uit:
                    axs[r, 0].axvline(uit, color=col, lw=1.8, ls="-", zorder=7,
                                      label=f"{solver} user convergence"
                                      if r == 0 else None)
                    k = min(int(np.searchsorted(it, uit)), len(tt) - 1)
                    if np.isfinite(tt[k]):
                        axs[r, 1].axvline(tt[k], color=col, lw=1.8, ls="-",
                                          zorder=7)
                    info.setdefault("user", {})[solver] = uit
                info.setdefault("stat", {})[solver] = (its, ev.get("W"),
                                                       ev.get("Cd_mean"),
                                                       ev.get("Cl_mean"))
            else:
                w = run_bench.stat_window(n)
                win = y[n - w:] if w else y
                info.setdefault("mon", {})[(solver, q)] = (float(np.mean(win)),
                                                           float(np.std(win)))
        k = f"{q}Mean"
        if k in extra:
            axs[r, 0].plot(cf["iter"], extra[k], color="k", lw=1.0, zorder=6,
                           label=f"coupledFoam {k} (solver)")
        axs[r, 0].set_ylabel(labels.get(q, q))
        # zoom: ignore the start-up transient in the y range
        vals = np.concatenate([s_[4][q][len(s_[4][q]) // 5:] for s_ in series
                               if q in s_[4] and len(s_[4][q]) > 5] or [np.array([0.0])])
        vals = vals[np.isfinite(vals)]
        if vals.size > 1:
            lo, hi = np.percentile(vals, [0.5, 99.5])
            pad = 0.25 * (hi - lo) + 1e-12
            for c in (0, 1):
                axs[r, c].set_ylim(lo - pad, hi + pad)
        bands(axs[r, 0], ph)
        iter_axis(axs[r, 0], [len(s_[2]) for s_ in series])
    axs[-1, 0].set_xlabel("outer iteration")
    axs[-1, 1].set_xlabel("wall-clock time [s] (fair comparison)")
    h, lab = axs[0, 0].get_legend_handles_labels()
    fig.suptitle(f"{title}: monitored quantities (dash-dot: iterations to "
                 "stationary window, D-042; shading: final averaging window)",
                 fontsize=7)
    fig.tight_layout(rect=(0, 0.06, 1, 0.97))
    fig.legend(h, lab, loc="lower center", ncol=2, fontsize=6)
    return save(fig, f"hist_{run}_loads", outdir)


def fig_residuals(run, title, cf, sf, ph, outdir):
    rows = [x for x in (("coupledFoam", cf), ("simpleFoam", sf)) if x[1] is not None]
    if not rows:
        return None
    fig, axs = plt.subplots(len(rows), 2, figsize=(6.5, 2.2 * len(rows) + 0.3),
                            squeeze=False)
    ylo, yhi = [], []
    for r, (solver, d) in enumerate(rows):
        if solver == "coupledFoam":
            curves = [(k, d[k]) for k in ("R", "rU", "rp") if k in d]
            x0, t = d["iter"], d["_t"]
        else:
            order = ["p", "Ux", "Uy", "Uz", "k", "omega", "epsilon", "nuTilda"]
            fields = sorted(d["res"], key=lambda f: order.index(f) if f in order else 99)
            curves = [(f, d["res"][f]) for f in fields]
            x0, t = np.arange(1, d["n"] + 1), d["t"]
        for j, (lab, y) in enumerate(curves):
            col = SERIES[j % len(SERIES)]
            for c, x in ((0, x0), (1, t)):
                x = x[:len(y)]
                ax = axs[r, c]
                if len(x) > NMAX:
                    raw(ax, x, y, col, log=True, alpha=0.3)
                    bx = decimate(x, x)[0]
                    by = np.exp(decimate(x, np.log(np.where(y > 0, y, np.nan)))[1])
                    ax.plot(bx, by, color=col, lw=1.2, zorder=4,
                            label=lab if c == 0 else None)
                else:
                    ax.plot(x, y, color=col, lw=1.2, zorder=4,
                            label=lab if c == 0 else None)
                ax.set_yscale("log")
            fy = y[np.isfinite(y) & (y > 0)]
            if fy.size:
                ylo.append(fy.min())
                yhi.append(fy.max())
        axs[r, 0].set_ylabel(f"{solver}\nresidual",
                             color=C_COUPLED if solver == "coupledFoam" else C_NATIVE)
        small_legend(axs[r, 0], ncol=4, loc="lower left")
        bands(axs[r, 0], ph if solver == "coupledFoam" else ([], {}))
        if solver == "coupledFoam":
            events(axs[r, 0], d["_events"])
    if ylo:
        for ax in axs.flat:
            ax.set_ylim(min(ylo) / 2, max(yhi) * 2)
    ns = [len(d["iter"]) if s_ == "coupledFoam" else d["n"] for s_, d in rows]
    for r in range(len(rows)):
        iter_axis(axs[r, 0], ns)
    # both solvers on the same wall-clock axis
    tmax = max((d["_t"][-1] if s_ == "coupledFoam" else d["t"][-1]) for s_, d in rows)
    for r in range(len(rows)):
        axs[r, 1].set_xlim(0, tmax * 1.02)
    axs[-1, 0].set_xlabel("outer iteration")
    axs[-1, 1].set_xlabel("wall-clock time [s] (same axis for both solvers)")
    fig.suptitle(f"{title}: residuals (coupledFoam: combined $R$, momentum "
                 "$r_U$, continuity $r_p$; simpleFoam: initial residuals; "
                 "long runs: band = per-bin min/max, line = bin mean). "
                 "Same y range; normalisations differ.", fontsize=6.5)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    return save(fig, f"hist_{run}_residuals", outdir)


def _timing(solver, d):
    if solver == "coupledFoam":
        t = d["_t"]
        it = d["iter"]
        per = (d["tIter"] * (t[-1] / np.nansum(d["tIter"])) if "tIter" in d
               else np.diff(t, prepend=0))
        summ = d["_summary"]
        wall, cpuh = summ.get("wallSeconds"), summ.get("cpuHours")
    else:
        t = d["t"]
        it = np.arange(1, d["n"] + 1, dtype=float)
        per = np.diff(t, prepend=0.0)
        ref = d["ref"]
        ta = ref.get("timingAllrun") or {}
        wall = (ref.get("wallSeconds") or ref.get("wallSecondsSolver")
                or ta.get("wallSeconds") or t[-1])
        cpuh = ref.get("cpuHours") or ta.get("cpuHours")
        if cpuh is None and wall:
            cpuh = wall * (ref.get("nProcs") or 1) / 3600.0
    return t, it, per, wall, cpuh


def fig_wall(run, title, cf, sf, ph, outdir, info):
    have = [x for x in (("simpleFoam", sf, C_NATIVE), ("coupledFoam", cf, C_COUPLED))
            if x[1] is not None]
    if not have:
        return None
    lens = {}
    for solver, d, col in have:
        lens[solver] = len(_timing(solver, d)[0])
    zoom = (cf is not None and sf is not None
            and lens["simpleFoam"] > 3 * lens["coupledFoam"])
    ncol = 2 if zoom else 1
    fig = plt.figure(figsize=(6.5, 6.8 if cf is not None else 4.6))
    nrow = 3 if cf is not None else 2
    gs = fig.add_gridspec(nrow, ncol)
    top = [fig.add_subplot(gs[0, c]) for c in range(ncol)]
    mid = [fig.add_subplot(gs[1, c]) for c in range(ncol)]
    for solver, d, col in have:
        t, it, per, wall, cpuh = _timing(solver, d)
        n = len(t)
        A = np.vstack([it, np.ones_like(it)]).T
        slope, icpt = np.linalg.lstsq(A, t, rcond=None)[0]
        resid = t - (slope * it + icpt)
        r2 = 1 - np.sum(resid ** 2) / max(np.sum((t - t.mean()) ** 2), 1e-30)
        q = max(1, n // 10)
        first, last = np.nanmedian(per[:q]), np.nanmedian(per[-q:])
        info.setdefault("wall", {})[solver] = dict(
            n=n, wall=wall, cpuh=cpuh, slope=slope, r2=r2, first=first,
            last=last)
        lab = (f"{solver}: {n} it., {wall:.4g} s, {cpuh:.3g} CPU-h"
               if wall and cpuh else solver)
        med = roll(per, max(5, n // 50))
        for c in range(ncol):
            main(top[c], it, t, col, label=lab if c == 0 else None, lw=1.8)
            main(top[c], it, slope * it + icpt, col, lw=1.0,
                 ls=":", zorder=5,
                 label=(f"{solver} linear fit {slope:.3g} s/it ($R^2$ {r2:.3f})"
                        if c == 0 else None))
            raw(mid[c], it, per, col, log=True, alpha=0.3)
            main(mid[c], it, med, col, lw=1.6,
                 label=f"{solver} rolling median" if c == 0 else None)
            mid[c].set_yscale("log")
    if zoom:
        nz = lens["coupledFoam"] * 1.05
        for ax in (top[1], mid[1]):
            ax.set_xlim(0, nz)
            ax.set_title(f"zoom: the {lens['coupledFoam']} coupledFoam iterations",
                         fontsize=6)
        # y range of the zoom: what both solvers reach within nz iterations
        ymax = 0.0
        for s_, d, _ in have:
            t, it = _timing(s_, d)[:2]
            ymax = max(ymax, float(np.nanmax(t[it <= nz])) if np.any(it <= nz) else 0.0)
        top[1].set_ylim(0, ymax * 1.05)
        top[0].set_title("all iterations (linear axis: linear growth looks linear)",
                         fontsize=6)
    top[0].set_ylabel("cumulative wall time [s]")
    mid[0].set_ylabel("wall time per iteration [s]")
    small_legend(top[0], ncol=1, loc="upper left")
    small_legend(mid[0], ncol=1, loc="upper right")
    axes = top + mid
    from matplotlib.ticker import MaxNLocator  # noqa: PLC0415
    for ax in axes:
        ax.xaxis.set_major_locator(MaxNLocator(5, integer=True))
    if cf is not None:
        ax3 = fig.add_subplot(gs[2, :])
        axes.append(ax3)
        it = cf["iter"]
        n = len(it)
        w = max(5, n // 50)
        parts = {k: roll(cf[k], w) for k in ("tAsm", "tSolve", "tTurb") if k in cf}
        if "tIter" in cf and parts:
            parts["other"] = np.clip(roll(cf["tIter"], w) - sum(parts.values()), 0, None)
            ax3.stackplot(it, *parts.values(), labels=list(parts),
                          colors=[SPLIT[k] for k in parts], lw=0)
            ax3.set_ylabel("coupledFoam time per\niteration [s] (rolling median)")
            small_legend(ax3, ncol=4, loc="upper left")
            events(ax3, cf["_events"])
        bands(ax3, ph)
    for ax in axes:
        ax.set_xlabel("outer iteration")
    fig.suptitle(f"{title}: wall-clock time (dotted: linear fit; light band: "
                 "per-iteration spread; vertical dotted lines: GAMG autoTune "
                 "events)", fontsize=7)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    return save(fig, f"hist_{run}_wall", outdir)


def fig_overview(group: str, runs: list[str], outdir: Path, log=print) -> str | None:
    """Main-text overview: Cd and Cl of both solvers over wall-clock time
    for a group of force cases, with the final D-042 window means."""
    rows = []
    for name in runs:
        cfn, sfn, title, monitor = RUNS[name]
        cfd, sfd = RUN / cfn, RUN / sfn
        cf = (read_cf(cfd) if cfd.is_dir() and not busy(cfd)
              and provenance.active().check_run(cfd) else None)
        sf = read_sf(sfd) if sfd.is_dir() and not busy(sfd) else None
        rows.append((title, load_series(cf, sf, cfd, sfd, "forces", None)))
    if not any(r[1] for r in rows):
        return None
    fig, axs = plt.subplots(len(rows), 2, figsize=(6.5, 1.75 * len(rows) + 0.5),
                            squeeze=False)
    for r, (title, series) in enumerate(rows):
        for c, q in enumerate(("Cd", "Cl")):
            ax = axs[r, c]
            if not series:
                ax.text(0.5, 0.5, "results pending", ha="center", va="center",
                        transform=ax.transAxes, color="#555555")
                ax.set_yticks([])
            vals = []
            for solver, col, it, tt, qs in series:
                if q not in qs:
                    continue
                y = qs[q]
                ok = np.isfinite(tt)
                raw(ax, tt[ok], y[ok], col, alpha=0.3)
                n = len(y)
                W = run_bench.stat_window(n)
                m = roll(y, W, np.nanmean)
                m[:max(W - 1, 0)] = np.nan
                main(ax, tt[ok], m[ok], col, lw=1.6, ls="--")
                mw = float(np.mean(y[n - W:])) if W and n >= W else np.nan
                if np.isfinite(mw):
                    ax.axhline(mw, color=col, lw=0.8, ls=":", zorder=3)
                vals.append(y[n // 5:])
            if vals:
                v = np.concatenate(vals)
                v = v[np.isfinite(v)]
                if v.size > 1:
                    lo, hi = np.percentile(v, [0.5, 99.5])
                    pad = 0.3 * (hi - lo) + 1e-9
                    ax.set_ylim(lo - pad, hi + pad)
            ax.set_ylabel(f"{title}\n" + ("$C_d$" if q == "Cd" else "$C_l$"),
                          fontsize=6)
            ax.tick_params(labelsize=6)
    for c in range(2):
        axs[-1, c].set_xlabel("wall-clock time [s]")
    from matplotlib.lines import Line2D  # noqa: PLC0415
    h = [Line2D([], [], color=C_NATIVE, lw=1.6), Line2D([], [], color=C_COUPLED, lw=1.6),
         Line2D([], [], color="k", lw=1.6, ls="--"), Line2D([], [], color="k", lw=0.8, ls=":")]
    fig.legend(h, ["simpleFoam", "coupledFoam", "running mean over the D-042 window",
                   "final window mean"], loc="lower center", ncol=4, fontsize=6)
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    stem = f"loads_overview_{group}"
    save(fig, stem, outdir)
    log(f"overview {group}: written")
    return stem


KNOWN = {"iter", "CFL", "omega", "cuts", "R", "rU", "rp", "linIters",
         "linRes", "tAsm", "tSolve", "tTurb", "tIter", "tWall", "nStat",
         "nDyn", "nLocLim", "nRollback", "nClamped", "eta", "rho", "nBound",
         "aa"}


def fig_linear(run, title, cf, sf, ph, outdir, info):
    if cf is None and sf is None:
        return None
    panels = []
    if cf is not None:
        it = cf["iter"]
        panels += [
            ("linear iterations\nper outer it.", [("linIters", cf.get("linIters"))], "lin", "step"),
            ("linear residual\n(final, rel.)", [("linRes", cf.get("linRes"))], "log", "line"),
            ("E-W forcing $\\eta$", [("eta", cf.get("eta"))], "log", "line"),
            ("precond. $\\rho$", [("rho", cf.get("rho"))], "lin", "line"),
            ("GAMG cycle", "cycle", None, None),
            ("nPostSweeps", "npost", None, None),
            ("CFL", [("CFL", cf.get("CFL"))], "log", "line"),
            ("line search $\\omega$", [("omega", cf.get("omega"))], "lin", "line"),
            ("CFL cuts, rollbacks\n(cumulative)", [
                ("cuts", np.nancumsum(cf["cuts"]) if "cuts" in cf else None),
                ("nRollback", cf.get("nRollback"))], "lin", "step"),
            ("cells", [("nLocLim", cf.get("nLocLim")), ("nStat", cf.get("nStat")),
                       ("nDyn", cf.get("nDyn")), ("nClamped", cf.get("nClamped"))],
             "symlog", "step"),
            ("turbulence bounding\nevents per it.", [("bounding events", cf.get("nBound"))], "lin", "step"),
        ]
        if "aa" in cf and np.nanmax(np.abs(cf["aa"])) > 0:
            panels.append(("Anderson", [("aa", cf["aa"])], "lin", "step"))
        for k in cf:
            if k.startswith("_") or k in KNOWN or re.match(r"^(Cd|Cl|Cm)", k):
                continue
            v = cf[k]
            if isinstance(v, np.ndarray) and v.size and np.nanmax(v) != np.nanmin(v):
                panels.append((k, [(k, v)], "lin", "line"))
    # panels whose every series is constant are listed in the text instead
    keep = []
    for p in panels:
        if isinstance(p[1], list):
            ser = [(nm, y) for nm, y in p[1] if y is not None and np.any(np.isfinite(y))]
            if ser and all(np.nanmax(y) == np.nanmin(y) for _, y in ser):
                info.setdefault("constant", []).extend(
                    f"{nm} = {np.nanmax(y):g}" for nm, y in ser)
                continue
        keep.append(p)
    panels = keep
    if sf is not None:
        panels.append(("simpleFoam linear\niterations per field", "sfnit", "log", None))
    fig, axs = plt.subplots(len(panels), 1, figsize=(6.5, 0.95 * len(panels) + 0.6),
                            squeeze=False)
    axs = axs[:, 0]
    # coupledFoam panels share the iteration axis; the simpleFoam panel has
    # its own (it may run many more iterations)
    for ax, p in zip(axs[1:], panels[1:]):
        if cf is not None and p[1] != "sfnit":
            ax.sharex(axs[0])
    for ax, (lab, data, scale, kind) in zip(axs, panels):
        if data == "cycle":
            seq = {"V": 0, "F": 1, "W": 2, "K": 3}
            cur = cf["_cycle0"]
            y = np.full(len(it), np.nan)
            evs = sorted((e for e in cf["_events"] if e[1] == "cycleType"),
                         key=lambda e: e[0])
            j = 0
            for i, x in enumerate(it):
                while j < len(evs) and evs[j][0] <= x:
                    cur = evs[j][3]
                    j += 1
                y[i] = seq.get(cur, np.nan) if cur else np.nan
            ax.step(it, y, where="post", color=C_COUPLED)
            ax.set_yticks(list(seq.values()))
            ax.set_yticklabels(list(seq))
            ax.set_ylim(-0.5, 3.5)
            if cur is None:
                ax.text(0.5, 0.5, "cycle type not in the log", transform=ax.transAxes,
                        ha="center", fontsize=6)
            for e in cf["_events"]:
                if e[1] in ("failure", "kcycle-saturated"):
                    ax.axvline(e[0], color="#d62728", lw=0.6, ls="--")
                    ax.text(e[0], 3.3, e[1], fontsize=5, color="#d62728",
                            rotation=0, va="top")
        elif data == "npost":
            cur = cf["_npost0"]
            y = np.full(len(it), np.nan)
            evs = sorted((e for e in cf["_events"] if e[1] == "nPostSweeps"),
                         key=lambda e: e[0])
            if cur is None and evs:
                cur = int(evs[0][2])
            j = 0
            for i, x in enumerate(it):
                while j < len(evs) and evs[j][0] <= x:
                    cur = int(evs[j][3])
                    j += 1
                y[i] = cur if cur is not None else np.nan
            ax.step(it, y, where="post", color=C_COUPLED)
            if np.all(np.isnan(y)):
                ax.text(0.5, 0.5, "nPostSweeps not in the log", transform=ax.transAxes,
                        ha="center", fontsize=6)
            else:
                ax.set_ylim(np.nanmin(y) - 0.5, np.nanmax(y) + 0.5)
        elif data == "sfnit":
            xs = np.arange(1, sf["n"] + 1)
            order = ["p", "Ux", "Uy", "Uz", "k", "omega", "epsilon", "nuTilda"]
            for j, f in enumerate(sorted(sf["nit"], key=lambda f: order.index(f) if f in order else 99)):
                bx, by = decimate(xs, np.maximum(sf["nit"][f], 0.8))
                ax.plot(bx, by, color=SERIES[j % len(SERIES)], lw=1.0, label=f)
            allv = np.concatenate([v[np.isfinite(v)] for v in sf["nit"].values()] or [np.array([1.0])])
            if allv.size and allv.max() > 10 * max(allv.min(), 1):
                ax.set_yscale("log")
            small_legend(ax, ncol=6, loc="upper right")
            ax.text(0.01, 0.85, "simpleFoam (x axis: its own iterations)",
                    transform=ax.transAxes, fontsize=5)
            ax.set_xlabel("simpleFoam iteration", fontsize=6)
        else:
            nplot = 0
            for j, (name, y) in enumerate(data):
                if y is None:
                    continue
                c = SERIES[j % len(SERIES)] if len(data) > 1 else C_COUPLED
                if kind == "step":
                    ax.step(it, y, where="mid", color=c, lw=0.7, label=name)
                else:
                    ax.plot(it, y, color=c, lw=0.7, label=name)
                nplot += 1
            ys = [y for _, y in data if y is not None]
            fin = np.concatenate([y[np.isfinite(y) & (y > 0)] for y in ys]) if ys else np.array([])
            if scale == "log" and fin.size and fin.max() > 10 * fin.min():
                ax.set_yscale("log")
            elif scale == "symlog":
                ax.set_yscale("symlog", linthresh=1)
            if nplot > 1:
                small_legend(ax, ncol=4, loc="upper right")
            if nplot == 0:
                ax.text(0.5, 0.5, "not in this log", transform=ax.transAxes,
                        ha="center", fontsize=6)
        ax.set_ylabel(lab, fontsize=6)
        ax.tick_params(labelsize=6)
        if cf is not None and data != "sfnit":
            bands(ax, ph)
            events(ax, cf["_events"])
    if cf is not None:
        last_cf = max(i for i, p in enumerate(panels) if p[1] != "sfnit")
        axs[last_cf].set_xlabel("outer iteration (coupledFoam)", fontsize=6)
        info["linTotal"] = float(np.nansum(cf.get("linIters", np.array([np.nan]))))
    fig.suptitle(f"{title}: linear solver and controls per outer iteration "
                 "(dotted: GAMG autoTune events)", fontsize=7)
    fig.tight_layout(h_pad=0.2, rect=(0, 0, 1, 0.985))
    return save(fig, f"hist_{run}_linear", outdir)


# --------------------------------------------------------------------------- #
# appendix text
# --------------------------------------------------------------------------- #

def _tex(s: str) -> str:
    return s.replace("_", r"\_").replace("%", r"\%")


def _fmt(v, f="{:.3g}"):
    return "n/a" if v is None or (isinstance(v, float) and not np.isfinite(v)) else f.format(v)


def tex_section(run, title, made, info, status) -> str:
    L = [rf"\subsection{{{_tex(title)}}}\label{{app:hist:{run}}}", ""]
    s = []
    w = info.get("wall", {})
    for solver in ("coupledFoam", "simpleFoam"):
        if solver in w:
            x = w[solver]
            s.append(f"{solver}: {x['n']} iterations, "
                     f"{_fmt(x['wall'], '{:.4g}')}\\,s wall-clock, "
                     f"{_fmt(x['cpuh'])} CPU-h; time per iteration "
                     f"{_fmt(x['first'])}\\,s (first 10\\,\\%) and "
                     f"{_fmt(x['last'])}\\,s (last 10\\,\\%), linear fit "
                     f"{_fmt(x['slope'])}\\,s/it ($R^2={_fmt(x['r2'], '{:.3f}')}$)")
        elif status.get(solver):
            s.append(f"{solver}: {status[solver]}")
    if info.get("linTotal") is not None:
        s.append(f"coupledFoam linear iterations in total: {info['linTotal']:.0f}")
    for solver, v in (info.get("stat") or {}).items():
        its, W, cd, cl = v
        s.append(f"{solver} D-042 window $W={_fmt(W, '{}')}$: "
                 f"$\\bar C_d={_fmt(cd, '{:.4f}')}$, $\\bar C_l={_fmt(cl, '{:.4f}')}$, "
                 f"stationary from iteration {_fmt(its, '{}')}")
    for (solver, q), (m, sd) in (info.get("mon") or {}).items():
        if q == "dp":
            s.append(f"{solver} window mean of $\\Delta p$: {m:.5g} "
                     f"$\\pm$ {sd:.2g}\\,m$^2$/s$^2$")
    if info.get("constant"):
        s.append("constant over the whole coupledFoam run (panel omitted): "
                 + ", ".join(sorted(set(info["constant"]))))
    for n_ in info.get("notes", []):
        s.append(n_)
    if s:
        L += [r"\begin{itemize}\setlength{\itemsep}{0pt}\small"]
        L += [rf"\item {x}." for x in s]
        L += [r"\end{itemize}", ""]
    caps = {
        "loads": "monitored quantities per outer iteration (left) and over "
                 "wall-clock time (right); thin: raw, dashed: running mean over "
                 "the D-042 window $W$ with $\\pm$RMS band, shaded: final "
                 "averaging window, dash-dot: iterations to a stationary window",
        "residuals": "residuals per outer iteration (left) and over wall-clock "
                     "time (right), coupledFoam top, simpleFoam bottom, same "
                     "y range",
        "wall": "cumulative wall time with linear fit (top), wall time per "
                "iteration with rolling median (middle), coupledFoam time "
                "split into assembly, linear solve, turbulence and other "
                "(bottom)",
        "linear": "linear solver and pseudo-time controls per outer "
                  "iteration; bottom: simpleFoam linear iterations per field",
    }
    for kind in ("loads", "residuals", "wall", "linear"):
        if kind == "loads" and kind not in made and info.get("noLoads"):
            continue
        L.append(rf"\cfhistfigure{{hist_{run}_{kind}}}{{{_tex(title)}: {caps[kind]}.}}"
                 rf"{{fig:hist:{run}:{kind}}}")
    L += [r"\clearpage", ""]
    return "\n".join(L)


# --------------------------------------------------------------------------- #
# driver
# --------------------------------------------------------------------------- #

# \cfhistfigure: \cffigure with [H] placement (package float), so every
# figure stays below its subsection heading
HIST_PREAMBLE = r"""\providecommand{\cfhistfigure}[3]{%
  \begin{figure}[H]
    \centering
    \IfFileExists{figures/#1.pdf}%
      {\includegraphics[width=0.95\linewidth,height=0.8\textheight,keepaspectratio]{figures/#1.pdf}}%
      {\pending{Figure \texttt{figures/\detokenize{#1}} is written by
        \texttt{bench/make\_report.py} once its data exist.}}
    \caption{#2}\label{#3}
  \end{figure}}
"""


def run(names: list[str] | None = None, outdir: Path = FIG,
        log=print) -> tuple[list[tuple[str, str]], list[str]]:
    outdir.mkdir(parents=True, exist_ok=True)
    figs, notes, secs = [], [], []
    for name, (cfn, sfn, title, monitor) in RUNS.items():
        if names and name not in names:
            continue
        cfd, sfd = RUN / cfn, RUN / sfn
        info: dict = {"notes": []}
        status = {}
        cf = sf = None
        for solver, d in (("coupledFoam", cfd), ("simpleFoam", sfd)):
            if not d.is_dir():
                status[solver] = "run not available yet (pending)"
            elif busy(d):
                status[solver] = "run in progress, skipped (pending)"
                notes.append(f"histories {name}: {solver} run {d.name} is being written, skipped")
            elif (solver == "coupledFoam"
                  and not provenance.active().check_run(d)):
                # staleness guard: run of another commit or without
                # provenance.json (listed in the report appendix)
                status[solver] = ("run not from the report's commit "
                                  "(pending)")
                notes.append(f"histories {name}: {d.name} stale, skipped")
        if "coupledFoam" not in status:
            cf = read_cf(cfd)
            if cf is None:
                status["coupledFoam"] = "no CF| lines in the log (pending)"
        if "simpleFoam" not in status:
            sf = read_sf(sfd)
            if sf is None:
                status["simpleFoam"] = "no simpleFoam log (pending)"
        made = []
        if cf is not None or sf is not None:
            ph = phases(cfd) if cf is not None else ([], {})
            if monitor is None:
                info["noLoads"] = True
                info["notes"].append("no monitored integral quantity in this "
                                     "case; convergence is judged by the "
                                     "residuals and the final centreline "
                                     "profiles")
            else:
                if monitor == "dp" and name.startswith("T2"):
                    info["notes"].append("the reattachment length is evaluated "
                                         "only at the end of the run (no "
                                         "history); the pressure difference is "
                                         "shown instead")
                st = fig_loads(name, title, cf, sf, cfd, sfd, monitor, ph, outdir, info)
                if st:
                    made.append("loads")
            for kind, fn in (("residuals", lambda: fig_residuals(name, title, cf, sf, ph, outdir)),
                             ("wall", lambda: fig_wall(name, title, cf, sf, ph, outdir, info)),
                             ("linear", lambda: fig_linear(name, title, cf, sf, ph, outdir, info))):
                try:
                    st = fn()
                except Exception as e:  # noqa: BLE001 - one broken panel must not stop the rest
                    notes.append(f"histories {name} {kind}: {type(e).__name__}: {e}")
                    st = None
                if st:
                    made.append(kind)
            if ph[0]:
                info["notes"].append("background colours: diagnostics phases "
                                     "(" + ", ".join(sorted({b[2] for b in ph[0]})) + ")")
            log(f"histories {name}: {', '.join(made) or 'nothing'}")
        for k in made:
            figs.append((f"hist_{name}_{k}", f"{title}: {k} history"))
        secs.append(tex_section(name, title, made, info, status))
    if not names:
        for group, members in (("wing", ["T3_kOmegaSST_np1", "T3_GEKO_np1"]),
                               ("motorbike", ["T4a_np10", "T4b_np10"]),
                               ("ahmed", [T5_RUN])):
            try:
                st = fig_overview(group, members, outdir, log)
            except Exception as e:  # noqa: BLE001
                notes.append(f"overview {group}: {type(e).__name__}: {e}")
                st = None
            if st:
                figs.append((st, f"Cd and Cl overview, {group}"))
            else:
                notes.append(f"overview {group}: no data yet (placeholder)")
    (outdir / "histories.tex").write_text(
        "% Generated by bench/plot_histories.py - do not edit\n" + HIST_PREAMBLE
        + "\n".join(secs))
    return figs, notes


def _cli(argv: list[str]) -> int:
    import argparse  # noqa: PLC0415
    ap = argparse.ArgumentParser()
    ap.add_argument("runs", nargs="*")
    ap.add_argument("--commit", default="HEAD")
    ap.add_argument("--allow-stale", action="store_true")
    a = ap.parse_args(argv)
    if not a.allow_stale:
        target = provenance.git_resolve(a.commit) or a.commit
        provenance.ACTIVE = provenance.Guard(target)
    f, n = run(a.runs or None)
    for x in n:
        print("note:", x)
    for r in provenance.active().rows():
        print(f"stale: {r['case']} {r['item']} found {r['found']}")
    return 0


if __name__ == "__main__":
    sys.exit(_cli(sys.argv[1:]))

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
                              at iters_to_stationary (D-042, T4/T5) or at
                              the 12.3(ii) point (T3); zoom panels on the
                              shorter run when the runs differ > 5x; a
                              simpleFoam reference continuation (mean
                              fields) is cut off, as in the tables
    hist_<run>_residuals.pdf  coupledFoam R, rU, rp and the simpleFoam
                              initial residuals of every solved field,
                              per iteration (left, own axis per solver) and
                              over wall time (right, common axis, inset
                              zoom on the short run), same y range
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
import os
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
import plotstyle as ps  # noqa: E402
import run_bench  # noqa: E402

BUSY_MINUTES = 10.0
# One colour pair everywhere (Okabe-Ito, colour-blind safe, bench/plotstyle):
# simpleFoam vermillion, coupledFoam blue
C_NATIVE = ps.C_NATIVE
C_COUPLED = ps.C_COUPLED
# per-field series (residuals, linear iterations), Okabe-Ito order
SERIES = ps.SERIES
NMAX = 1500      # points of a smooth curve (running means, fits)
SPLIT_RATIO = 5.0  # run lengths differing by more: zoom panels for the short run
SPLIT = {"tAsm": "#9ecae1", "tSolve": "#1f5fbf", "tTurb": "#fdae6b",
         "other": "#d9d9d9"}

def _t5_run() -> str:
    """T5 run the report shows, the same rule as bench/make_report.py
    _t5_run: CF_T5_MESH if set, else the variant (fine, then coarse) that
    has a test record, else the fine one (names: run_bench.t5_run_name)."""
    import os  # noqa: PLC0415
    env = os.environ.get("CF_T5_MESH")
    for v in ([env] if env else []) + ["fine", "coarse"]:
        n = run_bench.t5_run_name(run_bench.HEAVY_NP, v)
        if (REPO / "results" / "tests" / f"{n}.json").exists():
            return n
    return run_bench.t5_run_name(run_bench.HEAVY_NP, env or "fine")


# run name: coupledFoam dir, simpleFoam reference dir, title, monitor
T5_RUN = _t5_run()
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
# T5 (Ahmed body) is deferred by user decision (D-063): no appendix section
# and no overview figure unless CF_INCLUDE_T5=1
T5_DEFERRED = os.environ.get("CF_INCLUDE_T5") != "1"
if T5_DEFERRED:
    RUNS.pop(T5_RUN, None)

ps.apply()


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
                lst = store.get(f)
                if lst is None:      # first appearance: NaN before it
                    lst = store[f] = [np.nan] * (n - 1)
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


def _cut_reference(case: Path, h: dict | None) -> dict | None:
    """A continued simpleFoam reference is its original run only: samples
    beyond run_bench.reference_t_max(case) are dropped (M7, as in every
    evaluation of the reference)."""
    if h is None:
        return None
    t_max = run_bench.reference_t_max(case)
    if t_max is None:
        return h
    keep = h["iter"] <= t_max
    if not keep.any():
        return None
    return {k: v[keep] for k, v in h.items()}


def forces(case: Path) -> dict | None:
    """Cd, Cl, CmPitch per iteration, restarts concatenated (later start
    directories override earlier iterations); a continued reference is cut
    at its original budget."""
    return _cut_reference(case, _forces_all(case))


def _forces_all(case: Path) -> dict | None:
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
    """areaAverage(p) inlet - outlet per iteration (surfaceFieldValue); a
    continued reference is cut at its original budget."""
    return _cut_reference(case, _pressure_drop_all(case))


def _pressure_drop_all(case: Path) -> dict | None:
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
    """Diagnostics phases as light background colours."""
    b, col = ph
    for x0, x1, name in b:
        ax.axvspan(x0, x1, color=col.get(name, "#eeeeee"), alpha=0.35, lw=0,
                   zorder=0)


def events(ax, ev):
    """GAMG autoTune events: thin dotted vertical lines."""
    for it, what, a, b in ev:
        ax.axvline(it, color="k", lw=0.5, ls=":", alpha=0.7, zorder=1)


def small_legend(ax, **kw):
    """Legend of one panel, above the panel (outside the data)."""
    h, lab = ax.get_legend_handles_labels()
    if h:
        ax.legend(h, lab, ncol=kw.pop("ncol", 4), loc="lower left",
                  bbox_to_anchor=(0.0, 1.0), borderaxespad=0.1,
                  handlelength=1.5, columnspacing=1.0, **kw)


def decimate(x, y, nmax: int = NMAX):
    """Smooth curve with at most nmax points (bin means); exact below."""
    return ps.decimate(x, y, nmax)


def raw(ax, x, y, color, label=None, log=False, lw=ps.LW_RAW, alpha=0.9,
        z=None, nbins: int = 100):
    """Per-iteration data: a thin line when short (<= ps.RAW_MAX points),
    otherwise the per-bin min-max band with the bin median (plotstyle
    trace/envelope): thousands of raw points are never drawn as one
    polyline over another series."""
    ps.trace(ax, x, y, color, label=label, log=log, lw=lw, alpha=alpha,
             z=ps.Z_RAW if z is None else z, median_lw=0.7, nbins=nbins)


def main(ax, x, y, color, label=None, lw=ps.LW_MAIN, ls="-", zorder=ps.Z_MEAN):
    """Main curve (mean, median, cumulative, fit): strong, on top."""
    bx, by = decimate(x, y)
    ax.plot(bx, by, color=color, lw=lw, ls=ls, label=label, zorder=zorder)


def save(fig, stem: str, outdir: Path) -> str:
    ps.save(fig, outdir / f"{stem}.pdf")
    return stem


def _xmax(x) -> float:
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    return float(x.max()) if x.size else 0.0


def _clip_to(x, y, xmax):
    """The part of a history with x <= xmax (so an envelope in a zoom panel
    is binned over the visible range only)."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    k = np.isfinite(x) & (x <= xmax)
    return x[k], y[k]


# --------------------------------------------------------------------------- #
# figures
# --------------------------------------------------------------------------- #

def load_series(cf, sf, cfd, sfd, monitor, info):
    """[(solver, colour, iterations, wall time at each, {quantity: values},
    None)] (the last element is kept for the tuple layout of the callers).

    A continued simpleFoam reference (Allrun -restart for the mean fields,
    D-042 addendum) is its original run only: forces()/pressure_drop() cut
    it at run_bench.reference_t_max, as every evaluation of the reference
    (M7)."""
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
        t_max = run_bench.reference_t_max(case)
        if t_max is not None and info is not None:
            info["notes"].append(
                f"{solver}: the reference is its original run of "
                f"{int(t_max)} iterations, as in the tables; its "
                "continuation (mean fields, D-042 addendum) is neither "
                "shown nor evaluated")
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
        series.append((solver, col, it, tt, qs, None))
    return series


def _user_iteration(run: str, solver: str) -> int | None:
    """User-judged convergence iteration of a test run (D-060)."""
    try:
        import user_convergence as ucv  # noqa: PLC0415
        return ucv.lookup(ucv.load(), ucv.case_key(run), solver)
    except Exception:  # noqa: BLE001 - the plot must not fail on the file
        return None


def _case_of(run: str) -> str | None:
    """run_bench case key of a run name (T4a_np10, ref_T4a_np10 -> T4a;
    T3_kOmegaSST_np1 -> T3-SST; T5_coarse_np10 -> T5; T0 -> None). Every
    stat_window / iters_to_stationary / stationary_eval call of this module
    gets it, so both solvers of a wake case share one window (D-068)."""
    return run_bench.case_of_run(run)


def stat_window(n: int, case: str | None = None) -> int:
    return run_bench.stat_window(n, case)


def _oscillatory(run: str) -> bool:
    """D-042 wake case (T4, T5): stationary window mean; otherwise (T3)
    criterion 12.3(ii)."""
    return run_bench.is_oscillatory(run_bench.case_of_run(run) or "")


LOAD_LABELS = {"Cd": "$C_d$", "Cl": "$C_l$", "CmPitch": "$C_m$ (pitch)",
               "dp": r"$\Delta p$ inlet$-$outlet [m$^2$/s$^2$]"}


def _running(y, W):
    """Running mean and RMS over the D-042 window W (full windows only)."""
    m = roll(y, W, np.nanmean)
    sd = roll_std(y, W)
    m[:max(W - 1, 0)] = np.nan
    sd[:max(W - 1, 0)] = np.nan
    return m, sd


def _zoom_ranges(xs: list[float]) -> float | None:
    """x limit of a zoom on the shorter of two runs on a common axis, or
    None when the runs differ by less than SPLIT_RATIO."""
    xs = [x for x in xs if x and np.isfinite(x)]
    if len(xs) < 2 or max(xs) <= SPLIT_RATIO * min(xs):
        return None
    return 1.04 * min(xs)


def _legend_handles(series, marks: dict, extra_label=None, case=None):
    """Figure legend of the load plots: per solver the per-iteration data
    and the running mean with its RMS band, then the markers."""
    from matplotlib.lines import Line2D  # noqa: PLC0415
    from matplotlib.patches import Patch  # noqa: PLC0415
    h, lab = [], []
    for solver, col, it, tt, qs, cont in sorted(series, key=lambda s_: s_[0] != "simpleFoam"):
        n = len(it)
        W = stat_window(n, case)
        if n > ps.RAW_MAX:
            h.append(Patch(color=col, alpha=ps.BAND_ALPHA, lw=0))
            lab.append(f"{solver} per iteration (min-max per bin)")
        else:
            h.append(Line2D([], [], color=col, lw=ps.LW_RAW))
            lab.append(f"{solver} per iteration")
        h.append((Patch(color=col, alpha=0.18, lw=0),
                  Line2D([], [], color=col, lw=ps.LW_MAIN, ls="--")))
        lab.append(f"{solver} running mean $\\pm$RMS ($W$={W})")
    if marks.get("window"):
        h.append(Patch(color="#777777", alpha=0.18, lw=0))
        lab.append("final averaging window")
    if marks.get("stat"):
        h.append(Line2D([], [], color="#333333", lw=ps.LW_MARK, ls="-."))
        lab.append("iterations to a stationary window (D-042)")
    if marks.get("w123"):
        h.append(Line2D([], [], color="#333333", lw=ps.LW_MARK, ls=LS_W123))
        lab.append("criterion 12.3(ii) met (force window)")
    if marks.get("user"):
        h.append(Line2D([], [], color="#333333", lw=1.8, ls="-"))
        lab.append("user convergence point (D-060)")
    if extra_label:
        h.append(Line2D([], [], color="k", lw=1.0))
        lab.append(extra_label)
    return h, lab


LS_W123 = (0, (6, 1.5, 1, 1.5, 1, 1.5))   # dash-dot-dot: 12.3(ii) point


def _vline(ax, x, col, kind):
    if x is None or not np.isfinite(x):
        return
    if kind == "stat":
        ax.axvline(x, color=col, lw=ps.LW_MARK, ls="-.", zorder=ps.Z_MARK)
    elif kind == "w123":
        ax.axvline(x, color=col, lw=ps.LW_MARK, ls=LS_W123, zorder=ps.Z_MARK)
    else:
        ax.axvline(x, color=col, lw=1.8, ls="-", zorder=ps.Z_MARK + 0.5)


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
    # the longer run first, the shorter one is drawn on top of it
    series = sorted(series, key=lambda s_: -len(s_[2]))
    case = _case_of(run)
    osc = _oscillatory(run)
    # ---- evaluation: the averaging window W, the convergence point (D-042
    # stationary window on the wake cases T4/T5, criterion 12.3(ii) on the
    # other force cases), the user point (D-060)
    evals = {}
    for solver, col, it, tt, qs, cont in series:
        e = {"W": {}, "run": {}}
        for q in quants:
            if q not in qs:
                continue
            y = qs[q]
            n = len(y)
            W = stat_window(n, case)
            e["W"][q] = W
            e["run"][q] = _running(y, W)
            if q in ("Cd", "Cl") and monitor == "forces":
                if "its" not in e:
                    hist = {k: list(v) for k, v in qs.items() if k in ("Cd", "Cl")}
                    ev = {}
                    try:
                        if osc:
                            its = run_bench.iters_to_stationary(hist, case=case)
                            ev = run_bench.stationary_eval(hist, case=case)
                        else:
                            its = run_bench.iters_to_conv(hist)
                    except Exception:  # noqa: BLE001
                        its, ev = None, {}
                    if not osc:
                        # window means over the plotted averaging window
                        ev = {"W": W, "Cd_mean": float(np.mean(qs["Cd"][n - W:]))
                              if "Cd" in qs and W else None,
                              "Cl_mean": float(np.mean(qs["Cl"][n - W:]))
                              if "Cl" in qs and W else None}
                    e["its"] = its
                    e["kind"] = "stat" if osc else "w123"
                    uit = _user_iteration(run, solver)
                    e["user"] = uit
                    if uit:
                        info.setdefault("user", {})[solver] = uit
                    info.setdefault("stat", {})[solver] = (its, ev.get("W"),
                                                           ev.get("Cd_mean"),
                                                           ev.get("Cl_mean"),
                                                           "D-042" if osc else "12.3(ii)")
            else:
                w = stat_window(n, case)
                win = y[n - w:] if w else y
                info.setdefault("mon", {})[(solver, q)] = (float(np.mean(win)),
                                                           float(np.std(win)))
        evals[solver] = e

    def at_time(it, tt, x):
        """wall time of iteration x (for the markers in the right column)."""
        if x is None:
            return None
        k = min(int(np.searchsorted(it, x)), len(tt) - 1)
        return tt[k] if np.isfinite(tt[k]) else None

    # ---- layout. Every panel shows both solvers on a common axis (fair
    # comparison). Runs whose lengths differ by more than SPLIT_RATIO get a
    # zoom panel on the shorter run: one extra column when only one of the
    # two axes (iterations, wall time) needs it, otherwise a second row per
    # quantity.
    it_end = [_xmax(s_[2]) for s_ in series]
    t_end = [_xmax(s_[3]) for s_ in series]
    zit, zt = _zoom_ranges(it_end), _zoom_ranges(t_end)
    short_it = min(series, key=lambda s_: _xmax(s_[2]))
    short_t = min(series, key=lambda s_: _xmax(s_[3]) or np.inf)
    full_it = 1.02 * max(it_end)
    full_t = 1.02 * max(t_end) if max(t_end) > 0 else None
    if zit is not None and zt is not None:
        cols = [("it", None), ("t", None)]
        rows = [(q, kind) for q in quants for kind in ("full", "zoom")]
    else:
        cols = [("it", None)]
        if zit is not None:
            cols.append(("it", zit))
        cols.append(("t", None))
        if zt is not None:
            cols.append(("t", zt))
        rows = [(q, "full") for q in quants]
    nr, nc = len(rows), len(cols)
    hrow = max(1.1, min(2.3, 6.4 / nr)) if nc == 2 else max(1.35, min(2.2, 6.0 / nr))
    fig, axs = ps.subplots(nr, nc, height=hrow * nr + 0.8)
    marks = {"window": False, "stat": False, "w123": False, "user": False}
    for r, (q, kind) in enumerate(rows):
        for c, (axis, zlim) in enumerate(cols):
            ax = axs[r, c]
            if kind == "zoom":
                zlim = zit if axis == "it" else zt
            lim_series = []
            for k, (solver, col, it, tt, qs, cont) in enumerate(series):
                if q not in qs:
                    continue
                y = qs[q]
                n = len(y)
                W = evals[solver]["W"][q]
                m, sd = evals[solver]["run"][q]
                x = it if axis == "it" else tt
                ok = np.isfinite(x)
                xs, ys, ms, sds = x[ok], y[ok], m[ok], sd[ok]
                if zlim is not None:
                    sel = xs <= zlim
                    xs, ys, ms, sds = xs[sel], ys[sel], ms[sel], sds[sel]
                if not xs.size:
                    continue
                dz = 0.1 * k
                raw(ax, xs, ys, col, z=ps.Z_RAW + dz)
                bx, bm = decimate(xs, ms)
                _, bs = decimate(xs, sds)
                ps.fill(ax, bx, bm - bs, bm + bs, color=col, alpha=0.18,
                        lw=0, zorder=3 + dz)
                ax.plot(bx, bm, color=col, lw=ps.LW_MAIN, ls="--",
                        zorder=ps.Z_MEAN + dz)
                if W >= 2 and n >= W and np.isfinite(x[n - W]):
                    ax.axvspan(x[n - W], np.nanmax(x), color=col, alpha=0.08,
                               lw=0, zorder=ps.Z_SHADE)
                    marks["window"] = True
                e = evals[solver]
                if "its" in e:
                    if e["its"]:
                        _vline(ax, e["its"] if axis == "it" else at_time(it, tt, e["its"]),
                               col, e["kind"])
                        marks[e["kind"]] = True
                    if e.get("user"):
                        _vline(ax, e["user"] if axis == "it" else at_time(it, tt, e["user"]),
                               col, "user")
                        marks["user"] = True
                lim_series.append((xs, ys, col))
            if axis == "it" and f"{q}Mean" in extra:
                xe, ye = cf["iter"], extra[f"{q}Mean"]
                if zlim is not None:
                    xe, ye = _clip_to(xe, ye, zlim)
                main(ax, xe, ye, "k", lw=1.0, zorder=ps.Z_MARK)
            # y range from the post-transient part; the clipped start-up
            # transient is marked by triangles at the axis edge
            ps.robust_ylim(ax, lim_series,
                           skip_frac=0.25 if zlim is None else 0.4,
                           pct=(0.5, 99.5), margin=0.2)
            if zlim is not None:
                ax.set_xlim(0, zlim)
                s_ = short_it if axis == "it" else short_t
                if nc == 2 or r == 0:
                    ax.set_title(f"zoom: the {s_[0]} run" if nc > 2 else
                                 (f"zoom: the {len(s_[2])} {s_[0]} iterations"
                                  if axis == "it" else
                                  f"zoom: the {s_[0]} run ({_xmax(s_[3]):.3g} s)"),
                                 loc="left", fontsize=ps.FS_NOTE, color=ps.GREY)
            else:
                ax.set_xlim(0, full_it if axis == "it" else full_t)
            if axis == "it":
                bands(ax, ph)
        lab = LOAD_LABELS.get(q, q)
        axs[r, 0].set_ylabel(lab if kind == "full" else f"{lab} (zoom)")
    for c, (axis, zlim) in enumerate(cols):
        axs[-1, c].set_xlabel("outer iteration" if axis == "it" else
                              "wall-clock time [s]")
    h, lab = _legend_handles(series, marks,
                             "coupledFoam window mean (solver output)"
                             if any(f"{q}Mean" in extra for q in quants) else None,
                             case=case)
    ps.legend_below(fig, h, lab, ncol=2, h_pad=0.6)
    return save(fig, f"hist_{run}_loads", outdir)


def _residual_curves(solver, d):
    if solver == "coupledFoam":
        curves = [(k, d[k]) for k in ("R", "rU", "rp") if k in d]
        return curves, d["iter"], d["_t"]
    order = ["p", "Ux", "Uy", "Uz", "k", "omega", "epsilon", "nuTilda"]
    fields = sorted(d["res"], key=lambda f: order.index(f) if f in order else 99)
    return ([(f, d["res"][f]) for f in fields], np.arange(1, d["n"] + 1),
            d["t"])


def fig_residuals(run, title, cf, sf, ph, outdir):
    rows = [x for x in (("coupledFoam", cf), ("simpleFoam", sf)) if x[1] is not None]
    if not rows:
        return None
    data = [(s_, *_residual_curves(s_, d)) for s_, d in rows]
    tend = [_xmax(t) for _, _, _, t in data]
    tmax = max(tend)
    tzoom = _zoom_ranges(tend)
    fig, axs = ps.subplots(len(rows), 2, height=2.45 * len(rows) + 0.35)
    ylo, yhi = [], []
    for r, (solver, curves, x0, t) in enumerate(data):
        for j, (lab, y) in enumerate(curves):
            col = SERIES[j % len(SERIES)]
            for c, x in ((0, x0), (1, t)):
                x = np.asarray(x[:len(y)], float)
                ps.trace(axs[r, c], x, y[:len(x)], col, label=lab if c == 0 else None,
                         log=True, lw=1.0, band_alpha=0.14, median_lw=1.2,
                         z=ps.Z_RAW + 0.1 * j)
            fy = y[np.isfinite(y) & (y > 0)]
            if fy.size:
                ylo.append(fy.min())
                yhi.append(fy.max())
        for c in (0, 1):
            axs[r, c].set_yscale("log")
        col = C_COUPLED if solver == "coupledFoam" else C_NATIVE
        axs[r, 0].set_ylabel(f"{solver}\nresidual", color=col)
        small_legend(axs[r, 0], ncol=len(curves))
        axs[r, 0].set_xlim(0, 1.02 * _xmax(x0))
        # both solvers on the same wall-clock axis (fair comparison)
        axs[r, 1].set_xlim(0, 1.02 * tmax)
        bands(axs[r, 0], ph if solver == "coupledFoam" else ([], {}))
        if solver == "coupledFoam":
            events(axs[r, 0], rows[r][1]["_events"])
    if ylo:
        for ax in axs.flat:
            ax.set_ylim(min(ylo) / 2, max(yhi) * 2)
    # the short run on the common wall-clock axis: zoom inset
    if tzoom is not None:
        for r, (solver, curves, x0, t) in enumerate(data):
            if _xmax(t) * SPLIT_RATIO < tmax:
                ins = ps.zoom_inset(axs[r, 1], (0, tzoom), bounds=(0.30, 0.30, 0.67, 0.62),
                                    label=f"zoom: the {solver} run", ylog=True)
                for j, (lab, y) in enumerate(curves):
                    x = np.asarray(t[:len(y)], float)
                    ps.trace(ins, x, y[:len(x)], SERIES[j % len(SERIES)], log=True,
                             lw=1.0, band_alpha=0.14, median_lw=1.1, nbins=120)
                ins.set_ylim(axs[r, 1].get_ylim())
    for r in range(len(rows)):
        axs[r, 0].set_xlabel(f"{data[r][0]} outer iteration")
    axs[-1, 1].set_xlabel("wall-clock time [s] (common axis)")
    ps.tight(fig, h_pad=0.8)
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
    fig = ps.figure(6.9 if cf is not None else 4.6)
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
        z = ps.Z_MEAN + (0.2 if solver == "coupledFoam" else 0.0)
        for c in range(ncol):
            main(top[c], it, t, col, label=lab if c == 0 else None, lw=1.8, zorder=z)
            main(top[c], it, slope * it + icpt, col, lw=1.0,
                 ls=":", zorder=z + 1,
                 label=(f"{solver} linear fit {slope:.3g} s/it ($R^2$ {r2:.3f})"
                        if c == 0 else None))
            xx, pp = it, per
            if zoom and c == 1:
                xx, pp = _clip_to(it, per, lens["coupledFoam"] * 1.05)
            ps.trace(mid[c], xx, pp, col, log=True, lw=0.6, alpha=0.5,
                     z=ps.Z_RAW, band_alpha=0.18, line=False)
            main(mid[c], it, med, col, lw=1.5, zorder=z,
                 label=f"{solver} time per iteration, rolling median" if c == 0 else None)
            mid[c].set_yscale("log")
    if zoom:
        nz = lens["coupledFoam"] * 1.05
        for ax in (top[1], mid[1]):
            ax.set_xlim(0, nz)
        top[1].set_title(f"zoom: the {lens['coupledFoam']} coupledFoam iterations",
                         fontsize=ps.FS_NOTE)
        # y range of the zoom: what both solvers reach within nz iterations
        ymax = 0.0
        for s_, d, _ in have:
            t, it = _timing(s_, d)[:2]
            ymax = max(ymax, float(np.nanmax(t[it <= nz])) if np.any(it <= nz) else 0.0)
        top[1].set_ylim(0, ymax * 1.05)
        top[0].set_title("all iterations", fontsize=ps.FS_NOTE)
    for c in range(ncol):
        top[c].set_xlim(0, None)
        mid[c].set_xlim(top[c].get_xlim())
    top[0].set_ylabel("cumulative wall time [s]")
    mid[0].set_ylabel("wall time per\niteration [s]")
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
            bx = decimate(it, it)[0]
            ax3.stackplot(bx, *[decimate(it, v)[1] for v in parts.values()],
                          labels=[f"coupledFoam {k}" for k in parts],
                          colors=[SPLIT[k] for k in parts], lw=0)
            ax3.set_ylabel("coupledFoam time per\niteration [s]\n(rolling median)")
            events(ax3, cf["_events"])
            ax3.set_xlim(0, _xmax(it))
        bands(ax3, ph)
    for ax in axes:
        ax.set_xlabel("outer iteration")
    ps.legend_below(fig, ncol=2, h_pad=0.8)
    return save(fig, f"hist_{run}_wall", outdir)
def fig_overview(group: str, runs: list[str], outdir: Path, log=print) -> str | None:
    """Main-text overview: Cd and Cl of both solvers over wall-clock time
    for a group of force cases, with the final D-042 window means. A case
    whose two runs differ in wall time by more than SPLIT_RATIO gets a
    second row zoomed on the shorter run (same axes for both solvers)."""
    rows = []
    for name in runs:
        cfn, sfn, title, monitor = RUNS[name]
        cfd, sfd = RUN / cfn, RUN / sfn
        cf = (read_cf(cfd) if cfd.is_dir() and not busy(cfd)
              and provenance.active().check_run(cfd) else None)
        sf = read_sf(sfd) if sfd.is_dir() and not busy(sfd) else None
        rows.append((title, load_series(cf, sf, cfd, sfd, "forces", None),
                     _case_of(name)))
    if not any(r[1] for r in rows):
        return None
    layout = []
    for title, series, case in rows:
        series = sorted(series, key=lambda s_: -len(s_[2]))
        layout.append((title, series, None, case))
        tz = _zoom_ranges([_xmax(s_[3]) for s_ in series])
        if tz is not None:
            layout.append((title, series, tz, case))
    nr = len(layout)
    hrow = max(1.3, min(2.2, 6.0 / nr))
    fig, axs = ps.subplots(nr, 2, height=hrow * nr + 0.6)
    long_any = False
    for r, (title, series, tz, case) in enumerate(layout):
        for c, q in enumerate(("Cd", "Cl")):
            ax = axs[r, c]
            if not series:
                ps.note(ax, "results pending", loc="center")
                ax.set_yticks([])
            lim, inc, tmax = [], [], 0.0
            for k, (solver, col, it, tt, qs, cont) in enumerate(series):
                if q not in qs:
                    continue
                y = qs[q]
                n = len(y)
                W = stat_window(n, case)
                m = roll(y, W, np.nanmean)
                m[:max(W - 1, 0)] = np.nan
                ok = np.isfinite(tt)
                x, yy, mm = tt[ok], y[ok], m[ok]
                tmax = max(tmax, _xmax(x))
                if tz is not None:
                    sel = x <= tz
                    x, yy, mm = x[sel], yy[sel], mm[sel]
                if not x.size:
                    continue
                long_any |= len(x) > ps.RAW_MAX
                dz = 0.1 * k
                raw(ax, x, yy, col, z=ps.Z_RAW + dz)
                main(ax, x, mm, col, lw=ps.LW_MAIN, ls="--", zorder=ps.Z_MEAN + dz)
                mw = float(np.mean(y[n - W:])) if W and n >= W else np.nan
                if np.isfinite(mw):
                    ax.axhline(mw, color=col, lw=1.0, ls=":", zorder=ps.Z_MARK)
                    inc.append(mw)
                lim.append((x, yy, col))
            if lim:
                ps.robust_ylim(ax, lim, skip_frac=0.25, pct=(0.5, 99.5),
                               margin=0.2)
                ax.set_xlim(0, tz if tz is not None else 1.02 * tmax)
            ax.set_ylabel("$C_d$" if q == "Cd" else "$C_l$")
        head = title if tz is None else f"{title}: zoom on the shorter run"
        axs[r, 0].set_title(head, loc="left", fontsize=ps.FS)
    for c in range(2):
        axs[-1, c].set_xlabel("wall-clock time [s]")
    from matplotlib.lines import Line2D  # noqa: PLC0415
    from matplotlib.patches import Patch  # noqa: PLC0415
    h = [Line2D([], [], color=C_NATIVE, lw=ps.LW_MAIN),
         Line2D([], [], color=C_COUPLED, lw=ps.LW_MAIN),
         Line2D([], [], color="k", lw=ps.LW_MAIN, ls="--"),
         Line2D([], [], color="k", lw=1.0, ls=":")]
    lab = ["simpleFoam", "coupledFoam", "running mean over the D-042 window",
           "final window mean"]
    if long_any:
        h.append(Patch(color="#777777", alpha=ps.BAND_ALPHA, lw=0))
        lab.append("per iteration (band: min-max per bin)")
    ps.legend_below(fig, h, lab, ncol=3, h_pad=0.7)
    stem = f"loads_overview_{group}"
    save(fig, stem, outdir)
    log(f"overview {group}: written")
    return stem


KNOWN = {"iter", "CFL", "omega", "cuts", "R", "rU", "rp", "linIters",
         "linRes", "tAsm", "tSolve", "tTurb", "tIter", "tWall", "nStat",
         "nDyn", "nLocLim", "nRollback", "nClamped", "eta", "rho", "nBound",
         "aa"}
CYCLES = {"V": 0, "F": 1, "W": 2, "K": 3}


def _cycle_series(cf, it):
    """GAMG cycle type per outer iteration (0..3, NaN if unknown) from the
    start value in the log and the autoTune events."""
    cur = cf["_cycle0"]
    y = np.full(len(it), np.nan)
    evs = sorted((e for e in cf["_events"] if e[1] == "cycleType"), key=lambda e: e[0])
    j = 0
    for i, x in enumerate(it):
        while j < len(evs) and evs[j][0] <= x:
            cur = evs[j][3]
            j += 1
        y[i] = CYCLES.get(cur, np.nan) if cur else np.nan
    return y


def _npost_series(cf, it):
    cur = cf["_npost0"]
    y = np.full(len(it), np.nan)
    evs = sorted((e for e in cf["_events"] if e[1] == "nPostSweeps"), key=lambda e: e[0])
    if cur is None and evs:
        cur = int(evs[0][2])
    j = 0
    for i, x in enumerate(it):
        while j < len(evs) and evs[j][0] <= x:
            cur = int(evs[j][3])
            j += 1
        y[i] = cur if cur is not None else np.nan
    return y


def _draw_series(ax, x, y, col, name, kind, log):
    """One control history: exact step/line when short, envelope when long."""
    if len(x) > ps.RAW_MAX:
        ps.trace(ax, x, y, col, label=name, log=log, lw=0.9, band_alpha=0.2,
                 median_lw=1.1, nbins=200)
    elif kind == "step":
        ax.step(x, y, where="mid", color=col, lw=0.9, label=name)
    else:
        ax.plot(x, y, color=col, lw=0.9, label=name)


def fig_linear(run, title, cf, sf, ph, outdir, info):
    if cf is None and sf is None:
        return None
    panels = []
    if cf is not None:
        it = cf["iter"]
        cyc = _cycle_series(cf, it)
        npost = _npost_series(cf, it)
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
    # quantities constant over the whole run are listed in the text instead
    # of a panel (also single constant series of a multi-series panel)
    keep = []
    for p in panels:
        lab, data, scale, kind = p
        if data in ("cycle", "npost"):
            y = cyc if data == "cycle" else npost
            fails = [e for e in cf["_events"] if e[1] in ("failure", "kcycle-saturated")]
            if not np.any(np.isfinite(y)):
                info["notes"].append(("GAMG cycle type" if data == "cycle" else "nPostSweeps")
                                     + " not in the log (panel omitted)")
                continue
            if ps.is_constant(y) and not (data == "cycle" and fails):
                v = np.nanmax(y)
                info.setdefault("constant", []).append(
                    f"GAMG cycle = {[c for c, i in CYCLES.items() if i == v][0]}"
                    if data == "cycle" else f"nPostSweeps = {v:g}")
                continue
            keep.append(p)
            continue
        ser = [(nm, y) for nm, y in data if y is not None and np.any(np.isfinite(y))]
        const = [(nm, y) for nm, y in ser if ps.is_constant(y)]
        info.setdefault("constant", []).extend(
            f"{nm} = {np.nanmax(y):g}" for nm, y in const)
        var = [(nm, y) for nm, y in ser if not ps.is_constant(y)]
        if not var:
            continue
        keep.append((lab, var, scale, kind))
    panels = keep
    if sf is not None:
        panels.append(("simpleFoam linear\niterations per field", "sfnit", "log", None))
    if not panels:
        return None
    npan = len(panels)
    ncol = 2 if npan > 3 else 1
    nrow = int(np.ceil(npan / ncol))
    fig, axs = ps.subplots(nrow, ncol, height=1.2 * nrow + 0.45)
    # column-major: the coupledFoam panels of one column share the
    # iteration axis; the simpleFoam panel (last) has its own
    cells = [(i % nrow, i // nrow) for i in range(npan)]
    for i in range(npan, nrow * ncol):
        axs[i % nrow, i // nrow].axis("off")
    first_in_col = {}
    for (lab, data, scale, kind), (r, c) in zip(panels, cells):
        ax = axs[r, c]
        if data != "sfnit" and cf is not None:
            if c in first_in_col:
                ax.sharex(first_in_col[c])
            else:
                first_in_col[c] = ax
        if data == "cycle":
            ax.step(it, cyc, where="post", color=C_COUPLED, lw=1.0)
            ax.set_yticks(list(CYCLES.values()))
            ax.set_yticklabels(list(CYCLES))
            ax.set_ylim(-0.5, 3.5)
            for e in cf["_events"]:
                if e[1] in ("failure", "kcycle-saturated"):
                    ax.axvline(e[0], color="#d62728", lw=0.7, ls="--")
                    ax.text(e[0], 3.3, e[1], fontsize=ps.FS_NOTE, color="#d62728",
                            va="top")
        elif data == "npost":
            ax.step(it, npost, where="post", color=C_COUPLED, lw=1.0)
            ax.set_ylim(np.nanmin(npost) - 0.5, np.nanmax(npost) + 0.5)
        elif data == "sfnit":
            xs = np.arange(1, sf["n"] + 1)
            order = ["p", "Ux", "Uy", "Uz", "k", "omega", "epsilon", "nuTilda"]
            for j, f in enumerate(sorted(sf["nit"], key=lambda f: order.index(f) if f in order else 99)):
                ps.trace(ax, xs, np.maximum(sf["nit"][f], 0.8), SERIES[j % len(SERIES)],
                         label=f, lw=0.9, band_alpha=0.12, median_lw=1.1, nbins=200,
                         z=ps.Z_RAW + 0.1 * j)
            allv = np.concatenate([v[np.isfinite(v)] for v in sf["nit"].values()] or [np.array([1.0])])
            if allv.size and allv.max() > 10 * max(allv.min(), 1):
                ax.set_yscale("log")
            ps.inner_legend(ax, ncol=3)
            ax.set_xlim(0, sf["n"] * 1.02)
            ax.set_xlabel("simpleFoam iteration (own axis)")
        else:
            for j, (name, y) in enumerate(data):
                col = SERIES[j % len(SERIES)] if len(data) > 1 else C_COUPLED
                _draw_series(ax, it, y, col, name, kind, scale == "log")
            fin = np.concatenate([y[np.isfinite(y) & (y > 0)] for _, y in data])
            if scale == "log" and fin.size and fin.max() > 10 * fin.min():
                ax.set_yscale("log")
            elif scale == "symlog":
                ax.set_yscale("symlog", linthresh=1)
            if len(data) > 1:
                ps.inner_legend(ax, ncol=min(2, len(data)))
        ax.set_ylabel(lab)
        if cf is not None and data != "sfnit":
            bands(ax, ph)
            events(ax, cf["_events"])
            ax.set_xlim(0, _xmax(it) * 1.02)
    # x tick labels only on the last coupledFoam panel of each column (and
    # on the simpleFoam panel)
    for c in range(ncol):
        col_cells = [(k, p) for k, (p, cell) in enumerate(zip(panels, cells))
                     if cell[1] == c and p[1] != "sfnit"]
        for n_, (k, p) in enumerate(col_cells):
            ax = axs[cells[k]]
            if n_ < len(col_cells) - 1:
                ax.tick_params(labelbottom=False)
            else:
                ax.set_xlabel("coupledFoam outer iteration")
    if cf is not None:
        info["linTotal"] = float(np.nansum(cf.get("linIters", np.array([np.nan]))))
    ps.tight(fig, h_pad=0.4, w_pad=1.5)
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
        its, W, cd, cl, kind = v
        if kind == "D-042":
            s.append(f"{solver} D-042 window $W={_fmt(W, '{}')}$: "
                     f"$\\bar C_d={_fmt(cd, '{:.4f}')}$, $\\bar C_l={_fmt(cl, '{:.4f}')}$, "
                     f"stationary from iteration {_fmt(its, '{}')}")
        else:
            s.append(f"{solver} means over the last $W={_fmt(W, '{}')}$ "
                     f"iterations: $\\bar C_d={_fmt(cd, '{:.4f}')}$, "
                     f"$\\bar C_l={_fmt(cl, '{:.4f}')}$; criterion 12.3(ii) "
                     + (f"met at iteration {its}" if its else
                        "not met within the run"))
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
    forces_case = (RUNS.get(run) or (None,) * 4)[3] == "forces"
    marker = ("" if not forces_case else
              ", dash-dot: iterations to a stationary window (D-042)"
              if _oscillatory(run) else
              ", dash-dot-dot: criterion 12.3(ii) met")
    caps = {
        "loads": "monitored quantities per outer iteration (left) and over "
                 "wall-clock time (right), both solvers on the same axes; "
                 "thin line or light band: per-iteration values (long runs: "
                 "min--max per bin with the bin median), dashed with band: "
                 "running mean $\\pm$RMS over the averaging window $W$, "
                 f"shaded: final averaging window{marker}"
                 + (", solid: user convergence point (D-060, if set)" if forces_case else "")
                 + "; triangles at the axis "
                 "edge mark start-up values outside the plotted range; runs "
                 "of very different length get zoom panels on the shorter "
                 "run (an extra column, or a second row per quantity)",
        "residuals": "residuals per outer iteration (left, each solver on its "
                     "own iteration axis) and over wall-clock time (right, "
                     "same axis for both solvers; inset: zoom on the shorter "
                     "run), coupledFoam top, simpleFoam bottom, same y range; "
                     "long runs: light band = min--max per bin, line = bin "
                     "median; dotted vertical lines: GAMG autoTune events",
        "wall": "cumulative wall time with linear fit (top), wall time per "
                "iteration with rolling median (middle; light band: min--max "
                "per bin), coupledFoam time split into assembly, linear "
                "solve, turbulence and other (bottom); right column, if "
                "present: zoom on the coupledFoam iterations",
        "linear": "linear solver and pseudo-time controls per outer "
                  "iteration (the coupledFoam panels of one column share the "
                  "iteration axis; quantities constant over the run are "
                  "listed above instead of plotted; dotted vertical lines: "
                  "GAMG autoTune events); last panel: simpleFoam linear "
                  "iterations per field over its own iterations",
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
    ps.apply()      # other figure modules may have changed the rcParams
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
                               ("ahmed", [] if T5_DEFERRED else [T5_RUN])):
            if not members:          # deferred case (D-063)
                continue
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

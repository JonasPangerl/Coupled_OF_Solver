#!/usr/bin/env python3
"""Shared figure style of the report (paper.tex and paper_tutorial.tex).

One place for the print geometry, font sizes, line widths, colours and the
drawing helpers used by bench/make_report.py and bench/plot_histories.py
(and any other figure module):

    WIDTH, MAX_HEIGHT   figure size that \\cffigure / \\cfhistfigure print
                        at 1:1 (0.95 of the text width of paper.tex), so
                        the font sizes below are the sizes on paper
    apply()             matplotlib rcParams of the report
    C_COUPLED, C_NATIVE fixed solver colours (captions refer to them:
                        coupledFoam blue, simpleFoam orange/vermillion)
    SERIES              secondary colour-blind safe palette (Okabe-Ito)
    trace()             one history: a thin line when short, otherwise
                        envelope() (per-bin min-max band + bin median)
    envelope()          per-bin min-max band (light) and bin median (bold)
    decimate()          bin means of a smooth curve (running means, fits)
    robust_ylim()       y limits from post-transient percentiles with a
                        margin; values beyond the limits are marked by
                        small triangles at the axis edge
    zoom_inset()        inset axes for the short run on a common axis
    legend_below()      one figure legend under the panels
    save()              fixed-size PDF (and optional PNG), no tight bbox,
                        so the printed scale is exactly 1:1

Presentation only: nothing in here evaluates or changes a number.
"""

from __future__ import annotations

import re
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

REPO = Path(__file__).resolve().parents[1]
PAPER_TEX = REPO / "report" / "paper" / "paper.tex"


# --------------------------------------------------------------------------- #
# geometry
# --------------------------------------------------------------------------- #

def _geometry() -> tuple[float, float]:
    """Text width and height in inches from the geometry package options of
    paper.tex (a4paper, margin=...); A4 with 25 mm margins as fallback."""
    paper_w, paper_h, margin = 210.0, 297.0, 25.0
    try:
        tex = PAPER_TEX.read_text(errors="replace")
        m = re.search(r"\\usepackage\[([^\]]*)\]\{geometry\}", tex)
        if m:
            mm = re.search(r"margin\s*=\s*([\d.]+)\s*(mm|cm|in)", m.group(1))
            if mm:
                v = float(mm.group(1))
                margin = v * {"mm": 1.0, "cm": 10.0, "in": 25.4}[mm.group(2)]
        if re.search(r"\\documentclass\[[^\]]*letterpaper", tex):
            paper_w, paper_h = 215.9, 279.4
    except OSError:
        pass
    return (paper_w - 2 * margin) / 25.4, (paper_h - 2 * margin) / 25.4


TEXTWIDTH, TEXTHEIGHT = _geometry()          # 6.30 in x 9.72 in
# \cffigure and \cfhistfigure include every figure at 0.95\linewidth: a
# figure of WIDTH inches is printed at 1:1
WIDTH = 0.95 * TEXTWIDTH
# \cfhistfigure: height at most 0.8\textheight (keepaspectratio)
MAX_HEIGHT = 0.8 * TEXTHEIGHT

# font sizes in pt at print size (the figures are printed at 1:1)
FS = 8.0            # axis labels, panel titles
FS_TICK = 7.5
FS_LEGEND = 7.5
FS_NOTE = 7.0       # in-panel annotations

# line widths
LW_MAIN = 1.5       # running mean, median, fits
LW_RAW = 0.7        # short raw histories
LW_MARK = 1.1       # vertical/horizontal markers
BAND_ALPHA = 0.22   # envelope bands

# z-order (lowest first): background shading, envelopes, raw lines,
# means, markers
Z_SHADE, Z_ENV, Z_RAW, Z_MEAN, Z_MARK = 0.5, 2.0, 2.5, 4.0, 6.0

# --------------------------------------------------------------------------- #
# colours
# --------------------------------------------------------------------------- #

# fixed solver colours (Okabe-Ito); the captions refer to them
C_COUPLED = "#0072B2"      # coupledFoam: blue
C_NATIVE = "#D55E00"       # simpleFoam: vermillion ("orange" in captions)
C_COUPLED2 = "#56B4E9"     # second coupledFoam variant: sky blue
C_NATIVE2 = "#E69F00"      # second simpleFoam variant: orange
SOLVER_COLOR = {"coupledFoam": C_COUPLED, "simpleFoam": C_NATIVE}
# per-field series (residuals, linear iterations), colour-blind safe
SERIES = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9",
          "#000000", "#999999"]
GREY = "#555555"


def apply() -> None:
    """rcParams of the report figures."""
    plt.rcParams.update({
        "figure.dpi": 150, "savefig.dpi": 200,
        "font.size": FS, "axes.labelsize": FS, "axes.titlesize": FS,
        "xtick.labelsize": FS_TICK, "ytick.labelsize": FS_TICK,
        "legend.fontsize": FS_LEGEND, "legend.frameon": False,
        "legend.handlelength": 1.8, "legend.columnspacing": 1.2,
        "legend.borderaxespad": 0.3,
        "axes.grid": True, "grid.alpha": 0.3, "grid.linewidth": 0.5,
        "axes.linewidth": 0.6, "xtick.major.width": 0.6,
        "ytick.major.width": 0.6, "xtick.major.size": 2.5,
        "ytick.major.size": 2.5, "xtick.minor.size": 1.5,
        "ytick.minor.size": 1.5,
        "axes.titlepad": 3.0, "axes.labelpad": 2.0,
        "lines.linewidth": 1.0,
        "pdf.fonttype": 42, "ps.fonttype": 42,
        "axes.formatter.limits": (-4, 5),
        "axes.formatter.use_mathtext": True,
    })


def figure(height: float, width: float = WIDTH, **kw):
    """New figure of the print width (height capped at MAX_HEIGHT)."""
    return plt.figure(figsize=(width, min(height, MAX_HEIGHT)), **kw)


def subplots(nrows: int, ncols: int, height: float, width: float = WIDTH,
             **kw):
    kw.setdefault("squeeze", False)
    return plt.subplots(nrows, ncols, figsize=(width, min(height, MAX_HEIGHT)),
                        **kw)


# --------------------------------------------------------------------------- #
# data reduction
# --------------------------------------------------------------------------- #

RAW_MAX = 400        # a history up to this length is drawn as a thin line
NBINS = 250          # bins of an envelope over the full panel width


def _edges(n: int, nbins: int) -> np.ndarray:
    return np.unique(np.linspace(0, n, min(n, nbins) + 1).astype(int))


def decimate(x, y, nmax: int = 1500):
    """Bin means of a smooth curve with at most nmax points (exact below;
    all-NaN bins stay NaN, so gaps survive)."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    if len(x) <= nmax:
        return x, y
    starts = _edges(len(x), nmax)[:-1]

    def binmean(v):
        fin = np.isfinite(v)
        s = np.add.reduceat(np.where(fin, v, 0.0), starts)
        c = np.add.reduceat(fin.astype(float), starts)
        with np.errstate(invalid="ignore", divide="ignore"):
            return np.where(c > 0, s / np.maximum(c, 1.0), np.nan)
    return binmean(x), binmean(y)


def _edges_xlog(x: np.ndarray, nbins: int) -> np.ndarray:
    """Index bin edges that are uniform on a logarithmic x axis (x
    increasing): no bin covers decades at the start of a log axis."""
    pos = x[np.isfinite(x) & (x > 0)]
    if pos.size < 2:
        return _edges(len(x), nbins)
    x0 = max(pos.min(), pos.max() * 1e-5)
    g = np.geomspace(x0, pos.max(), nbins + 1)
    e = np.searchsorted(x, g, side="right")
    return np.unique(np.concatenate([[0], e, [len(x)]]))


def bin_stats(x, y, nbins: int = NBINS, log: bool = False,
              xlog: bool = False):
    """Per-bin (x mean, min, median, max) of y over index bins (xlog: bins
    uniform on a log x axis); y <= 0 is ignored on log axes."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    ok = np.isfinite(y) & np.isfinite(x)
    if log:
        ok &= y > 0
    e = _edges_xlog(x, nbins) if xlog else _edges(len(x), nbins)
    bx, lo, md, hi = [], [], [], []
    for a, b in zip(e[:-1], e[1:]):
        m = ok[a:b]
        if not m.any():
            continue
        yy = y[a:b][m]
        bx.append(x[a:b][m].mean())
        lo.append(yy.min())
        hi.append(yy.max())
        md.append(np.median(yy))
    return np.array(bx), np.array(lo), np.array(md), np.array(hi)


# --------------------------------------------------------------------------- #
# drawing
# --------------------------------------------------------------------------- #

def fill(ax, x, y1, y2, **kw):
    """fill_between that remembers its data, so mask_outside() can redraw
    it without the parts beyond the final y range."""
    coll = ax.fill_between(x, y1, y2, **kw)
    coll._ps_fill = (np.asarray(x, float), np.asarray(y1, float),
                     np.asarray(y2, float), dict(kw))
    return coll


def mask_outside(ax, lo: float, hi: float) -> None:
    """Remove what lies beyond the y range [lo, hi] instead of letting the
    renderer draw it clamped along the axis edge: data-coordinate lines get
    NaN there (the line is broken), bands drawn with fill() lose the parts
    that are entirely outside and are cut at the limits elsewhere. Axis-
    fraction artists (axvline, axhline, clip markers) are untouched."""
    for ln in list(ax.lines):
        if ln.get_transform() != ax.transData:
            continue
        y = np.asarray(ln.get_ydata(orig=True), float)
        if y.ndim != 1 or not y.size:
            continue
        out = (y > hi) | (y < lo)
        if out.any():
            ln.set_ydata(np.where(out, np.nan, y))
    for coll in list(ax.collections):
        data = getattr(coll, "_ps_fill", None)
        if data is None:
            continue
        x, y1, y2, kw = data
        a, b = np.minimum(y1, y2), np.maximum(y1, y2)
        gone = (a > hi) | (b < lo)
        if not gone.any() and a.min(initial=lo) >= lo and b.max(initial=hi) <= hi:
            continue
        a = np.where(gone, np.nan, np.clip(a, lo, hi))
        b = np.where(gone, np.nan, np.clip(b, lo, hi))
        coll.remove()
        fill(ax, x, a, b, **kw)


def envelope(ax, x, y, color, nbins: int = NBINS, log: bool = False,
             label=None, band_label=None, alpha: float = BAND_ALPHA,
             lw: float = LW_MAIN, ls: str = "-", z: float = Z_ENV,
             line: bool = True, xlog: bool = False):
    """Long or noisy history: per-bin min-max band (light) and the bin
    median (bold). Returns (bin x, bin median)."""
    bx, lo, md, hi = bin_stats(x, y, nbins, log, xlog)
    if not len(bx):
        return bx, md
    fill(ax, bx, lo, hi, color=color, alpha=alpha, lw=0, zorder=z,
         label=band_label)
    if line:
        ax.plot(bx, md, color=color, lw=lw, ls=ls, zorder=z + 0.2,
                label=label)
    return bx, md


def trace(ax, x, y, color, label=None, log: bool = False,
          lw: float = LW_RAW, alpha: float = 0.9, z: float = Z_RAW,
          nbins: int = NBINS, raw_max: int = RAW_MAX, ls: str = "-",
          band_alpha: float = BAND_ALPHA, median_lw: float | None = None,
          line: bool = True, xlog: bool = False):
    """One per-iteration history: a thin line up to raw_max points,
    otherwise the envelope (band + bin median; line=False: band only, when
    the caller draws its own smooth curve; xlog: bins uniform on a log x
    axis)."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    if len(x) <= raw_max:
        yy = np.where(y > 0, y, np.nan) if log else y
        ax.plot(x, yy, color=color, lw=lw, alpha=alpha, zorder=z, label=label,
                ls=ls)
        return
    envelope(ax, x, y, color, nbins=nbins, log=log, label=label,
             band_label=None if line else label, alpha=band_alpha,
             lw=median_lw or max(lw, 1.0), z=z - 0.5, ls=ls, line=line,
             xlog=xlog)


def robust_ylim(ax, series, skip_frac: float = 0.2, pct=(1.0, 99.0),
                margin: float = 0.15, log: bool = False, include=(),
                mark: bool = True, min_span: float = 0.0):
    """Set y limits from the post-transient part of every series.

    series: iterable of (x, y, colour) (colour None: no clip markers).
    The first skip_frac of each series (the start-up transient) and the
    outer (100 - pct) percentiles are ignored; values in `include` (window
    means, reference lines) are always inside. Points beyond the limits
    are shown as small triangles at the axis edge in the series colour;
    everything already drawn in data coordinates beyond the limits is
    masked (mask_outside), so call this after drawing the panel.
    Returns (lo, hi)."""
    vals = []
    for x, y, _ in series:
        y = np.asarray(y, float)
        k = int(len(y) * skip_frac)
        v = y[k:]
        v = v[np.isfinite(v)]
        if log:
            v = v[v > 0]
            v = np.log10(v)
        if v.size:
            vals.append(v)
    inc = np.asarray([v for v in include if v is not None and np.isfinite(v)
                      and (not log or v > 0)], float)
    if log and inc.size:
        inc = np.log10(inc)
    if not vals and not inc.size:
        return ax.get_ylim()
    allv = np.concatenate(vals) if vals else inc
    lo, hi = np.percentile(allv, pct)
    if inc.size:
        lo, hi = min(lo, inc.min()), max(hi, inc.max())
    span = max(hi - lo, min_span, 1e-12 * max(abs(lo), abs(hi), 1.0))
    lo, hi = lo - margin * span, hi + margin * span
    if log:
        lo, hi = 10 ** lo, 10 ** hi
    ax.set_ylim(lo, hi)
    mask_outside(ax, lo, hi)
    if mark:
        for x, y, col in series:
            if col is not None:
                clip_markers(ax, x, y, col, lo, hi)
    return lo, hi


def clip_markers(ax, x, y, color, lo, hi, nmax: int = 25) -> None:
    """Small triangles at the top/bottom axis edge where a series leaves the
    y range (the clipped start-up transient stays visible as a hint)."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    trans = ax.get_xaxis_transform()
    for sel, ypos, mk in ((y > hi, 0.985, "^"), (y < lo, 0.015, "v")):
        idx = np.flatnonzero(sel & np.isfinite(x))
        if not idx.size:
            continue
        # sparse markers: at least 4 % of the panel's data width apart (a
        # dense row of triangles would read as a line along the edge)
        xs = np.sort(x[idx])
        x0, x1 = ax.dataLim.intervalx
        dx = 0.04 * (x1 - x0) if np.isfinite(x1 - x0) and x1 > x0 else 0.0
        keep, last = [], -np.inf
        for v in xs:
            if v - last >= dx:
                keep.append(v)
                last = v
        xs = np.array(keep[:nmax])
        ax.plot(xs, np.full(xs.size, ypos), mk, color=color, ms=2.8,
                mec="none", alpha=0.8, transform=trans, clip_on=False,
                zorder=Z_MARK)


def zoom_inset(ax, xlim, bounds=(0.40, 0.40, 0.57, 0.55), label=None,
               ylog: bool = False):
    """Inset axes on ax for a short run on a long common axis; returns the
    inset (the caller draws into it and sets its y range). A light
    rectangle marks the zoomed range on the parent axis."""
    ins = ax.inset_axes(bounds)
    ins.set_xlim(*xlim)
    if ylog:
        ins.set_yscale("log")
    ins.tick_params(labelsize=FS_TICK - 0.5, pad=1.5)
    ins.grid(True, alpha=0.3)
    ins.set_facecolor("white")
    for s in ins.spines.values():
        s.set_edgecolor(GREY)
        s.set_linewidth(0.6)
    ax.axvspan(xlim[0], xlim[1], color=GREY, alpha=0.08, lw=0, zorder=Z_SHADE)
    if label:
        ins.set_title(label, fontsize=FS_NOTE, pad=2)
    return ins


def note(ax, text, loc="upper left", color=GREY):
    """Short in-panel annotation."""
    pos = {"upper left": (0.02, 0.96, "left", "top"),
           "upper right": (0.98, 0.96, "right", "top"),
           "lower left": (0.02, 0.04, "left", "bottom"),
           "lower right": (0.98, 0.04, "right", "bottom"),
           "center": (0.5, 0.5, "center", "center")}[loc]
    ax.text(pos[0], pos[1], text, transform=ax.transAxes, ha=pos[2],
            va=pos[3], fontsize=FS_NOTE, color=color, zorder=Z_MARK + 1,
            bbox=dict(boxstyle="square,pad=0.15", fc="white", ec="none",
                      alpha=0.75))


def unique_handles(axes):
    h, lab = [], []
    for ax in axes:
        for hh, ll in zip(*ax.get_legend_handles_labels()):
            if ll and not ll.startswith("_") and ll not in lab:
                h.append(hh)
                lab.append(ll)
    return h, lab


def legend_below(fig, handles=None, labels=None, ncol: int | None = None,
                 fontsize: float = FS_LEGEND, h_pad: float | None = None,
                 w_pad: float | None = None, top: float = 1.0):
    """One legend under all panels, then tight_layout above it (the legend
    height is measured, so nothing overlaps). Returns the legend."""
    if handles is None:
        handles, labels = unique_handles(fig.axes)
    kw = {}
    if h_pad is not None:
        kw["h_pad"] = h_pad
    if w_pad is not None:
        kw["w_pad"] = w_pad
    if not handles:
        fig.tight_layout(rect=(0, 0, 1, top), **kw)
        return None
    if ncol is None:
        ncol = min(len(handles), 3 if max(len(x) for x in labels) > 28 else 4)
    ncol = max(1, min(ncol, len(handles)))
    width = fig.get_figwidth() * fig.dpi
    while True:
        leg = fig.legend(handles, labels, loc="lower center",
                         bbox_to_anchor=(0.5, 0.0), ncol=ncol,
                         fontsize=fontsize, frameon=False)
        fig.canvas.draw()
        bb = leg.get_window_extent()
        if bb.width <= width - 4 or ncol == 1:
            break
        leg.remove()          # too wide for the page: one column less
        ncol -= 1
    frac = (bb.height + 4) / (fig.get_figheight() * fig.dpi)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        fig.tight_layout(rect=(0, frac, 1, top), **kw)
    return leg


def inner_legend(ax, ncol: int = 1, loc: str = "best"):
    """Legend inside a small panel: compact, on a semi-opaque white box so
    it stays readable where it has to cover data."""
    h, lab = ax.get_legend_handles_labels()
    if not h:
        return None
    return ax.legend(h, lab, loc=loc, ncol=ncol, fontsize=FS_LEGEND - 0.5,
                     frameon=True, framealpha=0.85, edgecolor="none",
                     facecolor="white", handlelength=1.3, borderpad=0.25,
                     labelspacing=0.2, columnspacing=0.8, handletextpad=0.4)


def tight(fig, **kw) -> None:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        fig.tight_layout(**kw)


def save(fig, pdf: Path, png: Path | None = None) -> None:
    """Fixed-size vector PDF (and PNG): the figure size is the print size,
    so no bbox cropping."""
    pdf = Path(pdf)
    fig.savefig(pdf.with_suffix(".pdf"))
    if png is not None:
        fig.savefig(Path(png).with_suffix(".png"))
    plt.close(fig)


def is_constant(y) -> bool:
    y = np.asarray(y, float)
    y = y[np.isfinite(y)]
    return y.size > 0 and np.nanmax(y) == np.nanmin(y)

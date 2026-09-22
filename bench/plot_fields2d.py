#!/usr/bin/env python3
"""Flow-field comparison figures for the 2D test cases (matplotlib).

For every 2D case with BOTH a coupledFoam run (run/<case>_np1) and a
simpleFoam reference on the same mesh (run/ref_<case>), this writes vector
PDFs into report/paper/figures/:

    fields_<case>_U.pdf    |U|/U_ref: coupledFoam, simpleFoam (same colour
                           scale, streamlines) and the delta
                           (|U|_cf - |U|_sf)/U_ref, diverging colour map
    fields_<case>_p.pdf    the same for C_p = (p - p_0)/(0.5 U_ref^2)
    profiles_<case>.pdf    case-specific line comparisons (T1/T2 velocity
                           profiles, T2 lower-wall skin friction, T3 surface
                           C_p); T0 has its centreline profiles already
                           (bench/make_report.py:fig_T0_profiles)

How the fields are read (documented choice): the VTK python module of the
venv (vtk 9.2, /home/jonas/OF/venv) contains vtkOpenFOAMReader, which reads
the binary OpenFOAM fields of the case directly, read-only - no foamToVTK,
no postProcess, no copies under run/. The reader interpolates the cell
values to the mesh points (boundary values included); the single-layer 2D
mesh is cut at mid-depth (vtkCutter) and triangulated, so the plots use the
real mesh geometry (holes such as the airfoil stay holes) and matplotlib's
tricontourf on the nodal values. Both solvers use the same mesh; if the
node sets differ, the simpleFoam field is interpolated linearly onto the
coupledFoam nodes. Streamlines are integrated on a regular grid
interpolated from the triangulation (points outside the mesh are masked).

The delta statistics (area-weighted RMS and max of the deltas) are returned
to bench/make_report.py, which quotes them in the paper, and are written to
report/paper/figures/fields2d_stats.json.

Usage: /home/jonas/OF/venv/bin/python bench/plot_fields2d.py [case ...]
       (no OpenFOAM environment needed; runs serially, nice it)
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.tri as mtri  # noqa: E402
import numpy as np  # noqa: E402

REPO = Path(__file__).resolve().parents[1]
RUN = REPO / "run"
FIG = REPO / "report" / "paper" / "figures"

C_NATIVE = "#7f7f7f"
C_COUPLED = "#1f5fbf"
CMAP_U = "viridis"
CMAP_P = "cividis"
CMAP_D = "RdBu_r"

# name: coupledFoam dir, simpleFoam dir, U_ref [m/s], pressure level,
# view (None: whole mesh; "patch:<name>": around that patch), extras
CASES = {
    "T0_Re100": dict(cf="T0_Re100_np1", sf="ref_T0_Re100", uref=1.0,
                     plevel="mean", view=None, L=0.1,
                     title="T0 lid-driven cavity, Re 100"),
    "T0_Re1000": dict(cf="T0_Re1000_np1", sf="ref_T0_Re1000", uref=1.0,
                      plevel="mean", view=None, L=0.1,
                      title="T0 lid-driven cavity, Re 1000"),
    "T1": dict(cf="T1_np1", sf="ref_T1", uref=10.0, plevel="outlet",
               view=None, h=0.0254, title="T1 pitzDaily",
               stations=(1, 3, 5, 7, 9, 11)),
    "T2": dict(cf="T2_np1", sf="ref_T2", uref=44.2, plevel="outlet",
               view=(-0.05, 0.30, None, None), h=0.0127,
               title="T2 backward-facing step",
               stations=(2, 4, 6, 8, 10, 14, 18), wallshear="lowerWallShear"),
    "T3_kOmegaSST": dict(cf="T3_kOmegaSST_np1", sf="ref_T3_kOmegaSST",
                         uref=26.0032, plevel="outlet", view="patch:walls",
                         title="T3 airFoil2D, k-omega SST", airfoil="walls"),
    "T3_GEKO": dict(cf="T3_GEKO_np1", sf="ref_T3_GEKO", uref=26.0032,
                    plevel="outlet", view="patch:walls",
                    title="T3 airFoil2D, GEKO", airfoil="walls"),
}

plt.rcParams.update({
    "figure.dpi": 150, "savefig.dpi": 300, "font.size": 8,
    "axes.grid": False, "legend.frameon": False,
    "pdf.fonttype": 42, "ps.fonttype": 42,
})


# --------------------------------------------------------------------------- #
# reading
# --------------------------------------------------------------------------- #

def _vtk():
    import vtk  # noqa: PLC0415
    vtk.vtkLogger.SetStderrVerbosity(vtk.vtkLogger.VERBOSITY_OFF)
    vtk.vtkObject.GlobalWarningDisplayOff()
    return vtk


def _blocks(ds, prefix=""):
    """(name, dataset) of all leaves of a vtkMultiBlockDataSet."""
    out = []
    for i in range(ds.GetNumberOfBlocks()):
        b = ds.GetBlock(i)
        md = ds.GetMetaData(i)
        name = md.Get(ds.NAME()) if md is not None and md.Has(ds.NAME()) else str(i)
        if b is None:
            continue
        if b.IsA("vtkMultiBlockDataSet"):
            out += _blocks(b, prefix + name + "/")
        else:
            out.append((prefix + name, b))
    return out


def latest_time(case: Path) -> float | None:
    ts = []
    for d in case.iterdir():
        try:
            t = float(d.name)
        except ValueError:
            continue
        if t > 0 and (d / "U").exists() and (d / "p").exists():
            ts.append(t)
    return max(ts) if ts else None


def read_case(case: Path, patches: tuple[str, ...] = ()) -> dict:
    """Mid-depth slice of the internal mesh (triangulated, nodal U and p)
    and mid-depth cuts of the listed patches (ordered polylines, nodal p)
    at the latest time of an OpenFOAM case."""
    vtk = _vtk()
    from vtk.util.numpy_support import vtk_to_numpy  # noqa: PLC0415

    t = latest_time(case)
    if t is None:
        raise FileNotFoundError(f"{case}: no time directory with U and p")
    r = vtk.vtkOpenFOAMReader()
    r.SetFileName(str(case / "system" / "controlDict"))
    r.UpdateInformation()
    r.DisableAllCellArrays()
    r.DisableAllPointArrays()
    for f in ("U", "p"):
        r.SetCellArrayStatus(f, 1)
    r.DisableAllPatchArrays()
    r.SetPatchArrayStatus("internalMesh", 1)
    for p in patches:
        r.SetPatchArrayStatus("patch/" + p, 1)
    r.SetCreateCellToPoint(1)
    r.UpdateTimeStep(t)
    blocks = dict(_blocks(r.GetOutput()))
    mesh = blocks.get("internalMesh")
    if mesh is None:
        raise RuntimeError(f"{case}: no internalMesh block")
    b = mesh.GetBounds()
    zmid = 0.5 * (b[4] + b[5])
    plane = vtk.vtkPlane()
    plane.SetOrigin(0, 0, zmid)
    plane.SetNormal(0, 0, 1)

    def cut(ds):
        c = vtk.vtkCutter()
        c.SetCutFunction(plane)
        c.SetInputData(ds)
        c.Update()
        return c.GetOutput()

    tri = vtk.vtkTriangleFilter()
    tri.SetInputData(cut(mesh))
    tri.Update()
    sl = tri.GetOutput()
    pts = vtk_to_numpy(sl.GetPoints().GetData()).astype(float)
    conn = vtk_to_numpy(sl.GetPolys().GetConnectivityArray()).reshape(-1, 3)
    # drop degenerate triangles of the cut (zero area), they make holes in
    # tricontourf
    x, y = pts[:, 0], pts[:, 1]
    area = 0.5 * np.abs((x[conn[:, 1]] - x[conn[:, 0]]) * (y[conn[:, 2]] - y[conn[:, 0]])
                        - (x[conn[:, 2]] - x[conn[:, 0]]) * (y[conn[:, 1]] - y[conn[:, 0]]))
    conn = conn[area > 1e-10 * np.median(area)]
    pd = sl.GetPointData()
    out = {
        "time": t, "x": pts[:, 0], "y": pts[:, 1], "tri": conn,
        "U": vtk_to_numpy(pd.GetArray("U")).astype(float),
        "p": vtk_to_numpy(pd.GetArray("p")).astype(float),
        "bounds": b, "nCells": mesh.GetNumberOfCells(), "patches": {},
    }
    for p in patches:
        blk = next((d for n, d in blocks.items() if n.endswith("/" + p)), None)
        if blk is None:
            continue
        strip = vtk.vtkStripper()
        strip.SetInputData(cut(blk))
        strip.JoinContiguousSegmentsOn()
        strip.Update()
        line = strip.GetOutput()
        lp = vtk_to_numpy(line.GetPoints().GetData()).astype(float)
        lpd = vtk_to_numpy(line.GetPointData().GetArray("p")).astype(float)
        conn_l = vtk_to_numpy(line.GetLines().GetConnectivityArray())
        offs = vtk_to_numpy(line.GetLines().GetOffsetsArray())
        polys = [conn_l[offs[i]:offs[i + 1]] for i in range(len(offs) - 1)]
        out["patches"][p] = {"x": lp[:, 0], "y": lp[:, 1], "p": lpd,
                             "polylines": polys,
                             "bounds": blk.GetBounds()}
    return out


def area_weights(x, y, tri) -> np.ndarray:
    """Nodal weights (one third of the adjacent triangle areas)."""
    a = 0.5 * np.abs((x[tri[:, 1]] - x[tri[:, 0]]) * (y[tri[:, 2]] - y[tri[:, 0]])
                     - (x[tri[:, 2]] - x[tri[:, 0]]) * (y[tri[:, 1]] - y[tri[:, 0]]))
    w = np.zeros_like(x)
    for k in range(3):
        np.add.at(w, tri[:, k], a / 3.0)
    return w


def on_same_nodes(cf: dict, sf: dict, key: str) -> np.ndarray:
    """simpleFoam nodal values of `key` on the coupledFoam nodes."""
    if len(cf["x"]) == len(sf["x"]):
        d = np.hypot(cf["x"] - sf["x"], cf["y"] - sf["y"]).max()
        scale = max(np.ptp(cf["x"]), np.ptp(cf["y"]))
        if d < 1e-9 * scale:
            return sf[key]
    trs = mtri.Triangulation(sf["x"], sf["y"], sf["tri"])
    v = sf[key]
    cols = [v] if v.ndim == 1 else [v[:, k] for k in range(v.shape[1])]
    res = [mtri.LinearTriInterpolator(trs, c)(cf["x"], cf["y"]).filled(np.nan)
           for c in cols]
    return res[0] if v.ndim == 1 else np.stack(res, axis=1)


# --------------------------------------------------------------------------- #
# plotting
# --------------------------------------------------------------------------- #

def _view(cfg: dict, d: dict) -> tuple[float, float, float, float]:
    b = d["bounds"]
    v = cfg.get("view")
    if isinstance(v, str) and v.startswith("patch:"):
        pb = d["patches"][v.split(":", 1)[1]]["bounds"]
        c = pb[1] - pb[0]
        return (pb[0] - 0.35 * c, pb[1] + 0.9 * c,
                0.5 * (pb[2] + pb[3]) - 0.45 * c, 0.5 * (pb[2] + pb[3]) + 0.45 * c)
    if isinstance(v, tuple):
        return tuple(b[i] if v[i] is None else v[i] for i in range(4))
    return (b[0], b[1], b[2], b[3])


def _stream(ax, tr, U, view, uref, color):
    """Streamlines of the in-plane velocity on a regular grid over view."""
    x0, x1, y0, y1 = view
    asp = (x1 - x0) / (y1 - y0)
    nx = 260 if asp >= 1 else max(40, int(260 * asp))
    ny = max(40, int(nx / asp)) if asp >= 1 else 260
    gx, gy = np.meshgrid(np.linspace(x0, x1, nx), np.linspace(y0, y1, ny))
    gu = mtri.LinearTriInterpolator(tr, U[:, 0])(gx, gy)
    gv = mtri.LinearTriInterpolator(tr, U[:, 1])(gx, gy)
    gu = np.ma.filled(gu, np.nan)
    gv = np.ma.filled(gv, np.nan)
    try:
        ax.streamplot(gx, gy, gu, gv, color=color, linewidth=0.35,
                      density=(1.4 * min(asp, 3.5), 1.4) if asp >= 1 else 1.4,
                      arrowsize=0.45, broken_streamlines=True)
    except (ValueError, IndexError):
        pass


def _levels(vals: list[np.ndarray], lo=0.5, hi=99.5, sym=False, zero_min=False):
    a = np.concatenate([v[np.isfinite(v)] for v in vals])
    if sym:
        m = np.percentile(np.abs(a), hi)
        m = m if m > 0 else 1e-12
        return -m, m
    vmin = 0.0 if zero_min else np.percentile(a, lo)
    vmax = np.percentile(a, hi)
    if vmax <= vmin:
        vmax = vmin + 1e-12
    return vmin, vmax


def _ticks(cb, n=5):
    from matplotlib.ticker import MaxNLocator  # noqa: PLC0415
    cb.locator = MaxNLocator(n, symmetric=False)
    cb.formatter.set_powerlimits((-2, 3))
    cb.update_ticks()


def _panel(ax, tr, v, vmin, vmax, cmap, view):
    lev = np.linspace(vmin, vmax, 33)
    # no clipping: values beyond the range go to the extend colours (a
    # clipped flat region exactly at a level boundary leaves holes)
    cs = ax.tricontourf(tr, np.nan_to_num(v, nan=0.5 * (vmin + vmax)),
                        levels=lev, cmap=cmap, extend="both")
    ax.set_xlim(view[0], view[1])
    ax.set_ylim(view[2], view[3])
    ax.set_aspect("equal")
    ax.tick_params(labelsize=6, length=2)
    return cs


def field_figure(name, cfg, cf, sf, quantity, stats, outdir) -> tuple[str, str]:
    tr = mtri.Triangulation(cf["x"], cf["y"], cf["tri"])
    uref = cfg["uref"]
    q = 0.5 * uref ** 2
    view = _view(cfg, cf)
    Usf = on_same_nodes(cf, sf, "U")
    psf = on_same_nodes(cf, sf, "p")
    w = area_weights(cf["x"], cf["y"], cf["tri"])
    if quantity == "U":
        a = np.linalg.norm(cf["U"][:, :2], axis=1) / uref
        b = np.linalg.norm(Usf[:, :2], axis=1) / uref
        label, dlabel = r"$|U|/U_{ref}$", r"$(|U|_{cf}-|U|_{sf})/U_{ref}$"
        vmin, vmax = _levels([a, b], zero_min=True)
        cmap = CMAP_U
    else:
        pa, pb = cf["p"].copy(), psf.copy()
        if cfg["plevel"] == "mean":
            ok = np.isfinite(pb)
            pa -= np.average(pa[ok], weights=w[ok])
            pb -= np.average(pb[ok], weights=w[ok])
        a, b = pa / q, pb / q
        label, dlabel = r"$C_p$", r"$C_{p,cf}-C_{p,sf}$"
        vmin, vmax = _levels([a, b], lo=1.0, hi=99.0)
        cmap = CMAP_P
    d = a - b
    ok = np.isfinite(d)
    rms = float(np.sqrt(np.average(d[ok] ** 2, weights=w[ok])))
    dmax = float(np.abs(d[ok]).max())
    stats[f"{name}_{quantity}_rms"] = rms
    stats[f"{name}_{quantity}_max"] = dmax
    dlim = _levels([d], hi=99.5, sym=True)
    dlim = (min(dlim[0], -1e-6), max(dlim[1], 1e-6))

    x0, x1, y0, y1 = view
    asp = (x1 - x0) / (y1 - y0)
    if asp > 1.8:
        fig, axs = plt.subplots(3, 1, figsize=(6.5, min(8.5, 3 * 6.0 / asp + 1.2)),
                                sharex=True)
    else:
        fig, axs = plt.subplots(1, 3, figsize=(6.5, 6.5 / 3 / asp + 0.9),
                                sharey=True)
    cs = _panel(axs[0], tr, a, vmin, vmax, cmap, view)
    _panel(axs[1], tr, b, vmin, vmax, cmap, view)
    ds = _panel(axs[2], tr, d, dlim[0], dlim[1], CMAP_D, view)
    if quantity == "U":
        trb = tr
        _stream(axs[0], trb, cf["U"], view, uref, "white")
        _stream(axs[1], trb, np.nan_to_num(Usf), view, uref, "white")
    axs[0].set_title("coupledFoam", fontsize=7, color=C_COUPLED)
    axs[1].set_title("simpleFoam", fontsize=7, color=C_NATIVE)
    axs[2].set_title(f"delta: RMS {rms:.2e}, max {dmax:.2e}", fontsize=7)
    if asp > 1.8:
        for ax in axs:
            ax.set_ylabel("y [m]", fontsize=6)
        axs[2].set_xlabel("x [m]", fontsize=6)
        for ax, m, lab in ((axs[0], cs, label), (axs[1], cs, label),
                           (axs[2], ds, dlabel)):
            cb = fig.colorbar(m, ax=ax, fraction=0.025, pad=0.01)
            _ticks(cb)
            cb.set_label(lab, fontsize=6)
            cb.ax.tick_params(labelsize=5)
    else:
        axs[0].set_ylabel("y [m]", fontsize=6)
        for ax in axs:
            ax.set_xlabel("x [m]", fontsize=6)
        cb = fig.colorbar(cs, ax=list(axs[:2]), orientation="horizontal",
                          fraction=0.06, pad=0.18, aspect=40)
        _ticks(cb)
        cb.set_label(label + " (same scale for both solvers)", fontsize=6)
        cb.ax.tick_params(labelsize=5)
        cb2 = fig.colorbar(ds, ax=axs[2], orientation="horizontal",
                           fraction=0.06, pad=0.18, aspect=20)
        _ticks(cb2, 3)
        cb2.set_label(dlabel, fontsize=6)
        cb2.ax.tick_params(labelsize=5)
    fig.suptitle(f"{cfg['title']}: {'velocity magnitude' if quantity == 'U' else 'pressure coefficient'}",
                 fontsize=8)
    stem = f"fields_{name}_{quantity}"
    fig.savefig(outdir / f"{stem}.pdf", bbox_inches="tight")
    plt.close(fig)
    what = ("velocity magnitude $|U|/U_{ref}$ with streamlines"
            if quantity == "U" else
            "pressure coefficient $C_p=(p-p_0)/(\\tfrac12 U_{ref}^2)$")
    cap = (f"{cfg['title']}: {what}; left/top coupledFoam, middle simpleFoam "
           f"(identical colour scale), right/bottom the difference "
           f"coupledFoam$-$simpleFoam (area-weighted RMS {rms:.1e}, "
           f"max {dmax:.1e}).")
    return stem, cap


def profiles_figure(name, cfg, cf, sf, outdir) -> tuple[str, str] | None:
    """Velocity profiles at x/h stations, skin friction, surface C_p."""
    panels = []
    if cfg.get("stations"):
        panels.append("stations")
    if cfg.get("wallshear"):
        panels.append("wallshear")
    if cfg.get("airfoil"):
        panels.append("airfoil")
    if not panels:
        return None
    uref = cfg["uref"]
    fig, axs = plt.subplots(len(panels), 1,
                            figsize=(6.5, 2.3 * len(panels)), squeeze=False)
    axs = axs[:, 0]
    trc = mtri.Triangulation(cf["x"], cf["y"], cf["tri"])
    trs = mtri.Triangulation(sf["x"], sf["y"], sf["tri"])
    for ax, kind in zip(axs, panels):
        if kind == "stations":
            h = cfg["h"]
            b = cf["bounds"]
            ys = np.linspace(b[2], b[3], 300)
            sc = 0.9 * float(np.min(np.diff(cfg["stations"])))
            for i, s in enumerate(cfg["stations"]):
                xs = np.full_like(ys, s * h)
                for tr, U, col, lw, ls, lab in (
                        (trs, sf["U"], C_NATIVE, 2.2, "-", "simpleFoam"),
                        (trc, cf["U"], C_COUPLED, 0.9, "--", "coupledFoam")):
                    u = mtri.LinearTriInterpolator(tr, U[:, 0])(xs, ys)
                    ax.plot(s + sc * np.ma.filled(u, np.nan) / uref, ys / h,
                            color=col, lw=lw, ls=ls,
                            label=lab if i == 0 else None)
                ax.axvline(s, color="k", lw=0.3, alpha=0.4)
            ax.set_xlabel(f"$x/h$ (station) $+ {sc}\\,u/U_{{ref}}$")
            ax.set_ylabel("$y/h$")
            ax.legend(fontsize=7, loc="upper right")
            ax.set_title("streamwise velocity profiles downstream of the step",
                         fontsize=8)
        elif kind == "wallshear":
            h = cfg["h"]
            q = 0.5 * uref ** 2
            for case, col, lw, ls, lab in ((RUN / cfg["sf"], C_NATIVE, 2.2, "-", "simpleFoam"),
                                           (RUN / cfg["cf"], C_COUPLED, 0.9, "--", "coupledFoam")):
                raw = _wallshear_raw(case, cfg["wallshear"])
                if raw is None:
                    continue
                o = np.argsort(raw[:, 0])
                x, tx = raw[o, 0], raw[o, 3]
                ax.plot(x / h, -tx / q, color=col, lw=lw, ls=ls, label=lab)
                # wallShearStress is the force on the wall per area with the
                # sign convention tau = -nu*dU/dn: C_f = -tau_x/q on the floor
            ax.axhline(0, color="k", lw=0.4)
            ax.set_xlim(-2, 25)
            ax.set_xlabel("$x/h$")
            ax.set_ylabel("$C_f$ (lower wall)")
            ax.legend(fontsize=7)
            ax.set_title("skin friction on the lower wall: the zero crossing "
                         "is the reattachment point", fontsize=8)
        elif kind == "airfoil":
            q = 0.5 * uref ** 2
            pb = cf["patches"][cfg["airfoil"]]["bounds"]
            c = pb[1] - pb[0]
            for d, col, lw, ls, lab in ((sf, C_NATIVE, 2.2, "-", "simpleFoam"),
                                        (cf, C_COUPLED, 0.9, "--", "coupledFoam")):
                pt = d["patches"].get(cfg["airfoil"])
                if pt is None:
                    continue
                for k, poly in enumerate(pt["polylines"]):
                    ax.plot((pt["x"][poly] - pb[0]) / c, -pt["p"][poly] / q,
                            color=col, lw=lw, ls=ls,
                            label=lab if k == 0 else None)
            ax.set_xlabel("$x/c$")
            ax.set_ylabel("$-C_p$")
            ax.legend(fontsize=7)
            ax.set_title("surface pressure coefficient on the airfoil "
                         "(upper curve: suction side)", fontsize=8)
        ax.grid(True, alpha=0.3)
    fig.suptitle(cfg["title"], fontsize=8)
    fig.tight_layout()
    stem = f"profiles_{name}"
    fig.savefig(outdir / f"{stem}.pdf", bbox_inches="tight")
    plt.close(fig)
    parts = {"stations": "streamwise velocity profiles at stations $x/h$",
             "wallshear": "lower-wall skin friction $C_f$",
             "airfoil": "surface pressure coefficient"}
    cap = (f"{cfg['title']}: " + ", ".join(parts[k] for k in panels)
           + "; simpleFoam thick grey, coupledFoam dashed blue.")
    return stem, cap


def _wallshear_raw(case: Path, fo: str) -> np.ndarray | None:
    d = case / "postProcessing" / fo
    if not d.is_dir():
        return None
    times = sorted((p for p in d.iterdir() if p.is_dir()),
                   key=lambda p: float(p.name) if p.name.replace(".", "").isdigit() else -1)
    for t in reversed(times):
        f = next(iter(sorted(t.glob("*.raw"))), None)
        if f is not None:
            try:
                return np.loadtxt(f, comments="#")
            except ValueError:
                return None
    return None


# --------------------------------------------------------------------------- #
# driver
# --------------------------------------------------------------------------- #

def run(names: list[str] | None = None, outdir: Path = FIG,
        log=print) -> tuple[list[tuple[str, str]], dict, list[str]]:
    """Make all figures. Returns (figures [(stem, caption)], stats, notes)."""
    outdir.mkdir(parents=True, exist_ok=True)
    figs, stats, notes = [], {}, []
    try:
        _vtk()
    except ImportError:
        notes.append("fields 2D: python module vtk missing, skipped")
        return figs, stats, notes
    for name, cfg in CASES.items():
        if names and name not in names:
            continue
        cfd, sfd = RUN / cfg["cf"], RUN / cfg["sf"]
        if not (cfd.is_dir() and sfd.is_dir()):
            notes.append(f"fields {name}: run/{cfg['cf']} or run/{cfg['sf']} missing")
            continue
        patches = tuple(p for p in (cfg.get("airfoil"),
                                    (cfg.get("view") or "")[6:]
                                    if isinstance(cfg.get("view"), str) else None)
                        if p)
        try:
            cf = read_case(cfd, patches)
            sf = read_case(sfd, patches)
        except (FileNotFoundError, RuntimeError, AttributeError) as e:
            notes.append(f"fields {name}: {e}")
            continue
        log(f"fields {name}: cf t={cf['time']:g} sf t={sf['time']:g} "
            f"nodes {len(cf['x'])}/{len(sf['x'])}")
        stats[f"{name}_cfTime"] = cf["time"]
        stats[f"{name}_sfTime"] = sf["time"]
        for qn in ("U", "p"):
            figs.append(field_figure(name, cfg, cf, sf, qn, stats, outdir))
        pr = profiles_figure(name, cfg, cf, sf, outdir)
        if pr:
            figs.append(pr)
    (outdir / "fields2d_stats.json").write_text(json.dumps(stats, indent=1))
    return figs, stats, notes


if __name__ == "__main__":
    f, s, n = run(sys.argv[1:] or None)
    for stem, _ in f:
        print("wrote", stem)
    for x in n:
        print("note:", x)

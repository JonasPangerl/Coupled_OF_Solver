#!/usr/bin/env pvbatch
"""ParaView renders of the 3D cases (T4a, T4b, T5) for the paper.

Run with ParaView's batch interpreter (5.11, offscreen):

    nice -n 19 pvbatch --force-offscreen-rendering bench/render_fields.py \
        [--cases T4a T4b T5] [--out report/paper/figures] [--busy-minutes 10]
        [--reference-only]

Inputs (read-only, decomposed cases read with the OpenFOAM reader in
decomposed mode, latest time; nothing is reconstructed or written under
run/):

    run/ref_<case>_np10   simpleFoam reference
    run/<case>_np10       coupledFoam run; after the heavy runs also UMean,
                          pMean and the delta fields UMeanDelta, pMeanDelta,
                          CpMeanDelta, magUMeanDeltaRel written by
                          applications/utilities/coupledFieldCompare

Outputs (PNG, 300 dpi at the paper's text width, report/paper/figures/):

    render_<C>_geometry.png    body surface with the mesh edges
    render_<C>_slice_U.png     |U|/U_inf on the mid-plane (y = y_c) and on
                               the wheel-height plane, simpleFoam | coupledFoam
    render_<C>_slice_p.png     the same for C_p = p / (0.5 U_inf^2)
    render_<C>_surface_Cp.png  C_p on the body, simpleFoam | coupledFoam
    render_<C>_mean_U.png      |UMean|/U_inf mid-plane, both solvers (if the
                               mean fields exist in both cases)
    render_<C>_delta.png       magUMeanDeltaRel on both planes and CpMeanDelta
                               on the body (if coupledFieldCompare has run)

Missing inputs are skipped with a note; a case whose logs (simpleFoam or coupledFoam
directory) were modified within --busy-minutes is treated as running and
skipped entirely (never read a case
while a solver writes it). A missing coupledFoam panel is drawn as a grey
"pending" panel, so the figure layout stays the same when the data arrive.

Reproducibility: camera positions, slice planes and colour ranges are
fixed below (VIEW, RANGES) and relative to the body bounding box; the same
colour range is used for both solvers. Nothing depends on the data range.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import sys
import time
from pathlib import Path

from paraview.simple import (  # noqa: E402
    AssignViewToLayout, Calculator, ColorBy, CreateLayout, CreateView,
    GetColorTransferFunction, GetOpacityTransferFunction, GetScalarBar,
    OpenFOAMReader,
    RemoveLayout, SaveScreenshot, Show, Slice, Text, Delete)

from vtkmodules.vtkCommonCore import vtkLogger  # noqa: E402

vtkLogger.SetStderrVerbosity(vtkLogger.VERBOSITY_ERROR)

REPO = Path(__file__).resolve().parents[1]
RUN = REPO / "run"

CASES = {
    "T4a": dict(cf="T4a_np10", sf="ref_T4a_np10", title="T4a motorBike (354k cells)"),
    "T4b": dict(cf="T4b_np10", sf="ref_T4b_np10", title="T4b motorBike (1.70M cells)"),
    "T5": dict(cf="T5_np10", sf="ref_T5_np10", title="T5 Ahmed body"),
}

# Colour ranges (fixed, identical for both solvers)
RANGES = {
    "magUrel": (0.0, 1.3),        # |U|/U_inf
    "Cp": (-1.5, 1.0),            # p/(0.5 U_inf^2), p kinematic, p_inf = 0
    "magUMeanDeltaRel": (0.0, 0.2),
    "CpMeanDelta": (-0.2, 0.2),
}
PRESETS = {"magUrel": "Viridis (matplotlib)", "Cp": "Cool to Warm",
           "magUMeanDeltaRel": "Inferno (matplotlib)",
           "CpMeanDelta": "Blue Orange (divergent)"}
BG = [1.0, 1.0, 1.0]
PANEL_PX = (1300, 620)     # one view; figures are 2 x 2 or 1 x 2 views

# Camera: direction from the body centre (x downstream, z up), distance in
# body lengths; slices: y = y_centre, z = z_min + WHEEL_Z * height
VIEW = {
    "iso": dict(direction=(-1.0, -1.25, 0.8), up=(0, 0, 1), dist=1.65),
    "side": dict(direction=(0.0, -1.0, 0.0), up=(0, 0, 1), dist=1.6),
    "top": dict(direction=(0.0, 0.0, 1.0), up=(0, 1, 0), dist=1.25),
}
WHEEL_Z = 0.2
SLICE_PAD = (0.6, 1.6, 0.9, 0.9)   # upstream, downstream, lateral, vertical

notes: list[str] = []
T_START = time.time()


# --------------------------------------------------------------------------- #
# case inspection (pure python, no ParaView)
# --------------------------------------------------------------------------- #

def is_busy(case: Path, minutes: float) -> bool:
    now = time.time()
    for f in case.glob("log.*"):
        if now - f.stat().st_mtime < minutes * 60:
            return True
    return False


def times_with(case: Path, fields: tuple[str, ...]) -> list[float]:
    """Times (> 0) of processor0 that contain all fields."""
    p0 = case / "processor0"
    out = []
    if not p0.is_dir():
        return out
    for d in p0.iterdir():
        try:
            t = float(d.name)
        except ValueError:
            continue
        if t > 0 and all((d / f).exists() or (d / (f + ".gz")).exists()
                         for f in fields):
            out.append(t)
    return sorted(out)


def body_patches(case: Path) -> list[str]:
    """Wall patches of the body: type wall, not the tunnel walls."""
    b = case / "processor0" / "constant" / "polyMesh" / "boundary"
    if not b.exists():
        b = case / "constant" / "polyMesh" / "boundary"
    txt = b.read_text(errors="replace")
    out = []
    for m in re.finditer(r"\n\s*([^\s{}()/;]+)\s*\n\s*\{([^}]*)\}", txt):
        name, body = m.group(1), m.group(2)
        if re.search(r"\btype\s+wall\s*;", body) and not re.search(
                r"(lowerWall|upperWall|ground|floor|tunnel|frontAndBack)",
                name, re.I):
            out.append(name)
    return out


def u_inf(case: Path) -> float:
    txt = (case / "system" / "controlDict").read_text(errors="replace")
    m = re.search(r"magUInf\s+([0-9.eE+-]+)\s*;", txt)
    return float(m.group(1)) if m else 1.0


# --------------------------------------------------------------------------- #
# ParaView helpers
# --------------------------------------------------------------------------- #

def reader(case: Path, regions: list[str], fields: list[str], t: float):
    r = OpenFOAMReader(FileName=str(case / "system" / "controlDict"))
    r.CaseType = "Decomposed Case"
    r.MeshRegions = regions
    r.CellArrays = fields
    r.Createcelltopointfiltereddata = 1
    r.Decomposepolyhedra = 1   # snappy polyhedra: slices need them split
    r.UpdatePipelineInformation()
    r.UpdatePipeline(t)
    return r


def with_norm(src, uinf: float, have: set[str]):
    """Calculators for |U|/U_inf and C_p (point data)."""
    out = src
    if "U" in have:
        out = Calculator(Input=out, ResultArrayName="magUrel",
                         Function=f"mag(U)/{uinf}")
    if "p" in have:
        out = Calculator(Input=out, ResultArrayName="Cp",
                         Function=f"p/{0.5 * uinf * uinf}")
    if "UMean" in have:
        out = Calculator(Input=out, ResultArrayName="magUMeanrel",
                         Function=f"mag(UMean)/{uinf}")
    return out


def lut(array: str, rng_key: str | None = None):
    lt = GetColorTransferFunction(array)
    lt.ApplyPreset(PRESETS.get(rng_key or array, "Viridis (matplotlib)"), True)
    lo, hi = RANGES[rng_key or array]
    lt.RescaleTransferFunction(lo, hi)
    lt.AutomaticRescaleRangeMode = "Never"
    GetOpacityTransferFunction(array).RescaleTransferFunction(lo, hi)
    return lt


def new_view():
    v = CreateView("RenderView")
    v.ViewSize = list(PANEL_PX)
    v.Background = BG
    try:
        v.UseColorPaletteForBackground = 0
    except AttributeError:
        pass
    v.OrientationAxesVisibility = 1
    v.OrientationAxesLabelColor = [0.1, 0.1, 0.1]
    return v


def label(view, text: str, size: int = 34):
    t = Text(Text=text)
    d = Show(t, view)
    d.FontSize = size
    d.Color = [0.1, 0.1, 0.1]
    d.WindowLocation = "Upper Center"
    return t


def colour(disp, view, array: str, rng_key: str, title: str, bar=True):
    ColorBy(disp, ("POINTS", array))
    lt = lut(array, rng_key)
    disp.LookupTable = lt
    disp.SetScalarBarVisibility(view, bar)
    if bar:
        sb = GetScalarBar(lt, view)
        sb.Title = title
        sb.ComponentTitle = ""
        sb.TitleColor = [0.1, 0.1, 0.1]
        sb.LabelColor = [0.1, 0.1, 0.1]
        sb.TitleFontSize = 30
        sb.LabelFontSize = 26
        sb.Orientation = "Vertical"
        sb.WindowLocation = "Any Location"
        sb.Position = [0.9, 0.12]
        sb.ScalarBarLength = 0.7
        sb.RangeLabelFormat = "%-#.2g"


def camera(view, bb, kind: str):
    c = [(bb[0] + bb[1]) / 2, (bb[2] + bb[3]) / 2, (bb[4] + bb[5]) / 2]
    L = max(bb[1] - bb[0], bb[3] - bb[2], bb[5] - bb[4])
    v = VIEW[kind]
    d = v["direction"]
    n = sum(x * x for x in d) ** 0.5
    view.CameraFocalPoint = c
    view.CameraPosition = [c[i] + d[i] / n * v["dist"] * L for i in range(3)]
    view.CameraViewUp = list(v["up"])
    view.CameraParallelProjection = 0
    view.CameraViewAngle = 30


def plane_camera(view, bb, normal: str, pad):
    """Parallel projection onto a slice plane, window = body box + pad."""
    L = bb[1] - bb[0]
    x0, x1 = bb[0] - pad[0] * L, bb[1] + pad[1] * L
    c = [(bb[0] + bb[1]) / 2, (bb[2] + bb[3]) / 2, (bb[4] + bb[5]) / 2]
    cx = 0.5 * (x0 + x1)
    view.CameraParallelProjection = 1
    if normal == "y":
        z0, z1 = bb[4], bb[5] + pad[3] * (bb[5] - bb[4])
        view.CameraFocalPoint = [cx, c[1], 0.5 * (z0 + z1)]
        view.CameraPosition = [cx, c[1] - 10 * L, 0.5 * (z0 + z1)]
        view.CameraViewUp = [0, 0, 1]
        h = z1 - z0
    else:
        view.CameraFocalPoint = [cx, c[1], c[2]]
        view.CameraPosition = [cx, c[1], c[2] + 10 * L]
        view.CameraViewUp = [0, 1, 0]
        h = (bb[3] - bb[2]) * (1 + 2 * pad[2])
    asp = PANEL_PX[0] / PANEL_PX[1]
    scale = 0.5 * max(h, (x1 - x0) / asp) * 1.02
    view.CameraParallelScale = scale
    if normal == "y":
        # put the ground (z_min of the body box) at the lower edge
        zc = bb[4] + 0.99 * scale
        view.CameraFocalPoint = [cx, c[1], zc]
        view.CameraPosition = [cx, c[1] - 10 * L, zc]


def save_grid(views: list[list], path: Path) -> None:
    """Save a rows x cols grid of views as one PNG."""
    rows, cols = len(views), len(views[0])
    lay = CreateLayout(path.stem)
    # build a balanced split tree: rows first, then columns
    cells = {}

    def split_rows(cell, r0, r1):
        if r1 - r0 == 1:
            split_cols(cell, r0, 0, cols)
            return
        m = r0 + (r1 - r0) // 2
        a = lay.SplitVertical(cell, (m - r0) / (r1 - r0))
        b = a + 1   # children of layout cell c are 2c+1 and 2c+2
        split_rows(a, r0, m)
        split_rows(b, m, r1)

    def split_cols(cell, r, c0, c1):
        if c1 - c0 == 1:
            cells[(r, c0)] = cell
            return
        m = c0 + (c1 - c0) // 2
        a = lay.SplitHorizontal(cell, (m - c0) / (c1 - c0))
        b = a + 1
        split_cols(a, r, c0, m)
        split_cols(b, r, m, c1)

    split_rows(0, 0, rows)
    for (r, c), cell in cells.items():
        AssignViewToLayout(view=views[r][c], layout=lay, hint=cell)
    lay.SetSize(PANEL_PX[0] * cols, PANEL_PX[1] * rows)
    SaveScreenshot(str(path), lay, ImageResolution=[PANEL_PX[0] * cols,
                                                    PANEL_PX[1] * rows],
                   TransparentBackground=0)
    RemoveLayout(lay)
    print("wrote", path.name, flush=True)


def pending_view(text: str):
    v = new_view()
    v.OrientationAxesVisibility = 0
    v.Background = [0.93, 0.93, 0.93]
    label(v, text, 34)
    t = Text(Text="results pending")
    d = Show(t, v)
    d.FontSize = 40
    d.Color = [0.4, 0.4, 0.4]
    d.WindowLocation = "Any Location"
    d.Position = [0.38, 0.45]
    return v


# --------------------------------------------------------------------------- #
# figures
# --------------------------------------------------------------------------- #

class Src:
    """One solver's case at one time: body surface and internal mesh."""

    def __init__(self, case: Path, t: float, fields: list[str], body: list[str],
                 uinf: float):
        self.case, self.t, self.uinf = case, t, uinf
        have = set(fields)
        self.vol = with_norm(reader(case, ["internalMesh"], fields, t), uinf, have)
        self.body = with_norm(reader(case, ["patch/" + b for b in body],
                                     fields, t), uinf, have)
        self.body.UpdatePipeline(t)
        info = self.body.GetDataInformation()
        self.bb = list(info.GetBounds())

    def slices(self, bb):
        c = [(bb[0] + bb[1]) / 2, (bb[2] + bb[3]) / 2]
        sy = Slice(Input=self.vol)
        sy.SliceType = "Plane"
        sy.SliceType.Origin = [c[0], c[1], 0]
        sy.SliceType.Normal = [0, 1, 0]
        sz = Slice(Input=self.vol)
        sz.SliceType = "Plane"
        sz.SliceType.Origin = [c[0], c[1], bb[4] + WHEEL_Z * (bb[5] - bb[4])]
        sz.SliceType.Normal = [0, 0, 1]
        try:
            sy.Triangulatetheslice = 0
            sz.Triangulatetheslice = 0
        except AttributeError:
            pass
        return {"y": sy, "z": sz}


def fig_geometry(name, src: Src, out: Path, title: str):
    v = new_view()
    d = Show(src.body, v)
    d.ColorArrayName = ["POINTS", ""]
    d.SetRepresentationType("Surface With Edges")
    d.AmbientColor = d.DiffuseColor = [0.78, 0.78, 0.80]
    d.EdgeColor = [0.15, 0.15, 0.25]
    camera(v, src.bb, "iso")
    label(v, f"{title}: body surface mesh", 40)
    v.ViewSize = [PANEL_PX[0] * 2, int(PANEL_PX[1] * 1.6)]
    SaveScreenshot(str(out / f"render_{name}_geometry.png"), v,
                   ImageResolution=[PANEL_PX[0] * 2, int(PANEL_PX[1] * 1.6)])
    print("wrote", f"render_{name}_geometry.png", flush=True)
    Delete(v)


def slice_panel(src: Src | None, bb, plane: str, array: str, rng: str,
                bar_title: str, head: str):
    if src is None:
        return pending_view(head)
    v = new_view()
    sl = src.slices(bb)[plane]
    d = Show(sl, v)
    colour(d, v, array, rng, bar_title)
    b = Show(src.body, v)
    b.ColorArrayName = ["POINTS", ""]
    b.AmbientColor = b.DiffuseColor = [0.55, 0.55, 0.58]
    plane_camera(v, bb, plane, SLICE_PAD)
    label(v, head)
    return v


def fig_slices(name, sf: Src, cf: Src | None, out: Path, quantity: str):
    arr, rng, bt = {"U": ("magUrel", "magUrel", "mag(U) / U_inf"),
                    "p": ("Cp", "Cp", "C_p")}[quantity]
    bb = sf.bb
    rows = []
    for plane, pname in (("y", "mid-plane y = y_c"),
                         ("z", f"plane z = z_min + {WHEEL_Z:g} h")):
        rows.append([
            slice_panel(sf, bb, plane, arr, rng, bt, f"simpleFoam, {pname}"),
            slice_panel(cf, bb, plane, arr, rng, bt, f"coupledFoam, {pname}"),
        ])
    save_grid(rows, out / f"render_{name}_slice_{quantity}.png")


def surface_panel(src: Src | None, bb, array, rng, bt, head):
    if src is None:
        return pending_view(head)
    v = new_view()
    d = Show(src.body, v)
    colour(d, v, array, rng, bt)
    camera(v, bb, "iso")
    label(v, head)
    return v


def fig_surface(name, sf: Src, cf: Src | None, out: Path):
    bb = sf.bb
    row = [surface_panel(sf, bb, "Cp", "Cp", "C_p", "simpleFoam: surface C_p"),
           surface_panel(cf, bb, "Cp", "Cp", "C_p", "coupledFoam: surface C_p")]
    save_grid([row], out / f"render_{name}_surface_Cp.png")


def fig_mean(name, sfm: Src, cfm: Src, out: Path):
    bb = sfm.bb
    row = [slice_panel(sfm, bb, "y", "magUMeanrel", "magUrel", "mag(UMean) / U_inf",
                       "simpleFoam, window mean"),
           slice_panel(cfm, bb, "y", "magUMeanrel", "magUrel", "mag(UMean) / U_inf",
                       "coupledFoam, window mean")]
    save_grid([row], out / f"render_{name}_mean_U.png")


def fig_delta(name, dl: Src, out: Path):
    bb = dl.bb
    v1 = slice_panel(dl, bb, "y", "magUMeanDeltaRel", "magUMeanDeltaRel",
                     "mag(dUMean) / U_inf", "mean-velocity delta, mid-plane")
    v2 = slice_panel(dl, bb, "z", "magUMeanDeltaRel", "magUMeanDeltaRel",
                     "mag(dUMean) / U_inf", "mean-velocity delta, wheel height")
    v3 = surface_panel(dl, bb, "CpMeanDelta", "CpMeanDelta", "dC_p (mean)",
                       "mean surface-pressure delta")
    v4 = new_view()
    d4 = Show(dl.body, v4)
    colour(d4, v4, "CpMeanDelta", "CpMeanDelta", "dC_p (mean)")
    camera(v4, bb, "top")
    label(v4, "mean surface-pressure delta, top view")
    save_grid([[v1, v2], [v3, v4]], out / f"render_{name}_delta.png")


# --------------------------------------------------------------------------- #
# driver
# --------------------------------------------------------------------------- #

def render_case(name: str, cfg: dict, out: Path, busy_minutes: float,
                reference_only: bool = False) -> None:
    sfd, cfd = RUN / cfg["sf"], RUN / cfg["cf"]
    if reference_only:
        # never touch the coupledFoam directory (e.g. while it runs): its
        # panels are drawn as pending
        cfd = RUN / f"__not_read__{cfg['cf']}"
    if not sfd.is_dir():
        notes.append(f"{name}: run/{cfg['sf']} missing, skipped")
        return
    if is_busy(sfd, busy_minutes) or (cfd.is_dir() and is_busy(cfd, busy_minutes)):
        # a live run (or its reference continuation) may start writing the
        # other directory at any moment: leave the whole case alone
        notes.append(f"{name}: run/{cfg['sf']} or run/{cfg['cf']} is being "
                     "written, case skipped (previous renders kept)")
        return
    t_sf = times_with(sfd, ("U", "p"))
    if not t_sf:
        notes.append(f"{name}: no simpleFoam result time, skipped")
        return
    body = body_patches(sfd)
    uinf = u_inf(sfd)
    sf = Src(sfd, t_sf[-1], ["U", "p"], body, uinf)
    print(f"{name}: simpleFoam t={t_sf[-1]:g}, {len(body)} body patches, "
          f"U_inf {uinf:g}, body bounds {[round(x, 3) for x in sf.bb]}",
          flush=True)
    cf = None
    cf_ok = cfd.is_dir()
    t_cf = times_with(cfd, ("U", "p")) if cf_ok else []
    if t_cf:
        cf = Src(cfd, t_cf[-1], ["U", "p"], body, uinf)
        print(f"{name}: coupledFoam t={t_cf[-1]:g}", flush=True)
    elif cf_ok:
        notes.append(f"{name}: no coupledFoam result yet, panels pending")

    fig_geometry(name, sf, out, cfg["title"])
    fig_slices(name, sf, cf, out, "U")
    fig_slices(name, sf, cf, out, "p")
    fig_surface(name, sf, cf, out)

    tm_sf = times_with(sfd, ("UMean", "pMean"))
    tm_cf = times_with(cfd, ("UMean", "pMean")) if cf_ok else []
    if tm_sf and tm_cf:
        f = ["UMean", "pMean"]
        fig_mean(name, Src(sfd, tm_sf[-1], f, body, uinf),
                 Src(cfd, tm_cf[-1], f, body, uinf), out)
    else:
        notes.append(f"{name}: mean fields UMean/pMean not in both cases, "
                     "mean figure skipped")
    delta = ("magUMeanDeltaRel", "CpMeanDelta")
    td = times_with(cfd, delta) if cf_ok else []
    if td:
        fig_delta(name, Src(cfd, td[-1], list(delta), body, uinf), out)
    else:
        notes.append(f"{name}: delta fields (coupledFieldCompare) missing, "
                     "delta figure skipped")


def main(argv: list[str]) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cases", nargs="*", default=list(CASES))
    ap.add_argument("--out", default=str(REPO / "report" / "paper" / "figures"))
    ap.add_argument("--busy-minutes", type=float, default=10.0)
    ap.add_argument("--reference-only", action="store_true",
                    help="render only the simpleFoam reference; the "
                    "coupledFoam directory is not read (pending panels)")
    a = ap.parse_args(argv)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    for name in a.cases:
        if name not in CASES:
            notes.append(f"unknown case {name}")
            continue
        try:
            render_case(name, CASES[name], out, a.busy_minutes,
                        a.reference_only)
        except Exception as e:  # noqa: BLE001 - one broken case must not stop the others
            notes.append(f"{name}: render failed: {type(e).__name__}: {e}")
    (out / "render_fields_notes.json").write_text(json.dumps(notes, indent=1))
    # pvbatch buffers python output in its own stream, which os._exit below
    # would drop: write the summary straight to the file descriptor
    written = sorted(p.name for p in out.glob("render_*.png")
                     if p.stat().st_mtime >= T_START)
    lines = [f"wrote {w}" for w in written] + [f"note: {n}" for n in notes]
    os.write(1, ("\n".join(lines) + "\n").encode())
    # skip the interpreter teardown: pvbatch under WSLg reports a harmless
    # GLXBadContext when the render windows are destroyed at exit
    os._exit(0)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

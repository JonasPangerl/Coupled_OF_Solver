#!/usr/bin/env python3
"""Ahmed body, 25 degree rear slant, half model (spec 13, T5).

Writes ahmed25.stl (ASCII, one solid named "ahmed") and ahmed25_stats.json
next to this script.

Geometry (Ahmed et al. 1984), metres, x downstream, z up, nose at x = 0:
    length L 1.044, width W 0.389, height H 0.288
    front edges rounded with R = 0.100: the two horizontal front edges (top,
    bottom) and the vertical front edge; the longitudinal edges are sharp
    rear slant: 0.222 long at 25 degrees to the roof, ending at the rear face
    ground clearance 0.050 (underside at z = 0.050), no stilts
    half model: symmetry plane y = 0, body from y = 0 to y = W/2

Construction: the body is the intersection of two extrusions, the side
profile (x-z: rounded top and bottom front edges, slant) extruded along y
and the plan profile (x-y: rounded vertical front edge) extruded along z.
Every cross-section x = const of that intersection is a rectangle, so the
surface is a sweep of rectangles along x, closed by flat nose and rear caps.
The front roundings are resolved with N_ARC stations (angle steps of
90/N_ARC degrees).

Symmetry plane: by default the surface is extended through the symmetry
plane to y = -Y_OVERLAP (5 mm) so that snappyHexMesh cuts it cleanly at the
domain boundary y = 0 instead of meeting a coincident face; the part at
y < 0 lies outside the domain and does not exist in the mesh. The
cross-section is constant in y there, so the meshed body is exactly the half
body y in [0, W/2]. --y-min 0 writes the exact closed half body instead.
All statistics refer to the half body y >= 0.

The output is deterministic (fixed station list, fixed float formatting).
The surface is checked before writing: closed and manifold (every edge
shared by exactly two triangles), consistently oriented (every directed
edge once), outward normals (positive volume), no degenerate triangles.

Usage:  python3 make_ahmed.py [--y-min Y] [--out DIR]
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from stl import mesh as stlmesh

# Dimensions (m)
L = 1.044
W = 0.389
H = 0.288
R = 0.100
SLANT_LEN = 0.222
SLANT_DEG = 25.0
CLEARANCE = 0.050

Y_OVERLAP = 0.005       # extension through the symmetry plane (m)
N_ARC = 24              # stations over each front rounding (3.75 deg steps)
DX_MAX = 0.02           # maximum station spacing elsewhere (m)

Z_BOT = CLEARANCE
Z_TOP = CLEARANCE + H
SLANT_DX = SLANT_LEN * math.cos(math.radians(SLANT_DEG))
SLANT_DZ = SLANT_LEN * math.sin(math.radians(SLANT_DEG))
X_SLANT = L - SLANT_DX          # start of the slant on the roof


def round_offset(x: float) -> float:
    """Inset of the rounded front edges at station x (0 for x >= R)."""
    if x >= R:
        return 0.0
    d = R - x
    return R - math.sqrt(max(R * R - d * d, 0.0))


def section(x: float) -> tuple[float, float, float]:
    """(y_max, z_min, z_max) of the rectangular cross-section at x."""
    e = round_offset(x)
    y_max = W / 2 - e
    z_min = Z_BOT + e
    z_max = Z_TOP - e
    if x > X_SLANT:
        z_max = Z_TOP - (x - X_SLANT) * math.tan(math.radians(SLANT_DEG))
    return y_max, z_min, z_max


def stations() -> list[float]:
    """Deterministic x stations: angle-uniform over the rounding, uniform
    (<= DX_MAX) on the flat part and on the slant, exact at the slant kink."""
    xs = [R * (1.0 - math.cos(0.5 * math.pi * i / N_ARC))
          for i in range(N_ARC + 1)]
    for a, b in ((R, X_SLANT), (X_SLANT, L)):
        n = max(1, math.ceil((b - a) / DX_MAX))
        xs += [a + (b - a) * i / n for i in range(1, n + 1)]
    return xs


def build(y_min: float) -> tuple[np.ndarray, np.ndarray]:
    """Vertices (n, 3) and triangles (m, 3), outward oriented."""
    xs = stations()
    verts = []
    for x in xs:
        y_max, z_min, z_max = section(x)
        # Loop, counter-clockwise seen from +x (looking upstream):
        # (y_min, z_min) -> (y_max, z_min) -> (y_max, z_max) -> (y_min, z_max)
        verts += [(x, y_min, z_min), (x, y_max, z_min),
                  (x, y_max, z_max), (x, y_min, z_max)]
    v = np.array(verts, dtype=float)

    tris = []
    n = len(xs)
    for i in range(n - 1):
        a, b = 4 * i, 4 * (i + 1)
        for k in range(4):
            k1 = (k + 1) % 4
            # Quad between loop points k, k1 of stations i, i+1; the loop
            # runs counter-clockwise seen from +x, so this order gives
            # outward normals
            tris.append((a + k, b + k1, b + k))
            tris.append((a + k, a + k1, b + k1))
    # Nose cap (normal -x) and rear cap (normal +x)
    tris.append((0, 2, 1))
    tris.append((0, 3, 2))
    e = 4 * (n - 1)
    tris.append((e, e + 1, e + 2))
    tris.append((e, e + 2, e + 3))
    return v, np.array(tris, dtype=np.int64)


def check_closed(v: np.ndarray, t: np.ndarray) -> dict:
    """Watertightness and orientation checks; raises on failure."""
    edges: dict[tuple[int, int], int] = {}
    directed: dict[tuple[int, int], int] = {}
    for tri in t:
        for j in range(3):
            p, q = int(tri[j]), int(tri[(j + 1) % 3])
            edges[(min(p, q), max(p, q))] = edges.get((min(p, q), max(p, q)), 0) + 1
            directed[(p, q)] = directed.get((p, q), 0) + 1
    bad = {e: c for e, c in edges.items() if c != 2}
    if bad:
        raise RuntimeError(f"not watertight: {len(bad)} edges not shared by "
                           f"exactly two triangles")
    dup = {e: c for e, c in directed.items() if c != 1}
    if dup:
        raise RuntimeError(f"inconsistent orientation: {len(dup)} directed "
                           f"edges used more than once")
    p0, p1, p2 = v[t[:, 0]], v[t[:, 1]], v[t[:, 2]]
    area2 = np.linalg.norm(np.cross(p1 - p0, p2 - p0), axis=1)
    if area2.min() <= 1e-14:
        raise RuntimeError("degenerate triangle")
    vol = float(np.einsum("ij,ij->i", p0, np.cross(p1, p2)).sum() / 6.0)
    if vol <= 0:
        raise RuntimeError("inward oriented surface (negative volume)")
    # Euler characteristic of a closed genus-0 surface: V - E + F = 2
    euler = len(v) - len(edges) + len(t)
    if euler != 2:
        raise RuntimeError(f"Euler characteristic {euler} != 2")
    return {"edges": len(edges), "euler": euler, "volume": vol,
            "area": float(area2.sum() / 2.0)}


def frontal_area(v: np.ndarray, t: np.ndarray) -> float:
    """Projected area onto the y-z plane: half the sum of |n_x| dA over the
    closed surface (the body is x-monotone in its silhouette)."""
    p0, p1, p2 = v[t[:, 0]], v[t[:, 1]], v[t[:, 2]]
    nx = np.cross(p1 - p0, p2 - p0)[:, 0]
    return float(np.abs(nx).sum() / 4.0)


def write_ascii_stl(path: Path, v: np.ndarray, t: np.ndarray,
                    name: str) -> None:
    """ASCII STL with fixed formatting (deterministic, numpy-stl for the
    normals)."""
    m = stlmesh.Mesh(np.zeros(len(t), dtype=stlmesh.Mesh.dtype))
    m.vectors[:] = v[t]
    m.update_normals()
    normals = m.normals / np.linalg.norm(m.normals, axis=1)[:, None]
    with open(path, "w", newline="\n") as fh:
        fh.write(f"solid {name}\n")
        for n, tri in zip(normals, m.vectors):
            fh.write(f"  facet normal {n[0]:.9e} {n[1]:.9e} {n[2]:.9e}\n")
            fh.write("    outer loop\n")
            for p in tri:
                fh.write(f"      vertex {p[0]:.9e} {p[1]:.9e} {p[2]:.9e}\n")
            fh.write("    endloop\n")
            fh.write("  endfacet\n")
        fh.write(f"endsolid {name}\n")


def main() -> int:
    here = Path(__file__).resolve().parent
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--y-min", type=float, default=-Y_OVERLAP,
                    help="inner y of the surface (default -%g: extended "
                    "through the symmetry plane; 0: exact half body)"
                    % Y_OVERLAP)
    ap.add_argument("--out", type=Path, default=here)
    a = ap.parse_args()

    v, t = build(a.y_min)
    chk = check_closed(v, t)

    # Statistics of the half body y >= 0 (the meshed body)
    vh, th = build(0.0)
    half = check_closed(vh, th)
    xs = stations()
    heights = [section(x)[2] - section(x)[1] for x in xs]
    # Area of the cut face on y = 0 (polygon; trapezoidal rule is exact)
    cap_y0 = sum(0.5 * (heights[i] + heights[i + 1]) * (xs[i + 1] - xs[i])
                 for i in range(len(xs) - 1))
    a_front = frontal_area(vh, th)
    stats = {
        "file": "ahmed25.stl",
        "solid": "ahmed",
        "units": "m",
        "length": L, "width": W, "height": H,
        "frontRadius": R, "slantLength": SLANT_LEN, "slantAngleDeg": SLANT_DEG,
        "groundClearance": CLEARANCE,
        "slantStartX": X_SLANT,
        "slantDrop": SLANT_DZ,
        "halfModel": True,
        "surfaceYmin": a.y_min,
        "bounds": {"min": v.min(axis=0).tolist(), "max": v.max(axis=0).tolist()},
        "nTriangles": int(len(t)),
        "nStations": len(xs),
        "watertight": True,
        "edges": chk["edges"],
        "eulerCharacteristic": chk["euler"],
        # half body y >= 0 (forceCoeffs Aref = frontalAreaHalf)
        "frontalAreaHalf": a_front,
        "frontalAreaFull": 2.0 * a_front,
        "volumeHalf": half["volume"],
        "volumeFull": 2.0 * half["volume"],
        "wettedAreaHalf": half["area"] - cap_y0,
        "centre": [L / 2, 0.0, Z_BOT + H / 2],
        "surfaceVolumeWritten": chk["volume"],
    }
    a.out.mkdir(parents=True, exist_ok=True)
    write_ascii_stl(a.out / "ahmed25.stl", v, t, "ahmed")
    (a.out / "ahmed25_stats.json").write_text(json.dumps(stats, indent=2) + "\n")
    print(f"ahmed25.stl: {len(t)} triangles, watertight, "
          f"half-body frontal area {stats['frontalAreaHalf']:.6f} m2, "
          f"volume {stats['volumeHalf']:.6f} m3")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

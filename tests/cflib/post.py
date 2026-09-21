"""postProcessing readers: sets (.xy), forceCoeffs, surfaceFieldValue."""

from __future__ import annotations

from pathlib import Path

import numpy as np


def _latest_subdir(d: Path) -> Path | None:
    subs = []
    for s in d.iterdir():
        try:
            subs.append((float(s.name), s))
        except ValueError:
            pass
    return max(subs)[1] if subs else None


def read_xy(path: Path) -> np.ndarray:
    return np.loadtxt(path, comments="#")


def sets_file(case: Path, fo: str, name: str, time: str | None = None) -> Path:
    base = case / "postProcessing" / fo
    d = base / time if time else _latest_subdir(base)
    if d is None:
        raise FileNotFoundError(f"no output of {fo} in {case}")
    matches = sorted(d.glob(f"{name}_*.xy"))
    if not matches:
        raise FileNotFoundError(f"no {name}_*.xy in {d}")
    return matches[0]


def line_profile(coord: np.ndarray, values: np.ndarray,
                 rel_tol: float = 1e-9) -> tuple[np.ndarray, np.ndarray]:
    """Sorted, de-duplicated line samples.

    A parallel `sets` run writes the sample points of every rank in rank
    order: points on a processor boundary appear twice (once per rank) and
    points that no rank claims are missing. Sort by the line coordinate and
    merge points whose coordinates agree to rel_tol of the line length
    (their values are averaged)."""
    order = np.argsort(coord, kind="stable")
    c = np.asarray(coord, dtype=float)[order]
    v = np.asarray(values, dtype=float)[order]
    span = float(c[-1] - c[0]) if c.size > 1 else 1.0
    tol = rel_tol * max(abs(span), 1.0e-300)
    out_c, out_v = [], []
    i = 0
    while i < c.size:
        j = i + 1
        while j < c.size and c[j] - c[i] <= tol:
            j += 1
        out_c.append(float(c[i:j].mean()))
        out_v.append(float(v[i:j].mean()))
        i = j
    return np.array(out_c), np.array(out_v)


def match_profiles(profile: tuple[np.ndarray, np.ndarray],
                   reference: tuple[np.ndarray, np.ndarray],
                   rel_tol: float = 1e-9
                   ) -> tuple[np.ndarray, np.ndarray, int]:
    """Values of two line profiles at the sample points present in both.

    Both profiles are cleaned with line_profile() (sorted, duplicates of a
    parallel run merged). A parallel `sets` run of a line that lies on
    cell faces (T0 centrelines) also drops whole runs of points: measured
    on T0 np4, up to 3 consecutive points of 129. Interpolating across such
    a gap is not an option - linear interpolation over 4 sample spacings
    of the lid boundary layer is off by 2.6e-3 while the sampled values
    agree to 1e-9 - so missing points are left out and counted instead.
    Returns (profile values, reference values, number of reference points
    missing in the profile)."""
    pc, pv = line_profile(*profile, rel_tol=rel_tol)
    rc, rv = line_profile(*reference, rel_tol=rel_tol)
    span = float(rc[-1] - rc[0]) if rc.size > 1 else 1.0
    tol = rel_tol * max(abs(span), 1.0e-300)
    idx = np.clip(np.searchsorted(pc, rc), 0, max(pc.size - 1, 0))
    best = idx.copy()
    left = np.maximum(idx - 1, 0)
    use_left = np.abs(pc[left] - rc) < np.abs(pc[idx] - rc)
    best[use_left] = left[use_left]
    present = np.abs(pc[best] - rc) <= tol
    return pv[best[present]], rv[present], int(np.count_nonzero(~present))


def l2rel(a: np.ndarray, b: np.ndarray) -> float:
    """||a - b||_2 / ||b||_2"""
    nb = np.linalg.norm(b)
    return float(np.linalg.norm(a - b) / nb) if nb > 0 else float(np.linalg.norm(a - b))


def read_dat(path: Path) -> dict[str, np.ndarray]:
    """Columnar .dat with a '# Time ...' header line (forceCoeffs,
    surfaceFieldValue)."""
    header = None
    rows = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                toks = line[1:].split()
                if toks and toks[0] == "Time":
                    header = toks
                continue
            if line.strip():
                rows.append([float(x) for x in line.split()])
    data = np.array(rows)
    if header is None:
        header = [f"c{i}" for i in range(data.shape[1])]
    return {h: data[:, i] for i, h in enumerate(header[: data.shape[1]])}


def function_object_file(case: Path, fo: str, fname: str) -> Path:
    """postProcessing/<fo>/<startTime>/<fname>, concatenating restarts is
    left to the caller; the latest start directory is returned."""
    base = case / "postProcessing" / fo
    d = _latest_subdir(base)
    if d is None:
        raise FileNotFoundError(f"no output of {fo} in {case}")
    f = d / fname
    if not f.exists():
        cands = sorted(d.glob("*.dat"))
        if not cands:
            raise FileNotFoundError(f"{f} not found")
        f = cands[0]
    return f


def all_starts(case: Path, fo: str, fname: str) -> list[Path]:
    base = case / "postProcessing" / fo
    out = []
    for s in sorted(base.iterdir(), key=lambda p: float(p.name) if p.name.replace('.', '', 1).isdigit() else -1):
        f = s / fname
        if f.exists():
            out.append(f)
    return out


def force_coeffs(case: Path, fo: str = "forceCoeffs") -> dict[str, np.ndarray]:
    """Cd/Cl history; restarts (several start directories) are merged."""
    files = all_starts(case, fo, "coefficient.dat")
    if not files:
        files = [function_object_file(case, fo, "coefficient.dat")]
    merged: dict[str, list] = {}
    for f in files:
        d = read_dat(f)
        for k, v in d.items():
            merged.setdefault(k, []).extend(v.tolist())
    out = {k: np.array(v) for k, v in merged.items()}
    # Keep the last value per time (restart overlap)
    if "Time" in out:
        _, idx = np.unique(out["Time"][::-1], return_index=True)
        keep = np.sort(len(out["Time"]) - 1 - idx)
        out = {k: v[keep] for k, v in out.items()}
    return out


def surface_value(case: Path, fo: str) -> dict[str, np.ndarray]:
    return read_dat(function_object_file(case, fo, "surfaceFieldValue.dat"))


def window_converged(h: np.ndarray, window: int, tol: float) -> int | None:
    """First index i (1-based count of samples) at which the last `window`
    samples satisfy max - min <= tol*|mean| (spec 12.3 ii)."""
    for i in range(window, len(h) + 1):
        w = h[i - window:i]
        if w.max() - w.min() <= tol * abs(w.mean()):
            return i
    return None

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

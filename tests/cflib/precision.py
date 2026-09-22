"""Single-precision (SP) harness: amendment D5 (geometry safeguards) and the
D11 record fields; DECISIONS.md D-065.

Selection: CF_PRECISION=sp (pytest also --precision sp). Precision is a
build choice, so pytest / bench/run_bench.py must then run in a shell with
the SP OpenFOAM environment sourced (WM_PRECISION_OPTION=SP, e.g.
`cfenv sp`) and the SP build of libcoupledFoam / coupledFoam first on
LD_LIBRARY_PATH / PATH (FOAM_USER_LIBBIN / FOAM_USER_APPBIN of the SP
private install). The DP steps run in a clean DP environment of their own
(CF_DP_BASHRC, default the system v2606 DP install).

Every SP case is prepared once, on its first cflib.case.allrun (not for
-restart):

 1. DP (D5.1): controlDict writeFormat ascii (the templates have
    writePrecision 12), `Allrun -mesh-only <args>`; a binary or compressed
    mesh (T3 polyMesh.orig) is converted with foamFormatConvert -constant;
    `checkMesh -allGeometry -allTopology -constant` -> log.checkMesh.dp.
 2. writeFormat restored (binary, D8.3).
 3. SP (D5.2): `transformPoints -translate '(-cx -cy -cz)'` with (cx cy cz)
    the centre of the "Overall domain bounding box" of the DP checkMesh;
    constant/meshShift records the vector (points_SP = points_DP + shift);
    every point-based setting (CofR of coupledForcesDict and forceCoeffs,
    sets start/end, probe locations, pRefPoint, origin/centre, ... see
    POINT_KEYS) in system/, constant/{MRFProperties,fvOptions} and 0.orig/
    is shifted by the same vector (text edit, full precision).
 4. SP (D5.3): `checkMesh -allGeometry -allTopology -constant`
    -> log.checkMesh.sp. The SP mesh must not fail a check that the DP mesh
    passes, and must have no negative-volume cell and no incorrectly
    oriented face; otherwise the case is SP-geometry-fail: SPGeometryFail
    is raised (a pytest skip; run_bench records and skips the run).
 5. The solver runs through the template's Allrun with -keep-mesh.

<case>/spHarness.json holds the shift, both checkMesh summaries, the diff,
the shifted entries and the timing of the preparation. annotate() copies
the D11 fields into a result record; the schema is bench/SCHEMA_precision.md.
Results with coordinates are mapped back to the DP frame with
post.unshift() (coordinate - shift) before they are compared.
"""

from __future__ import annotations

import json
import os
import re
import shlex
import time
from pathlib import Path

from . import env as cfenv

try:                                   # a skip inside pytest, else an error
    import pytest
    _SkipBase = pytest.skip.Exception
except ImportError:                    # pragma: no cover
    _SkipBase = Exception

PRECISIONS = ("dp", "sp")
SUFFIX = "_sp"
DP_BASHRC_DEFAULT = "/usr/lib/openfoam/openfoam2606/etc/bashrc"
MESH_SHIFT = "constant/meshShift"
HARNESS_JSON = "spHarness.json"
SP_GEOMETRY_FAIL = "SP-geometry-fail"
# D10: dp / Cl / Cd of an SP run within 0.3 % of the reference precision
# (DP here, D-062)
TOL_SP_VS_DP = 0.003

# Point-valued keywords (shifted with the mesh) and vector keywords that are
# directions (never shifted). A three-component entry with any other
# keyword outside 0.orig/ is listed in spHarness.json "unclassifiedVectors".
POINT_KEYS = ("CofR", "centreOfRotation", "origin", "point", "start", "end",
              "pRefPoint", "p1", "p2", "centre", "basePoint", "refPoint",
              "location", "position", "locationInMesh")
POINT_LIST_KEYS = ("probeLocations", "points")
DIRECTION_KEYS = ("dragDir", "liftDir", "pitchAxis", "axis", "normal",
                  "direction", "flowDir", "e1", "e2", "e3", "n",
                  "rotationAxis", "flowVelocity", "value", "inletValue",
                  "internalField", "Uinf", "freestreamValue", "g")
# Mesh-generation dictionaries: used in DP only (unshifted frame)
MESH_DICTS = ("blockMeshDict", "snappyHexMeshDict", "surfaceFeatureExtractDict",
              "surfaceFeaturesDict", "meshQualityDict", "decomposeParDict",
              "extrudeMeshDict", "foamyHexMeshDict", "cartesianMeshDict")

_NUM = r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?"
_VEC = rf"\(\s*({_NUM})\s+({_NUM})\s+({_NUM})\s*\)"
_POINT_RE = re.compile(
    rf"(?m)^([ \t]*({'|'.join(POINT_KEYS)})[ \t]+){_VEC}([ \t]*;)")
_VEC_NC = rf"\(\s*{_NUM}\s+{_NUM}\s+{_NUM}\s*\)"
# groups: 1 "key (", 2 key, 3 body of point triples, 4 ");"
_LIST_RE = re.compile(
    rf"(?s)(\b({'|'.join(POINT_LIST_KEYS)})\s*\()((?:\s*{_VEC_NC})*\s*)(\)\s*;)")
_ANYVEC_RE = re.compile(rf"(?m)^[ \t]*(\w+)[ \t]+(?:uniform[ \t]+)?{_VEC}[ \t]*;")
_TRIPLE_RE = re.compile(_VEC)


class SPGeometryFail(_SkipBase):
    """D5.3 gate failed: the SP benchmark/test of the case is skipped (not a
    solver failure)."""


# --------------------------------------------------------------------------- #
# selection
# --------------------------------------------------------------------------- #

def current() -> str:
    """'dp' (default) or 'sp' from CF_PRECISION."""
    p = os.environ.get("CF_PRECISION", "dp").strip().lower() or "dp"
    if p not in PRECISIONS:
        raise ValueError(f"CF_PRECISION={p!r}: use dp or sp")
    return p


def is_sp() -> bool:
    return current() == "sp"


def build_precision() -> str:
    """WM_PRECISION_OPTION of the sourced build (DP, SP, SPDP)."""
    return os.environ.get("WM_PRECISION_OPTION", "")


def label() -> str:
    """Value of the record field `precision`: 'dp', 'sp', or 'dp-shifted'
    for the harness self-test (the SP procedure run with a DP build,
    CF_SP_ALLOW_DP_BUILD=1: must reproduce the unshifted DP numbers)."""
    if not is_sp():
        return "dp"
    return "sp" if build_precision() == "SP" else "dp-shifted"


def check_environment() -> None:
    """CF_PRECISION must match the sourced build (WM_PRECISION_OPTION)."""
    wm = build_precision()
    if is_sp() and wm != "SP" and os.environ.get("CF_SP_ALLOW_DP_BUILD") != "1":
        raise RuntimeError(
            "CF_PRECISION=sp needs the SP OpenFOAM environment "
            f"(WM_PRECISION_OPTION=SP), found {wm!r} (harness self-test "
            "with a DP build: CF_SP_ALLOW_DP_BUILD=1)")
    if not is_sp() and wm == "SP":
        raise RuntimeError("an SP OpenFOAM environment is sourced: set "
                           "CF_PRECISION=sp (or source the DP build)")


def suffix() -> str:
    """Run/record name suffix of SP runs: _sp (CF_PRECISION_TAG overrides,
    e.g. _shiftdp for the DP self-test)."""
    return os.environ.get("CF_PRECISION_TAG", SUFFIX)


def tag(name: str) -> str:
    """Run / record name of the current precision: <name>_sp in SP."""
    s = suffix()
    return name + s if is_sp() and not name.endswith(s) else name


def dp_name(name: str) -> str:
    """DP counterpart of a (possibly tagged) name."""
    s = suffix()
    return name[:-len(s)] if name.endswith(s) else name


# --------------------------------------------------------------------------- #
# DP commands from an SP process
# --------------------------------------------------------------------------- #

_KEEP_ENV = ("HOME", "USER", "LOGNAME", "LANG", "LC_ALL", "TERM",
             "CF_MPI_BIND", "CF_MPI_CPUSET", "CF_NICE", "TMPDIR")


def dp_command(cmd: list[str], cwd: Path) -> list[str]:
    """cmd in a clean DP environment: env -i, source CF_DP_BASHRC (no
    positional arguments: etc/bashrc interprets them), check the precision,
    cd, run."""
    bashrc = os.environ.get("CF_DP_BASHRC", DP_BASHRC_DEFAULT)
    keep = [f"{k}={os.environ[k]}" for k in _KEEP_ENV if k in os.environ]
    script = ('set --; . "$CF_DP_BASHRC" > /dev/null 2>&1; '
              '[ "$WM_PRECISION_OPTION" = DP ] || '
              '{ echo "CF_DP_BASHRC is not a DP build" 1>&2; exit 96; }; '
              'cd "$CF_DP_CWD" || exit 98; eval "$CF_DP_CMD"')
    return (["env", "-i"] + keep
            + ["PATH=/usr/local/bin:/usr/bin:/bin",
               f"CF_DP_BASHRC={bashrc}", f"CF_DP_CWD={cwd}",
               f"CF_DP_CMD={shlex.join([str(c) for c in cmd])}",
               "bash", "--noprofile", "--norc", "-c", script])


def dp_run(cmd: list[str], cwd: Path, log: Path) -> int:
    return cfenv.run(dp_command(cmd, cwd), cwd=cwd, log=log)


# --------------------------------------------------------------------------- #
# checkMesh
# --------------------------------------------------------------------------- #

_BBOX = re.compile(rf"Overall domain bounding box {_VEC} {_VEC}")
_FAILED_N = re.compile(r"Failed (\d+) mesh checks")
_NEGVOL = re.compile(r"Number of negative volume cells: (\d+)")
_ORIENT = re.compile(r"(\d+) faces are incorrectly oriented")
_CELLS = re.compile(r"^\s*cells:\s+(\d+)", re.M)
_NUMTOK = re.compile(_NUM)


def normalise_check(line: str) -> str:
    """Failed-check message without its numbers (the kind of check)."""
    s = line.strip().lstrip("*").strip()
    s = _NUMTOK.sub("#", s)
    return re.sub(r"\s+", " ", s)


def parse_checkmesh(log: Path) -> dict:
    """Summary of a checkMesh log: bounding box, failed checks (lines
    starting with ***), negative-volume cells, incorrectly oriented faces."""
    text = log.read_text(errors="replace") if log.exists() else ""
    out: dict = {"log": log.name, "ok": "Mesh OK." in text}
    m = _BBOX.search(text)
    out["bbox"] = ([float(m.group(i)) for i in (1, 2, 3)],
                   [float(m.group(i)) for i in (4, 5, 6)]) if m else None
    failed = [ln.strip() for ln in text.splitlines()
              if ln.lstrip().startswith("***")]
    out["failed"] = failed
    out["failedKinds"] = sorted({normalise_check(f) for f in failed})
    m = _FAILED_N.search(text)
    out["nFailedChecks"] = int(m.group(1)) if m else 0
    out["negativeVolumeCells"] = sum(int(x) for x in _NEGVOL.findall(text))
    out["incorrectlyOrientedFaces"] = sum(int(x) for x in _ORIENT.findall(text))
    m = _CELLS.search(text)
    out["cells"] = int(m.group(1)) if m else None
    out["fatal"] = "FOAM FATAL" in text
    return out


def checkmesh_diff(dp: dict, sp: dict) -> list[dict]:
    """Failed-check kinds that differ between DP and SP:
    [{"check": <kind>, "only": "sp"|"dp"}]."""
    a, b = set(dp.get("failedKinds", [])), set(sp.get("failedKinds", []))
    return ([{"check": k, "only": "sp"} for k in sorted(b - a)]
            + [{"check": k, "only": "dp"} for k in sorted(a - b)])


def gate_overrides() -> list[str]:
    """CF_SP_GATE_OVERRIDE: comma-separated substrings of failed-check kinds
    that are reported but do not fail the D5.3 gate. Default empty (the
    gate as specified); an override is a user decision and is recorded
    (spGeometry "ok-overridden", gateOverridden)."""
    return [s.strip() for s in os.environ.get("CF_SP_GATE_OVERRIDE",
                                              "").split(",") if s.strip()]


def gate(dp: dict, sp: dict) -> tuple[bool, list[str]]:
    """D5.3: no SP-only failed check, 0 negative-volume cells and 0
    incorrectly oriented faces in SP. SP-only checks matching
    gate_overrides() are not counted (see gate_overridden)."""
    why = []
    if sp.get("fatal") or sp.get("bbox") is None:
        why.append("SP checkMesh did not complete")
    only_sp = [c for c in gate_overridden(dp, sp, invert=True)]
    if only_sp:
        why.append(f"SP-only failed checks: {only_sp}")
    if sp.get("negativeVolumeCells", 0):
        why.append(f"{sp['negativeVolumeCells']} negative-volume cells")
    if sp.get("incorrectlyOrientedFaces", 0):
        why.append(f"{sp['incorrectlyOrientedFaces']} incorrectly oriented faces")
    return not why, why


def gate_overridden(dp: dict, sp: dict, invert: bool = False) -> list[str]:
    """SP-only failed-check kinds that match CF_SP_GATE_OVERRIDE (invert:
    those that do not match)."""
    ov = gate_overrides()
    only_sp = [d["check"] for d in checkmesh_diff(dp, sp) if d["only"] == "sp"]
    hit = [c for c in only_sp if any(o in c for o in ov)]
    return [c for c in only_sp if c not in hit] if invert else hit


# --------------------------------------------------------------------------- #
# mesh shift
# --------------------------------------------------------------------------- #

def shift_from_bbox(bbox) -> list[float]:
    """-(bounding-box centre), rounded to 12 significant digits (the box is
    printed with 6 by checkMesh; this only removes the binary noise of
    (lo + hi)/2 so that the recorded and the applied vector read alike)."""
    lo, hi = bbox
    return [float(f"{-(a + b) / 2.0:.12g}") + 0.0 for a, b in zip(lo, hi)]


def _fmt(v: float) -> str:
    s = f"{v:.15g}"
    return "0" if s in ("-0", "0") else s


def fmt_vec(v) -> str:
    return "(" + " ".join(_fmt(float(x)) for x in v) + ")"


def write_mesh_shift(case: Path, shift, bbox_dp) -> None:
    """constant/meshShift: points_SP = points_DP + shift."""
    centre = [-s for s in shift]
    text = (
        "/*--------------------------------*- C++ -*----------------------------------*\\\n"
        "  Origin shift of the SP run (amendment D5.2, tests/cflib/precision.py):\n"
        "  points_SP = points_DP + shift; map a result coordinate back with\n"
        "  x_DP = x_SP - shift.\n"
        "\\*---------------------------------------------------------------------------*/\n"
        "FoamFile\n{\n    version     2.0;\n    format      ascii;\n"
        "    class       dictionary;\n    object      meshShift;\n}\n\n"
        f"shift           {fmt_vec(shift)};\n"
        f"boundingBoxCentre {fmt_vec(centre)};\n"
        f"boundingBoxMin  {fmt_vec(bbox_dp[0])};\n"
        f"boundingBoxMax  {fmt_vec(bbox_dp[1])};\n\n"
        "// ************************************************************************* //\n")
    (case / MESH_SHIFT).write_text(text)


def read_mesh_shift(case: Path) -> list[float] | None:
    f = Path(case) / MESH_SHIFT
    if not f.exists():
        return None
    m = re.search(rf"(?m)^shift\s+{_VEC}\s*;", f.read_text())
    return [float(m.group(i)) for i in (1, 2, 3)] if m else None


def _shift_triple(m: re.Match, shift, groups=(1, 2, 3)) -> list[float]:
    return [float(m.group(g)) + s for g, s in zip(groups, shift)]


def shift_text(text: str, shift) -> tuple[str, list[dict]]:
    """Shift every point-valued entry of a dictionary text by `shift`.
    Returns the new text and the list of changes."""
    changes: list[dict] = []

    def one(m: re.Match) -> str:
        new = _shift_triple(m, shift, (3, 4, 5))
        changes.append({"key": m.group(2),
                        "old": [float(m.group(i)) for i in (3, 4, 5)],
                        "new": new})
        return f"{m.group(1)}{fmt_vec(new)}{m.group(6)}"

    text = _POINT_RE.sub(one, text)

    def lst(m: re.Match) -> str:
        body = m.group(3)
        olds, news = [], []

        def tri(t: re.Match) -> str:
            new = _shift_triple(t, shift)
            olds.append([float(t.group(i)) for i in (1, 2, 3)])
            news.append(new)
            return fmt_vec(new)

        body = _TRIPLE_RE.sub(tri, body)
        changes.append({"key": m.group(2), "old": olds, "new": news})
        return f"{m.group(1)}{body}{m.group(4)}"

    text = _LIST_RE.sub(lst, text)
    return text, changes


def unclassified_vectors(text: str) -> list[str]:
    known = set(POINT_KEYS) | set(DIRECTION_KEYS) | set(POINT_LIST_KEYS)
    return sorted({m.group(1) for m in _ANYVEC_RE.finditer(text)
                   if m.group(1) not in known})


def shift_files(case: Path) -> list[Path]:
    """Dictionaries whose point settings belong to the solver run."""
    out = []
    sysdir = case / "system"
    for f in sorted(sysdir.iterdir()) if sysdir.is_dir() else []:
        if f.is_file() and f.name not in MESH_DICTS:
            out.append(f)
    for n in ("MRFProperties", "fvOptions", "dynamicMeshDict"):
        f = case / "constant" / n
        if f.is_file():
            out.append(f)
    orig = case / "0.orig"
    if orig.is_dir():
        out += sorted(p for p in orig.rglob("*") if p.is_file())
    return out


def shift_settings(case: Path, shift) -> dict:
    """Apply the shift to every point-based setting (D5.2). Returns
    {"shifted": [...], "unclassifiedVectors": [...]}."""
    shifted, unknown = [], []
    for f in shift_files(case):
        try:
            text = f.read_text()
        except UnicodeDecodeError:
            continue
        new, changes = shift_text(text, shift)
        rel = str(f.relative_to(case))
        if changes:
            f.write_text(new)
            shifted += [dict(c, file=rel) for c in changes]
        if not rel.startswith("0.orig"):
            unknown += [f"{rel}:{k}" for k in unclassified_vectors(new)]
    return {"shifted": shifted, "unclassifiedVectors": unknown}


# --------------------------------------------------------------------------- #
# SP solver settings (D-064 findings, lead 2026-09-23)
# --------------------------------------------------------------------------- #

# Precision-profile keywords set explicitly by the templates: removed in SP
# so that coupled.precisionProfile auto applies (D7, D-064 deviation 2)
PROFILE_ENTRIES = ("coupled.convergence.residualTol",
                   "solvers.coupled.tolerance",
                   "solvers.coupled.etaMin",
                   "coupled.bounds.omegaMin",
                   "coupled.bounds.kMin")
# Segregated (U|k|omega) solver tolerance in SP: 1e-10 hits the float floor
# (k solve 1000 sweeps per iteration, T1 4x slower; D-064)
SP_SEGREGATED_TOL = 1e-6
_SOLVER_BLOCK_RE = re.compile(
    r'(?s)(\n[ \t]*("[^"\n]*"|[A-Za-z_]\w*)[ \t]*\n[ \t]*\{)(.*?)(\n[ \t]*\})')


def _is_turb_or_u(key: str) -> bool:
    k = key.strip('"')
    try:
        return any(re.fullmatch(k, f) for f in ("U", "k", "omega"))
    except re.error:
        return False


def sp_solver_settings(case: Path) -> dict:
    """SP settings of system/fvSolution: (U|k|omega) solver tolerance raised
    to SP_SEGREGATED_TOL (text edit, only if smaller), the five profile
    keywords removed (foamDictionary -disableFunctionEntries -remove)."""
    f = case / "system" / "fvSolution"
    text = f.read_text()
    changed = []

    def blk(m: re.Match) -> str:
        if not _is_turb_or_u(m.group(2)):
            return m.group(0)

        def tol(t: re.Match) -> str:
            v = float(t.group(2))
            if v < SP_SEGREGATED_TOL:
                changed.append({"solver": m.group(2), "old": v,
                                "new": SP_SEGREGATED_TOL})
                return f"{t.group(1)}{SP_SEGREGATED_TOL:g};"
            return t.group(0)

        body = re.sub(rf"(\btolerance\s+)({_NUM})\s*;", tol, m.group(3))
        return m.group(1) + body + m.group(4)

    f.write_text(_SOLVER_BLOCK_RE.sub(blk, text))
    removed = []
    for e in PROFILE_ENTRIES:
        rc = cfenv.run(["foamDictionary", "-disableFunctionEntries", "-entry",
                        e, "-remove", "system/fvSolution"], cwd=case,
                       log=case / "log.foamDictionary.sp", nice=False)
        if rc == 0:
            removed.append(e)
    return {"segregatedTolerance": changed, "profileEntriesRemoved": removed}


# --------------------------------------------------------------------------- #
# controlDict text edits (no foamDictionary: it would re-quantise scalars)
# --------------------------------------------------------------------------- #

def _get_keyword(f: Path, key: str) -> str | None:
    m = re.search(rf"(?m)^{key}\s+([^;]+);", f.read_text())
    return m.group(1).strip() if m else None


def _set_keyword(f: Path, key: str, value: str) -> None:
    t = f.read_text()
    t2, n = re.subn(rf"(?m)^({key}\s+)[^;]+;", rf"\g<1>{value};", t)
    if n != 1:
        raise RuntimeError(f"{f}: top-level keyword {key} not found once")
    f.write_text(t2)


# --------------------------------------------------------------------------- #
# mesh helpers
# --------------------------------------------------------------------------- #

def _processor_dirs(case: Path) -> list[Path]:
    return sorted(p for p in case.glob("processor[0-9]*")
                  if (p / "constant" / "polyMesh").is_dir())


def _serial_mesh(case: Path) -> bool:
    poly = case / "constant" / "polyMesh"
    return any((poly / n).exists() for n in ("owner", "owner.gz"))


def _needs_ascii(case: Path) -> bool:
    """Serial mesh compressed or binary (needs foamFormatConvert in DP)."""
    poly = case / "constant" / "polyMesh"
    if any(poly.glob("*.gz")):
        return True
    pts = poly / "points"
    if not pts.exists():
        return False
    with open(pts, "rb") as fh:
        head = fh.read(2000).decode("latin-1")
    return re.search(r"format\s+binary", head) is not None


def _mpi(nprocs: int, cmd: list[str]) -> list[str]:
    return (cfenv.mpirun_prefix(nprocs) + cmd + ["-parallel"]) if nprocs > 1 \
        else cmd


# --------------------------------------------------------------------------- #
# case preparation
# --------------------------------------------------------------------------- #

def read_info(case: Path) -> dict:
    f = Path(case) / HARNESS_JSON
    return json.loads(f.read_text()) if f.exists() else {}


def _write_info(case: Path, info: dict) -> None:
    (case / HARNESS_JSON).write_text(json.dumps(info, indent=2) + "\n")


def prepare_case(case: Path, args: list[str]) -> list[str] | None:
    """SP preparation of a case before its Allrun (steps 1-4 of the module
    doc). Returns the Allrun arguments (+ -keep-mesh), or None for a
    -mesh-only request (done here, in DP). Raises SPGeometryFail."""
    args = list(args)
    if "-restart" in args:
        return args
    mesh_only = "-mesh-only" in args
    info = read_info(case)
    if info.get("status") == SP_GEOMETRY_FAIL:
        raise SPGeometryFail(f"{case.name}: {SP_GEOMETRY_FAIL}: {info.get('reasons')}")
    if not info.get("prepared"):
        info = _prepare(case, [a for a in args if a != "-mesh-only"])
    if info.get("status") == SP_GEOMETRY_FAIL:
        from . import results   # noqa: PLC0415
        results.write("tests", case.name, {
            "precision": label(), "spGeometry": SP_GEOMETRY_FAIL,
            "reasons": info.get("reasons"), "meshShift": info.get("shift"),
            "checkMeshDiff": info.get("checkMeshDiff"), "pass": None,
            "skipped": True})
        raise SPGeometryFail(f"{case.name}: {SP_GEOMETRY_FAIL}: {info.get('reasons')}")
    if mesh_only:
        return None
    return args + ["-keep-mesh"]


def _prepare(case: Path, args: list[str]) -> dict:
    t0 = time.perf_counter()
    info: dict = {"precision": label(), "prepared": False,
                  "dpBashrc": os.environ.get("CF_DP_BASHRC", DP_BASHRC_DEFAULT)}
    cdict = case / "system" / "controlDict"
    fmt = _get_keyword(cdict, "writeFormat") or "binary"
    info["writeFormatTemplate"] = fmt
    info["writePrecisionDP"] = _get_keyword(cdict, "writePrecision")

    # 1. mesh in DP, ascii
    _set_keyword(cdict, "writeFormat", "ascii")
    try:
        mesh_args = [a for a in args if a != "-keep-mesh"]
        rc = dp_run(["./Allrun", "-mesh-only"] + mesh_args, case,
                    case / "log.Allrun.meshDP")
        if rc != 0:
            raise RuntimeError(f"{case}: DP mesh generation failed (rc {rc}, "
                               "log.Allrun.meshDP)")
        procs = _processor_dirs(case)
        serial = _serial_mesh(case)
        nprocs = 1 if serial else len(procs)
        if serial and _needs_ascii(case):
            rc = dp_run(["foamFormatConvert", "-constant", "-noZero"], case,
                        case / "log.foamFormatConvert.dp")
            if rc != 0:
                raise RuntimeError(f"{case}: foamFormatConvert failed")
            poly = case / "constant" / "polyMesh"
            for gz in poly.glob("*.gz"):
                if (poly / gz.stem).exists():
                    gz.unlink()
        rc = dp_run(_mpi(nprocs, ["checkMesh", "-allGeometry", "-allTopology",
                                  "-constant"]),
                    case, case / "log.checkMesh.dp")
        dp = parse_checkmesh(case / "log.checkMesh.dp")
        dp["rc"] = rc
    finally:
        # 2. template write format for the SP run (binary, D8.3)
        _set_keyword(cdict, "writeFormat", fmt)
    if dp["bbox"] is None:
        raise RuntimeError(f"{case}: no bounding box in log.checkMesh.dp")
    t_dp = time.perf_counter() - t0

    # 3. origin shift in SP
    shift = shift_from_bbox(dp["bbox"])
    rc = cfenv.run(_mpi(nprocs, ["transformPoints", "-translate",
                                 fmt_vec(shift)]),
                   cwd=case, log=case / "log.transformPoints")
    if rc != 0:
        raise RuntimeError(f"{case}: transformPoints failed (rc {rc})")
    write_mesh_shift(case, shift, dp["bbox"])
    settings = shift_settings(case, shift)
    # SP solver settings (D-064): only with a real SP build
    sp_settings = (sp_solver_settings(case) if build_precision() == "SP"
                   else {})

    # 4. checkMesh gate in SP
    rc = cfenv.run(_mpi(nprocs, ["checkMesh", "-allGeometry", "-allTopology",
                                 "-constant"]),
                   cwd=case, log=case / "log.checkMesh.sp")
    sp = parse_checkmesh(case / "log.checkMesh.sp")
    sp["rc"] = rc
    for p in [case / "constant" / "polyMesh" / "sets"] + \
            [d / "constant" / "polyMesh" / "sets" for d in _processor_dirs(case)]:
        if p.is_dir():
            for f in p.iterdir():
                f.unlink()
            p.rmdir()
    ok, why = gate(dp, sp)
    overridden = gate_overridden(dp, sp)
    info.update({
        "prepared": True,
        "status": (SP_GEOMETRY_FAIL if not ok
                   else "ok-overridden" if overridden else "ok"),
        "reasons": why,
        "gateOverride": gate_overrides(),
        "gateOverridden": overridden,
        "nProcsMesh": nprocs,
        "shift": shift,
        "bboxCentreDP": [-s for s in shift],
        "checkMeshDP": dp, "checkMeshSP": sp,
        "checkMeshDiff": checkmesh_diff(dp, sp),
        "shiftedEntries": settings["shifted"],
        "spSolverSettings": sp_settings,
        "unclassifiedVectors": settings["unclassifiedVectors"],
        "wallSecondsDPMesh": t_dp,
        "wallSecondsPrepare": time.perf_counter() - t0,
    })
    _write_info(case, info)
    return info


# --------------------------------------------------------------------------- #
# D11 record fields
# --------------------------------------------------------------------------- #

def _first_nstat(log: Path) -> int | None:
    if not log.exists():
        return None
    from . import logs   # noqa: PLC0415
    rows = logs.parse_cf(log)
    v = rows[0].get("nStat") if rows else None
    return int(v) if isinstance(v, (int, float)) else None


def annotate(rec: dict, case: Path, metrics: dict | None = None,
             dp_rec: dict | None = None, dp_case: Path | None = None) -> dict:
    """Add the D11 fields to a result record (in place; returned).

    precision        "dp" | "sp"
    meshShift        [sx, sy, sz] (SP) | None
    checkMeshDiff    [{"check", "only"}] (SP) | None
    spGeometry       "ok" | "SP-geometry-fail" (SP) | None
    staticSetSize    size of the static remediation set of this run
    staticSetSizeDP, staticSetSizeDiff (SP - DP)
    <m>_rel_to_DP    (value - DP value)/|DP value| for every metric m
    spVsDpPass       all |<m>_rel_to_DP| <= TOL_SP_VS_DP (D10), None if n/a
    The DP counterpart is the record results/tests/<dp name>.json (dp_rec
    overrides) and the run directory <RUN_ROOT>/<dp name> (dp_case)."""
    metrics = metrics or {}
    case = Path(case)
    rec["precision"] = label()
    rec["buildPrecision"] = build_precision()
    if rec.get("staticSetSize") is None:
        rec["staticSetSize"] = _first_nstat(case / "log.coupledFoam")
    if not is_sp():
        rec.update({"meshShift": None, "checkMeshDiff": None,
                    "spGeometry": None, "staticSetSizeDP": None,
                    "staticSetSizeDiff": None, "spVsDpPass": None})
        for k in metrics:
            rec[f"{k}_rel_to_DP"] = None
        return rec

    info = read_info(case)
    rec["meshShift"] = info.get("shift")
    rec["checkMeshDiff"] = info.get("checkMeshDiff")
    rec["spGeometry"] = info.get("status")
    rec["spHarness"] = {k: info.get(k) for k in (
        "bboxCentreDP", "shiftedEntries", "unclassifiedVectors",
        "wallSecondsDPMesh", "wallSecondsPrepare", "reasons",
        "gateOverridden", "spSolverSettings")}
    rec["spHarness"]["checkMeshFailedDP"] = (info.get("checkMeshDP") or {}).get("failed")
    rec["spHarness"]["checkMeshFailedSP"] = (info.get("checkMeshSP") or {}).get("failed")

    name = dp_name(case.name)
    if dp_rec is None:
        from . import results   # noqa: PLC0415
        dp_rec = results.read_exact("tests", name) or {}
    if dp_case is None:
        dp_case = cfenv.RUN_ROOT / name
    rec["dpCounterpart"] = {"record": name if dp_rec else None,
                            "case": str(dp_case) if dp_case.exists() else None}

    dp_size = dp_rec.get("staticSetSize")
    if dp_size is None:
        dp_size = _first_nstat(dp_case / "log.coupledFoam")
    rec["staticSetSizeDP"] = dp_size
    sp_size = rec.get("staticSetSize")
    rec["staticSetSizeDiff"] = (sp_size - dp_size
                                if None not in (sp_size, dp_size) else None)

    rels = []
    for k, v in metrics.items():
        ref = dp_rec.get(k)
        if v is None or ref in (None, 0):
            rec[f"{k}_rel_to_DP"] = None
            rels.append(None)
            continue
        r = (float(v) - float(ref)) / abs(float(ref))
        rec[f"{k}_rel_to_DP"] = r
        rec[f"{k}_DP"] = ref
        rels.append(r)
    rec["spVsDpTol"] = TOL_SP_VS_DP
    rec["spVsDpPass"] = (all(abs(r) <= TOL_SP_VS_DP for r in rels)
                         if rels and None not in rels else None)
    return rec

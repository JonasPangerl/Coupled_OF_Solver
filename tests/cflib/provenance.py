"""Run provenance and the report staleness guard (OPUS_TASKS TASK 6).

Provenance
    cflib.case.allrun() writes <case>/provenance.json at the start of every
    run: git commit of the repository and a dirty flag, path, mtime, size
    and sha256 of the coupledFoam binary and of the libcoupledFoam.so it
    actually loads (`which coupledFoam`, `ldd`), a build id derived from
    both hashes, host, date, the OpenFOAM build (WM_OPTIONS) and every CF_*
    environment variable. A second run in the same directory (restart)
    keeps the earlier records under "previous".

Staleness guard
    bench/make_report.py (and the figure helpers it calls) accept a result
    or a run directory only if it was produced at the target commit
    (--commit, default `git rev-parse HEAD`) from a clean tree. Everything
    else is rejected ("pending" in the papers) and listed in the generated
    appendix "Missing or stale results" (tables/missing_results.tex).
    --allow-stale accepts it anyway, still listed, with a warning on the
    title page.

    simpleFoam reference directories (run/ref_*) are exempt from the run
    check: they are produced by the untouched system OpenFOAM, are cached
    on purpose and are not re-run in the final re-run (TASK 6 step 2).
    Their numbers enter the report through the test records, which are
    guarded.

    Build and record identity (review M9, D-068):
    - a coupledFoam run must have used ONE build (every start) equal to the
      report's build (Guard.build_id: $CF_REPORT_BUILD_ID, else the install
      on PATH, else the first run checked), and the binary / library at the
      recorded paths must not have been rebuilt since (binaries_changed);
    - a result record carries the fingerprint of its run directory
      (results.write: runFingerprint); check_record rejects the record when
      the directory was re-run or its solver log changed afterwards, or the
      run started at another commit than the record's; check_record_run
      does the same for a consumer that reads a record's run directory
      (older records: the directory was started after the record was
      written).

Standard library only: render_fields.py imports this under pvbatch.
"""

from __future__ import annotations

import functools
import hashlib
import json
import os
import re
import shutil
import socket
import subprocess
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
PROVENANCE = "provenance.json"

# Tracked paths whose modification does not change what a run computes:
# generated results and report outputs (written while the tests run)
DIRTY_IGNORE = ("results", "report")


# --------------------------------------------------------------------------- #
# git
# --------------------------------------------------------------------------- #

def _git(*args: str) -> str | None:
    try:
        # rstrip only: porcelain lines start with a status column (" M")
        return subprocess.run(["git", *args], cwd=REPO, capture_output=True,
                              text=True, check=True).stdout.rstrip()
    except (FileNotFoundError, subprocess.CalledProcessError):
        return None


def git_head() -> str | None:
    """Full commit hash of HEAD."""
    return _git("rev-parse", "HEAD")


def git_resolve(rev: str) -> str | None:
    """Full commit hash of a revision (sha prefix, tag, branch)."""
    return _git("rev-parse", "--verify", "--quiet", f"{rev}^{{commit}}")


def git_dirty_files() -> list[str]:
    """Tracked files with uncommitted changes, outside DIRTY_IGNORE."""
    out = _git("status", "--porcelain", "--untracked-files=no", "--", ".",
               *[f":(exclude){p}" for p in DIRTY_IGNORE])
    return [line[3:] for line in out.splitlines()] if out else []


def git_porcelain() -> list[str]:
    out = _git("status", "--porcelain")
    return out.splitlines() if out else []


# --------------------------------------------------------------------------- #
# binaries
# --------------------------------------------------------------------------- #

@functools.lru_cache(maxsize=None)
def _sha256(path: str, mtime: float, size: int) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def file_info(path: str | Path | None) -> dict | None:
    if not path:
        return None
    p = Path(path)
    if not p.exists():
        return {"path": str(p), "exists": False}
    real = p.resolve()
    st = real.stat()
    return {
        "path": str(p), "resolved": str(real),
        "mtime": time.strftime("%Y-%m-%dT%H:%M:%S",
                               time.localtime(st.st_mtime)),
        "mtimeEpoch": st.st_mtime, "size": st.st_size,
        "sha256": _sha256(str(real), st.st_mtime, st.st_size),
    }


def ldd_library(binary: str | Path, name: str) -> str | None:
    """Path of shared library `name` that `binary` loads (ldd)."""
    try:
        out = subprocess.run(["ldd", str(binary)], capture_output=True,
                             text=True, check=False).stdout
    except FileNotFoundError:
        return None
    for line in out.splitlines():
        m = re.match(rf"\s*{re.escape(name)}\S*\s+=>\s+(\S+)", line)
        if m and m.group(1) != "not":
            return m.group(1)
    return None


def binary_record(solver: str = "coupledFoam") -> dict:
    exe = shutil.which(solver)
    rec = {"which": exe, "binary": file_info(exe)}
    if exe and solver == "coupledFoam":
        rec["libcoupledFoam"] = file_info(ldd_library(exe, "libcoupledFoam.so"))
    return rec


def build_id(binaries: dict) -> str | None:
    """Short hash over the coupledFoam binary and library contents."""
    cf = binaries.get("coupledFoam") or {}
    parts = [((cf.get("binary") or {}).get("sha256")),
             ((cf.get("libcoupledFoam") or {}).get("sha256"))]
    if not any(parts):
        return None
    return hashlib.sha256("|".join(p or "-" for p in parts).encode()
                          ).hexdigest()[:16]


# --------------------------------------------------------------------------- #
# provenance.json
# --------------------------------------------------------------------------- #

def record(solver: str | None = None, args: list[str] | None = None) -> dict:
    """Provenance of a run started now from this process."""
    dirty = git_dirty_files()
    porcelain = git_porcelain()
    binaries = {"coupledFoam": binary_record("coupledFoam")}
    if solver and solver != "coupledFoam":
        binaries[solver] = binary_record(solver)
    head = git_head()
    return {
        "gitCommit": head,
        "gitCommitShort": head[:7] if head else None,
        "dirty": bool(dirty),
        "dirtyFiles": dirty[:100],
        "dirtyIgnores": list(DIRTY_IGNORE),
        "gitStatusPorcelain": porcelain[:200],
        "solver": solver, "args": list(args or []),
        "binaries": binaries,
        "buildId": build_id(binaries),
        "hostname": socket.gethostname(),
        "date": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "WM_OPTIONS": os.environ.get("WM_OPTIONS"),
        "WM_PROJECT_DIR": os.environ.get("WM_PROJECT_DIR"),
        "FOAM_USER_APPBIN": os.environ.get("FOAM_USER_APPBIN"),
        "FOAM_USER_LIBBIN": os.environ.get("FOAM_USER_LIBBIN"),
        "env": {k: v for k, v in sorted(os.environ.items())
                if k.startswith("CF_")},
    }


def solver_of(args: list[str] | None) -> str | None:
    args = list(args or [])
    if "-solver" in args:
        i = args.index("-solver")
        if i + 1 < len(args):
            return args[i + 1]
    return None


def write(case: Path, args: list[str] | None = None) -> Path:
    """Write <case>/provenance.json; an existing record (earlier run in the
    same directory, e.g. before a restart) moves to "previous"."""
    rec = record(solver_of(args), args)
    f = Path(case) / PROVENANCE
    old = read(case)
    if old:
        prev = old.pop("previous", [])
        rec["previous"] = (prev + [old])[-20:]
    f.write_text(json.dumps(rec, indent=2) + "\n")
    return f


def read(case: Path) -> dict | None:
    f = Path(case) / PROVENANCE
    try:
        return json.loads(f.read_text())
    except (OSError, json.JSONDecodeError):
        return None


# --------------------------------------------------------------------------- #
# binaries of a run vs the current install (M9)
# --------------------------------------------------------------------------- #

def current_build_id(solver: str = "coupledFoam") -> str | None:
    """Build id of the coupledFoam install on PATH now (None without an
    OpenFOAM environment)."""
    if not shutil.which(solver):
        return None
    return build_id({"coupledFoam": binary_record(solver)})


def binaries_changed(prov: dict) -> list[str]:
    """Binaries of a run's provenance whose file at the recorded (resolved)
    path now has another content (rebuilt since the run) or is gone."""
    out = []
    cf = (prov.get("binaries") or {}).get("coupledFoam") or {}
    for key in ("binary", "libcoupledFoam"):
        info = cf.get(key) or {}
        path = info.get("resolved") or info.get("path")
        if not path or not info.get("sha256"):
            continue
        now = file_info(path)
        if not now or not now.get("exists", True):
            out.append(f"{key} {path} is gone")
        elif now.get("sha256") != info["sha256"]:
            out.append(f"{key} {path} was rebuilt after the run")
    return out


# --------------------------------------------------------------------------- #
# record <-> run directory fingerprint (M9)
# --------------------------------------------------------------------------- #

FINGERPRINT_LOGS = ("log.coupledFoam", "log.simpleFoam", "log.potentialFoam")


def run_fingerprint(run_dir: Path | str | None) -> dict | None:
    """Identity of a run directory at the time a result record is written:
    start date, commit and build id of its provenance and the size of its
    solver logs. Stored in the record (results.write, runFingerprint); a
    later re-run in the same directory changes it."""
    if not run_dir:
        return None
    d = Path(run_dir)
    if not d.is_dir():
        return None
    prov = read(d) or {}
    logs_ = {}
    for name in FINGERPRINT_LOGS:
        f = d / name
        if f.exists():
            logs_[name] = {"size": f.stat().st_size}
    return {"runDir": str(d.resolve()), "provenanceDate": prov.get("date"),
            "provenanceCommit": prov.get("gitCommit"),
            "buildId": prov.get("buildId"), "logs": logs_}


def fingerprint_mismatch(fp: dict | None) -> str | None:
    """Why the run directory of a record no longer is the run the record
    was made from (None: it is, or the directory is gone)."""
    if not fp or not fp.get("runDir"):
        return None
    cur = run_fingerprint(fp["runDir"])
    if cur is None:
        return None
    if fp.get("provenanceDate") != cur.get("provenanceDate"):
        return (f"run directory re-run after the record was written (start "
                f"{fp.get('provenanceDate')} -> {cur.get('provenanceDate')})")
    for name, v in (fp.get("logs") or {}).items():
        c = cur["logs"].get(name)
        if c is None or c.get("size") != v.get("size"):
            return f"{name} changed after the record was written"
    return None


def record_run_mismatch(rec: dict, run_dir: Path | str) -> str | None:
    """fingerprint_mismatch of the record's own fingerprint; for records
    written before fingerprints: the run directory was started after the
    record was written (provenance date > record timestamp)."""
    fp = rec.get("runFingerprint")
    if fp:
        return fingerprint_mismatch(fp)
    if not run_dir or not Path(run_dir).is_dir():
        return None
    prov = read(Path(run_dir))
    ts, start = rec.get("timestamp"), (prov or {}).get("date")
    if ts and start and str(start) > str(ts):
        return (f"run directory started {start}, after the record was "
                f"written ({ts})")
    if ts:
        # no provenance either: a solver log written after the record
        # (60 s slack) belongs to a later run in the same directory
        try:
            t_rec = time.mktime(time.strptime(str(ts)[:19],
                                              "%Y-%m-%dT%H:%M:%S"))
        except ValueError:
            return None
        for name in FINGERPRINT_LOGS:
            f = Path(run_dir) / name
            if f.exists() and f.stat().st_mtime > t_rec + 60:
                return (f"{name} written "
                        + time.strftime("%Y-%m-%dT%H:%M:%S",
                                        time.localtime(f.stat().st_mtime))
                        + f", after the record ({ts})")
    return None


# --------------------------------------------------------------------------- #
# staleness guard
# --------------------------------------------------------------------------- #

def _short(c) -> str:
    return str(c)[:7] if c else "none"


class Guard:
    """Decides whether a result record or a run directory may enter the
    report and collects everything it rejected (or accepted as stale)."""

    def __init__(self, target: str | None, allow_stale: bool = False,
                 enabled: bool = True, build_id: str | None = "auto"):
        """build_id: the coupledFoam build every guarded coupledFoam run
        must have used (M9). "auto": $CF_REPORT_BUILD_ID, else the install
        on PATH now (current_build_id); if neither exists, the first
        coupledFoam run checked sets it (all runs of a report must share
        one build). None: no build-id comparison."""
        self.target = target
        self.allow_stale = allow_stale
        self.enabled = enabled and bool(target)
        if build_id == "auto":
            build_id = (os.environ.get("CF_REPORT_BUILD_ID")
                        or (current_build_id() if self.enabled else None)
                        or "first-run")
        self.build_id = build_id
        # (case, item) -> {"found", "date", "reason", "status"}
        self.items: dict[tuple[str, str], dict] = {}
        self._runs: dict[str, bool] = {}

    # -- matching ---------------------------------------------------------- #
    def matches(self, commit) -> bool:
        if not self.enabled:
            return True
        if not commit:
            return False
        c = str(commit)
        if c.endswith("-dirty"):
            return False
        c = c.split("-")[0].lower()
        t = self.target.lower()
        return len(c) >= 7 and (t.startswith(c) or c.startswith(t))

    def _reject(self, case: str, item: str, found, date, reason: str) -> bool:
        key = (case, item)
        if key not in self.items:
            self.items[key] = {
                "case": case, "item": item, "found": found or "none",
                "date": date or "", "reason": reason,
                "status": "used (stale)" if self.allow_stale else "pending",
            }
        return self.allow_stale

    def missing(self, case: str, item: str, reason: str = "no result") -> None:
        """An expected result that does not exist at all."""
        if self.enabled:
            key = (case, item)
            self.items.setdefault(key, {
                "case": case, "item": item, "found": "none", "date": "",
                "reason": reason, "status": "pending"})

    def historical(self, case: str, item: str, found, date,
                   reason: str) -> None:
        """Listed, but accepted on purpose (--allow-historical)."""
        self.items.setdefault((case, item), {
            "case": case, "item": item, "found": found or "none",
            "date": date or "", "reason": reason,
            "status": "used (historical)"})

    # -- checks ------------------------------------------------------------ #
    def check_record(self, case: str, rec: dict,
                     item: str = "result record") -> bool:
        """A results/*.json record: its gitCommit must be the target."""
        if not self.enabled:
            return True
        commit = rec.get("gitCommit")
        if self.matches(commit):
            # M9: the record's run directory must still be the run the
            # record was made from, started at the same commit
            fp = rec.get("runFingerprint") or {}
            why = fingerprint_mismatch(fp)
            if why is None and fp.get("provenanceCommit") and \
                    not self.matches(fp["provenanceCommit"]):
                why = (f"run started at {_short(fp['provenanceCommit'])}, "
                       "record written at another commit")
            if why is None:
                return True
            return self._reject(case, item, commit, rec.get("timestamp"), why)
        if not commit:
            reason = "record carries no gitCommit"
        elif str(commit).endswith("-dirty"):
            reason = "written from a modified tree"
        else:
            reason = "record from another commit"
        return self._reject(case, item, commit, rec.get("timestamp"), reason)

    def check_run(self, run_dir: Path, item: str = "run directory",
                  case: str | None = None) -> bool:
        """A run directory: provenance.json must name the target commit for
        every start in it (a restart continues an earlier run) and a clean
        tree."""
        if not self.enabled:
            return True
        run_dir = Path(run_dir)
        key = str(run_dir)
        if key in self._runs:
            ok = self._runs[key]
            return ok or self.allow_stale
        prov = read(run_dir)
        case = case or run_dir.name
        if prov is None:
            date = None
            logs_ = sorted(run_dir.glob("log.*"), key=lambda p: p.stat().st_mtime)
            if logs_:
                date = time.strftime("%Y-%m-%dT%H:%M:%S",
                                     time.localtime(logs_[-1].stat().st_mtime))
            self._runs[key] = False
            return self._reject(case, item, None, date,
                                "run has no provenance.json")
        starts = (prov.get("previous") or []) + [prov]
        commits = {s.get("gitCommit") for s in starts}
        dirty = any(s.get("dirty") for s in starts)
        ok = (len(commits) == 1 and self.matches(prov.get("gitCommit"))
              and not dirty)
        if ok:
            why = self.build_problem(starts)
            if why:
                self._runs[key] = False
                return self._reject(case, item, prov.get("buildId"),
                                    prov.get("date"), why)
        self._runs[key] = ok
        if ok:
            return True
        found = ",".join(_short(c) for c in sorted(commits, key=str))
        if dirty:
            found += "-dirty"
            reason = "run started from a modified tree"
        elif len(commits) > 1:
            reason = "run continued across commits"
        else:
            reason = "run from another commit"
        return self._reject(case, item, found, prov.get("date"), reason)

    def build_problem(self, starts: list[dict]) -> str | None:
        """M9: the coupledFoam binary of a run (every start of it) must be
        the report's build (self.build_id) and must not have been rebuilt
        since. None if fine or not a coupledFoam run."""
        cf_starts = [s for s in starts if s.get("solver") == "coupledFoam"]
        if not cf_starts:
            return None
        ids = {s.get("buildId") for s in cf_starts}
        if None in ids:
            return "run without a recorded coupledFoam build id"
        if len(ids) > 1:
            return "run continued across coupledFoam builds"
        bid = ids.pop()
        if self.build_id == "first-run":
            self.build_id = bid
        if self.build_id and bid != self.build_id:
            return (f"coupledFoam build {bid} is not the report's build "
                    f"{self.build_id} (stale or foreign binary)")
        changed = binaries_changed(cf_starts[-1])
        if changed:
            return "stale solver binary: " + "; ".join(changed)
        return None

    def check_record_run(self, case: str, rec: dict, run_dir: Path,
                         item: str = "run directory of the record") -> bool:
        """A result record together with the run directory a consumer is
        about to read for it: the directory must pass check_run (unless
        it is a simpleFoam reference, run/ref_*) and must still be the run
        the record was made from (record_run_mismatch)."""
        if not self.enabled:
            return True
        run_dir = Path(run_dir)
        if not run_dir.name.startswith("ref_") and \
                not self.check_run(run_dir, case=case):
            return self.allow_stale
        why = record_run_mismatch(rec, run_dir)
        if why is None:
            return True
        return self._reject(case, item, rec.get("gitCommit"),
                            rec.get("timestamp"), why)

    def run_commit(self, run_dir: Path) -> str | None:
        prov = read(run_dir)
        return (prov or {}).get("gitCommit")

    # -- output ------------------------------------------------------------ #
    def rows(self) -> list[dict]:
        return sorted(self.items.values(), key=lambda r: (r["case"], r["item"]))

    @property
    def n_rejected(self) -> int:
        return sum(1 for r in self.items.values()
                   if not r["status"].startswith("used (historical"))


# The guard of the current report run (set by bench/make_report.py, used by
# the figure helpers it imports); a permissive guard when none is set
ACTIVE: Guard | None = None


def active() -> Guard:
    return ACTIVE if ACTIVE is not None else Guard(None, enabled=False)

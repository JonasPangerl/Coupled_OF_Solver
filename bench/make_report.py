#!/usr/bin/env python3
"""Report generator (spec 15): everything from results/ and the run logs.

Outputs
    report/REPORT.md                 Markdown report, sections as spec 15
    report/figures/<name>.png        150 dpi raster figures for REPORT.md
    report/paper/figures/<name>.pdf  vector figures for the LaTeX paper
    report/paper/tables/<name>.tex   LaTeX tables (booktabs)
    report/paper/numbers.tex         every number quoted in the paper as a
                                     macro (\\cfNumber...), so the text never
                                     carries hand-copied values

Inputs: results/tests/*.json, results/bench/*.json (current records only,
see bench/run_bench.py:is_current), results/exploratory/*.json (written by
bench/exploratory_numbers.py), and read-only logs under run/.

Every performance number is given as wall-clock time AND CPU-hours (user
directive: the fastest solver, not the fewest iterations).

Figures whose data is missing are skipped with a note (the report lists
them), so the generator runs at every stage of the project.

Staleness guard (OPUS_TASKS TASK 6): every table row, figure and number
comes only from results whose gitCommit is the target commit (--commit,
default `git rev-parse HEAD`) and from run directories whose
provenance.json names that commit (tests/cflib/provenance.py; simpleFoam
reference directories are exempt, they are produced by the system
OpenFOAM). Everything else is left out - the papers show their "pending"
placeholders - and is listed in tables/missing_results.tex, the appendix
"Missing or stale results" of both papers. Before generating, the guard
removes every previously generated figure, table and numbers.tex, so no
output of an earlier run survives. --allow-stale restores the old
behaviour (use whatever exists, keep earlier outputs); the stale items are
still listed and both title pages carry a warning. results/exploratory
(historical studies of run/exp_*, no commit recorded) are stale unless
--allow-historical is given (then used and listed as historical).

Usage: bench/make_report.py [--commit SHA] [--allow-stale]
                            [--allow-historical]
       (no OpenFOAM environment needed)
"""

from __future__ import annotations

import argparse
import functools
import json
import subprocess
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "tests"))
sys.path.insert(0, str(REPO / "bench"))
from cflib import logs, post, provenance  # noqa: E402
import run_bench  # noqa: E402
import user_convergence as ucv  # noqa: E402

RESULTS = REPO / "results"
RUN = REPO / "run"
REPORT = REPO / "report"
FIG_PNG = REPORT / "figures"
PAPER = REPORT / "paper"
FIG_PDF = PAPER / "figures"
TABLES = PAPER / "tables"

# Colours: one pair everywhere (Okabe-Ito, colour-blind safe): simpleFoam
# vermillion, coupledFoam blue (was grey/blue; grey curves were too faint)
C_NATIVE = "#D55E00"
C_NATIVE2 = "#E69F00"
C_COUPLED = "#0072B2"
C_COUPLED2 = "#56B4E9"
CONFIG_COLOR = {"A": C_NATIVE2, "B": C_NATIVE, "C": C_COUPLED, "D": C_COUPLED2,
                "E": "#6baed6", "F": "#08306b", "G": "#d95f02", "H": "#9467bd"}
CONFIG_LABEL = {
    "A": "simpleFoam (tutorial)",
    "B": "simpleFoam SIMPLEC",
    "C": "coupledFoam",
    "D": "coupledFoam, blockDiagonal",
    "E": "coupledFoam, fixed V-cycle",
    "F": "coupledFoam, fixed relTol",
    "G": "coupledFoam, Anderson",
    "H": "coupledFoam, fixed K-cycle",
}

# Amendment B7 memory budget (spec constants, not measurements): total
# 105-120 GB at 45 M cells, worst case, i.e. per cell in bytes
B7_CELLS = 45e6
B7_TOTAL_GB = (105.0, 120.0)
B7_BYTES_PER_CELL = tuple(g * 1e9 / B7_CELLS for g in B7_TOTAL_GB)

plt.rcParams.update({
    "figure.dpi": 150, "savefig.dpi": 150, "font.size": 9,
    "axes.grid": True, "grid.alpha": 0.3, "legend.frameon": False,
    "pdf.fonttype": 42, "ps.fonttype": 42,
})

notes: list[str] = []
numbers: dict[str, str] = {}
figures: list[tuple[str, str]] = []

# --------------------------------------------------------------------------- #
# staleness guard (TASK 6)
# --------------------------------------------------------------------------- #

# The guard of this report run (set in main from --commit/--allow-stale)
GUARD = provenance.Guard(None, enabled=False)
ALLOW_HISTORICAL = False

# Result records the final re-run produces (TASK 6 step 3; fnmatch
# patterns over results/tests/*.json stems). A pattern without any record
# is listed as missing in the appendix.
EXPECTED_TESTS = [
    "T0_Re100_np1", "T0_Re100_np4", "T0_Re1000_np1", "T0_Re1000_np4",
    "T1_np1", "T1_np4", "T2_np1", "T2_np4",
    "T3_kOmegaSST_np1", "T3_kOmegaSST_np4", "T3_GEKO_np1", "T3_GEKO_np4",
    "T4a_np*", "T4b_np*", "T5*_np*",
    "T-restart_T1", "T-restart_T3-SST", "T-fpe_*", "T_scaling_T4a",
    "diagnostics_*", "Test-*", "test_env",
]
# results/exploratory records the report quotes
EXPLORATORY_USED = ("re1000", "linsolver", "symbol", "decisions")


def guard_tests(tests: dict) -> dict:
    """Test records at the target commit (all of them with --allow-stale);
    rejected and missing ones go to the appendix."""
    import fnmatch  # noqa: PLC0415
    ok = {n: d for n, d in tests.items() if GUARD.check_record(n, d)}
    for pat in EXPECTED_TESTS:
        if not any(fnmatch.fnmatch(n, pat) for n in tests):
            GUARD.missing(pat.replace("*", "<n>") if pat.endswith("_np*")
                          else pat, "result record")
    return ok


def guard_bench(rows: list[dict]) -> list[dict]:
    ok = [d for d in rows if GUARD.check_record(
        f"bench {d.get('case')} {d.get('config')} run {d.get('run', '?')}", d)]
    if not rows:
        GUARD.missing("benchmark A-H", "results/bench/*.json")
    return ok


def run_ok(case: Path) -> bool:
    """A coupledFoam run directory may be read (provenance at the target
    commit, or --allow-stale)."""
    return GUARD.check_run(case)


def purge_generated() -> int:
    """Remove every output of earlier report runs (figures, tables,
    numbers.tex, REPORT.md figures) so that nothing of another commit can
    survive: what is not regenerated becomes 'pending'. The ParaView
    renders are pruned separately (prune_renders), because pvbatch may be
    unavailable."""
    n = 0
    for d, pats in ((FIG_PDF, ("*.pdf", "*.png", "*.tex", "*.json")),
                    (FIG_PNG, ("*.png",)), (TABLES, ("*.tex",))):
        if not d.is_dir():
            continue
        for pat in pats:
            for f in d.glob(pat):
                if f.name.startswith("render_"):
                    continue
                f.unlink()
                n += 1
    for f in (PAPER / "numbers.tex",):
        if f.exists():
            f.unlink()
            n += 1
    return n


RENDER_PROVENANCE = "render_provenance.json"


def prune_renders() -> None:
    """ParaView renders: keep a render_*.png only if
    figures/render_provenance.json (written by bench/render_fields.py)
    says its coupledFoam panels came from a run at the target commit, or it
    has no coupledFoam panels. Without --allow-stale the others are
    deleted (pending); with it they are kept and listed."""
    pf = FIG_PDF / RENDER_PROVENANCE
    try:
        prov = json.loads(pf.read_text()) if pf.exists() else {}
    except json.JSONDecodeError:
        prov = {}
    for f in sorted(FIG_PDF.glob("render_*.png")):
        p = prov.get(f.name)
        if p is None:
            ok = GUARD.check_record(f.stem, {}, "ParaView render "
                                    "(no render provenance)")
        elif not p.get("cfUsed"):
            ok = True       # simpleFoam reference and geometry only
        else:
            ok = GUARD.check_record(
                f.stem, {"gitCommit": p.get("cfCommit"),
                         "timestamp": p.get("date")},
                f"ParaView render of run/{p.get('cfRun')}")
        if not ok:
            f.unlink()
            notes.append(f"render {f.name}: stale, removed (pending)")


def write_missing_table() -> None:
    """tables/missing_results.tex: the appendix 'Missing or stale results'
    of both papers (always written, also when nothing is missing)."""
    TABLES.mkdir(parents=True, exist_ok=True)
    tgt = (GUARD.target or "unknown")[:7]
    rows = GUARD.rows()
    mode = ("\\texttt{-{}-allow-stale}: the items marked ``used (stale)'' "
            "ARE INCLUDED in this document although they do not come from "
            "the target commit" if GUARD.allow_stale else
            "items not produced at the target commit are left out and "
            "shown as ``pending'' in the text")
    head = [
        "% Generated by bench/make_report.py - do not edit",
        rf"Target commit \texttt{{{tgt}}} (\texttt{{{tex_escape(GUARD.target or '')}}}); "
        + mode + ". A result counts when its record "
        r"(\texttt{results/*/<name>.json}, field \texttt{gitCommit}) or its "
        r"run directory (\texttt{provenance.json}, written by the test "
        r"harness at the start of every run) names this commit and a clean "
        r"source tree. The \code{simpleFoam} reference runs are exempt from "
        r"the run check (untouched system OpenFOAM, cached). "
        rf"Items listed: {len(rows)}.",
        "",
    ]
    if not rows:
        body = [r"\noindent All results used in this document were produced "
                rf"at commit \texttt{{{tgt}}}; nothing is missing.", ""]
    else:
        body = [
            r"{\small",
            r"\begin{longtable}{@{}p{0.21\linewidth}p{0.31\linewidth}"
            r"p{0.15\linewidth}p{0.1\linewidth}p{0.1\linewidth}@{}}",
            r"\caption{Missing or stale results (generated).}"
            r"\label{tab:missing}\\",
            r"\toprule case / record & item (reason) & found commit & date "
            r"& status\\ \midrule\endfirsthead",
            r"\toprule case / record & item (reason) & found commit & date "
            r"& status\\ \midrule\endhead",
            r"\bottomrule\endlastfoot",
        ]
        for r in rows:
            found = str(r["found"])
            found = found if found == "none" else found[:15]
            date = str(r["date"] or "").replace("T", " ")[:16]
            body.append(
                " & ".join([
                    r"\texttt{" + tex_escape(r["case"]) + "}",
                    tex_escape(r["item"]) + r" \newline{\footnotesize ("
                    + tex_escape(r["reason"]) + ")}",
                    r"{\footnotesize\texttt{" + tex_escape(found) + "}}",
                    tex_escape(date).replace(" ", r"\newline "),
                    tex_escape(r["status"])]) + r"\\")
        body += [r"\end{longtable}", "}", ""]
    (TABLES / "missing_results.tex").write_text("\n".join(head + body))


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #

def load_json(kind: str) -> dict[str, dict]:
    out = {}
    d = RESULTS / kind
    if d.is_dir():
        for f in sorted(d.glob("*.json")):
            try:
                out[f.stem] = json.loads(f.read_text())
            except json.JSONDecodeError:
                notes.append(f"unreadable {f}")
    return out


def load_bench() -> list[dict]:
    """Current benchmark records; stale ones (configuration changed since
    the run) are excluded and listed in the notes."""
    recs, stale = run_bench.load_current()
    for s in stale:
        notes.append(f"benchmark record {s} is stale (configHash mismatch), "
                     "excluded")
    # failed runs (M3): never a timing; listed as missing in the appendix
    for r in run_bench.load_failed():
        why = "; ".join(map(str, r.get("failure") or [f"rc {r.get('rc')}"]))
        tag = f"bench {r.get('case')} {r.get('config')} run {r.get('run', '?')}"
        notes.append(f"benchmark run {tag} FAILED ({why}), excluded")
        GUARD.missing(tag, "benchmark run", f"run failed: {why}"[:120])
    return recs


def save(fig, name: str, caption: str) -> None:
    FIG_PNG.mkdir(parents=True, exist_ok=True)
    FIG_PDF.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(FIG_PNG / f"{name}.png")
    fig.savefig(FIG_PDF / f"{name}.pdf")
    plt.close(fig)
    figures.append((name, caption))


def num(key: str, value, fmt: str = "{:.3g}") -> str:
    """Register a number for numbers.tex; returns the formatted string."""
    if value is None:
        s = "n/a"
    elif isinstance(value, str):
        s = value
    else:
        s = fmt.format(value)
    macro = "cf" + "".join(p[:1].upper() + p[1:] for p in
                           key.replace("-", " ").replace("_", " ").split())
    words = ["Zero", "One", "Two", "Three", "Four", "Five", "Six", "Seven",
             "Eight", "Nine"]
    macro = "".join(words[int(ch)] if ch.isdigit() else ch
                    for ch in macro if ch.isalnum())
    numbers[macro] = s
    return s


def git_commit() -> str:
    try:
        return subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=REPO,
                              capture_output=True, text=True,
                              check=True).stdout.strip()
    except (FileNotFoundError, subprocess.CalledProcessError):
        return "unknown"


def tex_escape(s: str) -> str:
    return (str(s).replace("_", r"\_").replace("%", r"\%")
            .replace("&", r"\&").replace("#", r"\#")
            .replace("±", r"$\pm$")
            .replace("†", r"\textsuperscript{\dag}"))


def write_table(name: str, header: list[str], rows: list[list], caption: str,
                label: str, resize: bool = False,
                note: str | None = None) -> str:
    TABLES.mkdir(parents=True, exist_ok=True)
    cols = "l" + "r" * (len(header) - 1)
    lines = [
        r"\begin{table}[tbp]", r"\centering", r"\small",
        rf"\caption{{{caption}}}", rf"\label{{{label}}}",
    ]
    if resize:
        lines.append(r"\resizebox{\linewidth}{!}{%")
    lines += [
        rf"\begin{{tabular}}{{{cols}}}", r"\toprule",
        " & ".join(tex_escape(h) for h in header) + r" \\", r"\midrule",
    ]
    for r in rows:
        lines.append(" & ".join(tex_escape(c) for c in r) + r" \\")
    lines += [r"\bottomrule", r"\end{tabular}"]
    if resize:
        lines.append("}")
    if note:
        lines += [r"\par\smallskip", r"{\footnotesize " + tex_escape(note)
                  + r"\par}"]
    lines += [r"\end{table}", ""]
    (TABLES / f"{name}.tex").write_text("\n".join(lines))
    # Markdown version for REPORT.md
    md = ["| " + " | ".join(header) + " |",
          "|" + "---|" * len(header)]
    md += ["| " + " | ".join(str(c) for c in r) + " |" for r in rows]
    if note:
        md += ["", note]
    return "\n".join(md)


def fmt(v, f="{:.3g}"):
    return "n/a" if v is None else (f.format(v) if not isinstance(v, str) else v)


def pct(v):
    return "n/a" if v is None else f"{100 * v:+.1f}%"


def flat(v) -> str:
    if not v:
        return "n/a"
    return " ".join(fmt(x, "{:.3g}") if isinstance(x, float) else str(x)
                    for x in v)


# --------------------------------------------------------------------------- #
# figures
# --------------------------------------------------------------------------- #

def table_convergence(tests: dict, bench: list[dict], uc: dict) -> None:
    """tables/convergence_choice.tex: automatic and user convergence
    iteration of every force-case run (D-060), and the skeleton
    report/user_convergence_template.json with the automatic iterations as
    hints for filling in report/user_convergence.json."""
    tmpl = REPO / "report" / "user_convergence_template.json"
    tmpl.write_text(json.dumps(ucv.template(tests, bench), indent=2) + "\n")
    rows = []

    def row(label, solver, r, auto_it, case_dir=None):
        u = r.get("user") or {}
        missing = ((r.get("wallToConv_s") or r.get("wall_to_conv_s")) is None
                   or (r.get("Cd") if r.get("Cd") is not None
                       else r.get("Cd_mean")) is None)
        if missing and auto_it and case_dir is not None and case_dir.is_dir():
            # not stored by the record (T3, references): time and the
            # coefficient means at the automatic point
            a = ucv.evaluate(case_dir, solver, auto_it, r)
            r = dict(r)
            for k, v in (("wallToConv_s", a.get("wall_to_conv_s")),
                         ("cpuHoursToConv", a.get("cpu_to_conv_h")),
                         ("Cd", a.get("Cd")), ("Cl", a.get("Cl"))):
                if r.get(k) is None:
                    r[k] = v
        used = "user" if u and not u.get("ignored") else "auto"
        rows.append([
            label, solver, fmt(u.get("iterationsRun") or r.get("iterationsRun")
                               or r.get("iterations"), "{}"),
            fmt(auto_it, "{}"),
            fmt(u.get("iters"), "{}") if u else "-",
            used + (" (user value ignored)" if r.get("userIgnored") else ""),
            fmt(r.get("wallToConv_s") or r.get("wall_to_conv_s")),
            fmt(r.get("cpuHoursToConv") or r.get("cpu_to_conv_h")),
            fmt(r.get("Cd") if r.get("Cd") is not None else r.get("Cd_mean"),
                "{:.4f}"),
            fmt(r.get("Cl") if r.get("Cl") is not None else r.get("Cl_mean"),
                "{:.4f}"),
        ])

    for name in sorted(tests):
        if ucv.case_key(name) is None:
            continue
        r = tests[name]
        row(name, "coupledFoam", r, r.get("autoIteration"), ucv.RUN / name)
        ref = r.get("reference") or {}
        row(name + " ref", "simpleFoam", ref, ref.get("autoIteration"),
            ucv.ref_dir(name))
    for r in sorted(bench, key=lambda d: (d.get("case", ""), d.get("config", ""),
                                          d.get("run", 0))):
        if r.get("case") in ucv.CASES and r.get("convergenceSource") == "user":
            row(f"bench {r['case']} run {r.get('run')}", r["config"], r,
                (r.get("auto") or {}).get("iters_to_conv"))
    if not rows:
        notes.append("convergence choice: no force-case results yet")
        return
    write_table(
        "convergence_choice",
        ["run", "solver / config", "iterations run",
         "automatic (criterion or solver stop)", "user",
         "used", "wall to conv. [s]", "CPU-h to conv.", "Cd", "Cl"],
        rows,
        "Convergence point of the force cases: automatic criterion "
        "(T3: 12.3(ii); T4, T5: D-042 stationary mean) and the iteration "
        "named by the user after inspecting the load histories (D-060). "
        "With a user iteration N, time and CPU-hours count up to N and "
        "Cd, Cl are the means over iterations N to the end of the run.",
        "tab:convchoice", resize=True,
        note=f"{len(uc)} case(s) with user entries. The pass/fail of the "
        "tests is unchanged and uses the automatic criteria. Benchmark rows "
        "appear only where a user value was given.")


def fig_T0_profiles() -> None:
    for re_ in (100, 1000):
        cp = RUN / f"T0_Re{re_}_np1"
        rp = RUN / f"ref_T0_Re{re_}"
        if cp.is_dir() and not run_ok(cp):
            notes.append(f"T0 Re{re_} profiles: run {cp.name} stale, pending")
            continue
        try:
            cv = post.read_xy(post.sets_file(cp, "centreLines", "vertical"))
            ch = post.read_xy(post.sets_file(cp, "centreLines", "horizontal"))
            rv = post.read_xy(post.sets_file(rp, "centreLines", "vertical"))
            rh = post.read_xy(post.sets_file(rp, "centreLines", "horizontal"))
        except (FileNotFoundError, OSError, StopIteration):
            notes.append(f"T0 Re{re_} profiles: data missing")
            continue
        fig, ax = plt.subplots(1, 2, figsize=(6.5, 2.8))
        ax[0].plot(rv[:, 2], rv[:, 0] / 0.1, color=C_NATIVE, lw=2.5,
                   label="simpleFoam")
        ax[0].plot(cv[:, 2], cv[:, 0] / 0.1, color=C_COUPLED, lw=1, ls="--",
                   label="coupledFoam")
        ax[0].set_xlabel("$u/U_{lid}$")
        ax[0].set_ylabel("$y/L$")
        ax[0].legend()
        ax[1].plot(rh[:, 0] / 0.1, rh[:, 3], color=C_NATIVE, lw=2.5)
        ax[1].plot(ch[:, 0] / 0.1, ch[:, 3], color=C_COUPLED, lw=1, ls="--")
        ax[1].set_xlabel("$x/L$")
        ax[1].set_ylabel("$v/U_{lid}$")
        save(fig, f"T0_Re{re_}_profiles",
             f"T0 lid-driven cavity, Re {re_}: centreline profiles, "
             "coupledFoam vs. simpleFoam.")


def fig_histories(tests: dict) -> None:
    for name, d in tests.items():
        h = d.get("history")
        if not h or not h.get("R"):
            continue
        it = np.arange(1, len(h["R"]) + 1)
        fig, ax = plt.subplots(1, 2, figsize=(6.5, 2.6))
        ax[0].semilogy(it, h["R"], color=C_COUPLED, label="coupledFoam $R$")
        case = RUN / name.replace("T0_", "ref_T0_").split("_np")[0] \
            if name.startswith("T0_") else None
        if case is not None and (case / "log.simpleFoam").exists():
            nat = logs.parse_native(case / "log.simpleFoam")
            for fld, ls in (("p", "-"), ("Ux", "--")):
                if fld in nat["residuals"]:
                    r = nat["residuals"][fld]
                    ax[0].semilogy(np.arange(1, len(r) + 1), r, color=C_NATIVE,
                                   ls=ls, lw=0.8, label=f"simpleFoam {fld}")
        ax[0].set_xlabel("iteration")
        ax[0].set_ylabel("residual")
        ax[0].legend(fontsize=7)
        ax[1].semilogy(it, h["CFL"], color=C_COUPLED, label="CFL")
        ax2 = ax[1].twinx()
        ax2.plot(it, h["omega"], color="#d95f02", lw=0.8, label=r"$\omega$")
        ax2.set_ylim(0, 1.05)
        ax2.set_ylabel(r"line-search $\omega$")
        ax[1].set_xlabel("iteration")
        ax[1].set_ylabel("CFL")
        save(fig, f"{name}_history",
             f"{name}: residual, CFL and line-search $\\omega$ histories.")


def _med(rows, c, cfg, key):
    v = [d.get(key) for d in rows
         if d["case"] == c and d["config"] == cfg and d.get(key) is not None]
    return float(np.median(v)) if v else None


def _bars(ax, rows, cases, cfgs, key, ylabel, log=True):
    width = 0.8 / max(1, len(cfgs))
    for j, cfg in enumerate(cfgs):
        vals = [_med(rows, c, cfg, key) or 0 for c in cases]
        ax.bar(np.arange(len(cases)) + j * width, vals, width,
               color=CONFIG_COLOR[cfg], label=CONFIG_LABEL[cfg])
    ax.set_xticks(np.arange(len(cases)) + width * (len(cfgs) - 1) / 2)
    ax.set_xticklabels(cases)
    ax.set_ylabel(ylabel)
    if log:
        ax.set_yscale("log")


def fig_bench(rows: list[dict]) -> str:
    """Wall and CPU time to convergence per case and configuration."""
    if not rows:
        notes.append("benchmark: no results yet")
        return ""
    cases = sorted({d["case"] for d in rows})
    cfgs = [c for c in run_bench.CONFIGS if any(d["config"] == c for d in rows)]

    for key, ylabel, fname in (
        ("wall_to_conv_s", "wall time to convergence [s]", "bench_wall"),
        ("cpu_to_conv_h", "CPU time to convergence [CPU-h]", "bench_cpu"),
        ("iters_to_conv", "iterations to convergence", "bench_iters"),
    ):
        fig, ax = plt.subplots(figsize=(6.5, 3.4))
        _bars(ax, rows, cases, cfgs, key, ylabel)
        ax.legend(fontsize=6, ncol=3, loc="upper center", bbox_to_anchor=(0.5, -0.12))
        save(fig, fname, f"Benchmark: {ylabel} (median of the repeats).")

    # Time per iteration breakdown of coupledFoam: wall and CPU
    br = [(c, _med(rows, c, "C", "t_assembly"), _med(rows, c, "C", "t_linsolve"),
           _med(rows, c, "C", "t_turb"), _med(rows, c, "C", "iters_to_conv"),
           _med(rows, c, "C", "time_per_iter_s"),
           _med(rows, c, "C", "cpu_per_iter_s")) for c in cases]
    br = [b for b in br if b[4] and None not in b[1:4]]
    if br:
        fig, ax = plt.subplots(1, 2, figsize=(6.5, 2.6))
        x = np.arange(len(br))
        a = np.array([b[1] / b[4] for b in br])
        s = np.array([b[2] / b[4] for b in br])
        t = np.array([b[3] / b[4] for b in br])
        tot = a + s + t
        # CPU per iteration split in the proportions of the wall breakdown
        cpu = np.array([b[6] or 0.0 for b in br])
        share = np.divide(cpu, tot, out=np.zeros_like(cpu), where=tot > 0)
        for k, f in enumerate((1.0, share)):
            ax[k].bar(x, a * f, color="#9ecae1", label="assembly")
            ax[k].bar(x, s * f, bottom=a * f, color=C_COUPLED,
                      label="linear solve")
            ax[k].bar(x, t * f, bottom=(a + s) * f, color="#fdae6b",
                      label="turbulence")
            ax[k].set_xticks(x)
            ax[k].set_xticklabels([b[0] for b in br])
        ax[0].set_ylabel("wall time per iteration [s]")
        ax[1].set_ylabel("CPU time per iteration [CPU-s]")
        ax[0].legend(fontsize=7)
        save(fig, "bench_breakdown",
             "coupledFoam time per iteration: assembly, linear solve, "
             "turbulence; wall (left) and CPU summed over ranks (right, "
             "split in the wall-time proportions).")

    fig_memory(rows, cases, cfgs)

    # Executive summary table
    header = ["case", "native best", "native wall [s]", "coupled wall [s]",
              "speed-up (wall)", "native CPU-h", "coupled CPU-h",
              "speed-up (CPU)", "memory ratio"]
    tab = []
    for c in cases:
        nat = [(_med(rows, c, k, "wall_to_conv_s"), k) for k in ("A", "B")
               if _med(rows, c, k, "wall_to_conv_s")]
        if not nat:
            continue
        wn, kbest = min(nat)
        wc = _med(rows, c, "C", "wall_to_conv_s")
        cn = _med(rows, c, kbest, "cpu_to_conv_h")
        cc = _med(rows, c, "C", "cpu_to_conv_h")
        mn = _med(rows, c, kbest, "peakRSS_GB_sum")
        mc = _med(rows, c, "C", "peakRSS_GB_sum")
        tab.append([c, kbest, fmt(wn), fmt(wc), fmt(wn / wc if wc else None),
                    fmt(cn), fmt(cc), fmt(cn / cc if cn and cc else None),
                    fmt(mc / mn if mn and mc else None)])
        num(f"speedup wall {c}", wn / wc if wc else None, "{:.2f}")
        num(f"speedup cpu {c}", cn / cc if cn and cc else None, "{:.2f}")
        num(f"pf wall {c}", _med(rows, c, "C", "preWallSeconds"), "{:.3g}")
        num(f"pf cpuh {c}", _med(rows, c, "C", "preCpuHours"), "{:.2g}")
    return write_table("executive_summary", header, tab,
                       "Executive summary: time to convergence (identical "
                       "criterion for all solvers: the force window, for the "
                       "wake cases T4 and T5 the first stationary window "
                       "mean, D-042; coupledFoam including its "
                       "potentialFoam initialisation), best native "
                       "configuration vs.\\ coupledFoam, wall-clock time and "
                       "CPU-hours.", "tab:summary")


def fig_memory(rows: list[dict], cases: list[str], cfgs: list[str]) -> None:
    """Peak RSS per case and configuration with the B7 budget scaled to the
    case's cell count (coupledFoam)."""
    fig, ax = plt.subplots(figsize=(6.5, 3.4))
    _bars(ax, rows, cases, cfgs, "peakRSS_GB_sum",
          "peak RSS, sum over ranks [GiB]", log=True)
    width = 0.8 / max(1, len(cfgs))
    first = True
    for i, c in enumerate(cases):
        ncell = _med(rows, c, "C", "nCells")
        if not ncell:
            continue
        lo, hi = (b * ncell / 1024**3 for b in B7_BYTES_PER_CELL)
        x0, x1 = i - width / 2, i + width * (len(cfgs) - 0.5)
        ax.fill_between([x0, x1], [lo, lo], [hi, hi], color="k", alpha=0.15,
                        lw=0, label="B7 budget per cell" if first else None)
        first = False
        mc = _med(rows, c, "C", "peakRSS_GB_sum")
        num(f"memory budget ratio {c}", mc / hi if mc else None, "{:.2f}")
    ax.legend(fontsize=6, ncol=3, loc="upper center", bbox_to_anchor=(0.5, -0.12))
    save(fig, "bench_memory",
         "Peak memory (sum of the per-rank maximum RSS) against the B7 "
         "budget (105-120 GB at 45 M cells) scaled per cell (grey band).")


def tests_table(tests: dict) -> str:
    header = ["test", "pass", "iterations", "final R", "wall [s]", "CPU-h",
              "key metric"]
    rows = []
    for name, d in sorted(tests.items()):
        key = ""
        for k in ("l2rel_u", "dpRelDiff", "xrRelDiff", "CdRelDiff",
                  "crossRankRelDiff", "maxInverseError", "relErrorDouble"):
            if k in d and d[k] is not None:
                key = f"{k} = {d[k]:.2e}"
                break
        rows.append([name, "yes" if d.get("pass") else "NO",
                     fmt(d.get("iterations"), "{}")
                     if not isinstance(d.get("iterations"), dict) else "",
                     fmt(d.get("finalR")),
                     fmt(d.get("wallSecondsSolver") or d.get("wallSeconds"))
                     if not isinstance(d.get("wallSeconds"), dict) else "",
                     fmt(d.get("cpuHoursSolver") or d.get("cpuHours"))
                     if not isinstance(d.get("cpuHours"), dict) else "",
                     key])
        num(f"test {name} pass", "yes" if d.get("pass") else "no")
    return write_table("tests", header, rows, "Test results (spec 13).",
                       "tab:tests")


# --------------------------------------------------------------------------- #
# amendment B (15.7) and remaining paper items
# --------------------------------------------------------------------------- #

def fig_eta_rho(tests: dict) -> None:
    """Eisenstat-Walker eta and preconditioner rho histories (15.7)."""
    cand = [(n, d) for n, d in tests.items()
            if (d.get("history") or {}).get("eta")]
    if not cand:
        notes.append("eta/rho history: no data")
        return
    fig, ax = plt.subplots(1, 2, figsize=(6.5, 2.6))
    for name, d in sorted(cand)[:4]:
        h = d["history"]
        eta = [e if e is not None else np.nan for e in h["eta"]]
        it = np.arange(1, len(eta) + 1)
        ax[0].semilogy(it, eta, lw=0.9, label=name)
        rho = [r if (r is not None and r >= 0) else np.nan
               for r in h.get("rho") or []]
        if rho:
            ax[1].plot(np.arange(1, len(rho) + 1), rho, lw=0.9, label=name)
        for ev in d.get("gamgTuneEvents") or []:
            ax[1].axvline(ev.get("iter", 0), color="k", lw=0.5, ls=":")
    ax[0].set_xlabel("iteration")
    ax[0].set_ylabel(r"$\eta_n$")
    ax[1].set_xlabel("iteration")
    ax[1].set_ylabel(r"$\rho$ (first preconditioner application)")
    ax[0].legend(fontsize=6)
    save(fig, "eta_rho_history",
         "Eisenstat-Walker forcing term and preconditioner efficiency; "
         "dotted lines: autoTune events.")


def fig_cycles(tests: dict) -> None:
    """Unit test Test-blockGAMG_cycles (B9): Krylov iterations and solve
    time per cycle type on the block-Poisson system. The JSON (written by
    tests/test_unit.py:_cycle_study) holds per-cycle dicts iterations{},
    wallSeconds{} (solve), processWallSeconds{}, processCpuSeconds{} and
    runs{} (the raw per-cycle records)."""
    for tname, mesh in (("Test-blockGAMG_cycles", "T0 cavity mesh"),
                        ("Test-blockGAMG_cycles_motorBike",
                         "motorBike mesh")):
        d = tests.get(tname)
        if not d:
            notes.append(f"cycle comparison: no {tname} result")
            continue
        runs = d.get("runs") or {}
        its = d.get("iterations") or {c: (runs.get(c) or {}).get("nIterations")
                                       for c in runs}
        wall = d.get("wallSeconds") or {c: (runs.get(c) or {}).get("wallSeconds")
                                        for c in runs}
        cpu = d.get("processCpuSeconds") or {}
        names = [c for c in ("V", "F", "W", "K") if its.get(c) is not None]
        if not names:
            notes.append(f"cycle comparison: {tname} has no iteration counts")
            continue
        npan = 3 if any(cpu.get(c) for c in names) else 2
        fig, ax = plt.subplots(1, npan, figsize=(6.5, 2.4))
        ax[0].bar(names, [its[c] or 0 for c in names], color=C_COUPLED)
        ax[0].set_ylabel("Krylov iterations")
        ax[1].bar(names, [wall.get(c) or 0 for c in names], color=C_COUPLED2)
        ax[1].set_ylabel("solve wall time [s]")
        if npan == 3:
            ax[2].bar(names, [(cpu.get(c) or 0) / 3600.0 for c in names],
                      color=C_NATIVE)
            ax[2].set_ylabel("process CPU [CPU-h]")
        else:
            notes.append(f"{tname}: no processCpuSeconds (older record); "
                         "CPU panel omitted")
        stem = "cycle_comparison" if tname.endswith("cycles") \
            else "cycle_comparison_motorBike"
        save(fig, stem,
             f"Unit test: block-GAMG cycle types on the block-Poisson "
             f"system ({mesh}, serial).")
        for c in names:
            num(f"cycle {stem} {c} iters", its[c], "{}")


def fig_bench_cycles(rows: list[dict]) -> None:
    """Benchmark C (K, autoTune) vs H (fixed K) vs E (fixed V): wall and
    CPU-hours to convergence on the B10 cycle cases."""
    sel = [d for d in rows if d["config"] in ("C", "E", "H")]
    cases = sorted({d["case"] for d in sel if d["config"] == "E"})
    if not cases:
        notes.append("benchmark cycle comparison: no configuration E results")
        return
    cfgs = [c for c in ("C", "H", "E") if any(d["config"] == c for d in sel)]
    fig, ax = plt.subplots(1, 2, figsize=(6.5, 2.6))
    _bars(ax[0], sel, cases, cfgs, "wall_to_conv_s",
          "wall time to convergence [s]", log=False)
    _bars(ax[1], sel, cases, cfgs, "cpu_to_conv_h",
          "CPU time to convergence [CPU-h]", log=False)
    ax[0].legend(fontsize=7)
    save(fig, "bench_cycles",
         "Benchmark: default K-cycle with autoTune (C), fixed K-cycle (H) "
         "and fixed V-cycle (E): wall time and CPU-hours to convergence.")


def fig_anderson(rows: list[dict]) -> None:
    sel = [d for d in rows if d.get("config") in ("C", "G")]
    cases = sorted({d["case"] for d in sel if d["config"] == "G"})
    if not cases:
        notes.append("Anderson comparison: no configuration G results")
        return
    fig, ax = plt.subplots(1, 3, figsize=(6.5, 2.4))
    x = np.arange(len(cases))
    for j, (cfg, col, lab) in enumerate((("C", C_COUPLED, "Anderson off"),
                                         ("G", "#d95f02", "Anderson on"))):
        for k, key in enumerate(("iters_to_conv", "wall_to_conv_s",
                                 "cpu_to_conv_h")):
            ax[k].bar(x + 0.4 * j, [_med(sel, c, cfg, key) or 0 for c in cases],
                      0.4, color=col, label=lab)
    for a in ax:
        a.set_xticks(x + 0.2)
        a.set_xticklabels(cases)
    ax[0].set_ylabel("iterations to convergence")
    ax[1].set_ylabel("wall time to conv. [s]")
    ax[2].set_ylabel("CPU time to conv. [CPU-h]")
    ax[0].legend(fontsize=7)
    save(fig, "anderson", "Anderson acceleration on and off (15.7).")


def fig_scaling(tests: dict) -> None:
    """T-scaling (tests/test_scaling.py -> results/tests/T_scaling_T4b.json):
    per solver timePerIter_s{ranks}, cpuHours{ranks} (whole run of the
    fixed iteration count), efficiency{ranks}."""
    d = next((v for k, v in tests.items()
              if k.startswith(("T_scaling", "T-scaling"))), None)
    if not d or not any(isinstance(d.get(s), dict) for s in
                        ("coupledFoam", "simpleFoam")):
        notes.append("strong scaling: no T-scaling result")
        return
    fig, ax = plt.subplots(1, 3, figsize=(6.5, 2.4))
    for solver, col in (("simpleFoam", C_NATIVE), ("coupledFoam", C_COUPLED)):
        s = d.get(solver) or {}
        tpi = {int(k): v for k, v in (s.get("timePerIter_s") or {}).items()
               if v}
        if not tpi:
            continue
        n = np.array(sorted(tpi), dtype=float)
        t = np.array([tpi[int(k)] for k in n])
        ax[0].loglog(n, t, "o-", color=col, label=solver)
        eff = t[0] * n[0] / (t * n)
        ax[1].plot(n, eff, "o-", color=col, label=solver)
        ch = {int(k): v for k, v in (s.get("cpuHours") or {}).items() if v}
        if ch:
            m = sorted(ch)
            ax[2].plot(m, [ch[k] for k in m], "o-", color=col, label=solver)
    ax[0].set_xlabel("ranks")
    ax[0].set_ylabel("wall time per iteration [s]")
    ax[1].set_xlabel("ranks")
    ax[1].set_ylabel("parallel efficiency")
    ax[2].set_xlabel("ranks")
    ax[2].set_ylabel(f"CPU-h per run ({d.get('iterations', '?')} it.)")
    ax[0].legend(fontsize=7)
    save(fig, "scaling", f"Strong scaling (T-scaling on {d.get('mesh', 'T4b')}).")
    num("scaling ranks max", d.get("efficiencyRanks"), "{}")
    num("scaling rel efficiency", d.get("relativeEfficiency"), "{:.2f}")


def fig_remediation(tests: dict) -> None:
    cand = [(n, d) for n, d in tests.items()
            if (d.get("history") or {}).get("nDyn")
            and n.startswith(("T4", "T5", "T-fpe", "T1", "T3"))]
    if not cand:
        notes.append("remediation history: no data")
        return
    fig, ax = plt.subplots(figsize=(6.5, 2.4))
    for name, d in sorted(cand)[:5]:
        nd = d["history"]["nDyn"]
        ax.plot(np.arange(1, len(nd) + 1), nd, lw=0.9, label=name)
    ax.set_xlabel("iteration")
    ax.set_ylabel("cells in the dynamic set")
    ax.legend(fontsize=6)
    save(fig, "remediation_history", "Dynamic remediation set size.")


def table_validation(tests: dict) -> str:
    header = ["case", "quantity", "coupledFoam", "simpleFoam", "rel. diff.",
              "tolerance", "field RMS dU / dp", "pass"]
    rows = []
    specs = [
        ("T1_np1", "dp", "dp", "dpRef", "dpRelDiff", "tolDp"),
        ("T2_np1", "x_r/h", "xr_over_h", "xrRef_over_h", "xrRelDiff", "tolXr"),
        ("T3_kOmegaSST_np1", "Cd", "Cd", "CdRef", "CdRelDiff", "tol"),
        ("T3_kOmegaSST_np1", "Cl", "Cl", "ClRef", "ClRelDiff", "tol"),
        ("T3_GEKO_np1", "Cd", "Cd", "CdRef", "CdRelDiff", "tol"),
        ("T3_GEKO_np1", "Cl", "Cl", "ClRef", "ClRelDiff", "tol"),
    ]
    for name, q, kc, kr, kd, kt in specs:
        d = tests.get(name)
        if not d:
            continue
        rows.append([name, q, fmt(d.get(kc), "{:.5g}"), fmt(d.get(kr), "{:.5g}"),
                     fmt(d.get(kd), "{:.2e}"), fmt(d.get(kt), "{:.3g}"), "",
                     "yes" if d.get("pass") else "NO"])
    osc_rows = False
    for name in sorted(tests):
        if not name.startswith(("T4", "T5")) or "Cd" not in tests[name]:
            continue
        d = tests[name]
        ok = "yes" if d.get("pass") else "NO"
        if str(d.get("criterion", "")).startswith("stationaryMean"):
            # wake case (D-042 and addendum): window mean +- std, tolerance
            # max(relative, absolute) of the reference mean; mean-field
            # delta (volume RMS of |dUMean|/U_inf and |dpMean|/p_ref)
            osc_rows = True
            ref = d.get("reference") or {}
            for q, (rel, ab) in run_bench.OSC_TOL.items():
                rows.append([
                    f"{name}†", f"{q} (mean ± std)",
                    _pm(d.get(f"{q}_mean"), d.get(f"{q}_std")),
                    _pm(ref.get(f"{q}_mean"), ref.get(f"{q}_std")),
                    fmt(d.get(f"{q}RelDiff"), "{:.2e}"),
                    f"max({rel:.0%}, {ab:g})",
                    _field_rms(d) if q == "Cd" else "", ok])
        else:
            rows.append([name, "Cd", fmt(d.get("Cd"), "{:.5g}"),
                         fmt(d.get("CdRef"), "{:.5g}"),
                         fmt(d.get("CdRelDiff"), "{:.2e}"), "0.01", "", ok])
    if not rows:
        notes.append("validation table: no T1-T5 results")
        return ""
    note = None
    if osc_rows:
        ft = run_bench.FIELD_TOL
        note = ("† Oscillating wake: stationary window mean over the "
                f"last W = max({run_bench.STAT_WINDOW_MIN}, "
                f"n/{run_bench.STAT_WINDOW_DIV}) iterations (capped at n), "
                "mean ± standard deviation; tolerance max(relative, "
                "absolute) of the simpleFoam mean (averaged force "
                "criterion, D-042 addendum). Field RMS dU / dp: "
                "volume-weighted RMS of the mean-velocity difference "
                "(magnitude, over U_inf) and of the mean-pressure "
                "difference (over the free-stream dynamic pressure) "
                "between the solvers, fields averaged over the same "
                "window (coupledFieldCompare); proposed limits "
                f"{ft['volRmsMagUDeltaRel']:.0%} / "
                f"{ft['volRmsPDeltaRel']:.0%}; n/a: not evaluated.")
    return write_table("validation", header, rows,
                       "Validation: integral quantities, coupledFoam vs. "
                       "simpleFoam on identical meshes and schemes; wake "
                       "cases T4 and T5 compared by window means and "
                       "mean-field deltas (D-042 addendum).",
                       "tab:validation", resize=True, note=note)


def _field_rms(d: dict) -> str:
    """'x % / y %' of the mean-field delta (D-042 addendum), n/a if the
    record predates it or the comparison failed."""
    u, p = d.get("fieldRmsU"), d.get("fieldRmsP")
    if u is None or p is None:
        return "n/a"
    return f"{100 * u:.2f}% / {100 * p:.2f}%"


def _pm(m, s) -> str:
    """mean +- std of a window (D-042 rows of the validation table)."""
    if m is None:
        return "n/a"
    return f"{m:.4f} ± {s:.4f}" if s is not None else f"{m:.4f}"


def _gamg_record(name: str, d: dict) -> dict:
    """GAMG statistics of a record, completed from the run's log (ranks per
    level, ratios, coarsest solver, cycle after autoTune)."""
    g = dict(d.get("gamg") or {})
    log = RUN / name / "log.coupledFoam"
    if log.exists() and run_ok(RUN / name):
        for k, v in logs.gamg_log_stats(log).items():
            if g.get(k) is None:
                g[k] = v
    cells = g.get("gamgCellsPerLevel") or []
    if not g.get("gamgRatios") and len(cells) > 1:
        g["gamgRatios"] = [cells[i] / cells[i + 1] for i in range(len(cells) - 1)
                           if cells[i + 1]]
    return g


def table_gamg_levels(tests: dict, bench: list[dict]) -> str:
    header = ["run", "cycle (start/end)", "post-sweeps", "levels",
              "mergeLevels", "C_op", "cells per level", "ranks per level",
              "coarsening ratios", "coarsest solver"]
    rows = []
    srcs = [(n, d) for n, d in sorted(tests.items())]
    seen = set()
    for d in sorted(bench, key=lambda r: (r["case"], r["config"], r["run"])):
        key = (d["case"], d["config"])
        if d["config"] == "C" and key not in seen:
            seen.add(key)
            srcs.append((f"bench_{d['case']}_C_{d['run']}", d))
    for name, d in srcs:
        g = _gamg_record(name, d)
        if not g.get("gamgLevels"):
            continue
        cyc = g.get("gamgCycleInitial")
        cyc_end = g.get("gamgCycleFinal")
        rows.append([
            name.replace("bench_", ""),
            f"{cyc or 'n/a'}/{cyc_end or 'n/a'}",
            f"{fmt(g.get('gamgNPostSweepsInitial'), '{}')}/"
            f"{fmt(g.get('gamgNPostSweepsFinal'), '{}')}",
            g.get("gamgLevels"), fmt(g.get("gamgMergeLevels"), "{}"),
            fmt(g.get("gamgCop"), "{:.3f}"),
            flat(g.get("gamgCellsPerLevel")),
            flat(g.get("gamgRanksPerLevel")),
            flat([round(r, 2) for r in g.get("gamgRatios") or []]),
            g.get("gamgCoarsestSolver") or "n/a"])
    if not rows:
        notes.append("GAMG level table: no data")
        return ""
    return write_table("gamg_levels", header, rows,
                       "Block-GAMG hierarchy per run: cycle type and "
                       "post-smoothing sweeps at start and end (autoTune), "
                       "levels, operator complexity, cells, ranks and "
                       "coarsening ratios per level, coarsest-level solver. "
                       "Ranks per level are n/a for logs written before "
                       "D-030.", "tab:gamglevels", resize=True)


def table_b10(rows: list[dict]) -> str:
    """Amendment B10 acceptance: C vs E/F/G/H in wall and CPU-hours."""
    b10 = run_bench.b10_evaluate(run_bench.summary_rows(rows)) if rows else []
    if not b10:
        notes.append("B10 acceptance: no E/F/G/H benchmark results")
        return ""
    header = ["case", "config", "wall X [s]", "wall C [s]", "d wall",
              "CPU-h X", "CPU-h C", "d CPU-h", "monitor rel. diff.", "status"]
    tab = []
    for r in b10:
        tab.append([r["case"], r["config"], fmt(r.get("wall_X")),
                    fmt(r.get("wall_C")), pct(r.get("dWall_X_vs_C")),
                    fmt(r.get("cpuh_X")), fmt(r.get("cpuh_C")),
                    pct(r.get("dCpu_X_vs_C")),
                    fmt(r.get("monitorRelDiff"), "{:.1e}"), r.get("status")])
        k = f"b10 {r['case']} {r['config']}"
        num(f"{k} dwall", pct(r.get("dWall_X_vs_C")).replace("%", r"\%"))
        num(f"{k} dcpu", pct(r.get("dCpu_X_vs_C")).replace("%", r"\%"))
        num(f"{k} status", r.get("status"))
    return write_table("b10_acceptance", header, tab,
                       "Amendment B10: configurations E (fixed V-cycle), "
                       "F (fixed relTol), G (Anderson) and H (fixed K-cycle) "
                       "against C; d = X/C$-$1 (positive: X slower). Pass "
                       "criterion only for F: C at most 5\\,\\% slower than F "
                       "in wall time and CPU-hours and $C_d$ or $\\Delta p$ "
                       "identical to $10^{-4}$.", "tab:btenacc", resize=True)


# --------------------------------------------------------------------------- #
# speed-up figures from the test records and the run logs
# --------------------------------------------------------------------------- #

# (label, test record, coupledFoam run dir, simpleFoam reference dir)
SPEED_CASES = [
    ("T0 Re100", "T0_Re100_np1", "T0_Re100_np1", "ref_T0_Re100"),
    ("T0 Re1000", "T0_Re1000_np1", "T0_Re1000_np1", "ref_T0_Re1000"),
    ("T1", "T1_np1", "T1_np1", "ref_T1"),
    ("T2", "T2_np1", "T2_np1", "ref_T2"),
    ("T3 SST", "T3_kOmegaSST_np1", "T3_kOmegaSST_np1", "ref_T3_kOmegaSST"),
    ("T3 GEKO", "T3_GEKO_np1", "T3_GEKO_np1", "ref_T3_GEKO"),
    ("T4a", "T4a_np10", "T4a_np10", "ref_T4a_np10"),
    ("T4b", "T4b_np10", "T4b_np10", "ref_T4b_np10"),
    ("T5", "T5_np10", "T5_np10", "ref_T5_np10"),
]
BUSY_MINUTES = 10.0   # a log written less than this ago belongs to a live run
_RES_LINE = logs._RES
_EXEC_LINE = logs._EXEC


def _busy(case: Path) -> bool:
    import time  # noqa: PLC0415
    return any(time.time() - f.stat().st_mtime < BUSY_MINUTES * 60
               for f in case.glob("log.*"))


@functools.lru_cache(maxsize=None)
def cf_timeline(case: Path) -> dict | None:
    """coupledFoam per-iteration R, tAsm/tSolve/tTurb and the wall time axis
    (cumulative tIter, scaled so that its end equals the solver wall time of
    postProcessing/coupledFoam/summary.json when it exists)."""
    log = case / "log.coupledFoam"
    if not log.exists() or _busy(case) or not run_ok(case):
        return None
    rows = [r for r in logs.parse_cf(log) if "R" in r and "tIter" in r]
    if not rows:
        return None
    t = np.cumsum([r["tIter"] for r in rows])
    summ = logs.coupled_summary(case)
    if not summ:
        return None      # unfinished or aborted run: no summary.json
    wall = summ.get("wallSeconds")
    if wall and t[-1] > 0:
        t = t * (wall / t[-1])
    return {
        "iter": np.array([r["iter"] for r in rows]),
        "R": np.array([r["R"] for r in rows]),
        "t": t,
        "tAsm": np.array([r.get("tAsm", np.nan) for r in rows]),
        "tSolve": np.array([r.get("tSolve", np.nan) for r in rows]),
        "tTurb": np.array([r.get("tTurb", np.nan) for r in rows]),
        "tIter": np.array([r["tIter"] for r in rows]),
        "summary": summ,
    }


@functools.lru_cache(maxsize=None)
def sf_timeline(case: Path) -> dict | None:
    """simpleFoam per-iteration initial residuals (first solve of each field)
    and the wall time axis: ExecutionTime per iteration, scaled so that its
    end equals the final ClockTime (ClockTime itself is in whole seconds)."""
    log = case / "log.simpleFoam"
    if not log.exists() or _busy(case):
        return None
    res: dict[str, list[float]] = {}
    ex, clk = [], []
    seen: set[str] = set()
    with open(log, errors="replace") as fh:
        for line in fh:
            if line.startswith("Time = "):
                seen = set()
                continue
            m = _RES_LINE.match(line)
            if m:
                f = m.group(2)
                if f not in seen:
                    res.setdefault(f, []).append(float(m.group(3)))
                    seen.add(f)
                continue
            m = _EXEC_LINE.search(line)
            if m:
                ex.append(float(m.group(1)))
                clk.append(float(m.group(2)))
    if not ex:
        return None
    t = np.array(ex)
    if clk[-1] > 0 and t[-1] > 0:
        t = t * (clk[-1] / t[-1])
    n = len(t)
    return {"t": t, "res": {k: np.array(v[:n]) for k, v in res.items()},
            "n": n}


# Residual target of the speed-up on the residual cases when the record does
# not carry one (T2: tests/test_T2 R_TARGET)
SPEED_R_DEFAULT = 1e-5


def speed_criterion(rec_name: str, d: dict) -> tuple[str, float | None]:
    """Criterion both solvers are timed to in the speed-up (M1, D-068):
    ("residual", R) on T0-T2 - the test's R target of coupledFoam's combined
    residual, applied to simpleFoam as "every initial residual of the log
    below R" (the D-024 definition of T2); ("forceWindow", None) on T3 -
    spec 12.3(ii), the 100-iteration Cd/Cl window; ("stationary", None) on
    the wake cases - the D-042 point under the case's common window
    (D-068)."""
    if str(d.get("criterion", "")).startswith("stationaryMean") or \
            run_bench.is_oscillatory(run_bench.case_of_run(rec_name) or ""):
        return "stationary", None
    if rec_name.startswith("T3"):
        return "forceWindow", None
    return "residual", float(d.get("Rtarget") or SPEED_R_DEFAULT)


def _sf_residual_iteration(sft: dict | None, target: float) -> int | None:
    """First simpleFoam iteration at which EVERY initial residual of the log
    (first solve of each field per iteration: p, Ux, Uy, [Uz], k, omega, ...)
    is below `target`."""
    if not sft or not sft["res"]:
        return None
    hist = list(sft["res"].values())
    n = min(len(h) for h in hist)
    for i in range(n):
        if all(h[i] < target for h in hist):
            return i + 1
    return None


def _frac(tl: dict | None, it: int | None, n: int | None) -> float | None:
    """Wall-clock fraction of a run spent up to iteration `it` (timeline of
    cf_timeline / sf_timeline; iteration share without one)."""
    if it is None:
        return None
    if tl is not None and len(tl["t"]) and it <= len(tl["t"]) \
            and tl["t"][-1] > 0:
        return float(tl["t"][it - 1] / tl["t"][-1])
    return it / n if n else None


def _speed_record(tests: dict, rec: str, cfd: str, sfd: str) -> dict | None:
    """Wall time and CPU-hours to convergence of both solvers for one case,
    both timed to the SAME criterion (speed_criterion; review M1, D-068):

    residual (T0-T2)  coupledFoam: first R < R_target (iterationsToR);
                      simpleFoam: first iteration with every initial
                      residual < R_target (its log). Times: solver only
                      (coupledFoam solver loop, simpleFoam ClockTime), the
                      wall-clock fraction of the run up to that iteration.
    forceWindow (T3)  12.3(ii) on both force histories (a D-060 user point
                      from the record takes precedence); times as above.
    stationary (T4/T5) the D-042 point under the case's common window
                      (D-068) from the test record (wallToConv_s,
                      cpuHoursToConv: rank timing incl. potentialFoam), or
                      recomputed from the run directories for records
                      written before D-068; the per-run-window values are
                      added as *_perRun. Flags of review M2
                      (run_bench.reference_timing_flags).

    A solver that never reaches the criterion is shown with its whole run
    (conv_* False; for simpleFoam the speed-up is then a lower bound).
    it_*_conv are the convergence iterations (None if not reached)."""
    d = tests.get(rec)
    if not d:
        return None
    # M9: the run directory read below must still be this record's run
    if (RUN / cfd).is_dir() and not GUARD.check_record_run(rec, d, RUN / cfd):
        return None
    ref = d.get("reference") or {}
    kind, target = speed_criterion(rec, d)
    out: dict = {"criterion": kind, "Rtarget": target,
                 "np_cf": d.get("nProcs") or 1,
                 "np_sf": ref.get("nProcs") or 1}
    if kind == "stationary":
        return _speed_record_stationary(d, ref, cfd, sfd, out)

    # ---- coupledFoam
    n_cf = d.get("iterations")
    wall_cf_run = d.get("wallSecondsSolver") or d.get("wallSeconds")
    cpu_cf_run = d.get("cpuHoursSolver") or d.get("cpuHours")
    if not (n_cf and wall_cf_run and cpu_cf_run):
        return None
    cft = cf_timeline(RUN / cfd)
    if kind == "residual":
        it_cf = d.get("iterationsToR") or d.get("iterationsToR_coupled")
        if it_cf is None and cft is not None:
            hit = np.nonzero(cft["R"] < target)[0]
            it_cf = int(hit[0]) + 1 if hit.size else None
    else:
        it_cf = d.get("itersToConv")
        if it_cf is None and run_ok(RUN / cfd):
            it_cf = run_bench.iters_to_conv(run_bench.force_history(RUN / cfd))
        if it_cf is None and cft is not None and \
                cft["summary"].get("converged"):
            it_cf = n_cf        # the solver's own stop (R below its tolerance)
    f_cf = _frac(cft, it_cf, n_cf)

    # ---- simpleFoam reference
    sft = sf_timeline(RUN / sfd)
    n_sf = ref.get("iterations") or (sft["n"] if sft else None)
    ta = ref.get("timingAllrun") or {}
    wall_sf_run = (float(sft["t"][-1]) if sft else None) \
        or ref.get("wallSecondsSolver") or ref.get("wallSeconds") \
        or ta.get("wallSeconds")
    cpu_sf_run = ref.get("cpuHours") or ta.get("cpuHours")
    if cpu_sf_run is None and wall_sf_run:
        cpu_sf_run = wall_sf_run * out["np_sf"] / 3600.0  # serial: CPU = wall
    if not wall_sf_run:
        return None
    if kind == "residual":
        it_sf = _sf_residual_iteration(sft, target)
    else:
        u = ref.get("user") or {}
        if u and not u.get("ignored"):
            it_sf = ref.get("itersToConv")       # D-060 user point
        elif (RUN / sfd).is_dir():
            it_sf = run_bench.iters_to_conv(run_bench.force_history(RUN / sfd))
        else:
            it_sf = ref.get("itersToConv")
    f_sf = _frac(sft, it_sf, n_sf)
    conv_cf, conv_sf = f_cf is not None, f_sf is not None
    out.update({
        "wall_cf": wall_cf_run * (f_cf if conv_cf else 1.0),
        "cpu_cf": cpu_cf_run * (f_cf if conv_cf else 1.0),
        "it_cf": it_cf if conv_cf else n_cf, "conv_cf": conv_cf,
        "it_cf_conv": it_cf if conv_cf else None, "n_cf": n_cf,
        "wall_sf": wall_sf_run * (f_sf if conv_sf else 1.0),
        "cpu_sf": cpu_sf_run * (f_sf if conv_sf else 1.0),
        "it_sf": it_sf if conv_sf else n_sf, "conv_sf": conv_sf,
        "it_sf_conv": it_sf if conv_sf else None, "n_sf": n_sf,
    })
    return out


def _speed_record_stationary(d: dict, ref: dict, cfd: str, sfd: str,
                             out: dict) -> dict | None:
    """Wake cases (T4, T5): see _speed_record."""
    case = run_bench.case_of_run(cfd)
    user = d.get("convergenceSource") == "user"
    common = str(d.get("windowRule", "")).startswith("common")
    pts = {}
    for side, r, rdir, solver in (("cf", d, RUN / cfd, "coupledFoam"),
                                  ("sf", ref, RUN / sfd, "simpleFoam")):
        p: dict = {}
        if (user or common) and r.get("wallToConv_s") is not None:
            p = {"iters": r.get("itersToConv"), "wall": r.get("wallToConv_s"),
                 "cpuh": r.get("cpuHoursToConv"),
                 "iters_perRun": r.get("iters_to_stationary_perRun"),
                 "wall_perRun": r.get("wallToConv_s_perRun"),
                 "cpuh_perRun": r.get("cpuHoursToConv_perRun"),
                 "n": r.get("iterationsRun") or r.get("iterations"),
                 "wallTotal": r.get("wallSeconds"),
                 "cpuhTotal": r.get("cpuHours")}
        elif rdir.is_dir() and (side == "sf" or run_ok(rdir)):
            # record from before D-068: recompute from the run directory
            p = run_bench.stationary_point(rdir, solver, case)
        if not p or not p.get("wallTotal"):
            return None
        pts[side] = p
    for s, p in pts.items():
        conv = p.get("wall") is not None and p.get("iters") is not None
        out.update({
            f"wall_{s}": p["wall"] if conv else p["wallTotal"],
            f"cpu_{s}": p["cpuh"] if conv else p["cpuhTotal"],
            f"it_{s}": p["iters"] if conv else p.get("n"),
            f"conv_{s}": conv, f"it_{s}_conv": p["iters"] if conv else None,
            f"n_{s}": p.get("n"),
            f"wall_{s}_perRun": p.get("wall_perRun"),
            f"cpu_{s}_perRun": p.get("cpuh_perRun"),
            f"it_{s}_perRun": p.get("iters_perRun"),
        })
    for k in ("wall", "cpu"):
        a, b = out.get(f"{k}_sf_perRun"), out.get(f"{k}_cf_perRun")
        out[f"speedup_{k}_perRun"] = a / b if a and b else None
    out["W"] = case and run_bench.case_window(case)
    out.update(run_bench.reference_timing_flags(RUN / sfd, ref))
    return out


def fig_speed(tests: dict) -> None:
    """Speed-up section: residual vs wall time, time and CPU-hours to
    convergence with speed-up factors, time per iteration breakdown and the
    iterations/cost trade-off."""
    # ---- (a) residual against wall-clock time
    avail = []
    for lab, rec, cfd, sfd in SPEED_CASES:
        cft, sft = cf_timeline(RUN / cfd), sf_timeline(RUN / sfd)
        if cft is not None or sft is not None:
            avail.append((lab, rec, cft, sft))
        elif (RUN / cfd).is_dir() and _busy(RUN / cfd):
            notes.append(f"speed-up: {lab} is running, skipped")
    if avail:
        ncol = 3
        nrow = int(np.ceil(len(avail) / ncol))
        fig, axs = plt.subplots(nrow, ncol, figsize=(6.5, 2.1 * nrow),
                                squeeze=False)
        for ax, (lab, rec, cft, sft) in zip(axs.flat, avail):
            if sft is not None:
                for f, ls in (("p", "-"), ("Ux", "--")):
                    r = sft["res"].get(f)
                    if r is not None and len(r):
                        ax.semilogy(sft["t"][:len(r)], r, color=C_NATIVE,
                                    ls=ls, lw=1.2, zorder=3,
                                    label=f"simpleFoam {f}")
            if cft is not None:
                ax.semilogy(cft["t"], cft["R"], color=C_COUPLED, lw=1.4,
                            zorder=4, label="coupledFoam $R$")
            else:
                ax.text(0.5, 0.5, "coupledFoam\npending", ha="center",
                        va="center", transform=ax.transAxes, fontsize=7,
                        color=C_COUPLED)
            ax.set_xscale("symlog", linthresh=1.0)
            ax.set_title(lab, fontsize=8)
            ax.set_xlabel("wall-clock time [s]", fontsize=7)
            ax.tick_params(labelsize=6)
        for ax in list(axs.flat)[len(avail):]:
            ax.axis("off")
        for r in range(nrow):
            axs[r, 0].set_ylabel("residual", fontsize=7)
        axs.flat[0].legend(fontsize=6, loc="lower left")
        save(fig, "speed_residual_wall",
             "Residual against wall-clock time on the same axes: "
             "coupledFoam combined residual $R$ (blue) and simpleFoam "
             "initial residuals of $p$ and $U_x$ (orange). The two residual "
             "normalisations differ (Section 3.2); the time axis is the "
             "comparable quantity.")

    # ---- (b) wall time and CPU-hours to convergence with speed-ups
    recs = [(lab, _speed_record(tests, rec, cfd, sfd))
            for lab, rec, cfd, sfd in SPEED_CASES]
    recs = [(lab, r) for lab, r in recs if r]
    if recs:
        fig, axs = plt.subplots(1, 2, figsize=(6.5, 3.0))
        x = np.arange(len(recs))
        w = 0.38
        for ax, key, ylab in ((axs[0], "wall", "wall-clock time [s]"),
                              (axs[1], "cpu", "CPU time [CPU-h]")):
            for j, (s, col, name) in enumerate((("sf", C_NATIVE, "simpleFoam"),
                                                ("cf", C_COUPLED, "coupledFoam"))):
                vals = [r[f"{key}_{s}"] for _, r in recs]
                hatch = ["" if r[f"conv_{s}"] else "////" for _, r in recs]
                bars = ax.bar(x + (j - 0.5) * w, vals, w, color=col,
                              label=name, edgecolor="white", linewidth=0.5)
                for b, h in zip(bars, hatch):
                    b.set_hatch(h)
            for i, (_, r) in enumerate(recs):
                a, b = r[f"{key}_sf"], r[f"{key}_cf"]
                if not r["conv_cf"]:
                    txt = "not conv."
                else:
                    txt = f"{'≥' if not r['conv_sf'] else ''}{a / b:.1f}×"
                ax.text(i, max(a, b) * 1.25, txt, ha="center", va="bottom",
                        fontsize=6)
            ax.set_yscale("log")
            ax.set_xticks(x)
            ax.set_xticklabels([lab for lab, _ in recs], fontsize=6,
                               rotation=30)
            ax.set_ylabel(ylab, fontsize=7)
            ax.tick_params(labelsize=6)
            lo, hi = ax.get_ylim()
            ax.set_ylim(lo, hi * 3)
        axs[0].legend(fontsize=6, loc="upper left")
        save(fig, "speed_time_to_conv",
             "Wall-clock time (left) and CPU-hours (right) to convergence, "
             "simpleFoam (orange) and coupledFoam (blue), serial runs; "
             "numbers: speed-up simpleFoam/coupledFoam. Hatched: the run "
             "did not reach its criterion, the bar is the whole run "
             "(for simpleFoam a lower bound, so the speed-up is marked "
             "$\\geq$).")
        for lab, r in recs:
            k = lab.replace(" ", "")
            if r["conv_cf"]:
                num(f"speed wall {k}", r["wall_sf"] / r["wall_cf"], "{:.1f}")
                num(f"speed cpu {k}", r["cpu_sf"] / r["cpu_cf"], "{:.1f}")
            num(f"speed wall cf {k}", r["wall_cf"], "{:.4g}")
            num(f"speed wall sf {k}", r["wall_sf"], "{:.4g}")
            num(f"speed cpuh cf {k}", r["cpu_cf"], "{:.2g}")
            num(f"speed cpuh sf {k}", r["cpu_sf"], "{:.2g}")
            num(f"speed iters cf {k}", r["it_cf"], "{}")
            num(f"speed iters sf {k}", r["it_sf"], "{}")

    # ---- (c) time per outer iteration: breakdown and cost vs simpleFoam
    br = []
    for lab, rec, cfd, sfd in SPEED_CASES:
        cft = cf_timeline(RUN / cfd)
        sft = sf_timeline(RUN / sfd)
        if cft is None:
            continue
        per_sf = (sft["t"][-1] / sft["n"]) if sft and sft["n"] else None
        br.append((lab, np.nanmean(cft["tAsm"]), np.nanmean(cft["tSolve"]),
                   np.nanmean(cft["tTurb"]), np.nanmean(cft["tIter"]), per_sf))
    if br:
        fig, axs = plt.subplots(1, 2, figsize=(6.5, 2.7))
        x = np.arange(len(br))
        a = np.array([b[1] for b in br])
        s = np.array([b[2] for b in br])
        t = np.array([b[3] for b in br])
        tot = np.array([b[4] for b in br])
        other = np.clip(tot - a - s - t, 0, None)
        axs[0].bar(x, a / tot, color="#9ecae1", label="assembly $t_{asm}$",
                   edgecolor="white", linewidth=0.5)
        axs[0].bar(x, s / tot, bottom=a / tot, color=C_COUPLED,
                   label="linear solve $t_{solve}$", edgecolor="white",
                   linewidth=0.5)
        axs[0].bar(x, t / tot, bottom=(a + s) / tot, color="#fdae6b",
                   label="turbulence $t_{turb}$", edgecolor="white",
                   linewidth=0.5)
        axs[0].bar(x, other / tot, bottom=(a + s + t) / tot, color="#d9d9d9",
                   label="other", edgecolor="white", linewidth=0.5)
        axs[0].set_ylabel("share of the iteration time", fontsize=7)
        axs[0].set_ylim(0, 1.0)
        axs[0].legend(fontsize=6, loc="upper center",
                      bbox_to_anchor=(0.5, -0.25), ncol=2)
        w = 0.38
        per_sf = np.array([b[5] if b[5] else np.nan for b in br])
        axs[1].bar(x - w / 2, per_sf, w, color=C_NATIVE, label="simpleFoam")
        axs[1].bar(x + w / 2, tot, w, color=C_COUPLED, label="coupledFoam")
        for i in range(len(br)):
            if np.isfinite(per_sf[i]) and per_sf[i] > 0:
                axs[1].text(i, max(per_sf[i], tot[i]) * 1.3,
                            f"{tot[i] / per_sf[i]:.1f}×", ha="center",
                            fontsize=6)
        axs[1].set_yscale("log")
        lo, hi = axs[1].get_ylim()
        axs[1].set_ylim(lo, hi * 3)
        axs[1].set_ylabel("wall time per outer iteration [s]", fontsize=7)
        axs[1].legend(fontsize=6, loc="upper left")
        for ax in axs:
            ax.set_xticks(x)
            ax.set_xticklabels([b[0] for b in br], fontsize=6, rotation=30)
            ax.tick_params(labelsize=6)
        save(fig, "speed_iteration_cost",
             "coupledFoam time per outer iteration: share of assembly, "
             "linear solve and turbulence (left, mean over the run); wall "
             "time per iteration of both solvers with the cost ratio "
             "coupledFoam/simpleFoam (right).")
        for b in br:
            k = b[0].replace(" ", "")
            num(f"speed iter cost ratio {k}",
                b[4] / b[5] if b[5] else None, "{:.1f}")
            num(f"speed solve share {k}", 100 * b[2] / b[4], "{:.0f}")

    # ---- (d) iterations vs wall time trade-off
    if recs:
        fig, ax = plt.subplots(figsize=(6.5, 3.2))
        xs = [r["it_sf"] for _, r in recs] + [r["it_cf"] for _, r in recs]
        ys = [r["wall_sf"] for _, r in recs] + [r["wall_cf"] for _, r in recs]
        lx = np.logspace(np.log10(min(xs) / 2), np.log10(max(xs) * 2), 10)
        ytop = max(ys) * 2
        for per in (1e-3, 1e-2, 1e-1, 1, 10):
            ax.plot(lx, per * lx, color="k", lw=0.4, alpha=0.25)
            # label where the diagonal leaves the plot (right or top edge)
            xm = min(lx[-1], ytop / per) / 1.15
            if lx[0] < xm and min(ys) / 2 < per * xm:
                ax.text(xm, per * xm, f"{per:g} s/it", fontsize=5,
                        alpha=0.6, ha="right", va="bottom")
        for lab, r in recs:
            ax.annotate("", xy=(r["it_cf"], r["wall_cf"]),
                        xytext=(r["it_sf"], r["wall_sf"]),
                        arrowprops=dict(arrowstyle="->", color="#9e9e9e",
                                        lw=0.7))
            ax.plot(r["it_sf"], r["wall_sf"], "o", color=C_NATIVE, ms=5,
                    mfc="white" if not r["conv_sf"] else C_NATIVE)
            ax.plot(r["it_cf"], r["wall_cf"], "s", color=C_COUPLED, ms=5,
                    mfc="white" if not r["conv_cf"] else C_COUPLED)
            ax.text(r["it_cf"] * 0.8, r["wall_cf"], lab, fontsize=6,
                    ha="right", va="center")
        ax.plot([], [], "o", color=C_NATIVE, label="simpleFoam")
        ax.plot([], [], "s", color=C_COUPLED, label="coupledFoam")
        ax.plot([], [], "o", color="k", mfc="white",
                label="criterion not reached (whole run)")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_ylim(min(ys) / 2, max(ys) * 2)
        ax.set_xlim(min(xs) / 2, max(xs) * 2)
        ax.set_xlabel("outer iterations (to convergence or run length)",
                      fontsize=7)
        ax.set_ylabel("wall-clock time [s]", fontsize=7)
        ax.tick_params(labelsize=6)
        ax.legend(fontsize=6, loc="upper left")
        save(fig, "speed_tradeoff",
             "Iterations against wall-clock time: each arrow goes from "
             "simpleFoam to coupledFoam for one case. Diagonal lines are "
             "constant time per iteration; coupledFoam moves far to the "
             "left (fewer iterations) and up to a more expensive diagonal. "
             "It is faster where the arrow points down.")


def fig_fields() -> None:
    """2D flow-field comparisons (bench/plot_fields2d.py, always) and the
    ParaView renders of the 3D cases (bench/render_fields.py, if pvbatch
    exists and CF_NO_RENDER is not set)."""
    import os  # noqa: PLC0415
    import shutil  # noqa: PLC0415
    try:
        import plot_fields2d  # noqa: PLC0415
    except ImportError as e:
        notes.append(f"fields 2D: {e}")
        return
    figs, stats, nts = plot_fields2d.run(log=lambda s: print(s, flush=True))
    notes.extend(nts)
    for stem, cap in figs:
        figures_pdf_only.append((stem, cap))
    for k, v in stats.items():
        if k.endswith(("_rms", "_max")):
            num(f"fld {k}", v, "{:.1e}")
    pv = shutil.which("pvbatch")
    if os.environ.get("CF_NO_RENDER"):
        notes.append("3D renders: CF_NO_RENDER set, skipped")
    elif not pv:
        notes.append("3D renders: pvbatch not found, skipped (figures keep "
                     "their previous state or show the placeholder)")
    else:
        cmd = ["nice", "-n", "19", pv, "--force-offscreen-rendering",
               str(REPO / "bench" / "render_fields.py"),
               "--out", str(FIG_PDF)]
        if GUARD.enabled:
            cmd += ["--commit", GUARD.target]
            if GUARD.allow_stale:
                cmd.append("--allow-stale")
        try:
            p = subprocess.run(cmd, cwd=REPO, capture_output=True, text=True,
                               timeout=3600)
            # pvbatch may route python output to stderr: scan both
            for line in (p.stdout + "\n" + p.stderr).splitlines():
                if line.startswith("wrote"):
                    print("render:", line)
            nf = FIG_PDF / "render_fields_notes.json"
            if nf.exists():
                for n in json.loads(nf.read_text()):
                    notes.append("3D renders: " + n)
        except (OSError, subprocess.TimeoutExpired) as e:
            notes.append(f"3D renders failed: {e}")
    if GUARD.enabled:
        prune_renders()


def fig_iteration_histories() -> None:
    """Per-case iteration histories and the appendix text
    (bench/plot_histories.py -> figures/hist_*.pdf, figures/histories.tex)."""
    try:
        import plot_histories  # noqa: PLC0415
    except ImportError as e:
        notes.append(f"histories: {e}")
        return
    figs, nts = plot_histories.run(log=lambda s: print(s, flush=True))
    notes.extend(nts)
    figures_pdf_only.extend(figs)


# figures written by the helper modules (PDF/PNG in report/paper/figures
# only, not in report/figures); listed in REPORT.md without an image
figures_pdf_only: list[tuple[str, str]] = []


# --------------------------------------------------------------------------- #
# exploratory numbers (bench/exploratory_numbers.py)
# --------------------------------------------------------------------------- #

def exploratory_numbers() -> None:
    ex = load_json("exploratory")
    if not ex:
        notes.append("exploratory numbers: run bench/exploratory_numbers.py")
        return
    # staleness guard: these records describe historical studies of
    # run/exp_* and carry no gitCommit; used only with --allow-historical
    for k in EXPLORATORY_USED:
        d = ex.get(k)
        if d is None:
            GUARD.missing(f"exploratory {k}", f"results/exploratory/{k}.json")
            continue
        if GUARD.matches(d.get("gitCommit")):
            continue
        if ALLOW_HISTORICAL:
            GUARD.historical(f"exploratory {k}",
                             f"results/exploratory/{k}.json",
                             d.get("gitCommit"), d.get("timestamp"),
                             "historical study (run/exp_*), accepted by "
                             "--allow-historical")
        elif not GUARD.check_record(f"exploratory {k}", d,
                                    f"results/exploratory/{k}.json"):
            ex.pop(k)
    re1 = (ex.get("re1000") or {}).get("runs") or {}
    for run, keys in (
        ("default", (("cflCutsTotal", "{}"), ("iterations", "{}"),
                     ("finalR", "{:.1e}"), ("minR", "{:.1e}"),
                     ("cflP05AfterStartup", "{:.2g}"),
                     ("cflP95AfterStartup", "{:.2g}"),
                     ("cflMinAfterStartup", "{:.2g}"),
                     ("cflMaxAfterStartup", "{:.2g}"),
                     ("linItersTotal", "{}"), ("tSolveTotal_s", "{:.0f}"))),
        ("upwind", (("firstIterBelowTarget", "{}"), ("cflCutsTotal", "{}"),
                    ("firstIterAtCFLmax", "{}"))),
        ("ilu0", (("firstIterBelowTarget", "{}"), ("cflCutsTotal", "{}"),
                  ("firstIterAtCFLmax", "{}"), ("tSolveTotal_s", "{:.0f}"))),
    ):
        r = re1.get(run)
        if r is None:
            notes.append(f"exploratory re1000/{run}: missing")
            continue
        for k, f in keys:
            num(f"exp re1000 {run} {k}", r.get(k), f)
        if run == "default":
            s = r.get("settings") or {}
            num("exp re1000 default relTol", s.get("relTol"), "{:g}")
    lin = (ex.get("linsolver") or {}).get("runs") or {}
    gs = lin.get("blockBiCGStab_blockGaussSeidel")
    il = lin.get("blockBiCGStab_blockILU0")
    gm = lin.get("blockGMRES_blockGaussSeidel")
    for tag, r in (("bicgGs", gs), ("bicgIlu", il), ("gmresGs", gm)):
        if r is None:
            notes.append(f"exploratory linsolver/{tag}: missing")
            continue
        num(f"exp lin {tag} iterations", r.get("iterations"), "{}")
        num(f"exp lin {tag} linIters", r.get("linItersTotal"), "{}")
        num(f"exp lin {tag} maxed", r.get("linSolvesAtMaxIter"), "{}")
        num(f"exp lin {tag} tSolve", r.get("tSolveTotal_s"), "{:.0f}")
    if gs and il and il.get("tSolveTotal_s"):
        num("exp lin solve ratio gs ilu",
            gs["tSolveTotal_s"] / il["tSolveTotal_s"], "{:.1f}")
    sym = ex.get("symbol") or {}
    if sym:
        num("exp sym pe low", sym.get("peLow"), "{:.2g}")
        num("exp sym pe high", sym.get("peHigh"), "{:.2g}")
        for k in ("rhoInfPeLow", "rhoInfPeHigh", "rhoCflFiveHundredPeLow"):
            num(f"exp sym {k}", (sym.get(k) or {}).get("rho"), "{:.2f}")
    dec = ex.get("decisions") or {}
    for k, f in (("dTwoTwoFrozenR", "{:.1e}"),
                 ("dTwoThreeRelDiffTolEightMinus", "{:.1e}"),
                 ("dTwoThreeRelDiffTolNineMinus", "{:.1e}"),
                 ("dTwoThreeRelDiffFloor", "{:.0e}")):
        num(f"exp {k}", (dec.get(k) or {}).get("value"), f)


# --------------------------------------------------------------------------- #
# main
# --------------------------------------------------------------------------- #

def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--commit", default=None,
                    help="target commit (sha, tag, branch); default: "
                    "git rev-parse HEAD")
    ap.add_argument("--allow-stale", action="store_true",
                    help="use results of other commits too (old behaviour); "
                    "they are listed and the title pages carry a warning")
    ap.add_argument("--allow-historical", action="store_true",
                    help="use results/exploratory (historical studies "
                    "without a commit) and list them as historical")
    return ap.parse_args(argv)


def setup_guard(a: argparse.Namespace) -> None:
    global GUARD, ALLOW_HISTORICAL
    rev = a.commit or "HEAD"
    target = provenance.git_resolve(rev)
    if target is None:
        # a commit not in this clone (e.g. a short sha of another clone)
        if a.commit and len(a.commit) >= 7:
            target = a.commit
        else:
            sys.exit(f"make_report: cannot resolve commit {rev!r}")
    GUARD = provenance.Guard(target, allow_stale=a.allow_stale)
    provenance.ACTIVE = GUARD
    ALLOW_HISTORICAL = a.allow_historical
    dirty = provenance.git_dirty_files()
    if dirty and target == provenance.git_head():
        print(f"make_report: note: working tree has modified tracked files "
              f"({len(dirty)}); results are matched against the commit "
              f"{target[:7]} only", flush=True)
    if a.allow_stale:
        print("#" * 72 + "\n# make_report: --allow-stale: results of OTHER "
              "commits are used; the PDFs\n# carry a warning on the title "
              "page\n" + "#" * 72, flush=True)
    else:
        n = purge_generated()
        print(f"make_report: guard on commit {target[:7]}: removed {n} "
              "previously generated outputs", flush=True)


def main(argv: list[str] | None = None) -> int:
    setup_guard(parse_args(argv))
    tests = guard_tests(load_json("tests"))
    bench = guard_bench(load_bench())
    # user-judged convergence iterations (D-060) replace the automatic
    # ones in every number derived from the convergence point
    uc = ucv.load()
    tests = {n: ucv.apply_test(n, d, uc) for n, d in tests.items()}
    bench = [ucv.apply_bench(r, uc) for r in bench]
    table_convergence(tests, bench, uc)

    fig_T0_profiles()
    fig_histories(tests)
    summary_md = fig_bench(bench)
    tests_md = tests_table(tests)
    fig_eta_rho(tests)
    fig_cycles(tests)
    fig_bench_cycles(bench)
    fig_anderson(bench)
    fig_scaling(tests)
    fig_remediation(tests)
    validation_md = table_validation(tests)
    levels_md = table_gamg_levels(tests, bench)
    b10_md = table_b10(bench)
    exploratory_numbers()
    fig_speed(tests)
    fig_fields()
    fig_iteration_histories()

    t0 = tests.get("T0_Re100_np1")
    if t0 is not None:      # no record: the macros stay undefined (pending)
        num("T0 Re100 iterations", t0.get("iterations"), "{}")
        num("T0 Re100 l2u", t0.get("l2rel_u"), "{:.1e}")
        num("T0 Re100 l2v", t0.get("l2rel_v"), "{:.1e}")
        num("T0 Re100 native iterations",
            (t0.get("reference") or {}).get("iterations"), "{}")
    num("commit", GUARD.target[:7] if not GUARD.allow_stale
        else git_commit())
    num("target commit", GUARD.target[:7])
    num("b7 anderson GB", 320 * B7_CELLS / 1e9, "{:.1f}")

    # staleness guard: appendix table and title-page notes
    write_missing_table()
    nrej = GUARD.n_rejected
    num("missing count", nrej, "{}")
    if GUARD.allow_stale and nrej:
        num("stale warning",
            rf"WARNING: built with \texttt{{-{{}}-allow-stale}}. {nrej} results "
            r"from other commits or without provenance are INCLUDED; they "
            r"are listed in Appendix~\ref{app:missing}.")
    elif nrej:
        num("missing note",
            rf"{nrej} results are missing or stale for commit "
            rf"\texttt{{{GUARD.target[:7]}}} and are shown as pending "
            r"(Appendix~\ref{app:missing}).")

    PAPER.mkdir(parents=True, exist_ok=True)
    (PAPER / "numbers.tex").write_text(
        "% Generated by bench/make_report.py - do not edit\n"
        + "".join(f"\\newcommand{{\\{k}}}{{{v}}}\n"
                  for k, v in sorted(numbers.items())))

    md = [
        "# coupledFoam - Report",
        "",
        f"Generated by `bench/make_report.py` from `results/` at commit "
        f"`{git_commit()}`. The LaTeX paper is `report/paper/paper.tex`.",
        "",
        "## Executive summary",
        "",
        summary_md or "_No benchmark results yet._",
        "",
        "## Method",
        "",
        "Block-coupled pressure-velocity solution (4x4 blocks u, v, w, p per "
        "cell) with Rhie-Chow interpolation, pseudo-transient continuation "
        "with adaptive CFL and physicality line search, two-tier remediation "
        "cell sets, single-precision block linear algebra with double-precision "
        "residuals and reductions, block-GAMG preconditioned Krylov solvers. "
        "Full description: `SPEC_coupledFoam.md`, deviations: `DECISIONS.md`, "
        "paper Section 2.",
        "",
        "## Correctness and tests",
        "",
        tests_md or "_No test results yet._",
        "",
        validation_md or "",
        "",
        "## Linear solver",
        "",
        levels_md or "_No block-GAMG statistics yet._",
        "",
        "## Amendment B10 acceptance",
        "",
        b10_md or "_No E/F/G/H benchmark results yet._",
        "",
        "## Figures",
        "",
    ]
    for name, cap in figures:
        md += [f"![{name}](figures/{name}.png)", "", f"*{cap}*", ""]
    for name, cap in figures_pdf_only:
        md += [f"- `paper/figures/{name}.pdf`: {cap}"]
    if figures_pdf_only:
        md += ["", "3D renders (ParaView, `bench/render_fields.py`): "
               "`paper/figures/render_*.png`.", ""]
    md += ["## Limitations and Phase 2", "",
           "See `DECISIONS.md` and the paper, Section 6.", "",
           "## Reproduction", "",
           "```", "source <openfoam2606>/etc/bashrc", "./Allwmake -j 8",
           "~/OF/venv/bin/pytest tests/", "bench/run_bench.py",
           "bench/exploratory_numbers.py", "bench/make_report.py",
           "cd report/paper && latexmk -pdf paper.tex",
           "```", ""]
    grows = GUARD.rows()
    md += ["## Missing or stale results", "",
           f"Target commit `{GUARD.target}`"
           + (" (**--allow-stale: stale items are INCLUDED**)"
              if GUARD.allow_stale else "") + ".", ""]
    if grows:
        md += ["| case | item | found commit | date | status | reason |",
               "|---|---|---|---|---|---|"]
        md += [f"| {r['case']} | {r['item']} | {r['found']} | {r['date']} "
               f"| {r['status']} | {r['reason']} |" for r in grows]
    else:
        md += ["Nothing missing."]
    md += [""]
    if notes:
        md += ["## Generator notes", ""] + [f"- {n}" for n in notes] + [""]
    (REPORT / "REPORT.md").write_text("\n".join(md))
    print(f"figures: {len(figures)}, numbers: {len(numbers)}, notes: {len(notes)}")
    for n in notes:
        print("note:", n)
    print(f"staleness guard: target {GUARD.target[:7]}, "
          f"{GUARD.n_rejected} missing/stale item(s)"
          + (" INCLUDED (--allow-stale)" if GUARD.allow_stale else
             " rendered as pending"))
    for r in grows:
        print(f"  {r['status']:>17}  {r['case']:<34} {r['item']:<40} "
              f"found {r['found']}")
    if GUARD.allow_stale and GUARD.n_rejected:
        print("#" * 72 + "\n# WARNING: --allow-stale: the report contains "
              f"{GUARD.n_rejected} stale item(s)\n" + "#" * 72)
    return 0


if __name__ == "__main__":
    sys.exit(main())

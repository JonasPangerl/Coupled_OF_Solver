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
# Campaign scope D-063: T5 (Ahmed body) is not run by user decision and is
# not expected (the papers state it as deferred); T4b is run only after T4a
# has been accepted, so a missing T4b record is listed with that reason.
EXPECTED_TESTS = [
    "T0_Re100_np1", "T0_Re100_np4", "T0_Re1000_np1", "T0_Re1000_np4",
    "T1_np1", "T1_np4", "T2_np1", "T2_np4",
    "T3_kOmegaSST_np1", "T3_kOmegaSST_np4", "T3_GEKO_np1", "T3_GEKO_np4",
    "T4a_np*", "T4b_np*",
    "T-restart_T1", "T-restart_T3-SST", "T-fpe_*", "T_scaling_T4a",
    "diagnostics_*", "Test-*", "test_env",
]
# reason shown in the appendix for an expected record that is missing
EXPECTED_REASON = {
    "T4b_np*": "not run yet (D-063: T4b only after T4a is accepted)",
}
# cases deferred by user decision (D-063): never listed as missing
DEFERRED_CASES = ("T5",)
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
                          else pat, "result record",
                          EXPECTED_REASON.get(pat, "no result"))
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
    flagged: dict[str, str] = {}
    for f in sorted(FIG_PDF.glob("render_*.png")):
        p = prov.get(f.name)
        stale = False
        if p is None:
            ok = GUARD.check_record(f.stem, {}, "ParaView render "
                                    "(no render provenance)")
            stale = True
            why = "render without provenance record"
        elif not p.get("cfUsed"):
            ok = True       # simpleFoam reference and geometry only
        else:
            ok = GUARD.check_record(
                f.stem, {"gitCommit": p.get("cfCommit"),
                         "timestamp": p.get("date")},
                f"ParaView render of run/{p.get('cfRun')}")
            stale = not GUARD.matches(p.get("cfCommit"))
            why = (f"coupledFoam run run/{p.get('cfRun')} of commit "
                   f"{str(p.get('cfCommit') or 'unknown')[:15]}")
        if not ok:
            f.unlink()
            notes.append(f"render {f.name}: stale, removed (pending)")
        elif stale:
            # kept with --allow-stale: flag it in the caption
            # (\cfrenderflag{<case>} in both papers)
            case = f.stem.split("_")[1] if f.stem.count("_") >= 2 else f.stem
            flagged.setdefault(case, why)
    for case, why in flagged.items():
        num(f"render flag {case}",
            r"\textbf{Stale render:} " + tex_escape(why)
            + r", not the target commit (Appendix~\ref{app:missing}).")


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
                       "Amendment B10: configurations F (fixed relTol), "
                       "G (Anderson) and H (fixed K-cycle) against C (the "
                       "defaults, V-cycle since D-043; E is identical to C "
                       "and not evaluated, D-068); d = X/C$-$1 (positive: "
                       "X slower). Pass "
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


def _speed_record(tests: dict, rec: str, cfd: str, sfd: str) -> dict | None:
    """Wall time and CPU-hours to convergence of both solvers for one case,
    with flags for 'did not reach the criterion' (then: whole run, a lower
    bound for the time to convergence)."""
    d = tests.get(rec)
    if not d:
        return None
    ref = d.get("reference") or {}
    n_cf = d.get("iterations")
    it_cf = d.get("iterationsToR") or d.get("iterationsToR_coupled")
    wall_cf_run = d.get("wallSecondsSolver") or d.get("wallSeconds")
    cpu_cf_run = d.get("cpuHoursSolver") or d.get("cpuHours")
    if not (n_cf and wall_cf_run and cpu_cf_run):
        return None
    # records without an iteration count to the residual target (T3: the
    # residual target was not reached in the 3000 iterations) count as not
    # converged; the bar is then the whole run
    conv_cf = bool(it_cf)
    cft = cf_timeline(RUN / cfd)
    if conv_cf and cft is not None and it_cf <= len(cft["t"]):
        frac = cft["t"][it_cf - 1] / cft["t"][-1]
    elif conv_cf:
        frac = it_cf / n_cf
    else:
        frac = 1.0
    it_sf = ref.get("convergedAt")
    n_sf = ref.get("iterations")
    conv_sf = it_sf is not None
    ta = ref.get("timingAllrun") or {}
    wall_sf = ref.get("wallSecondsSolver") or ref.get("wallSeconds") \
        or ta.get("wallSeconds")
    cpu_sf = ta.get("cpuHours") or ref.get("cpuHours")
    nproc_sf = ref.get("nProcs") or 1
    if cpu_sf is None and wall_sf:
        cpu_sf = wall_sf * nproc_sf / 3600.0   # serial reference: CPU = wall
    if not wall_sf:
        return None
    return {
        "wall_cf": wall_cf_run * frac, "cpu_cf": cpu_cf_run * frac,
        "it_cf": it_cf if conv_cf else n_cf, "conv_cf": conv_cf,
        "wall_sf": wall_sf, "cpu_sf": cpu_sf,
        "it_sf": it_sf if conv_sf else n_sf, "conv_sf": conv_sf,
        "np_cf": d.get("nProcs") or 1, "np_sf": nproc_sf,
    }


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
# single precision (amendment D11; D-062: the reference precision is DP)
# --------------------------------------------------------------------------- #
#
# Records (schema of branch amend-d-forces, bench/SCHEMA_precision.md): an SP
# test record is named like its DP counterpart with the suffix "_sp"
# (results/tests/T1_np1_sp.json) or carries "precision": "SP" and names its
# DP record in "dpRecord". Fields used when present: precision, meshShift
# (vector or {"vector": ...}), checkMeshDiff (list of checks that differ
# from DP; empty = gate passed), staticSetSizeDiff, Cd_rel_to_DP, and the
# usual timing/memory fields. Benchmark configurations F1 (simpleFoam SP)
# and F2 (coupledFoam SP) are compared with their DP counterparts
# (record field "dpConfig", default A for F1 and C for F2).

SP_SUFFIX = "_sp"
SP_TOL = 0.005          # D11 verdict: monitored quantity within 0.5 % of DP
SP_BENCH_DP = {"F1": "A", "F2": "C"}
# monitored integral quantity per case family: (label, record keys)
SP_QUANTITY = (("C_d", ("Cd", "Cd_mean")), ("dp", ("dp",)),
               ("x_r/h", ("xr_over_h",)))


def _sp_pairs(tests: dict) -> list[tuple[str, dict, str, dict | None]]:
    """(SP name, SP record, DP name, DP record or None), sorted."""
    out = []
    for n, d in sorted(tests.items()):
        is_sp = n.endswith(SP_SUFFIX) or str(d.get("precision", "")).upper() == "SP"
        if not is_sp:
            continue
        dpn = d.get("dpRecord") or (n[:-len(SP_SUFFIX)] if n.endswith(SP_SUFFIX)
                                    else None)
        out.append((n, d, dpn or "?", tests.get(dpn) if dpn else None))
    return out


def _sp_quantity(sp: dict, dp: dict | None) -> tuple[str, float | None,
                                                      float | None, float | None]:
    """(label, DP value, SP value, relative deviation SP vs DP)."""
    for lab, keys in SP_QUANTITY:
        vs = next((sp[k] for k in keys if sp.get(k) is not None), None)
        if vs is None:
            continue
        vd = next((dp[k] for k in keys if dp and dp.get(k) is not None), None)
        rel = sp.get("Cd_rel_to_DP") if lab == "C_d" else None
        if rel is None and vd:
            rel = (vs - vd) / abs(vd)
        return lab, vd, vs, rel
    return "-", None, None, sp.get("Cd_rel_to_DP")


def _sp_gate(sp: dict) -> tuple[bool | None, str]:
    """checkMesh gate of D5.3: (passed, text). None: not recorded."""
    if str(sp.get("status", "")).lower() == "sp-geometry-fail" \
            or sp.get("spGeometryFail"):
        return False, "SP-geometry-fail"
    diff = sp.get("checkMeshDiff")
    if diff is None:
        return None, "not recorded"
    if isinstance(diff, dict):
        diff = [k for k, v in diff.items() if v]
    if diff:
        return False, "differs: " + ", ".join(str(x) for x in diff)[:60]
    return True, "passed"


def _shift_norm(v) -> str:
    if isinstance(v, dict):
        v = v.get("vector") or v.get("translate")
    try:
        return f"{float(np.linalg.norm(np.asarray(v, dtype=float))):.3g} m"
    except (TypeError, ValueError):
        return "n/a"


def _per_iter(d: dict, run_dir: Path | None) -> dict:
    """coupledFoam time per iteration split (s): assembly, linear solve,
    turbulence, other (incl. I/O; or t_io if recorded), total. From the
    run log if readable, else from the record totals."""
    cft = cf_timeline(run_dir) if run_dir is not None and run_dir.is_dir() \
        else None
    if cft is not None:
        a, s, t = (float(np.nanmean(cft[k])) for k in ("tAsm", "tSolve", "tTurb"))
        tot = float(np.nanmean(cft["tIter"]))
    else:
        n = d.get("iterations") or d.get("iterationsRun")
        if not n or d.get("t_assembly") is None:
            return {}
        a, s, t = (d.get(k, 0.0) / n for k in ("t_assembly", "t_linsolve",
                                                "t_turb"))
        w = d.get("wallSecondsSolver") or d.get("wallSeconds")
        tot = w / n if w else a + s + t
    io = d.get("t_io")
    n = d.get("iterations") or d.get("iterationsRun") or 1
    other = io / n if io is not None else max(tot - a - s - t, 0.0)
    return {"asm": a, "solve": s, "turb": t, "other": other, "total": tot,
            "ioRecorded": io is not None}


def sp_section(tests: dict, bench: list[dict]) -> None:
    """Section 'Single precision' (D11): performance table SP vs DP for both
    solvers, per-iteration breakdown, verdict table, convergence-floor
    figure. Nothing is written without SP records (the papers show their
    pending lines)."""
    pairs = _sp_pairs(tests)
    brows = [d for d in bench if d.get("config") in SP_BENCH_DP]
    num("sp n cases", len(pairs), "{}")
    if not pairs and not brows:
        notes.append("single precision: no SP records (*_sp, F1/F2) yet")
        return
    # ---- performance: wall, CPU-h, peak RSS, both solvers
    perf = []

    def ratio(a, b):
        return fmt(a / b, "{:.2f}") if a and b else "n/a"

    for n, sp, dpn, dp in pairs:
        dp = dp or {}
        for solver, s, d in (("coupledFoam", sp, dp),
                             ("simpleFoam", sp.get("reference") or {},
                              dp.get("reference") or {})):
            ws = s.get("wallSecondsSolver") or s.get("wallSeconds")
            wd = d.get("wallSecondsSolver") or d.get("wallSeconds")
            cs = s.get("cpuHoursSolver") or s.get("cpuHours")
            cd = d.get("cpuHoursSolver") or d.get("cpuHours")
            ms, md = s.get("peakRSS_GB_sum"), d.get("peakRSS_GB_sum")
            if not any((ws, cs, ms)):
                continue
            perf.append([dpn, solver, fmt(wd), fmt(ws), ratio(wd, ws),
                         fmt(cd, "{:.3g}"), fmt(cs, "{:.3g}"), ratio(cd, cs),
                         fmt(md, "{:.3g}"), fmt(ms, "{:.3g}"), ratio(md, ms)])
            if solver == "coupledFoam":
                num(f"sp speedup wall {dpn}", wd / ws if wd and ws else None,
                    "{:.2f}")
    for cfg, dpcfg0 in SP_BENCH_DP.items():
        for c in sorted({d["case"] for d in brows if d["config"] == cfg}):
            dpcfg = next((d.get("dpConfig") for d in brows
                          if d["config"] == cfg and d.get("dpConfig")), dpcfg0)
            ws, wd = _med(bench, c, cfg, "wall_to_conv_s"), \
                _med(bench, c, dpcfg, "wall_to_conv_s")
            cs, cd = _med(bench, c, cfg, "cpu_to_conv_h"), \
                _med(bench, c, dpcfg, "cpu_to_conv_h")
            ms, md = _med(bench, c, cfg, "peakRSS_GB_sum"), \
                _med(bench, c, dpcfg, "peakRSS_GB_sum")
            perf.append([f"bench {c} ({cfg} vs {dpcfg})",
                         "simpleFoam" if cfg == "F1" else "coupledFoam",
                         fmt(wd), fmt(ws), ratio(wd, ws), fmt(cd, "{:.3g}"),
                         fmt(cs, "{:.3g}"), ratio(cd, cs), fmt(md, "{:.3g}"),
                         fmt(ms, "{:.3g}"), ratio(md, ms)])
    if perf:
        write_table(
            "sp_performance",
            ["case", "solver", "wall DP [s]", "wall SP [s]", "DP/SP",
             "CPU-h DP", "CPU-h SP", "DP/SP", "RSS DP [GB]", "RSS SP [GB]",
             "DP/SP"], perf,
            "Single precision (SP) against double precision (DP, the "
            "reference build, D-062): solver wall-clock time, CPU-hours and "
            "peak memory (sum over ranks) of the same test for both solvers. "
            "Ratios DP/SP above one: SP faster or smaller. Test rows: whole "
            "run of the test; benchmark rows: time to convergence (median).",
            "tab:spperf", resize=True)
    # ---- per-iteration breakdown (coupledFoam)
    brk = []
    for n, sp, dpn, dp in pairs:
        for prec, d, rd in (("DP", dp, RUN / dpn), ("SP", sp, RUN / n)):
            if not d:
                continue
            b = _per_iter(d, rd if run_ok(rd) else None)
            if not b:
                continue
            brk.append([dpn, prec] + [fmt(1e3 * b[k], "{:.3g}") for k in
                                      ("asm", "solve", "turb", "other", "total")])
    if brk:
        write_table(
            "sp_breakdown",
            ["case", "precision", "assembly [ms/it]", "linear solve [ms/it]",
             "turbulence [ms/it]", "other incl. I/O [ms/it]", "total [ms/it]"],
            brk,
            "coupledFoam wall time per outer iteration by component, DP and "
            "SP. ``Other'' is the remainder of the iteration (residual "
            "evaluation, line search, bookkeeping and field output); it is "
            "the recorded I/O time where the record carries one (t_io).",
            "tab:spbreak", resize=True)
    # ---- verdict (D11)
    ver = []
    n_ok = 0
    for n, sp, dpn, dp in pairs:
        lab, vd, vs, rel = _sp_quantity(sp, dp)
        gate, gtxt = _sp_gate(sp)
        reasons = []
        if gate is False:
            reasons.append(f"checkMesh gate {gtxt}")
        if rel is None:
            reasons.append("no DP value to compare")
        elif abs(rel) > SP_TOL:
            reasons.append(f"{lab} deviates {100 * rel:+.2f}% from DP")
        if sp.get("rc") not in (None, 0) or sp.get("fpeTrap"):
            reasons.append("run failed")
        if reasons and (gate is False or (rel is not None and abs(rel) > SP_TOL)
                        or "run failed" in reasons):
            verdict = "SP not usable: " + "; ".join(reasons)
        elif gate is None or rel is None:
            verdict = "undetermined: " + ("; ".join(reasons) if reasons
                                          else "checkMesh gate not recorded")
        else:
            verdict = "SP usable"
            n_ok += 1
        ver.append([dpn, lab, fmt(vd, "{:.5g}"), fmt(vs, "{:.5g}"),
                    fmt(100 * rel if rel is not None else None, "{:+.3f}") +
                    ("%" if rel is not None else ""),
                    gtxt, fmt(sp.get("staticSetSizeDiff"), "{}"),
                    _shift_norm(sp.get("meshShift")), verdict])
    if ver:
        write_table(
            "sp_verdict",
            ["case", "quantity", "DP", "SP", "SP vs DP", "checkMesh gate",
             "static set diff. [cells]", "origin shift", "verdict"], ver,
            "Single-precision verdict per case (amendment D11 with DP as the "
            "reference, D-062): ``SP usable'' if the monitored quantity "
            "($C_d$; $\\Delta p$ for T1, $x_r/h$ for T2) of coupledFoam in SP "
            "is within 0.5\\,\\% of DP and the checkMesh gate of the shifted SP "
            "mesh passed (no check failing in SP that passes in DP, no "
            "negative volumes), else ``SP not usable'' with the reason.",
            "tab:spverdict", resize=True)
        num("sp n usable", n_ok, "{}")
    # ---- convergence floor: R_n in SP against DP
    hist = [(dpn, (dp or {}).get("history", {}).get("R"),
             sp.get("history", {}).get("R")) for _, sp, dpn, dp in pairs]
    hist = [h for h in hist if h[1] or h[2]]
    if hist:
        ncol = min(3, len(hist))
        nrow = int(np.ceil(len(hist) / ncol))
        fig, axs = plt.subplots(nrow, ncol, figsize=(6.5, 2.3 * nrow),
                                squeeze=False)
        for ax, (lab, rd, rs) in zip(axs.flat, hist):
            for r, col, ls, name in ((rd, C_COUPLED, "-", "DP"),
                                     (rs, "#E69F00", "--", "SP")):
                if r:
                    ax.semilogy(np.arange(1, len(r) + 1), r, color=col, ls=ls,
                                lw=1.3, label=f"coupledFoam {name}")
            ax.axhline(1e-5, color="k", lw=0.6, ls=":",
                       label="SP residual target (D7)")
            ax.set_title(lab, fontsize=8)
            ax.set_xlabel("outer iteration", fontsize=8)
            ax.tick_params(labelsize=7)
        for ax in list(axs.flat)[len(hist):]:
            ax.axis("off")
        for r in range(nrow):
            axs[r, 0].set_ylabel("combined residual $R_n$", fontsize=8)
        axs.flat[0].legend(fontsize=7, loc="upper right")
        save(fig, "sp_convergence_floor",
             "Convergence floor: combined residual $R_n$ of coupledFoam in "
             "double (solid) and single precision (dashed) on the same case; "
             "dotted: the SP residual target of the precision profile (D7).")


# --------------------------------------------------------------------------- #
# metrics used by other solvers' publications (cost per cell, memory per
# cell, residual reduction, scaling, grid) - computed from existing records
# --------------------------------------------------------------------------- #

def _cells_of(d: dict, run_dir: Path | None) -> int | None:
    """Cell count: record field, else the owner-file headers of the run
    (serial or the sum over processor*/)."""
    for k in ("meshCells", "nCells", "cells"):
        if isinstance(d.get(k), (int, float)) and d.get(k):
            return int(d[k])
    fc = d.get("fieldCompare") or {}
    if fc.get("nCells"):
        return int(fc["nCells"])
    if run_dir is None or not run_dir.is_dir():
        return None
    import gzip  # noqa: PLC0415
    import re  # noqa: PLC0415
    owners = []
    for pat in ("processor*/constant/polyMesh/owner",
                "processor*/constant/polyMesh/owner.gz",
                "constant/polyMesh/owner", "constant/polyMesh/owner.gz"):
        owners = sorted(run_dir.glob(pat))
        if owners:
            break
    tot = 0
    for f in owners:
        try:
            opener = gzip.open if f.suffix == ".gz" else open
            with opener(f, "rb") as fh:
                head = fh.read(4096).decode("latin-1")
        except OSError:
            return None
        m = re.search(r"nCells:\s*(\d+)", head)
        if not m:
            return None
        tot += int(m.group(1))
    return tot or None


def metrics_section(tests: dict) -> None:
    """Tables metrics_cost, metrics_residual, metrics_scaling,
    metrics_grid (see the section 'Metrics for comparison with other
    solvers' of the papers)."""
    rows, rres = [], []
    for lab, rec, cfd, sfd in SPEED_CASES:
        if lab.split()[0] in DEFERRED_CASES:
            continue
        d = tests.get(rec)
        if not d:
            continue
        ref = d.get("reference") or {}
        cells = _cells_of(d, RUN / cfd) or _cells_of(ref, RUN / sfd)
        for solver, r, conv_it, conv_cpu in (
                ("coupledFoam", d,
                 d.get("itersToConv") or d.get("iterationsToR"),
                 d.get("cpuHoursToConv")),
                ("simpleFoam", ref, ref.get("itersToConv") or ref.get("convergedAt"),
                 ref.get("cpuHoursToConv"))):
            n = r.get("iterations") or r.get("iterationsRun")
            w = r.get("wallSecondsSolver") or r.get("wallSeconds")
            c = r.get("cpuHoursSolver") or r.get("cpuHours")
            m = r.get("peakRSS_GB_sum")
            if m is None and r.get("peakRSS_MB_sum") is not None:
                m = r["peakRSS_MB_sum"] / 1024.0
            if solver == "simpleFoam" and c is None and w:
                c = w * (r.get("nProcs") or 1) / 3600.0
            if not (n and w):
                continue
            mc = cells / 1e6 if cells else None
            npr = r.get("nProcs") or d.get("nProcs") or 1
            rows.append([
                lab, solver, fmt(cells, "{:,}").replace(",", r"\,") if cells else "n/a",
                fmt(npr, "{}"), fmt(n, "{}"), fmt(w / n, "{:.3g}"),
                fmt(c * 3600 / n / mc if c and mc else None, "{:.3g}"),
                fmt(mc * 1e6 * n / (c * 3600) / 1e3 if c and mc else None,
                    "{:.3g}"),
                fmt(m / mc if m and mc else None, "{:.3g}"),
                fmt(conv_it, "{}"),
                fmt(conv_cpu / mc if conv_cpu and mc else None, "{:.3g}")])
        # residual reduction (orders of magnitude below the first value)
        cft, sft = cf_timeline(RUN / cfd), sf_timeline(RUN / sfd)
        for solver, R, t in (
                ("coupledFoam R", cft["R"] if cft else None,
                 cft["t"] if cft else None),
                ("simpleFoam p", sft["res"].get("p") if sft else None,
                 sft["t"] if sft else None)):
            if R is None or len(R) < 2 or not R[0]:
                continue
            cells_ = []
            for k in (2, 3, 4):
                hit = np.nonzero(np.asarray(R) <= R[0] * 10.0 ** (-k))[0]
                cells_.append(f"{hit[0] + 1} / {t[hit[0]]:.3g}" if len(hit)
                              else "not reached")
            rres.append([lab, solver, f"{R[0]:.2e}"] + cells_)
    if rows:
        write_table(
            "metrics_cost",
            ["case", "solver", "cells", "ranks", "iterations", "wall/it [s]",
             "CPU-s/it/Mcell", "k cell-it per CPU-s", "RSS [GB/Mcell]",
             "it. to conv.", "CPU-h to conv./Mcell"], rows,
            "Cost and memory normalised by the mesh size, as reported in "
            "solver validation reports: CPU time per outer iteration and "
            "million cells, throughput in thousand cell-iterations per "
            "CPU-second, peak memory (sum over ranks) per million cells, "
            "and CPU-hours to convergence per million cells. Test records "
            "(whole run of the test, shared machine, D-007). On meshes below "
            "about 10$^5$ cells the fixed per-process memory of OpenFOAM "
            "dominates the memory per cell.", "tab:metcost", resize=True)
    if rres:
        write_table(
            "metrics_residual",
            ["case", "residual", "first value", "to 1e-2 x first",
             "to 1e-3 x first", "to 1e-4 x first"], rres,
            "Residual reduction: iteration / wall-clock time [s] at which the "
            "residual first falls two, three and four orders of magnitude "
            "below its first value (coupledFoam combined residual $R$, "
            "simpleFoam initial residual of $p$). The two residuals are "
            "normalised differently (Section on residual norms), so only the "
            "reductions of each solver, not the levels, are comparable.",
            "tab:metres", resize=True)
    # ---- scaling table (T-scaling record)
    d = next((v for k, v in tests.items()
              if k.startswith(("T_scaling", "T-scaling"))), None)
    if d:
        srows = []
        per = {s: {int(k): v for k, v in ((d.get(s) or {}).get("timePerIter_s")
                                          or {}).items() if v}
               for s in ("coupledFoam", "simpleFoam")}
        ranks = sorted(set(per["coupledFoam"]) | set(per["simpleFoam"]))
        for nr in ranks:
            row = [f"{nr}"]
            for s in ("coupledFoam", "simpleFoam"):
                t = per[s].get(nr)
                n0 = min(per[s]) if per[s] else None
                eff = (per[s][n0] * n0 / (t * nr)) if t and n0 else None
                row += [fmt(t, "{:.3g}"), fmt(eff, "{:.2f}")]
            cf_e = row[2] if row[2] != "n/a" else None
            sf_e = row[4] if row[4] != "n/a" else None
            row.append(fmt(float(cf_e) / float(sf_e), "{:.2f}")
                       if cf_e and sf_e and float(sf_e) else "n/a")
            srows.append(row)
        if srows:
            write_table(
                "metrics_scaling",
                ["ranks", "cF wall/it [s]", "cF efficiency", "sF wall/it [s]",
                 "sF efficiency", "efficiency ratio cF/sF"], srows,
                "Strong scaling on "
                f"{tex_escape(str(d.get('mesh', 'T4a')))}: wall time per "
                "iteration and parallel efficiency $E(n)=t(n_0)n_0/(t(n)n)$ "
                "of coupledFoam (cF) and simpleFoam (sF); pass criterion: "
                "efficiency ratio at the largest rank count at least 0.8.",
                "tab:metscal")
    # ---- grid: T4a against T4b
    a, b = tests.get("T4a_np10"), tests.get("T4b_np10")
    if a and b and a.get("Cd") is not None and b.get("Cd") is not None:
        grows = []
        for lab, d in (("T4a", a), ("T4b", b)):
            ref = d.get("reference") or {}
            grows.append([lab, fmt(_cells_of(d, RUN / f"{lab}_np10"), "{}"),
                          fmt(d.get("Cd"), "{:.4f}"), fmt(ref.get("Cd"), "{:.4f}"),
                          fmt(d.get("itersToConv"), "{}"),
                          fmt(ref.get("itersToConv"), "{}"),
                          fmt(d.get("cpuHoursToConv"), "{:.3g}"),
                          fmt(ref.get("cpuHoursToConv"), "{:.3g}")])
        write_table(
            "metrics_grid",
            ["mesh", "cells", "Cd cF", "Cd sF", "it. to conv. cF",
             "it. to conv. sF", "CPU-h to conv. cF", "CPU-h to conv. sF"],
            grows,
            "Mesh dependence on the two motorBike meshes: window-mean $C_d$ "
            "and the cost to the stationary window mean (D-042) of "
            "coupledFoam (cF) and simpleFoam (sF). A solver whose iteration "
            "count to convergence grows little with the mesh size keeps its "
            "advantage on production meshes.", "tab:metgrid", resize=True)


# --------------------------------------------------------------------------- #
# speed-up text and the motorbike speed-up against the cached references
# --------------------------------------------------------------------------- #

def speed_notes(tests: dict) -> None:
    """\\cfSpeedNotConvNote: which test runs of the speed-up figure did not
    reach their criterion (generated instead of a hand-written claim)."""
    nc = []
    for lab, rec, cfd, sfd in SPEED_CASES:
        r = _speed_record(tests, rec, cfd, sfd)
        if not r:
            continue
        who = [s for s, k in (("coupledFoam", "conv_cf"), ("simpleFoam", "conv_sf"))
               if not r[k]]
        if who:
            nc.append(f"{lab} ({' and '.join(who)})")
    if nc:
        txt = ("The following runs did not reach their criterion within their "
               "iteration limit, so their bars compare the run lengths: "
               + ", ".join(nc) + ".")
    else:
        txt = "All runs of the figure reached their criterion."
    num("speed not conv note", tex_escape(txt))


def _flag(d: dict, key: str) -> str:
    v = d.get(key)
    if v is None:
        v = (d.get("reference") or {}).get(key)
    return "not recorded" if v is None else ("yes" if v else "no")


def table_wake_speedup(tests: dict) -> None:
    """tables/wake_speedup.tex: T4a/T4b coupledFoam test run against the
    cached simpleFoam reference (D-059), with the window rule used and the
    fairness flags of the reference (harness-fix: referenceNoPotentialStart,
    referenceTimingConditionsUnknown)."""
    rows = []
    for name in sorted(tests):
        if not name.startswith(("T4a_np", "T4b_np")):
            continue
        d = tests[name]
        ref = d.get("reference") or {}
        wc, ws = d.get("wallToConv_s"), ref.get("wallToConv_s")
        cc, cs = d.get("cpuHoursToConv"), ref.get("cpuHoursToConv")
        if wc is None and ws is None:
            continue
        spw = d.get("speedupWall") or (ws / wc if ws and wc else None)
        spc = d.get("speedupCpu") or (cs / cc if cs and cc else None)
        W, Wr = d.get("W"), ref.get("W")
        common = d.get("commonWindow") or (W is not None and W == Wr)
        rule = (f"common W = {W}" if common else
                f"per-run W = {W} / {Wr} (earlier rule)")
        sens = (d.get("sensitivity") or {}).get("perRunWindow") \
            or d.get("perRunWindow") or {}
        sens_txt = (f"{sens['speedupWall']:.2f} / {sens['speedupCpu']:.2f}"
                    if sens.get("speedupWall") and sens.get("speedupCpu")
                    else ("= left" if not common else "n/a"))
        rows.append([name, fmt(d.get("itersToConv"), "{}"),
                     fmt(ref.get("itersToConv"), "{}"), fmt(wc), fmt(ws),
                     fmt(cc, "{:.3g}"), fmt(cs, "{:.3g}"),
                     fmt(spw, "{:.2f}"), fmt(spc, "{:.2f}"), rule, sens_txt,
                     _flag(d, "referenceNoPotentialStart"),
                     _flag(d, "referenceTimingConditionsUnknown")])
        k = name.split("_")[0]
        num(f"wake speedup wall {k}", spw, "{:.2f}")
        num(f"wake speedup cpu {k}", spc, "{:.2f}")
    if not rows:
        notes.append("wake speed-up: no T4 records")
        return
    write_table(
        "wake_speedup",
        ["run", "it. cF", "it. sF", "wall cF [s]", "wall sF [s]",
         "CPU-h cF", "CPU-h sF", "speed-up wall", "speed-up CPU",
         "window", "per-run-window speed-up (wall / CPU)",
         "ref. without potentialFoam start", "ref. timing conditions unknown"],
        rows,
        "Motorbike: coupledFoam (cF) test run against the cached simpleFoam "
        "(sF) reference, iterations, wall-clock time and CPU-hours to the "
        "stationary window mean (D-042, common window D-068), and the "
        "speed-up sF/cF. The window of the earlier rule was derived from each "
        "run's own budget and fixed the earliest possible convergence "
        "of the longer reference run; its speed-up is a sensitivity value. "
        "The last two columns disclose the conditions of the reference run.",
        "tab:wakespeed", resize=True)


# --------------------------------------------------------------------------- #
# remediation per category (rem-cat: meshQuality, badMesh, processor, wall)
# --------------------------------------------------------------------------- #

# per-category cell counts in a record: a dict under one of these keys,
# values int or {"cells"|"nCells"|"size": int, ...}
REM_CAT_KEYS = ("remediationCategories", "staticCategories", "remediation")
REM_CAT_ORDER = ("meshQuality", "badMesh", "processor", "wall")


def _rem_categories(d: dict) -> dict[str, int]:
    for k in REM_CAT_KEYS:
        v = d.get(k)
        if isinstance(v, dict) and v:
            out = {}
            for cat, x in v.items():
                if isinstance(x, dict):
                    x = x.get("cells", x.get("nCells", x.get("size")))
                if isinstance(x, (int, float)):
                    out[cat] = int(x)
            if out:
                return out
    # flat fields staticCells_<cat> / nCells_<cat>
    out = {}
    for k, x in d.items():
        for pre in ("staticCells_", "remCells_"):
            if k.startswith(pre) and isinstance(x, (int, float)):
                out[k[len(pre):]] = int(x)
    return out


def table_remediation_categories(tests: dict) -> None:
    """tables/remediation_categories.tex: cells in the remediation sets per
    test run. Works with the old records (static set nStat/staticCells,
    dynamic set nDyn) and with per-category records (rem-cat)."""
    recs = []
    cats: list[str] = []
    for n, d in sorted(tests.items()):
        h = d.get("history") or {}
        c = _rem_categories(d)
        if not (c or d.get("staticCells") is not None or h.get("nStat")
                or h.get("nDyn")):
            continue
        for k in c:
            if k not in cats:
                cats.append(k)
        recs.append((n, d, c, h))
    if not recs:
        notes.append("remediation categories: no data")
        return
    cats.sort(key=lambda k: (REM_CAT_ORDER.index(k) if k in REM_CAT_ORDER
                             else len(REM_CAT_ORDER), k))

    def last(v):
        v = [x for x in (v or []) if x is not None]
        return v[-1] if v else None

    def vmax(v):
        v = [x for x in (v or []) if x is not None]
        return max(v) if v else None

    rows = []
    for n, d, c, h in recs:
        stat = d.get("staticCells")
        if stat is None:
            stat = last(h.get("nStat"))
        cells = _cells_of(d, RUN / n)
        rows.append([n, fmt(cells, "{:,}").replace(",", r"\,")
                     if cells else "n/a"]
                    + [fmt(c.get(k), "{}") for k in cats]
                    + [fmt(stat, "{}"),
                       fmt(100 * stat / cells if stat is not None and cells
                           else None, "{:.2f}"),
                       fmt(vmax(h.get("nDyn")), "{}"),
                       fmt(last(h.get("nDyn")), "{}")])
    write_table(
        "remediation_categories",
        ["run", "cells"] + [f"{k}" for k in cats]
        + ["static total", "static [%]", "dynamic max", "dynamic final"],
        rows,
        "Remediation cells per test run: the pre-selected cells per category"
        + (" (" + ", ".join(cats) + ")" if cats else
           " (records of this commit carry no per-category counts yet)")
        + ", the pre-selected total and its share of the mesh, and the "
        "largest and final size of the dynamic set.",
        "tab:remcat", resize=True)


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
    # report-d additions: single precision (D11), metrics of other solvers'
    # validation reports, remediation per category
    sp_section(tests, bench)
    metrics_section(tests)
    table_remediation_categories(tests)
    speed_notes(tests)
    table_wake_speedup(tests)

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

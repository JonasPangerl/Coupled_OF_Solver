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

Figures whose data is missing are skipped with a note (the report lists
them), so the generator runs at every stage of the project.

Usage: bench/make_report.py            (no OpenFOAM environment needed)
"""

from __future__ import annotations

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
from cflib import logs, post  # noqa: E402

RESULTS = REPO / "results"
RUN = REPO / "run"
REPORT = REPO / "report"
FIG_PNG = REPORT / "figures"
PAPER = REPORT / "paper"
FIG_PDF = PAPER / "figures"
TABLES = PAPER / "tables"

# Colours (spec 15): native grey, coupled blue; configuration shades
C_NATIVE = "#7f7f7f"
C_NATIVE2 = "#b0b0b0"
C_COUPLED = "#1f5fbf"
C_COUPLED2 = "#7fa6e0"
CONFIG_COLOR = {"A": C_NATIVE2, "B": C_NATIVE, "C": C_COUPLED, "D": C_COUPLED2,
                "E": "#6baed6", "F": "#08306b", "G": "#d95f02"}
CONFIG_LABEL = {
    "A": "simpleFoam (tutorial)",
    "B": "simpleFoam SIMPLEC",
    "C": "coupledFoam",
    "D": "coupledFoam, blockDiagonal",
    "E": "coupledFoam, V-cycle",
    "F": "coupledFoam, fixed relTol",
    "G": "coupledFoam, Anderson",
}

plt.rcParams.update({
    "figure.dpi": 150, "savefig.dpi": 150, "font.size": 9,
    "axes.grid": True, "grid.alpha": 0.3, "legend.frameon": False,
    "pdf.fonttype": 42, "ps.fonttype": 42,
})

notes: list[str] = []
numbers: dict[str, str] = {}
figures: list[tuple[str, str]] = []


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
            .replace("&", r"\&").replace("#", r"\#"))


def write_table(name: str, header: list[str], rows: list[list], caption: str,
                label: str) -> str:
    TABLES.mkdir(parents=True, exist_ok=True)
    cols = "l" + "r" * (len(header) - 1)
    lines = [
        r"\begin{table}[t]", r"\centering", r"\small",
        rf"\caption{{{caption}}}", rf"\label{{{label}}}",
        rf"\begin{{tabular}}{{{cols}}}", r"\toprule",
        " & ".join(tex_escape(h) for h in header) + r" \\", r"\midrule",
    ]
    for r in rows:
        lines.append(" & ".join(tex_escape(c) for c in r) + r" \\")
    lines += [r"\bottomrule", r"\end{tabular}", r"\end{table}", ""]
    (TABLES / f"{name}.tex").write_text("\n".join(lines))
    # Markdown version for REPORT.md
    md = ["| " + " | ".join(header) + " |",
          "|" + "---|" * len(header)]
    md += ["| " + " | ".join(str(c) for c in r) + " |" for r in rows]
    return "\n".join(md)


def fmt(v, f="{:.3g}"):
    return "n/a" if v is None else (f.format(v) if not isinstance(v, str) else v)


# --------------------------------------------------------------------------- #
# figures
# --------------------------------------------------------------------------- #

def fig_T0_profiles() -> None:
    for re_ in (100, 1000):
        cp = RUN / f"T0_Re{re_}_np1"
        rp = RUN / f"ref_T0_Re{re_}"
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
        ref = d.get("reference", {})
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
        _ = ref


def fig_bench(bench: dict) -> str:
    """Wall and CPU time to convergence per case and configuration."""
    rows = [d for d in bench.values() if "case" in d]
    if not rows:
        notes.append("benchmark: no results yet")
        return ""
    cases = sorted({d["case"] for d in rows})
    cfgs = [c for c in "ABCDEFG" if any(d["config"] == c for d in rows)]

    def med(c, cfg, key):
        v = [d.get(key) for d in rows
             if d["case"] == c and d["config"] == cfg and d.get(key)]
        return float(np.median(v)) if v else None

    for key, ylabel, fname in (
        ("wall_to_conv_s", "wall time to convergence [s]", "bench_wall"),
        ("cpu_to_conv_h", "CPU time to convergence [CPU-h]", "bench_cpu"),
        ("iters_to_conv", "iterations to convergence", "bench_iters"),
    ):
        fig, ax = plt.subplots(figsize=(6.5, 2.8))
        width = 0.8 / max(1, len(cfgs))
        for j, cfg in enumerate(cfgs):
            vals = [med(c, cfg, key) or 0 for c in cases]
            ax.bar(np.arange(len(cases)) + j * width, vals, width,
                   color=CONFIG_COLOR[cfg], label=CONFIG_LABEL[cfg])
        ax.set_xticks(np.arange(len(cases)) + width * (len(cfgs) - 1) / 2)
        ax.set_xticklabels(cases)
        ax.set_ylabel(ylabel)
        ax.set_yscale("log")
        ax.legend(fontsize=7, ncol=2)
        save(fig, fname, f"Benchmark: {ylabel} (median of the repeats).")

    # Time per iteration breakdown of coupledFoam
    br = [(c, med(c, "C", "t_assembly"), med(c, "C", "t_linsolve"),
           med(c, "C", "t_turb"), med(c, "C", "iters_to_conv")) for c in cases]
    br = [b for b in br if b[4]]
    if br:
        fig, ax = plt.subplots(figsize=(6.5, 2.6))
        x = np.arange(len(br))
        a = [b[1] / b[4] for b in br]
        s = [b[2] / b[4] for b in br]
        t = [b[3] / b[4] for b in br]
        ax.bar(x, a, color="#9ecae1", label="assembly")
        ax.bar(x, s, bottom=a, color=C_COUPLED, label="linear solve")
        ax.bar(x, t, bottom=np.add(a, s), color="#fdae6b", label="turbulence")
        ax.set_xticks(x)
        ax.set_xticklabels([b[0] for b in br])
        ax.set_ylabel("time per iteration [s]")
        ax.legend(fontsize=7)
        save(fig, "bench_breakdown",
             "coupledFoam time per iteration: assembly, linear solve, "
             "turbulence.")

    # Memory
    fig, ax = plt.subplots(figsize=(6.5, 2.6))
    for j, cfg in enumerate(cfgs):
        vals = [med(c, cfg, "peakRSS_GB_sum") or 0 for c in cases]
        ax.bar(np.arange(len(cases)) + j * 0.8 / len(cfgs), vals,
               0.8 / len(cfgs), color=CONFIG_COLOR[cfg], label=CONFIG_LABEL[cfg])
    ax.set_xticks(np.arange(len(cases)) + 0.4 - 0.4 / len(cfgs))
    ax.set_xticklabels(cases)
    ax.set_ylabel("peak RSS, sum over ranks [GB]")
    ax.legend(fontsize=7, ncol=2)
    save(fig, "bench_memory", "Peak memory (sum of the per-rank maximum RSS).")

    # Executive summary table
    header = ["case", "native best wall [s]", "coupled wall [s]",
              "speed-up (wall)", "native best CPU-h", "coupled CPU-h",
              "speed-up (CPU)", "memory ratio"]
    tab = []
    for c in cases:
        nat = [(med(c, k, "wall_to_conv_s"), k) for k in ("A", "B")
               if med(c, k, "wall_to_conv_s")]
        if not nat:
            continue
        wn, kbest = min(nat)
        wc = med(c, "C", "wall_to_conv_s")
        cn = med(c, kbest, "cpu_to_conv_h")
        cc = med(c, "C", "cpu_to_conv_h")
        mn = med(c, kbest, "peakRSS_GB_sum")
        mc = med(c, "C", "peakRSS_GB_sum")
        tab.append([c, fmt(wn), fmt(wc), fmt(wn / wc if wc else None),
                    fmt(cn), fmt(cc), fmt(cn / cc if cn and cc else None),
                    fmt(mc / mn if mn and mc else None)])
        num(f"speedup wall {c}", wn / wc if wc else None, "{:.2f}")
        num(f"speedup cpu {c}", cn / cc if cn and cc else None, "{:.2f}")
    return write_table("executive_summary", header, tab,
                       "Executive summary: time to convergence (identical "
                       "window criterion), best native configuration vs. "
                       "coupledFoam.", "tab:summary")


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
                     fmt(d.get("iterations"), "{}"), fmt(d.get("finalR")),
                     fmt(d.get("wallSecondsSolver") or d.get("wallSeconds")),
                     fmt(d.get("cpuHoursSolver")), key])
        num(f"test {name} pass", "yes" if d.get("pass") else "no")
    return write_table("tests", header, rows, "Test results (spec 13).",
                       "tab:tests")


# --------------------------------------------------------------------------- #
# amendment B (15.7) and remaining paper items
# --------------------------------------------------------------------------- #

def fig_eta_rho(tests: dict) -> None:
    """Eisenstat-Walker eta and preconditioner rho histories (15.7)."""
    cand = [(n, d) for n, d in tests.items()
            if d.get("history", {}).get("eta")]
    if not cand:
        notes.append("eta/rho history: no data")
        return
    fig, ax = plt.subplots(1, 2, figsize=(6.5, 2.6))
    for name, d in sorted(cand)[:4]:
        h = d["history"]
        it = np.arange(1, len(h["eta"]) + 1)
        ax[0].semilogy(it, h["eta"], lw=0.9, label=name)
        rho = [r if (r is not None and r >= 0) else np.nan for r in h.get("rho", [])]
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
    d = tests.get("Test-blockGAMG_cycles")
    if not d:
        notes.append("cycle comparison: no Test-blockGAMG_cycles result")
        return
    per = d.get("cycles") or d.get("perCycle") or {}
    if isinstance(per, list):
        per = {x.get("cycleType"): x for x in per}
    names = [c for c in ("V", "F", "W", "K") if c in per]
    if not names:
        notes.append("cycle comparison: unexpected JSON layout")
        return
    its = [per[c].get("nIterations") or 0 for c in names]
    wall = [per[c].get("wallSeconds") or per[c].get("solveWallSeconds") or 0
            for c in names]
    fig, ax = plt.subplots(1, 2, figsize=(6.5, 2.4))
    ax[0].bar(names, its, color=C_COUPLED)
    ax[0].set_ylabel("Krylov iterations")
    ax[1].bar(names, wall, color=C_COUPLED2)
    ax[1].set_ylabel("solve wall time [s]")
    save(fig, "cycle_comparison",
         "Block-GAMG cycle types on the block-Poisson system (T0 mesh).")


def fig_anderson(bench: dict) -> None:
    rows = [d for d in bench.values() if d.get("config") in ("C", "G")]
    cases = sorted({d["case"] for d in rows})
    if not rows or not any(d["config"] == "G" for d in rows):
        notes.append("Anderson comparison: no configuration G results")
        return

    def med(c, cfg, key):
        v = [d.get(key) for d in rows
             if d["case"] == c and d["config"] == cfg and d.get(key)]
        return float(np.median(v)) if v else 0.0

    fig, ax = plt.subplots(1, 2, figsize=(6.5, 2.4))
    x = np.arange(len(cases))
    for j, (cfg, col, lab) in enumerate((("C", C_COUPLED, "Anderson off"),
                                         ("G", "#d95f02", "Anderson on"))):
        ax[0].bar(x + 0.4 * j, [med(c, cfg, "iters_to_conv") for c in cases],
                  0.4, color=col, label=lab)
        ax[1].bar(x + 0.4 * j, [med(c, cfg, "wall_to_conv_s") for c in cases],
                  0.4, color=col, label=lab)
    for a in ax:
        a.set_xticks(x + 0.2)
        a.set_xticklabels(cases)
    ax[0].set_ylabel("iterations to convergence")
    ax[1].set_ylabel("wall time to convergence [s]")
    ax[0].legend(fontsize=7)
    save(fig, "anderson", "Anderson acceleration on and off (15.7).")


def fig_scaling(tests: dict) -> None:
    d = next((v for k, v in tests.items() if k.startswith("T-scaling")), None)
    pts = (d or {}).get("points") or (d or {}).get("scaling")
    if not pts:
        notes.append("strong scaling: no T-scaling result")
        return
    fig, ax = plt.subplots(1, 2, figsize=(6.5, 2.6))
    for solver, col in (("simpleFoam", C_NATIVE), ("coupledFoam", C_COUPLED)):
        p = sorted((q for q in pts if q.get("solver") == solver),
                   key=lambda q: q["nProcs"])
        if not p:
            continue
        n = np.array([q["nProcs"] for q in p], dtype=float)
        t = np.array([q["timePerIter"] for q in p], dtype=float)
        ax[0].loglog(n, t, "o-", color=col, label=solver)
        eff = t[0] * n[0] / (t * n)
        ax[1].plot(n, eff, "o-", color=col, label=solver)
    ax[0].set_xlabel("ranks")
    ax[0].set_ylabel("time per iteration [s]")
    ax[1].set_xlabel("ranks")
    ax[1].set_ylabel("parallel efficiency")
    ax[0].legend(fontsize=7)
    save(fig, "scaling", "Strong scaling (T-scaling on T4b).")


def fig_remediation(tests: dict) -> None:
    cand = [(n, d) for n, d in tests.items()
            if d.get("history", {}).get("nDyn")
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
              "tolerance", "pass"]
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
                     fmt(d.get(kd), "{:.2e}"), fmt(d.get(kt), "{:.3g}"),
                     "yes" if d.get("pass") else "NO"])
    for name in sorted(tests):
        if name.startswith(("T4", "T5")) and "Cd" in tests[name]:
            d = tests[name]
            rows.append([name, "Cd", fmt(d.get("Cd"), "{:.5g}"),
                         fmt(d.get("CdRef"), "{:.5g}"),
                         fmt(d.get("CdRelDiff"), "{:.2e}"), "0.01",
                         "yes" if d.get("pass") else "NO"])
    if not rows:
        notes.append("validation table: no T1-T5 results")
        return ""
    return write_table("validation", header, rows,
                       "Validation: integral quantities, coupledFoam vs. "
                       "simpleFoam on identical meshes and schemes.",
                       "tab:validation")


def table_gamg_levels(tests: dict) -> str:
    header = ["run", "levels", "mergeLevels", "C_op", "cells per level"]
    rows = []
    for name, d in sorted(tests.items()):
        g = d.get("gamg") or {}
        if g.get("gamgLevels"):
            rows.append([name, g.get("gamgLevels"), g.get("gamgMergeLevels"),
                         fmt(g.get("gamgCop"), "{:.3f}"),
                         " ".join(str(c) for c in (g.get("gamgCellsPerLevel") or []))])
    if not rows:
        notes.append("GAMG level table: no data")
        return ""
    return write_table("gamg_levels", header, rows,
                       "Block-GAMG hierarchy per run.", "tab:gamglevels")


# --------------------------------------------------------------------------- #
# main
# --------------------------------------------------------------------------- #

def main() -> int:
    tests = load_json("tests")
    bench = load_json("bench")

    fig_T0_profiles()
    fig_histories(tests)
    summary_md = fig_bench(bench)
    tests_md = tests_table(tests)
    fig_eta_rho(tests)
    fig_cycles(tests)
    fig_anderson(bench)
    fig_scaling(tests)
    fig_remediation(tests)
    validation_md = table_validation(tests)
    levels_md = table_gamg_levels(tests)

    t0 = tests.get("T0_Re100_np1", {})
    num("T0 Re100 iterations", t0.get("iterations"), "{}")
    num("T0 Re100 l2u", t0.get("l2rel_u"), "{:.1e}")
    num("T0 Re100 l2v", t0.get("l2rel_v"), "{:.1e}")
    num("T0 Re100 native iterations",
        (t0.get("reference") or {}).get("iterations"), "{}")
    num("commit", git_commit())

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
        "## Figures",
        "",
    ]
    for name, cap in figures:
        md += [f"![{name}](figures/{name}.png)", "", f"*{cap}*", ""]
    md += ["## Limitations and Phase 2", "",
           "See `DECISIONS.md` and the paper, Section 6.", "",
           "## Reproduction", "",
           "```", "source <openfoam2606>/etc/bashrc", "./Allwmake -j 8",
           "~/OF/venv/bin/pytest tests/", "bench/run_bench.py",
           "bench/make_report.py", "cd report/paper && latexmk -pdf paper.tex",
           "```", ""]
    if notes:
        md += ["## Generator notes", ""] + [f"- {n}" for n in notes] + [""]
    (REPORT / "REPORT.md").write_text("\n".join(md))
    print(f"figures: {len(figures)}, numbers: {len(numbers)}, notes: {len(notes)}")
    for n in notes:
        print("note:", n)
    return 0


if __name__ == "__main__":
    sys.exit(main())

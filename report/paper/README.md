# coupledFoam paper

`paper.tex` is the LaTeX paper on coupledFoam. It never contains
hand-copied results: every number, table and figure of the verification and
performance sections is generated from the JSON records in `results/` by
`bench/make_report.py`.

## Two versions (D-032)

Two documents are built from the same generated data:

| Source | PDF | Reader |
|---|---|---|
| `paper.tex` | `paper.pdf` | CFD numerics specialists: the professional paper |
| `paper_tutorial.tex` | `paper_tutorial.pdf` | engineers without a numerics background: every concept explained from intuition, with "Intuition", "Why it matters", "Pitfall" and "Building block" boxes, TikZ sketches, a guide to reading each figure and table, and a glossary |

Both `\input{numbers.tex}` and include the same `tables/*.tex` and
`figures/*.pdf` through the same `\cffigure`/`\cftable`/`\cfnum`/`\cfsci`
helpers (defined identically in both preambles), so neither can drift from
the measured data. Neither file contains a hand-typed measured number. When
a technical statement changes in `paper.tex`, check the corresponding
section of `paper_tutorial.tex` as well. The tutorial additionally needs the
LaTeX packages `tcolorbox` (libraries `breakable`, `skins`), `tikz` and
`array`.

## Files

| File | Written by | Content |
|---|---|---|
| `paper.tex` | hand | text, equations, algorithms |
| `references.bib` | hand | bibliography (BibTeX, `plainnat`) |
| `numbers.tex` | `bench/make_report.py` | `\newcommand{\cf...}{value}` for every quoted number |
| `tables/*.tex` | `bench/make_report.py` | booktabs tables (`tests.tex`, `executive_summary.tex`, ...) |
| `figures/*.pdf` | `bench/make_report.py` | vector figures |

Every generated file is optional. `paper.tex` loads each one through
`\IfFileExists`, so the paper compiles at every stage of the project; a
missing table or figure is shown as a boxed "Results pending" placeholder
and a missing number as a small boxed "pending" mark.

## Regenerate numbers, tables and figures

From the repository root (no OpenFOAM environment needed; Python with
numpy and matplotlib, e.g. the project venv):

```
bench/make_report.py            # or: ~/OF/venv/bin/python bench/make_report.py
```

or from this directory `make numbers PYTHON=~/OF/venv/bin/python`, which
also rebuilds the PDF.

## Build the PDF

Requires a TeX distribution with `latexmk`, `pdflatex` and `bibtex`
(on Ubuntu: `texlive-latex-recommended texlive-latex-extra
texlive-science latexmk`). Packages used: amsmath, amssymb, graphicx,
booktabs, natbib, hyperref, geometry, fontenc; siunitx and lmodern are used
if installed and replaced by simple fallbacks otherwise.

```
cd report/paper
make            # builds paper.pdf and paper_tutorial.pdf (latexmk -pdf)
make paper_tutorial.pdf   # only one of them
make clean      # remove intermediate files of both
make distclean  # also remove both PDFs
```

## Contract with the report generator

Numbers are referenced with `\cfnum{Name}` (plain) or `\cfsci{Name}`
(through siunitx `\num`, e.g. `2.9e-06`), where `Name` is the generated
macro without the `cf` prefix, e.g. `\cfnum{TZeroReOneZeroZeroIterations}`
for `\cfTZeroReOneZeroZeroIterations`. Macro names follow `num()` in
`bench/make_report.py` (digits spelled out, words capitalised).

Generated files the paper includes (stem, label):

| Kind | Stem | Label | Status of generator |
|---|---|---|---|
| table | `tests` | `tab:tests` | generated |
| table | `executive_summary` | `tab:summary` | generated once benchmark JSON exist |
| table | `validation` | `tab:validation` (must be set in the file) | to be added: T1-T5 integral quantities vs. simpleFoam |
| table | `gamg_levels` | `tab:gamglevels` (must be set in the file) | to be added: levels, cells/ranks per level, ratios, C_op |
| figure | `T0_Re100_profiles`, `T0_Re1000_profiles` | set by the paper | generated when the run data exist |
| figure | `T0_Re100_np1_history` | set by the paper | generated |
| figure | `bench_wall`, `bench_cpu`, `bench_iters`, `bench_breakdown`, `bench_memory` | set by the paper | generated once benchmark JSON exist |
| figure | `scaling` | set by the paper | to be added (T-scaling) |
| figure | `cycle_comparison` | set by the paper | to be added (configuration E) |
| figure | `eta_rho_history` | set by the paper | to be added (eta and rho histories, tune events) |
| figure | `anderson` | set by the paper | to be added (Anderson on/off) |
| figure | `remediation_history` | set by the paper | to be added (T4/T5 set sizes, rollbacks) |

A generated table carries its own `\caption` and `\label`; the paper uses
the label given above. Figures get caption and label from the paper.

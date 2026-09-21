# coupledFoam paper - build with latexmk (GNU make).
#
#   make           build paper.pdf
#   make numbers   regenerate numbers.tex, tables/*.tex, figures/*.pdf
#                  from results/ (runs bench/make_report.py), then build
#   make clean     remove LaTeX intermediate files (keeps paper.pdf)
#   make distclean also remove paper.pdf
#
# Recipes use '>' as prefix instead of a TAB (GNU make >= 3.82).
.RECIPEPREFIX = >

LATEXMK ?= latexmk
PYTHON  ?= python3
REPO    := ../..

GENERATED := $(wildcard numbers.tex tables/*.tex figures/*.pdf)

.PHONY: all numbers clean distclean

all: paper.pdf

paper.pdf: paper.tex references.bib $(GENERATED)
> $(LATEXMK) -pdf -interaction=nonstopmode -halt-on-error paper.tex

numbers:
> cd $(REPO) && $(PYTHON) bench/make_report.py
> $(MAKE) all

clean:
> $(LATEXMK) -c paper.tex
> rm -f paper.bbl paper.run.xml paper.bcf

distclean: clean
> $(LATEXMK) -C paper.tex
> rm -f paper.pdf

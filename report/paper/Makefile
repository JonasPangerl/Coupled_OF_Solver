# coupledFoam paper - build with latexmk (GNU make).
#
# Two documents are built from the same generated data (D-032):
#   paper.pdf           the professional paper
#   paper_tutorial.pdf  the tutorial (plain-language) version
#
#   make           build paper.pdf and paper_tutorial.pdf
#   make paper.pdf, make paper_tutorial.pdf
#                  build only one of them
#   make numbers   regenerate numbers.tex, tables/*.tex, figures/*.pdf
#                  from results/ (runs bench/make_report.py), then build both
#   make clean     remove LaTeX intermediate files (keeps the PDFs)
#   make distclean also remove the PDFs
#
# Recipes use '>' as prefix instead of a TAB (GNU make >= 3.82).
.RECIPEPREFIX = >

LATEXMK ?= latexmk
PYTHON  ?= python3
REPO    := ../..

DOCS      := paper paper_tutorial
PDFS      := $(addsuffix .pdf,$(DOCS))
GENERATED := $(wildcard numbers.tex tables/*.tex figures/*.pdf)

.PHONY: all numbers clean distclean

all: $(PDFS)

paper.pdf: paper.tex references.bib $(GENERATED)
> $(LATEXMK) -pdf -interaction=nonstopmode -halt-on-error paper.tex

paper_tutorial.pdf: paper_tutorial.tex references.bib $(GENERATED)
> $(LATEXMK) -pdf -interaction=nonstopmode -halt-on-error paper_tutorial.tex

numbers:
> cd $(REPO) && $(PYTHON) bench/make_report.py
> $(MAKE) all

clean:
> for d in $(DOCS); do $(LATEXMK) -c $$d.tex; rm -f $$d.bbl $$d.run.xml $$d.bcf; done

distclean: clean
> for d in $(DOCS); do $(LATEXMK) -C $$d.tex; done
> rm -f $(PDFS)

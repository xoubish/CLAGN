.PHONY: all figures clean
PYTHON ?= python
export PYTHONDONTWRITEBYTECODE = 1

all:
	latexmk -pdf -interaction=nonstopmode -halt-on-error proposal.tex

figures:
	$(PYTHON) make_timing_audit.py
	$(PYTHON) make_selection.py
	$(PYTHON) make_fig1_connected.py

clean:
	latexmk -c proposal.tex

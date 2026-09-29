.PHONY: all figures clean
PYTHON ?= python
export PYTHONDONTWRITEBYTECODE = 1

all:
	latexmk -pdf -interaction=nonstopmode -halt-on-error proposal_reframed.tex

figures:
	$(PYTHON) make_selection.py
	$(PYTHON) make_fig1_five.py stack
	$(PYTHON) make_fig2_memory.py
	$(PYTHON) make_fig2_forecast.py

clean:
	latexmk -c proposal_reframed.tex

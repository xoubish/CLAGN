.PHONY: all figures single-example-preview bundle clean

PYTHON ?= python

all:
	latexmk -pdf -interaction=nonstopmode -halt-on-error proposal.tex

figures:
	$(PYTHON) make_fig1_histories.py
	$(PYTHON) make_fig2_model.py

single-example-preview:
	$(PYTHON) make_fig1_single_example.py
	$(PYTHON) make_fig1_single_example.py --target P2190

bundle: all
	zip -r clagn_jwst_cycle6_draft.zip proposal.tex proposal.pdf jwstproposaltemplate_v6.sty JWST_proposal_template_cy6.tex apt_title_abstract.txt fig1_histories.pdf fig1_histories.png fig2_model.pdf fig2_model.png make_fig1_histories.py make_fig2_model.py inputs README.md requirements.txt Makefile .gitignore

clean:
	latexmk -c proposal.tex

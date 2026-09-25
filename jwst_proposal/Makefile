.PHONY: all figures single-example-preview bundle clean

PYTHON ?= python
TECTONIC ?= tectonic

all:
	@if command -v latexmk >/dev/null 2>&1; then \
		latexmk -pdf -interaction=nonstopmode -halt-on-error proposal.tex; \
	else \
		$(TECTONIC) --untrusted --keep-logs proposal.tex; \
	fi

figures:
	$(PYTHON) make_fig1_single_example.py
	$(PYTHON) make_fig2_model.py

single-example-preview:
	$(PYTHON) make_fig1_single_example.py
	$(PYTHON) make_fig1_single_example.py --target P2190
	$(PYTHON) make_fig1_single_example.py --target P1823

bundle: all
	zip -r clagn_jwst_cycle6_draft.zip proposal.tex proposal.pdf jwstproposaltemplate_v6.sty JWST_proposal_template_cy6.tex apt_title_abstract.txt fig1_single_example_preview.pdf fig1_single_example_preview.png fig2_model.pdf fig2_model.png make_fig1_single_example.py spherex_line_labels.py make_fig2_model.py spaxel_scryer_*_cleaned.csv inputs README.md revision_notes_v0.4.md requirements.txt Makefile .gitignore

clean:
	rm -f proposal.aux proposal.log proposal.out proposal.fdb_latexmk proposal.fls proposal.synctex.gz

# JWST Cycle 6 proposal

This is the active, Git-visible proposal folder. Copy this directory into another repository or open it in another checkout; building the proposal and regenerating its figures require no files from the ignored `observing/` or `data/` trees.

The current scientific text and figure designs are v0.3. This relocation preserves them. The previous working copies, historical drafts and review notes remain in the original ignored `observing/jwst_cycle6_20260925/proposal_draft/` directory.

## Edit and build

Edit `proposal.tex`, then run from this directory:

```sh
latexmk -pdf -interaction=nonstopmode -halt-on-error proposal.tex
```

Alternatively, run `make`. The existing figure PDFs are included, so Python is unnecessary for ordinary text editing and compilation. A TeX distribution with `latexmk`, `pdflatex`, `mathptmx`, and `graphicx` is sufficient.

## Files

| File | Purpose |
| --- | --- |
| `proposal.tex` / `proposal.pdf` | Main proposal and compiled preview |
| `apt_title_abstract.txt` | Separate title and abstract for APT |
| `jwstproposaltemplate_v6.sty` | Unmodified official Cycle 6 style |
| `JWST_proposal_template_cy6.tex` | Original official template for reference |
| `fig1_histories.pdf` / `.png` | Archival/NGPS spectra, ZTF and WISE histories |
| `fig2_model.pdf` / `.png` | Illustrative dust-model comparison |
| `make_fig1_histories.py`, `make_fig2_model.py` | Figure generators |
| `inputs/figure_targets.json` | Only the three targets and data fields used by these generators |
| `inputs/provenance.json` | Snapshot origin and checksums |
| `inputs/allwise_w1w4_candidates.csv` | Existing AllWISE planning table |

The figure snapshot preserves the local collaboration epochs already used in v0.3; it does not switch to the smaller public-page spectral set. The AllWISE table is copied unchanged; its error-column units need checking before use in calculations.

## Regenerate figures

```sh
python -m pip install -r requirements.txt
python make_fig1_histories.py
python make_fig2_model.py
```

Or run `make figures`. Both scripts load `inputs/figure_targets.json` relative to their own location. They do not download data or access the original workspace. Times New Roman/Times fonts reproduce the original figure appearance when installed.

## Git and Overleaf

From the repository root, include `jwst_proposal/`, `.gitignore`, and the root `README.md` in your commit. Proposal/figure PDFs are explicitly allowed by the repository's ignore rules; compilation scratch files and ZIPs are ignored.

For Overleaf, run `make bundle`, upload `clagn_jwst_cycle6_draft.zip`, and select `proposal.tex` as the main document. The folder can also be copied directly to another repository; check that repository's PDF ignore rules.

## Working status

This is an editable proposal draft, not a submitted program. The relocation does not validate the model assumptions, figure typography, target states, ETC sensitivities, or APT charged times. Those scientific and submission checks remain part of proposal development. Keep the official style unchanged and update the AI disclosure as the proposal is revised.

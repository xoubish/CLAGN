# Current Figure 2: spectrum, space, and time

[Preview](fig2_three_panel_draft.png) · [Overleaf figure PDF](fig2_three_panel_draft.pdf) · [Caption](caption.tex)

This is the current two-row design: a freshly plotted vector spectrum spans
the top row, with spatial and temporal examples side by side below it. The
lower panels retain their source aspect ratios and are about 1.5 times larger
at proposal width than in the original three-panel row.
Titles use descriptive names, avoiding confusion with MRS sub-band letters.
The PDF is installed as `../../fig2_miri.pdf`, and the caption is mirrored in
`../../proposal.tex`. The previous mock figure and generator are archived in
`../../older/before_literature_figure2_20260929.zip`.
Do not compile the proposal PDF during routine figure edits.

- **Top — Separating hot and warm dust:** Donnan et al. (2024), Figure 8,
  NGC 7469, digitized from the native 4026 × 3355 image embedded in the paper.
  Five wavelength ticks and three flux ticks calibrate the logarithmic axes
  to maximum residuals of 0.12 and 0.33 pixels. The extracted traces are saved
  in `../../inputs/fig2_ngc7469_digitized.csv`, with provenance in the adjacent
  `fig2_ngc7469_digitization.json`. The observed trace uses column-median black
  symbol positions; these are display values, not original spectral samples.
  Hot/warm curves join dash gaps and smooth pixel steps; no extrapolation or
  new fit is made. The full model and negligible star-forming dust component
  are omitted. Dark, unconnected dots display every second extracted spectrum
  column. Shading marks the 8–13 µm integration band used to measure dust
  luminosity after PAH/line subtraction, as explained in the caption.
  SPHEREx and JWST MIRI/MRS headings and coverage bars are aligned; four
  colours identify the MRS channels without individual channel labels.
  The bars use NGC 7469's redshift of 0.0163, with observed instrument limits
  divided by 1 + z to match the rest-wavelength axis. The plotted spectrum is
  from NIRSpec IFU and MIRI MRS (Donnan et al., Section 4); the SPHEREx bar
  shows comparison access only. The observed channel limits follow the STScI
  MRS documentation linked in the provenance. The galaxy label includes its
  redshift, and the spectrum legend names NIRSpec + MRS.
- **Lower left — Isolating the central emission:** González-Martín et al. (2025), Figure 4,
  NGC 6552 (downloaded 2026 revision), with spectra and three 11.6-micron maps.
- **Lower right — Reading the dust’s memory:** Lyu & Rieke (2021), Figure 15, NGC 4151.
  The year axis is moved below the curves. Complete 20–24 and 34–37 micron
  photometric series and labels are omitted. The original fit-parameter labels
  are replaced by a fitted delay of approximately 8 years and a smoothing window
  of approximately 14 years. These are literature-model values, not predictions
  or directly measurable lags for our single-visit sample.

The caption records these adaptations and connects the examples to our experiment.
The top panel uses rest wavelengths; the lower-left panel retains observed wavelengths.

The top panel retains vector axes, curves and searchable labels. The lower
source figure regions are rendered at 600 dpi and embedded as losslessly
compressed crops, using `../pdf_figure_crop.py`. Full source pages, off-crop
text, and tables are not embedded.

The lower-left panel uses [arXiv v2](https://arxiv.org/abs/2504.01103v2), submitted in April
2026, of the [2025 paper](https://doi.org/10.1093/mnras/staf573). The 2026
preprint header does not identify a different publication. The published
[2026 correction](https://doi.org/10.1093/mnras/stag1045) concerns an omitted
acknowledgement. Reference [28] records both the publication and revised source.

## Reproduction

From the repository root, with the saved literature PDFs available:

```sh
/opt/anaconda3/bin/python jwst_proposal/make_fig2_miri.py
```

This runs `../make_fig2_spectrum.py` and `../make_fig2_three_panel.py`
sequentially and installs the final PDF. `panel_a/` contains the vector spectrum;
`build/` contains intermediate PDF layers. Both are ignored by Git. Final PDF,
PNG, caption, scripts, and checks remain eligible for version control.
Source papers and download manifests are in `../miri_literature/`.
The component generators verify that the proposal source and installed figure
remain unchanged while rendering; the entry point then installs the result.
It does not rebuild the proposal PDF or overwrite caption edits in the source.

Previous two-panel, six-panel, spatial-pair and stacked-layout attempts are
preserved in `../../older/figure2_trials_20260929.zip`; original paths and verified
checksums are in `../cleanup_20260929.json`. The current generator no longer
requires an archived spatial-pair PDF.

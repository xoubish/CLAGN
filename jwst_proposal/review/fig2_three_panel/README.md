# Current Figure 2: spectrum, space, and time

[Preview](fig2_three_panel_draft.png) · [Overleaf figure PDF](fig2_three_panel_draft.pdf) · [Caption](caption.tex)

This is the accepted three-panel design. The panels have equal content heights,
proportional widths, centered titles, and galaxy names and citations inside.
The PDF is installed as `../../fig2_miri.pdf`, and the caption is mirrored in
`../../proposal.tex`. The previous mock figure and generator are archived in
`../../older/before_literature_figure2_20260929.zip`.
Do not compile the proposal PDF during routine figure edits.

- **A — Separating hot and warm dust:** Donnan et al. (2024), Figure 8,
  NGC 7469, with nominal SPHEREx/MRS coverage, PAH positions and silicate markers.
- **B — Isolating the central emission:** González-Martín et al. (2025), Figure 4,
  NGC 6552 (downloaded 2026 revision), with spectra and three 11.6-micron maps.
- **C — Reading the dust’s memory:** Lyu & Rieke (2021), Figure 15, NGC 4151.
  The year axis is moved below the curves. Complete 20–24 and 34–37 micron
  photometric series and labels are omitted. The original fit-parameter labels
  are replaced by a fitted delay of approximately 8 years and a smoothing window
  of approximately 14 years. These are literature-model values, not predictions
  or directly measurable lags for our single-visit sample.

The caption records these adaptations and connects the examples to our experiment.
Panel A uses rest wavelengths; panel B retains observed wavelengths.

Source figure regions are rendered at 600 dpi and embedded as losslessly
compressed cropped images before assembly. Full source pages, off-crop text,
and tables are not embedded. Final titles and panel B/C identifications remain
searchable text. This is implemented in `../pdf_figure_crop.py`.

Panel B uses [arXiv v2](https://arxiv.org/abs/2504.01103v2), submitted in April
2026, of the [2025 paper](https://doi.org/10.1093/mnras/staf573). The 2026
preprint header does not identify a different publication. The published
[2026 correction](https://doi.org/10.1093/mnras/stag1045) concerns an omitted
acknowledgement. Reference [28] records both the publication and revised source.

## Reproduction

From the repository root, with the saved literature PDFs available:

```sh
/opt/anaconda3/bin/python jwst_proposal/make_fig2_miri.py
```

This runs both component generators sequentially and installs the final PDF. `panel_a/` contains the generated panel-A dependency;
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

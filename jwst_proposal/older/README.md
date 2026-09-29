# Superseded proposal material

The active proposal source is [../proposal.tex](../proposal.tex); no compiled
proposal PDF is maintained.

- `before_literature_figure2_20260929.zip`: previous installed mock Figure 2,
  its generator and provenance, and the proposal source and README before
  installing the accepted three-panel literature figure. Also retains the
  previous validation metadata.

- `figure2_trials_20260929.zip`: local-only archive of the superseded two-panel,
  six-panel, spatial-pair and stacked layouts, their trial generators, and the
  previous mock-figure preview. Original paths are relative to `jwst_proposal/`;
  extract into a separate directory for review. All archive members were verified
  against SHA-256 hashes before removing the loose files. See
  [the manifest](../review/cleanup_20260929.json). The accepted current design is
  [here](../review/fig2_three_panel/README.md).

- `before_readable_checkpoint/`: the workflow figure and proposal before returning
  to the measured R01 spectrum plus conditional MIRI scenarios; includes both
  figure PDFs, generator, style file and prior validation.
- `before_measurement_workflow/`: the last spectrum-plus-model-comparison
  figure, source, PDF, generator, numerical inputs and metadata before replacing
  it with the measurement-to-science workflow.
- `before_spectrum_comparison_cleanup/`: the spectrum-led version with both
  model spectra overlaid, before moving the comparison entirely to panel C.
- `before_spectrum_figure/`: the two-panel combined figure and proposal before
  adding the dominant model-spectrum panel; includes its generator and curves.
- `before_combined_figure/`: the last reframed draft, both former Figures 2–3,
  Figure 1, style file, build instructions, and generator snapshots before the
  rename and combined science figure. The saved source compiles in that folder.
- `previous_main/`: the incoming `proposal.tex` and PDF, with the figures and
  style file needed to compile that saved draft.
- `reframed_before_revision/`: the incoming reframed draft, its figures,
  relevant generator snapshots, forecast tables and prior validation metadata.
  The saved LaTeX can be compiled within that folder. Generator copies are
  historical records; their data paths refer to the original working folder.
- `fig1_connected.pdf`: previous single-target P9694 figure.
- `fig3_dust_evolution.pdf`: previous idealized boundary-evolution figure.
- `fig1_stack_preview.pdf` / `.png`: incoming three-target layout preview.
- `reframe.diff`: the earlier narrative comparison.

Shared numerical code, observations and input data remain in the working folder.
Nothing in this directory is part of the active build.

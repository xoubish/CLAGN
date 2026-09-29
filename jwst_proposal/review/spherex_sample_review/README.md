# Full-sample exploratory review

Open `index.html` to inspect all 24 targets and the three candidate sample plots.
The atlas contains four pages, six targets per page. Figure PDFs are standalone
review artifacts; no proposal PDF is built and no proposal input is replaced.

## What the first pass shows

- Most targets have useful measurements in two or three time groups, but some
  groups cover only a small wavelength interval. A time group is not necessarily
  a complete spectrum.
- The history-versus-color diagram shows a broad spread of continuum shapes at
  similar W1 history ratios. It does not establish a clean separation or a
  significant correlation. These are total-light colors, not hot-dust temperatures.
- The visit-change diagram is worth inspecting alongside the atlas. F05, for
  example, has a larger long-wavelength difference between its two visits than
  several other targets. This is a candidate change, pending extraction and
  calibration checks; error bars alone do not establish intrinsic variability.
- R10 and R11 lack adequate coverage in the chosen short continuum window in
  the initial snapshot. R12 has no accepted extraction. All remain visible in
  the inventory, coverage plot and atlas.
- Three rejected R12 cutouts were inspected: two have only 13 and 0 usable
  background pixels (the extractor requires 15), despite acceptable central
  pixels. A third has the source outside the permitted cutout margin. These
  examples explain specific extraction failures, not every R12 rejection;
  they are not evidence of source nondetection. See `rejection_checks/diagnostics.json`.
- Some downloader results labelled accepted contain non-finite photometry.
  The atlas excludes those values and records their counts as `invalid_accepted`.
  Earlier download-status counts therefore exceed the finite points available
  for scientific use. The underlying per-image records are preserved.

## Definitions and limits

Visits are split at observation-time gaps exceeding 45 days. Colors compare
rest-frame 3.6–3.9 µm and 2.1–2.4 µm continuum windows. Each window has three
fixed wavelength cells; each cell must contain at least one valid point. We
take the inverse-variance mean inside each cell, then average the three cells
equally. Both window means must exceed three times their approximate errors
to enter a ratio plot. All other targets and visits remain in the CSVs and atlas.
No interpolation fills missing cells, and no selection cuts are made on
individual positive flux or individual signal-to-noise.

Ratio errors use independent propagation of supplied errors. The supplied R01
uncertainties and new extraction uncertainties have not been homogenized, and
correlated calibration errors have not been propagated. The two coordinates
of the visit-change diagram share a flux measurement and are correlated.
The plotted changes use earliest and latest **eligible** visits, not necessarily
the same dates across targets. They are exploratory diagnostics, not a fitted
temperature distribution, dust-response recovery, or variability significance test.

The WISE history values come from the existing
`review/data_anchored_scenarios/sample_measurements.csv`. The horizontal
history statistic is median flux of the last three visits divided by the first
three; amplitude is the 95th-to-5th-percentile ratio. Neither specifies a unique
history class. Inspect the actual light curves before assigning flare,
recovery, sustained rise or decline labels.

Only a subset of archival/NGPS spectra is currently indexed locally by the
existing scripts. The coverage figure therefore shows WISE and SPHEREx only;
this local inventory does not establish that optical spectra are unavailable
elsewhere or that scheduled spectra have already been obtained.

## Refresh after further downloads

From the project root:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 /opt/anaconda3/bin/python jwst_proposal/review/review_spherex_sample.py
```

This reads saved per-image JSON results directly, so it can be run while the
background downloader is active. It records the snapshot date, each input
checksum, per-target completeness and exclusions. `*_snapshot.csv` preserve
the exact finite measurements used. `summary.json`, `inventory.csv`,
`visit_measurements.csv`, and `visit_changes.csv` provide machine-readable results.

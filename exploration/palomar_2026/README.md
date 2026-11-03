# Palomar 2026 discovery exploration

**Start with [the executed census notebook](notebooks/00_start_here_census_and_quality.ipynb), or [its browser-readable HTML](outputs/html/00_start_here_census_and_quality.html).**

This is the first reproducible exploration of the completed three-night campaign,
using the Palomar science notes and current JWST proposal, as agreed. It includes
executed notebooks, all-source plots, measurement tables and a concrete science
plan. Snapshot date: **2026-11-03**. Original reductions, caches, proposal and
observing products are read only.

## What is here

| Notebook | Question and products |
|---|---|
| [00 — census and quality](notebooks/00_start_here_census_and_quality.ipynb) · [HTML](outputs/html/00_start_here_census_and_quality.html) | What was observed, what histories exist, and which measurements are usable? |
| [01 — lines and continuum](notebooks/01_emission_lines_and_continuum.ipynb) · [HTML](outputs/html/01_emission_lines_and_continuum.html) | Hβ, Hα and Mg II wing indices, continuum slopes, archival comparisons, five repeat sources, and calibration/continuum sensitivity. |
| [02 — photometry and timing](notebooks/02_photometry_timing_and_cross_band.ipynb) · [HTML](outputs/html/02_photometry_timing_and_cross_band.html) | Optical/IR histories, actual temporal overlap, amplitude diagnostics and matched g−r colors. |
| [03 — individual sources](notebooks/03_individual_sources_and_atlas.ipynb) · [HTML](outputs/html/03_individual_sources_and_atlas.html) | Editable source inspection, eight highlighted cases, and the complete 99-source atlas. |
| [04 — population experiments](notebooks/04_population_and_selection_experiments.ipynb) · [HTML](outputs/html/04_population_and_selection_experiments.html) | Frozen selection groups, measurable denominators, descriptive source-level bootstrap, and requirements for a valid enrichment test. |
| [05 — proposal and next steps](notebooks/05_proposal_mapping_and_next_measurements.ipynb) · [HTML](outputs/html/05_proposal_mapping_and_next_measurements.html) | Actual JWST coverage, proposal-to-measurement map, and per-source work queue. |

**All targets:** [99-page PDF atlas](outputs/all_sources_atlas.pdf) ·
[individual-source index](outputs/source_index.html) ·
[source inventory](outputs/objects.csv).

**Plan and implementation:** [science plan](SCIENCE_PLAN.md) ·
[measurement code](analysis.py) · [plotting code](plots.py) ·
[notebook source builder](build_notebooks.py) ·
[measurement validation](test_measurements.py).

## Findings from this pass

There are **104 delivered visits to 99 distinct targets**: 18 September visits,
44 October 26 visits including the 16 Cyg stellar observation, and 42 October 27
visits. Five sources were repeated: P7837, P8544, P9506, P9694 and P10381.
The census separates 95 AGN-program sources, two stellar/special targets and two
ancillary nearby galaxies; these are working analysis strata, not new classifications.

**85 targets have retained archival spectra.** The measurement catalog contains
365 epochs: 104 NGPS, 213 SDSS with recovered native pixel errors, and 48 DESI
spectra whose current compact representation lacks propagated errors. Twenty-six
retained epochs are flagged as collaboration data. All rejected and duplicate
cache records remain traceable in the audit. Eighty sources have at least one
usable Hβ comparison under the conservative coverage screen.

The most useful first inspections are:

| Source | Measured exploratory result | Interpretation to test |
|---|---|---|
| **P21422** | Hβ wing EW is 0.54–0.62 of the 2022 archive across three continuum choices; Hα is about 0.67–0.68. Hβ wing/[O III] is about 0.51 of the reference. | A coherent weakening lead with fading optical/IR history. Broad lines remain visible. Literature status and calibration need adjudication. |
| **P24930** | Hβ wing EW ratio 0.58–0.61 relative to the 2003 reference; wing/[O III] ratio about 0.66. | Long-baseline weakening; the long interval leaves the transition date poorly constrained. |
| **P7947** | Hβ wing EW ratio 1.83–2.31 relative to the 2001 reference; wing/[O III] ratio about 1.60. | Long-baseline strengthening; confirm broad-profile/host separation and check novelty. |
| **P17881** | Hβ wing EW ratio 2.17–3.27 relative to February 2025; Hα 1.58–1.65. Hβ wing/[O III] increases only about 1.32×. | Known-CLAGN/JWST follow-up, with noisy broad Hβ and a substantial normalization dependence. Prioritize joint fitting over a simple luminosity claim. |
| **P9694** | October Hβ wing EW is about 1.40× the 2021 reference, but October/September is only 1.00–1.03. | Persistence of the September strengthening lead; no large new month-scale EW change established. |
| **P8544** | October/September Hβ EW ratio ranges from 0.99 to 1.17 with continuum choice; narrow-line normalization differs. | The simple interpretation of continuing weakening is not established by these repeats. |
| **P10381** | September fails the primary 98%-coverage Hβ cut. On explicitly shared pixels (~92% coverage), October/September is 1.17–1.19. | A partial-coverage sensitivity result worth fitting, not a primary complete-window measurement or formal significance. |

Ranges above span continuum choices; **they are not confidence intervals**.
The P10381 shared-support experiment uses resampling and a different support
definition; it is not mixed into the primary catalog.

Eight additional Mg II comparisons trigger the adopted large-change screen:
P10252, P11010, P16145, P16218, P21841, P2921, P6658 and P7344. Fe II and
pseudo-continuum sensitivity make these fitting targets, not demonstrated Mg II
transitions. Across all lines, 12 distinct sources trigger this exploratory
threshold. Neither this count nor its fraction is a confirmed discovery yield.

P11530 and P1823 remain useful time-domain cases even though their new Hβ indices
are close to the latest accepted archival measurements. The atlas also retains
H5MUSES-327's stellar optical spectrum and the ancillary targets as separate
discovery opportunities.

## The most consequential gaps

* **Current optical photometry:** 97 sources have an optical cache, but the strict
  ±30-day match finds only **one Palomar visit** with accepted unbinned optical
  data: P10381 in September, in both g and r. None of the October visits passes.
  These are cache-availability statements, not absence of observations or variability.
  Median survey-photometry gaps are about one year; the IR histories are older.
* **Joint spectral fitting:** indices do not separate host, Fe II, broad/narrow
  lines or instrumental resolution. Slit/cloud losses and color calibration are
  not included in the statistical errors. Widths, black-hole masses, intrinsic
  responsivities and CL classifications have deliberately not been inferred.
* **JWST completeness:** **23 of 24** proposal targets have a delivered Palomar
  spectrum in this snapshot. **P24976 is missing.** P11530 has its September
  spectrum. Proposed MIRI measurements are future observations, not data here.
* **Selection inference:** actual observations, known-source status, line
  detectability, archival intervals and selection decisions must enter the
  population analysis. A targeted follow-up fraction is not an AGN incidence rate.

The [work queue](outputs/work_queue.csv) names tasks for every source. The
[science plan](SCIENCE_PLAN.md) gives the order, deliverables and acceptance criteria.

## Measurement contracts

`line_indices.csv` contains signed rest-frame wing pseudo-EWs in Å and integrated
wing flux proxies in **10⁻¹⁷ erg s⁻¹ cm⁻²**. `continuum` and the continuum-window
levels are observed fλ evaluated at rest-defined windows, in
**10⁻¹⁷ erg s⁻¹ cm⁻² Å⁻¹**; they are not luminosities. The pseudo-continuum slope
uses fλ ∝ λα, includes host light and is not an AGN-only spectral index.

The line windows, continuum variants and exact integration implementation are
visible in notebook 01. Missing errors stay missing. Error columns explicitly
refer to diagonal statistical propagation and omit covariance, sensitivity-fit,
slit, cloud, host, Fe II and residual telluric uncertainties. Alternative-standard
spectra are systematic comparisons of the same observation, not new epochs.

`comparisons.csv` distinguishes latest archive, first archive and previous NGPS.
The latest reference is the latest **measurable** accepted epoch for that line,
so it may predate a newer spectrum with incomplete coverage. There is no automatic
state label. Ratio screens require both wing EWs >2 Å, full coverage, and the
same direction across three sideband choices. This is a review heuristic, not
a line-detection significance or a transition definition.

`photometry.csv` retains band, magnitude system, stream, source and object ID.
NEOWISE/unWISE and ZTF survey/corrected-alert data are never silently merged.
Photometric plots include total source light, and display-binned data are not
used to claim precise ±30-day matches. The census counts retained spectral
epochs/intervals, not statistically independent transitions.

## Reproduce

The existing `/opt/anaconda3/bin/python` environment has the required packages.
From the project root:

```sh
/opt/anaconda3/bin/python exploration/palomar_2026/analysis.py
/opt/anaconda3/bin/python -m unittest discover -s exploration/palomar_2026 -p 'test_*.py' -v
/opt/anaconda3/bin/python exploration/palomar_2026/build_notebooks.py
/opt/anaconda3/bin/python exploration/palomar_2026/execute_notebooks.py
```

Jupyter needs local loopback sockets. The sandboxed first execution was blocked
by that restriction; the permitted local Jupyter run completed the notebooks.
No source retrieval, external upload or publication is performed by these commands.
For just notebooks 01 and 05, append `01 05` to both notebook commands.

The builder overwrites generated notebooks. Preserve personal exploratory edits
under a new filename before rebuilding. Rerunning notebook 03 rebuilds the atlas;
set its `REBUILD_ATLAS=False` when you only want to inspect existing figures.

`outputs/provenance.json` hashes the actual inputs; notebook cell counts and
execution status are in `outputs/notebook_execution.json`. All six notebooks
were executed successfully. Numerical validation covers units/redshift factors,
grey-calibration invariance, masked gaps, common support, missing errors and
sideband error propagation against noise realizations.

Executed notebooks, figures and tables are ignored by Git because they can
contain collaboration data. The code, this guide and the science plan are kept
as small reviewable source files. No commit or publication was made.

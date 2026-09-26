# JWST Cycle 6 proposal

The current draft is [proposal.pdf](proposal.pdf), edited in `proposal.tex`.
The science case tests the dust response to luminosity changes in strongly
varying quasars and other AGN, including changing-look objects. Broad-line
appearance or disappearance is not a selection requirement. Figure 1 shows
P1823 in a high optical state, alongside its earlier infrared decline and
recovery. The text beside Figure 1 explains how these measurements and the
requested JWST spectra test the dust response.
The optical spectrum spans 0–0.6 mJy. Figure 1 has no source heading or outer box;
the caption identifies the source.

This folder contains only the active proposal, its figure, and the files
needed to rebuild them. All figure inputs are bundled under `inputs/`.

## Build

```sh
make
```

To regenerate the timing audit, selection and figure, install `requirements.txt`,
then run `make figures PYTHON=/path/to/python` followed by `make`.
Run `make clean` to remove compilation files while retaining the PDFs.
For Overleaf, upload this folder and select `proposal.tex` as the main document.

## Files

| File | Purpose |
| --- | --- |
| `proposal.tex` / `proposal.pdf` | Current text and compiled proposal |
| `apt_title_abstract.txt` | Title and abstract for APT |
| `jwstproposaltemplate_v6.sty` | Unmodified official style |
| `fig1_connected.pdf` / `.png` | P1823 light curves and spectra with emission-line labels |
| `make_fig1_connected.py` | Figure generator |
| `make_selection.py` | Reproduce the provisional sample and documented membership decisions |
| `make_timing_audit.py` | Reproduce the descriptive history audit from the bundled snapshot |
| `spherex_data.py` | SPHEREx input loading and epoch grouping |
| `inputs/` | Spectra, light curves, cutouts, selection and provenance |
| `validation.json` | Latest document and reproducibility checks |

## Working status

The compiled draft has six pages: four core pages, supplemental information
on page 5, references on pages 5–6, and the disclosure on page 6.
The working design has 11 targets and a 37–51 h planning allocation, scaled
from the previous per-target allowances rather than recomputed with ETC/APT.
ETC/APT will settle the [program category](https://jwst-docs.stsci.edu/jwst-opportunities-and-policies/jwst-call-for-proposals-for-cycle-6/jwst-proposal-types-and-categories/jwst-general-observer-go-proposals); a Medium request remains possible.
Final event-date and dust-delay audits, target-specific ETC/APT calculations, scheduling
and duplication checks remain outstanding.

Figure 1 contains only P1823, with its measured optical and infrared histories.
There is no target table in the PDF. P9694 is retained as a variable science
candidate because its optical and infrared histories show brightening; P9584
is removed from the working request. There are no dedicated controls. P1823's
18-micron feature lies beyond MRS, so its test uses the continuum and accessible
9.7-micron feature. SPHEREx availability is not asserted for every target.

Figure 1 preserves the existing P330E calibration and the previously requested
exclusion of one P1823 point; the CSV is intact. The input provenance manifest
records the extraction and calibration metadata, historical sources and current
file checksums. The plotted spectra do not establish significant SPHEREx
variability or quantitative separation of competing dust models.

SPHEREx already covers 0.75–5 microns. MIRI adds the longer-wavelength warm
continuum and accessible silicate features. NIRSpec adds spectral/spatial detail
in the overlapping near-IR range near the MIRI epoch. The test first asks whether
a fixed dust distribution under changing illumination explains the dated data,
then whether allowing the inner dust distribution to change is required.

Gold dashed lines mark expected redshifted emission-line positions, not fitted
detections. Vacuum wavelengths follow the [SDSS optical line table](https://classic.sdss.org/dr6/algorithms/linestable.php)
and [STScI infrared line lists](https://www.stsci.edu/instruments/nicmos/documents/handbooks/instrument/v5/Appendix_26.html).
The figure provenance records every plotted line position; a few labels are
offset horizontally to keep neighboring lines readable.

Date keys above the optical and SPHEREx panels use the same colors as their
spectra and the time markers on both light curves. Archival optical spectra
are grouped by year in the key; the latest spectrum is labelled 24 Sep 2026.
Each of the 20 optical epochs since 2015 has a dashed marker at its saved MJD.
The labelled 2002 spectrum is outside the displayed window. SPHEREx keys show
each visit's month range; dotted lines mark median observation dates, with
shading covering the full visit. Exact dates and colors are recorded in
`inputs/fig1_connected_provenance.json`.

Both ZTF and WISE panels start at 2015. W1 is shown as purple circles and W2 as
brown squares, on the same mJy scale without offsets. The original W1 data are
retained (filled circles: saved unWISE series; open circles: NEOWISE visits).
W2 uses the bundled target-only NEOWISE exposure snapshot in
`inputs/figure_neowise_exposures.csv`. The figure generator applies the existing
W1 frame-quality and six-month grouping rules, additionally requiring finite,
positive W2 uncertainties, and keeps visits with at least three measurements.
Visit magnitudes are medians; their approximate statistical errors are
`1.2533 * sample_std / sqrt(n)`. Conversion to mJy uses a W2 Vega zero point of
171.787 Jy from the [WISE calibration table](https://irsa.ipac.caltech.edu/data/WISE/docs/release/All-Sky/expsup/sec4_4h.html).
These plotted errors exclude absolute calibration uncertainty; no color
correction is applied. Figure provenance records all W2 visits and quality cuts.
Earlier measurements remain in the inputs, outside the displayed time window.

## Working selection provenance

The 11-object list retains the earlier seven-source core and three literature
anchors, plus P9694 as a variable candidate. It has not yet passed a selection
by event age relative to dust delay. Final membership and sample size require
event intervals, nuclear luminosities, lag distributions and model predictions
at the allowed JWST dates. The manifold records how the original pool was
assembled; it does not establish timing coverage.

The bundled pool contains the 18 NGPS-observed nuclei and 108 low-redshift
catalogue CLAGN with saved manifold projections. The reference map contains
1,960 AGN, 167 labelled as CLAGN. Selection uses the original coordinates:

1. Require a known-CLAGN fraction of at least 0.20 among 50 neighbors and at
   least 0.16 for both 25 and 75 neighbors, within 0.75 units of a reference
   point. Select the connected passing region containing the most NGPS targets.
2. Retain its seven NGPS targets: P1823, P7281, P7837, P8548, P11113, P11530
   and P10381. Require z<0.3, with P1823 retained for its optical/SPHEREx history.
3. Add all three qualifying catalogue anchors in that region: P22470, P2759
   and P16663. They have z<0.3, W3/W4 S/N>=3, AllWISE ccf=0000 and ext_flag=0.
4. The earlier selection added P9694 and P9584 outside the enriched region.
   The timing review retains P9694 as a variable science candidate and removes
   P9584, which has no demonstrated matched-control role.

`inputs/jwst_sample.csv` records the resulting list and scores;
`inputs/selection_summary.json` records the rules and the audit of all 18 NGPS
objects. `inputs/selection_pool.csv` is the compact input snapshot. These are
descriptive neighborhood fractions, not calibrated transition probabilities.
P8548's weak W4 baseline (S/N=2.5) is described in the observations section.

## Timing evidence and limits

`inputs/timing_history_snapshot.json` bundles the saved displayed optical
series, NEOWISE visits, spectral dates and prior review summaries for the
eleven candidates and the removed control. It records original source hashes.
`inputs/timing_audit.json` gives reproducible calendar-year medians, date
coverage and the selection decisions. These are descriptive summaries without
host subtraction, not fitted event times or lag measurements. Figure 1 and
its original measurements are unchanged.

The histories include P1823's 2019–2023 rise, P11530's 2018–2021 decline and
partial 2022–2023 recovery, and P9694's 2021–2022 rise followed by a brighter
state through 2025. These trends support comparing different histories, but do
not establish that the sample spans early, late and recovering dust-response
phases. That requires rest-frame event ages divided by luminosity-dependent
lag estimates, evaluated at each allowed JWST date. Two literature anchors
lack full optical display series in the bundled snapshot.

The baseline remains one visit per instrument per source. The provisional
90-day observer-frame inter-instrument limit requires source-specific checks;
shorter necessary intervals must be justified and entered in APT before
submission. No repeat visits are included without predicted detectable change.
Cycle 6 nominally runs from July 2027 to June 2028, so current optical coverage
must bridge the saved light curves to the JWST epoch. Timing sources are cited
in proposal references 3, 7, 8, 26 and 27.

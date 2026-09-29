# A data-anchored MIRI scenario figure

[Preview PDF](data_anchored_preview.pdf) · [Preview PNG](data_anchored_preview.png)

This is a scientific development preview, not a replacement for the active
proposal figure. It implements the idea of using actual sample measurements
as anchors and then showing contrasting possible MIRI outcomes.

## What was measured in this pass

- Optical g/r and W1/W2 visit summaries for all 24 targets, from the saved
  quality-selected history table: dates, number of bins, 95th/5th percentile
  flux ratio, median last-three/first-three flux ratio, and WISE last-three/
  2010 catalogue flux ratio.
- 2010 W1–W4 catalogue fluxes and statistical errors from saved photometry.
- Nine SPHEREx epoch summaries for the three locally extracted targets
  F04, F06 and R01: native-point counts and median total flux densities in
  rest 2.23–2.45 and 3.55–3.90 micron windows, grouped by gaps over 45 days.
- Rest 5050–5150 Angstrom continuum medians in 37 available optical spectra
  of five targets, including the five local NGPS reductions. These are
  total aperture continuum measurements, not host-subtracted AGN luminosities.

The inventory finds SPHEREx extractions for three targets and calibrated NGPS
products for five in the current working inputs. This is an availability
statement about those inputs, not a claim that the remaining observations do
not exist elsewhere or have not been obtained.

SPHEREx window ratios use medians in finite wavelength intervals, with point
errors recorded separately; no formal uncertainty on those medians or claim
of temperature evolution is assigned. Some windows contain only one or two
samples, as explicitly recorded. The optical medians have no new absolute
calibration or aperture correction. WISE ratios use observed total flux and
can be diluted by host light. No bolometric corrections or UV luminosities
are inferred in this pass.

The catalogue flux conversions use the existing project's WISE zero points.
They are catalogue-equivalent flux densities, without new color corrections.
The [WISE calibration documentation](https://irsa.ipac.caltech.edu/data/WISE/docs/release/All-Sky/expsup/sec4_4h.html)
explains why a final W3/W4 comparison must fit through the bandpasses and account
for spectral shape. These are historical 2010 broadband measurements, not
contemporaneous MIRI spectral points.

## Contrasting MIRI outcomes

The first candidate is R01/P9694, with real optical snapshots, WISE histories,
and 289 native SPHEREx measurements. Its latest three W1/W2 bins have median
fluxes 1.1825 and 1.0775 times the 2010 catalogue values, respectively.

The saved exploratory predictivity audit contains fits to its optical, WISE
and SPHEREx data. From that archive, this preview selects the weakest and
strongest integrated 8–13 micron outcomes among **fixed-dust** profiles with
exactly the same forecast date (2028.0) and future illumination factor (1).
It does not independently maximize every wavelength or splice spectra.
These are two intact, previously calculated spectra. Their rest-12-micron
observed flux densities are 2.10 and 23.00 mJy; their integrated 8–13 micron
fluxes differ by a factor 14.52. The spectrum includes the model's stellar/disc
terms; it is not a uniquely isolated nuclear measurement.

These selected extremes are conditional on the original finite shell-model
exploration, additional independent 5% historical scatter, energy/sublimation
constraints and an upper bound from total 2010 W3/W4 flux. The original fits
fail adequacy with the primary error model. The archive's ranges are neither
posterior credible intervals nor validated population forecasts. See the
[existing predictivity audit](../p9694_predictivity/README.md).

The 8–13 micron integration inserts exact band endpoints and integrates Fnu
over frequency. Flux densities are observed mJy at the labelled rest-frame
wavelength; the reported energy flux includes the 1/(1+z) Jacobian converting
the rest-frequency integration to the observed frequency interval.

## What this establishes and what the next physical test requires

The actual data plus conditional spectra can make a concrete figure of the
warm spectral information MIRI supplies. These two endpoints **cannot** be
labelled fixed versus evolving dust: both belong to the fixed family.
The previous audit also provides fixed/evolving spectra that are nearly
indistinguishable while meeting its conditional historical constraints.

To illustrate an evolution-specific MIRI outcome, inject a full candidate
spectrum and refit the allowed delayed fixed-distribution models, including
host and illumination alternatives. A discriminating example must require
coherent continuum/silicate differences that these models cannot reproduce.
A contrasting successful fixed-dust fit then illustrates the other scientific
outcome. First establish an adequate fit to the existing data; do not obtain
apparently separate curves by freezing unconstrained nuisance parameters or
adjusting the plotted amplitudes. No such validated mechanism classification
is claimed by this first preview.

## Reproduce

From `jwst_proposal/`:

```sh
python review/data_anchored_scenarios/measure_and_plot.py
```

Outputs: `sample_measurements.csv`, `spherex_epoch_measurements.csv`,
`optical_continuum_measurements.csv`, `r01_miri_scenarios.csv`, the preview,
and `summary.json` with exact selected curve indices and source hashes.
The proposal, Figure 1, active Figure 2, input data and existing fits are unchanged.

## SED plus physical diagnostic planes

The [new diagnostic preview](sed_and_diagnostics.pdf) retains the measured SED
and the same two MIRI spectra, replaces the redundant light-curve panel with
warm spectral ratios and a paired-silicate contrast plot, and places both
selected outcomes among 27 saved fixed/evolving conditional scenarios. All
27 use the same 2028 epoch and unit future illumination. Gray symbols are
model scenarios, not observations or probability samples; circles are fixed
dust and triangles are evolving inner dust mass. The colored outcomes are
both fixed-dust cases.

The warm ratios are exact interpolated spectral flux-density ratios at rest
4, 14 and 22 microns. They are not independently decomposed continuum
components. Their interpretation as temperature or emitting-area changes
requires the joint fit; the near-IR comparison is at the same model epoch.
The 14/4 ratios of the selected outcomes are 0.164 and 0.921, despite nearly
equal model 4-micron fluxes (10.25 and 10.75 mJy).

The feature coordinates are illustrative local contrasts
S = ln(Fnu / Fnu,local), using a power law through 6 and 14 microns for S9.7,
and 14 and 22 microns for S18. They are explicitly labelled local-continuum
indices and are not a substitute for full host/continuum/silicate decomposition.
The selected coordinates are (-0.015, -0.038) and (1.948, 0.392). These panels
show the different observable information MIRI supplies; overlap of response
families is retained and no mechanism classification is asserted.

Reproduce with `python review/data_anchored_scenarios/plot_diagnostics.py`.
`r01_diagnostic_scenarios.csv` and `diagnostics_provenance.json` record all
coordinates, exact definitions and curve indices. The active proposal is unchanged.

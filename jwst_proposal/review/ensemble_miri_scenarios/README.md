# All 24 targets: conditional MIRI spectral spans

- [Ensemble figure](ensemble_scenarios.pdf): nested spectral bands and shaded
  diagnostic planes for weaker and stronger warm endpoints.
- [All 24 targets](all_24_scenarios.pdf): each individual spectrum pair with
  its complete model-grid span in gray.
- [Scenario table](scenario_grid.csv), [selected endpoints](selected_endpoints.csv),
  [spectral quantiles](spectral_bands.csv), [provenance](provenance.json).

## What this calculation uses

The same saved AllWISE/NEOWISE W1 histories and AllWISE W1–W2 colors are available
for all 24 targets. This calculation applies the existing `forecast_warm_response.py`
physics consistently to that common dataset. It does **not** repeat the full R01
optical/WISE/SPHEREx fit for every target. The current extracted SPHEREx inputs
cover F04, F06 and R01; the provenance file states that spectra exist for all
24 but are not bundled here. The other locations were requested from the user.

The heating proxy subtracts the estimated constant W1 host contribution and
inverts the existing model's hot-response exponent. The grid varies the W1
host fraction by −0.15, 0, +0.15 (clipped to 0–0.9, duplicate values removed)
and the dust light-crossing scale by factors 0.5, 1, 2. Fixed dust with delayed
heating and an immediately adjusting sublimation boundary are both evaluated.
All forecasts use epoch 2028.0 and constant illumination outside the recorded
W1 history. The radial dust mass profile, grain mix and optical properties are
those of the existing shell model. No new dust-mass profile is fitted.

The result contains **384 spectra for 24 targets**. Each is normalized by its
own model rest-4-micron flux at the same epoch, exposing differences in spectral
shape. Normalization does not assert an observed 2028 hot-dust measurement.
The W3/W4 and SPHEREx fluxes are not newly fitted in this calculation. This is
a WISE-conditioned sensitivity grid, not a posterior or validated end-to-end
forecast of the sample's nuclear spectra.

## How the colors and shades are defined

For each object, integrate the rest 8–13 micron model spectrum over frequency
and divide by the model nu Fnu at 4 microns. Select the minimum and maximum of
this ratio over that object's host/size/response grid. This yields one weaker
blue endpoint and one stronger orange endpoint per target, 24 of each.
These labels are relative within each object's conditional grid. They do not
mean absolute weak/strong populations or fading/rising cohorts.

The top panel shows pointwise medians and nested 25–75% / 10–90% spans of those
24 endpoint spectra. Each target contributes once to each band. The shades
represent object-to-object scenario spread, not confidence intervals. A median
curve or envelope is a summary and need not be a realizable single spectrum.

The bottom panels plot all 48 endpoint coordinates. Translucent polygons are
convex hulls of each set of 24 points. They show the geometric span of the
selected scenarios, not density contours, likelihoods, or decision boundaries.
Interiors need not correspond to individually computed models. No outlier is
removed from these hulls. Circle/triangle markers identify fixed/evolving
boundary cases independently of the weak/strong color coding. Overlap is shown.

The selected endpoints include 8 fixed and 16 evolving models on the weaker
side, and 17 fixed and 7 evolving models on the stronger side. Thus color is
not an evolution classification. The median within-target strong/weak warm
ratio is about 1.135 (range 1.018–1.559) in this restricted grid.

## Spectral diagnostics

The continuum-shape plot uses Fnu(14)/Fnu(4) and Fnu(21)/Fnu(14). These are
spectral ratios, not components extracted through a full continuum/feature fit.
21 microns is used instead of 22 so both MIR anchors lie within the nominal
MRS coverage common to all 24 targets. The 4-micron normalization is a model
reference and is not within MRS coverage for every target.

The feature coordinates are illustrative contrasts S = ln(Fnu/Fnu,local).
The local power law goes through 6 and 14 microns for S9.7 and through 14 and
21 microns for S18. Full nuclear/host/silicate decomposition would replace
these simple indices in the science analysis. The shifted 21-micron anchor
means the S18 values are not exactly the same as the previous R01-only preview.

## Limits that matter for interpreting the bands

The earlier R01 preview varied fitted shell weights and consequently admitted
a much wider warm range. This uniform sample grid holds that distribution
fixed while varying host, delay scale and response prescription. Its narrower
bands must not be described as the full allowed spectral range or evidence
that adding the other targets removes R01's degeneracy. They illustrate a
consistent subset of the experiment. Equivalent data-constrained envelopes
require the remaining extracted spectra and comparable joint fits.

The inherited W1 heating proxy floors the inferred nuclear flux ratio at 0.05;
18 model curves use a history with at least one floored point. Each affected
row is flagged in `scenario_grid.csv`. This clipping is part of the existing
conditional model, not a measurement of an extremely faint nucleus.

## Numerical checks and reproduction

The temperature integral is accelerated by redistributing each cell's weight
linearly onto a fine log-temperature grid. Tests against the original direct
integrator at rest 4, 9.7 and 18 microns, for both response hypotheses and all
24 baseline objects, give maximum fractional error 1.18e−5. All 144 comparisons
are recorded in `provenance.json`. The model physics itself is unchanged.

From `jwst_proposal/`:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python review/ensemble_miri_scenarios/build_ensemble.py
python review/ensemble_miri_scenarios/plot_ensemble.py
```

The active proposal and the earlier R01 figures are unchanged. These are
ensemble development figures awaiting the equivalent full-spectral fits.

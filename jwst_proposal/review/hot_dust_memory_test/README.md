# Conditional hot-dust memory test

The requested last-year linear ZTF extrapolation has been tested. Some delayed
models improve held-out prediction, but the result is not robust enough to
claim a detected hot-dust memory signal or a measured response time.

The initial catalogue-only pass and a separate pass incorporating available
saved alert photometry are both retained. The latter adds measured points for
F02, F04, F06, F07 and R01, after matching their magnitude zero points to the
catalogue light curves. Matching offsets require at least five near-simultaneous
points. This avoids extrapolating over saved measurements; no new alert query
was made. Alert-reference calibration uncertainties are not fully modeled.

## What was tested

- The final 365.25 observed days of each optical series were fit in **flux**,
  using 14-day bins. At least three bins spanning 90 days were required. The
  slope is extended continuously from the final measured bin. Flat continuation
  is a sensitivity test. Negative predicted illumination is excluded, not clipped.
- Spectra were fit jointly per target with a constant host, variable disc, and
  one blackbody per sufficiently covered visit. Ell5 and S0 host templates were
  tested. The thermal quantity is a model-dependent bolometric blackbody proxy.
- Observed changes between first and last adequately covered SPHEREx visits
  were compared with changes in the instantaneous or delayed optical driver.
  A constant source-specific reprocessing normalization cancels in these ratios.
- The response is a uniform delay distribution from zero to twice its mean.
  Fixed mean lags of 0, 0.1, 0.25, 0.5, 1, 2 and 3 rest years were explored.
  The 0–1 year and 0–3 year ranges are reported separately because requiring
  prehistory coverage for the longer grid removes additional targets.
- A second family uses luminosity-scaled delays. The reference follows
  log10(tau/day) = -2.11 - 0.2 M_V (Koshida et al. 2014), with approximate V-band
  luminosities inferred from the median optical flux, a fixed Fnu ∝ nu^(1/3)
  continuum, and Planck18 distances. AB/Vega differences, extinction and host
  luminosity make this an approximate prior, not a measured lag. Factors
  0, 0.5, 1, 2 and 4 times that reference were tested.
- Both unit response and one globally fitted response exponent in [0,3] were
  evaluated. Each target was left out when selecting the lag and, where free,
  fitting that exponent. The current-state comparison has the same freedom.
  Targets within any one comparison are identical across all tested lags.
- Optical host sensitivity subtracts either nothing or 50% of the minimum
  observed optical flux. The latter is an assumption, not a measured host.

## Interpretation

With saved alerts included, the simple common-lag/unit-response model worsens
held-out prediction by about 12% in g and 6% in r for the baseline host template.
The luminosity-scaled/fitted-response family improves it by about 7% in g and
18% in r. These families use different common-coverage cohorts, so compare
current and delayed models **within** each row, not raw RMSE between rows.

The improvements depend on optical host subtraction, response freedom,
continuation and host template. Some gains reverse under the host-subtraction
sensitivity case. The trial families were explored in this session, including
the shorter lag range after the initial broad-grid pass; this is exploratory,
not a preregistered significance test. All comparisons are in
`all_comparisons.csv`, including less favorable results.

The main limitation is spectral model adequacy: **none of these simple joint
spectral fits passes a p > 0.01 goodness-of-fit check with the supplied independent
errors**. Reduced chi-squared ranges from about 1.6 to 62. Scaling local parameter
covariances by reduced chi-squared does not fix those physical/model residuals.
Thus even the more promising predictive improvements are conditional on
unvalidated dust decompositions. The `screened` cohort only removes boundary,
rank, extreme-uncertainty and predominantly nonphysical-driver cases; it is
not a statistically adequate-fit cohort. The `adequate` cohort is empty.

Errors propagate local spectral covariance (including shared host), an assumed
3% independent scale error per visit, and Monte Carlo optical/trend errors.
The model-selection score uses conditional mean predictions, not a fully
marginalized optical-history likelihood. It is equal-target log-change RMSE.
Bootstrap intervals resample targets only. They are not detection significances.

Saved spectra/optical points, exclusions, slope errors, extrapolation horizons,
all parameter-grid predictions and full comparison reports are retained.
SPHEREx visits are represented by median times, and unmodeled intra-visit
evolution, extinction, complex dust SEDs and extraction/calibration effects remain.
The experiment tests an ordinary echo. It does not measure a dust formation or
destruction timescale and does not rule out structural memory.

## Reproduce

From the project root:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 /opt/anaconda3/bin/python jwst_proposal/review/test_hot_dust_memory.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 /opt/anaconda3/bin/python jwst_proposal/review/test_hot_dust_memory.py --use-saved-alerts
/opt/anaconda3/bin/python jwst_proposal/review/summarize_hot_dust_memory.py
```

No proposal source or proposal PDF is modified.

Reference: https://arxiv.org/abs/1406.2078

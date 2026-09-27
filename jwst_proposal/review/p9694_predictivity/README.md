# P9694: exploratory MIRI predictivity audit

This calculation checks whether the existing P9694 data independently determine
distinct MIRI spectra for fixed and changing dust distributions. It does **not**
establish that the proposal can classify these mechanisms from one MIRI epoch.

The explored fits fail with the primary error model. With an explicitly assumed
additional 5% independent scatter, both fixed and changing distributions can fit
the historical data and predict MIRI spectra differing by less than 1% across the
sampled wavelength range. Thus the former Figure 2 (now Figure 3)'s difference between two selected models
is not a forecast of discrimination between all admissible models.

This result neither excludes AGN memory nor demonstrates that dust is fixed.
Both response families already contain memory through light-travel delays.
See [the AGN science framing](agn_memory_framing.md) for the distinction and
the implications for the proposal.

## Results and scope

| Inner mass response | Primary χ² | χ² with additional 5% scatter and physical constraints |
| --- | ---: | ---: |
| Fixed (η = 0) | 1012.93 | 132.24 |
| Mass ∝ illumination⁻⁰·⁵ | 1032.09 | 135.07 |
| Mass ∝ illumination⁻¹ | 1250.93 | 152.64 |

There are 258 measurements, 12 nonnegative shell weights and seven nonlinear
parameters. A nominal 239 degrees of freedom gives a 99% adequacy threshold of
χ² = 292.78. Nonnegative boundaries, fitted histories, approximate errors and
unmeasured covariance make this an exploratory selection criterion, not a
calibrated likelihood-ratio test. The primary results are the best solutions
found in the finite search, reevaluated with finer angular quadrature; they do
not prove that no better solution exists. One primary nonlinear optimization
hit its evaluation limit. All three final extra-scatter optimizations converged.

The [forecast audit](forecast_audit.pdf) shows primary WISE residuals and the
overlapping conditional forecast ranges. These ranges vary selected shell
weights, observing dates and future illumination; they are **not posterior
credible intervals** and do not explore every geometry or nonlinear parameter.

The [whole-spectrum counterexample](counterexample.pdf) is a stronger check
than overlapping intervals at individual wavelengths. A fixed model is refitted
to each changing-model forecast while retaining the historical-data threshold,
energy, archival W3/W4 and sublimation constraints. Both predictions use the
same date (2028.0) and future multiplier (1); no arbitrary final flux rescaling
is applied. The fixed models retain historical χ² of 206.59 and 189.50. An
independent 64-node calculation gives maximum spectral differences of 0.845%
and 0.436% for η = −0.5 and −1, respectively. These are conditional mathematical
counterexamples, not claims of subpercent physical prediction accuracy.

## Measurements and errors

- 209 line-masked native SPHEREx measurements from July 2025, December
  2025–January 2026, and June–July 2026. The same masks as the hot-dust fits
  exclude bands intersecting major lines and the 3.3 micron aromatic feature.
- 24 W1 visits: 22 NEOWISE visits plus two pre-2011.5 archival visits;
  22 W2 NEOWISE visits. Exposures require positive quality/error, clean W1/W2
  contamination flags, separation below 3 arcsec, and at least three exposures
  in each rounded 180-day group. Fluxes use the median magnitude; magnitude
  errors use 1.2533 times the sample standard deviation divided by √N.
- Three optical continuum medians over observed 0.86–0.90 micron, from 2018,
  2021 and 2026. A 15% absolute uncertainty is **assumed**, since measured
  cross-epoch calibration/slit-loss errors are unavailable.
- Archival 2010 W3/W4 total-source values are 23.52/74.81 mJy. In the constrained
  analysis the predicted nuclear-plus-model-host flux cannot exceed 1.2 times
  these values. The 20% allowance is calibration/colour-conversion slack,
  **not a measured statistical uncertainty** or a known nuclear fraction.

The primary covariance uses supplied SPHEREx errors (already including its
2% term), WISE statistical errors, and fully correlated scale terms of 2.4%
within W1 and 2.8% within W2. Cross-band covariance is not available.
The sensitivity analysis adds `(0.05 × observed flux)²` independently to every
diagonal element. This is an assumption about additional scatter, not evidence
that the physical model is adequate. Its low reduced χ² values do not validate it.

## Forward model and exploration

The model uses 0.1 micron Laor–Draine silicate and graphite grains. Six shells
per species have reference temperatures of 200/350/550/750/1000/1200 K
(silicate) and 250/450/700/1000/1300/1600 K (graphite). Shell radii follow their
radiative-equilibrium cooling integrals, with one free reference radius for
1500 K graphite. Fluxes use Qabs × Bν and nonnegative shell masses.
An Ell5 stellar template and power-law disc provide the other continuum terms.
This is an optically thin shell basis, **not** a clumpy-torus or disc-wind fit.

The illumination proxy is the host-subtracted, 60-day median ZTF g light curve.
Unknown pre-2018 illumination is parameterized with pre-2010 and 2014 knots.
The 2026 endpoint uses the NGPS/2021 optical continuum ratio, with a free relative
calibration factor. This optical proxy does not measure the full UV heating
history. Delays are integrated over 0 to a free extent times radius/c, including
redshift time dilation. Temperatures follow retarded illumination.

For the changing models, shells with reference silicate temperature ≥750 K or
graphite temperature ≥1000 K have mass multiplied by illumination^η. This
instantaneous local prescription has no independent formation time and differs
from Figure 2's moving inner boundary. η = 0 is the fixed-distribution family.

The seven parameters are optical host fraction [0,0.75], pre-2010 relative
flux [0.5,1.5], 2014 relative flux [0.5,1.3], **log10(reference radius/light-year)**
[−1.3,0.4], disc Fν slope [−1.5,1/3], NGPS relative calibration [0.8,1.2],
and delay extent [1,2]. The original archive's fourth parameter was labelled
`hot_radius_lightyears`; its stored values are logarithms. The source and archive
metadata now correct that label without changing any fitted values.

An initial 128-point Sobol exploration for each η is followed by Powell
refinement and nonnegative least squares. Final constrained fits also enforce
an upper bound on dust power from the assumed heating luminosity and omit
components exceeding illustrative sublimation limits of 1500 K for silicate
or 2000 K for graphite. These checks are necessary but not sufficient for
full radiative-transfer consistency.

SPHEREx predictions use top-hat band averages. WISE predictions integrate the
official photon response curves, with their constant-Fλ reference normalization.
Final integrations use 32 wavelength and 32 angular nodes. Doubling angular
nodes changes best-fit historical fluxes by at most 0.31%; independent
counterexample checks remain below 1%. Constant unit illumination gives identical
fixed and changing spectral bases, as required.

For each constrained optimum, shell weights are profiled at rest 6, 9.7, 12
and 18 micron within the nominal adequacy threshold. Predictions sample
2027 July 1, 2028 January 1 and 2028 June 30, with future illumination approaching
0.5, 1 or 2 times the 2026 proxy at 2028.5. Curves failing future energy or
sublimation bounds are discarded. Wavelengths are rest-frame; plotted flux
densities are observed mJy. Neither future history nor a unique warm-dust
distribution is determined by the present data.

## Files and reproduction

Run from the CLAGN root:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python jwst_proposal/test_p9694_miri_predictivity.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python jwst_proposal/review/p9694_predictivity/analyse_forecasts.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python jwst_proposal/review/p9694_predictivity/check_and_plot.py
```

`profile_fits.json` stores the initial search; `adequacy_and_forecasts.json`
stores the final conditional results and checks. `conditional_curves.npz` has
forecast arrays whose labels are in the JSON. `counterexample_spectra.csv`
contains the η = −0.5 example at 32 angular nodes. `provenance.json` records
input and source hashes. Original run-time source hashes are retained; the
later parameter-label and plot-title corrections do not alter fitted values.
No proposal spectrum or target measurement was replaced by these predictions.

Responses and calibration conventions come from the
[WISE explanatory supplement](https://irsa.ipac.caltech.edu/data/WISE/docs/release/All-Sky/expsup/sec4_4h.html).
The future-date window follows the
[Cycle 6 call](https://jwst-docs.stsci.edu/jwst-opportunities-and-policies/jwst-call-for-proposals-for-cycle-6).
Grain inputs are the [Draine optical-property tables](https://www.astro.princeton.edu/~draine/dust/dust.diel.html);
the host is from the [SWIRE templates](https://www.iasf-milano.inaf.it/~polletta/templates/swire_templates.html).

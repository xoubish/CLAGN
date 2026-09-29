# Figure 2: a conditional measurement experiment

2026-09-29. The active figure now connects R01's measured light curves to a mock
MIRI spectrum, refits competing shell responses, and maps three complementary
constraints to the scientific questions. Continuum and silicate precision
comes from the saved simulation; PAH and spatial constraints are not simulated.
It replaces the two extreme weak/strong curves. The
proposal's old 4.5σ/3.7σ ensemble claims were removed: they were not calibrated
model-recovery results. No proposal PDF was compiled.

## Result

An evolving inner-mass model (η = −0.5) supplies the mock truth. Both it and a
fixed-mass model can fit the same realization after refitting nonnegative shell
weights and accounting for correlated nuisance terms. Their MIRI quadratic
statistics are 152.96 and 152.70 for 174 bins, respectively. These are **not**
Bayesian evidences or calibrated model-selection statistics. Historical
constraints remain satisfied under the adopted additional 5% scatter; the fixed
model reaches the historical adequacy boundary. The third η = −1 response was
also fitted and is retained in the output, not hidden from the audit.

For 4000 independent noise realizations of this spectrum:

| Observable | Input value | Statistical standard deviation | With assumed correlated errors |
|---|---:|---:|---:|
| Observed 8–13 µm band flux | 2.318 × 10⁻¹² erg s⁻¹ cm⁻² | 0.31% | 5.38% |
| Local S9.7 | 0.945 | 0.023 | 0.082 |
| Local S18 | 0.094 | 0.037 | 0.041 |

This precision concerns the extracted spectrum and the stated error model. It
does not imply 5% knowledge of nuclear luminosity, dust mass, covering factor,
geometry, or historical illumination. A band flux includes broad dust features;
it is not a separately decomposed feature-free continuum integral.

## Assumptions and display

- Rest-frame wavelength, observed mJy. Band integration includes exact rest
  8 and 13 µm boundaries and the redshift conversion of the observed frequency
  interval. Nominal observed MRS coverage is 4.9–27.9 µm, with 174 contiguous
  bins and eight quadrature nodes per bin.
- Statistical S/N is **assumed at the requested goal**, 50 through rest 12 µm
  and 30 at/above 18 µm, interpolated between. This is not a new full-band ETC
  run for R01 and does not establish sensitivity for all 24 targets.
- Correlated errors: a 3% normalization coefficient, a 3% coefficient on
  ln(λ/12)/ln(22/6), and one smooth residual host mode normalized to 10% of
  the flux at 12 µm. The host shape is a 200 K modified blackbody with β=1.5.
  These are sensitivity assumptions, not measured calibration/host priors.
- Covariance is diagonal statistical noise plus outer products of those three
  modes. Refits use the full covariance. Displayed fits include the conditional
  best nuisance offsets; the intrinsic curves are saved separately. The
  residual shading shows marginal pointwise errors, which are correlated.
- The local indices use log-linear continua through rest 6/14 and 14/21 µm.
  They are not the published spline-based silicate strengths or a full physical
  decomposition. Narrow lines and PAHs are absent from the mock; a real analysis
  must fit or mask them.
- Grey ranges come from 27 archived conditional curves at a common 2028 epoch
  and future multiplier of one, with fixed/evolving shell weights explored.
  The envelope is not a posterior or complete model family. The former panel C
  grey-to-orange ranges have been removed to avoid suggesting a measured
  before/after posterior reduction. Panel C now explains warm emission,
  silicate profiles and host contribution, retaining the conditional precision.
- Panel A normalizes each measured band to its own 2018–2020 median and leaves
  the unobserved future blank. The archival W2 point is placed at the W1 mean
  AllWISE epoch as a plotting approximation; the input export lacks W2's mean
  epoch. The forecast assumes constant illumination after the 2026 optical
  anchor; this is not inferred from the plotted optical light curve. Open g
  points use Figure 1's overlap-corrected R01 alert photometry. No NGPS marker
  has been restored, and no extrapolated optical observations are plotted.

### Added measured epochs

The former SPHEREx date shading has been replaced by three synthetic W1 points
computed from the measured visit spectra, sharing the WISE W1 normalization.
They are 7.595 ± 0.154, 7.889 ± 0.160 and 8.029 ± 0.163 mJy. The photon-response
integration retains more than 99.99999% of the W1 response; the supplied
calibration term is propagated coherently within each visit. Piecewise-linear
interpolation is used between native spectral measurements; interpolation
uncertainty and additional inter-instrument aperture differences are not included.

The NGPS point and its 2018 reference marker were removed at the user's request.
The extraction recipe and numerical values are saved in
`../../inputs/fig2_epoch_points.json`; regenerate with the figure script.

## What was checked

All mock fits satisfy the adopted historical, energy and sublimation constraints.
The marginalized-covariance statistic equals the explicit Gaussian nuisance-fit
statistic. Doubling angular quadrature from 64 to 128 nodes changes the injected
spectrum by at most 0.115%. The mock was not chosen by searching random seeds for
large separation: the seed is fixed at 96942026. Original figure/source backups
are in `older/before_recovery_figure`.

An additional three-start instantaneous-response fit is saved in
`../instantaneous_check.json`. Its best historical statistic is approximately
377.8, above the previous 292.8 screening threshold even with 5% extra scatter.
Thus that particular instantaneous shell prescription already struggles with
the existing record; rejecting it in a joint fit would not demonstrate new
discrimination supplied by MIRI. This finite search does not rule out flexible
equilibrium torus models or establish a physical memory detection.

The nonlinear geometry/illumination parameters are held at the separate
historical optima during these mock refits. Full torus-library fitting and
sample-wide null/recovery calibration remain outstanding. This revision is a
defensible measurement illustration, not the requested ultimate 24-object
mechanism-discrimination forecast.

## Reproduction

From the project root:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 /opt/anaconda3/bin/python jwst_proposal/review/validate_fig2_recovery.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 /opt/anaconda3/bin/python jwst_proposal/make_fig2_miri.py
```

`summary.json` records assumptions and refits; `mock_spectrum.csv` separates
injected truth, noisy data, intrinsic models and fitted observed models;
`measurement_draws.npz` saves covariance and diagnostic realizations.

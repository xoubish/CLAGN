# A testable hot-dust history model

This is a model specification, not a fit or a detected correlation. The previous
raw continuum-color plots do not test the models below. No proposal text or
figures are changed by this note.

## 1. Start with a delayed illumination model

Let L_i(t) denote the intrinsic illuminating luminosity, with host-subtracted
optical continuum as an imperfect proxy. Use rest-frame times. Define

    L_eff,i(t) = integral Psi_i(u) L_i(t-u) du
    integral Psi_i(u) du = 1, u >= 0
    M_i(t) = ln[L_eff,i(t) / L_i(t)]

Start with one fixed kernel family, e.g. a uniform distribution of delays from
0 to 2 tau_i (a thin spherical shell's geometrical response), as an economical
approximation rather than an assertion of spherical torus geometry. Compare
with a zero-delay model and one other prespecified kernel as a sensitivity test.
Use luminosity-informed lag priors with scatter, not a separate lag selected
to maximize correlation in each object. Reverberation measurements support an
approximate luminosity^0.5 scaling; W1 is itself reprocessed emission and must
not be substituted for the optical/UV driver without an explicit transfer model.

For a *bolometric* reprocessed component with fixed interception fraction A_i,
the simplest approximation is D_i(t) = A_i L_eff,i(t). Then

    ln[D_i(t)/(A_i L_i(t))] = M_i(t).

A finite near-IR band is not bolometric: its fraction changes with temperature.
For that band, synthesize the observed flux from the temperature-dependent SED,
or use a local responsivity model for small perturbations,

    delta ln D_lambda(t) = eta_lambda integral Psi_lambda(u) delta ln L(t-u) du.

The general spectral calculation should integrate contributions from different
radii/delays before measuring a fitted color; a single effective luminosity and
temperature are only an initial low-dimensional approximation.

For slowly varying luminosity, Taylor expansion gives

    M(t) approximately -tau_mean d ln L(t)/dt
         = -tau_mean [dL(t)/dt] / L(t).

This is the physically scaled version of "slope / today". A fading driver gives
positive excess relative to instantaneous equilibrium; a rising driver gives
a deficit, conditional on the fixed-dust approximation. Flares/recoveries need
the full convolution: a zero present slope can coexist with delayed emission.
The approximation fails for large curvature/rapid changes and is not a new
empirical estimator of lag.

## 2. A separate hypothesis: slow evolution of the emitting radius

After fitting the echo-only baseline, ask whether a slowly adjusting inner
dust distribution is needed. An illustrative one-radius toy model is

    R_inner(t)^2 = K L_structure(t)
    dL_structure/dt = [L(t) - L_structure(t)] / t_adjust
    T(t)^p proportional to L_heat(t) / R_inner(t)^2,

where L_heat is the appropriately delayed illuminating state and p = 4 for
grey grains; an emissivity-dependent exponent is a model assumption, not a
measured universal constant. Relative to an equilibrium sublimation radius,

    (T/T_sub)^p approximately L_heat / L_structure.

A nucleus that faded can retain a larger inner radius and cooler dust until
dust repopulates smaller radii. A brightening nucleus can heat surviving dust
toward sublimation, after which destruction changes the emitting area. Enforce
sublimation limits; the toy relation cannot predict unlimited temperature rises.
Destruction and formation may need different timescales, but do not add both
freely before the data support even a common adjustment time.

This is structural memory, distinct from a light echo and from thermal energy
being stored for years by individual grains. There is no universal predicted
sign for raw 3.7/2.2 micron color: wavelength-dependent lags, responsivities,
host light, disc light and multiple dust temperatures all contribute.

## 3. Fit the measurements, then display the history diagnostic

1. Fit each source's SPHEREx visits jointly with a constant host component,
   variable disc and dust SEDs, source extinction/nuisance calibration, actual
   passbands and observation dates. Visits span days to weeks and are not
   simultaneous full spectra. Check detector-boundary steps and extraction
   effects before reporting temperature changes. Use more than one plausible
   host template to test identifiability.
2. Reconstruct the optical driver from original photometry, maintaining host,
   extinction, calibration and interpolation uncertainties. Fit WISE as another
   delayed response. Optical observations from the same SPHEREx visits can
   help anchor current state but share calibration/decomposition uncertainty.
3. Compare an instantaneous response with the prespecified echo model on the
   same data and nuisance terms. Use source-specific, time-independent dust
   normalizations and partial pooling of lag/response parameters. With only
   two or three visits per source, a freely adjustable lag, covering factor,
   temperature and host for every visit would not be identifiable.
4. Only then test whether an adjustment-time parameter improves prediction
   beyond the echo-only model. Assess held-out visits where possible, simulated
   null data using the actual cadence, and sensitivity to individual targets.
   Do not choose timescales or wavelength windows to maximize the final plot.
5. Display modelled history contrast M versus the dust excess relative to
   current illumination, with multiple visits connected per target. Source
   offsets absorb persistent covering-factor differences; changes within a
   target are especially useful. Shared L(t) on both axes creates correlated
   errors, so fit the original fluxes jointly, not a naive ratio regression.

For structural memory, an additional diagnostic is fitted temperature against
L_heat/L_structure, after testing host/temperature degeneracies. It is secondary
to the echo baseline, not an interchangeable way to search for any correlation.

## 4. What the existing local data can presently support

The first atlas has finite data for 23 targets, usable paired continuum colors
for 16, and only two or three observation groups per target. These support an
initial low-dimensional comparison, not guaranteed individual lag recovery.

`driver_coverage.csv` audits the existing 90-day optical bins at 60 SPHEREx visit
midpoints. For trial mean rest-frame delays 0.25, 0.5 and 1 year, only 10 visits
have the g-band record bracketing the full [t-2tau, t] interval, and 13 have
r-band bracketing. All 23 targets' latest SPHEREx visits are later than their
last retained g-band bin. This is **not** a claim that newer raw observations
do not exist: existing alert photometry for R01, for example, extends beyond
the binned audit. The next data step is therefore to assemble recent optical
photometry, rather than extrapolate the old endpoints as measured current states.
Bracketing alone also does not ensure sufficiently dense or accurate sampling.

The existing long-term WISE slopes describe 2014–2024 behavior. They should not
be called the contemporaneous driver slope at a 2025–2026 SPHEREx visit.

## Primary literature

- [Koshida et al. 2014](https://arxiv.org/abs/1406.2078): measured dust delays
  in 17 Seyferts and the lag–luminosity relation.
- [Hönig & Kishimoto 2011](https://arxiv.org/abs/1109.3465): wavelength-dependent
  smoothing, dust temperature response, and transfer-model predictions.
- [Kishimoto et al. 2013](https://arxiv.org/abs/1308.6517): interpretation of
  NGC 4151 radius evolution in terms of earlier luminosity and multi-year
  structural response. This is a hypothesis to test, not a universal prior.
- [Schnülle et al. 2015](https://www.aanda.org/articles/aa/abs/2015/06/aa25733-15/aa25733-15.html):
  temperature changes without inferred sublimation-driven radius growth over
  their monitoring interval.
- [Lyu & Rieke 2021](https://arxiv.org/abs/2011.07638): multiple dust components,
  distinct hot and warm delays, and a roughly stable hottest inner edge over
  a longer NGC 4151 record. It cautions against treating slow radius adjustment
  as already established by any temperature or flux variation.

# P9694 exploratory hot-dust continuum fits

Update: the subsequent [pooled fit and figure revision](mean_fit.md) adds an
illustrative mean curve and rounded ~1200 K colour-temperature label. The report
below records the earlier separate-visit analysis and its original figure decision.

The fits give an approximate hot-dust colour temperature, but **do not support
precise temperature labels in the proposal or a robust temperature-variability
claim**. The fixed-disc baseline gives about 1160–1180 K; a freely fitted red
power law gives about 1330–1390 K. Both leave structured residuals larger than
the supplied errors. Figure 1 and the proposal PDF were therefore not changed.

| Model | Jul 2025 T (K) | Dec 2025–Jan 2026 T (K) | Jun–Jul 2026 T (K) | χ² / nominal dof |
| --- | ---: | ---: | ---: | ---: |
| Ell5 host, fixed disc slope 1/3 | 1162 | 1157 | 1182 | 3.17 |
| Ell5 host, free common power-law slope | 1383 | 1333 | 1393 | 2.21 |
| Ell2 host, fixed disc slope 1/3 | 1228 | 1212 | 1236 | 3.85 |
| S0 host, fixed disc slope 1/3 | 1151 | 1147 | 1172 | 3.06 |
| Ell5 host, modified blackbody β=1 | 947 | 941 | 959 | 4.99 |
| Ell5 host, wider line masks | 1157 | 1148 | 1185 | 2.93 |

These alternatives are sensitivity checks, **not a confidence interval**. In
particular β=1 changes the assumed emissivity as well as the inferred temperature
and provides a worse fit. The free power-law index is α =
-1.136, where Fν ∝ ν^α. This red slope is a
phenomenological continuum; it is not evidence for a physical thin-disc spectrum.
The shared host amplitude falls from 2.60 mJy in the
baseline to 0.40 mJy with the free
slope (both at rest 1.6 µm), demonstrating a substantial decomposition degeneracy.
Dust amplitudes and any extrapolated dust luminosities are correspondingly
model dependent; no physical bolometric luminosity is reported.

## Inputs and method

- Source P9694, J004607.98+090720.9, z = 0.2377224; 289 native, unbinned measurements.
  Rest wavelength coverage is approximately 0.60–4.03 µm.
- Simultaneous fit to three visits: constant host normalization and common disc
  slope; independent disc amplitude, dust amplitude and dust temperature per visit.
  All component amplitudes are nonnegative. Temperatures are bounded at 600–2400 K.
- Baseline host: SWIRE Ell5; Ell2 and S0 are sensitivity alternatives. The library
  Fλ spectra are converted to Fν before fitting. Ell13 was downloaded but not used.
- Each predicted Fν is averaged uniformly in wavelength over the exported centre
  ± half-width, using 32-point Gauss–Legendre quadrature. This is a top-hat
  approximation: the full instrument response curves were not supplied.
- Fixed masks exclude any band touching ±2% in wavelength around Hα, Paδ,
  He I/Paγ, Paβ, Paα, Brγ, Brβ or Brα, and rest 3.20–3.35 µm around PAH 3.3.
  The wider-mask check uses ±3%. No residual-based clipping is applied.
  The baseline uses 209 points (67, 71, 71
  per visit); all supplied points remain visible in the diagnostic plots.
- The supplied total flux errors are used without rescaling. Their 2% calibration
  term is retained, but covariance is unavailable and the fit treats errors as
  independent. Nominal χ² probabilities and profile intervals are conditional on
  this approximation; active zero-amplitude boundaries further limit their interpretation.
- No foreground or intrinsic extinction correction is applied. Each visit spans
  roughly two to three weeks; within-visit source variability is not modelled.
  Noncontemporaneous optical/NGPS spectra are not used as rigid flux constraints.
- Model components are evaluated in the rest frame, with amplitudes in observer
  mJy; the common redshift/distance normalization is absorbed in the amplitudes.

## Formal errors and visit comparison

For reproducibility, the baseline profile-likelihood intervals (Δχ²=1, all other
parameters refitted) are recorded below. **These are not credible total errors**:
model dependence and structured residuals dominate this apparent precision.

| Visit | Baseline T (K) | Conditional 68% interval (K) |
| --- | ---: | ---: |
| Jul 2025 | 1162 | 1155–1169 |
| Dec 2025–Jan 2026 | 1157 | 1150–1165 |
| Jun–Jul 2026 | 1182 | 1176–1189 |

- baseline: shared T = 1169 K; allowing separate visit temperatures improves χ² by 8.64 for two extra parameters (nominal p = 0.0133).
- free_disc_slope: shared T = 1364 K; allowing separate visit temperatures improves χ² by 7.79 for two extra parameters (nominal p = 0.0204).

These nominal comparisons assume an adequate spectral model and independent
errors, which the residuals do not establish. They must not be presented as a
robust detection of heating or cooling. An extinction-aware continuum model,
better host constraints, instrumental response/covariance and a check of
within-visit variability are needed before making that claim.

## Files and verification

- [Baseline components and residuals](continuum_fits.pdf)
- [Free-slope components and residuals](free_slope_fits.pdf)
- [Results, profiles, assumptions and input/script hashes](results.json)
- [Every measurement, mask decision and baseline residual](measurements.csv)

Reproduce from the project root with
`python jwst_proposal/fit_spherex_hot_dust.py`.
Noiseless synthetic data recover their input temperatures with maximum error
0 K.
Increasing the quadrature to 64 nodes changes predicted band fluxes by at most
0.034%.
Multistart optimization was used for all six model variants. The diagnostics
show that the main limitation is model adequacy, not numerical convergence.

Host templates and their units: [SWIRE library, Polletta et al. (2007)](https://www.iasf-milano.inaf.it/~polletta/templates/swire_templates.html).
The use of near-IR spectral decomposition to infer dust colour temperature, and
its dependence on emissivity assumptions, has precedent in
[Landt et al. (2019)](https://arxiv.org/abs/1908.01627).
Template download provenance and checksums are in `../../inputs/host_templates/provenance.json`.

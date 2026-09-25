# P1823 SPHEREx continuum variability check

Regenerate with `python check_a3_spherex_variability.py --target P1823`.

The three epochs are Jul–Aug 2025, Jan–Feb 2026, Jul–Aug 2026. Fit an independent quadratic continuum to each epoch at its actual wavelengths, then average the fitted F_nu over identical wavelength intervals. This accounts for unequal wavelength sampling without imposing the same continuum shape.

| Observed interval (µm) | Mean fluxes in chronological order (mJy) | First-to-last change | Nominal significance | Shared 2% scenario | Shared 2% + excess scatter |
|---|---|---|---|---|---|
| 0.85–1 | 0.182, 0.193, 0.213 | +17.3% | 1.1σ | 1.1σ | 1.1σ |
| 1.15–1.5 | 0.236, 0.235, 0.279 | +18.2% | 2.3σ | 2.2σ | 1.6σ |
| 2.2–2.6 | 0.438, 0.416, 0.425 | -2.9% | -0.7σ | -0.6σ | -0.5σ |
| 2.6–3.5 | 0.620, 0.645, 0.620 | -0.0% | -0.0σ | 0.1σ | 0.1σ |
| 3.5–4 | 0.859, 0.855, 0.920 | +7.1% | 2.2σ | 1.7σ | 1.5σ |
| 4–4.5 | 1.128, 1.106, 1.178 | +4.4% | 1.2σ | 1.0σ | 0.9σ |
| 4.5–5 | 1.218, 1.256, 1.282 | +5.3% | 1.4σ | 1.1σ | 1.0σ |
| 4–5 | 1.179, 1.175, 1.227 | +4.0% | 1.6σ | 1.1σ | 1.0σ |

The nominal column treats the supplied total errors as independent. The CSV includes a 2% calibration term. In the sensitivity scenario this same term is shared within each epoch, independently between epochs, while retaining the supplied measurement variance; it is not added twice. A scale error common to all epochs would instead largely cancel. The CSV alone cannot determine the correct covariance.

The final column additionally inflates measurement variances by max(1, reduced chi-square of the epoch-specific continuum fits), keeping the shared calibration term fixed. This is a residual-scatter sensitivity check, not a measured noise model.

Linear, quadratic and cubic fits give first-to-last changes of 4.0–4.2% for 4–5 µm. The P1823 data do not provide a robust detection of continuum variability. The shorter 1.15–1.5 micron interval has a possible increase, but its significance also falls when excess residual scatter is allowed. The isolated 4.06996 micron point is excluded in the primary calculation and restored with --include-excluded for sensitivity analysis; the first-to-last comparison is unaffected because the point belongs to the middle epoch.

## Limits

- Exploratory, uncorrected local significances; aggregate 4-5 um overlaps two narrower windows and is not independent of them.
- Smooth fits use wavelength centres, not full spectral response curves. Band means are model estimates, not direct synthetic photometry.
- CSV does not specify cross-measurement or cross-epoch calibration covariance; the correlated scenario is illustrative, not established.
- Background, aperture, processing-version and other systematics have not been independently validated with control sources.
- Nominal significances assume the fitted smooth continuum is adequate; inspect fit chi-square and degree sensitivity.

## Joint test of a constant continuum spectrum

Use the seven non-overlapping intervals and all three epochs, with one constant mean per interval (14 degrees of freedom). The overlapping 4–5 micron aggregate is omitted to avoid double counting. The shared 2% scenario correlates the epoch scale across wavelengths as well as within each interval.

- independent_errors: chi-square = 24.56, dof = 14, p = 0.0392.
- epoch_correlated_calibration: chi-square = 19.06, dof = 14, p = 0.1625.
- excess_scatter_and_correlated_calibration: chi-square = 11.97, dof = 14, p = 0.6085.

The nominal result is marginal and does not survive the shared-calibration or excess-scatter checks. Individual-window results are exploratory and are not corrected for the multiple windows/epoch pairs inspected. The nominal smallest window p-value (about 0.03) also does not survive an eight-window Bonferroni correction.


Instrument context: [IRSA spectral calibration products](https://irsa.ipac.caltech.edu/data/SPHEREx/docs/spherex_spectral_calibrations.html) provide the spectral response curves needed for a full response-based comparison. The numerical findings above come from the supplied CSV, not from that documentation.

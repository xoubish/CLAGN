# A3 SPHEREx continuum variability check

Regenerate with `python check_a3_spherex_variability.py --target A3`.

The three epochs are Jul 2025, Dec 2025–Jan 2026, Jun–Jul 2026. Fit an independent quadratic continuum to each epoch at its actual wavelengths, then average the fitted F_nu over identical wavelength intervals. This accounts for unequal wavelength sampling without imposing the same continuum shape.

| Observed interval (µm) | Mean fluxes in chronological order (mJy) | First-to-last change | Nominal significance | Shared 2% scenario |
|---|---|---|---|---|
| 2.6–3.5 | 6.822, 7.123, 7.246 | +6.2% | 5.4σ | 2.1σ |
| 3.5–4 | 8.395, 8.917, 8.986 | +7.0% | 5.1σ | 2.3σ |
| 4–4.5 | 9.672, 9.860, 10.243 | +5.9% | 4.7σ | 1.9σ |
| 4.5–5 | 10.536, 10.877, 10.862 | +3.1% | 2.4σ | 1.0σ |
| 4–5 | 10.090, 10.367, 10.550 | +4.6% | 5.2σ | 1.5σ |

The nominal column treats the supplied total errors as independent. The CSV includes a 2% calibration term. In the sensitivity scenario this same term is shared within each epoch, independently between epochs, while retaining the supplied measurement variance; it is not added twice. A scale error common to all epochs would instead largely cancel. The CSV alone cannot determine the correct covariance.

Linear, quadratic and cubic fits give first-to-last changes of 4.4–4.6% for 4–5 µm. The 4.5–5 µm epoch ordering is not monotonic. These results support a candidate few-percent continuum change, not a secure claim of monotonic brightening. Proposed wording: “Three SPHEREx epochs reveal a rising near-infrared continuum with indications of continued evolution, providing the short-wavelength context for the proposed MIRI spectroscopy.”

## Limits

- Exploratory, uncorrected local significances; aggregate 4-5 um overlaps two narrower windows and is not independent of them.
- Smooth fits use wavelength centres, not full spectral response curves. Band means are model estimates, not direct synthetic photometry.
- CSV does not specify cross-measurement or cross-epoch calibration covariance; the correlated scenario is illustrative, not established.
- Background, aperture, processing-version and other systematics have not been independently validated with control sources.
- Nominal significances assume the fitted smooth continuum is adequate; inspect fit chi-square and degree sensitivity.

Instrument context: [IRSA spectral calibration products](https://irsa.ipac.caltech.edu/data/SPHEREx/docs/spherex_spectral_calibrations.html) provide the spectral response curves needed for a full response-based comparison. The numerical findings above come from the supplied CSV, not from that documentation.

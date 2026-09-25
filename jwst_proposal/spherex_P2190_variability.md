# P2190 SPHEREx continuum variability check

Regenerate with `python check_a3_spherex_variability.py --target P2190`.

The three epochs are Jun–Aug 2025, Jan–Feb 2026, Jul–Aug 2026. Fit an independent quadratic continuum to each epoch at its actual wavelengths, then average the fitted F_nu over identical wavelength intervals. This accounts for unequal wavelength sampling without imposing the same continuum shape.

| Observed interval (µm) | Mean fluxes in chronological order (mJy) | First-to-last change | Nominal significance | Shared 2% scenario |
|---|---|---|---|---|
| 2.6–3.5 | 0.558, 0.563, 0.581 | +4.1% | 1.7σ | 1.1σ |
| 3.5–4 | 0.753, 0.733, 0.746 | -0.9% | -0.3σ | -0.2σ |
| 4–4.5 | 0.856, 0.925, 0.934 | +9.1% | 2.2σ | 1.9σ |
| 4.5–5 | 0.916, 1.014, 1.031 | +12.6% | 2.9σ | 2.4σ |
| 4–5 | 0.890, 0.972, 0.982 | +10.3% | 3.5σ | 2.5σ |

The nominal column treats the supplied total errors as independent. The CSV includes a 2% calibration term. In the sensitivity scenario this same term is shared within each epoch, independently between epochs, while retaining the supplied measurement variance; it is not added twice. A scale error common to all epochs would instead largely cancel. The CSV alone cannot determine the correct covariance.

Linear, quadratic and cubic fits give first-to-last changes of 9.9–10.6% for 4–5 µm. P2190 shows a candidate long-wavelength brightening concentrated between summer 2025 and early 2026; the last two epochs are consistent with each other. There is no comparable significant change at 3.5–4 µm. Proposed wording: “Three SPHEREx epochs reveal a rising near-infrared continuum with tentative long-wavelength brightening, providing the short-wavelength context for MIRI spectroscopy.”

## Limits

- Exploratory, uncorrected local significances; aggregate 4-5 um overlaps two narrower windows and is not independent of them.
- Smooth fits use wavelength centres, not full spectral response curves. Band means are model estimates, not direct synthetic photometry.
- CSV does not specify cross-measurement or cross-epoch calibration covariance; the correlated scenario is illustrative, not established.
- Background, aperture, processing-version and other systematics have not been independently validated with control sources.
- Nominal significances assume the fitted smooth continuum is adequate; inspect fit chi-square and degree sensitivity.

Instrument context: [IRSA spectral calibration products](https://irsa.ipac.caltech.edu/data/SPHEREx/docs/spherex_spectral_calibrations.html) provide the spectral response curves needed for a full response-based comparison. The numerical findings above come from the supplied CSV, not from that documentation.

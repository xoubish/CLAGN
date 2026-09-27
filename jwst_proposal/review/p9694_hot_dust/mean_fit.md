# P9694 pooled SPHEREx continuum

The figure now shows one dashed **total host + disc + blackbody** curve fitted
to the continuum measurements from all three visits, with a rounded,
model-dependent hot-dust colour temperature of **1200 K**. The exact baseline
optimum is 1167.0 K. All three component amplitudes and the
temperature are shared; this is not an average of the earlier fitted temperatures.
The fixed-disc baseline uses Fν ∝ ν^(1/3) and the SWIRE Ell5 host template.

The visits span **July 2025–July 2026**, approximately one year. Pooling gives a
useful descriptive average; it does not establish simultaneity or remove
variability, host/disc degeneracy, extinction or calibration systematics.
All 289 measurements are displayed unchanged; 209 continuum points enter the
fit with their supplied errors and band widths. The earlier flux-independent
emission-line/PAH masks are retained. No binning or inter-visit rescaling is used.

| Pooled model | Colour T (K) | χ² / nominal dof |
| --- | ---: | ---: |
| baseline | 1167 | 3.44 |
| free_disc_slope | 1379 | 2.50 |
| young_host | 1225 | 4.02 |
| S0_host | 1159 | 3.33 |
| modified_blackbody_beta1 | 949 | 5.14 |
| wider_line_masks | 1160 | 3.28 |

The baseline has χ² = 704.2 for 205 nominal degrees of freedom,
compared with 631.2 for the earlier model allowing
separate visit disc amplitudes, dust amplitudes and temperatures (six additional
parameters). These residuals exceed the supplied errors. The figure's 1200 K
label is consequently an illustrative, conditional colour temperature without
a formal precision claim. The alternative models in this table are sensitivity
checks, not a confidence interval. Pooling does not resolve the earlier model
dependence, so no temperature-change or bolometric luminosity claim is made.

Reproduce with `python jwst_proposal/fit_spherex_mean.py`, then regenerate Figure 1.
[Machine-readable fit and curve provenance](mean_fit.json).
The [earlier fit report](README.md) documents shared assumptions, line masks,
template sources, uncertainty limitations and the original separate-visit fits.
Numerical checks recover a known synthetic temperature to better than 1 K;
64-node quadrature changes predicted band fluxes by at most 0.035%.

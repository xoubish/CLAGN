# R06 first ETC sensitivity trial

Input and output are from the user-supplied archive
`../downloads/web_exports/wb325241_c7_2026-09-26_18.58.14.tar`, not a local Pandeia run.
The archive checksum is in `analysis.json`. The original archive is unchanged.

## Verified in the archive

- Point source at scene center, no emission lines or extinction.
- Power law in fnu, wavelength exponent +2; WISE W3 normalization 0.95 mJy.
  This is an assumed 50%-of-archival-total nuclear-flux scenario, not a fit.
- Saved redshift is 0.211, rather than the sample value 0.2111. Correct this
  for later template work. For this featureless power law normalized after
  redshifting, the difference does not change the normalized continuum shape.
- MIRI/MRS channel 3, short; FULL, FASTR1, 27 groups, 1 integration, nexp=4.
- Source-centered circular extraction radius 1 arcsec; off-scene sky strategy.
- Background at 12.81 microns is 37.45 MJy/sr, changed from the original
  example's 31.18 MJy/sr. Input JSON records the background spectrum but not
  coordinates or percentile; the requested position was
  08:38:26.50 +37:19:06.7 with Medium background.
- The saved spectrum reproduces the reported native S/N: 5.3725 at 12.81 um.
- The output grid covers 11.5316–13.4781 um; this does not cover the target's
  entire observed warm band (9.69–15.74 um), or rest 12 um at observed 14.53 um.

## Timing from the user's copied ETC scalar report

On-source exposure: 299.70 s (four exposures of 74.93 s, rounded).
Total required for the on-source plus off-source strategy: 599.41 s.
The extra sky exposure time is not included in the displayed 299.70 s.
These are not APT charged times and exclude the full overhead accounting.
Reported maximum saturation fraction: 0.00104.
The archive does not include the complete scalar report or warnings.

## Approximate R=100 calculation

`../bin_archive.py` reads the exported extracted flux and noise on the
907-point output wavelength grid. It sums whole spectral pixels with centers
inside lambda +/- lambda/200 and divides by the quadrature sum of their noise.
The exact resolving power changes slightly with rounding to whole pixels;
`r100_estimates.csv` records it. The coarser 42-point input/model wavelength
grid is not used as the output pixel grid. Flux/noise was checked against the
exported S/N at every output pixel.

Approximate S/N at bin centers 11.70, 12.00, 12.50, 12.81, 13.00, 13.30 um is
37.04, 32.81, 37.79, 41.10, 39.13, 35.05, respectively. These estimates assume
independent spectral-pixel errors. The exported 1-D arrays supply no spectral
covariance; this calculation cannot establish the uncertainty after real cube
reconstruction, calibration or host subtraction. Host photon noise is also
absent from this point-source-only scene. These are planning estimates, not
achieved or guaranteed precision. No proposal sensitivity claim was updated.

Next comparison: retain the exposure setup and change only the extraction
radius from 1.0 to 0.6 arcsec. Compare the signal loss and noise reduction
before increasing exposure time; this is a trial radius, not an optimum.

References checked:
- https://jwst-docs.stsci.edu/jwst-exposure-time-calculator-overview/jwst-etc-outputs-overview/jwst-etc-downloads
- https://jwst-docs.stsci.edu/jwst-exposure-time-calculator-overview/jwst-etc-calculations-page-overview/jwst-etc-strategies/jwst-etc-ifu-strategies
- https://jwst-docs.stsci.edu/known-issues/miri-known-issues/miri-mrs-known-issues

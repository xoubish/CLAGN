# Single-galaxy Figure 1 preview

Example: **A3 / P9694 / J004607.98+090720.9**, z = 0.2377224.

Open `fig1_single_example_preview.pdf` or `.png`. Regenerate with:

```sh
python make_fig1_single_example.py
```

The active `proposal.tex` and its existing Figure 1 are unchanged. This alternative gives the
example more space without removing the six-target table. Figure 2 still models
the faded anchor A1, not this brightening example; its numerical predictions
must not be attributed to A3.

## Suggested caption

**One nucleus connects the optical history to the infrared continuum.**
A3 (J004607.98+090720.9, z = 0.238): (a) SDSS gri image, 15 arcsec across,
with a 5-arcsec scale bar.
(b) ZTF g and r light curves and WISE W1 flux history; filled symbols are the
bundled unWISE measurements and open symbols the NEOWISE measurements. Dashed
lines mark the three optical spectral epochs, with matching colours in panel
(c). Coloured shaded intervals mark the SPHEREx observations, with
colours matching panel (d). (c) All three optical epochs, from
2018 November 2, 2021 October 14 and 2026 September 24, with labelled positions
of [O II], [Ne III], H-delta, H-gamma, H-beta, [O III], [O I], H-alpha/[N II]
and [S II]. The optical and SPHEREx panels both use observed wavelength in
microns and F-nu in mJy, with linear axes and independent axis ranges.
No offsets, inter-dataset normalisation or host subtraction are applied.
The NGPS spectrum is smoothed for display with a 6-Angstrom FWHM
kernel in observed wavelength, separately within each valid arm segment.
(d) Cleaned SPHEREx aperture measurements show a continuum rising toward longer
wavelengths over 2–5 microns, providing a short-wavelength complement to the
proposed MIRI spectra. The three observing periods are 2025 July 1–14,
2025 December 19–2026 January 6, and 2026 June 27–July 15 (UTC).
Vertical bars are the supplied flux uncertainties and horizontal bars are
the exported spectral half-widths, not wavelength uncertainties. All 289
measurements are plotted individually, without smoothing, rebinning or
averaging across observing periods. Dashed markers locate expected emission
lines at the target redshift; they are not fitted detections.

The standalone `optical_spherex_A3_observed.pdf` / `.png` contains a wider
version of the combined spectrum, retained as supporting material with a
logarithmic wavelength axis. Optical input fluxes in units of
1e-17 erg s^-1 cm^-2 Angstrom^-1 are converted using
`F_nu[mJy] = F_lambda[input] * 1e-17 * lambda_obs[Angstrom]^2 / 2.99792458e18 / 1e-26`.
No factor of (1+z) is applied to flux: the input wavelengths and fluxes are
already observed-frame. This conversion was checked against Astropy's
spectral-density unit equivalency. The plotted spectra have different
apertures, dates and spectral resolutions; peak heights in the overlap
are not directly comparable integrated line fluxes.

## Inputs and remaining scientific checks

- `inputs/figure_targets.json`: existing proposal snapshot, unchanged; all three
  available A3 optical epochs, ZTF g/r, unWISE W1 and NEOWISE W1. The preview
  reads the saved P330E NGPS flux array from the FITS file below.
- `inputs/A3_ngps_flux_standards.fits`: unchanged copy of
  `docs/observed/sep23_p330e/P9694.fits`. The preview uses `FLUX`,
  `WAVE_VAC_HELIO_A` and `MASK`, preserving native channel boundaries.
  `inputs/A3_ngps_flux_standards.json` records its checksum and provenance.
- `inputs/A3_sdss.jpg` and `.json`: SDSS SkyServer DR18 cutout and retrieval URL;
  400 pixels at 0.1 arcsec per pixel (a resampled display image, not a claim of
  0.1-arcsec instrumental resolution). The preview displays its central
  15-by-15-arcsec region; the original cutout is unchanged.
- `spaxel_scryer_P9694_ra11p5333_dec9p1225_cleaned.csv`: user-supplied cleaned
  measurements; `wavelength_um`, `flux_mjy`, `flux_err_mjy`,
  `wavelength_half_width_um` and `mjds` are plotted directly.
- `inputs/A3_spherex_csv_provenance.json`: source checksum, plotted row count,
  exact observing intervals, extraction metadata and boolean flag summary.
- `inputs/A3_spherex_inspector.png` and `.json`: earlier screenshot retained
  for provenance; no longer used by the A3 figure.

Chronological gaps exceeding 45 days separate the three observing periods
(95, 100 and 94 measurements). No additional rows are rejected from the
cleaned export. All four boolean columns for incomplete aperture, incomplete
annulus, asymmetric annulus and saturation/nonlinearity are false. The CSV
records a local-background extraction with an aperture radius of 10.67–10.68
arcsec; this is not a host-subtracted nuclear spectrum. The supplied total
flux-error variance includes a calibration contribution. No independent
assessment of calibration or covariance between different spectral samples
has been performed. An exploratory comparison is recorded in
`spherex_A3_variability.md` and its reproducible script; it finds approximately
4.6% first-to-last brightening over 4–5 microns, with a significance that
depends strongly on the unprovided covariance. No secure variability
detection is claimed. Technical calibration provenance remains in these
notes rather than on the proposal figure.

The SPHEREx panel now also marks expected redshifted emission-line positions;
these are not confirmed detections. The larger `spherex_A3_lines.pdf` makes the
labels easier to read. See `spherex_line_identifications.md` for wavelengths,
sources, and the companion P2190 figure.

The NGPS FITS header identifies P330E as `FLUXSTD` and BD284211 as `ALTSTD`.
The current figure selects P330E (`FLUX`) for both targets, replacing the
previous BD+28 4211 selection. Both arrays remain available in the unchanged FITS files.
`FLUX_ALTSTD` uses the same extraction and telluric correction, with a
wavelength-dependent change to the mean BD+28 4211 sensitivity function
implemented in `scripts/59_ngps_deliver.py`. This is a saved alternative flux
calibration, not an independent reduction of the raw exposures. Median
BD/P330E ratios in the U/G/R/I continuum comparison windows are
0.650/0.670/0.704/0.691, respectively: approximately 30–35% lower. This
comparison alone does not establish which absolute calibration is correct.
The archival spectra retain their original calibrations. The NGPS snapshot
notes that absolute slit losses are not included in its statistical
uncertainties. Check aperture/calibration differences before
claiming a numerical change in broad-line flux from this panel. A host and
continuum decomposition or a justified calibration comparison would strengthen
that claim. Do not infer dust re-formation from a rising continuum alone.

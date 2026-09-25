# SPHEREx line labels and the second example

Both example figures now show expected line positions in their SPHEREx panels.
Both targets are plotted directly from their supplied cleaned CSVs, with
uncertainties and three dated observing periods each.
The wavelengths are calculated from **lambda_observed = lambda_rest × (1 + z)**,
using z = 0.2377224 for A3 / J004607.98+090720.9 and z = 0.739785 for
P2190 / J160101.84+365619.7. The latter coordinate-based name differs by 0.1
arcsec in the final declination digit from the older AllWISE planning table;
the supplied screenshot and saved observer-page target refer to the same source.

## Outputs

- `fig1_single_example_preview.pdf` / `.png`: updated A3 composite.
- `fig1_P2190_example_preview.pdf` / `.png`: companion P2190 composite.
- `spherex_A3_lines.pdf` / `.png`: larger annotated A3 spectrum view.
- `spherex_P2190_lines.pdf` / `.png`: larger annotated P2190 spectrum view.
- `optical_spherex_A3_observed.pdf` / `.png` and
  `optical_spherex_P2190_observed.pdf` / `.png`: wide, combined optical and
  SPHEREx spectra in observed microns and F-nu (mJy).
- `spherex_A3_line_positions.csv` and `spherex_P2190_line_positions.csv`:
  rest and observed vacuum wavelengths, coverage flags and source references.

Run `make single-example-preview` to regenerate everything. Alternatively,
run `python make_fig1_single_example.py --target A3` or `--target P2190`.

## Interpretation

Gold dashed markers locate H-alpha, Pa-delta, He I / Pa-gamma, Pa-beta,
Pa-alpha, Br-gamma and Br-beta where they fall within the nominal 0.75–5-micron
range. P2190 also includes H-beta and the [O III] doublet, which are below the
SPHEREx wavelength range for A3. He I and Pa-gamma are grouped because their
separation is smaller than a nominal SPHEREx resolution element at these
wavelengths. The H-beta / [O III] label groups nearby markers for readability;
it does not assert that all three lines form a single unresolved blend.

The prominent narrow-looking spike in each spectrum is consistent with
H-alpha at its expected observed wavelength. A3 has an apparent excess near
the expected Pa-alpha position, and P2190 has suggestive structure near
Pa-beta and Pa-alpha. These are visual candidate identifications, not fitted
line measurements. No line detection significance, flux or width is measured
by these labels. Exploratory continuum variability comparisons are recorded
separately in `spherex_A3_variability.md` and `spherex_P2190_variability.md`. Hydrogen labels mean
Balmer (H), Paschen (Pa), or Brackett (Br) transitions; He I is neutral helium
and [O III] is forbidden emission from doubly ionized oxygen.

Br-delta is retained in the tables but omitted from the compact plots to reduce
crowding near Pa-alpha. Br-alpha falls just outside the nominal SPHEREx range
for A3 and well outside it for P2190, and is therefore not drawn. No markers
identify the broad continuum rise as an emission line.

P2190 now uses all 364 individual cleaned CSV measurements: 134 from
2025 June 29–August 4, 125 from 2026 January 18–February 10, and 105 from
2026 July 1–August 6. Its supplied errors and spectral half-widths are drawn
directly, with no rebinning or rescaling. The CSV records local-background
aperture photometry with a radius of approximately 9.15 arcsec; all four
boolean aperture/annulus/saturation flags are false. The spectrum is not
host-subtracted. Its earlier screenshot is retained for provenance only.
A3 now uses all 289 individual CSV measurements, separated by observation
date into July 2025, December 2025–January 2026 and June–July 2026. The plot
uses exported flux errors and spectral half-widths, with no rebinning. Its
larger panel also includes a rest-wavelength axis. See
`fig1_single_example_notes.md` for its updated caption and extraction metadata.

## Wavelength references

- [SDSS vacuum spectral-line table](https://classic.sdss.org/dr6/algorithms/linestable.php):
  H-beta 0.486268, [O III] 0.4960295 and 0.5008240, H-alpha 0.656461 microns.
- [STScI NICMOS infrared line lists](https://www.stsci.edu/instruments/nicmos/documents/handbooks/instrument/v5/Appendix_26.html):
  vacuum Paschen, He I and Brackett wavelengths. The old Br-alpha wavelength
  entry is transposed: its listed wavenumber 2467.765 cm^-1 gives
  10000 / 2467.765 = 4.05225 microns, which is the value used here, not 4.5225.
- [IRSA SPHEREx mission characteristics](https://irsa.ipac.caltech.edu/Missions/spherex.html):
  nominal 0.75–5-micron coverage and wavelength-dependent resolving power.

## P2190 companion-figure caption

**J1601+3656 links long-term optical and infrared histories to its current spectra.**
(a) SDSS gri cutout, 15 arcsec across, with a 5-arcsec scale bar.
(b) ZTF g/r and WISE W1 histories;
filled symbols denote the bundled unWISE measurements and open symbols denote
NEOWISE. The shaded interval spans the 2021 June 10–2022 May 29 spectral coadd,
with a dashed line at its representative date. The 2026 NGPS epoch is marked
in orange; the 2004 optical spectrum predates the displayed light curves.
The narrow coloured shaded intervals mark the three SPHEREx observing periods.
(c) All three optical epochs: 2004 June 17, the 2021–2022 coadd and
2026 September 24 NGPS. The optical and SPHEREx panels use observed-frame
microns and F-nu in mJy, with linear axes and independent ranges,
without offsets or host subtraction. Optical
F-lambda is converted to F-nu without renormalising the spectra to SPHEREx.
NGPS is smoothed by 6 Angstrom FWHM in observed wavelength for
display, within each valid arm segment. (d) The three SPHEREx epochs retain
individual flux uncertainties and spectral half-widths, and
expected redshifted line positions indicated by dashed gold markers;
these markers do not assert detections.
See `fig1_single_example_notes.md` for the common unit-conversion formula
and interpretation of spectra with different apertures and resolutions.

P2190 data and cutout were extracted only for this target from `docs/index.html`;
source checksums and scope are recorded in `inputs/P2190_figure_target_provenance.json`.
The numerical SPHEREx provenance is in `inputs/P2190_spherex_csv_provenance.json`.
The unchanged `inputs/P2190_ngps_flux_standards.fits` and its JSON record the
saved P330E calibration (`FLUX`) used for the NGPS preview, as for A3;
these technical details are omitted from the proposal figure.
NGPS slit losses and inter-instrument calibration differences still need
checking before quantifying changes in optical line flux.

P2190 is an illustrative alternative, not an added target in the six-object
proposal table. At its redshift the rest-frame 18-micron silicate feature lies
at 31.32 microns, outside MIRI/MRS. The active proposal and its target list
have not been changed.

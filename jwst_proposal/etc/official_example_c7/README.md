# Official ETC example: workbook 325241, calculation 7

User-downloaded archive: `wb325241_c7_2026-09-26_18.26.09.tar` (retained unchanged in `../downloads/web_exports/`).
SHA-256: `114b2c819e415854b90d4f27f385582961de1a9f84a9e96bee68ee4ac1327b19`

`input.json` is copied byte-for-byte from the archive. Other products remain
inside the archive to avoid duplicating the FITS cubes.

The example is from the Cassiopeia A science-program workbook #26.
Inspected configuration:
- MIRI/MRS, channel 3, short setting.
- FASTR1, full array; 27 groups, 1 integration, 4 exposures.
- An extended Gaussian source with an uploaded Cassiopeia A spectrum.
- Off-scene background strategy (`ifunodoffscene`), aperture size 1 arcsec,
  reference wavelength 12.81 microns.
- A supplied numerical background spectrum, not the background for R06.

This is a configuration reference only. Its exposure result and S/N do not
establish AGN feasibility. Replace the source spectrum, geometry, target
background and extraction strategy before calculating any AGN sensitivity.
No engine/reference-data version is recorded in input.json or the inspected
FITS headers, so the download date alone does not establish the calculation
version. Recalculate with the current ETC before using results.

Next case: R06 / P17881, J083826.50+371906.7 (z=0.2111), following the
provisional scenarios in ../../completion_notes.md. Source photometry,
host contribution and spectral shape still need to be established for that run.

# Selected MIRI/MRS sky fields

All 24 sky fields have been screened against infrared catalogs and inspected in archival images. The selection replaces 23 of the original arbitrary north offsets; R07 retains its 60-arcsec north position. The selected coordinates are in `selected_sky_fields.csv`; image panels are in `sky_field_review.pdf`.

## Screening and review

- Queried AllWISE and 2MASS point-source catalogs through IRSA around each AGN. Catalog entries were retained regardless of quality flags so potentially real neighbors were not filtered out.
- Retrieved WISE W1/W3/W4 FITS images and 2MASS J cutouts with WCS metadata. SDSS cutouts provide an additional optical check where returned successfully; absent SDSS images are explicitly labeled, not treated as empty sky.
- Searched offsets of 60–120 arcsec and retained the original north position only if it passed. The crowded R12 field required an expanded search; the adopted offset is 135 arcsec at PA 175 degrees, measured east of north.
- Required at least 25 arcsec from every catalogued AllWISE and 2MASS point source. The 20-arcsec inspection circle provides a surrounding region larger than the central MRS field. A conservative bright-source screen excludes W1<10 sources within 90 arcsec for the other fields; R12 has a 63.7-arcsec bright-source clearance and was individually inspected.
- Checked the image pixels within the 20-arcsec circle against a local 30–60-arcsec annulus to reject conspicuous positive structure, then visually reviewed all 24 fields in WISE and 2MASS. This image statistic is a ranking/screening aid, not a calibrated detection probability or a diffuse-emission upper limit.
- Individual choices, old-position metrics, candidate rankings and review notes are saved in the per-target JSON files. All downloaded data have adjacent JSON files recording the actual request URL. Candidate rankings, failures and download caches are retained locally and ignored by Git; accepted per-target selections and the assembled review PDF remain versioned; the final `selections.json` contains all 24 accepted fields.

## Notes for reduction and final observing review

The selection avoids catalogued and visually obvious compact infrared contamination at the surveys' resolution and depth; it does not claim zero emission at JWST sensitivity. R12 lies in a crowded region with structured diffuse W3 emission, so its foreground residuals require attention during reduction. The red linear trail in R03's SDSS image appears to be an optical imaging artifact and is absent in the infrared images. F08 has a faint optical object outside the central MRS footprint; it is not an infrared catalog detection within the inspection region.

The sky targets have `Extended=YES`, no acquisition, `PrimaryChannel=ALL`, and four-point `EXTENDED SOURCE` dithers, with the same FASTR1 readout, groups and integration counts as the science observations. The source/sky sequence links are retained. Science target positions and science exposure parameters were not changed.

Scripts: `../sky_review.py` fetches, screens and plots; `../apply_sky_fields.py` applies visually accepted coordinates to the APT file and invalidates old scheduling caches. `sky_review.py --plot-only` regenerates the assembled panels from reviewed per-target selections without downloading or resetting their review status.

References: [STScI MIRI dedicated sky guidance](https://jwst-docs.stsci.edu/jwst-mid-infrared-instrument/miri-operations/miri-dithering/miri-mrs-dedicated-sky-observations), [AllWISE image API](https://irsa.ipac.caltech.edu/ibe/docs/wise/allwise/p3am_cdd/), [2MASS image API](https://irsa.ipac.caltech.edu/ibe/docs/twomass/allsky/allsky/).

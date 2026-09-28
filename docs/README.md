# Public observing website

The single [observer page](index.html) has six tabs: **All · Sep23-targets · Observed · Oct26 · Oct27 · About**. The first five each contain exactly one target-list CSV button, with uppercase headers. About explains the science, selection rationale, observing setup, September experience and JWST coverage.

- **All:** 218 unique catalogue targets, including all 24 JWST proposal objects. The catalogue CSV uses default settings and is not an observing sequence.
- **Sep23-targets:** the 14 planned science targets, with their original exposure settings. Standards and backup notes remain within the tab.
- **Observed:** the 18 targets actually observed, in exposure order; the CSV is the actual observation log. Spectrum FITS downloads remain on the cards.
- **Oct26 / Oct27:** 40 science targets per night, with one NGPS science CSV each. Standards and alternatives appear as notes rather than extra CSV buttons.

The five current tab exports are in [targets/](targets/). Search filters the All table; its CSV always contains the complete catalogue. All six tabs stay on the same webpage.

The orange **Observed · 18** tab ([direct view](index.html#observed)) shows the
September 23 science targets in actual observing order, including backups. Each
card has its P330E-calibrated NGPS 1D spectrum, selectable archival comparisons,
actual exposure details, light curves, imaging and FITS downloads. The 18
requested NGPS products are packaged in `observed/sep23_p330e/`; private archival
spectra remain in local data files.

- [Observer page](index.html) — the single page for the run, public copy: About this run (the science question, who is in the pool, what a spectrum tells us, how the sequence was chosen, then the setting, standards and collapsed calibration and procedure notes), the September 23 sequence with timeline, table and science CSV, and the candidate catalogue with visibility charts, summary table and one card per target (geometry, light curves, public spectra, image, manifold position).

Generator: `scripts/51_observer_page.py`, template `web/observer_page_template.html`, data from `scripts/16_candidate_webpage.py` (payload JSON) and the September packet from `scripts/40_september_observing_packet.py`. The page includes up to two public spectroscopic quasar backups per primary within the September tab, with strict Moon >40°, X <1.5 and r <19 cuts. Each backup matches the primary's full visit and exposure count; literature/Zeltyn regions have first priority. Brightness uses the latest cached ZTF r median, falling back to archival r, and is not a forecast. This is the only observer webpage; private collaboration inputs remain in local data files. Regenerating does not commit or push.

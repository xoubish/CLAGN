# Public observing website

The orange **Observed · 18** tab ([direct view](index.html#observed)) shows the
September 23 science targets in actual observing order, including backups. Each
card has its P330E-calibrated NGPS 1D spectrum, selectable archival comparisons,
actual exposure details, light curves, imaging and CSV/FITS downloads. The 18
requested NGPS products are packaged in `observed/sep23_p330e/`; private archival
spectra remain confined to the local collaboration copy.

- [Observer page](index.html) — the single page for the run, public copy: About this run (the science question, who is in the pool, what a spectrum tells us, how the sequence was chosen, then the setting, standards and collapsed calibration and procedure notes), the September 23 sequence with timeline, table and primary CSV, and the candidate pool for all three nights with visibility charts, summary table and one card per target (geometry, light curves, public spectra, image, manifold position).

Generator: `scripts/51_observer_page.py`, template `web/observer_page_template.html`, data from `scripts/16_candidate_webpage.py` (payload JSON) and the September packet from `scripts/40_september_observing_packet.py`. Both pages include up to two public spectroscopic quasar backups per primary and their CSV, with strict Moon >40°, X <1.5 and r <19 cuts. Each backup matches the primary's full visit and exposure count; literature/Zeltyn regions have first priority. Brightness uses the latest cached ZTF r median, falling back to archival r, and is not a forecast. The local collaboration copy with complete spectra is `data/reselection_2026-09-20/observer_page_local.html`. Regenerating does not commit or push.

# Public observing website

The **October sequences + CSVs** section at the top of the [observer page](index.html#october) contains 40 targets for each of October 26 and 27, with JWST labels, visit times, coordinates, cached brightness and geometry. Downloads include science/standard/backup NGPS lists, detailed timing, coordinates, all 80 targets and the JWST coverage audit. Downloadable files are also in [october/](october/). All CSV column headers and displayed page headings are uppercase; coordinate-only exports intentionally have no header. October uses 12-minute visits including the measured two-minute overhead. Public exports omit private reference-dependent measurements.

The orange **Observed · 18** tab ([direct view](index.html#observed)) shows the
September 23 science targets in actual observing order, including backups. Each
card has its P330E-calibrated NGPS 1D spectrum, selectable archival comparisons,
actual exposure details, light curves, imaging and CSV/FITS downloads. The 18
requested NGPS products are packaged in `observed/sep23_p330e/`; private archival
spectra remain in local data files.

- [Observer page](index.html) — the single page for the run, public copy: About this run (the science question, who is in the pool, what a spectrum tells us, how the sequence was chosen, then the setting, standards and collapsed calibration and procedure notes), the September 23 sequence with timeline, table and primary CSV, and the candidate pool for all three nights with visibility charts, summary table and one card per target (geometry, light curves, public spectra, image, manifold position).

Generator: `scripts/51_observer_page.py`, template `web/observer_page_template.html`, data from `scripts/16_candidate_webpage.py` (payload JSON) and the September packet from `scripts/40_september_observing_packet.py`. The page includes up to two public spectroscopic quasar backups per primary and their CSV, with strict Moon >40°, X <1.5 and r <19 cuts. Each backup matches the primary's full visit and exposure count; literature/Zeltyn regions have first priority. Brightness uses the latest cached ZTF r median, falling back to archival r, and is not a forecast. This is the only observer webpage; private collaboration inputs remain in local data files. Regenerating does not commit or push.

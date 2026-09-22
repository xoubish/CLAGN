# Public observing website

- [Observer page](index.html) — the single page for the run, public copy: About this run (setting, calibrations, standards, procedure, science aim), the September 23 sequence with timeline, table and primary CSV, and the candidate pool for all three nights with visibility charts, summary table and one card per target (geometry, light curves, public spectra, image, manifold position).

Generator: `scripts/51_observer_page.py`, template `web/observer_page_template.html`, data from `scripts/16_candidate_webpage.py` (payload JSON) and the September packet from `scripts/40_september_observing_packet.py`. The local collaboration copy with complete spectra and the backup CSV is `data/reselection_2026-09-20/observer_page_local.html`. Regenerating does not commit or push.

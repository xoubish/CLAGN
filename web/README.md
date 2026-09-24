# Observer page sources

`observer_page_template.html` and `candidate_spectra.js` render the current page.
Run `scripts/16_candidate_webpage.py` to assemble the two candidate payloads,
then `scripts/51_observer_page.py` to build `docs/index.html` (public) and
`data/reselection_2026-09-20/observer_page_local.html` (local collaboration copy).

The public page retains all public primary identities and telescope settings;
private-reference S/N and reference-dependent science fields stay local. Public
primary downloads use sanitized comments. Neither generator publishes.

The orange Observed tab shows the 18 September 23 targets in actual exposure
order, with P330E-only calibrated NGPS spectra, archival overlays and CSV/FITS
downloads. These specific science products are included in both copies at the
PI's request. `scripts/observed_ngps.py` packages them from
`sep23_data/reduction_20260924/products_p330e` into
`docs/observed/sep23_p330e`; its manifest also permits rebuilding without the raw
data. Open either page with `#observed` to go directly to this view.

`candidate_review_template.html` is retained for the local legacy explorer and
upstream packet input. The separate primary-page templates and copies are retired.

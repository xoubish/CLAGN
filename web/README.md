# Observer page sources

`observer_page_template.html` and `candidate_spectra.js` render the current page.
Run `scripts/16_candidate_webpage.py` to assemble the two candidate payloads,
then `scripts/51_observer_page.py` to build `docs/index.html` (public) and
`data/reselection_2026-09-20/observer_page_local.html` (local collaboration copy).

The public page retains all public primary identities and telescope settings;
private-reference S/N, reference-dependent science fields and NGPS spectra stay
local. Public primary downloads use sanitized comments. Neither generator publishes.

`candidate_review_template.html` is retained for the local legacy explorer and
upstream packet input. The separate primary-page templates and copies are retired.

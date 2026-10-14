# Observer page sources

`observer_page_template.html` and `candidate_spectra.js` render the current page.
Run `scripts/16_candidate_webpage.py` to assemble the two candidate payloads,
then `scripts/51_observer_page.py` to build the only observer webpage, `docs/index.html`.
The JSON payloads remain separate data inputs; neither creates another webpage.


`observer_tabs.js` controls five tabs: About (default), Sep23-targets,
Sep23-observed, Oct26 and Oct27. About begins with the full available target
manifold and a collapsible searchable catalogue, then the science notes.
Clicking a point or target name renders one detailed card; never render the
entire catalogue's spectra and charts at once. The catalogue includes targets
without saved manifold projections, which are explicitly counted beside the map.
`scripts/observer_tabs.py` supplies the five dataset CSVs under `docs/targets/`:
the complete catalogue under About and one export per observing tab. Standards
and backups are notes within their night tabs. Old `#all`/`#pool` links redirect
to About. `tests/observer_navigation.cjs` checks the generated page's navigation
and guards against restoring eager rendering of every target card.

The October section is supplied by `scripts/october_web.py` and rendered by
`observer_october.js`. Rebuild its observing packet with
`scripts/53_october_observing_packet.py`, then run `scripts/51_observer_page.py`.
The page shows the 80 selected October targets and CSV download buttons;
the build also writes sanitized files to `docs/october/`. CSV headers
and displayed page headings are uppercase. Coordinate-only files are headerless.
The October 14 review uses 14-minute default visits and a 20-minute, 3×300 s
P8544 follow-up on October 27. Exposure counts, visit lengths, target counts and
spare time are rendered from the packet; the historical September setup remains
separate. `scripts/october_review.py` records the timing and science decisions.

The public page retains all public primary identities and telescope settings;
private-reference S/N and reference-dependent science fields stay local. Public
primary downloads use sanitized comments. Neither generator publishes.

The orange Observed tab shows the 18 September 23 targets in actual exposure
order, with P330E-only calibrated NGPS spectra, archival overlays and CSV/FITS
downloads. These specific science products are included on the page at the
PI's request. `scripts/observed_ngps.py` packages them from
`sep23_data/reduction_20260924/products_p330e` into
`docs/observed/sep23_p330e`; its manifest also permits rebuilding without the raw
data. Open the page with `#observed` to go directly to this view.

The Observed view opens with the supplied Palomar night photograph,
`docs/observed/palomar-sep23.jpg`. The desktop banner links to the uncropped image;
mobile displays its full frame. The page uses relative asset paths.

NGPS traces in archival overlays default to a display-only Gaussian smoothing
with 6 Å FWHM (selectable: off, 3, 6 or 10 Å). Smoothing uses native wavelength
spacing and stops at masked gaps and arm boundaries. It is applied only while
an archival trace is also selected; it is not an instrumental resolution match.
Standalone plots, the embedded native flux arrays and downloadable files retain
their original resolution.

`candidate_review_template.html` is a retained historical source template.
Current builds do not render it; upstream tools read JSON payloads directly.
The separate October review, local observer copy, and legacy explorers are retired.

Public exports omit proposal affiliations, sample IDs and coverage audits.
The source packet retains internal planning metadata; public rows retain target
identities, coordinates, sequence, exposures and operational limits.

`observer_visibility.js` renders the queue-ordered start-window charts for
September 23 and both October nights. `scripts/observer_visibility.py` computes
full-visit maximum airmass and minimum topocentric Moon separation at Palomar,
including overhead and excluding fixed standard visits. The interactive time
cursor compares each target at the same proposed start; it does not reserve
other science visits. Per-target limits come from the September plan and the
October NGPS files. Boundaries are rounded inward on a 30-second grid and checked
independently at 20-second cadence in `tests/test_observer_visibility.py`.

`scripts/observer_names.py` applies selected J-prefixed display names and CSV
identities after geometry generation. Internal keys, old deep links and spectral
asset filenames stay stable; both name forms are searchable in About.

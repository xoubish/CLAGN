# Public observing website

The single [observer page](index.html) opens on **About**, followed by **Sep23-targets · Sep23-observed · Oct26 · Oct27**. About starts with the Palomar telescope photograph, followed by the complete available target manifold, a collapsible searchable catalogue and its CSV, then the science and observing notes. Each of the four night/history tabs contains one target-list CSV button, with uppercase headers.

- **About:** the manifold shows all 199 targets with saved projections; the catalogue includes all 224 targets. Click a point or name to load one detailed card. The catalogue CSV uses default settings and is not an observing sequence.
- **Sep23-targets:** the 14 planned science targets, with their original exposure settings. Standards and backup notes remain within the tab.
- **Sep23-observed:** the 18 targets actually observed, in exposure order; the CSV is the actual observation log. Spectrum FITS downloads remain on the cards.
- **Oct26 / Oct27:** 40 science targets per night, with one NGPS science CSV each. Standards and alternatives appear as notes rather than extra CSV buttons.

The **October 22 update** retains 40 targets per night and all protected observations.
H5MUSES-327 is scheduled on October 26. Four requested galaxies
are scheduled; HUGC-03097 remains in the catalogue because the Moon blocks both
nights. All six requested identities have H prefixes. HIRAS-F23290-1252 uses
4×300 s, 26 minutes and X≤1.5; P8544 retains 3×300 s and 20 minutes. The other
78 visits use 2×300 s and 14 minutes. Galaxy brightness is central optical
aperture photometry, with explicitly approximate exposure scenarios on each card.

The five dataset exports are in [targets/](targets/). Search filters the catalogue under About; its CSV always contains the full catalogue, and the manifold keeps all available target positions visible. All five tabs stay on the same webpage. Legacy `#all` and `#pool` links open About; target links still open the selected card.

Sep23-targets, Oct26 and Oct27 include queue-ordered visibility charts. Yellow
marks full-visit starts at X≤1.5, blue X≤1.8, and purple X≤2 only when permitted
for that target. Every colored start retains Moon separation ≥40° throughout
the visit, including overhead, and avoids the reserved standards. Outlines show
the planned visits. Move the time slider or click the chart to compare swaps;
the right column reports full-visit maximum airmass and minimum Moon separation.
These windows do not reserve other science visits, so swaps still require a
queue-conflict check. Geometry uses a 30-second grid at Palomar.

Selected targets use a J prefix on their existing names (for example JP9506).
The page and CSV downloads use the same names; catalogue search also accepts
the original aliases, and existing target links and spectral filenames remain valid.

The October 22 addition **H5MUSES-327** is a requested comparison target in About,
with its saved-model W1 projection, ZTF/WISE curves, SDSS imaging, October
visibility and archive-review links. Gaia astrometry favors a stellar optical
counterpart; its redshift and AGN identity are unconfirmed. It is scheduled on October 26. Its archive payload and the five galaxy reviews are loaded by
`scripts/additional_targets.py` from the local `data/additional_targets/` cache.

The orange **Sep23-observed · 18** tab ([direct view](index.html#observed)) shows the
September 23 science targets in actual observing order, including backups. Each
card has its P330E-calibrated NGPS 1D spectrum, selectable archival comparisons,
actual exposure details, light curves, imaging and FITS downloads. The 18
requested NGPS products are packaged in `observed/sep23_p330e/`; private archival
spectra remain in local data files.

- [Observer page](index.html) — the single page for the run, public copy: About this run (the science question, who is in the pool, what a spectrum tells us, how the sequence was chosen, then the setting, standards and collapsed calibration and procedure notes), the September 23 sequence with timeline, table and science CSV, and the candidate catalogue with visibility charts, summary table and one card per target (geometry, light curves, public spectra, image, manifold position).

Generator: `scripts/51_observer_page.py`, template `web/observer_page_template.html`, data from `scripts/16_candidate_webpage.py` (payload JSON) and the September packet from `scripts/40_september_observing_packet.py`. The page includes up to two public spectroscopic quasar backups per primary within the September tab, with strict Moon >40°, X <1.5 and r <19 cuts. Each backup matches the primary's full visit and exposure count; literature/Zeltyn regions have first priority. Brightness uses the latest cached ZTF r median, falling back to archival r, and is not a forecast. This is the only observer webpage; private collaboration inputs remain in local data files. Regenerating does not commit or push.

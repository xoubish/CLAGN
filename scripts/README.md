# Pipeline scripts

Run scripts from the project root, using `python scripts/<filename>.py`. Numbered module names are retained because scripts import one another. Root/data paths and subprocess entry points have been updated for this directory.

## Current observing workflow

| Scripts | Purpose |
| --- | --- |
| `15a`–`15c`, `32`, `34`, `36`, `37` | Broad parent, W1 acquisition/projection, rolling spectral identities, and completion pipeline. |
| `05`, `23`, `24`, `33`, `39` | Geometry, neighbour screening, prepared lists, and airmass windows with per-tier S/N for the adopted 16-minute visit (`39`). |
| `15f`, `17`–`19`, `28`, `29`, `35`, `38` | Spectral inventories, archive retrieval, imaging, ZTF, and acquisition audits. |
| `22`, `27`, `30`, `31` | NGPS ETC, Hβ sensitivity (legacy 2×600 s columns feed the frozen review score; `snr300_*` columns give S/N per Å for the adopted 1.5″ / 2×3 / 2×300 s setting), dated evidence, and three-night review. |
| `16_candidate_webpage.py` | Build the single observer page (public `docs/index.html` and the local collaboration copy): candidate pool, September sequence, run information; embeds the packet from `40`. |
| `40_september_observing_packet.py` | Build the September sequence packet: NGPS CSVs, backups, timing, README and run information (pages retired 2026-09-21; run `16` afterwards). |
| `41_september_snr5_plan.py` | One instrument setting for all science targets (1.5″ slit, 2×3 binning, 2×300 s); ETC screen for continuum S/N≥5 per Å with airmass-scaled seeing and a moonlit-sky model per slot; integer-program order that maximises slot S/N (airmass, Moon distance, visibility); protects the eleven previous primaries and fills 16-minute visits from the reviewed pool. |
| `42_science_reselection.py` | Separate projected-parent inventory and dated evidence review; no protected primaries and no telescope-packet writes. |
| `43_public_alert_review.py` | Retrieve/cache public ALeRCE detections; reference-corrected alert photometry remains separate from catalog photometry. |
| `44_science_discussion.py` | Local diagnostic page, qualitative review notes, and provisional matched-control pairs. |
| `45_sep23_decision_page.py` | On-demand decision board (not an observer page): every candidate with a feasible visit in observability order, with geometry, per-slot S/N, manifold, images, light curves and spectra; used to choose the September primaries. |

Typical commands from the project root:

```sh
# Rebuild the explorer from prepared tables.
/opt/anaconda3/bin/python scripts/16_candidate_webpage.py

# Deliberately recalculate the September S/N screen, promotions and geometry-optimised order.
/opt/anaconda3/bin/python scripts/41_september_snr5_plan.py

# Regenerate the packet from that plan, including backup choices.
/opt/anaconda3/bin/python scripts/40_september_observing_packet.py

# Regenerate the packet and fetch its SDSS images when needed.
/opt/anaconda3/bin/python scripts/40_september_observing_packet.py --fetch-images

# Continue the broad-parent search after acquisition; do not start a second copy.
/opt/anaconda3/bin/python scripts/37_finish_three_night_search.py

# Check observing geometry.
/opt/anaconda3/bin/python -m unittest discover -s tests -v
```

Packet generation writes to `observing/sep23/`, plus the public primary pages under `docs/` and `web/`. The completion pipeline refreshes the broader candidate explorer; it does not regenerate the frozen September packet. These commands do not commit or push.

The separate science review is rebuilt with `42_science_reselection.py inventory`, then `42_science_reselection.py evidence`, then `44_science_discussion.py`. The inventory stage replaces its parent snapshot; omit it when reproducing evidence against the existing snapshot. Public alerts are optional: run `43_public_alert_review.py <targets.csv>` with columns `name,ra,dec` before building the discussion page. The current discussion builder expects its alert summary. All review products stay under ignored `observing/science_review/`; they may include private spectra and are not publication artifacts. Spectral screening indices are not line fits or CLAGN classifications. Matching calipers are exploratory and still require state/host/cadence review.

## Earlier work retained for reproducibility

Scripts `00`–`14` contain the original manifold, calibration, selection, and enrichment work. Some remain shared dependencies: for example, `07_make_webpage.py` provides light-curve code to the current pages, and `03d_fetch_spectra.py` provides spectral parsing. Earlier selection/report scripts such as `15d`, `15e`, `20`, `21`, `25`, and `26` are also retained. Their main routines can recreate older selections or reports; use the current workflow above for observing products.

The six historical shell drivers are archived under `archive/2026-09-20/legacy_drivers/`. They preserve the old workflow verbatim and have not been ported to this layout. `13c_window_manifold.py` also retains its original scratch-output convention; it is a research prototype, not an observing entry point.

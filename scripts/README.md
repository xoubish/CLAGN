# Pipeline scripts

Run scripts from the project root, using `python scripts/<filename>.py`. Numbered module names are retained because scripts import one another. Root/data paths and subprocess entry points have been updated for this directory.

## Current observing workflow

| Scripts | Purpose |
| --- | --- |
| `15a`–`15c`, `32`, `34`, `36`, `37` | Broad parent, W1 acquisition/projection, rolling spectral identities, and completion pipeline. |
| `05`, `23`, `24`, `33`, `39` | Geometry, neighbour screening, prepared lists, and airmass/exposure alternatives. |
| `15f`, `17`–`19`, `28`, `29`, `35`, `38` | Spectral inventories, archive retrieval, imaging, ZTF, and acquisition audits. |
| `22`, `27`, `30`, `31` | NGPS ETC, Hβ sensitivity, dated evidence, and three-night review. |
| `16_candidate_webpage.py` | Build public and private candidate explorers. |
| `40_september_observing_packet.py` | Build the curated September sequence, backups, NGPS CSVs, and dark primary pages. |

Typical commands from the project root:

```sh
# Rebuild the explorer from prepared tables.
/opt/anaconda3/bin/python scripts/16_candidate_webpage.py

# Deliberately regenerate the September packet, including backup choices.
/opt/anaconda3/bin/python scripts/40_september_observing_packet.py

# Regenerate the packet and fetch its SDSS images when needed.
/opt/anaconda3/bin/python scripts/40_september_observing_packet.py --fetch-images

# Continue the broad-parent search after acquisition; do not start a second copy.
/opt/anaconda3/bin/python scripts/37_finish_three_night_search.py

# Check observing geometry.
/opt/anaconda3/bin/python -m unittest discover -s tests -v
```

Packet generation writes to `observing/sep23/`, plus the public primary pages under `docs/` and `web/`. The completion pipeline refreshes the broader candidate explorer; it does not regenerate the frozen September packet. These commands do not commit or push.

## Earlier work retained for reproducibility

Scripts `00`–`14` contain the original manifold, calibration, selection, and enrichment work. Some remain shared dependencies: for example, `07_make_webpage.py` provides light-curve code to the current pages, and `03d_fetch_spectra.py` provides spectral parsing. Earlier selection/report scripts such as `15d`, `15e`, `20`, `21`, `25`, and `26` are also retained. Their main routines can recreate older selections or reports; use the current workflow above for observing products.

The six historical shell drivers are archived under `archive/2026-09-20/legacy_drivers/`. They preserve the old workflow verbatim and have not been ported to this layout. `13c_window_manifold.py` also retains its original scratch-output convention; it is a research prototype, not an observing entry point.

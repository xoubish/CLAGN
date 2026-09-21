# Observing files

## September 23, 2026 — first half-night

- **[NGPS primary CSV](sep23/sep23_primaries_ngps.csv)** — eight science targets plus two standards, in observing order.
- **[NGPS backup CSV](sep23/sep23_backups_ngps.csv)** — 24 alternatives, three per primary. Substitute a backup when needed; do not append the entire backup list to the observing sequence.
- [Detailed timing CSV](sep23/sep23_sequence.csv)
- [Dark primary page](sep23/sep23_primaries_local.html) — ordered targets, 40″ images, light curves, and all available spectral dates.
- [Packet notes](sep23/README.md) — windows, settings, caveats, and checks.
- [Full local candidate explorer](../data/reselection_2026-09-20/candidate_review_local.html)

The science sequence uses 2 × 600 s, a 1″ slit, and 2×2 binning. Its existing timing and target choices were preserved during cleanup. Standard exposures are initial settings requiring a saturation check on the night.

## Visibility for all three nights

| Night | Image | PDF |
| --- | --- | --- |
| September 23 | [PNG](visibility/sep23_visibility.png) | [PDF](visibility/sep23_visibility.pdf) |
| October 26 | [PNG](visibility/oct26_visibility.png) | [PDF](visibility/oct26_visibility.pdf) |
| October 27 | [PNG](visibility/oct27_visibility.png) | [PDF](visibility/oct27_visibility.pdf) |

Blue shades show the preferred X≤1.5 windows, extensions to X≤1.8, and fallback windows to X≤2.0. These are availability plots; the September primary packet gives the actual sequence. October observing sequences have not been frozen into equivalent packets.

The files in `visibility/` link to the live generated plots, so rebuilding them updates these views without duplicating files. The packet and these links are local and git-ignored. Public versions of the pages remain under `docs/`.

The September CSVs previously lived in `data/reselection_2026-09-20/sep23_packet/`; their new location is `observing/sep23/`.

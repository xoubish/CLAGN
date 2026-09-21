# Observing files

## Science reassessment — discussion, not a replacement queue

- [Local science review](science_review/index.html) — 38 comparison targets with optical/IR light curves and dated Balmer profiles, including available SDSS-V spectra.
- [Assessment and proposed experiment](science_review/README.md) — findings, selection limits, and provisional allocation of science goals.
- [Discussion CSV](science_review/proposed_science_discussion.csv) — evidence and flags; not an NGPS upload file.

This review reopens the projected parent without protecting previous primaries. It contains collaboration data and stays local. The September telescope files below remain the operational packet pending a science-driven replacement.

## September 23, 2026 — first half-night

- **[NGPS primary CSV](sep23_primaries_ngps.csv)** — 11 science targets plus two standards, in observing order.
- **[NGPS backup CSV](sep23_backups_ngps.csv)** — three replacement choices per primary slot. Names can recur for different slots; use the matching row and skip objects already observed. Do not run the whole backup list.
- [Detailed timing CSV](sep23/sep23_sequence.csv)
- [Dark primary page](sep23/sep23_primaries_local.html) — ordered targets, 40″ images, light curves, and all available spectral dates.
- [Packet notes](sep23/README.md) — windows, settings, caveats, and checks.
- [Full local candidate explorer](../data/reselection_2026-09-20/candidate_review_local.html)

The science sequence uses **2×300 seconds for every primary and backup**, a 1″ slit, and 2×2 binning. Each visit allows ten minutes of integration and ten minutes for acquisition/readout. The original eight primaries remain and three backups have been promoted. ETC screening requires nominal combined continuum S/N ≥5 near Hβ under the archival-brightness/bright-sky assumptions. The buffer runs 23:36–00:08 PDT; P12457 follows at 00:08–00:28 on September 24, at lower airmass. For P8548 and P12457, host light makes the total-continuum estimates overstate AGN-only S/N.

Standards use initial settings of 2×30 s for P330E and 2×5 s for BD+28 4211. Inspect the first exposure in all four channels before repeating and adjust as needed. The [NGPS manual](https://caltechopticalobservatories.github.io/NGPS/users-manual/data-reduction-and-telemetry.html) recommends at least about 1,000 counts, preferably about 10,000, and below roughly 40,000 in the Quicklook spatial-profile display. [Archival V-band brightness comparison](sep23/standards_comparison.csv): the science targets span V≈16.9–18.1; P330E (V=13.01) is about 36–106 times brighter, and BD+28 4211 (V=10.58) about 338–993 times brighter. These are same-band integrated-flux ratios, not detector saturation predictions.

## Visibility for all three nights

| Night | Image | PDF |
| --- | --- | --- |
| September 23 | [PNG](visibility/sep23_visibility.png) | [PDF](visibility/sep23_visibility.pdf) |
| October 26 | [PNG](visibility/oct26_visibility.png) | [PDF](visibility/oct26_visibility.pdf) |
| October 27 | [PNG](visibility/oct27_visibility.png) | [PDF](visibility/oct27_visibility.pdf) |

Blue shades show the preferred X≤1.5 windows, extensions to X≤1.8, and fallback windows to X≤2.0. These are availability plots; the September primary packet gives the actual sequence. October observing sequences have not been frozen into equivalent packets.

The files in `visibility/` link to the live generated plots, so rebuilding them updates these views without duplicating files. The packet and these links are local and git-ignored. Public versions of the pages remain under `docs/`.

The two NGPS CSVs are directly accessible in `observing/`. They link to the canonical files in `observing/sep23/`, so rebuilding the packet keeps them current.

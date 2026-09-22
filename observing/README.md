# Observing files

## Science reassessment — discussion, not a replacement queue

- [Local science review](science_review/index.html) — 38 comparison targets with optical/IR light curves and dated Balmer profiles, including available SDSS-V spectra.
- [Assessment and proposed experiment](science_review/README.md) — findings, selection limits, and provisional allocation of science goals.
- [Discussion CSV](science_review/proposed_science_discussion.csv) — evidence and flags; not an NGPS upload file.

This review reopens the projected parent without protecting previous primaries. It contains collaboration data and stays local. The September telescope files below remain the operational packet pending a science-driven replacement.

## September 23, 2026 — first half-night

- **[Decision board](sep23/sep23_decision_local.html)** — all 85 candidates that fit a full 16-minute visit on the night, numbered P01–P85 in observability order, with the visibility chart, a comparison table, and per-target airmass, Moon, brightness, per-slot S/N, manifold, image, light curves and spectra. Companion tables: [summary CSV](sep23/sep23_decision_table.csv), [per-slot S/N CSV](sep23/snr5_slot_table.csv). Local only; it contains collaboration spectra.
- **[NGPS primary CSV](sep23_primaries_ngps.csv)** — 14 science targets plus two standards, in observing order.
- **[NGPS backup CSV](sep23_backups_ngps.csv)** — up to three replacement choices per primary slot. Names can recur for different slots; use the matching row and skip objects already observed. Do not run the whole backup list.
- [Detailed timing CSV](sep23/sep23_sequence.csv) · [Per-slot S/N table](sep23/snr5_slot_table.csv)
- [Dark primary page](sep23/sep23_primaries_local.html) — ordered targets, 40″ images, light curves, and all available spectral dates.
- [Packet notes](sep23/README.md) — windows, settings, caveats, and checks.
- [Full local candidate explorer](../data/reselection_2026-09-20/candidate_review_local.html)

One instrument setting for every science target and both standards: **1.5″ slit, 2×3 binning (spatial × spectral), 300 s sub-exposures**, 2 per target except the two faint favourites P1823 (4×300 s) and P8544 (3×300 s). Decided on 2026-09-21 from an ETC slit scan under the 93 % Moon: at the 1.5–1.8″ seeing expected at the slit, 1.5″ recovers 30–37 % of exposure time relative to the earlier 1.0″ slit, R≈1,650 is ample for Balmer-line work, and 2×3 is the documented binning for that slit. A standard visit is 16 minutes (10 min integration plus 6 min for slew, acquisition and the readout between exposures); the deeper visits are 27 and 22 minutes. Unscheduled time inside the science block is 11 minutes in total, 8 of them before the closing standard (00:20–00:28 PDT); the sequence is packed, so any delay beyond that eats into the second half-night.

The 14 primaries are the PI's choice from the [decision board](sep23/sep23_decision_local.html) on 2026-09-21 (board numbers P3, P1, P5, P8, P11, P22, P32, P26, P38, P37, P46, P52, recorded in [user_selection.json](sep23/user_selection.json)) plus P9227 and P9506, filled automatically. Three of the picks are comparison reserves observed as controls (P9584, P9694, P9599) and two sit below the S/N floor even with extra exposures (P1823, predicted 2.9 per Å at its slot; P8544, 3.7). The previous primaries P3642, P8548, P11113, P11010, P10596 and P12457 are now backups. The order is set by an integer program that maximises each target's predicted continuum S/N per Å near observed Hβ in its slot, with the five setting RA 16h targets forced into the first 90 minutes by their windows; the S/N screen uses the official NGPS ETC with seeing scaled to the target airmass and a Krisciunas & Schaefer moonlight sky at each mid-visit.

Standards use initial settings of 2×60 s for P330E and 2×10 s for BD+28 4211 with the same slit and binning (sized for the 1.5″ slit to reach roughly 5,000–8,000 peak counts per binned pixel), so the afternoon calibrations are one set: 3 ThAr and 3 FeAr arcs and 7 biases per channel at 2×3 binning, and at least 5 dome flats per channel (7–10 for U and G) at 1.5″ + 2×3. Inspect the first exposure in all four channels before repeating and adjust as needed. The [NGPS manual](https://caltechopticalobservatories.github.io/NGPS/users-manual/data-reduction-and-telemetry.html) recommends at least about 1,000 counts, preferably about 10,000, and below roughly 40,000 in the Quicklook spatial-profile display. [Archival V-band brightness comparison](sep23/standards_comparison.csv): the science targets span V≈16.9–18.1; P330E (V=13.01) is about 36–106 times brighter, and BD+28 4211 (V=10.58) about 338–993 times brighter. These are same-band integrated-flux ratios, not detector saturation predictions.

## Visibility for all three nights

| Night | Image | PDF |
| --- | --- | --- |
| September 23 | [PNG](visibility/sep23_visibility.png) | [PDF](visibility/sep23_visibility.pdf) |
| October 26 | [PNG](visibility/oct26_visibility.png) | [PDF](visibility/oct26_visibility.pdf) |
| October 27 | [PNG](visibility/oct27_visibility.png) | [PDF](visibility/oct27_visibility.pdf) |

Blue shades show the preferred X≤1.5 windows, extensions to X≤1.8, and fallback windows to X≤2.0. These are availability plots; the September primary packet gives the actual sequence. October observing sequences have not been frozen into equivalent packets.

The files in `visibility/` link to the live generated plots, so rebuilding them updates these views without duplicating files. The packet and these links are local and git-ignored. Public versions of the pages remain under `docs/`.

The two NGPS CSVs are directly accessible in `observing/`. They link to the canonical files in `observing/sep23/`, so rebuilding the packet keeps them current.

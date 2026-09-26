# 24-AGN MIRI program — 26 September 2026

Open **clagn24_miri_draft.aptx** in APT 2026.5.1 or later. APT 2026.5.1 has successfully loaded this file, run STScI's online Visit Planner, exported timing reports and saved the scheduling results. It contains 24 AGNs (12 fading, 12 rising), plus 24 background fields, in 48 observations linked as 24 consecutive, noninterruptible science/sky pairs. It is a draft, not a submitted or fully validated proposal.

## Current budget and configuration

The final online planning run gives **67.33 charged hours** in the timing JSON.
APT rounds its saved proposal-summary request upward to **67.4 hours**, matching
the proposal text. The program remains Medium. All **48 visits** have current,
nonempty scheduling windows from constraint generator 19.0.1.

Source-plus-sky photon collection totals **17.29 hours**; APT's science-duration
accounting gives 17.33 hours. Acquisition and operational overheads are included
in the charged total. All four MRS channels are simultaneous.

| Setting | Default groups | Integrations per dither | Dithers | On-source time |
| --- | ---: | ---: | ---: | ---: |
| A | 27 | 1 | 4 | 299.7 s |
| B | 60 | 1 | 4 | 666.0 s |
| C | 27 | 1 | 4 | 299.7 s |

- R06 uses **80 groups in B (888.0 s)**; R01/P1823 uses **75 groups in A (832.5 s)**.
  Each matching sky uses the same exposure sequence.
- MIRI/MRS FULL/FASTR1 on both detectors, all channels, four-point point-source
  science dithers, no simultaneous imaging. Background targets are extended,
  use four-point extended-source dithers, and have no acquisition.
- Science acquisitions use F560W/FAST/4 groups for four sources,
  FND/FAST/10 groups for eight, and FND/FASTGRPAVG/10 groups for twelve.
  Recipes and tested flux brackets are in `../review/acquisition_recipes.csv`
  and `../review/README.md`. FASTGRPAVG coadds four frames per group.
- Selected sky positions and sample coordinates are unchanged by this review.
  See `sky_fields/README.md` for the completed catalog/image checks.
- Remaining source-model assumptions and the nine extended-host acquisition
  centering checks are explicit in the technical review. Individual nuclear
  fluxes and sample-wide sensitivity have not all been measured/validated.
- Administrative/PI fields are blank. The proposal has not been submitted.

The preceding sky-reviewed budget was 66.59 hours. The technical fixes add
0.74 charged hours. Earlier versions remain under `work/before_technical_review/`,
`work/before_sky_review/` and `work/pre_visit_planner/`.

## Files

- `clagn24_miri_draft.aptx`: APT-processed review file, schema 71, APT 2026.5.1, bundled PRD `PRDAPTSOC-074-002`.
- `clagn24_miri_draft.xml`: matching XML extracted from that file.
- `clagn24_miri_draft.times`: APT's per-observation exposure and overhead report.
- `clagn24_miri_draft.timing.json`: APT's timing export, including the program total.
- `clagn24_miri_draft.pointing`: APT's pointing and dither report after sky selection.
- `sky_fields/sky_field_review.pdf`: annotated archival images of all selected sky fields.
- `sky_fields/selected_sky_fields.csv`: coordinates, offsets and catalog clearances.
- `observation_inventory.csv`: compact source/sky and exposure inventory.
- `build_draft.py`: reproducible input generator. Defaults to `work/generated/` so it does not overwrite the APT-processed review file or its timing reports.
- `runtime/`: locally unpacked official APT installer; the older application in `/Applications` was not replaced. Download SHA-256: `28411dcd58f3560cb7a35a73d849fbd0817d2172479f2c07fe1a078d5c311b28`.

The user explicitly approved transmitting the draft and target coordinates to STScI's planning services. The online Visit Planner then completed successfully. **The proposal has not been submitted.**

Sources: [official APT installer](https://apst.stsci.edu/apt/external/downloads/APT-2026.5.1/Web/install.html), [Cycle 6 GO size categories](https://jwst-docs.stsci.edu/jwst-opportunities-and-policies/jwst-call-for-proposals-for-cycle-6/jwst-proposal-types-and-categories/jwst-general-observer-go-proposals). The XML representation was informed by public MIRI program 1523, with no investigator information, program constraints or cached timings inherited.

# JWST duplication screening — 2026 September 26

All 24 science coordinates returned zero JWST observation matches in a
120-arcsec search. All queries completed successfully. No potential MIRI/MRS
duplications were identified in this archive snapshot.

The search follows the [STScI MAST duplication notebook](https://spacetelescope.github.io/mast_notebooks/notebooks/JWST/duplication_checking/duplication_checking.html),
using `astroquery.mast.Observations.query_criteria` with `obs_collection='JWST'`.
There are no instrument, date, public-data, or calibration-level restrictions;
planned records (`calib_level=-1`) are included by the query. Coordinates come
from `ra_hms` and `dec_dms` in `../inputs/jwst_sample_cycle6.csv`, matching the
science positions used in APT. The search radius exceeds the MRS field size.

- `target_screening.csv`: per-target coordinates, query times, counts and results.
- `screening_summary.json`: aggregate result, input checksum and review status.
- `matches.json`: returned science-target records (empty).
- `cache/`: original ECSV responses and per-target query metadata, retained locally
  and ignored by Git. The compact dated audit tables above remain versioned.
- `positive_control.json`: the same search at NGC 7469 returned 71 observations
  from program 1328, including MIRI/IFU, confirming that the service returned
  known JWST observations. This control is separate from the 24-target sample.

With no spatial matches, no program-by-program configuration or exposure-time
comparison was needed. This records the MAST contents on the query date;
approved footprints can be approximate, and new programs can be added later.

To reproduce the saved result, run
`/opt/anaconda3/bin/python jwst_proposal/duplication/check_targets.py` from the
project root. The script reuses cached per-target responses. For a new archive
snapshot, preserve this audit directory and remove the 24 science-target ECSV
cache files before rerunning; review the new output before updating the text.

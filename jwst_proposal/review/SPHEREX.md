# SPHEREx retrieval and local data

All 24 positions were queried against the public QR2 collection on 2026-09-29.
Retrieval for the 21 targets beyond the three existing Figure 1 examples is
complete: **6081 image outcomes, 5194 accepted and 887 quality-rejected, with
zero remaining download errors**. Usable spectra exist for 23 targets overall;
R12 is crowded and has no accepted extraction. Rejection is not a nondetection.

## Storage

The Git-ignored cache is `jwst_proposal/local_data/spherex/sample/`, containing
query XML, image-index CSVs, per-image JSON results, extracted tables, logs,
progress and background-worker state. It was moved intact from
`review/spherex_sample/` during the September 29 cleanup. URLs and per-image
hashes are retained. Temporary cutouts are normally deleted after extraction;
the three downloaded diagnostic FITS now live in
`jwst_proposal/local_data/spherex/rejection_checks/`.

Curated Figure 1 inputs remain in `inputs/spherex_extracted/` and the supplied
R01 CSV. Review outputs remain in [spherex_sample_review/](spherex_sample_review/README.md).
The original three spectra were not replaced. The compact R12 diagnostic report
is `spherex_sample_review/rejection_checks/diagnostics.json`.

## Status and restart

From the repository root:

```sh
/opt/anaconda3/bin/python jwst_proposal/review/spherex_background.py status
/opt/anaconda3/bin/python jwst_proposal/review/spherex_background.py start
```

There is no need to restart the completed retrieval. If it is needed later,
`start` launches a detached worker with automatic retries and a duplicate-worker
lock. It reuses accepted and quality-rejected per-image records and retries
failed or truncated records. `stop` saves the current batch before stopping.
After reboot, use `start` again; no login service is installed. Computer sleep
pauses downloads. Do not run the direct extractor alongside the worker.

The direct equivalent is:

```sh
/opt/anaconda3/bin/python jwst_proposal/review/retrieve_spherex_sample.py --ids F01 F02 F03 F05 F07 F08 F09 F10 F11 F12 R02 R03 R04 R05 R06 R07 R08 R09 R10 R11 R12
```

`progress.json`, `extraction_status.csv`, and individual JSON outcomes record
completion; image-index rows alone do not establish usable spectral coverage.
The extractor reproduces the existing F04 measurements for the eight trial
images checked. QR2 was selected for consistency with Figure 1; QR3 would need
separate extraction validation.

[IRSA cutout documentation](https://irsa.ipac.caltech.edu/data/SPHEREx/docs/cutout_tool.html).

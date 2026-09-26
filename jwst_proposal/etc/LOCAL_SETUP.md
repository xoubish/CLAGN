# Local MIRI/MRS ETC

Installed and validated on 2026-09-26. Current technical results are in `../review/README.md`.

## Installed components

- Isolated Python 3.12 environment: `.venv/`, with `pandeia.engine==2026.7`.
  Exact dependencies are recorded in `requirements-local.lock`.
- Full JWST instrument reference data, version 2026.7:
  `reference_data/pandeia_data-2026.7-jwst/`.
- All 104 MIRI/MRS channel 1–4 PSFs from the official release, plus metadata:
  `reference_data/mrs_psfs_2026.7/`. The archive's internal directory name
  includes `2026.7rc1`, but its VERSION_PSF explicitly identifies 2026.7.
- Official Synphot WISE W1–W4 bandpasses in
  `reference_data/synphot/comp/nonhst/`, downloaded individually. The Vega
  spectrum from Pandeia's reference data is in the local Synphot calspec folder.

This supports our MIRI/MRS scenarios, not every JWST mode or template library.
An additional 58 complete MIRI imager PSFs now support the dedicated target-acquisition
runner in `../review/select_acquisition.py`; their checksums and FITS verification
are in `../review/ta_psf_manifest.json`. Other modes, photometric systems and
external spectral templates require their own reference files. Full PSF and Synphot atlas
downloads were stopped after obtaining the required components; incomplete
downloads were explicitly named `.partial`; those unused partial archives were
deleted during cleanup after verifying the installed reference-file checksums. The complete MRS files were read
from the beginning of the official PSF archive, before the next aperture's
files. No truncated FITS files are used. All 109 installed PSF/bandpass/Vega
FITS files passed verification; hashes are in `installed_reference_manifest.json`.

Environment, downloads, reference files and regenerable runs are ignored by Git.

## Running a calculation

From the workspace root, select a new output directory for each trial:

```sh
jwst_proposal/etc/.venv/bin/python jwst_proposal/etc/run_local.py \
  jwst_proposal/etc/r06_ch4b_60groups_r09/input.json \
  jwst_proposal/etc/local_runs/my_trial
```

The runner sets process-local paths automatically; no shell-profile edits
are needed:

| Variable | Directory relative to this file |
| --- | --- |
| `pandeia_refdata` | `reference_data/pandeia_data-2026.7-jwst` |
| `PSF_DIR` | `reference_data/mrs_psfs_2026.7` |
| `PYSYN_CDBS` | `reference_data/synphot` |

It saves inputs, scalar results, warnings, versions, 1-D FITS products and a
`local_export.tar` compatible with `bin_archive.py`. It refuses non-MRS modes
and existing output folders. It adds a primary HDU to Pandeia's table-only
FITS HDU lists without changing numerical data. It configures the local Vega
file explicitly, avoiding an irrelevant missing-default-Vega import warning.

## Web validation

The unchanged input from `downloads/web_exports/wb325241_c7_2026-09-26_19.17.43.tar` reproduces
the 60-group channel 4B web calculation. Results are in `local_validation.json`;
regenerable products are in `local_runs/validation_ch4b_60groups/`.

- All 916 spectral samples match for extracted flux, noise and S/N, with
  maximum relative difference below 8e-16.
- Native scalar S/N: 6.80178; web table: 6.80.
- On-source time: 666.0096 s; strategy time including sky: 1332.0192 s.
- Calculation warnings: none.
- Exported local FITS products pass the binning script; R~100 S/N agrees
  with the downloaded web output to below 2e-14 in the comparison bins.

The scalar uses the nearest pixel, 21.80115 um. Interpolating at exactly
21.80 um gives a slightly different native S/N; this is not a local/web
discrepancy. Optional argument `--compare-web ARCHIVE` validates another run
against identical web inputs and full signal/noise/SNR arrays.

Engine validation does not establish the assumed source spectrum, host noise,
spectral covariance, calibration errors or final observing times. Existing
numerical backgrounds apply to this source; other targets need appropriate
background spectra of their own.

## Official sources

- https://outerspace.stsci.edu/spaces/PEN/pages/77530136/Pandeia%2BEngine%2BInstallation
- https://stsci.box.com/v/pandeia-data-v2026p7-jwst
- https://stsci.box.com/v/pandeia-psfs-v2026p7-jwst
- https://ssb.stsci.edu/trds/comp/nonhst/wise_w3_001_syn.fits
  (W1, W2 and W4 use the same URL pattern.)
- Full Synphot base atlas:
  https://archive.stsci.edu/hlsps/reference-atlases/hlsp_reference-atlases_hst_multi_everything_multi_v16_sed.tar

## Local archives and Git

Original web exports are retained in ignored `downloads/web_exports/`. Generated
ETC arrays and caches stay local; compact inputs, results and validation manifests
remain eligible for Git. Duplicate `local_export.tar` bundles were removed only
after comparing every member with its retained loose file. A fresh `run_local.py`
run writes a new export when needed.

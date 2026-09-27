# Focused technical and editorial review — 26 September 2026

The proposal and saved APT file now request **67.0 hours**; the timing export
contains **67.00 charged hours**. APT rounds the proposal summary upward to a tenth.
The program has 24 MIRI/MRS targets (12 fading, 12 rising), 24 matching sky visits,
and current nonempty scheduling windows for all 48 visits. Source-plus-sky photon
collection is 17.00 hours; APT science-duration accounting is 17.03 hours.
The official style is unchanged; the final PDF has five core pages plus two
supplement/reference pages.

## Changes justified by the calculations

1. R06/J083826.50+371906.7: increased setting B from 60 to 80 groups in both
   science and sky observations. Rest-5-micron S/N rises from 37.9 to 54.8.
2. R01/P9694 replaces P1823 and uses the default 27/60/27 sequence. Its new sky field, acquisition recipe and online scheduling have been checked.
3. Replaced uniform F560W/FAST/10-group acquisitions, which saturate in the bright
   model cases, with flux-dependent recipes in `acquisition_recipes.csv`.
4. Updated text, APT inventory, reproducible generator and time request; updated
   R01 sample membership, sky coordinates and Figure 1 to P9694.

## MRS sensitivity and saturation

Local Pandeia engine/reference data/PSFs are version 2026.7. The engine previously
reproduced the supplied web calculation to floating-point precision. Inputs,
reports, spectral arrays and R=100 calculations are under the Git-ignored local
`technical_runs/` directory. Scripts and compact results remain versioned;
`sensitivity_summary.json` retains the quoted sensitivity results.

| Source / case | Rest wavelength | Groups | Radius | R=100 statistical S/N |
| --- | ---: | ---: | ---: | ---: |
| R06 | 5 microns | B, 80 | 0.3 arcsec | 54.8 |
| R06 | 12 microns | B, 80 | 0.6 arcsec | 225.2 |
| R06 | 18 microns | B, 80 | 0.9 arcsec | 61.5 |
| P9694 | 5 microns | B, 60 | 0.3 arcsec | 336.0 |
| P9694 | 12 microns | B, 60 | 0.6 arcsec | 1343.6 |
| P9694 | 18 microns | B, 60 | 0.9 arcsec | 472.2 |

R06 uses z=0.2111, a point source with Fnu proportional to wavelength squared,
and W3=0.95 mJy (half its archival total). P9694 uses z=0.2377224,
W3=11.76 mJy (half its archival total) and the same wavelength-squared shape.
Both use the adopted numerical sky spectrum from the R06 web calculation,
with equal-duration off-source subtraction. These are representative-background
estimates. `p9694_sensitivity.json` retains full compact inputs and anchor results;
its maximum tested saturation fraction is 0.022, with no warnings. Historical
P1823 calculations remain in `sensitivity_summary.json` but no longer justify
an exposure or science claim for the current sample.

The bright-envelope B-setting checks use W3=48.4 mJy, flat Fnu for channels 1–2
and wavelength-squared spectra for channels 3–4. Maximum saturation fractions
are 0.100/0.125/0.086/0.031, with no warnings. This is a focused test of the
longer default setting, not a calculation of all twelve subbands for every target.

Binning sums signal and adds sample variances, without spectral covariance.
Nuclear fractions and shapes are planning assumptions; extended-host photon
noise and calibration systematics are not included in these ETC estimates.
The science analysis propagates host and calibration uncertainties separately.
S/N near the specified continuum wavelengths is not a guarantee throughout each
feature, or a detection significance for model evolution. Source histories and
individual nuclear flux fractions remain to be reconciled with the full data.

## Target acquisition

The tested lower-flux model is half each source's archival W3 with Fnu proportional
to wavelength squared; the upper model is its full W3 with flat Fnu. These are
assumed planning brackets, not measured bounds on JWST-epoch nuclear flux.
The saved recipe passes both cases without saturation warnings and with S/N>50:

- Three targets: F560W, FAST, 4 groups (11.10 s).
- Nine targets: FND, FAST, 10 groups (27.75 s).
- Twelve targets: FND, FASTGRPAVG, 10 groups (111.00 s; four coadded frames/group).

Across the selected recipes, the minimum lower-case S/N is 51.79 and maximum
upper-case saturation fraction is 0.567. Full inputs/results are in
`technical_runs/ta_grid/` and `acquisition_recipes.csv`. The 58 complete imager
PSFs used were extracted from the available official archive, read and FITS-verified;
`ta_psf_manifest.json` records their hashes. No truncated FITS files are used.

**Remaining centering check:** AllWISE marks F08, F09, F10, F11, F12, R08, R09,
R10 and R12 as extended. Host extension does not rule out acquisition on a compact
nucleus, but higher-resolution morphology or nuclear dominance should be reviewed
before final submission. Catalog positions agree within 0.68 arcsec for all 24.
The only listed neighbour within 5 arcsec is at 4.98 arcsec for R08, outside the
centered 5-by-5-arcsec acquisition region; unresolved structure is not excluded.
See `acquisition_morphology_screen.json`. Recipe S/N tests do not establish
centroid suitability.

## References, anonymity and format

All 35 numbered entries were checked against primary publication metadata or
STScI/IPAC documentation, with links retained in `reference_audit.csv`. Minezaki
et al. (2019) now has its published journal reference. Ansky is described as
recently brightened, consistent with evidence for prior nuclear activity.
This bibliographic check is not a full replication of the cited science.

The PDF text, metadata, figures and APT title/abstract were screened for identity
leaks and NGPS observing-status/month/program details. The public wording remains
“2026 Palomar/NGPS spectra.” Published author lists in references are retained. No
investigator names, affiliations or email addresses were introduced into the
anonymous attachment or abstract. The APT administrative fields still need filling.

Restored the official template-version statement and added a factual AI-use
disclosure in the references, as required by the Cycle 6 PDF instructions. Both
are outside the five-page required sections. The style file was not modified.

## Reproduction and audit

- `run_focused_etc.py`: focused MRS inputs and runs.
- `select_acquisition.py`: model bracket tests and acquisition selection.
- `apply_technical_updates.py`: applies the selected recipes; strips cached tool
  data, so rerun APT afterwards if using this script again.
- `audit_consistency.py` / `consistency.json`: saved APT/XML, inventory, recipe,
  scheduling-window, abstract and layout comparisons.
- `apt_recompute.log`: completed APT 2026.5.1 online run; constraint generator 19.0.1.
- Earlier APT files are preserved under `../../archive/jwst_cleanup_20260926.zip::apt/work/before_technical_review/`.

Official guidance consulted:
- [MRS target acquisition](https://jwst-docs.stsci.edu/jwst-mid-infrared-instrument/miri-operations/miri-target-acquisition/miri-mrs-target-acquisition)
- [ETC MIRI target acquisition](https://jwst-docs.stsci.edu/jwst-exposure-time-calculator-overview/jwst-etc-calculations-page-overview/jwst-etc-target-acquisition/jwst-etc-miri-target-acquisition)
- [Cycle 6 PDF requirements](https://jwst-docs.stsci.edu/jwst-opportunities-and-policies/jwst-call-for-proposals-for-cycle-6/jwst-preparing-the-proposal-pdf-attachment)
- [Dual-anonymous review](https://jwst-docs.stsci.edu/jwst-opportunities-and-policies/jwst-call-for-proposals-for-cycle-6/jwst-key-policies/jwst-dual-anonymous-peer-review)

No proposal has been submitted. Focused checks passed under explicit assumptions;
this review does not claim full-sample, every-wavelength ETC or astrophysical validation.

## Workspace cleanup

`cleanup_manifest.json` records removed duplicate exports/unused downloads and
relocated web exports/automatic backups. Original web ETC archives now live in
`../etc/downloads/web_exports/`; APT snapshots live in `../../archive/jwst_cleanup_20260926.zip::apt/work/automatic_backups/`.
These are local-only. Generated ETC grids and raw archive-query caches are
ignored by Git, while recipes, scripts and summary results are retained.

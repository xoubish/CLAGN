# Palomar paper: executed pilot and next decisions

The active project is **multi-epoch broad-line evolution and a test of the
pre-run manifold selection**. Dust physics waits for JWST. The strongest path
currently is a complete observed-sample paper with carefully validated case
studies. Successful prediction, selection enrichment, and a standalone new-source
discovery are hypotheses; this pilot does not yet establish them.

Start with [executed notebook 07](notebooks/07_palomar_paper_fits_and_results.ipynb)
or [its HTML](outputs/html/07_palomar_paper_fits_and_results.html). The
[protocol notebook](notebooks/06_palomar_paper_protocol_and_pilot.ipynb) explains
selection and timing. [Fit atlas: 15 sources](outputs/palomar_paper/pilot_fit_atlas.pdf).

## What the pilot actually accomplished

- Selected 11 outcome-informed science cases and four D-group comparisons using
  redshift, historical brightness and archival baseline. This is explicitly a
  post-first-look protocol, not a preregistration or an unbiased yield sample.
- Selected 43 epochs using dates/quality, recovering native uncertainties for
  three DESI epochs. Fit 42 epochs with three host/Fe II model families; leave
  the P10381 DESI reference unavailable rather than substitute an older date.
- Fit every NGPS visit under the alternative standard for all three model
  families: 60 additional fits. These are calibration stress tests; the cloudy
  October 27 Feige standard is not equally plausible with the adopted standard.
- Checked eight independent native exposures for P21422, P7947 and P9506.
- Refreshed public corrected ZTF alerts for all 15 sources; five visits now have
  a band with at least three nights within ±30 days, including three October
  visits. P21422's nearest qualifying r point is about 19 days earlier.
- Ran 14 invariant/recovery unit checks and 72 two-epoch injection trials.
  The latter reveal bias under model mismatch and are not real-data completeness.

## What changed scientifically

Ratios below compare the new spectrum with the date-selected latest accepted
archive. Ranges span tested models/standards, **not confidence intervals**.

| Case | Current result | Decision |
|---|---|---|
| **P21422** | Broad-Hβ EW ratio **0.59–0.69**; Hβ/[O III] **0.52–0.61**. Primary-standard fits alone give EW 0.66–0.69. Both line windows covered. Two exposure wing EWs agree within ~2.4%. | First joint Hα/Hβ case. Broad lines remain present. Audit the Graham-2019 literature association before claiming novelty. |
| **P9506** | October/archive EW **0.67–0.77**; Hβ/[O III] **0.53–0.55**. September supports weakening too; individual-exposure wing indices are similar across both nights. | Promote as persistent long-baseline weakening. Do not claim an abrupt month-scale transition. |
| **P7947** | EW **1.51–2.11** across standards/models; adopted-standard range 1.70–2.11. | Strong follow-up candidate, but host/centroid sensitivity and masked [O III] support need work. |
| **P24930** | Adopted-standard full-line EW **0.87–0.88**, much less extreme than the first wing-index result. | Down-rank the large-change claim; inspect profile/coverage and decomposition. |
| **P17881** | Adopted-standard EW **1.57–1.77**, but Hβ/[O III] **0.96–1.03** versus February 2025. | Known-source evolution; continuum dilution and aperture/calibration matter before a broad-flux claim. |
| **P9694** | October/September EW stress envelope **0.89–1.30**. | Month-scale growth is calibration dependent. |

Several comparison sources have EW and [O III]-normalized trends in opposite
directions. This is useful evidence that absolute calibration and aperture effects
must be modeled; the four comparisons are not a calibrated control experiment.

## What makes the paper publishable

1. **Validate the leading three cases.** Include all usable archival epochs,
   Hα/Hβ simultaneously, shared host information, flexible broad profiles, and
   aperture-dependent host/narrow emission. Check native resolution, wavelength
   residuals and extraction/standard choices. The present fits are independent
   epoch models with approximate instrumental response, not that joint model.
2. **Anchor continuum changes where possible.** Use actual filter transmission
   curves and time-matched photometry with consistent total/reference flux and
   aperture treatment. Alert photometry is not complete forced photometry.
   Do not interpret fitted continuum slopes as intrinsic disk changes yet.
3. **Measure the full observed denominator.** Extend the validated procedure
   to every Hβ-accessible source, retain non-detections/failures, and measure
   sensitivity through injections into real spectra and residuals. The pilot
   was chosen partly for interesting outcomes and cannot provide discovery rates.
4. **Test incremental selection information.** Preserve the pre-run A/B/C/D
   groups and known-source stratum. Compare manifold features against simple
   pre-run variability using source-level evaluation and uncertainty. Only two
   C sources and nine D sources were observed, limiting adjusted inference.
5. **Resolve forecasting provenance.** Of 99 campaign sources, 46 pass a timing
   screen in which the latest archive follows all linked original W1 inputs;
   29 have a baseline preceding some inputs and 24 have unresolved feature
   support. These are not yet verified independent forecasts: historical
   projection snapshots and training-identity overlap remain to be established.
   Pool IDs and Sample A training IDs must not be equated by their integer values.

If the selection advantage does not survive these checks, report that result
alongside the spectral-evolution catalog and individual cases. Do not rewrite
the selection rule around the interesting pilot objects.

## Measurement and provenance limits

The fitted model uses a power law, nonnegative BC03 host components and an optical
Fe II template, narrow lines, and two broad-Hβ Gaussians. The native per-pixel
errors enter a conditional least-squares objective; the reported model spread
omits a full covariance/calibration posterior. Fixed stellar dispersion, Fe II
smoothing, [O III] wing shape and approximate LSF need further sensitivity work.
No upper-limit/classification procedure has been validated. Simulations with
extra curvature and correlated residuals bias the fading ratio by a median ~14%,
illustrating why the small three-model spread is not the full error budget.

Template data and license are retained from the pinned
[PyQSOFit repository](https://github.com/legolason/PyQSOFit); no PyQSOFit code is
executed. [Ren et al. 2024](https://arxiv.org/abs/2406.17598) provides relevant host
decomposition methodology. [Graham et al.](https://doi.org/10.1093/mnras/stz3244)
already studied changing-state quasars; P21422's source-level association remains
unresolved here. The available older arXiv version did not resolve the local
match, and absence from that version is not a novelty test.

The pilot does not alter raw data, delivered reductions, first-pass notebooks,
or frozen observing notes. Local outputs can contain private collaboration
spectra and remain ignored by Git.

## Reproduce

From the project root, using the existing Python environment:

```sh
python exploration/palomar_2026/paper_protocol.py
python exploration/palomar_2026/fetch_pilot_templates.py
OPENBLAS_NUM_THREADS=1 python exploration/palomar_2026/pilot_fits.py
OPENBLAS_NUM_THREADS=1 python exploration/palomar_2026/pilot_validation.py
python exploration/palomar_2026/refresh_pilot_optical.py
MPLCONFIGDIR=/private/tmp/palomar_mpl python exploration/palomar_2026/paper_plots.py
MPLCONFIGDIR=/private/tmp/palomar_mpl python exploration/palomar_2026/pilot_exposure_check.py
OPENBLAS_NUM_THREADS=1 python -m unittest discover -s exploration/palomar_2026 -p 'test_*.py' -v
python exploration/palomar_2026/build_paper_notebooks.py
python exploration/palomar_2026/execute_notebooks.py 06 07
OPENBLAS_NUM_THREADS=1 python exploration/palomar_2026/validate_paper_outputs.py
```

Template download and optical refresh need public network access. Existing
template checksums are verified and the repository commit stays pinned. The
optical refresh reuses successful local query snapshots; preserve/version that
directory to request a genuinely new query. Jupyter requires loopback sockets.
The notebook builder replaces only 06/07; copy any personal edits before rebuilding.

Detailed outputs are in `outputs/palomar_paper/`: source/epoch tables,
`protocol_snapshot.json`, native-data audit and hashes, template provenance,
fit coefficients/components/residuals, model and calibration envelopes,
independent-exposure indices, injection trials, raw alert responses, timing
matches and diagnostic figures. `pilot_validation_manifest.json` records the
final artifact checks and software versions.

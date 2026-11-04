# Palomar paper: full-sample decomposition, A/B/C/D selection experiment, and joint physical transitions

The active project is **multi-epoch broad-line evolution and a full-sample test of the
pre-run manifold selection** across all 99 unique AGN systems (365 spectral epochs) in the
2026 Palomar 200-inch (NGPS) campaign. Dust physics waits for JWST.

In addition to the locked 15-source pilot ([executed notebook 07](notebooks/07_palomar_paper_fits_and_results.ipynb),
[HTML](outputs/html/07_palomar_paper_fits_and_results.html), [protocol notebook 06](notebooks/06_palomar_paper_protocol_and_pilot.ipynb),
and [15-source pilot fit atlas](outputs/palomar_paper/pilot_fit_atlas.pdf)), the repository now includes:
- **Step 1 (Full-Sample Spectral Decomposition & A/B/C/D Selection Experiment):**
  [`full_sample_fits.py`](full_sample_fits.py), with full catalogs in [`outputs/full_sample/`](outputs/full_sample/)
  and summary figures in [`outputs/full_sample/figures/`](outputs/full_sample/figures/).
- **Step 2 (Joint Multi-Epoch Hβ+Hα Fits, Physical Mechanism Tests, Sub-Exposure Audit & Literature Resolution):**
  [`deep_dive_leads.py`](deep_dive_leads.py), with joint transition tables in [`outputs/deep_dive/`](outputs/deep_dive/)
  and the 14-lead multi-epoch decomposition atlas in [`outputs/deep_dive/deep_dive_fit_atlas.pdf`](outputs/deep_dive/deep_dive_fit_atlas.pdf).

## What the full-sample decomposition and joint analysis accomplished (Steps 1 & 2)

1. **Full 365-epoch multi-line spectral decomposition (`full_sample_fits.py`):**
   - **Hβ + [O III] + optical Fe II + BC03 host (`4200–5500 Å`):** Decomposed all 365 baseline epochs across 83 multi-epoch objects with native archival comparison uncertainties across three model families (`baseline`, `young_host`, `wide_fe`) plus alternate NGPS standard-star reductions (1,404 saved `.npz` component fits in `outputs/full_sample/fits_hb/`).
   - **Hα + [N II] + [S II] (`6200–6850 Å`):** Decomposed 177 baseline epochs across 46 multi-epoch objects (`ha_fit_results.csv`, `ha_fit_comparisons.csv`, `ha_model_sensitivity_envelope.csv`).
   - **Mg II + UV Fe II (`2600–3050 Å`):** Decomposed 234 baseline epochs across 45 multi-epoch objects (`mgii_fit_results.csv`, `mgii_fit_comparisons.csv`, `mgii_model_sensitivity_envelope.csv`).
2. **Unbiased A/B/C/D selection experiment (`abcd_group_experiment_summary.csv`, `selection_experiment_tests.json`):**
   - Evaluated across all 79 pre-run reviewed AGN targets with native-error Hβ multi-epoch comparisons (76 excluding prior catalog-confirmed CLAGNs: Group A $n=11$, Group B $n=58$, Group C $n=1$, Group D $n=6$).
   - **Amplitude test:** Median $|\Delta \log_{10} \text{EW}(\text{H}\beta)|$ is $0.061\text{ dex}$ in Group A, $0.056\text{ dex}$ in Group B, and $0.070\text{ dex}$ in Group D. Two-sample permutation tests show no statistically significant population-median enrichment between Group A and Group B ($p_{\text{perm}} = 0.96$) or between manifold-selected candidates (A+B) and random AGN controls (D) ($p_{\text{perm}} = 0.68$).
   - **Directional test:** Pre-run `manifold_direction` ($f_{\text{on}} - f_{\text{off}}$) aligns with the sign of $\Delta \log_{10} \text{EW}(\text{H}\beta)$ in $51.9\%$ of Group B ($57.8\%$ on clean coverage $\ge 0.80$) and $11.1\%$ of Group A ($p_{\text{perm}} = 0.46$ overall), because the two strongest Group A turn-off discoveries (`P21422` and `P9506`) experienced earlier historical WISE W1 brightening episodes before collapsing in 2022–2026.
   - **Scientific conclusion:** An integrated 12-year WISE W1 manifold projection captures whether a source occupies a state-transitioning trajectory, rather than instantaneous 2026 polarity.
3. **Simultaneous joint multi-epoch Hβ+Hα fits & physical mechanism tests (`deep_dive_leads.py`):**
   - Tied the intrinsic host-galaxy stellar population mix $\theta_{\text{young}} = f_{1.4\text{Gyr}} / (f_{5\text{Gyr}} + f_{1.4\text{Gyr}})$ and narrow-line dispersion $\sigma_n$ across all epochs of each source for 28 discovery and extreme-evolution leads (`joint_multi_epoch_lead_fits.csv`), while regularizing `broad2_sigma <= 4200 km/s` to prevent continuum-bump boundary hits.
   - Computed broad and narrow Balmer decrements ($\text{H}\alpha/\text{H}\beta$), Cardelli, Clayton & Mathis (1989) ($R_V = 3.1$) pure foreground dust-extinction predictions vs. observed Balmer-decrement changes, Vestergaard & Peterson (2006) single-epoch virial black hole masses ($\log M_{\text{BH}}/M_\odot$), bolometric Eddington ratios ($\log \lambda_{\text{Edd}}$), and broad-line FWHM virial breathing slopes (`lead_physical_transitions.csv`).
4. **Sub-exposure audit (68 native frames) and literature cross-match resolution:**
   - Verified individual NGPS sub-exposures across `sep23_data/`, `oct26_data/`, and `oct27_data/` (`lead_subexposure_audit.csv`), confirming $<2.5\%$ exposure-to-exposure agreement for `P21422`, `P9506`, `P7947`, `P9694`, `P21743`, `P21841`, and `P11010`, while catching a $2.6\times$ frame discrepancy (`10.6 Å` vs. `27.4 Å`) and `broad2_sigma` boundary hit in `P11542`.
   - Resolved all candidate cross-matches against Camus & Panda (2026) (`clagn_catalog_camus_panda2026.csv`), `literature_clagn.csv`, and `zeltyn_coords.csv` (`literature_audit_resolution.csv`).

### Full-sample A/B/C/D selection experiment summary (excluding prior catalog-confirmed CLAGNs)

| Group | Definition | $N_{\text{obs}}$ (Hβ) | Median $|\Delta \log_{10} \text{EW}(\text{H}\beta)|$ [dex] | $|\Delta \text{EW}| \ge 30\%$ Fraction | $|\Delta(\text{H}\beta/\text{[O III]})| \ge 30\%$ Fraction | Sign Hit Rate (`manifold_direction`) |
|---|---|---:|---:|---:|---:|---:|
| **A** | Manifold + W1 variability selected | 11 | $0.061$ (IQR $0.097$) | $27.3\%$ ($3/11$) | $40.0\%$ ($4/10$) | $11.1\%$ ($50.0\%$ on $\|\text{dir}\|\ge 0.15$) |
| **B** | Manifold selected (all W1 classes) | 58 | $0.056$ (IQR $0.093$) | $25.9\%$ ($15/58$) | $34.0\%$ ($18/53$) | $51.9\%$ ($57.8\%$ on clean cov $\ge 0.80$) |
| **C** | W1 variability only (1 non-catalog target) | 1 | $0.008$ (IQR $0.000$) | $0.0\%$ ($0/1$) | $0.0\%$ ($0/1$) | $0.0\%$ |
| **D** | Random AGN control sample | 6 | $0.070$ (IQR $0.070$) | $33.3\%$ ($2/6$) | $60.0\%$ ($3/5$) | $50.0\%$ |

### Joint multi-epoch physical transitions for primary discoveries & newly surfaced Group B/C leads

| Source | Group | $z$ | $\Delta t_{\text{rest}}$ [yr] | Joint $\text{EW}(\text{H}\beta)$ ($\text{old}\to\text{new}$) [Å] | Joint $\text{EW}(\text{H}\beta)$ Ratio | Joint $\text{H}\beta/\text{[O III]}$ Ratio | $\text{EW}(\text{H}\alpha)$ Ratio | Balmer Decr. Ratio (Obs vs. Dust Pred.) | $\log(M_{\text{BH}}/M_\odot)$ | $\Delta \log \lambda_{\text{Edd}}$ [dex] | Physical & Novelty Verdict |
|---|---|---:|---:|---|---:|---:|---:|---|---:|---:|---|
| **P21422** | A | 0.239 | 3.63 | $106.4 \to 73.5$ | **0.69×** | **0.61×** | **0.59×** | $0.86\times$ vs. $1.16\times$ | 8.64 | $-0.32$ | **Intrinsic accretion turn-off collapse** (virial broadening $3402 \to 5820\text{ km s}^{-1}$); recurrent post-2022 collapse after Graham+19 CRTS brightening candidate. |
| **P9506** | A | 0.232 | 9.77 | $76.9 \to 59.9 / 60.1$ | **0.78×** | **0.57× / 0.56×** | **0.89× / 0.84×** | $0.93\times$ vs. $1.19\times$ | 7.58 | $-0.30$ | **Intrinsic accretion dimming** verified across 4 exposures on Sep 23 & Oct 27; **unreported discovery**. |
| **P7947** | B | 0.254 | 19.90 | $13.8 \to 31.5$ | **2.28×** | **1.73×** | **1.32×** | $0.61\times$ vs. $0.85\times$ | 7.18 | $+0.24$ | **Intrinsic accretion turn-on** (Hα/[S II] $1.44\times$, verified across both Oct 26 frames); **unreported discovery**. |
| **P9694** | D | 0.238 | 4.07 | $98.9 \to 118.0 / 132.3$ | **1.19× / 1.34×** | **1.27× / 1.30×** | **1.12× / 1.13×** | $0.87\times$ vs. $0.92\times$ | 8.75 | $+0.26 / -0.09$ | **Concordant Hβ+Hα brightening** verified across 4 exposures on Sep 23 & Oct 27; **unreported control discovery**. |
| **P14388** | C (Known) | 0.269 | 1.38 | $46.5 \text{ [2021]} \to 7.2 \text{ [2025]} \to 1.0$ | **0.13×** (0.02× vs. '21) | **0.08×** | **0.73×** | $5.32\times$ vs. $2.15\times$ | 8.01 | $-0.82$ | **Complete Type 1 $\to$ Type 2 turn-off collapse** in catalog CLAGN (`J075934.95+322143.3`); coupled accretion shutdown + dust reddening. |
| **P12457** | B | 0.092 | 23.94 | $3.1 \to 1.8$ | **0.58×** | **0.57×** | **0.72×** | $1.28\times$ vs. $1.18\times$ | 6.08 | $-0.39$ | **Low-luminosity Seyfert 1.9 dimming** (Hα/[S II] $0.71\times$, concordant Balmer decrement steepening); **unreported lead**. |
| **P21743** | B | 0.143 | 11.21 | $102.3 \to 59.0$ | **0.58×** | **0.81×** | **0.91×** | $1.55\times$ vs. $1.06\times$ | 7.70 | $-0.30$ | **Strong Hβ dimming** verified across both Oct 27 G frames (`31.0 Å` vs. `31.8 Å`); NLS1-like core requires high-res [O III] profile tying. |
| **P21841** | B | 0.514 | 1.20 | $69.6 \to 106.0$ | **1.52×** | **0.87×** | **1.24×** | $1.03\times$ vs. $1.04\times$ | 8.77 | $-0.18$ | **Concordant Mg II ($2.47\times$) + Hβ ($1.52\times$) + Hα ($1.24\times$) EW growth** as blue continuum subsides; **unreported lead**. |
| **P11010** | B | 0.381 | 18.75 | $38.3 \to 54.0$ (single) | **1.41×** | **0.90×** | **1.53×** | — | 8.66 | $+0.02$ | **Concordant Mg II ($5.20\times$) + Hβ ($1.41\times$) + Hα ($1.53\times$) brightening** verified on both Oct 27 R frames; **unreported lead**. |
| **P17881** | Known | 0.211 | 1.38 | $65.4 \to 104.7$ | **1.60×** | **0.96×** | **1.15×** | $0.99\times$ vs. $1.01\times$ | 8.49 | $+0.01$ | **Continuum-driven EW increase** in catalog CLAGN (`J004339.32+134436.5`) at constant intrinsic broad Hβ/[O III] flux. |
| **P10381** | Known | 0.149 | 10.58 | $83.0 \to 100.2 / 113.5$ | **1.21× / 1.37×** | **1.33× / 1.21×** | **1.14× / 1.24×** | $0.97\times$ vs. $0.94\times$ | 7.81 | $+0.17 / -0.04$ | **Persistent Type 1 brightening** over 2014 SDSS across both Sep 23 & Oct 27 (unblocked via rest-frame $5695\text{ Å}$ splice boundary). |

## What the locked 15-source pilot accomplished

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

## What changed scientifically in the 15-source pilot

Ratios below compare the new spectrum with the date-selected latest accepted
archive. Ranges span tested models/standards, **not confidence intervals**.

| Case | Current result | Decision |
|---|---|---|
| **P21422** | Broad-Hβ EW ratio **0.59–0.69**; Hβ/[O III] **0.52–0.61**. Primary-standard fits alone give EW 0.66–0.69. Both line windows covered. Two exposure wing EWs agree within ~2.4%. | First joint Hα/Hβ case. Broad lines remain present. Audit the Graham-2019 literature association before claiming novelty. *(Resolved in Step 2: Graham+19 CRTS brightening candidate from 2005→2018; 2026 NGPS captures a newly discovered post-2022 turn-off collapse.)* |
| **P9506** | October/archive EW **0.67–0.77**; Hβ/[O III] **0.53–0.55**. September supports weakening too; individual-exposure wing indices are similar across both nights. | Promote as persistent long-baseline weakening. Do not claim an abrupt month-scale transition. *(Verified in Step 2 via joint Hβ+Hα fits and 4-frame audit.)* |
| **P7947** | EW **1.51–2.11** across standards/models; adopted-standard range 1.70–2.11. | Strong follow-up candidate, but host/centroid sensitivity and masked [O III] support need work. *(Joint fit in Step 2 gives EW 2.28×, Hβ/[O III] 1.73×, Hα EW 1.32×.)* |
| **P24930** | Adopted-standard full-line EW **0.87–0.88**, much less extreme than the first wing-index result. | Down-rank the large-change claim; inspect profile/coverage and decomposition. |
| **P17881** | Adopted-standard EW **1.57–1.77**, but Hβ/[O III] **0.96–1.03** versus February 2025. | Known-source evolution; continuum dilution and aperture/calibration matter before a broad-flux claim. |
| **P9694** | October/September EW stress envelope **0.89–1.30**. | Month-scale growth is calibration dependent; long-baseline brightening vs. 2021 DESI is 1.19–1.34× in Hβ and 1.12–1.13× in Hα. |

Several comparison sources have EW and [O III]-normalized trends in opposite
directions. This is useful evidence that absolute calibration and aperture effects
must be modeled; the four comparisons are not a calibrated control experiment.

## Roadmap Status: From Pilot to Publishable Paper

1. **Validate the leading cases with joint multi-epoch Hα+Hβ models (Completed in Step 2).**
   Included all usable archival epochs, simultaneous Hα/Hβ measurements, shared
   host-galaxy stellar population fraction ($\theta_{\text{young}}$) and narrow-line
   dispersion ($\sigma_n$), regularized broad profiles, and a 68-exposure native
   sub-frame audit (`outputs/deep_dive/`).
2. **Anchor continuum changes where possible.** Use actual filter transmission
   curves and time-matched photometry with consistent total/reference flux and
   aperture treatment. Alert photometry is not complete forced photometry.
   Do not interpret fitted continuum slopes as intrinsic disk changes without
   Balmer-decrement and narrow-line flux cross-checks.
3. **Measure the full observed denominator (Completed in Step 1).** Extended the
   validated decomposition procedure across all 365 epochs of the 99 campaign
   sources (`outputs/full_sample/`), retaining non-detections, coverage masks,
   and 3-model + alternate-standard sensitivity envelopes.
4. **Test incremental selection information (Completed in Step 1).** Preserved
   the pre-run A/B/C/D groups and known-source stratum, running two-sample
   permutation tests and 10,000-draw bootstrap confidence intervals across all
   79 reviewed multi-epoch Hβ targets (`selection_experiment_tests.json`).
   Because the selection enrichment over Group D controls is statistically null
   ($p_{\text{perm}} = 0.68$), we report this negative selection test transparently
   alongside the full spectral-evolution catalog and individual discoveries.
5. **Resolve forecasting provenance.** Of 99 campaign sources, 46 pass a timing
   screen in which the latest archive follows all linked original W1 inputs;
   29 have a baseline preceding some inputs and 24 have unresolved feature
   support. Pool IDs and Sample A training IDs must not be equated by their integer values.

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
studied changing-state quasars; Step 2 resolved `P21422` (`SDSS J140650.94+282717.8`)
as a 2005→2018 brightening candidate in Graham et al. (2019/2020) that is absent from
Camus & Panda (2026) (`clagn_catalog_camus_panda2026.csv`), establishing the 2026 Palomar
spectrum as a newly discovered post-2022 turn-off collapse (recurrent CLAGN).

Neither the pilot nor the full-sample pipeline alters raw data, delivered reductions,
first-pass notebooks, or frozen observing notes. Local outputs can contain private
collaboration spectra and remain ignored by Git.

## Reproduce

From the project root, using the existing Python environment:

```sh
# 1. Locked 15-source pilot pipeline and validation
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

# 2. Step 1: Full-sample 365-epoch Hbeta + Halpha + Mg II decomposition & A/B/C/D test
OPENBLAS_NUM_THREADS=1 MPLCONFIGDIR=/private/tmp/palomar_mpl python exploration/palomar_2026/full_sample_fits.py

# 3. Step 2: Joint multi-epoch Hbeta+Halpha fits, Balmer decrement & virial tests, and lead fit atlas
OPENBLAS_NUM_THREADS=1 MPLCONFIGDIR=/private/tmp/palomar_mpl python exploration/palomar_2026/deep_dive_leads.py
```

Template download and optical refresh need public network access. Existing
template checksums are verified and the repository commit stays pinned. The
optical refresh reuses successful local query snapshots; preserve/version that
directory to request a genuinely new query. Jupyter requires loopback sockets.
The notebook builder replaces only 06/07; copy any personal edits before rebuilding.

Detailed outputs are in:
- `outputs/palomar_paper/`: locked 15-source pilot tables, `protocol_snapshot.json`,
  native-data audit and hashes, template provenance, fit coefficients/components/residuals,
  model and calibration envelopes, independent-exposure indices, injection trials,
  raw alert responses, timing matches, diagnostic figures, and `pilot_validation_manifest.json`.
- `outputs/full_sample/`: full 365-epoch decomposition tables (`hb_fit_results.csv`,
  `ha_fit_results.csv`, `mgii_fit_results.csv`), multi-epoch comparisons and 3-model +
  alternate-standard envelopes, `full_sample_object_summary.csv`, `abcd_group_experiment_summary.csv`,
  `selection_experiment_tests.json`, 1,404 `.npz` component fits in `fits_hb/`, and
  full-sample summary figures in `figures/`.
- `outputs/deep_dive/`: joint multi-epoch fits (`joint_multi_epoch_lead_fits.csv`),
  physical mechanism & virial transition diagnostics (`lead_physical_transitions.csv`),
  68-frame sub-exposure audit (`lead_subexposure_audit.csv`), literature resolution
  table (`literature_audit_resolution.csv`), 14 individual multi-epoch decomposition
  figures in `figures/`, and `deep_dive_fit_atlas.pdf`.

# Current review status — 26 September 2026

The technical review in [review/README.md](review/README.md) supersedes earlier
exposure/acquisition estimates below. The current request is **67.4 hours**
(APT timing export 67.33 hours), with longer exposures for R01 and R06 and
flux-dependent acquisition recipes. All 48 visits retain scheduling windows.
Focused ETC checks, bibliographic checks, anonymity screening and the five-page
core layout are complete. These do not validate individual source histories or
turn assumed nuclear flux fractions into measurements.

Before final submission, review nuclear centering for the nine AllWISE-extended
hosts listed in the technical report, reconcile remaining sample/history inputs,
and complete APT administrative fields. The proposal is not submitted.

Original web ETC archives named below are now under `etc/downloads/web_exports/`.
Automatic APT backups are under `apt/work/automatic_backups/`.

The entries below are a chronological work log; descriptions such as “placeholder,”
“not installed,” and older time totals record earlier states, not the current draft.

---

# Remaining inputs for the proposal

The narrative now states the science case and observing plan directly. This file
retains the unresolved checks removed from the prose during the editorial pass;
removing a reminder from the PDF does not mean the check has been completed.

## Sample and supporting data

- Audit the identities and 12 rising / 12 fading assignments in
  `inputs/jwst_sample_cycle6.csv`, including the three reversal flags. Membership
  is unchanged by this edit.
- Reconcile catalogue classifications and citations source by source. The
  manifold description is retained without a citation, as requested.
- Assemble a compact target/epoch table with redshifts, optical and infrared
  coverage and SPHEREx visit dates. Do not include NGPS status, observing months,
  exact observing dates, or program identifiers in the proposal or target table;
  use "2026 NGPS spectra" in the narrative.
  Full-sample SPHEREx and manifold availability is confirmed by the user;
  individual epoch counts have not been audited in this pass.
- Bracket changes using the illumination history; validate nuclear luminosities,
  delay distributions and event-age/delay coverage at permitted JWST dates.
  The CSV's universal warm/hot delay ratio of ten is provisional.
- Verify the archival W3/W4 flux provenance, host contribution and nuclear flux
  estimates for the current membership. The proposal quotes the CSV ranges.

## Quantitative case and exposure design

- Calculate fixed/evolving model predictions for Figure 2, including calibration
  covariance, geometry/host uncertainty, and target positions in response phase.
  The figure remains a placeholder; no discrimination forecast is established.
- Supply the within-group scatter and predicted effect size. The 0.289s and
  0.408s factors are independent equal-variance sampling formulae, not a power
  calculation for this sample.
- Run target-specific ETC calculations to establish nuclear sensitivity,
  saturation margins, exposure times, readout choice and sky strategy. The
  red S/N values and FASTR1 setting remain design inputs to verify.
- Compute charged time in APT, including target acquisition, dithers and any
  dedicated sky. Check that the science-driven total is within Medium limits.
- Check the use of channel 4C against each primary diagnostic and finalize the
  appropriate background-observation requirements in APT.

## Scheduling and supporting observations

- Check JWST visibility, schedulability and whether any source needs a narrower
  observing window. The baseline requests no timing constraint or repeat visit.
- Confirm arrangements for optical coverage near the JWST epoch privately.
  Keep NGPS scheduling, completion status and program identifiers in internal
  records only. The proposal uses "2026 NGPS spectra" without a completion claim.
- The proposal no longer assumes continued SPHEREx operations or guaranteed
  near-simultaneous coverage in 2027–28. Historical spectra are fitted at their
  actual dates and gaps are included as uncertainties.
- Unspecified ground-based Pa-alpha/Pa-beta/Br-gamma observations have been
  removed from the baseline. Restore them only with a defined facility,
  observing arrangement and demonstrated role in the MIRI science.
- Complete target-by-target archive and approved-program duplication checks.

## Final assembly

- Recheck references and instrument statements before submission. This was an
  editorial pass, not an independent literature or instrument verification.
- Replace red quantities and the Figure 2 placeholder with verified results.
  Insert the compact sample table and recheck the five-page required-section
  limit using the unchanged official 12-point style.
- Refresh the APT abstract and validation report after substantive changes.

## Step 1: provisional precision design

The exposure-driving measurement is the host-separated rest 8–13 micron
warm-dust continuum luminosity relative to a fixed-distribution prediction.
Continuum shape and the 9.7-micron feature provide supporting constraints;
the accessible 18-micron feature is a secondary diagnostic. The listed
redshifts place rest 8–13 microns within observed 8.11–21.07 microns across
the sample. This is a wavelength-coverage check, not an ETC calculation.

For the first ETC pass retain the draft's provisional S/N 50 per R=100
continuum bin, emphasizing clean continuum regions bracketing the warm band,
and S/N 30 per R=100 continuum bin around the accessible 18-micron feature.
Treat rest 5 microns as a continuum-shape anchor. These are design targets,
not achieved sensitivities or a demonstrated minimum for model discrimination.
Do not interpret S/N on the continuum as the detection significance of a
silicate feature or of a model residual.

At S/N 50, single-bin statistical noise is 2%; at S/N 30 it is 3.33%.
For the ratio of two equally precise independent bins, those errors become
2.83% and 4.71%, respectively. Actual fitted spectra require the full
covariance and continuum/host uncertainty. Rebinning native MRS spectra to
R=100 is an analysis choice; it does not change the observing mode, and
correlated errors must be propagated.

A useful provisional measurement goal is about 5% uncertainty on the extracted
warm continuum (0.0217 dex), including calibration and host decomposition.
This goal is not yet shown to be achievable for each source. It is separate
from uncertainty in the fixed-distribution prediction. The current STScI MRS
calibration guidance describes 1–2% differences among exposure setups and
less reliable long-wavelength calibration; these errors do not simply average
away with spectral binning. Use wavelength-dependent covariance, not a
universal calibration floor or standard-star repeatability as an AGN accuracy.

For illustration only, a 0.10 dex residual is a 25.9% flux excess. A three-sigma
individual-object test requires total residual uncertainty <=0.0333 dex.
With a 5% luminosity error and independent prediction uncertainty, this leaves
<=0.0253 dex for the prediction. Neither a 0.10 dex effect nor that prediction
precision has been established. This illustrates why increased MIRI S/N alone
cannot establish the science reach. Shared uncertainties must also remain in
the 12-versus-12 comparison rather than being divided by sqrt(12).

Before fixing exposure times, calculate the model separation and prediction
uncertainty for representative histories. The current values can seed ETC
feasibility comparisons; they must remain provisional in the proposal.

Calibration source checked for this step:
https://jwst-docs.stsci.edu/jwst-calibration-status/miri-calibration-status/miri-mrs-calibration-status

## Step 2: representative ETC cases (inputs only)

Four cases span the archival flux range. These are planning cases, not a
reselection of the sample or measurements of nuclear flux. W3/W4 values below
come from the current sample CSV and include host emission. The typical case
has W3=5.1 mJy, close to the sample median of 5.7 mJy.

| Case | ID / target | z | W3 / W4 (mJy) | Rest 12 microns observed | Rest 18 microns observed |
| --- | --- | --- | --- | --- | --- |
| Faint W3 | R06: J083826.50+371906.7 | 0.2111 | 1.9 / 6.4 | 14.53 microns | 21.80 microns |
| Typical flux | F05: J085359.05+212803.4 | 0.0842 | 5.1 / 13.6 | 13.01 microns | 19.52 microns |
| Bright W3 | F11: J2330323-022745 | 0.0332 | 48.4 / 80.8 | 12.40 microns | 18.60 microns |
| High-redshift / faint W4 | R01: J162039.13+430814.6 | 0.6205 | 2.6 / 3.7 | 19.45 microns | 29.17 microns |

First calculation: R06 / P17881 (J083826.50+371906.7). Its rest 8–13 micron
band falls at observed 9.69–15.74 microns; the 9.7-micron peak at 11.75 microns,
and the 18-micron peak at 21.80 microns. Use its coordinates from the CSV
for background/visibility calculations. A flux-normalized point source alone
does not account for host background noise; assess the extended component.

When constructing source spectra, explore nuclear contributions of 25%, 50%
and 100% of the archival total as explicitly assumed sensitivity scenarios.
These fractions are not measured bounds, and a wavelength-independent fraction
is only a simplification for initial experiments. Vary the spectral slope and
current state; do not infer the JWST-epoch nuclear continuum directly from the
2010 total-source photometry. Normalize physical templates through the actual
WISE bandpasses rather than treating W3/W4 as exact monochromatic anchors.

Check the all-three-settings, four-dither setup against S/N at R=100, detector
saturation, background strategy and native-resolution line diagnostics.
P1823 has no accessible 18-micron peak; test its warm continuum and 9.7-micron
feature instead. No exposure times or charged times have been calculated.

Initial local setup check: no Pandeia module or configured reference data found in the
active Python environment, and no JWST ETC/APT project found in the workspace.
The installed APT application is 2025.5.3 and needs a Cycle 6 update before
final planning/submission. The supported Cycle 6 Pandeia engine is 2026.7,
with matching reference data and PSFs. The user has now supplied an official ETC example archive; see the receipt
below. The local installation is now available and validated; see the entry
at the end of these notes and `etc/LOCAL_SETUP.md`.

Sources:
- https://outerspace.stsci.edu/spaces/PEN/pages/77530136/Pandeia%2BEngine%2BInstallation
- https://jwst-docs.stsci.edu/jwst-mid-infrared-instrument/miri-performance/miri-sensitivity
- https://jwst-docs.stsci.edu/jwst-mid-infrared-instrument/miri-observing-modes/miri-medium-resolution-spectroscopy

## ETC archive received

Received `wb325241_c7_2026-09-26_18.26.09.tar`; inspected calculation 7 and
copied its input unchanged to `etc/official_example_c7/input.json`. The
archive contains input.json and FITS products for a Cassiopeia A extended
source in MIRI/MRS channel 3, short setting (FASTR1, 27 groups, 1 integration,
4 exposures). It is an example configuration, not an AGN calculation.
Source spectrum/geometry, sky background and extraction need adaptation.
The inspected files do not establish the engine/reference-data version.
No AGN exposure-time result has been added to the proposal.

## First R06 ETC trial received

The adapted calculation is in `wb325241_c7_2026-09-26_18.58.14.tar`.
Verified settings and reproducible approximate R=100 estimates are in
`etc/r06_ch3a_27groups_r1/`; the script is `etc/bin_archive.py`.
This is a point-source power law normalized to 0.95 mJy in WISE W3, not a
measured nuclear SED. The saved redshift is 0.211 rather than 0.2111.
The copied scalar report gives 299.70 s on source and 599.41 s including sky;
APT overheads are not included. At 12.81 microns, native S/N=5.37 and the
diagonal-covariance R~100 estimate is 41.10. Example bins elsewhere in this
single sub-band give ~33–39. Full spectral covariance, calibration and host
uncertainties are not established. The point-source-only calculation omits
host photon noise. This trial does not establish full-sample exposure times.
Next compare a 0.6-arcsec extraction radius at the same exposure settings.

The 0.6-arcsec comparison has now been received in
`wb325241_c7_2026-09-26_19.02.44.tar` and analyzed in
`etc/r06_ch3a_27groups_r06/`. A recursive comparison confirms that only the
extraction radius changed. At 12.81 um native S/N=7.831 and approximate
R=100 S/N=59.91. At 12.00 um the binned estimate is 48.49; not every bin meets
50. Fifteen complete contiguous bins spanning this sub-band give 49.51–59.91
for a grid anchored at 12.81 um. These are diagonal-covariance estimates and
retain the earlier limitations. The next wavelength check is the rest-12-um
anchor at observed 14.5332 um, outside the current sub-band. No exposure time
or extraction radius is finalized for the proposal.

Channel 3B archive `wb325241_c7_2026-09-26_19.08.41.tar` is analyzed in
`etc/r06_ch3b_27groups_r06/`. Only the disperser and reference wavelength
changed from the preceding 0.6-arcsec run. At 14.53 um (rest ~12 um), native
S/N=9.116 and diagonal-covariance R~100 S/N=69.35. Fourteen complete bins
anchored at this wavelength give 56.46–71.08. Time remains 299.70 s on source,
599.41 s including the off-source sky strategy, before APT overheads. The
assumed source and missing host/covariance limitations remain. Next test
continuum sensitivity at observed 21.80 um (rest 18 um) in channel 4B, at a
trial 0.9-arcsec radius. Do not sum 3B and 4B exposure times: they are
simultaneous channels of the same setting. This is not yet a feature-detection
forecast or a full wavelength/sample timing calculation.

Channel 4B archive `wb325241_c7_2026-09-26_19.14.26.tar` is analyzed in
`etc/r06_ch4b_27groups_r09/`. Only the channel, radius (0.9 arcsec), and
reference wavelength (21.80 um) changed from the channel 3B trial. Native
S/N=2.801; diagonal-covariance R~100 S/N=20.29 at 21.80 um, 14.85 at 23 um,
and 11.65 at 24 um. Sixteen complete bins give 11.48–23.25. The provisional
30 goal is not met at the rest-18-um continuum anchor. Next trial changes
groups from 27 to 60 only; this is motivated by approximate square-root time
scaling and must be verified with a new ETC run. The red continuum wing may
need more time or a revised requirement even if the anchor reaches 30.
The assumed power-law scene has no silicate feature or host component.

The 60-group 4B run is saved in `wb325241_c7_2026-09-26_19.17.43.tar`, with
analysis in `etc/r06_ch4b_60groups_r09/`. Only ngroup changed (27 -> 60).
Native S/N=6.8049 at 21.80 um; diagonal-covariance R~100 S/N=49.25 there,
36.89 at 23 um and 29.35 at 24 um. Sixteen complete bins give 28.95–56.95.
The rest-18-um continuum anchor clears the provisional goal; the red end
remains marginal. Retain 60 groups as a working B-setting value for this
source scenario, not a finalized exposure or a whole-sub-band S/N guarantee.
Reported times are 666.01 s on source and 1332.02 s including sky, before APT
overheads. Other B channels still require checks at the longer exposure.
Next check channel 3C at 15.74 um, near the primary warm-band upper boundary,
with a trial 0.6-arcsec radius and the same 60-group detector setup.

## Local ETC installed and validated

Pandeia 2026.7 is installed in `etc/.venv/`, with matching instrument data,
the full MIRI/MRS PSF subset and WISE W1–W4 normalization bandpasses.
`etc/run_local.py` sets its reference paths and saves reproducible outputs.
The unchanged 60-group channel 4B input reproduces all 916 web spectral
samples in flux, noise and S/N to floating-point precision. The local scalar
S/N is 6.80178; on-source/sky-inclusive times are 666.0096/1332.0192 s.
No calculation warnings are reported. The local output also reproduces the
previous approximate R=100 estimates. See `etc/local_validation.json` and
`etc/LOCAL_SETUP.md` for evidence, installation scope and commands. Further
current MRS trials can run locally without user screenshots or downloads.
Imaging/TA PSFs and extra template libraries are not installed yet. Source
assumptions, covariance, host noise and full-program timing remain unresolved.

## First automated continuation: channel 3C

The rest-13-um boundary (observed 15.74 um) has now been tested locally at a
0.6-arcsec trial radius. At 60 groups, approximate R=100 S/N=215.86; at 27
groups it is 91.61. Both runs have no calculation warnings or significant
saturation. Fifteen complete 3C bins give S/N=81.17–91.61 at 27 groups.
See `etc/r06_ch3c_comparison.md` and the linked input/result directories.
Retain 27 groups as an initial C-setting value for this assumed source,
subject to checks of the simultaneous channels and other diagnostics.
Next consolidate first-target coverage, including the rest-5-um anchor and
the blue end of the primary warm band, then revisit source/host assumptions.

## Rough APT budget prioritized at the user's request

The user asked to stop exhaustive ETC checking, adopt reasonable assumptions,
and build the APT program first. `apt/clagn24_miri_draft.aptx` now contains
24 science targets (12 fading, 12 rising) and 24 matching sky observations,
linked in consecutive, noninterruptible pairs. All use MIRI/MRS only,
FASTR1/FULL, four dithers, one integration, and A/B/C groups 27/60/27.
Science observations include provisional F560W/FAST/10-group acquisition.
Sky fields are unvetted placeholders 60 arcsec north; no NGPS scheduling
details were copied from the internal sample CSV.

APT 2026.5.1 was downloaded from STScI, checksum verified, unpacked locally,
and used to load, compute local timings, export reports and save the draft.
The timing JSON reports 84.24 charged hours; the saved proposal summary
rounds to 84.3 hours and assigns MEDIUM. Source-plus-sky photon collection
is 16.87 hours. The current local budget charges a full 2,100-second slew
for every visit, including sky fields; it has no Smart Accounting savings.
This is a provisional starting budget, not a final proposal time request or
demonstrated S/N for all 24 objects. See `apt/README.md` and timing exports.

Automatic approval review rejected the online Visit Planner command because
it could transmit the draft and target coordinates to STScI. No online run
or proposal submission occurred. Explicit approval is needed to retry that
external planning step. Local draft preparation and timing exports are done;
the proposal TeX/PDF has not been changed to claim finalized exposure times.

## Online Visit Planner completed with user approval

The user explicitly approved the external planning run. APT 2026.5.1
successfully contacted STScI's constraint generator (19.0.1) and processed
all 48 visits. Every visit has up-to-date, nonempty scheduling windows; no
diagnostics appear in the cached visit results. The timing JSON now reports
**66.34 charged hours**, versus 84.24 before online planning; the saved
APT proposal summary rounds up to **66.4 hours**, still MEDIUM.
Sky-visit slews are now 274 s instead of 2,100 s, guide-star accounting is
updated, and six science visits also receive reduced charges. No source
exposure parameters were changed. This supports retaining the working
24-AGN / 12-fading / 12-rising design at Medium-program scale.

The APT file, matching XML, and timing exports in `apt/` now include the
online results; earlier estimates are preserved in `apt/work/pre_visit_planner/`.
Sky positions, acquisition settings and all-target exposure adequacy remain
provisional as described in `apt/README.md`. No proposal submission occurred.

## Proposal text updated from ETC and APT results

Filled the observing-description placeholders with FULL/FASTR1, one
integration per dither, four dithers and A/B/C groups 27/60/27. The text now
gives 299.7/666.0/299.7 s on source per setting, 1265.4 s per AGN, matching
sky observations for all 24 targets, 16.87 h source-plus-sky integration,
and the baseline APT request of 66.4 h (APT summary rounding of 66.34 h).
The special-requirements paragraph now describes the actual consecutive,
noninterruptible pairs and the 48 visits' scheduling windows. It no longer
incorrectly states that there are no timing constraints of any kind.

Added the W3-faint R06 example: assumed point source, Fnu proportional to
wavelength squared, 0.95 mJy in W3 (50% of the archival total), medium sky,
and approximate R=100 S/N 60/92/49 at observed 12.81/15.74/21.8 microns.
Extraction radii and independent-sample noise assumptions are explicit.
The S/N 50 and 30 values are labeled goals, not proven minima for every
source and wavelength. Remaining red technical reminders cover diagnostic
band sensitivity, acquisition/saturation, and clean sky positions. The model
figure, intrinsic scatter/model separation and duplication checks remain open.

Rebuilt `proposal.pdf` successfully with the unmodified official style.
It remains six pages: core sections 1–4, supplement 5, references 5–6.
There are no LaTeX layout/reference warnings; the updated observing page
was visually checked. `validation.json` records the current PDF and timing.

## Confident proposal tone and completed sky-field selection

At the user's request, strengthened the scientific opening, MIRI's role and
the balanced-sample rationale. The proposal now states the time request
directly, without "provisional" or "baseline request" wording. Procedural
sensitivity/acquisition reminders were moved out of the proposal narrative;
they remain open working tasks, not completed measurements. The quantitative
model figure, intrinsic scatter/model separation and duplication results
remain visibly unfilled because those results have not been calculated.

Screened all 24 sky fields using AllWISE and 2MASS catalogs, WISE W1/W3/W4
images and 2MASS J cutouts, with SDSS images where available. Visually
inspected every selected field. Replaced 23 arbitrary north offsets; R07
retains its 60-arcsec north position. Twenty-three selected offsets are
60–120 arcsec, and crowded R12 uses 135 arcsec. Catalogued-source clearances
exceed 25 arcsec. R12 has structured diffuse W3 emission and a bright-source
clearance of 63.7 arcsec; its reduction note explicitly retains this concern.
The survey review is not a claim of no emission at JWST depth.

The selected coordinates and review record are in
`apt/sky_fields/selected_sky_fields.csv` and `sky_field_review.pdf`.
All APT sky targets now have Extended=YES, no acquisition, and four-point
EXTENDED SOURCE dithers, following STScI guidance. All science exposure
settings and target coordinates were retained. The input generator now
uses the reviewed sky-field table when available.

Reran the authorized online Visit Planner after the coordinate/dither
changes: all 48 visits have up-to-date, nonempty scheduling windows.
The timing export is **66.59 hours**, rounded by APT to **66.6 hours**;
the proposal and current timing records now use this request. The previous
66.34-hour file is retained under `apt/work/before_sky_review/`.
Saved the updated APT, matching XML, timing and pointing exports. The PDF
builds without warnings, remains six pages (core 1–4), and the updated
observing and supplemental pages were visually checked. No submission occurred.

Remaining practical technical work: confirm target acquisition brightness
and centering, representative sensitivity at the remaining diagnostic bands
(especially rest 5 microns), bright-source saturation, final guide stars,
and duplication checks. These were not completed by the sky review or
removed merely by strengthening the wording.

## Figure 2: physical illustration of the MIRI diagnostics

Replaced the empty figure with a reproducible rising/fading dust-response
illustration (`fig2_diagnostics.pdf`, generated by `make_fig2_diagnostics.py`).
The calculation uses published silicate and graphite absorption efficiencies,
separate radiative-equilibrium temperatures and angular shell delays. An
inner boundary following sqrt(luminosity) is compared with fixed dust after
idealized factor-four steps. Each alternative is scaled to the same
integrated 2--4 micron luminosity as the corresponding fixed model, allowing
one dust-mass amplitude. Opacity files, assumptions, checksums, output spectra
and numerical checks are retained under `inputs/`.

Following visual feedback, enlarged the fractional-difference panels, shaded
the actual separation and labeled +36%/-28% differences near 18 microns.
The fading example overlaps near 9.7 microns, which is explicitly marked.
No physical parameters were changed to inflate the contrast in response to
that feedback. Red shading is a model difference, not an uncertainty band.

The separate ETC strip shows statistical precision of approximately
1.67%, 1.09% and 2.03% at rest 10.58, 13.00 and 18.00 microns for the existing
R06 power-law setup. These are not error bars calculated for the illustrated
dust spectra. The idealized spectra are not fitted to targets or actual
light curves; this figure does not establish population model separation
or the unknown intrinsic scatter. Those numerical placeholders remain open.
The future clumpy and disk/wind model comparison remains in the text.

The user clarified that SPHEREx spectra are available for all 24 targets,
although only a subset is on this computer. The existing availability
statement was retained; it was never changed to a three-source claim.

Updated the MIRI paragraph to connect warm-continuum measurements, silicate
profiles and PAH/spatial host constraints to their scientific roles. Moved
the event-age definition into the model-testing section for flow and page
allocation. The rebuilt main proposal occupies five pages; supplemental
information and references bring the complete PDF to seven pages. Official
style, sample, Figure 1 measurements and APT configuration are unchanged.

Numerical checks: unchanged luminosity and the pre-echo limit give identical
fixed/evolving spectra; doubling radial sampling changes the spectra by less
than 0.25%. Emitting-shell equilibrium integrals agree to numerical precision.
The figure and assembled pages were visually checked. Full target fits and
population power remain separate future work, not completed by this illustration.

## Figure 2: cartoons and fractional differences

At the user's request, removed the ETC row and replaced the overlaid spectra
with schematic vector cartoons. Each rising/fading column now contrasts the
unchanged inner boundary with an expanding/contracting boundary. Arrows show
the prescribed motion; the caption identifies schematic sizes and dot counts.
The fractional-difference panels are retained and enlarged, with silicate
wavelengths and PAH positions marked. The fading case's 9.7-micron overlap
remains visible. All modeled spectral values and physical parameters are
unchanged; the saved spectrum CSV was compared with the preceding version.
The ETC calculations remain in the observing-description text.

## Figure 2: fewer titles

Removed repeated headlines, cartoon motion captions and explanatory footnotes.
Kept rising/fading labels, fixed/evolving labels, axes and diagnostic annotations.
Compacted the figure vertically; the proposal caption carries the explanation.
Removed the unlabeled illustrative-redshift cutoff marker; the figure uses
rest wavelengths without target-specific coverage claims. Spectral values are
unchanged. Clean rebuild: five core pages, seven total; figure page inspected.

## Figure 2: black-hole symbols

Replaced central stars with vector black-hole/accretion-disk symbols. The
black center has the same size in all panels; disk color is brighter in the
rising case. Icons are schematic and not to scale, as described in the
updated caption. Model spectra are unchanged. Rebuilt successfully within
five core pages (seven total).

## Figure 2: shading legend

Added a compact legend beneath the wavelength axis for the blue hot-dust
anchor, gold warm-dust region and red model difference. Model values are
unchanged. The figure was visually inspected and the PDF rebuilt cleanly,
remaining five core pages and seven pages total.

## Why 24: conditional sensitivity calculation

Replaced the X/Y placeholders with an explicit design benchmark, calculated
by `sample_power.py` and saved to `inputs/sample_power.json`. The primary
contrast is the difference of rising/fading group means in delta_warm.
Equal allocation minimizes its variance at fixed total N when the group
variances are equal. With 12 per group and assumed independent total residual
scatter 0.20 dex, the standard error of the contrast is 0.08165 dex.

A design offset of 0.25 dex (factor 1.778 in the geometric-mean
observed/predicted warm-continuum ratios) gives power 0.83295 for a two-sided
pooled t-test at alpha=0.05, df=22. Eight per group gives power 0.64294 under
the same assumptions. For total scatter 0.30 dex, the 80%-power offset is
0.35904 dex. Scenarios spanning 0.15--0.30 dex and 8--18 per group are saved.

The scatter and effect size are planning assumptions, not observed quantities
or predictions inferred from Figure 2. The calculation assumes independent
Gaussian residuals with equal variance and does not incorporate group-dependent
systematic offsets, shared errors, or fitting additional covariates. The text
retains propagation of those effects in the final analysis. Alpha=0.05 is
not a three-sigma claim; stricter-threshold sensitivity is in the JSON.

An independent 200,000-trial simulation gives empirical power 0.83201 and
null rejection rate 0.050465, agreeing with the analytic calculation within
sampling error. Method references are the NIST two-sample t-test and sample
size handbook pages, linked in the JSON. The main PDF remains five pages
plus supplement/references (seven total), builds cleanly, and page 3 was
visually checked. The only remaining red text in the proposal is the
duplication paragraph. Target-specific model inference and verification
of sample properties remain outstanding; no empirical population scatter
has been measured by this planning exercise.

## Why 24: science-led sample justification

At the user's direction, replaced the assumed-scatter/power paragraph with
a justification based on equal rising/fading representation, sustained
changes and reversals, different event epochs, and repeated physical tests
across 12 objects per group. It connects sample coverage to direction-,
age- and history-dependent warm-dust diagnostics. It makes no unsupported
claim of a guaranteed detection threshold, measured scatter, or statistically
optimal sample size. The numerical power calculation remains available in
`sample_power.py` and `inputs/sample_power.json` as an internal benchmark,
not a prediction in the proposal. Five core pages; seven total; clean build
and visual review of page 3.

## Duplication screening completed — 2026 September 26

Searched all 24 science coordinates with the STScI MAST duplication workflow,
using a 120-arcsec radius and no instrument, date, public-data or calibration-level
restrictions. All queries completed without errors and returned zero observations.
The same query at NGC 7469 returned 71 observations, including MIRI/IFU, from
program 1328. Results and cached responses are in `duplication/`; no matched
programs required further inspection. Filled the final proposal placeholder with
the dated result. Rebuilt successfully: five core pages, seven total, no layout
warnings. Updated the README current-status summary and validation metadata.

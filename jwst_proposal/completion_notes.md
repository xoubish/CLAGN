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

Local setup check: no Pandeia module or configured reference data found in the
active Python environment, and no JWST ETC/APT project found in the workspace.
The installed APT application is 2025.5.3 and needs a Cycle 6 update before
final planning/submission. The supported Cycle 6 Pandeia engine is 2026.7,
with matching reference data and PSFs. An existing ETC workbook/setup has been
requested from the user before choosing the calculation route.

Sources:
- https://outerspace.stsci.edu/spaces/PEN/pages/77530136/Pandeia%2BEngine%2BInstallation
- https://jwst-docs.stsci.edu/jwst-mid-infrared-instrument/miri-performance/miri-sensitivity
- https://jwst-docs.stsci.edu/jwst-mid-infrared-instrument/miri-observing-modes/miri-medium-resolution-spectroscopy

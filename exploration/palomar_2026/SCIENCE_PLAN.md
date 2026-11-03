# Scientific program after the Palomar campaign

This plan is grounded in the [Palomar science notes](../../observing/science_review/README.md)
and the [current JWST proposal](../../jwst_proposal/proposal.tex), as requested.
It is a discovery and measurement program, not a claim that every listed analysis
has already been performed. Executed measurements and initial results are in the
[notebooks](README.md). The original Palomar proposal was not used.

## 1. Deliver the promised selection experiment

The local [Hemmati et al. 2026 paper](../../reference/papers/Hemmati_2026_ApJ_998_130.pdf)
motivates using a frozen light-curve manifold to select unusual AGN states.
Palomar now provides a prospective spectral endpoint. The essential result is
whether that representation adds predictive information beyond simple light-curve
and observational variables, not whether a few interesting sources can be found.

Preserve the pre-observation A/B/C/D groups from the science notes: manifold with
a dated trigger, manifold without one, variability outside the selected regions,
and comparison sources without a dated trigger. Separately label catalog-confirmed
CLAGN, literature candidates, ambiguous associations and ancillary targets. The
present observed counts are 11, 59, 2 and 9 in these pre-review groups, with 18
targets outside that review. These labels are not the result of the new spectra.

The primary outcome should be reviewed broad-Hβ state change relative to a
usable archival baseline. Separate ordinary line variability, extreme variability,
a credible appearance/disappearance, a possible recurrence and unusable evidence.
Do not turn a wing-EW ratio threshold into a changing-look definition. The
multi-epoch approach of [Zeltyn et al. 2024](https://arxiv.org/abs/2401.01933)
is relevant to the distinction between CL and other extreme-variable systems.

**Required implementation:**

1. Freeze object identities, pre-run features, archival reference policy and the
   actual observed denominator. Reconcile primaries, backups, requested galaxies,
   protected JWST sources and weather/time substitutions.
2. Validate object-level outcomes with native spectra, common instrumental
   resolution, calibration/host modeling and independent review where possible.
3. Inject broad lines of different widths and fluxes into each spectrum, using
   the real noise, masks and continuum/host mix. Recover them with the same fit
   and review rules. Record non-detection limits and fit failures.
4. Compare groups at common redshift, archival brightness, host fraction,
   line detectability, time baseline and photometric coverage. Report covariate
   overlap and effective sample size before adjusting for many predictors.
5. Model object-level outcomes, allowing repeated visits within an object.
   Bootstrap by source; do not resample all epoch pairs as independent objects.
   Include measurement uncertainty, ambiguous labels and censoring.
6. Compare manifold features against simpler predictors, such as IR amplitude,
   slope and optical change, using object-level held-out evaluation if the sample
   supports it. Reusing training identities or refitting the manifold with the
   new outcome labels invalidates an independent validation claim.

**Deliverable:** observed-sample enrichment and incremental predictive value,
with uncertainty and sensitivity to matching/censoring. A population incidence
rate requires a reconstructable parent selection function and inclusion
probabilities; these are not supplied by the observed sample alone. With only
two C-group sources, complex adjusted subgroup inference may simply be too weak.

## 2. Make the multi-epoch spectral measurements reliable

The strongest reusable asset will be a catalog of consistent components and
their uncertainty across all epochs, rather than a set of visually interesting
overlays. Current notebooks deliberately stop at signed wing indices and
total-continuum colors.

**Calibration and extraction first:** inspect raw-exposure repeatability, sky
residuals, telluric-corrected versus uncorrected profiles, standard choices,
channel overlap and slit/acquisition geometry. September's standard disagreement
and October 27's cloud losses must enter as nuisance parameters or sensitivity
families, not invented Gaussian uncertainties. Isolated peaks should survive
individual-exposure checks. Features at joins should be measured in native
channels. Verify wavelength frames and flexure before interpreting centroids.

Use synthetic photometry through actual passbands only where throughput coverage
is complete. A continuum value at a filter pivot is not a synthetic magnitude.
Match measured total photometry in time and account for aperture/host light.
Treat narrow lines as useful calibration checks only after testing their spatial
extent and aperture sensitivity. The combination of narrow-line,
pseudo-photometric and visual checks in the
[DESI CLAGN study](https://arxiv.org/abs/2408.00402) motivates this redundant approach.

**Joint fit specification:**

* Correct Galactic extinction using a documented map/law; maintain the original
  spectra and correction provenance. Use vacuum rest wavelengths consistently.
* Fit a shared host spectral shape where justified, with aperture-dependent
  amplitudes; add an AGN continuum and appropriate optical/UV Fe II templates.
  Low states help constrain host light, but are not assumed to be pure host.
* Fit narrow Balmer, [O III], [N II], [S II] and other covered narrow lines with
  justified kinematic/doublet constraints. Allow extended narrow emission to
  change its observed contribution across apertures.
* Fit broad Hβ, Hα and Mg II flexibly enough to capture asymmetry without forcing
  a component in an off state. Compare a broad-component and no-broad-component
  model with calibration-aware, injection-validated decision rules. Boundary
  problems make a naive likelihood-ratio significance unreliable.
* Convolve the intrinsic model with each instrument's wavelength-dependent
  resolution. Match grids through the model, retain resampling covariance or
  an empirically justified correlation model, and test sensitivity to LSF errors.
* Fit grey/color calibration terms and exposure/night effects jointly with
  physical components. Anchor them only with valid contemporaneous information.

Prior-informed host decomposition is a useful starting point, as described by
[Ren et al. 2024](https://arxiv.org/abs/2406.17598), but template and prior sensitivity
must be evaluated on these sources rather than assumed away.

**Save:** broad/narrow fluxes and EWs, FWHM and line dispersion, centroids and
asymmetry, continuum luminosity/slope, host fraction, Fe II strength, covariance
or posterior draws, upper limits, model components, residual plots, fit flags,
LSF/calibration assumptions and provenance. Avoid black-hole-mass/Eddington-ratio
claims until the luminosity, width and host corrections are credible.

**Acceptance:** repeat-exposure consistency, injection/recovery completeness,
stable qualitative conclusions across supported calibration/host choices,
reasonable residuals and an adjudicated fit-quality flag. Missing coverage is
not a non-detection. A weak line should yield a meaningful upper limit.

## 3. Push individual-source discovery beyond binary state changes

The [atlas](outputs/all_sources_atlas.pdf) gives every observed object a place
in the review. A review queue based only on extreme ratios would discard useful
recurrent, slowly evolving or spectrally unusual systems.

| Track | Candidates or data subset | Measurement and discriminant |
|---|---|---|
| Coherent Balmer weakening | P21422 first; P24930 | Joint Hβ/Hα fits, broad/narrow behavior, host and optical-state consistency; verify literature status. |
| Broad-profile strengthening | P7947; known-source P17881 | Separate nuclear broad flux from continuum/host dilution; the P17881 [O III] normalization disagreement is a key check. |
| Month-scale evolution | All five repeated sources | Joint native-spectrum fits with night-specific calibration; compare changes over ~33–34 observed days in the source rest frame. |
| Persistence vs new transition | P9694, P11530, P1823 | Distinguish continued old behavior from a newly changed spectral state; retain all usable archival epochs. |
| Differential Mg II response | Eight index-screen leads plus sources with both lines | Optical/UV Fe II decomposition and response relative to Hβ at the same epochs. |
| Broad-profile asymmetry/centroid shifts | All high-S/N multi-epoch sources | Instrumental/wavelength checks, velocity-resolved difference spectra, reverberation/profile variability alternatives; no binary-black-hole inference from one shift. |
| Balmer-decrement evolution | Low-z sources with both broad Balmer lines | Broad-only ratio with cross-arm/color calibration and host absorption; compare accretion/extinction hypotheses. |
| He II/Bowen/coronal changes | Covered high-S/N sources | A systematic blind line search with sky/blend masks and multiple-testing control; transient/outflow explanations need follow-up. |
| Dusty-star/blend anomaly | H5MUSES-327 | Spatial association and optical stellar subtraction; keep dusty-star and background-source possibilities open. |
| Ancillary galaxies | HIC-1953, HIRAS-F23290-1252 | Stellar continuum and narrow-line fitting, aperture checks, censored diagnostic line ratios; do not force into CLAGN incidence statistics. |

These are proposed tracks, not completed measurements. For each candidate,
record the physical alternatives and the one observation or fit that would
discriminate between them. Literature searches must use resolved coordinates,
aliases and a specific event epoch; “no match in checked catalogs” is not proof
of novelty.

## 4. Exploit time-domain structure without overinterpreting cadence

**Optical endpoint recovery is the highest-impact acquisition task.** Local survey
curves commonly stop in 2025. Saved alerts extend some histories into 2026, but
only one Palomar visit passes the ±30-day unbinned-data match. Retrieve all
available total-flux measurements around each visit, with field/object-ID,
reference-flux and association audits. Retain signed low-S/N fluxes and upper
limits. Do not infer stability from absent alerts, and do not combine reference
systems without a measured transformation.

With reliable measurements, prioritize:

* **Line responsivity:** within-source broad-line flux versus AGN continuum,
  accounting for temporal offsets, measurement covariance and luminosity-state
  dependence. A trend in EW versus continuum is partly algebraic and insufficient.
* **Breathing:** line width/dispersion versus AGN continuum using matched
  instrumental resolution and joint posterior draws. Separate profile-component
  changes from a smooth width response.
* **Color and hysteresis:** matched g−r and host-corrected spectral slope at the
  same luminosity on rising and falling branches. Use actual time order and
  explicit host/line-contamination alternatives.
* **Recurrence and timescales:** reconstruct a spectral state sequence. Treat
  transition times as intervals between bracketing observations. Prefer a small
  set of secure recurrences over a large count inferred from photometric peaks.
* **Optical–IR mismatch:** compare dated optical nuclear states with delayed IR
  histories. A mismatch can be useful without claiming a unique lag or geometry.
* **Lag recovery:** use irregular-cadence transfer models only on suitable
  overlapping data. Inject known lags and red-noise drivers into the observed
  cadence, test seasonal aliases and held-out seasons, and report nonidentifiable
  lags explicitly. WISE's cadence does not support arbitrary short-timescale claims.

For ensemble trends, use a hierarchy that separates within-object changes from
between-object luminosity, host and redshift differences. Avoid correlations
between noncontemporaneous extrema from different surveys. A matrix of hypotheses,
data eligibility and null tests should precede an unrestricted feature search.

## 5. Use Palomar to strengthen the JWST experiment

The proposal asks how long warm-dust memory persists and whether fading/recovering
sources need dust replenishment beyond a fixed distribution responding to the
recorded illumination. Palomar supplies a new optical spectral endpoint and
constraints on nuclear versus host emission. It does not by itself measure
MIRI silicate evolution or demonstrate changes in dust geometry.

Update every proposal target with the actual 2026 state, broad-line and continuum
posterior, host fraction, last reliable optical driver point, and remaining
calibration/association uncertainty. Keep the proposal's historical rise/fade
family separate from the new optical-state measurement. There are delivered
spectra for **23/24** targets; identify an alternative current-state measurement
for **P24976** or obtain a replacement observation. P11530 already has September.

Build a forward comparison of fixed-dust illumination response and evolving
dust geometry. Marginalize over unobserved driver intervals, host/PAH contribution,
lag/kernel shape, line-of-sight extinction and spectral calibration. Use the
existing SPHEREx extraction-quality flags and variable-epoch uncertainties;
do not promote a scenario curve into observed data. Fit actual MIRI continuum,
silicates and PAHs when available, with spatial host information where measurable.

**Deliverable now:** an updated 24-source optical/IR endpoint package and explicit
forecast identifiability limits. **Deliverable later:** a comparison of dust-response
hypotheses using the IR spectra, with model degeneracies and selection effects.

## 6. Execution order and review products

| Priority | Concrete work | Completion criterion |
|---|---|---|
| P0 | Review P21422/P24930/P7947/P17881, all five repeats, and Mg II leads with extraction/calibration diagnostics | Each lead has a stable interpretation or a named calibration/coverage failure. |
| P0 | Refresh optical data at all unmatched NGPS visits; locate P24976 endpoint | A dated match/upper limit or explicit acquisition failure for each visit/band. |
| P0 | Recover native DESI spectra/errors/LSFs and missing JWST archival spectra | Native provenance and usable-error/coverage flags; compact caches never given invented errors. |
| P1 | Joint fits on the priority set, then all eligible sources | Components, posteriors, limits, QA and successful injection/recovery checks. |
| P1 | Literature/identity audit and state adjudication | Separate known objects, new events, plausible new objects and ambiguous cases. |
| P1 | Frozen selection/denominator and completeness audit | Actual observation decisions and measurable comparison groups documented. |
| P2 | Response, breathing, recurrence and rare-line studies | Claims pass cadence/calibration null tests and multiplicity controls. |
| P2 | Population enrichment and JWST endpoint release | Source-level uncertainties, reproducible tables and explicitly limited scope. |

The per-source [work queue](outputs/work_queue.csv) implements this ordering at
screening level. Final priorities should be updated after visual adjudication;
the threshold is not a physical ranking score. Preserve an adjudication log with
reviewer, date, baseline, evidence, alternatives, state label and next action.

The first scientific outputs should be a validated epoch-measurement catalog,
a reviewed variable-source catalog and a measured test of the frozen selection.
These support individual case studies and the JWST dust-history analysis without
forcing every source into a changing-look narrative.

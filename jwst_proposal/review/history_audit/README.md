# Evidence from the 24 histories — 27 September 2026

## Assessment

The histories justify a sample of infrared-variable AGN with useful rises,
declines, and recoveries. They do **not yet establish a detectable 0.24-dex
rising–fading warm-continuum separation at the future JWST epoch**. A full-history
analysis is better motivated than treating the original 12/12 labels as current
nuclear states. Neither the measurements nor the response experiment below
establish that the programme will fail; the missing nuclear decomposition and
warm-response constraints prevent a reliable detection forecast in either direction.

The proposal text, target membership, exposure settings, and original sample
timing fields were not changed by this audit.

After the audit, an authorized editorial update added the measured monitoring
baseline and W1 percentile-ratio range to the proposal, clarified that the
rising/fading labels describe selection histories, and made the group offset a
summary comparison within the full-history analysis. No response-experiment
contrast was inserted as a detection forecast.

## Measured evidence

All 24 sources now have inspected ZTF g/r and NEOWISE W1/W2 photometry. Public
ZTF DR24 queries supplied eight previously uncached targets; four especially
sparse curves were also refreshed. The refreshed curves remain sparse under
the selection below. This is a photometric audit, not a spectral-state audit.

- Each source has 20–22 W1 and W2 observing seasons across 9.4–10.5 observed
  years (roughly 2014–2024).
- W1 95th-to-5th-percentile flux ratios span **1.30–2.29**, with a median of
  **1.65**. W2 ratios span **1.23–2.91**, median **1.77**. These are total-source
  ratios including host light, not bolometric nuclear amplitudes.
- Matched W1/W2 seasonal magnitudes have Spearman correlations **0.84–0.99**.
  Coherent variations across both bands strengthen the evidence for infrared
  variability, although this is not a systematic-error or origin-of-variability test.
- A descriptive linear slope over each source's final three observed W1 years
  is opposite to its selection label by more than 0.05 mag/year in **five** cases:
  F06, F11, R02, R04, and R06. This threshold is a diagnostic convention, not
  a fitted transition criterion or a statistical significance threshold. A final
  three-year average is also not the instantaneous direction at the last point.
- Optical curves provide additional counterexamples to persistent labels:
  **F04/P10381** brightened substantially after its earlier infrared decline;
  **R04/J082942.66+415436.8** faded after its 2019–2022 optical high state.
  Recoveries and recurrent flares must be retained in a forward model.
- Four targets have fewer than ten retained 90-day bins in both optical bands:
  **F01, R03, R06, R08**. In particular, the retained R03 optical bins stop in
  2023, despite later individual detections. Their current nuclear states require
  better photometry or spectroscopy; interpolated lines do not fix these gaps.

See [all 24 curves](history_atlas.pdf), [per-target measurements](history_summary.csv),
and [optical data inventory](optical_inventory.csv).

## What the data can say about future memory

`response_stress_test.py` calculates an intentionally simple sensitivity
experiment, not a fitted torus model. It forms

`history contrast = log10(filtered flux history / assumed current flux)`.

The inputs are the measured g, r, or W1 flux curves. Each curve is interpolated
linearly in flux, continued constantly before its first point and after its last
point, and smoothed with either an exponential or uniform delay kernel. The grid
uses mean rest-frame response times 0.5, 1, 2, 3, 5, 10, and 20 years, three
illustrative dates from July 2027 to June 2028, and two endpoint definitions.
The 2028 date is a scenario, not a scheduled visit.

Host sensitivity is explored by subtracting 0%, 50%, or 80% of each curve's
minimum observed flux as a constant. **None of those fractions is measured.**
Response times are not assigned to targets from their provisional sample-table
lags. W1 is already reprocessed infrared emission: further filtering of it is
only a descriptive history contrast, not a reconstruction of the optical/UV driver.

At 2028 January 1, the r-band experiment gives the following mean contrast
difference (original fading group minus original rising group), using an
exponential kernel and a last-400-day median endpoint:

| Assumed mean response time (rest years) | Host included | Subtract 50% of minimum flux | Subtract 80% of minimum flux |
| --- | ---: | ---: | ---: |
| 3 | 0.018 dex | 0.031 dex | 0.058 dex |
| 5 | 0.038 dex | 0.063 dex | 0.112 dex |
| 10 | 0.068 dex | 0.114 dex | 0.196 dex |

For the five-year, 50%-subtraction example, changing the date, endpoint, and
kernel gives 0.038–0.083 dex. Larger separations are possible in the explored
grid: a ten-year uniform kernel with 80% subtraction gives 0.248 dex. That
example has only about 32% of its kernel weight within the observed optical
time span; it relies strongly on the unobserved history. It is **not** evidence
for a measured 0.248-dex warm-dust signal.

The general problem appears at both ends of the response-time range. Short
responses can equilibrate during the gap before JWST if the source remains
steady. Long responses retain information from before the optical record
began. At five years the exponential example has median kernel weights of
about 49% inside the optical time span, 18% before it, and 34% after it. Even
the within-span weight includes interpolated gaps. These individually computed
medians need not sum to exactly one.

The proposal's 0.24-dex threshold is a conditional planning calculation for
its **warm-continuum observable**, assuming 0.20-dex independent residual
scatter. The filter contrasts are **different observables**. Comparing scales
shows why variability alone does not demonstrate the claimed sensitivity; it
does not calculate the programme's statistical power, set a lower or upper
bound on the real warm signal, or prove that MIRI cannot detect memory.

Missing ingredients include measured optical host fractions, optical-to-UV
heating conversion, nonlinear wavelength-dependent dust response, geometry,
covering-factor scatter, the warm equilibrium calibration, and the future
nuclear state. The user confirmed that the remaining SPHEREx spectra are not
on this computer. They were not assumed or simulated as measurements.

See [sensitivity plot](response_stress_test.pdf),
[full assumptions](response_assumptions.json), and
[group results](conditional_group_contrasts.csv).

## Consequences for the proposal

1. **The sample has an empirical justification.** Every target has sustained
   multiyear infrared variability, with coherent W1/W2 behaviour. The selection
   spans declines, rises, plateaus, and recoveries. This is stronger than an
   assertion based only on catalogue labels or manifold position.
2. **The simple two-group detection claim is not established.** Keep 12/12 as
   the historical selection balance, but use the reconstructed history and
   response stage as the scientific predictor. Recovering objects can dilute
   or reverse a binary contrast while remaining useful tests of history dependence.
3. **Sensitivity is not an expected effect size.** Do not present 0.24 dex as a
   forecast from these histories. Do not replace it with a favourable number
   chosen from the response grid. The assumed 0.20-dex scatter is also unmeasured
   for this sample and this analysis.
4. **Prioritize constraints, not additional citations.** Obtain nuclear optical
   decompositions and the available SPHEREx spectra, constrain the response
   families jointly with W1/W2, and compare history-dependent and equilibrium
   predictions with the same nuisance freedoms. Establish a shared, justified
   covering-factor/geometry distribution or an external comparison: freely
   fitting a separate warm normalization for each target can absorb the signal.
   MIRI-epoch monitoring is essential for measuring the current denominator.
5. **The earlier P9694 result remains relevant but limited.** Its admissible
   fixed and changing dust models can give nearly identical MIRI spectra under
   an assumed extra-scatter error model. That is a mechanism degeneracy, not a
   failure to detect history, and does not constitute a 24-object power forecast.

Candidate wording supported by the measurements:

> The sample is selected for measured changes in activity rather than a uniform
> spectroscopic label. Each target has approximately a decade of W1/W2 monitoring;
> the W1 95th-to-5th-percentile flux ratios span 1.3–2.3 before host subtraction.
> The histories include sustained changes and recoveries, providing tests of how
> infrared emission depends on the path to the current nuclear state. The rising
> and fading designations describe the selection histories; predictions will use
> the full multiband records, including reversals, rather than assume that every
> source remains in its original phase.

This paragraph justifies the selection. It does not resolve the quantitative
detectability question, which should remain explicit until the missing
decompositions and response constraints are available.

## Methods, limits, and reproduction

NEOWISE single exposures are matched within 3 arcsec, require positive frame
quality, clean artifact flags in the band used, and finite positive magnitude
errors. Seasons are separated by gaps over 60 days and require at least three
exposures. Seasonal magnitudes are medians. The reported standard error is
1.2533 times the within-season sample standard deviation divided by sqrt(N);
this does not include a calibration floor or covariance. No variability
significance claim is derived from that error alone. Separating seasons by gaps
avoids splitting a visit with an arbitrary 180-day boundary.

ZTF detections require finite positive errors, clean bit 32768, and separation
under 1.5 arcsec. One object ID per band is chosen by number of populated
90-day bins, then baseline, then detection count. Retaining a single object ID
avoids mixing field-dependent offsets but can discard useful overlapping-field
measurements. Bins require three detections. No forced-photometry upper limits
are available in these queries; faint states may be incompletely sampled.
No optical/infrared constant host component is measured in this audit.

The period-averaged metrics use the first/last 400 days, percentile amplitudes,
and an unweighted last-three-year slope. They do not identify spectroscopic
changing-look dates. The original sample's factor-of-ten warm/hot delay rule is
explicitly provisional in `inputs/selection_summary.json` and was not used as
a measured timing constraint.

Run from the project root:

```sh
MPLCONFIGDIR=/tmp/clagn_mpl python jwst_proposal/review/history_audit/audit_histories.py
MPLCONFIGDIR=/tmp/clagn_mpl python jwst_proposal/review/history_audit/response_stress_test.py
```

`provenance.json` records input hashes. The response code checks its constant-input
limit and its piecewise-linear exponential solution against numerical integration.

Relevant primary literature:

- [Hönig & Kishimoto 2011](https://arxiv.org/abs/1109.3465): response times and
  smoothing depend on the radial brightness distribution and wavelength; this
  does not assign a universal warm/hot lag ratio.
- [Lyu, Rieke & Smith 2019](https://arxiv.org/abs/1909.11101): W1/W2 reverberation
  and much weaker longer-wavelength variability in their quasar sample; their
  results are not a measurement of this sample's warm signal.
- [Almeyda et al. 2020](https://arxiv.org/abs/2002.12823): geometry-dependent
  dust response functions.

## Target-by-target reading

The notes below describe observed photometric shapes. They are not new
spectroscopic classifications or fitted event dates.

| ID | Target | W1 flux ratio (95th/5th) | Reading and limitation |
| --- | --- | ---: | --- |
| F01 | P14388 | 1.45 | Clear later decline in optical and W1/W2; sparse optical sampling limits the current-state constraint. |
| F02 | P9506 | 1.53 | Long infrared decline; later optical fluctuations and a renewed 2025 rise complicate a single fading state. |
| F03 | P24976 | 1.77 | Recurrent behaviour: infrared high state around 2022 followed by a decline; not a single monotonic event. |
| F04 | P10381 | 1.59 | Earlier infrared decline, followed by substantial optical recovery through 2025. Persistent fading is inappropriate. |
| F05 | CLAGN_0816 | 2.10 | Strong infrared decline after 2018; recent optical rise warrants a current-state update. |
| F06 | P11530 | 1.90 | 2018–2019 optical flare, decline, and 2022–2023 recovery; full history is essential. |
| F07 | P7837 | 1.30 | Most sustained infrared decline precedes the ZTF record; later optical contrast is relatively small. |
| F08 | CLAGN_0826 | 1.79 | Infrared decline and subsequent recovery; recurrent structure rather than a persistent monotonic fall. |
| F09 | CLAGN_0987 | 1.80 | NGC 2617: optical/infrared decline into the 2022–2025 low state; useful fading history. |
| F10 | J085337.26+014303.6 | 1.57 | Earlier infrared flares with comparatively weak recent optical change; pre-ZTF driver is missing. |
| F11 | CLAGN_0795 | 1.72 | Optical and infrared recur on different apparent timescales; optical high state around 2022 and later decline. |
| F12 | CLAGN_0801 | 1.48 | Much of the infrared decline is early; optical fluctuations and recent recovery complicate the label. |
| R01 | P9694 | 1.38 | P9694: optical rise in 2021–2022 with infrared brightening; useful rising history and available 2026 spectrum. |
| R02 | CLAGN_0215 | 1.45 | Rise to a 2020–2021 infrared high state followed by later fading; not a persistent rise. |
| R03 | P24541 | 1.94 | Strong infrared rise after the 2020 trough; promising contrast but sparse optical record and late-state gap. |
| R04 | CLAGN_0030 | 2.29 | Optical high state in 2019–2022, then fading; the original rising label can predict the wrong sign. |
| R05 | CLAGN_0824 | 1.51 | Overall optical/infrared brightening with substantial intervening optical fluctuations. |
| R06 | P17881 | 1.59 | Rise and renewed infrared decline; sparse optical bins make a simple current-state label unreliable. |
| R07 | CLAGN_0045 | 2.26 | Optical rise toward 2021 followed by W1/W2 rise toward 2023; strong example of a delayed infrared response, without a fitted lag. |
| R08 | CLAGN_0818 | 1.34 | Infrared rise after 2018; optical sampling is sparse and the two bands have different late coverage. |
| R09 | CLAGN_0969 | 1.73 | Rise into a brighter state around 2020, then fluctuations and a late infrared decline. |
| R10 | CLAGN_0972 | 1.45 | Infrared recovery from a trough with weak net early-to-late change; recent optical amplitude is modest. |
| R11 | CLAGN_0808 | 2.04 | Large infrared event peaks before the ZTF record; later optical recovery does not reconstruct the missing early driver. |
| R12 | CLAGN_0973 | 2.16 | Strong decade-long infrared rise; comparatively modest later optical changes and recurrent dips. |

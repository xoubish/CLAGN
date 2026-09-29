# Sample reassessment — 29 September 2026

## Recommendation

Do not replace the current 12+12 count with 10+10, 12+12 or 15+15 yet. The
available evidence supports a **first modelling shortlist of 12 objects: seven
fading and five rising histories**, consisting of four current targets and eight
possible replacements. This is a validation worklist, not a proposed observing
sample or a claim that all twelve will retain measurable warm-dust memory.

Among the existing 24, prioritize four for detailed modelling, hold seven for
additional checks, and exclude thirteen from a *simple sustained-change core*.
Exclusion here does not mean scientifically uninteresting or unsuitable for a
separately motivated recovery/flaring experiment. The current source CSV, APT,
figures and proposal text have not been changed by this assessment. No proposal
PDF was generated.

A cleaner sample is plausible, but it cannot be obtained by deleting only the
obvious 2019 flare and keeping the other labels. The top and middle examples in
Figure 1 both belong outside the proposed clean core. The bottom example has a
retained rise but a recent optical downturn.

## Current targets to model first

| ID | Object | Interpretation | Main outstanding issue |
| --- | --- | --- | --- |
| F01 | P14388 / J075934.95+322143.3 | Sustained optical decline after 2018–2019 and later infrared fading | Optical sampling is sparse; date the transition and establish current state |
| F09 | NGC 2617 | Optical and infrared decline into the 2022–2025 low state | Earlier fluctuations and strong host sensitivity in the old shell forecast |
| R01 | P9694 / J004607.98+090720.9 | Retained rise, with optical decline from the 2025 high into 2026 | Do not describe as monotonically or currently rising; use the full history |
| R07 | J101152.98+544206.4 | Optical rise to a high plateau with delayed W1/W2 brightening | Determine whether useful warm memory survives after the plateau |

These four are the first objects to inspect and model, not four confirmed
monotonic sources. F01 and F09 have broad earlier high states; the distinction
between a long event and a persistent transition is observationally limited.

**Hold:** F02, F05, R02, R03, R05, R08, R12. Reasons include recent recovery,
an infrared downturn, sparse optical sampling, and uncertain nuclear contrast.

**Exclude from the simple sustained-change core:** F03, F04, F06, F07, F08,
F10, F11, F12, R04, R06, R09, R10, R11. These histories are recurrent, have
reversed, or lack the optical record of the main transition.

The complete source-by-source reasons and required checks are in
[current_target_assessment.csv](current_target_assessment.csv). Descriptive
measurements are in [current_metrics.csv](current_metrics.csv). An endpoint
infrared cut alone is inadequate: F02 and F03 pass that cut despite their
optical recovery/recurrent behaviour, whereas R01 fails the conservative
amplitude cut despite a useful recorded rise.

## Replacement search

The search examined the **2,376 objects in the saved W3-matched forecast pool**,
using its coordinates and photometric metadata, not its forecast significance
to choose targets. This pool was previously selected for z<0.3, infrared
variability, observing coverage and an AllWISE match. It is a subset of the
larger parent catalogue, not an exhaustive search of all nearby AGN.

Of the objects outside the current sample:

- 428 passed the initial two-band infrared screen: 275 fading, 153 rising.
- 34 also passed the optical-direction check: 20 fading, 14 rising.
- 65 had an optical conflict or reversal under the screening criteria.
- 329 had weak, insufficiently recent, or missing local optical evidence.
- Visual inspection of all 34 supported candidates prioritized **five fading
  and three rising replacements**. Another 23 remain reserves; three have
  clear flare/recurrent histories unsuitable for the intended clean core.

The full screening tables are [parent_wise_screen.csv](parent_wise_screen.csv)
and [replacement_candidates.csv](replacement_candidates.csv). Visual decisions
are recorded for all 34 in [replacement_assessment.csv](replacement_assessment.csv).
These categories are judgement-based triage, not probabilities or complete
scientific classifications.

| Candidate | Direction | z | W1 late-minus-early (mag) | W3 (mJy) | Catalogue spectral epochs |
| --- | --- | ---: | ---: | ---: | ---: |
| P20778 | Fading | 0.2471 | +1.337 | 2.97 | 2 |
| P14318 | Fading | 0.1359 | +0.589 | 2.05 | 2 |
| P18762 | Fading | 0.1540 | +0.709 | 1.24 | 3 |
| P16096 | Fading | 0.1052 | +0.451 | 5.50 | 2 |
| P16524 | Fading | 0.2267 | +0.361 | 2.20 | 1 |
| P12673 | Rising | 0.1897 | −0.728 | 4.26 | 5 |
| P7750 | Rising | 0.1571 | −0.842 | 12.16 | 1 |
| P19631 | Rising | 0.2776 | −0.333 | 3.46 | 2 |

Coordinates and all measured screening fields are in
[priority_replacements.csv](priority_replacements.csv). Positive magnitude
change means fading. Endpoint values are medians of the first and last 400
days of available NEOWISE coverage, approximately 2014–2024, and include host
light. Spectral counts come from the saved master catalogue; they do not
establish usable broad-line changes or independently calibrated epochs.
W3 values are archival total-source fluxes, not future nuclear predictions.

See the actual four-band curves:

- [Priority candidates, page 1](priority_histories_1.png)
- [Priority candidates, page 2](priority_histories_2.png)
- [All 34 reviewed candidates](replacement_visual_review.csv), with the first
  24 plotted in `replacement_histories_1.png` through `replacement_histories_4.png`
  and the remaining ten in `additional_candidates_1.png` and `additional_candidates_2.png`.

The strongest illustrative replacements include P20778's prolonged decline,
P12673's rise into a retained high state, and P19631's substantial optical rise
with continuing infrared brightening. Normal fluctuations are allowed: none is
asserted to be mathematically monotonic, and P7750 has a later optical dip.

### Follow-up implications

None of the 34 supported replacements matches the checked proposed October
NGPS queue within 3 arcsec; none has a matching name among the checked local
September reduced NGPS spectra. This is a statement about the checked local
products, not proof that no other spectroscopy exists. SPHEREx completeness,
extraction quality and repeat-visit availability have not been verified for
these replacements. P16524 and P7750 have only one epoch in the saved master
spectral inventory. The original claim that every target has the full rich
dataset must not be transferred automatically to a replacement list.

Some candidates are fainter in W3 than the present faintest target. New nuclear
flux estimates, ETC calculations, sky selection, duplication screening and
APT timing are needed after membership is selected; current 67 h cannot simply
be inherited. October visibility must be assessed for replacements rather than
assuming that a coordinate absent from the queue is an available night target.

## Does removing flares solve the timing problem?

**No. It improves interpretability, but does not guarantee a surviving signal.**
An additional sensitivity calculation applies simple exponential and uniform
delay kernels to the measured optical flux histories, using rest-frame mean
response times of 1, 3, 5 and 10 years and illustrative dates from July 2027 to
June 2028. These dates are not scheduled observations.

The full results are [timing_sensitivity.csv](timing_sensitivity.csv), with a
[summary for the twelve priorities](priority_timing_summary.csv). For example,
in the January 2028 scenario, the signed r-band history contrast for P20778 is
about 0.013–0.025 dex for a 3-year response and 0.053–0.066 dex for a 5-year
response. For P19631 the corresponding values are about 0.065–0.087 and
0.109–0.150 dex. Short responses can produce negligible contrast even for
the visually clear histories; large assumed delays increasingly depend on
unrecorded illumination before the optical survey.

**These numbers are not MIRI flux predictions, detection significance, measured
lags, or fitted nuclear responses.** They use total-source optical flux, no
host subtraction or optical-to-UV conversion, a unit-gain linear response, and
constant illumination before/after the record. A last-400-day endpoint is used.
They do not model silicate spectra, composition, geometry, or covariance. The
purpose is to expose timing dependence rather than rank targets by a favourable
assumed dust timescale. Visual classifications were assigned independently of
these contrasts.

The updated alert photometry available for F02, F04, F06, F07 and R01 was
inspected in [current_recent_updates.png](current_recent_updates.png). For the
timing calculation these five use 90-day bins from the Figure 1 IRSA+alert
assembly, retaining bins with at least three nights. Other current targets use
the prior audit's single-ZTF-object bins. Sparse alerts do not necessarily form
a new bin; `last_observed_year` records the actual usable endpoint. R01 uses g
in the summary because its well-sampled extension is more recent than r. The
alert assembly's object handling differs from the original audit; these remain
sensitivity calculations, not final scientific light-curve fits.

The existing shell-model predictions for the current sample are also retained
in the assessment CSV as signed delayed-minus-instantaneous contrasts at
0.5, 1 and 2 times the old dust scale. They assume a W1-derived heating proxy
and constant future illumination; the raw contrasts, not the old
`memory_signal` weighting, are used. None of the four prioritized current
objects reaches 3 sigma individually for that memory comparison across those
saved scale choices and quoted errors. This does not rule out an ensemble or
full-spectrum test, but rules out claiming that their appearance alone proves
detectability.

## Selection decision and next analysis

1. Start comparable joint history/SED fits with the four current and eight
   replacement priorities; verify the seven current holds before excluding
   them definitively. Maintain reserves instead of filling a numerical quota.
2. Fit current illumination, host fractions, plausible dust geometries and
   response scales using the available optical, infrared and spectral data.
   Marginalize missing prehistory and plausible illumination up to JWST.
3. Evaluate the full continuum and silicate predictions at candidate observing
   dates. Compare equilibrium and delayed models, then fixed and evolving dust,
   under the same nuisance freedoms. A broad-band endpoint difference is not
   the complete science criterion.
4. Choose the final count using information gained per observing hour, group
   comparability (luminosity, redshift, host fraction, amplitude and event age),
   and actual follow-up readiness. Equal counts are optional. Do not quote the
   original 24-object ensemble significance for a revised subset.

**A clean 10+10 programme is a design to investigate, not a result of this
audit.** The first worklist has 7 fading and 5 rising objects, and even those
require the checks above. The broader pool offers reserves, but committing to
20 or 30 targets now would exchange the original history ambiguity for
unverified data completeness and detectability.

## Reproduction and scope

Run from the repository root:

```sh
/opt/anaconda3/bin/python jwst_proposal/review/sample_reassessment/screen.py
/opt/anaconda3/bin/python jwst_proposal/review/sample_reassessment/assess.py
```

The infrared screen requires at least 15 visits in each band, same-sign
endpoint changes of at least 0.30 mag in W1 and W2, at least 65% retention of
the smoothed excursion, terminal retreat at most 0.15 mag, and no opposing
last-three-year slope exceeding 0.05 mag/year. Smoothing uses a three-visit
rolling median; endpoint medians can slightly exceed that smoothed span.
Thresholds are descriptive choices, not significance cuts. Their variation is
tabulated in [screen_threshold_sensitivity.csv](screen_threshold_sensitivity.csv).

For the optical screen, use cached detections within 1.5 arcsec, remove flagged
photometry, choose a single object ID per filter by coverage, deduplicate
exposures, and require three detections per 90-day bin. At least one band must
have ten bins, reach 2025, change by at least 0.20 mag in the selected direction,
and have terminal retreat no larger than 0.20 mag. An informative band changing
oppositely by more than 0.15 mag or with an opposing recent slope above 0.05
mag/year triggers review rather than support. Visual review then checks the
full paths, including peaks that endpoint cuts miss.

No new remote photometry or spectra were acquired. The parent WISE visit
tables were used as saved; their underlying per-exposure filtering was not
repeated for every candidate. Cached blend flags were recorded but do not
replace image inspection. The remaining 329 candidates without adequate
optical support were not classified as physically unsuitable. Input hashes
are in [provenance.json](provenance.json) and
[recent_alert_provenance.json](recent_alert_provenance.json). Numerical and
bookkeeping checks are in [assessment_checks.json](assessment_checks.json);
their passing status does not validate the astrophysical assumptions.

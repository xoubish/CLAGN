# Figure 2: scientific audit and replacement requirements

Checked against the active `make_fig2_miri.py`, its saved scenarios, the P9694
predictivity audit and the ensemble forecast on 2026-09-29. No new physical fits
were performed. The active figure and proposal source were not changed.

## What the current figure establishes

The selected blue and red spectra are both fixed-distribution models, at the
same illustrative 2028 epoch and future driver factor. Their integrated 8–13
micron fluxes differ by a factor of 14.52. This demonstrates that these
exploratory historical fits leave substantial warm-spectrum freedom. It does
not demonstrate that this freedom is caused by variability, nor separate
current-power equilibrium, delayed fixed dust and evolving dust.

The lower panels re-express those spectra as flux ratios and local silicate
contrasts. Grey symbols are selected conditional solutions, not probability
distributions or confidence regions. Fixed and evolving predictions overlap.
There are no forecast measurement errors or recovery results in these panels.

The highlighted solutions come from an optically thin shell basis whose
historical fits fail the primary error model. They become acceptable under an
explicit additional 5% independent-scatter assumption. The extreme silicate
contrast should not be presented as a validated radiative-transfer prediction.

## What the existing checks warn against

The saved P9694 audit contains fixed-distribution fits that reproduce two
evolving-distribution forecast spectra while retaining the adopted historical
constraints. Independent angular-quadrature checks give maximum spectral
differences of 0.845% and 0.436%. These are conditional counterexamples to
unique mechanism identification, not a forecast that all sources or models
will be indistinguishable, or a claim of subpercent physical accuracy.

The ensemble's 4.5 and 3.7 sigma numbers divide mean model offsets by an
assumed measurement/host error. They are not results of fitting competing
families to mock spectra with nuisance parameters varied. The baseline has
zero individual >3 sigma delayed-versus-instant or evolving-versus-delayed
separations; the five objects quoted in the text exceed 3 sigma for the widest
pair of the three predictions. Halving the response scale changes the
conditional ensemble memory estimate from 4.5 to 2.7 sigma. Neither statistic
includes a general exploration of geometry, unobserved illumination or
correlated modelling errors.

## Replacement figure: the experiment and its demonstrated reach

1. **Observed history and its timing.** Use an actual target's measured history,
   identify the observed epochs and gaps, and show an explicitly illustrative
   MIRI epoch. Connect this history to the response being tested rather than
   assuming a monotonic rise/fall label describes the entire object.
2. **Predictions that existing data allow.** Show current-power equilibrium and
   delayed fixed-distribution spectral ranges fitted under the same historical
   constraints, allowing each family its nuisance parameters. Add evolving dust
   only if its additional information is demonstrated. Distinguish measured
   points from simulated MIRI data; show statistical precision and relevant
   host/calibration covariance separately. Any range must say how it was made;
   a finite grid envelope is not a posterior interval.
3. **What can be recovered across 24 histories.** Fit mock observations with
   competing families, including null injections. Report recovery and false
   positive performance, or the constraints attainable where models overlap.
   Account for missing/future illumination, host uncertainty, response scale
   and geometry. A selected pair of separated curves cannot replace this test.

The strongest defensible lead test is whether incorporating the recorded
history improves on an adequately flexible current-power equilibrium model.
The dust-evolution test is secondary and must permit an inconclusive result.
Geometry and historical illumination can imitate spectral memory; the figure
must show how the available constraints limit that freedom.

If this validation does not demonstrate discrimination, the figure should
instead quantify the attainable dust/host constraints and the science claims
must be revised accordingly. Better graphics alone cannot establish the
missing inference. No replacement detection significance is justified yet.

## Inputs checked

- `make_fig2_miri.py`
- `review/data_anchored_scenarios/summary.json`
- `review/p9694_predictivity/README.md`
- `review/p9694_predictivity/adequacy_and_forecasts.json`
- `inputs/warm_response_forecast_summary.json`
- `forecast_warm_response.py`

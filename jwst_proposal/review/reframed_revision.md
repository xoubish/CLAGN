# Reframed proposal revision — 28 September 2026

The active source/PDF is `proposal_reframed.tex` / `proposal_reframed.pdf`.
Prior main and reframed drafts are preserved under `../older/`.

## Scientific argument

The two tests are now distinct: measure history-dependent warm emission, then
compare the spectrum with delayed heating of fixed dust to test structural
evolution. The historical reference is a joint fit of 2010 W3/W4 and MRS through
their instrument bandpasses, rather than treating W3 as an observed rest-frame
8–13-micron luminosity. The monochromatic observed-frame 12-micron forecast is
explicitly a planning proxy for W3.

The immediate-adjustment boundary forecast has no fitted condensation timescale.
Finite adjustment times are an intended science-analysis extension. The text
no longer equates an instantaneous-response residual with condensation time or
dust mass. It also removes the universal months-long sublimation statement and
claims of no previous re-formation constraints. References 18–19 already discuss
year-scale boundary adjustment/replenishment.

The obscuration/excitation paragraph is shortened: silicates and ionic lines
are complementary diagnostics, not a binary emission/absorption test or a
universal few-year clock. The nominal 27.9-micron coverage includes [Ne V]
24.3 microns for 13 sample targets and [O IV] 25.9 microns for nine; coverage
is not a guaranteed line detection. The unsupported S/N>10 line claim is removed.

## Forecast scope

The saved baseline gives sign-weighted differences of 0.06203 and 0.05084 dex,
with adopted mean error 0.01377 dex (4.50 and 3.69 sigma). These are conditional
on the shell model and error assumptions; no marginalized model-discrimination
significance has been established. W1 is a host-subtracted heating proxy with
an assumed response exponent, not an independently recovered optical/UV driver.

Photometric/decomposition errors alone are 0.03068–0.04445 dex; adding the
current W1 host-fraction term gives 0.03148–0.18778 dex (median 0.04347).
Twelve objects have a maximum pairwise model separation exceeding three times
the measurement-only error; five do so with the host term. This is not a count
of detections. The error model does not include the full distribution of dust
geometries, heating conversions, future histories or correlated uncertainties.
The text states that recovery simulations will address these effects.

The saved signs are 13 negative and 11 positive. A common additive log-flux
shift contributes -1/12 of that shift to the unweighted sign-weighted mean;
it is suppressed, not exactly cancelled. Family labels remain 12 rising and
12 fading and are not the same thing as the forecast signs.

## Figures

Figure 1 retains F06/F04/R01 and the existing measurement and epoch selections,
with smaller cutouts, larger type, shorter legends and six selected expected
line identifiers. It writes its final PDF directly and saves a display inventory.

Figure 2 removes the inherited 0.20-dex-scatter precision band, which was not the
error model used by the revised forecast. Its physical response curves stay unchanged.

Figure 3 now plots actual measured total W1 points instead of labelling the
transformed heating proxy as measurements. Its example spectrum uses the same
saved heating exponent as the ensemble forecast (previously the imported module
left this unset). The hypothetical MIRI point is labelled Forecast. Generic 2%
and 5% bands are removed; the primary warm band is highlighted instead. The
ensemble retains its adopted measurement and host-fraction bars, with a shorter
heading and an explicit observed-frame W3-proxy axis label.

## Validation

Build/layout and APT consistency results are recorded in `consistency.json`.
The official style, sample and observing configuration are unchanged. Full
astrophysical forecast validation remains separate from these editorial checks.

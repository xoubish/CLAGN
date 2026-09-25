# Revision v0.4 — 25 September 2026

The rewrite retains the six-target MIRI/MRS program and makes the central test explicit: can a fixed dust distribution, responding to an observationally constrained luminosity history, explain the nuclear spectrum and infrared time series together?

## Scientific changes

- Replaced the claim that NEOWISE's end leaves no infrared monitoring with a wavelength/diagnostic argument. SPHEREx supplies short-wavelength information; MIRI adds the warm-dust continuum, both silicates, and spectral constraints on contamination. No unsupported claim that CLAGN lack all prior mid-infrared spectroscopy remains.
- Adopted A3 as the main observational example, preserving the requested 15-arcsec cutout, separate observed-wavelength optical/SPHEREx panels, common flux units and alternate-standard NGPS reduction. P2190 remains supporting material: its redshift places the 18-micron silicate peak outside MRS.
- Incorporated the transferable diagnostic strategy from ProtoSB: fit aromatic bands and ionic lines alongside nuclear dust emission. High-ionization lines are measurements or upper limits, not promised detections or instantaneous accretion tracers. No new instrument or resolved-torus claim was added.
- Replaced deterministic silicate/obscuration claims with a joint test of profiles, continuum, geometry and histories. The comparisons are not assumed perfectly static.
- Removed unvalidated nuclear fade factors, broad-line percentage changes, all-target SPHEREx coverage, achieved S/N and a unique dust-mass/formation-timescale promise.
- Specified a positive result and a useful null result, with host emission, uncertain heating history, epoch differences, calibration covariance and model diversity included.

## Figure 2

The regenerated model is explicitly illustrative, not fitted to A1. It integrates the emitted spectra over each optically thin spherical shell's light-travel response rather than treating an infrared echo as a directly measured illuminating history. Newly formed shells have no pre-fade emission inside the original dust-free cavity. Both models use the same outer density normalization and are then normalized to their own rest 13–14 micron continuum.

Assumptions: factor-six step decline in 2019.5, observation in 2027.5, A1 redshift, old inner radius 0.25 light-years, outer radius 75 light-years, emissivity index 1.5 and sublimation temperature 1500 K. The added inner component follows the same density law down to the new sublimation radius. This is a limiting response, not a grain-nucleation calculation; opacity features are schematic. Density indices 1 and 1.5 give rest 4.5–5.5/13–14 micron colour changes of 2.63% and 12.25%. Panel (b) uses index 1.5; the shaded ratio spans both. Assumptions and numerical outputs are in `inputs/illustrative_model_summary.json`.

## Required technical closure before submission

1. **ETC and APT:** replace the planning 60-group/33-minute exposures and 14.1 charged hours with exported target-specific calculations. Statistical S/N of 50 near rest 5 and 13–14 microns and 30 across the 18-micron feature are design goals, not verified achievements. Use conservative nuclear/host spectra, the current channel-4 throughput, and check saturation and acquisition. Do not substitute the schematic normalized model for a measured flux prediction.
2. **Backgrounds:** the existing allocation includes dedicated skies for A1, A2 and C3. STScI recommends dedicated backgrounds. Inspect whether A3, C1 and C2 provide genuinely source-free local sky and whether it is adequate; otherwise add matched skies and update the charged time. Brightness alone is not a reason to omit a sky. Simultaneous imaging during dedicated sky is an optional later optimization requiring field/roll/saturation/APT checks; it has not been promised.
3. **Quantitative discrimination:** obtain posterior predictive spectra for each anchor using its own optical history, propagate host and calibration uncertainty, and inject/recover fixed/evolving cases at the final ETC precision. The smallest illustrative colour changes may not be distinguishable. A one-epoch MIRI spectrum plus sparse histories need not uniquely identify restructuring; demonstrate which parameter combinations can be constrained.
4. **Ancillary data:** A3 has a verified extracted SPHEREx CSV with three intervals. Extract/validate other sample targets before promising their availability. Treat epochs independently, account for aperture and correlated errors, and update optical photometry near the JWST epoch. The bundled ZTF snapshot ends in October 2025, whereas the latest optical spectra are September 2026.
5. **Sample and planning values:** quantify nuclear transition amplitudes and comparison matching with consistent host subtraction. Validate the AllWISE error-column units before using them in likelihoods; the tabulated flux densities retain the prior draft's values. Finish coordinates, visibility/duplication checks and the other required APT/team fields using the final configuration.

## Document checks

The official Cycle 6 style file is unchanged. Explicit T1 font encoding allows Tectonic to render the template's Times bold/italic correctly without changing font size or margins. The final PDF has four core pages and one page containing supplemental information, references and disclosure. All core pages were rendered and inspected; the Makefile build succeeds, with no overfull boxes or font-substitution warnings. One reference paragraph has a harmless underfull-box warning. The source, APT abstract, active figures, build instructions and portable bundle recipe have been updated together. The previous v0.3 proposal and model remain in `previous_v0.3/`.

Relevant official guidance: [Cycle 6 PDF requirements](https://jwst-docs.stsci.edu/jwst-opportunities-and-policies/jwst-call-for-proposals-for-cycle-6/jwst-preparing-the-proposal-pdf-attachment), [MRS strategies](https://jwst-docs.stsci.edu/jwst-mid-infrared-instrument/miri-observing-strategies/miri-mrs-recommended-strategies), [dedicated sky](https://jwst-docs.stsci.edu/jwst-mid-infrared-instrument/miri-operations/miri-dithering/miri-mrs-dedicated-sky-observations), and [calibration status](https://jwst-docs.stsci.edu/jwst-calibration-status/miri-calibration-status/miri-mrs-calibration-status).

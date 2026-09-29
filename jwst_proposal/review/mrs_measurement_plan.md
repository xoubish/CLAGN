# From MRS data to the proposal's answers

Figure 2 is organized around the analysis of the proposed data. It is not a
claim that two particular illustrative curves establish model discrimination.
The primary products are spectra, continuum luminosities and shapes, spectral
features, line fluxes or limits, and spatial profiles. Dust radii, temperatures,
covering factors and adjustment times are model inferences, not direct MRS
measurements of an unresolved torus.

| Data product / measurement | Comparison or analysis | Scientific question / result |
| --- | --- | --- |
| Wavelength-dependent unresolved nuclear spectrum, extended host spectrum, nuclear fraction and spatial profiles | Joint point-source and extended-emission fitting to the cube; surrounding spaxels and PAH emission constrain the host | Is the inferred warm response nuclear? Establish host-subtraction bounds for both main tests |
| Nuclear rest 8–13 micron continuum luminosity, fitted continuum levels and slopes near 5, 12 and 18 microns | Joint fit with historical W3/W4 through their bandpasses and with optical/WISE/SPHEREx at their actual epochs | How much has warm emission changed since 2010? How does the warm SED relate to the illumination and hot-dust history? |
| Warm continuum residual relative to a current-power equilibrium prediction | Compare across recorded direction, amplitude, event age and recoveries, allowing for geometry and current nuclear state | Strength and persistence of spectral memory; whether the path to the current state matters |
| Silicate strengths S9.7 and S18 and their full profiles (peak, width, emission/absorption) | Joint continuum and feature decomposition; radiative-transfer models constrain geometry, optical depth and composition | Which physical dust distributions explain the spectrum? Can obscuration explain changes when combined with the optical history? |
| Continuum and silicate residuals relative to delayed fixed-distribution models driven by each history | Test allowed fixed geometries and compare with evolving distributions, carrying host and calibration uncertainty through the fit | Does delayed heating suffice, or are changes in dust amount/location required? If fixed dust suffices, constrain the permitted departure rather than assert evolution |
| PAH fluxes, equivalent widths and spatial distributions; extended continuum | Joint spatial/spectral host decomposition allowing AGN processing of PAHs | Host contribution and circumnuclear star-forming emission; test its effect on continuum and silicate results |
| [Ne V]/[Ne VI] fluxes or limits and line ratios where measured; spatial extent | Excitation diagnostics combined with optical state and silicate constraints | Characterize highly ionized gas; support the intrinsic-change versus obscuration assessment, not an instantaneous nuclear-luminosity measurement |
| H2 rotational-line fluxes or limits and spatial distributions | Characterize warm molecular emission alongside PAHs and ionic lines | Host gas/excitation context for the two cohorts; no new promised feedback or outflow experiment |

S_sil = ln(F_nu,feature / F_nu,continuum); use a consistently fitted continuum
and report the profiles, rather than treating either strength as a unique
optical-depth or composition measurement. Continuum slopes refer to the fitted
continuum, not raw points that include lines and features. No torus boundary
motion is spatially resolved by this experiment. Comparing one spectrum from
each object across known histories is an ensemble response test, not a directly
measured multi-epoch MIRI lag.

The warm reference is a forward fit to the archival bandpasses. MRS does not
supply the entire redshifted W4 band for every target; extrapolation beyond
coverage is part of that model comparison. SPHEREx and WISE constrain their
own observed epochs; optical monitoring anchors the nuclear state near JWST.
The figure does not assume a contemporaneous SPHEREx observation.

## Figure design

1. Panel A: R01's measured optical and WISE histories, corrected optical alerts,
   and three synthetic-W1 SPHEREx points. The MIRI date is illustrative.
2. Panel B: measured SPHEREx and simulated MIRI spectrum, with 8–13 micron
   shading, both silicate features and PAH positions. The residual strip shows
   that both refitted dust responses remain viable under the error assumptions.
3. Panel C: warm emission to illumination history; silicate profiles to optical
   depth and geometry; PAHs/spatial profiles to bounds on host emission. Only
   continuum and silicate precision is quantified in the mock. No posterior
   reduction or model-discrimination significance is implied.

## Source checks

The diagnostics retain the proposal's existing references. In particular:

- [Sirocky et al. 2008](https://uknowledge.uky.edu/physastron_facpub/199/)
  discusses the paired silicate features, continuum placement and constraints
  on dust geometry/chemistry. It supports jointly fitting both profiles rather
  than interpreting a feature strength in isolation.
- [García-Bernete et al. 2024](https://arxiv.org/abs/2409.05686)
  measures nuclear/circumnuclear PAH properties and their modification in AGN
  environments; PAHs are constraints on the host, not a guaranteed fixed
  conversion from continuum to star formation.
- [Hermosa Muñoz et al. 2025](https://www.aanda.org/articles/aa/pdf/2025/01/aa52437-24.pdf)
  supplies the MRS ionized/molecular-gas context already cited in the proposal.
- Proposal references 13, 14, 28 and 49 support the existing joint spatial and
  spectral decomposition plan. No new instrument-performance claim is added.

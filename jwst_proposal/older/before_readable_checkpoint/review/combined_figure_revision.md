# Combined science figure and main-file rename

The active draft is now `proposal.tex` / `proposal.pdf`. The build, consistency
audit and current documentation point to those names. The immediately preceding
draft and both old science figures are saved in `older/before_combined_figure/`.
Figure 1 is byte-for-byte unchanged, retaining F06, F04 and R01.

Figure 2 now gives the two proposal questions large panels on a shared rest-frame
wavelength axis. Rising/fading colors are consistent throughout. The top strip
identifies the additional spectral range supplied by MIRI; the lower annotations
connect continuum/silicate measurements to the science and spatial/PAH data to
host separation. No instrumental uncertainty band or detection claim is added.
The old per-target forecast plot is archived; its calculations and quantitative
paragraph remain. Figure references and the caption have been updated.

The curves use the unchanged shared shell calculation in
`make_fig2_diagnostics.py`. Both panels use age 10 initial silicate Rin/c after
fourfold changes. Panel A divides delayed fixed dust by its own current-power
equilibrium. Panel B compares evolving with delayed fixed dust after fitting one
amplitude to equal integrated 2–4 micron luminosity. This does not imply identical
hot spectral shapes. The boundary responds immediately when illumination arrives;
residual signs depend on that prescription and normalization. This is an
illustration of discriminating observables, not a target fit or a population
power calculation. Exact plotted curves and provenance are saved under `inputs/`.

The compiled PDF has five core pages and seven pages total. The existing APT,
anonymity, page-limit, reference and layout audit passes. Page 4 was rendered
and visually reviewed; the official proposal style is unchanged.

## Spectrum-led revision

The preceding two-panel version is saved in `older/before_spectrum_figure/`.
The active figure now has a dominant model-spectrum panel (A), a smaller
memory panel (B), and a boundary-evolution panel (C). A/C use exactly the same
rising example and hot-band normalization. C shows the MRS-range fractional
separation of the two spectra in A; this equality was checked numerically.
The former response-curve CSV and Figure 1 are byte-for-byte unchanged.

PAH and high-ionization markers locate complementary diagnostics; they are
not synthesized features or forecasts of line detections. Their positions
are inherited from the earlier diagnostic figure and proposal. The coverage
strip shows the intersection across the actual sample redshifts. The
caption defines the physical illustration, normalization and panel baselines.

The final page 4 was rendered and visually reviewed. The current PDF remains
seven pages overall with five core pages, and the consistency audit passes.

## Clarifying the model comparison

The two-model overlay made the common spectral shape dominate the main panel.
Panel A now displays one illustrative nuclear spectrum; panel C retains the
quantitative model separation and labels the computed differences at 12 and
18 microns (+20% and +33%, rounded). No curve amplitudes, assumptions or
normalizations changed. The source spectra remain in the saved CSV. These
values describe this illustrative calculation, not a detection significance.
The caption and provenance distinguish the single spectrum from the model
comparison. The preceding version is in
`older/before_spectrum_comparison_cleanup/`. Page 4 was visually reviewed;
the five-core-page consistency audit passed.

## Measurement-to-science workflow

Figure 2 is now organized around the analysis of the proposed cubes. The full
measurement inventory and scientific mapping are in `mrs_measurement_plan.md`.
A labelled spatial schematic shows nuclear/host decomposition; a nuclear
spectrum locates continuum and spectral diagnostics. Observable cards feed
each target's historical record and two science branches: current-power
spectral memory and delayed fixed-dust versus evolving-distribution fits.
The supporting obscuration/excitation and host-gas tests are explicit.

No toy model separation is used as evidence of detectability in this figure.
The previous version is preserved in `older/before_measurement_workflow/`.
The first analysis paragraph and caption now describe the same measurements.
The history-driven forecast paragraph is unchanged. Figure 1 is byte-for-byte
unchanged. The spectrum is numerically identical to the preceding fixed-dust
illustration; the spatial profile is a dimensionless schematic, not an empirical
PSF or a simulation of a target.

Pages 3 and 4 were rendered and reviewed. The build and consistency audit pass
with five core pages and seven total pages. The official style is unchanged.

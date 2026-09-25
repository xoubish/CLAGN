# Ideas from ProtoSB for the CLAGN proposal

Reviewed 2026-09-25. Source: `/Users/shemmati/Desktop/ProtoSB.pdf`,
Cycle 5 GO proposal 12259, “Intensely Radio Quiet Galaxies: Proto-Starbursts
or Proto-AGN?” All 11 PDF pages were read; Figures 1 and 2 on PDF pages 8
and 9 were inspected. Page numbers below refer to PDF pages, not the
scientific-justification page numbering. This is an assessment, not an
implementation of changes to the science program or observing setup.

The transferable strength is the connection between each physical alternative,
the diagnostic that discriminates it, and the required observation. Our central
question should remain whether the dated luminosity history can explain the
dust spectrum with a fixed dust distribution. ProtoSB's radio deficiency and
very young starburst interpretation are not established for our targets.

## Highest-priority additions

1. **Specify what MIRI contributes beyond continuum photometry.** ProtoSB
   pages 7–8 name PAHs, molecular hydrogen, recombination lines, and
   high-ionization lines as distinct physical diagnostics. In our proposal,
   add [Ne VI] 7.652 and [Ne V] 14.32 micron measurements or upper limits as
   constraints on the hard ionizing radiation, and fit the PAH bands and
   lower-ionization lines jointly with the continuum to constrain the host
   contribution. These support the dust-history test by reducing plausible
   contamination explanations. Do not promise detections without line-flux
   sensitivity calculations, or treat narrow-line emission as an instantaneous
   luminosity measurement. PAHs and [Ne II] are not perfectly pure star-formation
   tracers in AGN. A relevant demonstration of joint PAH, nuclear-continuum,
   ionization and H2 analysis is [Donnan et al.](https://arxiv.org/abs/2210.04647).

2. **Show a predicted MIRI spectrum with the decisive measurements marked.**
   ProtoSB Figure 2 (page 9) makes the gain in spectral information visible.
   Our new Figure 1 already establishes the target histories and SPHEREx
   coverage; the complementary prediction figure should show the competing
   dust models, the MIRI wavelength range, both silicate features, and the
   continuum intervals used in the model comparison. Use actual SPHEREx
   measurements and WISE photometry where available, and clearly distinguish
   observations from predictions. If the model remains for A1, label it as
   A1 rather than connecting its numerical forecast to the A3 example.

3. **Set exposure requirements from the discriminating measurement.** ProtoSB
   page 9 derives MRS time from the faint continuum in the silicate trough.
   For us, derive a target-specific ETC requirement for continuum ratios and
   silicate-profile measurements that separate the models, with nuclear/host
   components and wavelength-dependent correlated errors included. A continuum
   colour cancels a common multiplicative scale error, not all spectrophotometric
   errors. Our current ETC/APT numbers are explicitly provisional; this review
   does not validate them. Use the SPHEREx spectra to constrain the short-wave
   model, not as measurements of the unobserved 10–25 micron continuum.

## Useful observational option

ProtoSB pages 4 and 8 put the target on the MIRI imager while MRS obtains a
dedicated sky exposure. Our draft already includes backgrounds for A1, A2 and
C3, making this an option worth evaluating for those visits. It could supply
contemporaneous imaging of the nucleus and host. The imager and MRS have
different fields: ordinary simultaneous imaging during MRS science does not
image the same nucleus. Roll angle, empty MRS sky, imaging background,
saturation and APT timing must be checked. This is an opportunity, not a claim
that imaging has zero overhead or that backgrounds should be added to every
target just to obtain it. See [STScI dedicated-sky guidance](https://jwst-docs.stsci.edu/jwst-mid-infrared-instrument/miri-operations/miri-dithering/miri-mrs-dedicated-sky-observations)
and [simultaneous-imaging guidance](https://jwst-docs.stsci.edu/jwst-mid-infrared-instrument/miri-operations/miri-mrs-simultaneous-imaging).
Choose filters at our targets' redshifted features rather than copying ProtoSB's
nearly zero-redshift filter choices.

## Coverage for the two illustrated objects

Observed microns, calculated as rest wavelength times (1+z), using
z=0.2377224 for A3 and z=0.739785 for P2190. Coverage assumes MRS 4.9–27.9 microns;
being in range does not establish detectability.

| Diagnostic | Rest µm | A3 observed µm | P2190 observed µm |
|---|---:|---:|---:|
| PAH | 6.2 | 7.67 | 10.79 |
| [Ne VI] | 7.652 | 9.47 | 13.31 |
| Silicate | 9.7 | 12.01 | 16.88 |
| PAH | 11.3 | 13.99 | 19.66 |
| [Ne II] | 12.81 | 15.86 | 22.29 |
| [Ne V] | 14.32 | 17.72 | 24.91 |
| [Ne III] | 15.56 | 19.26 | 27.07 |
| H2 S(1) | 17.035 | 21.08 | 29.64 — outside |
| Silicate | 18 | 22.28 | 31.32 — outside |
| [O IV] | 25.89 | 32.04 — outside | 45.04 — outside |

This reinforces A3's suitability for the existing two-silicate-feature science
case. In particular, do not insert [O IV] 25.89 microns as a promised diagnostic
for our z=0.11–0.28 sample: it lies beyond MRS throughout that redshift range.

## Interpretive limits to keep explicit in the science argument

- ProtoSB's dust-geometry emphasis is useful, but our statement that silicates
  “settle the obscuration alternative” is too categorical. Fit both profiles
  and their intervening continuum jointly with the time-domain evidence,
  host contribution and dust geometry. A single-epoch silicate spectrum does
  not identify the cause of the optical transition on its own. See
  [Sirocky et al.](https://arxiv.org/abs/0801.4776) and
  [Nenkova et al.](https://arxiv.org/abs/0806.0512).
- ProtoSB's approximately 10-pc spatial goals do not transfer to our distances.
  For example, with Astropy Planck18, 0.3 arcsec corresponds to about 1.17 kpc
  at A3 and 2.25 kpc at P2190. The torus is unresolved; use PSF-based nuclear
  extraction and constrain extended host emission where measurable.
- H2 excitation/kinematics can be supporting science, but a line detection
  alone is not proof of shocks, inflow or outflow. Keep this secondary to the
  dust-history experiment.
- Additional NIRSpec/NIRCam observations need a specific unmet requirement;
  adopting all three instruments would add scope and overhead without yet
  establishing a benefit for our core test.

## Candidate paragraph for a later revision

“MIRI will connect the measured near-infrared histories to the warm-dust
continuum and both silicate features. Joint modelling of PAH emission and
ionic lines will constrain host-galaxy contamination and the ionizing source,
while the silicate profiles constrain dust opacity and geometry. These
measurements test whether a fixed dust distribution illuminated by the
observed luminosity history can reproduce the nuclear infrared spectrum,
or whether an evolving dust distribution is required.”

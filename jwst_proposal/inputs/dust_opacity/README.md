# Optical properties used in Figure 2

Public tabulations by B. T. Draine, from Draine & Lee (1984) and
Laor & Draine (1993, ApJ, 402, 441):

- https://www.astro.princeton.edu/~draine/dust/dust.diel.html
- https://www.astro.princeton.edu/~draine/dust/diel/Sil_21.gz
- https://www.astro.princeton.edu/~draine/dust/diel/Gra_21.gz

Downloaded unchanged on 2026-09-26. Figure 2 selects the 0.1-micron radius
block and interpolates positive absorption efficiencies in log wavelength
and log efficiency. File checksums are recorded in
`../fig2_diagnostics_provenance.json`.

Run `MPLCONFIGDIR=/tmp/clagn_mpl /opt/anaconda3/bin/python
jwst_proposal/make_fig2_diagnostics.py` from the workspace root to reproduce
the figure, numerical checks, spectrum CSV and provenance JSON.

The calculation integrates separate silicate/graphite radiative-equilibrium
temperatures through optically thin, isotropic shells. It uses idealized
factor-four steps and an instantaneous luminosity-dependent inner boundary,
including the full angular light-travel delay of each shell. Both hypotheses
are scaled to the same integrated rest 2--4 micron luminosity, representing
a fitted dust-mass amplitude. This does not establish agreement with actual
SPHEREx colors or optical/WISE histories.

These spectra are illustrations, not clumpy-torus radiative-transfer fits,
measurements of any target, or a forecast of population discrimination.
The current figure shows schematic vector cartoons of the boundary motion
above the computed fractional spectral differences. The earlier spectrum
overlays and ETC precision strip were removed following user feedback.
No sensitivity has been calculated for these illustrative spectral shapes;
the existing ETC calculations remain in the observing description.

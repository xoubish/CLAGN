# R06 aperture comparison: 0.6 arcsec

Archive: `../downloads/web_exports/wb325241_c7_2026-09-26_19.02.44.tar`.
The download completed before analysis. `analysis.json` records its SHA-256.
Compared recursively with `../r06_ch3a_27groups_r1/input.json`, the only input
change is strategy.aperture_size: 1.0 to 0.6 arcsec. Source spectrum,
normalization, numerical sky background and detector configuration match.

The exported flux/noise arrays reproduce native S/N=7.8311 at 12.81 microns.
The user's scalar report gives 299.70 s on source and 599.41 s including the
off-source sky strategy, without APT overheads. The reported maximum
saturation fraction remains 0.00104. The scalar report is not in the archive.

| Observed bin center (um) | R~100 S/N, radius 1.0 arcsec | R~100 S/N, radius 0.6 arcsec |
| --- | --- | --- |
| 11.70 | 37.04 | 54.91 |
| 12.00 | 32.81 | 48.49 |
| 12.50 | 37.79 | 55.36 |
| 12.81 | 41.10 | 59.91 |
| 13.00 | 39.13 | 56.95 |
| 13.30 | 35.05 | 50.93 |

`r100_full_subband.csv` also covers this sub-band using 15 complete contiguous
nominal R=100 bins anchored at 12.81 um. Their S/N range is 49.51–59.91.
The six comparison bins above use specific centers and are a separate check;
the bin at 12.00 um gives 48.49. Values depend slightly on bin placement and
rounding to whole output spectral pixels. Incomplete edge bins are omitted.
Do not claim that every bin meets S/N=50.

These are diagonal-covariance estimates, calculated using ../bin_archive.py,
not a Pandeia rerun. The spectra contain no full spectral covariance and the
scene omits host photon noise. Calibration, cube reconstruction, aperture
correction and host-decomposition errors are not established. This is still
the assumed point-source power-law scenario, not a fitted nuclear spectrum.
The smaller aperture improves statistical sensitivity in this scenario but
has not been established as optimal or as the final extraction choice.

Next inspect the channel/setting controls for the rest-12-micron anchor,
observed 14.5332 microns at sample z=0.2111. The current 11.53–13.48 um
sub-band does not reach that anchor. Continue wavelength checks before
setting the final exposure time. No scientific-proposal quantities updated.

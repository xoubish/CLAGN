# R06 channel 4B: continuum near rest 18 microns

Archive: `../downloads/web_exports/wb325241_c7_2026-09-26_19.14.26.tar`.
The checksum, copied input and bin statistics are recorded in this directory.
Comparison with the channel 3B run finds exactly three input changes: channel
3 -> 4, reference wavelength 14.53 -> 21.80 um, and aperture radius 0.6 ->
0.9 arcsec. Source, sky spectrum and detector settings are unchanged.

The exported arrays reproduce native S/N=2.8010 at 21.80 um. The output covers
20.5723–24.3967 um. The source is still an assumed point-source power law with
WISE W3 normalization 0.95 mJy; it has no silicate feature or extended host.
The 0.9-arcsec radius is a trial extraction, not an established optimum.

| Observed bin center (um) | Approximate R=100 S/N |
| --- | --- |
| 21.00 | 21.60 |
| 21.80 | 20.29 |
| 22.50 | 17.22 |
| 23.00 | 14.85 |
| 24.00 | 11.65 |

Sixteen complete contiguous bins anchored at 21.80 um give S/N=11.48–23.25.
This calculation assumes independent spectral-pixel errors and sums the
exported extracted signal and noise in quadrature. Host photon noise,
host-decomposition errors, spectral covariance and additional calibration/
cube-reconstruction uncertainties are not included. These are continuum
sensitivity estimates, not silicate-feature detection significances.

The user's scalar report gives 299.70 s on source, 599.41 s including sky,
maximum saturation fraction 0.00100, and background 313.51 MJy/sr at the
reference wavelength. APT charged time is not established by these values.

At 21.80 um the estimate of 20.29 is below the provisional goal of 30.
Pure square-root-of-time scaling suggests an exposure multiplier
(30/20.289862)^2 = 2.186, but changing the ramp length changes read-noise
performance and must be recalculated. Next trial: change only groups per
integration from 27 to 60, retaining 1 integration, 4 dithers, FULL/FASTR1,
the sky strategy and 0.9-arcsec extraction radius. This is a trial, not a final
exposure or a claim that the whole 4B sub-band will meet S/N=30.

Longer-wavelength bins are weaker; even if the 21.80-um anchor reaches 30,
the red continuum wing must be assessed separately against the science need.
Channels 3B and 4B are simultaneous: do not add their times. A longer B-setting
exposure should later be checked for sensitivity and saturation in all four
channels. No full-sample or Medium-program timing claim is established.

Reproduce from the workspace root:

```sh
python jwst_proposal/etc/bin_archive.py jwst_proposal/etc/downloads/web_exports/wb325241_c7_2026-09-26_19.14.26.tar jwst_proposal/etc/r06_ch4b_27groups_r09 --centers 21 21.8 22.5 23 24
```

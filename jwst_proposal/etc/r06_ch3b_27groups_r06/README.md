# R06 channel 3B: rest-12-micron continuum anchor

Archive: `../downloads/web_exports/wb325241_c7_2026-09-26_19.08.41.tar`.
The checksum, copied input and calculated bin statistics are saved here.
Recursive comparison with the 0.6-arcsec channel 3A run confirms only two
changes: disperser short -> medium and reference wavelength 12.81 -> 14.53 um.
The source, WISE W3 normalization, sky spectrum, aperture and detector setup
are unchanged. The point-source power law remains an assumed sensitivity
scenario, with no extended host in the scene.

The saved arrays reproduce native S/N=9.1160 at 14.53 um. Output wavelengths
span 13.3712–15.6298 um. This includes the rest-12-um anchor at 14.5332 um for
sample z=0.2111 (the rounded saved source redshift remains 0.211).

| Observed bin center (um) | Approximate R=100 S/N |
| --- | --- |
| 13.50 | 54.27 |
| 14.00 | 63.45 |
| 14.53 | 69.35 |
| 15.00 | 70.98 |
| 15.50 | 69.47 |

Fourteen complete contiguous nominal R=100 bins anchored at 14.53 um give
S/N=56.46–71.08. The separately centered comparison bin at 13.50 um gives
54.27. Incomplete band-edge bins are omitted. In this assumed source model,
these bins exceed the provisional statistical S/N goal of 50. No claim is
made for covariance, calibration, aperture correction or host-decomposition
accuracy; the exported one-dimensional noise does not supply full spectral
covariance and host photon noise is absent.

From the user's copied scalar report: on-source exposure 299.70 s; strategy
time including sky 599.41 s; saturation fraction 0.00121. These scalar values
are not included in the archive. Charged time still requires APT overheads.

Reproduce from the workspace root:

```sh
python jwst_proposal/etc/bin_archive.py jwst_proposal/etc/downloads/web_exports/wb325241_c7_2026-09-26_19.08.41.tar jwst_proposal/etc/r06_ch3b_27groups_r06 --centers 13.5 14 14.53 15 15.5
```

The script now accepts explicit --centers; without them it checks the archive's
reference wavelength. Earlier 3A comparison points can be reproduced using
--centers 11.7 12 12.5 12.81 13 13.3.

Next trial: channel 4, Medium (B), reference wavelength 21.80 um (rest 18 um),
same detector settings. Use a trial radius of 0.9 arcsec, approximately the
3B radius scaled by 21.80/14.53. This is not an optimized aperture or a
validated encircled-energy match. The power law has no silicate feature;
this trial measures continuum sensitivity near the feature, not feature
detection significance. Channels 3B and 4B are simultaneous and their exposure
times must not be added as separate spectral-setting observations.

Official wavelength coverage and simultaneous-channel documentation:
https://jwst-docs.stsci.edu/jwst-mid-infrared-instrument/miri-instrumentation/miri-filters-and-dispersers

No proposal sensitivity, total time, or program-size claim was updated.

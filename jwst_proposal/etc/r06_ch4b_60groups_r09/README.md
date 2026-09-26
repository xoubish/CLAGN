# R06 channel 4B: 60-group comparison

Archive: `../downloads/web_exports/wb325241_c7_2026-09-26_19.17.43.tar`.
Its checksum, copied input and approximate R=100 estimates are saved here.
A recursive input comparison with the preceding 27-group channel 4B run
confirms that only configuration.detector.ngroup changed, from 27 to 60.

The saved arrays reproduce native S/N=6.8049 at 21.80 um. Output coverage
remains 20.5723–24.3967 um.

| Observed bin center (um) | R~100 S/N, 27 groups | R~100 S/N, 60 groups |
| --- | --- | --- |
| 21.00 | 21.60 | 53.17 |
| 21.80 | 20.29 | 49.25 |
| 22.50 | 17.22 | 42.28 |
| 23.00 | 14.85 | 36.89 |
| 24.00 | 11.65 | 29.35 |

Sixteen complete contiguous bins anchored at 21.80 um give S/N=28.95–56.95.
The weakest complete bin is centered at 24.09275 um. This result clears the
provisional continuum S/N=30 goal at the rest-18-um anchor and through the
sampled 23-um wing, but not throughout the whole sub-band. The 24-um region
is marginal. Retain 60 groups as a working B-setting exposure for this assumed
source, not a final requirement or a guarantee of S/N=30 at all wavelengths.

All binned values assume independent spectral-pixel errors. The point-source
power law normalized to 0.95 mJy in WISE W3 is an assumed nuclear-flux scenario.
It includes neither an extended host nor a silicate feature. The estimates
do not establish feature significance, host photon noise, calibration,
aperture-correction accuracy, cube-reconstruction covariance or decomposition
uncertainty. The extraction radius remains an unoptimized 0.9 arcsec.

The user's scalar report gives 666.01 s on source, 1332.02 s including sky,
and maximum saturation fraction 0.00223. These scalar values are not included
in the archive. Charged time requires APT overheads. Channels 3B and 4B are
simultaneous and their times must not be added. The longer B-setting exposure
still needs sensitivity/saturation checks in the other simultaneous channels.

Next test the upper boundary of the primary rest 8–13 um band: observed
15.7443 um at z=0.2111, using channel 3, Long (C), reference 15.74 um. Start
with a 0.6-arcsec trial radius and retain 60 groups, 1 integration and 4 dithers
for the initial C-setting calculation. This does not finalize C-setting time
or cover every remaining continuum diagnostic.

Reproduce from the workspace root:

```sh
python jwst_proposal/etc/bin_archive.py jwst_proposal/etc/downloads/web_exports/wb325241_c7_2026-09-26_19.17.43.tar jwst_proposal/etc/r06_ch4b_60groups_r09 --centers 21 21.8 22.5 23 24
```

No proposal quantities, total hours or program-size claims were updated.

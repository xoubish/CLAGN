# First automated continuation: R06 channel 3C

Pandeia 2026.7 local runs, using the installation validated against the web
ETC. This checks observed 15.74 um, near rest 13 um at sample z=0.2111.

Starting from `r06_ch4b_60groups_r09/input.json`, change the channel to ch3,
disperser to long, reference wavelength to 15.74 um and extraction radius to
0.6 arcsec. Preserve the same sky spectrum and assumed point-source power law
normalized to 0.95 mJy in WISE W3. The source redshift remains the inherited
0.211; for this renormalized featureless power law it does not affect shape.
The second trial changes only ngroup from 60 to 27.

| Groups | On-source time (s) | Time including sky (s) | Native S/N at reference | Approximate R=100 S/N at 15.74 um | Maximum saturation fraction |
| --- | --- | --- | --- | --- | --- |
| 60 | 666.0096 | 1332.0192 | 29.1844 | 215.86 | 0.00382589 |
| 27 | 299.70432 | 599.40864 | 12.3952 | 91.61 | 0.00172165 |

Neither calculation reports warnings. Fifteen complete R~100 bins in channel
3C have S/N=194.97–215.86 for 60 groups and 81.17–91.61 for 27 groups.
The 27-group trial clears the provisional statistical goal for this part of
the assumed continuum. It is a starting value for C, not a final setting time:
other simultaneous channels, relevant long-wavelength diagnostics and source
uncertainties still need assessment.

The corresponding `r06_ch3c_60groups_r06/` and `r06_ch3c_27groups_r06/`
directories contain input, scalar report, bin estimates and provenance.
Regenerable FITS outputs and local archives are under `local_runs/` with
matching directory names. Both use run_local.py and bin_archive.py with
--centers 15.6 15.74 16 17 17.8.

Binned estimates assume diagonal spectral covariance. The unresolved power-law
scene omits extended-host photon noise and does not establish calibration,
aperture correction, decomposition error or silicate-feature significance.
No APT overheads, full-program time or proposal quantities were updated.

Next consolidate the first target's three-setting feasibility calculation,
including the rest-5-um continuum anchor, the blue end of rest 8–13 um, and
sensitivity/saturation in the other simultaneous channels. Reassess source
SED and nuclear-flux scenarios before treating these as exposure requirements.

"""Evaluate H-beta S/N=10 exposure scenarios without changing target selection.

All outputs stay in the ignored research directory because normalization can
use internal spectra. This is continuum S/N per binned spectral pixel, not a
measurement of integrated broad-line significance.
"""
from pathlib import Path
import importlib
import json
import re
import hashlib

import numpy as np
import pandas as pd
import astropy.units as u
from astropy.time import Time

ROOT = Path(__file__).resolve().parent
OUT = ROOT / 'data/reselection_2026-09-20'
DEST = OUT / 'hbeta_snr10'
MODEL = importlib.import_module('22_september_etc')


def main():
    DEST.mkdir(exist_ok=True)
    page = (OUT / 'candidate_review_local.html').read_text()
    data = json.loads(re.search(r'<script id="candidate-data" type="application/json">(.*?)</script>', page, re.S).group(1))
    targets = [t for t in data['targets'] if any(w['night'] == 'sep23' for w in t['nights'])]
    for channel, noise, scale in zip(MODEL.CFG.channels, [2.8, 7.8, 3.7, 4.6], [.193, .189, .186, .186]):
        MODEL.CFG.readnoise[channel] = noise * u.count / u.pix
        MODEL.CFG.platescale[channel] = scale * u.arcsec / u.pix
    # These bracket sky sensitivity; neither is claimed as the predicted lunar sky.
    scenarios = [
        ('sky19_X1p3', 19., 1.3, 1.3, 0.),
        ('sky18p5_X1p3', 18.5, 1.3, 1.3, 0.),
        ('sky18_X1p3', 18., 1.3, 1.3, 0.),
        ('sky18_X1p5', 18., 1.5, 1.3, 0.),
        ('sky18_X1p5_fainter0p5', 18., 1.5, 1.3, .5),
        ('sky18_X1p5_seeing1p8', 18., 1.5, 1.8, 0.),
    ]
    rows = []
    channel_checks = []
    for target in targets:
        window = next(w for w in target['nights'] if w['night'] == 'sep23')
        base = dict(name=target['name'], pool=target['pool_role'], field_status=target['field_status'],
                    z=target['z'], historical_r=target['rmag'], longest_window_min=window['longest_minutes'],
                    windows_pdt='; '.join(r['start']+' → '+r['end'] for r in window['ranges']))
        reference = MODEL.reference(target)
        if reference is None:
            rows.append(base | dict(status='No accepted continuum reference'))
            continue
        mjd, mag, private, continuum = reference
        wavelength = 4862.7 * (1 + target['z'])  # vacuum wavelength for SDSS spectra
        low, high = wavelength/10 - 4, wavelength/10 + 4
        channels = [ch for ch in ['R', 'I', 'G', 'U']
                    if MODEL.CFG.channelRange[ch][0].to_value(u.nm) <= low
                    and MODEL.CFG.channelRange[ch][1].to_value(u.nm) >= high]
        if not channels:
            rows.append(base | dict(status='No full H-beta window in one channel'))
            continue
        options = [dict(name=target['name'], channel=ch,
                        **MODEL.calculate(ch, low, high, mag, 10., 1.3, 18.5, airmass=1.3, noslicer=True))
                   for ch in channels]
        channel = min(options, key=lambda option: option['seconds'])['channel']
        channel_checks.extend(options)
        base |= dict(reference_mjd=mjd, reference_date=Time(mjd, format='mjd').strftime('%Y-%m-%d'),
                     reference_private=private, continuum_flux_1e17=continuum, continuum_AB=mag,
                     hbeta_observed_A=wavelength, low_nm=low, high_nm=high, channel=channel)
        for label, sky, airmass, seeing, extra in scenarios:
            calculation = MODEL.calculate(channel, low, high, mag+extra, 10., seeing, sky, airmass=airmass, noslicer=True)
            # Forward-evaluate the returned sub-exposure to check the solver and
            # the read-noise-aware equal-exposure combination.
            cmd = [channel, str(low), str(high), 'EXPTIME', str(calculation['each_seconds']),
                   '-slit', 'SET', '1.0', '-binspect', '2', '-binspat', '2',
                   '-seeing', str(seeing), '500', '-airmass', str(airmass),
                   '-skymag', str(sky), '-mag', str(mag+extra), '-magsystem', 'AB', '-magfilter', 'match', '-noslicer']
            args = MODEL.ETC.parser.parse_args(cmd)
            MODEL.ETC.check_inputs_add_units(args)
            actual = MODEL.ETC.main(args, quiet=True)['SNR'].value * np.sqrt(calculation['nexp'])
            assert abs(actual-10) < .15, (target['name'], label, actual)
            assert calculation['each_seconds'] <= 900.01
            result = base | dict(status='Scenario calculation', scenario=label, goal_snr=10.,
                                  achieved_snr=actual, sky_V=sky, airmass=airmass,
                                  seeing_zenith_500nm=seeing, extra_mag=extra) | calculation
            result['integration_min'] = result['seconds']/60
            # A placeholder budget, not a measured NGPS 2x2 readout time.
            result['illustrative_slot_min'] = result['integration_min'] + 5 + result['nexp']
            result['slot_assumption'] = '5 min acquisition/slew + 1 min per read; excludes standards/calibrations'
            rows.append(result)
        print(target['name'], 'done', flush=True)
    table = pd.DataFrame(rows)
    pd.DataFrame(channel_checks).to_csv(DEST / 'channel_comparison_sky18p5.csv', index=False)
    table.to_csv(DEST / 'all_september_scenarios.csv', index=False)
    valid = table[table.status.eq('Scenario calculation')]
    summary = valid.pivot(index='name', columns='scenario', values='integration_min')
    info = valid.drop_duplicates('name').set_index('name')
    summary = info[['pool','field_status','historical_r','z','continuum_AB','reference_date',
                    'reference_private','hbeta_observed_A','channel','longest_window_min','windows_pdt']].join(summary)
    summary.to_csv(DEST / 'integration_minutes_by_target.csv')
    screened = summary[summary.field_status.eq('clear')]
    screened.to_csv(DEST / 'screened_integration_minutes.csv')
    assumptions = dict(
        snapshot=data['generated'], target_count=len(targets), computed_targets=len(summary),
        configuration='2 spatial × 2 spectral on-chip binning, 1.0 arcsec per slice; central slice only (No Slicer checked)',
        goal='Mean continuum S/N=10 per binned spectral pixel over ±40 observed Angstrom about vacuum H-beta',
        continuum='Latest accepted loaded spectrum; interpolate median rest-frame 4750–4790 and 5100–5140 Angstrom continuum windows to H-beta. Flat f_nu within the ETC calculation window. No broad-line flux added.',
        brightness='Archival spectrophotometric normalization, no unverified extrapolation or cross-band ZTF correction. The +0.5 mag case shows sensitivity to dimming/aperture loss, not a forecast.',
        source_geometry='Point-source optimal extraction. Host-dominated aperture flux can overestimate nuclear flux admitted to the extraction.',
        sky='Official ETC Gemini template uniformly scaled to Johnson V Vega mag/arcsec^2. Values 19, 18.5 and 18 are sensitivity scenarios, not computed from Moon angle. OH and scattered lunar continuum need not scale together.',
        seeing='Zenith FWHM at 500 nm; ETC applies its wavelength/airmass scaling.',
        channel_choice='Compare every channel covering the full 80-Angstrom window; choose the fastest at sky V=18.5 and X=1.3, then keep that channel fixed for all sensitivity cases. All channels expose simultaneously.',
        runtime_overrides=dict(readnoise_e=dict(zip(MODEL.CFG.channels,[2.8,7.8,3.7,4.6])),
                               platescale_arcsec=dict(zip(MODEL.CFG.channels,[.193,.189,.186,.186]))),
        slicer_guidance='https://caltechopticalobservatories.github.io/NGPS/users-manual/exposure-time-calculator.html: central slice recommended while current three-slice model is inaccurate.',
        retained_defaults='Repository throughput, dark-current assumptions, LSF, dispersion and slicer model retained; ETC source files unchanged.',
        instrument_source='https://caltechopticalobservatories.github.io/NGPS/technical-specifications.html',
        repo_sha=json.loads((OUT/'ngps_etc/github_tree.json').read_text())['sha'],
        source_hashes={name:hashlib.sha256((OUT/'ngps_etc'/name).read_bytes()).hexdigest()
                       for name in ['ETC_main.py','ETC_config.py','ETC_import.py','ETC_arguments.py']},
        split='Equal exposures <=900 s; solve separately at 10/sqrt(N) to include each read. Forward S/N checked for every scenario.',
        slot='Integration plus an illustrative 5-minute acquisition/slew allowance and 1 minute per read. These are budget assumptions; actual 2x2 overhead has not been verified.',
        scientific_limit='Does not provide integrated broad-H-beta S/N or a host/narrow-line-subtracted broad-line upper limit.',
    )
    (DEST / 'assumptions.json').write_text(json.dumps(assumptions, indent=2))
    columns = ['historical_r','hbeta_observed_A','continuum_AB','sky19_X1p3','sky18p5_X1p3','sky18_X1p3','sky18_X1p5']
    early = summary.loc[summary.index.intersection(['P1736','P3584','P3813','P7457','P7281','P7837'])]
    report = f'''# H-beta exposure calculation: S/N 10

Selection, webpage and visibility thresholds were not changed by this calculation.

**Corrected to central-slice-only extraction**, following the NGPS ETC manual recommendation to check No Slicer because the current three-slice hardware model is inaccurate. Earlier results are preserved in `superseded_three_slice/` and should not be used for scheduling. The ETC assumes optimal PSF extraction; Quicklook uses box extraction, so its measured S/N need not match the model exactly.

Source: https://caltechopticalobservatories.github.io/NGPS/users-manual/exposure-time-calculator.html

The local official NGPS ETC was run for **{len(summary)} September targets** using 2×2 spatial/spectral binning and a **1 arcsec slice width**. The goal is **continuum S/N 10 per binned spectral pixel near redshifted H-beta**, averaged over an 80-Angstrom observed-frame interval. This is not integrated broad-line significance.

Normalization uses the latest accepted loaded spectrum, including internal data where available. Brightness is archival. Sky V=19, 18.5 and 18 mag/arcsec² are explicit sensitivity cases, not a prediction for a 40-degree lunar separation. Seeing is 1.3 arcsec at zenith and 500 nm for the table below; the model scales it with airmass and wavelength. Airmass is 1.3 except in the final column.

## Early and bridging candidates: integration minutes

{early[columns].round(2).to_markdown()}

## All fields without catalogue flags: integration minutes

{screened[columns].round(2).to_markdown()}

## Interpreting a required visibility window

The visibility requirement should contain the **integration + acquisition/readout overhead + scheduling margin**, all within acceptable airmass and Moon separation. A two-hour window is not required for a 15–20-minute integration. Conversely, short geometric visibility alone does not guarantee that the chosen S/N can be reached. No new target cut or schedule is adopted here.

The CSV also includes brighter-sky calculations at airmass 1.5, a source 0.5 mag fainter, and 1.8 arcsec zenith seeing. Exposures longer than 900 seconds were split with repeated read noise included. Every solution was forward-checked to return combined S/N within 0.15 of 10. The optional slot column assumes five minutes for acquisition/slew and one minute per read; these are illustrative allowances, not measured NGPS overheads.

The native ETC sky spectrum is uniformly rescaled and does not predict the changing mixture of lunar continuum and airglow. Host/aperture differences and source variability introduce additional uncertainty. Confirming or excluding a weak broad component ultimately needs a line/continuum/host model. A continuum S/N threshold is only an exposure-planning baseline.

This report includes internal spectral normalizations and stays local. Full reproducibility settings are in `assumptions.json`.
'''
    (DEST / 'H_BETA_SNR10.md').write_text(report)
    print('\nEARLY TARGETS\n'+early[columns].round(2).to_string(), flush=True)
    print('\nSCREENED MEDIANS\n'+screened[['sky19_X1p3','sky18p5_X1p3','sky18_X1p3','sky18_X1p5']].median().round(2).to_string(), flush=True)


if __name__ == '__main__':
    main()

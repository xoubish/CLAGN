"""Build requested-galaxy archive overlays from cached public query responses.

Raw responses and downloaded files live beside target_payload.json. This separate
overlay survives regeneration of the geometry/brightness payload and changes no
observing slots or slit-aperture brightness estimates. Run before script 51.
"""
import csv
import importlib.util
import json
from pathlib import Path
from urllib.parse import quote

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.time import Time
import astropy.units as u
from bs4 import BeautifulSoup

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / 'data/additional_targets'
NAMES = ['HNGC-1377', 'HIC-1953', 'HUGC-03097', 'HIRAS-F23290-1252', 'HNGC-0814']


def read(path):
    return json.loads(path.read_text())


def ztf(target, directory):
    center = SkyCoord(target['ra']*u.deg, target['dec']*u.deg)
    candidates = {}
    for oid, obj in read(directory/'ztf_snad.json').items():
        c = obj['meta']['coord']
        sep = center.separation(SkyCoord(c['ra']*u.deg, c['dec']*u.deg)).arcsec
        band = obj['meta']['filter'].removeprefix('z')
        points = sorted(set((p['mjd'], p['mag'], p['magerr']) for p in obj['lc']
                            if np.isfinite([p['mjd'], p['mag'], p['magerr']]).all() and p['magerr'] > 0))
        if sep <= 2 and points and band in ('g', 'r', 'i'):
            candidates.setdefault(band, []).append((len(points), oid, sep, points))
    series, matches, native = {}, {}, []
    for band, options in candidates.items():
        # Do not concatenate different fields/PSF associations with unknown offsets.
        count, oid, sep, points = max(options)
        arr = np.array(points)
        groups = np.floor(arr[:, 0]/7)
        series[band] = []
        for group in np.unique(groups):
            p = arr[groups == group]
            series[band].append([round(float(np.median(p[:, 0])), 5),
                                 round(float(np.median(p[:, 1])), 5),
                                 round(float(max(np.median(p[:, 2])/np.sqrt(len(p)),
                                                 1.4826*np.median(abs(p[:, 1]-np.median(p[:, 1])))/np.sqrt(len(p)))), 5)])
        matches[band] = dict(oid=oid, separation_arcsec=float(sep), count=count,
                             first_date=Time(arr[0, 0], format='mjd').iso[:10],
                             last_date=Time(arr[-1, 0], format='mjd').iso[:10],
                             median_mag=float(np.median(arr[:, 1])))
        native.extend(dict(band=band, oid=oid, mjd=p[0], mag=p[1], magerr=p[2]) for p in points)
    assert series.get('g') and series.get('r'), target['name']
    with (directory/'ztf_selected.csv').open('w') as f:
        writer = csv.DictWriter(f, fieldnames=['band', 'oid', 'mjd', 'mag', 'magerr'])
        writer.writeheader(); writer.writerows(native)
    return series, matches


def epoch(wave, flux, label, instrument, mjd=None, **kwargs):
    date = Time(mjd, format='mjd').iso[:10] if mjd is not None else None
    return dict(wave=[round(float(w), 2) for w in wave],
                flux=[round(float(f), 5) if f is not None and np.isfinite(f) else None for f in flux],
                label=label, instrument=instrument, date=date, mjd=mjd,
                epoch_day=int(mjd) if mjd is not None else None, coadd=False, **kwargs)


def spectra(name, directory):
    epochs, dates = [], []
    if name == 'HNGC-0814':
        spec = importlib.util.spec_from_file_location('sdss_parser', ROOT/'scripts/03d_fetch_spectra.py')
        parser = importlib.util.module_from_spec(spec); spec.loader.exec_module(parser)
        data = parser.parse(directory/'daily-spec-105657-60652-63050395812350306.fits')
        assert abs(data['meta']['z']-.00539) < .001 and data['meta']['zwarning'] == 0
        epochs.append(epoch(data['wave'], data['flux'], 'SDSS-V · 2024-12-08', 'SDSS-V', 60652,
                            quality_note='Position offset 1.97 arcsec; galaxy z=0.005348 agrees with the target.', meta=data['meta']))
    if name == 'HNGC-1377':
        with fits.open(directory/'6dFJ0336391-205407:S:VR:j2009.fits.gz') as h:
            d = h['SPECTRUM VR'].data
            flux = d[0].astype(float)
            flux[~np.isfinite(d[1]) | (d[1] <= 0)] = np.nan
            epochs.append(epoch(d[3], flux, '6dFGS · 2003-11-30 · relative flux', '6dFGS', h[0].header['MJDOBS_V'],
                                flux_unit='relative', quality_note='No absolute flux calibration; 6.7 arcsec fiber.'))
        # NED provides these as ASCII even where its FITS column says N/A.
        for filename, label, relative in [
            ('NGC_1377_nuc:S:Opt:mk2006_NED.txt', 'Nuclear · Moustakas 2006 · date unverified', False),
            ('NGC_1377:S:Opt:mk2006_NED.txt', 'Integrated · Moustakas 2006 · date unverified', False),
            ('NGC_1377:S:B:ksv1995_NED.txt', 'DBSP blue · Kim 1995 · date unverified', True),
            ('NGC_1377:S:R:ksv1995_NED.txt', 'DBSP red · Kim 1995 · date unverified', True),
        ]:
            rows = [line.split('|')[0].split() for line in (directory/filename).read_text().splitlines()
                    if line.strip() and not line.startswith('#')]
            values = np.array([[float(row[0]), float(row[1])] for row in rows])
            epochs.append(epoch(values[:, 0], values[:, 1]/1e-17, label, 'NED optical',
                                flux_unit='relative' if relative else '1e-17',
                                quality_note='Publication year is not an observing date. Different apertures cannot establish variability.'))
        dates.append(dict(src='6dFGS', mjd=52969.61, date='2003-11-26', file_available=False))
    dates.extend(dict(src=e['instrument'], mjd=e['mjd'], date=e['date'], file_available=True)
                 for e in epochs if e['mjd'] is not None)
    return epochs, dates


def main():
    for name in NAMES:
        directory = BASE/name
        target = read(directory/'target_payload.json')
        series, matches = ztf(target, directory)
        epochs, dates = spectra(name, directory)
        # Require actual successful responses, not absent cache files, for negatives.
        assert read(directory/'sparcl_search.json')['status'] == 'success'
        for f in ['sdss_spectra_cone.json', 'sdss_dr20_allspec.json']:
            assert any(t.get('TableName') == 'Table1' for t in read(directory/f))
        html = (directory/'ned_spectra.html').read_text()
        assert 'Spectral data in NED archive' in html or 'There are no spectral data found' in html
        ned = read(directory/'ned_spectra_urls.json')['query']
        ned_url = ned if isinstance(ned, str) else 'https://ned.ipac.caltech.edu/cgi-bin/NEDspectra?objname='+quote(target['jname'])
        links = [dict(label='NED spectra (all wavelengths)', url=ned_url)]
        links += [dict(label=f'ZTF {band} · {m["count"]} detections', url=f'https://ztf.snad.space/dr24/view/{m["oid"]}') for band, m in matches.items()]
        if name == 'HNGC-0814':
            links.append(dict(label='SDSS-V spectrum FITS', url='https://data.sdss.org/sas/dr20/spectro/boss/redux/v6_2_1/spectra/daily/lite/105XXX/105657/60652/spec-105657-60652-63050395812350306.fits'))
        if name == 'HNGC-1377':
            soup = BeautifulSoup(html, 'html.parser')
            links += [dict(label='Optical spectrum data · '+a['href'].split('/')[-1].split(':')[0], url='https://ned.ipac.caltech.edu'+a['href'])
                      for a in soup.find_all('a', href=True) if a['href'].endswith('_NED.txt') and any(s in a['href'] for s in ['mk2006', 'ksv1995', 'j2009'])]
        status = {
            'HNGC-0814': 'SDSS-V optical spectrum loaded (2024-12-08; z=0.005348). Daily/epoch/allepoch products represent the same observing date.',
            'HNGC-1377': 'Optical spectra loaded: 6dFGS (2003-11-30), nuclear and integrated Moustakas 2006 spectra, and DBSP blue/red Kim 1995 spectra. The latter publication years are not observing dates. An additional 6dFGS date (2003-11-26) is identified but its file is not loaded.',
            'HIC-1953': 'No optical spectrum found in the searched archives. The 6dFGS target entry has no exposures; NED lists infrared and H I spectra.',
            'HUGC-03097': 'No optical spectrum found in the searched archives. NED lists an H I radio spectrum.',
            'HIRAS-F23290-1252': 'No matching optical spectrum found. A nearby DESI spectrum at z=0.0885 is inconsistent with this target at z=0.0212 and was rejected.',
        }[name]
        counts = ', '.join(f'{b}: {m["count"]} detections, offset {m["separation_arcsec"]:.2f} arcsec' for b, m in matches.items())
        summary = f'Archive search 2026-10-22: {status} ZTF DR24 light curves loaded ({counts}); plots use 7-day medians.'
        note = ('Searched SDSS DR18 (10 arcsec), SDSS DR20 allspec (5 arcsec), SPARCL including DESI (10 arcsec coordinate box), 6dFGS (10 arcsec), and NED by name. '
                'No match means none found in these searches, not proof that no spectrum exists. ZTF uses one object per band within 2 arcsec; other field IDs are not combined. '
                'PSF photometry of resolved host light is not a clean nuclear variability measurement; intrinsic variability has not been established.')
        if name in ['HNGC-1377', 'HIC-1953']:
            note += ' ZTF associations are offset 1–2 arcsec: inspect the source position before attributing changes to the nucleus.'
        overlay = dict(ztf=series, ztf_archive_matches=matches, spec=dict(epochs=epochs) if epochs else None,
                       epochs=dates, n_spec=len({d['date'] for d in dates}),
                       n_spec_available_dates=len({d['date'] for d in dates if d['file_available']}),
                       optical_archive_status=status, archive_summary=summary, archive_note=note, archive_sources=links)
        (directory/'archive_payload.json').write_text(json.dumps(overlay, allow_nan=False))
        print(name, counts, 'spectra:', len(epochs))


if __name__ == '__main__':
    main()

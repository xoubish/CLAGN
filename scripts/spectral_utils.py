"""Shared wavelength centers and archival reduction/reference selection.

Legacy six-Angstrom caches were labeled with left edges. Correct their labels
on read, without resampling flux or rewriting the downloaded cache.
"""
import json
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def bin_indices(centers, wave):
    centers = np.asarray(centers, float)
    if len(centers) < 2 or not np.all(np.diff(centers) > 0):
        raise ValueError('Need increasing bin centers')
    edges = np.r_[centers[0] - (centers[1]-centers[0])/2,
                  (centers[:-1]+centers[1:])/2,
                  centers[-1] + (centers[-1]-centers[-2])/2]
    return np.searchsorted(edges, wave, side='right') - 1


def centered_record(record):
    r = dict(record)
    w = np.asarray(r.get('wave', []), float)
    if (not r.get('grid_version') and len(w) > 2 and w[0] in (3000., 3600.)
            and np.allclose(np.diff(w), 6.)):
        r['wave'] = (w + 3.).tolist()
        r['grid_version'] = 2
        r['legacy_centers_corrected'] = True
    return r


def quality(r):
    sn = r.get('sn_median_all', r.get('meta', {}).get('sn_median_all'))
    return (r.get('metadata_quality_ok') is not False,
            r.get('archive_version') == 'master', r.get('run2d') == 'v6_2_1',
            float(sn) if sn is not None and np.isfinite(sn) else -1)


def selected_records(records):
    best = {}
    for raw in records:
        r = centered_record(raw)
        if r.get('coadd') and r.get('source', 'SDSS') != 'DESI':
            continue
        w, f = np.asarray(r.get('wave', []), float), np.asarray(r.get('flux', []), float)
        if len(w) != len(f) or np.isfinite(f).sum() < 20:
            continue
        mjd = r.get('mjd')
        day = int(np.floor(mjd)) if mjd is not None and np.isfinite(mjd) and mjd > 40000 else None
        key = (r.get('source', 'SDSS'), day)
        if r.get('source') == 'DESI' and r.get('date_verified'):
            key += (str(r.get('meta', {}).get('specid')), r.get('survey'), r.get('program'))
        if key not in best or quality(r) > quality(best[key]):
            best[key] = r
    return sorted(best.values(), key=lambda r: (r.get('mjd') is None, r.get('mjd') or 0, r.get('source', 'SDSS')))


def archival_records(name):
    records = []
    for path in [ROOT/'data/spectra_dl'/f'{name}.json',
                 ROOT/'data/reselection_2026-09-20/sdssv_spectra'/f'{name}.json']:
        if path.exists():
            records.extend(json.loads(path.read_text()))
    return records


def accepted_reference(target, records=None):
    accepted = []
    for r in selected_records(archival_records(target['name']) if records is None else records):
        mjd = r.get('mjd')
        if mjd is None or not np.isfinite(mjd) or mjd <= 40000 or r.get('coadd') or r.get('metadata_quality_ok') is False:
            continue
        warning = r.get('meta', {}).get('zwarning')
        sn = r.get('sn_median_all', r.get('meta', {}).get('sn_median_all'))
        if warning is not None and (not np.isfinite(warning) or warning != 0):
            continue
        if sn is not None and (not np.isfinite(sn) or sn < 5):
            continue
        rest = np.asarray(r['wave'], float)/(1+target['z'])
        f = np.asarray(r['flux'], float)
        a = f[(rest >= 4750) & (rest <= 4790) & np.isfinite(f)]
        b = f[(rest >= 5100) & (rest <= 5140) & np.isfinite(f)]
        if len(a) < 4 or len(b) < 4:
            continue
        continuum = float(np.interp(4862.68, [4770, 5120], [np.median(a), np.median(b)]))
        if continuum <= 0:
            continue
        fnu = continuum*1e-17*(4862.68*(1+target['z']))**2/2.99792458e18
        accepted.append((r, float(-2.5*np.log10(fnu)-48.6), continuum))
    return max(accepted, key=lambda v: (v[0]['mjd'], quality(v[0]))) if accepted else None

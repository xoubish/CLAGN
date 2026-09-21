"""14c_v2_epochs.py -- build data/spectra_epochs_v2.csv for the v2 targets from local knowledge (DR16 plate-mjd-fiber URLs, cached files,
proprietary SDSS-V epochs from 03f), so 03d/07 work even when SkyServer is down.  Usage: 14c_v2_epochs.py"""
import os, importlib.util, numpy as np, pandas as pd
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); DATA = os.path.join(HERE, 'data')
spec = importlib.util.spec_from_file_location('f13b', os.path.join(HERE, 'scripts', '13b_calibration_lines.py')); f13b = importlib.util.module_from_spec(spec); spec.loader.exec_module(f13b)
t = pd.concat([pd.read_csv(os.path.join(DATA, f'targets_{n}_v2.csv')) for n in ['sep23', 'oct26', 'oct27']]).drop_duplicates('name'); names = set(t.name)
t[['name', 'ra', 'dec']].to_csv(os.path.join(DATA, 'v2_targets_positions.csv'), index=False)
def pick_url(name, plate, mjd, fiber, survey):
    d = os.path.join(DATA, 'spectra_cache', name); cached = os.listdir(d) if os.path.isdir(d) else []
    urls = f13b.dr16_urls(plate, mjd, fiber, survey); return next((u for u in urls if os.path.basename(u) in cached), urls[0])
rows = []
pool = pd.read_csv(os.path.join(DATA, 'parent_pool_scored.csv')); pool['name'] = 'P' + pool.poolid.astype(str); pool = pool[pool.name.isin(names)]
for r in pool.itertuples():
    rows.append({'name': r.name, 'source': 'SDSS', 'sas_url': pick_url(r.name, r.plate, r.mjd, r.fiberid, r.survey), 'mjd': float(r.mjd), 'sdss_phase': 4 if r.plate >= 3500 else 2, 'programname': 'eboss' if r.plate >= 3500 else 'legacy',
                 'is_coadd': False, 'proprietary': False, 'ra': r.ra, 'dec': r.dec, 'plate_or_fps_field': r.plate, 'fiberid': r.fiberid, 'run2d': 'v5_13_2' if r.plate >= 3500 else '26', 'class': 'QSO', 'subclass': str(r.subclass), 'z': r.z, 'zwarning': 0, 'sn_median_all': r.sn_median})
sv = pd.read_csv(os.path.join(DATA, 'sdssv_internal_epochs.csv')).rename(columns={'class': 'cls'}); sv = sv[sv.name.isin(names)]
for r in sv.itertuples():
    rows.append({'name': r.name, 'source': 'SDSS', 'sas_url': r.sas_url, 'mjd': float(r.mjd), 'sdss_phase': 5, 'programname': str(r.programname), 'is_coadd': bool(r.is_coadd), 'proprietary': True, 'ra': r.ra, 'dec': r.dec,
                 'plate_or_fps_field': r.field, 'fiberid': np.nan, 'run2d': r.run2d, 'class': r.cls, 'subclass': str(r.subclass), 'z': r.z, 'zwarning': r.zwarning, 'sn_median_all': r.sn_median_all})
gal = pd.read_csv(os.path.join(DATA, 'turnon_parent_galaxies.csv')); gal = gal[gal.name.isin(names)]
for r in gal.itertuples():
    rows.append({'name': r.name, 'source': 'SDSS', 'sas_url': pick_url(r.name, r.plate, r.mjd, r.fiberid, 'sdss'), 'mjd': float(r.mjd), 'sdss_phase': 2, 'programname': 'legacy', 'is_coadd': False, 'proprietary': False, 'ra': r.ra, 'dec': r.dec,
                 'plate_or_fps_field': r.plate, 'fiberid': r.fiberid, 'run2d': '26', 'class': 'GALAXY', 'subclass': str(r.subclass), 'z': r.z, 'zwarning': 0, 'sn_median_all': r.snMedian})
E = pd.DataFrame(rows); E.to_csv(os.path.join(DATA, 'spectra_epochs_v2.csv'), index=False)
print(f'{len(E)} epochs for {E.name.nunique()} v2 targets ({int(E.proprietary.sum())} proprietary SDSS-V epochs) -> data/spectra_epochs_v2.csv')

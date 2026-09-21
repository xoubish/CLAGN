"""
13b_calibration_lines.py  --  line measurements for the SDSS-V calibration set (13_calibrate_sdssv.py): download the DR16
spectrum and the latest good SDSS-V spectrum of every calibration object and measure Hbeta / Halpha with the 03d code,
so the "dimmed" outcome can be defined by broad-line change instead of the spectrophotometric continuum.

Writes ONLY data/sdssv_calibration_lines.csv (git-ignored, proprietary-derived); never touches spectra_lines*.csv or spectra_dl/
(the scoring inputs).  FITS files are cached in data/spectra_cache/<name>/ like 03d.  Order: changed objects first, so partial
output is useful early.  Usage: /opt/anaconda3/bin/python 13b_calibration_lines.py [calibration csv] [--threads 6]
"""
import os, sys, time, argparse, importlib.util
import numpy as np, pandas as pd
from concurrent.futures import ThreadPoolExecutor
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); DATA = os.path.join(HERE, 'data')
spec = importlib.util.spec_from_file_location('f03d', os.path.join(HERE, 'scripts', '03d_fetch_spectra.py')); f03d = importlib.util.module_from_spec(spec); spec.loader.exec_module(f03d)
DR17 = 'https://data.sdss.org/sas/dr17'


def dr16_urls(plate, mjd, fiber, survey):
    """candidate DR17 spec-lite paths for the DR16 spectrum (eBOSS v5_13_2, legacy 26/103/104); 404s are skipped by 03d.fetch"""
    p, f = int(plate), int(fiber); m = int(mjd)
    eboss = f'{DR17}/eboss/spectro/redux/v5_13_2/spectra/lite/{p}/spec-{p}-{m}-{f:04d}.fits'
    legacy = [f'{DR17}/sdss/spectro/redux/{r}/spectra/lite/{p:04d}/spec-{p:04d}-{m}-{f:04d}.fits' for r in ('26', '103', '104')]
    return ([eboss] + legacy) if (str(survey).lower() in ('boss', 'eboss') or p >= 3500) else (legacy + [eboss])


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('csv', nargs='?', default=os.path.join(DATA, 'sdssv_calibration_pool.csv')); ap.add_argument('--threads', type=int, default=6)
    a = ap.parse_args()
    c = pd.read_csv(a.csv)
    changed = c.dimmed_any.astype(bool) | c.bright07.astype(bool)
    c = pd.concat([c[changed], c[~changed]])
    sv = pd.read_csv(os.path.join(DATA, 'sdssv_internal_epochs.csv'))
    sv = sv[(sv.zwarning == 0) & (sv.sn_median_all >= 2) & (sv.sep_arcsec < 1.0) & sv['class'].isin(['QSO', 'GALAXY']) & (sv.objtype == 'science')].sort_values('mjd')
    jobs = []
    for r in c.itertuples():
        ep = [dict(sas_url=u, mjd=float(r.mjd), sdss_phase=np.nan, programname='DR16', is_coadd=False, proprietary=False) for u in dr16_urls(r.plate, r.mjd, r.fiberid, r.survey)]
        g = sv[sv.name == r.name]
        if len(g):
            keep = pd.concat([g.tail(1), g.loc[[g.sn_median_all.idxmax()]]]).drop_duplicates('sas_url')       # latest good epoch + best-S/N epoch
            ep += [dict(sas_url=e.sas_url, mjd=float(e.mjd), sdss_phase=5, programname=str(e.programname), is_coadd=False, proprietary=True) for e in keep.itertuples()]
        jobs.append((r.name, pd.DataFrame(ep)))
    out = os.path.join(DATA, 'sdssv_calibration_lines.csv'); done = set()
    if os.path.exists(out):
        done = set(pd.read_csv(out).name); jobs = [j for j in jobs if j[0] not in done]
    print(f'{len(c)} calibration objects, {len(done)} already measured, {len(jobs)} to do ({int(changed.sum())} changed first)', flush=True)
    t0 = time.time(); rows = []; header = not os.path.exists(out)
    def flush():
        nonlocal rows, header
        if rows:
            pd.DataFrame(rows).to_csv(out, mode='a', header=header, index=False); header = False; rows = []
    with ThreadPoolExecutor(max_workers=a.threads) as ex:
        for k, (name, recs) in enumerate(ex.map(lambda j: f03d.one(j[0], j[1]), jobs)):
            for r in recs:
                L = r['lines']; g = lambda ln, key: L.get(ln, {}).get(key, np.nan)
                rows.append(dict(name=name, proprietary=bool(r.get('proprietary', False)), mjd=r['mjd'], program=r['program'], cls=r['meta'].get('class'), subclass=r['meta'].get('subclass'),
                                 z=r['meta'].get('z'), sn=r['meta'].get('sn_median_all'), EW_Hb_rest=r.get('ew', {}).get('Hb', np.nan), EW_Ha_rest=r.get('ew', {}).get('Ha', np.nan),
                                 Hb_area=g('H_beta', 'area'), Hb_area_err=g('H_beta', 'area_err'), Hb_ew=g('H_beta', 'ew'), Hb_sigma=g('H_beta', 'sigma'), OIII_area=g('OIII_5007', 'area'),
                                 OIII_ew=g('OIII_5007', 'ew'), Ha_area=g('H_alpha', 'area'), Ha_area_err=g('H_alpha', 'area_err'), Ha_ew=g('H_alpha', 'ew'), Ha_sigma=g('H_alpha', 'sigma'),
                                 MgII_area=g('MgII', 'area'), MgII_ew=g('MgII', 'ew'), cont5100=g('H_beta', 'cont')))
            if not recs:
                rows.append(dict(name=name, proprietary=False, mjd=np.nan, program='NONE'))
            if (k + 1) % 50 == 0:
                flush(); print(f'[{time.time()-t0:5.0f}s] {k+1}/{len(jobs)}', flush=True)
    flush(); print(f'done in {time.time()-t0:.0f}s -> {out}', flush=True)


if __name__ == '__main__':
    main()

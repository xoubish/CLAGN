"""
03d_fetch_spectra.py  --  download the archival SDSS spectra (all epochs, SDSS-I..V) for the listed targets and extract
a compact spectral history.

Sources: the `sas_url` of every SDSS epoch in data/spectra_epochs_<tag>.csv (DR17 lite files for SDSS-I..IV, DR19 v6_1_3
lite files incl. the 'allepoch' coadd for SDSS-V).  Files are cached in data/spectra_cache/<name>/ (~200 kB each).
DESI DR1 spectra are not downloaded here (per-healpix coadd files of 100+ MB); see SPARCL note in FOLLOWUP_PLAN.md.

Usage: /opt/anaconda3/bin/python 03d_fetch_spectra.py [targets csv ...]   (default: data/targets_*.csv; primaries + backups)

Outputs
  data/spectra_dl/<name>.json      per target: list of epochs {mjd, phase, program, url, wave (Å, 6 Å bins 3600-10300),
                                   flux (1e-17 erg/s/cm2/Å, median in bin), ivar_ok fraction, class, subclass, z, lines{...}}
  data/spectra_lines.csv           per epoch: pipeline (SPZLINE) single-Gaussian line fits for Hβ, [O III] 5007, Hα, Mg II:
                                   area (1e-17 erg/s/cm2), EW (Å), sigma (km/s), continuum; plus ratios useful for CLAGN:
                                   Hb_area/OIII_area, Ha_area/OIII_area
"""
import io, os, sys, glob, json, time
import numpy as np
import pandas as pd
import requests
from astropy.io import fits
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__)); DATA = os.path.join(HERE, 'data')
CACHE = os.path.join(DATA, 'spectra_cache'); OUT = os.path.join(DATA, 'spectra_dl')
GRID = np.arange(3600.0, 10300.0, 6.0)
LINES = {'H_beta': 4862.7, 'OIII_5007': 5008.2, 'H_alpha': 6564.6, 'MgII': 2800.3}


def fetch(url, path, tries=3):
    if os.path.exists(path) and os.path.getsize(path) > 10000:
        return path
    for attempt in range(tries):
        try:
            r = requests.get(url, timeout=180)
            if r.ok and len(r.content) > 10000:
                open(path, 'wb').write(r.content); return path
            if r.status_code == 404:
                return None
        except Exception:
            time.sleep(5 * (attempt + 1))
    return None


def parse(path):
    h = fits.open(path)
    d = h[1].data
    lam = 10 ** d['loglam']; flux = d['flux'].astype(float); ivar = d['ivar'].astype(float)
    good = ivar > 0
    # median-rebin onto a fixed 6 Å grid (compact, enough for broad lines)
    idx = np.searchsorted(GRID, lam) - 1
    fb = np.full(len(GRID), np.nan); ok = np.zeros(len(GRID))
    for i in np.unique(idx[(idx >= 0) & (idx < len(GRID))]):
        m = (idx == i) & good
        if m.sum():
            fb[i] = np.median(flux[m]); ok[i] = m.mean()
    # metadata: legacy lite files carry SPECOBJ; eBOSS carry SPALL; SDSS-V carry SPALL (exposure info) + ZALL (class, z)
    meta = {}
    for hdu in ['SPECOBJ', 'SPALL', 'ZALL']:
        if hdu in h and h[hdu].data is not None and len(h[hdu].data):
            row = h[hdu].data[0]; names = h[hdu].columns.names
            for k in ['CLASS', 'SUBCLASS', 'Z', 'ZWARNING', 'PLATE', 'MJD', 'FIBERID', 'FIELD', 'CATALOGID', 'SN_MEDIAN_ALL', 'PROGRAMNAME']:
                if k in names and k.lower() not in meta:
                    v = row[k]
                    if isinstance(v, (str, bytes, np.str_, np.bytes_)):
                        meta[k.lower()] = (v.decode() if isinstance(v, (bytes, np.bytes_)) else str(v)).strip()
                    elif np.ndim(v) == 0:
                        meta[k.lower()] = float(v)
    z = meta.get('z')
    # uniform, model-free line indices in the rest frame (EW in Å): H-beta 4800-4930 vs continua 4700-4790 / 5090-5150,
    # H-alpha 6480-6650 vs 6350-6450 / 6700-6800.  Positive EW = emission.  Uses the raw pixels, not the rebinned grid.
    ew = {}
    if z is not None and np.isfinite(z) and z >= 0:
        rest = lam / (1 + z)
        for nm, (l0, l1, c0, c1, c2, c3) in {'Hb': (4800, 4930, 4700, 4790, 5090, 5150), 'Ha': (6480, 6650, 6350, 6450, 6700, 6800)}.items():
            ml = (rest > l0) & (rest < l1) & good; mc = (((rest > c0) & (rest < c1)) | ((rest > c2) & (rest < c3))) & good
            if ml.sum() > 10 and mc.sum() > 10:
                cont = np.polyfit(rest[mc], flux[mc], 1); cfit = np.polyval(cont, rest[ml])
                with np.errstate(divide='ignore', invalid='ignore'):
                    ew[nm] = float(np.nansum((flux[ml] / cfit - 1) * np.gradient(rest[ml]))) if np.all(cfit > 0) else np.nan
    lines = {}
    zl_name = 'SPZLINE' if 'SPZLINE' in h else ('ZLINE' if 'ZLINE' in h else None)
    if zl_name and h[zl_name].data is not None and len(h[zl_name].data):
        zl = h[zl_name].data; cols = zl.columns.names
        def col(c, j, default=np.nan):
            return float(zl[c][j]) if c in cols else default
        for nm, w0 in LINES.items():
            j = int(np.argmin(np.abs(zl['LINEWAVE'] - w0)))
            if abs(zl['LINEWAVE'][j] - w0) < 3:
                lines[nm] = dict(area=col('LINEAREA', j), area_err=col('LINEAREA_ERR', j), ew=col('LINEEW', j), ew_err=col('LINEEW_ERR', j),
                                 sigma=col('LINESIGMA', j), cont=col('LINECONTLEVEL', j), npix=int(col('LINENPIX', j, 0)))
    return dict(wave=[round(float(x), 1) for x in GRID], flux=[None if np.isnan(v) else round(float(v), 3) for v in fb],
                ok=[round(float(v), 2) for v in ok], meta=meta, lines=lines, ew=ew)


def desi_sparcl(name, ra, dec):
    """DESI DR1 (and BOSS/SDSS DR17 as a cross-check) spectra from SPARCL, rebinned like the SDSS files."""
    try:
        from sparcl.client import SparclClient
    except Exception:
        return []
    out = []
    try:
        c = SparclClient()
        found = c.find(outfields=['sparcl_id', 'ra', 'dec', 'redshift', 'spectype', 'data_release', 'specid'],
                       constraints={'data_release': ['DESI-DR1'], 'ra': [ra - 0.0006, ra + 0.0006], 'dec': [dec - 0.0006, dec + 0.0006]}, limit=10)
        if not found.records:
            return []
        got = c.retrieve(uuid_list=found.ids, include=['sparcl_id', 'specid', 'data_release', 'redshift', 'spectype', 'flux', 'wavelength', 'ivar', 'ra', 'dec'])
        for r in got.records:
            lam = np.asarray(r.wavelength); flux = np.asarray(r.flux, float); ivar = np.asarray(r.ivar, float); good = ivar > 0
            idx = np.searchsorted(GRID, lam) - 1
            fb = np.full(len(GRID), np.nan); ok = np.zeros(len(GRID))
            for i in np.unique(idx[(idx >= 0) & (idx < len(GRID))]):
                m = (idx == i) & good
                if m.sum():
                    fb[i] = np.median(flux[m]); ok[i] = m.mean()
            out.append(dict(wave=[round(float(x), 1) for x in GRID], flux=[None if np.isnan(v) else round(float(v), 3) for v in fb],
                            ok=[round(float(v), 2) for v in ok], meta=dict(**{'class': str(r.spectype), 'z': float(r.redshift), 'specid': str(r.specid)}),
                            lines={}, mjd=np.nan, phase=None, program=f'DESI {r.data_release}', coadd=True, url=f'sparcl:{r.sparcl_id}', source='DESI'))
    except Exception as ex:
        print(f'   {name}: SPARCL failed {str(ex)[:60]}', flush=True)
    return out


def one(name, epochs, ra=None, dec=None):
    os.makedirs(os.path.join(CACHE, name), exist_ok=True)
    out = []
    epochs = epochs.drop_duplicates('sas_url')
    for e in epochs.sort_values('mjd').itertuples():
        url = str(e.sas_url)
        if not url.startswith('http'):
            continue
        path = os.path.join(CACHE, name, os.path.basename(url))
        p = fetch(url, path)
        if p is None:
            continue
        try:
            rec = parse(p)
        except Exception as ex:
            print(f'   {name} {os.path.basename(url)}: parse failed {str(ex)[:60]}', flush=True); continue
        rec.update(mjd=float(e.mjd), phase=(int(e.sdss_phase) if pd.notna(e.sdss_phase) else None), program=str(getattr(e, 'programname', '')),
                   coadd=bool(getattr(e, 'is_coadd', False)), url=url, source='SDSS')
        out.append(rec)
    # several allspec rows can point at different reductions of the same night (daily / epoch / allepoch coadds);
    # keep one spectrum per (rounded MJD), the one with the highest S/N, but always keep the allepoch coadd too
    best = {}
    for r in out:
        key = ('coadd' if r['coadd'] else round(r['mjd']))
        sn = r['meta'].get('sn_median_all') or 0
        if key not in best or sn > (best[key]['meta'].get('sn_median_all') or 0):
            best[key] = r
    out = sorted(best.values(), key=lambda r: r['mjd'])
    if ra is not None:
        desi = desi_sparcl(name, ra, dec)
        seen = set(); dd = []
        for d in desi:
            if d['meta'].get('specid') in seen:
                continue
            seen.add(d['meta'].get('specid')); dd.append(d)
        dm = epochs_desi_mjd.get(name)
        for d in dd:
            d['mjd'] = float(dm) if dm is not None else np.nan
            # same model-free EW indices for DESI, from the rebinned grid
            zz = d['meta'].get('z'); w = np.asarray(d['wave']); f = np.asarray([np.nan if v is None else v for v in d['flux']], float)
            if zz is not None and np.isfinite(zz):
                rest = w / (1 + zz); d['ew'] = {}
                for nm, (l0, l1, c0, c1, c2, c3) in {'Hb': (4800, 4930, 4700, 4790, 5090, 5150), 'Ha': (6480, 6650, 6350, 6450, 6700, 6800)}.items():
                    ml = (rest > l0) & (rest < l1) & np.isfinite(f); mc = (((rest > c0) & (rest < c1)) | ((rest > c2) & (rest < c3))) & np.isfinite(f)
                    if ml.sum() > 5 and mc.sum() > 5:
                        cont = np.polyfit(rest[mc], f[mc], 1); cfit = np.polyval(cont, rest[ml])
                        d['ew'][nm] = float(np.nansum((f[ml] / cfit - 1) * np.gradient(rest[ml]))) if np.all(cfit > 0) else np.nan
        out += dd
    return name, out


epochs_desi_mjd = {}


def main():
    files = sys.argv[1:] or sorted(glob.glob(os.path.join(DATA, 'targets_*.csv')))
    t = pd.concat([pd.read_csv(f) for f in files]).drop_duplicates('name')
    epall = pd.concat([pd.read_csv(p, low_memory=False) for p in glob.glob(os.path.join(DATA, 'spectra_epochs_*.csv'))], ignore_index=True)
    epall = epall[epall.name.isin(t.name)]
    for n, g in epall[epall.source == 'DESI'].groupby('name'):
        epochs_desi_mjd[n] = float(g.mjd.median())
    ep = epall[(epall.source == 'SDSS') & epall.sas_url.notna()]
    pos = t.set_index('name')
    os.makedirs(OUT, exist_ok=True); t0 = time.time(); rows = []
    names = sorted(set(ep.name) | set(epochs_desi_mjd))
    jobs = [(n, ep[ep.name == n], float(pos.at[n, 'ra']), float(pos.at[n, 'dec'])) for n in names]
    print(f'{len(t)} targets, {len(jobs)} with archival spectra ({len(ep.drop_duplicates("sas_url"))} SDSS files, {len(epochs_desi_mjd)} with DESI)', flush=True)
    with ThreadPoolExecutor(max_workers=6) as ex:
        for k, (name, recs) in enumerate(ex.map(lambda a: one(*a), jobs)):
            json.dump(recs, open(os.path.join(OUT, f'{name}.json'), 'w'), separators=(',', ':'))
            for r in recs:
                L = r['lines']; g = lambda ln, key: L.get(ln, {}).get(key, np.nan)
                rows.append(dict(name=name, source=r.get('source', 'SDSS'), mjd=r['mjd'], phase=r['phase'], program=r['program'], coadd=r['coadd'], cls=r['meta'].get('class'),
                                 subclass=r['meta'].get('subclass'), z=r['meta'].get('z'), sn=r['meta'].get('sn_median_all'),
                                 EW_Hb_rest=r.get('ew', {}).get('Hb', np.nan), EW_Ha_rest=r.get('ew', {}).get('Ha', np.nan),
                                 Hb_area=g('H_beta', 'area'), Hb_area_err=g('H_beta', 'area_err'), Hb_ew=g('H_beta', 'ew'), Hb_sigma=g('H_beta', 'sigma'),
                                 OIII_area=g('OIII_5007', 'area'), OIII_ew=g('OIII_5007', 'ew'),
                                 Ha_area=g('H_alpha', 'area'), Ha_area_err=g('H_alpha', 'area_err'), Ha_ew=g('H_alpha', 'ew'), Ha_sigma=g('H_alpha', 'sigma'),
                                 MgII_area=g('MgII', 'area'), MgII_ew=g('MgII', 'ew'), cont5100=g('H_beta', 'cont')))
            if (k + 1) % 20 == 0:
                print(f'[{time.time()-t0:5.0f}s] {k+1}/{len(jobs)}', flush=True)
    L = pd.DataFrame(rows)
    if len(L):
        L['Hb_over_OIII'] = L.Hb_area / L.OIII_area.where(L.OIII_area > 0)
        L['Ha_over_OIII'] = L.Ha_area / L.OIII_area.where(L.OIII_area > 0)
    L.to_csv(os.path.join(DATA, 'spectra_lines.csv'), index=False)
    print(f'done: {len(L)} epoch spectra for {L.name.nunique() if len(L) else 0} targets in {time.time()-t0:.0f}s -> data/spectra_dl/*.json, data/spectra_lines.csv')


if __name__ == '__main__':
    main()

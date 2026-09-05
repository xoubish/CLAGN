"""
03f_sdssv_internal.py  --  proprietary SDSS-V BOSS epochs for our targets from the collaboration's Science Archive Server.

Public DR19 (used by 03_spectra_inventory.py) ends at MJD 60130 (2023-06).  The internal spAll of the newest pipeline
version carries every BHM epoch since then.  This script:
  1. lists https://data.sdss5.org/sas/sdsswork/bhm/boss/spectro/redux/ and picks the newest tagged version (or --version),
  2. downloads spAll-lite-<version>.fits.gz once into data/sdssv_internal/ (a few hundred MB, resumable cache),
  3. matches every row within 2" of a master-list position (data/master_list_scored.csv, or the CSVs given),
  4. writes  data/sdssv_internal_epochs.csv        -> read by 04_score_tiers.py (recency, class change, r at the new epoch)
             data/spectra_epochs_sdssvint.csv      -> read by 03d_fetch_spectra.py, which then downloads the spec-lite files
     Rows already public in DR19 (MJD <= 60130 for the same object) are dropped, so only genuinely new epochs are added.

Credentials: ~/.netrc with
    machine data.sdss5.org
        login <sdss-v username>
        password <sdss-v password>
(chmod 600).  They are handed out by phone or in person only; the IPAC list is in SDSSV-Data-Access.pdf.

Everything this script writes is PROPRIETARY collaboration data: data/sdssv_internal/ and both CSVs are git-ignored,
03d marks the spectra `proprietary`, and 07_make_webpage.py leaves them out of the public docs/index.html copy.

Usage: /opt/anaconda3/bin/python 03f_sdssv_internal.py [--version v6_2_1|master] [--coadd daily|epoch] [--dry-run] [targets.csv ...]
"""
import os, re, sys, glob, time, netrc, argparse
import numpy as np
import pandas as pd
import requests
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, 'data')
CACHE = os.path.join(DATA, 'sdssv_internal')
HOST = 'data.sdss5.org'
REDUX = f'https://{HOST}/sas/sdsswork/bhm/boss/spectro/redux'
DR19_LAST_MJD = 60130          # public DR19 BOSS cutoff; internal epochs up to here are already in spectra_epochs_*.csv
MATCH_ARCSEC = 2.0


def auth():
    try:
        rc = netrc.netrc(os.path.expanduser('~/.netrc'))
        a = rc.authenticators(HOST)
    except (FileNotFoundError, netrc.NetrcParseError):
        a = None
    if a is None:
        sys.exit(f'no credentials: add a "machine {HOST}" entry to ~/.netrc (see the docstring), chmod 600, and rerun')
    return (a[0], a[2])


def list_versions(s):
    r = s.get(REDUX + '/', timeout=60)
    if r.status_code == 401:
        sys.exit('HTTP 401 from the SAS: the .netrc credentials were rejected')
    r.raise_for_status()
    vers = sorted(set(re.findall(r'href="(v\d+_\d+_\d+)/"', r.text)), key=lambda v: [int(x) for x in v[1:].split('_')])
    return vers, ('master' in r.text)


def download(s, url, path, tries=3):
    """Resumable download with a progress line; returns path (cached if complete)."""
    done = os.path.join(path + '.ok')
    if os.path.exists(done) and os.path.exists(path):
        return path
    os.makedirs(os.path.dirname(path), exist_ok=True)
    for attempt in range(tries):
        have = os.path.getsize(path) if os.path.exists(path) else 0
        hdr = {'Range': f'bytes={have}-'} if have else {}
        with s.get(url, headers=hdr, stream=True, timeout=300) as r:
            if r.status_code == 416:                    # already complete
                break
            if r.status_code == 200 and have:           # server ignored Range: start over
                have = 0
            r.raise_for_status()
            total = int(r.headers.get('Content-Length', 0)) + have
            t0 = time.time()
            with open(path, 'ab' if have else 'wb') as f:
                for chunk in r.iter_content(1 << 20):
                    f.write(chunk); have += len(chunk)
                    if int(time.time() - t0) % 15 == 0:
                        print(f'\r   {os.path.basename(path)}: {have/1e6:.0f}/{total/1e6:.0f} MB', end='', flush=True)
            print()
            if not total or have >= total:
                break
    open(done, 'w').write('ok')
    return path


def spall_table(path):
    """Columns we need from spAll-lite, robust to the column renames between v6_0 and v6_2."""
    with fits.open(path, memmap=False) as h:
        d = h[1].data; names = set(d.columns.names)
        def col(*cands, default=np.nan):
            for c in cands:
                if c in names:
                    return d[c]
            return np.full(len(d), default)
        t = pd.DataFrame(dict(
            ra=col('FIBER_RA', 'PLUG_RA', 'RACAT').astype(float), dec=col('FIBER_DEC', 'PLUG_DEC', 'DECCAT').astype(float),
            field=col('FIELD', 'PLATE').astype(int), mjd=col('MJD').astype(int), catalogid=col('CATALOGID').astype(np.int64),
            sdss_id=col('SDSS_ID', default=-1).astype(np.int64),
            cls=[str(x).strip() for x in col('CLASS', default='')], subclass=[str(x).strip() for x in col('SUBCLASS', default='')],
            z=col('Z').astype(float), zwarning=col('ZWARNING').astype(int), sn_median_all=col('SN_MEDIAN_ALL').astype(float),
            firstcarton=[str(x).strip() for x in col('FIRSTCARTON', default='')], programname=[str(x).strip() for x in col('PROGRAMNAME', default='')],
            objtype=[str(x).strip() for x in col('OBJTYPE', default='')], nexp=col('NEXP', default=0).astype(int), exptime=col('EXPTIME', default=0).astype(float),
        ))
        sf = np.asarray(col('SPECTROFLUX'), dtype=float)          # FITS arrays are big-endian: cast before pandas sees them
        t['spectroflux_g'] = sf[:, 1] if sf.ndim == 2 else np.nan
        t['spectroflux_r'] = sf[:, 2] if sf.ndim == 2 else np.nan
        t['spectroflux_i'] = sf[:, 3] if sf.ndim == 2 else np.nan
    return t


def spec_url(version, coadd, field, mjd, catalogid):
    """spec-lite path conventions of idlspec2d v6_1+ (identical to the public DR19 tree, host and root aside)."""
    if coadd == 'allepoch':
        return f'{REDUX}/{version}/spectra/lite/allepoch/{mjd:5d}/spec-allepoch-{mjd:5d}-{catalogid}.fits'
    sub = 'epoch/spectra' if coadd == 'epoch' else 'spectra'
    return f'{REDUX}/{version}/{sub}/lite/{field:06d}/{mjd:5d}/spec-{field:06d}-{mjd:5d}-{catalogid}.fits'


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--version', default=None, help='pipeline version directory, e.g. v6_2_1 or master (default: newest tagged)')
    ap.add_argument('--coadd', default='daily', choices=['daily', 'epoch'], help='which spAll to use (daily = one row per visit)')
    ap.add_argument('--dry-run', action='store_true', help='only list the versions available on the SAS')
    ap.add_argument('--keep-public', action='store_true', help='also keep epochs already covered by DR19 (MJD <= %d)' % DR19_LAST_MJD)
    ap.add_argument('--public-dr19', action='store_true', help='TEST MODE: run the same code on the public DR19 spAll-lite (no credentials); writes to data/_sdssv_test_*.csv')
    ap.add_argument('targets', nargs='*', help='CSV(s) with name, ra, dec (default: data/master_list_scored.csv)')
    a = ap.parse_args()

    global REDUX
    s = requests.Session()
    if a.public_dr19:
        REDUX = 'https://data.sdss.org/sas/dr19/spectro/boss/redux'; a.version = 'v6_1_3'; a.keep_public = True
        vers, has_master = ['v6_1_3'], False
    else:
        s.auth = auth()
        vers, has_master = list_versions(s)
    print(f'SAS redux versions: {", ".join(vers)}{" + master" if has_master else ""}')
    version = a.version or (vers[-1] if vers else 'master')
    print(f'using {version} ({a.coadd} coadd)')
    if a.dry_run:
        return

    fname = f'spAll-lite-{version}.fits.gz'
    url = f'{REDUX}/{version}/' + ('epoch/' if a.coadd == 'epoch' else '') + fname
    path = download(s, url, os.path.join(CACHE, 'public_dr19' if a.public_dr19 else a.coadd, fname))
    t0 = time.time(); sp = spall_table(path); print(f'spAll rows: {len(sp)}  (read in {time.time()-t0:.0f}s)')
    sp = sp[np.isfinite(sp.ra) & np.isfinite(sp.dec)]

    files = a.targets or [os.path.join(DATA, 'master_list_scored.csv')]
    tg = pd.concat([pd.read_csv(f, low_memory=False, usecols=lambda c: c in ('name', 'ra', 'dec')) for f in files]).drop_duplicates('name')
    ct = SkyCoord(tg.ra.values * u.deg, tg.dec.values * u.deg); cs = SkyCoord(sp.ra.values * u.deg, sp.dec.values * u.deg)
    # every target within the radius, not only the nearest: the master list can hold one object under two names
    # (a Zeltyn J-name and a pool P-name), and both inventories must see the epoch
    isp, it, sep, _ = ct.search_around_sky(cs, MATCH_ARCSEC * u.arcsec)      # returns (idx into the argument, idx into self)
    m = sp.iloc[isp].copy(); m['name'] = tg.name.values[it]; m['sep_arcsec'] = sep.arcsec
    print(f'{len(m)} spAll rows within {MATCH_ARCSEC}" of {m.name.nunique()} of our {len(tg)} targets')
    if not a.keep_public:
        m = m[m.mjd > DR19_LAST_MJD]
        print(f'{len(m)} epochs newer than DR19 (MJD > {DR19_LAST_MJD}) for {m.name.nunique()} targets')

    # 04_score_tiers.py schema (+ provenance)
    m['run2d'] = version; m['coadd'] = a.coadd; m['sdss_phase'] = 5; m['source'] = 'SDSS'; m['proprietary'] = True
    m['is_coadd'] = False
    m['sas_url'] = [spec_url(version, a.coadd, f, j, c) for f, j, c in zip(m.field, m.mjd, m.catalogid)]
    m = m.rename(columns={'cls': 'class'}).sort_values(['name', 'mjd'])
    out04 = os.path.join(DATA, '_sdssv_test_epochs.csv' if a.public_dr19 else 'sdssv_internal_epochs.csv'); m.to_csv(out04, index=False)
    m['proprietary'] = not a.public_dr19
    # 03d_fetch_spectra.py schema (its glob spectra_epochs_*.csv picks this up; file is git-ignored)
    cols03 = ['name', 'ra', 'dec', 'mjd', 'sdss_phase', 'run2d', 'coadd', 'programname', 'firstcarton', 'class', 'subclass', 'z', 'zwarning',
              'sn_median_all', 'spectroflux_g', 'spectroflux_r', 'spectroflux_i', 'catalogid', 'field', 'sas_url', 'is_coadd', 'source', 'proprietary']
    m[cols03].rename(columns={'field': 'plate_or_fps_field'}).to_csv(os.path.join(DATA, '_sdssv_test_epochs03d.csv' if a.public_dr19 else 'spectra_epochs_sdssvint.csv'), index=False)

    prim = set(pd.concat([pd.read_csv(f) for f in glob.glob(os.path.join(DATA, 'targets_*.csv'))]).query('rank > 0').name) if glob.glob(os.path.join(DATA, 'targets_*.csv')) else set()
    print(f'wrote {out04} ({len(m)} rows); primaries with a new internal epoch: {len(prim & set(m.name))}')
    if len(m):
        last = m.groupby('name').agg(n=('mjd', 'size'), last_mjd=('mjd', 'max'), cls=('class', 'last'), carton=('firstcarton', 'first'))
        print(last[last.index.isin(prim)].to_string() if len(prim) else last.head(20).to_string())
    print('next: ./rescore.sh   (04 uses the epochs, 03d downloads the spectra with the same .netrc, 07 keeps them out of docs/)')


if __name__ == '__main__':
    main()

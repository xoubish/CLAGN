"""
03f_sdssv_internal.py  --  proprietary SDSS-V BOSS epochs for our targets from the collaboration's Science Archive Server.

Public DR19 (used by 03_spectra_inventory.py) ends at MJD 60130 (2023-06).  The internal spAll of the newest pipeline
version carries every BHM epoch since then.  This script:
  1. lists https://data.sdss5.org/sas/sdsswork/bhm/boss/spectro/redux/ and picks the newest tagged version (or --version),
  2. downloads <version>/summary/<coadd>/spAll-lite-<version>.fits.gz once into data/sdssv_internal/<coadd>/ (v6_2_1 daily:
     4.6 GB gzipped, 15.8 M rows x 1193 B = 18.8 GB uncompressed; resumable, parallel-range download, cached),
     then STREAMS it: rows are decoded chunk by chunk from the gzip and only those in an (RA, Dec) cell touched by a
     target survive, so the table is never held in memory,
  3. matches every surviving row within 2" of a master-list position (data/master_list_scored.csv, or the CSVs given),
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
import gzip
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


def download(s, url, path, tries=3, workers=8, seg=64 << 20):
    """Resumable download; returns path (cached once <path>.ok exists).
    The SAS delivers only ~0.5 MB/s per connection but accepts many, so the file is fetched as `workers` parallel byte
    ranges of `seg` bytes written in place (os.pwrite); finished segments are listed in <path>.parts, so an interrupted
    run resumes.  A plain sequential download is the fallback when the server does not advertise byte ranges."""
    done = path + '.ok'
    if os.path.exists(done) and os.path.exists(path):
        return path
    os.makedirs(os.path.dirname(path), exist_ok=True)
    h = s.head(url, timeout=120, allow_redirects=True)
    if h.status_code == 401:
        sys.exit('HTTP 401 from the SAS: the credentials were rejected')
    h.raise_for_status()
    total = int(h.headers.get('Content-Length', 0))
    if total and h.headers.get('Accept-Ranges') == 'bytes':
        import json, threading
        from concurrent.futures import ThreadPoolExecutor
        parts = path + '.parts'
        finished = set(json.load(open(parts))) if os.path.exists(parts) else set()
        if not finished and os.path.exists(path):                      # partial file from the sequential fallback: its prefix is good
            finished = set(range(os.path.getsize(path) // seg))
        with open(path, 'ab') as f:
            if f.tell() < total:
                f.truncate(total)                                      # sparse preallocation
        nseg = -(-total // seg); todo = [i for i in range(nseg) if i not in finished]
        fd = os.open(path, os.O_RDWR); lock = threading.Lock(); t0 = time.time(); got = [len(finished) * seg]; new = [0]
        print(f'   {os.path.basename(path)}: {len(finished)}/{nseg} segments cached, fetching {len(todo)} with {workers} connections')
        def fetch(i):
            lo, hi = i * seg, min(total, (i + 1) * seg) - 1
            for attempt in range(tries):
                try:
                    with s.get(url, headers={'Range': f'bytes={lo}-{hi}'}, stream=True, timeout=300) as r:
                        if r.status_code != 206:
                            raise IOError(f'HTTP {r.status_code} for bytes {lo}-{hi}')
                        off = lo
                        for chunk in r.iter_content(1 << 20):
                            os.pwrite(fd, chunk, off); off += len(chunk)
                        if off != hi + 1:
                            raise IOError(f'short range: got {off - lo} of {hi + 1 - lo} bytes')
                    with lock:
                        finished.add(i); json.dump(sorted(finished), open(parts, 'w')); got[0] += hi + 1 - lo; new[0] += hi + 1 - lo
                        rate = new[0] / 1e6 / max(time.time() - t0, 1)
                        print(f'\r   {os.path.basename(path)}: {min(got[0], total)/1e6:.0f}/{total/1e6:.0f} MB  {rate:.1f} MB/s, ~{(total - got[0]) / 1e6 / max(rate, 0.01) / 60:.0f} min left', end='', flush=True)
                    return
                except Exception as e:
                    err = e; time.sleep(5 * (attempt + 1))
            raise err
        with ThreadPoolExecutor(workers) as ex:
            list(ex.map(fetch, todo))
        os.close(fd); print()
        if os.path.exists(parts):
            os.remove(parts)
        open(done, 'w').write('ok')
        return path
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


def spall_url(version, coadd, public=False):
    """Where the spAll-lite summary table lives.
    Internal SAS (idlspec2d v6_2+):  <version>/summary/<coadd>/spAll-lite-<version>[-epoch|-allepoch].fits.gz
                                     (v6_2_1: daily 4.6 GB, epoch 5.3 GB, allepoch 0.5 GB; master = rolling daily reduction)
    Public DR19 tree (v6_1_3):       <version>/[epoch/]spAll-lite-<version>.fits.gz"""
    if public:
        return f'{REDUX}/{version}/' + ('epoch/' if coadd == 'epoch' else '') + f'spAll-lite-{version}.fits.gz'
    return f'{REDUX}/{version}/summary/{coadd}/spAll-lite-{version}{"" if coadd == "daily" else "-" + coadd}.fits.gz'


FITS_DTYPES = {'L': 'S1', 'B': 'u1', 'I': '>i2', 'J': '>i4', 'K': '>i8', 'E': '>f4', 'D': '>f8'}


def _fits_row_dtype(hdr):
    """numpy dtype of one binary-table row from its header (fixed-width columns only, which is all spAll has)."""
    fields = []
    for i in range(1, hdr['TFIELDS'] + 1):
        name, form = hdr[f'TTYPE{i}'], hdr[f'TFORM{i}'].strip()
        mo = re.fullmatch(r'(\d*)([A-Z])(.*)', form)
        n, code = int(mo.group(1) or 1), mo.group(2)
        if code == 'A':
            fields.append((name, f'S{n}'))
        elif code in FITS_DTYPES:
            fields.append((name, FITS_DTYPES[code]) if n == 1 else (name, FITS_DTYPES[code], (n,)))
        else:
            raise ValueError(f'column {name}: unsupported TFORM {form}')
    dt = np.dtype(fields)
    if dt.itemsize != hdr['NAXIS1']:
        raise ValueError(f'row size mismatch: dtype {dt.itemsize} B vs NAXIS1 {hdr["NAXIS1"]} B')
    return dt


def _read_header(f):
    """Read 2880-byte blocks from a stream until the END card; returns the astropy Header."""
    raw = b''
    while True:
        blk = f.read(2880)
        if len(blk) < 2880:
            raise IOError('truncated FITS header')
        raw += blk
        if any(blk[j:j + 80] == b'END'.ljust(80) for j in range(0, 2880, 80)):
            return fits.Header.fromstring(raw)


def _skip_data(f, hdr):
    n = abs(hdr.get('BITPIX', 8)) // 8
    for i in range(1, hdr.get('NAXIS', 0) + 1):
        n *= hdr[f'NAXIS{i}']
    if hdr.get('NAXIS', 0) and n:
        f.read(n + (-n) % 2880)


def radec(arr):
    """(ra, dec) as native floats from a structured spAll chunk, whatever the column generation calls them."""
    ra = next(c for c in ('FIBER_RA', 'PLUG_RA', 'RACAT') if c in arr.dtype.names)
    dec = next(c for c in ('FIBER_DEC', 'PLUG_DEC', 'DECCAT') if c in arr.dtype.names)
    return arr[ra].astype(float), arr[dec].astype(float)


def stream_spall(path, select, convert=None, chunk_rows=250_000):
    """Stream the gzipped spAll-lite chunk by chunk (~300 MB in flight, the 18.8 GB table is never held in memory).
    select(arr)  -> boolean mask over a structured-array chunk (which rows to keep);
    convert(arr) -> optional reduction of the kept rows (e.g. spall_table, or a compact recarray) applied per chunk.
    Returns the concatenation of the kept (converted) chunks: a structured array, or whatever convert returns
    (DataFrames are concatenated with pandas)."""
    keep = []; nkept = 0; t0 = time.time()
    with gzip.open(path, 'rb') as f:
        h0 = _read_header(f); _skip_data(f, h0)
        h1 = _read_header(f); dt = _fits_row_dtype(h1); nrows = h1['NAXIS2']
        print(f'   spAll-lite: {nrows:,} rows x {dt.itemsize} B ({nrows * dt.itemsize / 1e9:.1f} GB), {len(dt.names)} columns; streaming ...')
        done = 0
        while done < nrows:
            n = min(chunk_rows, nrows - done)
            buf = f.read(n * dt.itemsize)
            if len(buf) < n * dt.itemsize:
                raise IOError(f'truncated spAll data after {done + len(buf) // dt.itemsize:,} rows: delete {path} and its .ok file and rerun')
            arr = np.frombuffer(buf, dtype=dt, count=n)
            sel = select(arr)
            if sel.any():
                part = arr[sel].copy(); nkept += len(part)
                keep.append(convert(part) if convert else part)
            done += n
            if done % (chunk_rows * 8) == 0 or done == nrows:
                print(f'\r   {done / 1e6:5.1f} M rows, {nkept:,} kept, {time.time() - t0:.0f}s', end='', flush=True)
    print()
    if not keep:
        return np.zeros(0, dtype=dt) if convert is None else convert(np.zeros(0, dtype=dt))
    return pd.concat(keep, ignore_index=True) if isinstance(keep[0], pd.DataFrame) else np.concatenate(keep)


def spall_rows_near(path, ra_t, dec_t, cell_deg=0.01, chunk_rows=250_000):
    """Structured array of the spAll rows that fall in an (RA, Dec) cell touched by a target (cell = 36", plus the 8
    neighbours, so every row within 2" of a target is kept for the exact match in main()); a few thousand rows survive."""
    ncell = int(round(360 / cell_deg))
    def key(ra, dec):
        return (np.floor(ra / cell_deg).astype(np.int64) % ncell) * 100_000 + np.floor((dec + 90) / cell_deg).astype(np.int64)
    ira = np.floor(ra_t / cell_deg).astype(np.int64); idec = np.floor((dec_t + 90) / cell_deg).astype(np.int64)
    keys = np.unique(np.concatenate([((ira + di) % ncell) * 100_000 + idec + dj for di in (-1, 0, 1) for dj in (-1, 0, 1)]))
    def select(arr):
        ra, dec = radec(arr)
        ok = np.isfinite(ra) & np.isfinite(dec) & (np.abs(dec) <= 90)
        return ok & np.isin(key(np.where(ok, ra, 0.0), np.where(ok, dec, 0.0)), keys)
    return stream_spall(path, select, chunk_rows=chunk_rows)


def spall_table(d):
    """Columns we need, from a structured array of spAll-lite rows; robust to the column renames between v6_0 and v6_2."""
    names = set(d.dtype.names)
    def col(*cands, default=np.nan):
        for c in cands:
            if c in names:
                return d[c]
        return np.full(len(d), default)
    def s(*cands):
        v = col(*cands, default=b'')
        return [x.decode('ascii', 'replace').strip() if isinstance(x, bytes) else str(x).strip() for x in v]
    t = pd.DataFrame(dict(
        ra=col('FIBER_RA', 'PLUG_RA', 'RACAT').astype(float), dec=col('FIBER_DEC', 'PLUG_DEC', 'DECCAT').astype(float),
        field=col('FIELD', 'PLATE').astype(int), mjd=col('MJD').astype(int), catalogid=col('CATALOGID').astype(np.int64),
        sdss_id=col('SDSS_ID', default=-1).astype(np.int64), spec_file=s('SPEC_FILE'), obs=s('OBS'), run2d_row=s('RUN2D'),
        cls=s('CLASS'), subclass=s('SUBCLASS'),
        z=col('Z').astype(float), zwarning=col('ZWARNING').astype(int), sn_median_all=col('SN_MEDIAN_ALL').astype(float),
        firstcarton=s('FIRSTCARTON'), programname=s('PROGRAMNAME'), objtype=s('OBJTYPE'),
        nexp=col('NEXP', default=0).astype(int), exptime=col('EXPTIME', default=0).astype(float),
        fieldquality=s('FIELDQUALITY'), specprimary=col('SPECPRIMARY', default=-1).astype(int), nspecobs=col('NSPECOBS', default=-1).astype(int),
    ))
    sf = np.asarray(col('SPECTROFLUX'), dtype=float)          # FITS arrays are big-endian: cast before pandas sees them
    t['spectroflux_g'] = sf[:, 1] if sf.ndim == 2 else np.nan
    t['spectroflux_r'] = sf[:, 2] if sf.ndim == 2 else np.nan
    t['spectroflux_i'] = sf[:, 3] if sf.ndim == 2 else np.nan
    return t


def spec_url(version, coadd, field, mjd, catalogid, spec_file='', public=False):
    """Where one spec-lite file lives (the spAll row's SPEC_FILE name when present, else the spec-<field>-<mjd>-<catalogid> form).
    Internal SAS (v6_2+):  <version>/spectra/<coadd>/lite/<field//1000>XXX/<field>/<mjd>/<file>
                           e.g. v6_2_1/spectra/daily/lite/030XXX/030073/60188/spec-030073-60188-27021599486824279.fits
    Public DR19 (v6_1_3):  <version>/[epoch/]spectra/lite/<field>/<mjd>/<file>"""
    name = spec_file or f'spec-{field:06d}-{mjd:5d}-{catalogid}.fits'
    if public:
        return f'{REDUX}/{version}/{"epoch/spectra" if coadd == "epoch" else "spectra"}/lite/{field:06d}/{mjd:5d}/{name}'
    return f'{REDUX}/{version}/spectra/{coadd}/lite/{field // 1000:03d}XXX/{field:06d}/{mjd:5d}/{name}'


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

    files = a.targets or [os.path.join(DATA, 'master_list_scored.csv')]
    tg = pd.concat([pd.read_csv(f, low_memory=False, usecols=lambda c: c in ('name', 'ra', 'dec')) for f in files]).drop_duplicates('name')
    tg = tg[np.isfinite(tg.ra) & np.isfinite(tg.dec)].reset_index(drop=True)

    url = spall_url(version, a.coadd, public=a.public_dr19)
    path = download(s, url, os.path.join(CACHE, 'public_dr19' if a.public_dr19 else a.coadd, os.path.basename(url)))
    t0 = time.time(); sp = spall_table(spall_rows_near(path, tg.ra.values, tg.dec.values))
    print(f'{len(sp)} spAll rows in cells around our {len(tg)} targets  (streamed in {time.time()-t0:.0f}s)')

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
    m['sas_url'] = [spec_url(version, a.coadd, f, j, c, sf, public=a.public_dr19) for f, j, c, sf in zip(m.field, m.mjd, m.catalogid, m.spec_file)]
    m = m.rename(columns={'cls': 'class'}).sort_values(['name', 'mjd'])
    if len(m):                                       # one HEAD request: catches a changed spec-lite layout before 03d silently skips every file
        r = s.head(m.sas_url.iloc[0], timeout=60, allow_redirects=True)
        print(f'spec-lite path check: HTTP {r.status_code} for {m.sas_url.iloc[0]}' + ('' if r.ok else '   <-- LAYOUT CHANGED? fix spec_url()'))
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

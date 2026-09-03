"""
07b_cutouts.py  --  image cutouts for the night-sheet cards, cached as small PNG/JPEG files in data/cutouts/<name>_<kind>.

  sdss : SDSS DR18 gri colour composite JPEG from SkyServer ImgCutout (64" x 64", 0.4"/pix -> 160 px)
  ztf_g, ztf_r : ZTF reference-image (deep stack) cutouts from the IRSA IBE service, 64" x 64", asinh-stretched PNG

Usage: /opt/anaconda3/bin/python 07b_cutouts.py [targets csv ...]   (default: data/targets_*.csv)
The generator 07_make_webpage.py embeds whatever is in data/cutouts/ as data URIs, so nothing external is loaded by the page.
"""
import io, os, sys, glob, time
import numpy as np
import pandas as pd
import requests
from concurrent.futures import ThreadPoolExecutor
from astropy.io import fits
from astropy.visualization import AsinhStretch, PercentileInterval, ImageNormalize
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, 'data')
CUT = os.path.join(DATA, 'cutouts')
SIZE_ARCSEC = 64
IBE_SEARCH = 'https://irsa.ipac.caltech.edu/ibe/search/ztf/products/ref'
IBE_DATA = 'https://irsa.ipac.caltech.edu/ibe/data/ztf/products/ref'


def sdss_jpeg(ra, dec, path):
    if os.path.exists(path):
        return True
    r = requests.get('https://skyserver.sdss.org/dr18/SkyServerWS/ImgCutout/getjpeg',
                     params={'ra': ra, 'dec': dec, 'scale': 0.4, 'width': 160, 'height': 160}, timeout=90)
    if r.ok and r.headers.get('content-type', '').startswith('image'):
        open(path, 'wb').write(r.content); return True
    return False


def ztf_refs(ra, dec):
    """IBE metadata rows for reference images covering the position; best row per filter = most frames."""
    q = requests.get(IBE_SEARCH, params={'POS': f'{ra},{dec}', 'ct': 'csv'}, timeout=180)
    lines = [l for l in q.text.splitlines() if l and not l.startswith('\\')]
    if len(lines) < 2:
        return {}
    df = pd.read_csv(io.StringIO('\n'.join(lines)))
    best = {}
    for fc, g in df.groupby('filtercode'):
        best[fc] = [row for _, row in g.sort_values('nframes', ascending=False).iterrows()]   # candidates, deepest first
    return best


def ztf_png(ra, dec, row, path):
    if os.path.exists(path):
        return True
    fld, cc, q, fc = int(row.field), int(row.ccdid), int(row.qid), row.filtercode
    f6 = f'{fld:06d}'
    url = f'{IBE_DATA}/{f6[:3]}/field{f6}/{fc}/ccd{cc:02d}/q{q}/ztf_{f6}_{fc}_c{cc:02d}_q{q}_refimg.fits'
    r = requests.get(url, params={'center': f'{ra},{dec}', 'size': f'{SIZE_ARCSEC}arcsec', 'gzip': 'false'}, timeout=240)
    if not r.ok:
        return False
    hdu = fits.open(io.BytesIO(r.content))[0]
    if hdu.data is None or hdu.data.ndim != 2 or min(hdu.data.shape) < 8:
        return False                      # position falls off the edge of this reference image; caller tries another field
    img = hdu.data.astype(float)
    # orient North up / East left from the WCS: with origin='lower', north-up needs CD2_2 > 0 and east-left CD1_1 < 0
    # (ZTF reference images are stored with CD1_1 > 0 and CD2_2 < 0, i.e. flipped both ways)
    h = hdu.header
    cd11 = h.get('CD1_1', -1.0); cd22 = h.get('CD2_2', 1.0)
    if cd22 < 0:
        img = img[::-1, :]
    if cd11 > 0:
        img = img[:, ::-1]
    norm = ImageNormalize(img, interval=PercentileInterval(99.3), stretch=AsinhStretch(0.1))
    fig = plt.figure(figsize=(1.6, 1.6), dpi=100); ax = fig.add_axes([0, 0, 1, 1]); ax.axis('off')
    ax.imshow(img, origin='lower', cmap='gray_r', norm=norm, interpolation='nearest')
    fig.savefig(path, dpi=100); plt.close(fig)
    return True


def one(row):
    name, ra, dec = row['name'], float(row['ra']), float(row['dec'])
    got = {}
    try:
        got['sdss'] = sdss_jpeg(ra, dec, os.path.join(CUT, f'{name}_sdss.jpg'))
    except Exception as e:
        got['sdss'] = f'ERR {str(e)[:40]}'
    need = [fc for fc in ('zg', 'zr') if not os.path.exists(os.path.join(CUT, f'{name}_ztf_{fc[1]}.png'))]
    if need:
        try:
            refs = ztf_refs(ra, dec)
            for fc in need:
                ok = False
                for row in refs.get(fc, [])[:3]:          # fall back to the next field if the cutout is empty at the edge
                    try:
                        ok = ztf_png(ra, dec, row, os.path.join(CUT, f'{name}_ztf_{fc[1]}.png'))
                    except Exception:
                        ok = False
                    if ok:
                        break
                got[fc] = ok if fc in refs else 'no ref image'
        except Exception as e:
            got['ztf'] = f'ERR {str(e)[:40]}'
    return name, got


def main():
    files = sys.argv[1:] or sorted(glob.glob(os.path.join(DATA, 'targets_*.csv')))
    t = pd.concat([pd.read_csv(f) for f in files]).drop_duplicates('name')
    os.makedirs(CUT, exist_ok=True)
    t0 = time.time(); bad = 0
    with ThreadPoolExecutor(max_workers=4) as ex:
        for k, (name, got) in enumerate(ex.map(one, [r for _, r in t.iterrows()])):
            if any(v is not True for v in got.values()):
                bad += 1; print(f'   {name}: {got}', flush=True)
            if (k + 1) % 20 == 0:
                print(f'[{time.time()-t0:5.0f}s] {k+1}/{len(t)}', flush=True)
    n = len(glob.glob(os.path.join(CUT, '*')))
    print(f'done: {len(t)} targets, {n} cutout files in data/cutouts, {bad} with problems, {time.time()-t0:.0f}s')


if __name__ == '__main__':
    main()

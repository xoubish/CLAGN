"""
06_finder_charts.py  --  Pan-STARRS1 finder charts + an NGPS-style target list for a night file.

Usage: /opt/anaconda3/bin/python 06_finder_charts.py data/targets_<night>.csv [size_arcmin]
Writes finders/<night>/<rank>_<name>.png (PS1 r-band 3' cutout, N up E left, target circle, 30" scale bar)
and finders/<night>/targetlist_<night>.txt  (name, RA/Dec sexagesimal, epoch, magnitude, priority, comment).
Only rows with rank > 0 (the primary list) get charts; backups go in the text list flagged 'backup'.
"""
import io, os, sys, time
import glob
import numpy as np
import pandas as pd
import requests
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.visualization import ZScaleInterval

HERE = os.path.dirname(os.path.abspath(__file__))


def ps1_cutout(ra, dec, size_pix, band='r', tries=3):
    """Return a PS1 stack cutout as a 2-D array (0.25 arcsec/pix) or None if outside PS1 / on failure."""
    for attempt in range(tries):
        try:
            tab = requests.get('https://ps1images.stsci.edu/cgi-bin/ps1filenames.py',
                               params={'ra': ra, 'dec': dec, 'filters': band, 'type': 'stack'}, timeout=120).text
            lines = [l for l in tab.splitlines() if l and not l.startswith('projcell')]
            if not lines:
                return None
            fname = lines[0].split()[7]
            url = ('https://ps1images.stsci.edu/cgi-bin/fitscut.cgi?ra=%f&dec=%f&size=%d&format=fits&red=%s' % (ra, dec, size_pix, fname))
            r = requests.get(url, timeout=180); r.raise_for_status()
            return fits.open(io.BytesIO(r.content))[0].data
        except Exception as e:
            print(f'   cutout attempt {attempt+1} failed: {str(e)[:80]}', flush=True); time.sleep(5)
    return None


def chart(row, outpng, size_arcmin=3.0):
    size_pix = int(size_arcmin * 60 / 0.25)
    img = ps1_cutout(row.ra, row.dec, size_pix)
    fig, ax = plt.subplots(figsize=(6, 6.4))
    if img is None:
        ax.text(0.5, 0.5, 'no PS1 image', ha='center', va='center', transform=ax.transAxes)
    else:
        img = np.nan_to_num(img)
        lo, hi = ZScaleInterval().get_limits(img)
        ax.imshow(img, origin='lower', cmap='gray_r', vmin=lo, vmax=hi)   # PS1 cutouts: N up, E left
    c = size_pix / 2
    ax.add_patch(plt.Circle((c, c), 5 / 0.25, fill=False, color='red', lw=1.5))
    # scale bar lower-left; compass lower-right (PS1 cutouts with origin='lower': North up, East LEFT)
    ax.plot([size_pix * 0.06, size_pix * 0.06 + 30 / 0.25], [size_pix * 0.06] * 2, color='red', lw=2)
    ax.text(size_pix * 0.06 + 15 / 0.25, size_pix * 0.085, '30"', color='red', ha='center')
    ox, oy = size_pix * 0.90, size_pix * 0.08
    ax.annotate('N', xy=(ox, oy), xytext=(ox, oy + size_pix * 0.12), ha='center', color='red',
                arrowprops=dict(arrowstyle='<-', color='red'))
    ax.annotate('E', xy=(ox, oy), xytext=(ox - size_pix * 0.12, oy), va='center', ha='right', color='red',
                arrowprops=dict(arrowstyle='<-', color='red'))
    sc = SkyCoord(row.ra * u.deg, row.dec * u.deg)
    ax.set_title(f"{row['name']}  {row.tier}  rank {int(row['rank'])}\n{sc.to_string('hmsdms', precision=2)}   z={row.z:.3f}  r={row.r_mag:.1f}",
                 fontsize=10)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_xlabel(f"{size_arcmin:.0f}' PS1 r   |   {row.get('trend', '')}   priority {row.priority:.2f}", fontsize=9)
    fig.tight_layout(); fig.savefig(outpng, dpi=110); plt.close(fig)


def main():
    inp = sys.argv[1]; size = float(sys.argv[2]) if len(sys.argv) > 2 else 3.0
    night = os.path.basename(inp).replace('targets_', '').replace('.csv', '')
    df = pd.read_csv(inp)
    outdir = os.path.join(HERE, 'finders', night); os.makedirs(outdir, exist_ok=True)
    sc = SkyCoord(df.ra.values * u.deg, df.dec.values * u.deg)
    with open(os.path.join(outdir, f'targetlist_{night}.txt'), 'w') as f:
        f.write('# name                   RA(J2000)     Dec(J2000)    z      r    tier  exposure    status  comment\n')
        for (_, row), c in zip(df.iterrows(), sc):
            status = 'primary' if row['rank'] > 0 else 'backup'
            f.write(f"{row['name']:22s} {c.ra.to_string(u.hour, sep=':', precision=2, pad=True):12s} "
                    f"{c.dec.to_string(u.deg, sep=':', precision=1, alwayssign=True, pad=True):12s} {row.z:5.3f} {row.r_mag:5.1f} "
                    f"{row.tier:4s} {str(row.get('exp_plan', '')):11s} {status:7s} {row.get('trend', '')}, last spec {row.years_since_last_spec:.1f} yr ago; "
                    f"{row.get('notes', '')}\n")
    prim = df[df['rank'] > 0]
    for _, row in prim.iterrows():
        outpng = os.path.join(outdir, f"{int(row['rank']):02d}_{row['name']}.png")
        if not os.path.exists(outpng):                       # same target under an old rank prefix: rename instead of refetching
            for old in glob.glob(os.path.join(outdir, f"[0-9][0-9]_{glob.escape(row['name'])}.png")):
                os.replace(old, outpng); break
        if not os.path.exists(outpng):
            chart(row, outpng, size)
            print(f'   {os.path.basename(outpng)}', flush=True)
    print(f'{len(prim)} charts + target list in {outdir}')


if __name__ == '__main__':
    main()

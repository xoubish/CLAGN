"""SDSS field of H5MUSES-327 with the Oct 26 NGPS slit, Gaia, WISE/Spitzer positions and IR beam sizes.

Run with /opt/anaconda3/bin/python. Read-only against SkyServer DR18 and IRSA.
Writes data/additional_targets/H5MUSES-327/field_slit.{png,pdf} and field_sources.csv.
"""
from pathlib import Path
import io

import numpy as np
import pandas as pd
import requests
from astropy.io import ascii, fits
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Polygon

ROOT = Path(__file__).resolve().parents[1]
DEST = ROOT / 'data/additional_targets/H5MUSES-327'
RA, DEC = 261.137, 59.44636111111111       # 5MUSES catalogue position (planning coordinates)
RAW = Path('/Users/shemmati/Desktop/Palomar_27/20261027/ngps_261027_0049.fits')
NEIGH = ROOT / 'data/reselection_2026-09-20/neighbour_cache/H5MUSES-327_gaia.csv'
# FWHM / aperture sizes in arcsec (WISE Explanatory Supplement; Spitzer IRAC/MIPS/IRS handbooks)
BEAMS = [('IRAC 3.6 µm', 1.66), ('MIPS 24 µm (5MUSES selection)', 5.9),
         ('WISE W1 / W2', 6.1), ('WISE W4 (22 µm)', 12.0)]


def offsets(ra, dec):
    """Tangent-plane offsets in arcsec, +x = east."""
    return (np.asarray(ra)-RA)*np.cos(np.radians(DEC))*3600, (np.asarray(dec)-DEC)*3600


def fetch_sdss(field, scale):
    n = int(round(field/scale))
    r = requests.get('https://skyserver.sdss.org/dr18/SkyServerWS/ImgCutout/getjpeg',
                     params=dict(ra=RA, dec=DEC, scale=scale, width=n, height=n), timeout=60)
    r.raise_for_status()
    return plt.imread(io.BytesIO(r.content), format='jpeg')


def seip():
    q = ("SELECT ra,dec,i1_f_ap1,i2_f_ap1,m1_f_psf FROM slphotdr4 WHERE "
         f"CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{RA},{DEC},0.0042))=1")
    r = requests.get('https://irsa.ipac.caltech.edu/TAP/sync', params=dict(QUERY=q, FORMAT='csv'), timeout=60)
    r.raise_for_status()
    return pd.read_csv(io.StringIO(r.text))


def main():
    h = fits.getheader(RAW)
    slit_pa, slit_w = float(h['SLITA']) % 360, float(h['SLITW'])
    slit_len = 50.                                      # NGPS slice length
    gaia = pd.read_csv(NEIGH)
    wise = ascii.read(DEST/'allwise.tbl', format='ipac').to_pandas().iloc[0]
    irac = seip()
    rows = [dict(source='5MUSES-327 catalogue', ra=RA, dec=DEC, note='planning coordinates')]
    rows += [dict(source=f'Gaia DR3 {g.source_id}', ra=g.ra, dec=g.dec, note=f'G={g.phot_g_mean_mag:.2f}')
             for g in gaia.itertuples()]
    rows.append(dict(source=f'AllWISE {wise.designation}', ra=wise.ra, dec=wise.dec,
                     note=f'W1={wise.w1mpro:.2f} W2={wise.w2mpro:.2f} W3={wise.w3mpro:.2f} W4={wise.w4mpro:.2f}'))
    rows += [dict(source='Spitzer SEIP', ra=s.ra, dec=s.dec,
                  note=f'3.6µm={s.i1_f_ap1:.0f} µJy, 4.5µm={s.i2_f_ap1:.0f} µJy'
                       + (f', 24µm={s.m1_f_psf/1000:.1f} mJy' if np.isfinite(s.m1_f_psf) else ''))
             for s in irac.itertuples()]
    table = pd.DataFrame(rows)
    x, y = offsets(table.ra, table.dec)
    table['sep_arcsec'], table['pa_deg'] = np.hypot(x, y), np.degrees(np.arctan2(x, y)) % 360
    # perpendicular / along-slit offsets relative to the slit centre (the target)
    pa = np.radians(slit_pa)
    table['along_slit'] = x*np.sin(pa)+y*np.cos(pa)
    table['across_slit'] = x*np.cos(pa)-y*np.sin(pa)
    table['in_slit'] = (np.abs(table.across_slit) <= slit_w/2) & (np.abs(table.along_slit) <= slit_len/2)
    table = table.sort_values('sep_arcsec')
    table.to_csv(DEST/'field_sources.csv', index=False, float_format='%.6f')
    print(table.to_string(index=False))

    is_gaia = table.source.str.startswith('Gaia')
    star = table[is_gaia].iloc[0]
    wpos, ipos = table[table.source.str.startswith('AllWISE')].iloc[0], table[table.source.eq('Spitzer SEIP')].iloc[0]

    fig, axes = plt.subplots(1, 2, figsize=(14, 7.2), constrained_layout=True)
    for ax, field, scale in [(axes[0], 80., 0.2), (axes[1], 20., 0.05)]:
        img = fetch_sdss(field, scale)
        half = field/2
        ax.imshow(img, extent=[half, -half, -half, half], origin='upper', interpolation='nearest')
        # slit (1.5" x 50" at the sky PA in SLITA)
        u = np.array([np.sin(pa), np.cos(pa)]); v = np.array([np.cos(pa), -np.sin(pa)])
        corners = [s*slit_len/2*u+t*slit_w/2*v for s, t in [(1, 1), (1, -1), (-1, -1), (-1, 1)]]
        ax.add_patch(Polygon(corners, closed=True, fill=False, ec='#ffd23f', lw=1.6,
                             label=f'NGPS slit {slit_w:.1f}″ × {slit_len:.0f}″, PA {slit_pa:.1f}°'))
        gx, gy = offsets(table[is_gaia].ra, table[is_gaia].dec)
        ax.scatter(gx, gy, s=260 if field < 30 else 160, facecolors='none', edgecolors='#3fd0ff', lw=1.6,
                   label='Gaia DR3 sources')
        ax.plot(0, 0, '+', color='#ff4f4f', ms=22 if field < 30 else 14, mew=1.8,
                label=f'H5MUSES-327 (5MUSES) {RA:.5f}, {DEC:+.5f}')
        wx, wy = offsets(wpos.ra, wpos.dec); ix, iy = offsets(ipos.ra, ipos.dec)
        ax.plot(wx, wy, 'x', color='#ff9cf0', ms=10, mew=1.6, label='AllWISE position')
        ax.plot(ix, iy, 'D', mfc='none', mec='#9cff6e', ms=8, mew=1.4, label='Spitzer IRAC/MIPS position')
        if field < 30:
            for (label, fwhm), ls in zip(BEAMS, [':', '--', '-.', '-']):
                ax.add_patch(Circle((ix, iy), fwhm/2, fill=False, ec='white', ls=ls, lw=1.1,
                                    label=f'{label} FWHM {fwhm}″'))
        ax.set_xlim(half, -half); ax.set_ylim(-half, half)
        ax.set_xlabel('ΔRA (arcsec, east left)'); ax.set_ylabel('ΔDec (arcsec)')
    others = table[is_gaia].iloc[1:]
    axes[0].set_title(f'SDSS DR18 gri, 80″ field. Nearest other Gaia source: '
                      f'{others.sep_arcsec.iloc[0]:.1f}″ (G={gaia.set_index(gaia.source_id.astype(str)).loc[others.source.iloc[0].split()[-1]].phot_g_mean_mag:.1f})',
                      fontsize=10)
    axes[1].set_title(f'20″ zoom: Gaia DR3 {star.source.split()[-1]} is {star.sep_arcsec:.2f}″ from the 5MUSES position\n'
                      f'AllWISE {wpos.sep_arcsec:.2f}″, Spitzer {ipos.sep_arcsec:.2f}″; circles = IR beam FWHM', fontsize=10)
    axes[0].legend(loc='lower left', fontsize=8, framealpha=.85)
    axes[1].legend(loc='lower left', fontsize=7.5, framealpha=.85)
    fig.suptitle('H5MUSES-327 — NGPS 2026-10-27 UT, slit at parallactic angle', fontsize=12)
    for suffix in ['png', 'pdf']:
        fig.savefig(DEST/f'field_slit.{suffix}', dpi=160)


if __name__ == '__main__':
    main()

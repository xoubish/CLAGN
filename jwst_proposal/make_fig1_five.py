"""Figure 1: three NGPS-observed targets with different recorded histories.

Columns: SDSS image; ZTF (g, r; nightly medians, alert-extended) with AllWISE+NEOWISE W1/W2
on a right axis; optical spectra (archival SDSS epochs + 2026 NGPS) joined to SPHEREx on one
log-wavelength axis. Epoch colours match between the light-curve ticks and the spectra.
Inputs: inputs/lightcurves/, inputs/archival_spectra/, sep23 NGPS P330E products,
inputs/spaxel_scryer_*.csv (SPHEREx, where available), data/cutouts/<P>_sdss.jpg.
"""
import json, os
from pathlib import Path
os.environ.setdefault('MPLCONFIGDIR', '/tmp/clagn-matplotlib')
import numpy as np, pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FixedLocator, FuncFormatter, NullFormatter
from astropy.time import Time
from scipy.ndimage import gaussian_filter1d

HERE = Path(__file__).resolve().parent; ROOT = HERE.parent
LC = HERE/'inputs/lightcurves'; ARCH = HERE/'inputs/archival_spectra'
NGPS = ROOT/'sep23_data/reduction_20260924/products_p330e/spectra'
IDS = ['F06', 'F04', 'R01']
SPHEREX = {'R01': HERE/'inputs/spaxel_scryer_P9694_ra11p5333_dec9p1225_cleaned.csv'}
EXTRACTED = HERE/'inputs/spherex_extracted'   # spherex_extract.py outputs, <id>.csv
W_ZERO = {'W1': 309.540, 'W2': 171.787}
green, pink, purple, w2c, orange, teal = '#298b71', '#b74680', '#7654a2', '#a85a16', '#e36a30', '#1f7a8c'
ARCH_COLORS = ['#2878a5', '#7654a2', '#4c9a6a']
C_AA = 2.99792458e18
# Vacuum rest wavelengths (um); same table as make_fig1_connected.py
LINES = [('Mg II', [.2799117]), ('[Ne V]', [.3426850]), ('[O II]', [.3727092]), ('[Ne III]', [.3869860]),
         ('Hγ', [.434168]), ('He II', [.4687020]), ('Hβ', [.486268]), ('[O III]', [.4960295, .5008240]),
         ('He I', [.5877290]), ('Hα', [.656461]), ('[S II]', [.671829, .673267]), ('[S III]', [.9071100, .9533200]),
         ('He I, Paγ', [1.0833, 1.0941]), ('Paβ', [1.2822]), ('Paα', [1.8756]), ('Brγ', [2.1661])]
GOLD, GOLD_TEXT = '#a68143', '#735625'
SELECTED_LINES = {'Mg II', 'Hβ', '[O III]', 'Hα', 'Paβ', 'Paα'}
HISTORY_LABELS = {'F06': 'Flare → fade', 'F04': 'IR decline', 'R01': 'Rise'}


def mark_lines(ax, z, fontsize=9.5):
    xlo, xhi = ax.get_xlim(); placed = []
    for label, rest in LINES:
        if label not in SELECTED_LINES: continue
        waves = [w*(1+z) for w in rest]
        if not all(xlo < w < xhi for w in waves): continue
        for w in waves: ax.axvline(w, ymax=.9, color=GOLD, lw=.5, ls=(0, (3, 3)), alpha=.65, zorder=1)
        x = float(np.mean(waves))
        # push the label right if it would overlap the previous one (log spacing)
        if placed and np.log10(x/placed[-1]) < .033: x = placed[-1]*10**.033
        placed.append(x)
        ax.text(x, .975, label, rotation=90, va='top', ha='center', fontsize=fontsize, color=GOLD_TEXT,
                transform=ax.get_xaxis_transform(), bbox=dict(facecolor='white', edgecolor='none', pad=.15, alpha=.85), zorder=5)
        if abs(np.log10(x/np.mean(waves))) > 1e-6:
            ax.plot([x, np.mean(waves)], [.74, .70], transform=ax.get_xaxis_transform(), color=GOLD, lw=.5, zorder=1)
year = lambda m: Time(np.asarray(m, float), format='mjd').decimalyear
plt.rcParams.update({'font.family': 'serif', 'font.serif': ['Times New Roman'], 'font.size': 11,
                     'axes.labelsize': 11, 'xtick.labelsize': 10, 'ytick.labelsize': 10, 'pdf.fonttype': 42,
                     'axes.spines.top': False, 'axes.spines.right': False, 'axes.linewidth': .5})


def flam_to_mjy(wave_A, flux_1e17):
    return flux_1e17*1e-17*wave_A**2/C_AA/1e-26


def ztf(tid):
    dr = pd.read_csv(LC/f'{tid}_ztf_irsa.csv'); dr = dr[dr.catflags.eq(0)]
    al = pd.read_csv(LC/f'{tid}_alerce.csv').drop_duplicates('candid')
    al = al[al.corrected.astype(bool) & ~al.dubious.astype(bool) & np.isfinite(al.magpsf_corr)]
    out = {}
    for band, code, fid in [('g', 'zg', 1), ('r', 'zr', 2)]:
        d = dr[dr.filtercode.eq(code)]
        nights = d.groupby(np.floor(d.mjd)).agg(mjd=('mjd', 'median'), mag=('mag', 'median')).to_numpy()
        a = al[al.fid.eq(fid)]
        near = np.abs(d.mjd.values[None, :]-a.mjd.values[:, None]).argmin(1) if len(d) and len(a) else np.array([], int)
        matched = np.abs(d.mjd.values[near]-a.mjd.values) < .02 if len(near) else np.array([], bool)
        off = float(np.median(a.magpsf_corr.values[matched]-d.mag.values[near][matched])) if matched.sum() >= 5 else 0.
        post = a[a.mjd > d.mjd.max()] if len(d) else a
        post = post.groupby(np.floor(post.mjd)).agg(mjd=('mjd', 'median'), mag=('magpsf_corr', 'median')).to_numpy()
        if len(post): post[:, 1] -= off
        out[band] = (nights, post)
    return out


def wise(tid, internal):
    mep = pd.read_csv(LC/f'{tid}_allwise_mep.csv', dtype={'cc_flags': str, 'moon_masked': str})
    mep = mep[(mep.qi_fact > 0) & mep.cc_flags.str.zfill(4).str[:2].eq('00')]
    neo = pd.read_csv(HERE/'review/history_audit/binned_measurements.csv')
    neo = neo[neo.id.eq(tid)]
    out = {}
    for band in ['W1', 'W2']:
        col = band.lower()+'mpro_ep'
        g = mep.dropna(subset=[col]); g = g.assign(visit=np.round(g.mjd/180))
        early = [(float(v.mjd.median()), float(W_ZERO[band]*1e3*10**(-.4*v[col].median()))) for _, v in g.groupby('visit') if len(v) >= 3]
        n = neo[neo.band.eq(band)].sort_values('mjd')
        late = [(float(m), float(W_ZERO[band]*1e3*10**(-.4*mag))) for m, mag in zip(n.mjd, n.mag)]
        out[band] = np.array(early+late)
    return out


def spectra(tid, z):
    idx = pd.read_csv(ARCH/'index.csv'); idx = idx[idx.id.eq(tid)].sort_values('mjd').drop_duplicates('mjd')
    # earliest, latest public, and one in between if the record is long
    picks = [idx.iloc[0]]
    yr = lambda e: Time(e.mjd, format='mjd').datetime.year
    if len(idx) > 2:
        mid = idx.iloc[len(idx)//2]
        if yr(mid) not in (yr(idx.iloc[0]), yr(idx.iloc[-1])): picks.append(mid)
    if len(idx) > 1 and yr(idx.iloc[-1]) != yr(idx.iloc[0]): picks.append(idx.iloc[-1])
    epochs = []
    for k, e in enumerate(picks):
        d = pd.read_csv(HERE/e.file); ok = (d.ivar > 0) & np.isfinite(d.flux)
        w, f = d.wave_A.values[ok], d.flux.values[ok]
        f = gaussian_filter1d(f, 2.)
        epochs.append(dict(label=Time(e.mjd, format='mjd').datetime.strftime('%Y'), mjd=float(e.mjd),
                           wave_um=w/1e4, mjy=flam_to_mjy(w, f), color=ARCH_COLORS[k % 3], ngps=False))
    internal = pd.read_csv(HERE/'inputs/jwst_sample_cycle6.csv').set_index('id').loc[tid, 'internal_id']
    n = pd.read_csv(NGPS/f'{internal}.csv')
    w, f = n.WAVE_VAC_HELIO_A.values, np.where(n.MASK.values, n.FLUX.values, np.nan)
    good = np.isfinite(f)
    fs = f.copy(); fs[good] = gaussian_filter1d(f[good], 4.)
    from astropy.io import fits
    mjd = fits.getheader(NGPS/f'{internal}.fits', 1).get('MJD', 61307.)
    epochs.append(dict(label='2026 NGPS', mjd=float(mjd), wave_um=w/1e4, mjy=flam_to_mjy(w, fs), color=orange, ngps=True))
    return epochs


def spherex(tid):
    if tid in SPHEREX and SPHEREX[tid].exists():
        d = pd.read_csv(SPHEREX[tid])
        return d[d.spectrum_type.eq('cleaned')] if 'spectrum_type' in d else d
    if (EXTRACTED/f'{tid}.csv').exists():
        d = pd.read_csv(EXTRACTED/f'{tid}.csv')
        return d[np.isfinite(d.flux_mjy) & (d.flux_err_mjy > 0) & (d.flux_mjy > 3*d.flux_err_mjy - 1e9)]
    return None


def main(layout='stack'):
    sample = pd.read_csv(HERE/'inputs/jwst_sample_cycle6.csv').set_index('id')
    n = len(IDS)
    stack = layout == 'stack'
    fig = plt.figure(figsize=((9.0 if stack else 7.2), (1.86*n+.4) if stack else (1.3*n+.55)), facecolor='white')
    left, top, gap = .045, .975 if stack else .96, .024 if stack else .015
    rowh = (top-(.06 if stack else .085))/n
    inventory = []
    for i, tid in enumerate(IDS):
        r = sample.loc[tid]; z = float(r.z)
        y0 = top-(i+1)*rowh+gap
        h = rowh-2*gap
        if stack:
            usable = h-.055                # light curves : spectrum = 1 : 1.6
            lch, sph = usable/2.6, 1.6*usable/2.6
            imh = .12*fig.get_figwidth()/fig.get_figheight()
            ax_im = fig.add_axes([.015, y0+h-imh, .12, imh])
            x0 = .20
            wl = (.99-x0-.065)/2.
            ax_lc = fig.add_axes([x0, y0+sph+.055, wl, lch])
            ax_w = fig.add_axes([x0+wl+.065, y0+sph+.055, wl, lch])
            ax_sp = fig.add_axes([x0, y0, .99-x0, sph])
            # light box around the block: image + light curves + spectrum
            from matplotlib.patches import FancyBboxPatch
            fig.patches.append(FancyBboxPatch((.004, y0-.028), .992, h+.046, boxstyle='round,pad=0,rounding_size=.006',
                                              transform=fig.transFigure, facecolor='#f6f7f8', edgecolor='#c9cfd4', lw=.6, zorder=-5))
        else:
            ax_im = fig.add_axes([.012, y0, h*fig.get_figheight()/fig.get_figwidth(), h])
            ax_lc = fig.add_axes([.225, y0, .155, h])
            ax_w = fig.add_axes([.43, y0, .155, h])
            ax_sp = fig.add_axes([.635, y0, .355, h])
        # image: 40" field, show the central 24"
        hires = HERE/f'inputs/cutouts_hires/{tid}_sdss_24arcsec.jpg'
        if hires.exists():
            ax_im.imshow(plt.imread(hires), extent=[-12, 12, -12, 12])
        else:
            ax_im.imshow(plt.imread(ROOT/f'data/cutouts/{r.internal_id}_sdss.jpg'), extent=[-20, 20, -20, 20])
        ax_im.set(xlim=(-12, 12), ylim=(-12, 12)); ax_im.set_axis_off()
        ax_im.plot([5.5, 10.5], [-9.8, -9.8], color='white', lw=1.1); ax_im.text(8, -8.9, '5″', color='white', ha='center', fontsize=10)
        if stack:
            fig.text(.022,y0+.018,f'{tid}   z = {z:.2f}\n{HISTORY_LABELS[tid]}',
                     fontsize=10.5,ha='left',va='bottom',linespacing=1.3)
        else:
            ax_im.text(-11,10.6,f'{tid}\nz = {z:.2f}',color='white',fontsize=9,va='top')
        # light curve
        lc = ztf(tid)
        for band, c in [('g', green), ('r', pink)]:
            nights, post = lc[band]
            ax_lc.plot(year(nights[:, 0]), nights[:, 1], '.', ms=1.2, color=c, alpha=.6, rasterized=True)
            if len(post): ax_lc.plot(year(post[:, 0]), post[:, 1], 'o', ms=1.8, mfc='white', mec=c, mew=.4)
        ax_lc.invert_yaxis(); ax_lc.set(xlim=(2017.4, 2027.4), xticks=[2019, 2022, 2025])
        ax_lc.tick_params(pad=1, length=2); ax_lc.yaxis.set_major_locator(plt.MaxNLocator(3, prune='both'))
        ax_w.set(xlim=(2009.6, 2027.4), xticks=[2010, 2015, 2020, 2025])
        ws = wise(tid, r.internal_id)
        allw = np.concatenate([ws[b][:, 1] for b in ['W1', 'W2'] if len(ws[b])])
        ax_w.set_ylim(allw.min()*.93, allw.max()*1.07)
        for band, c, mk in [('W1', purple, 'o'), ('W2', w2c, 's')]:
            a = ws[band]
            if len(a):
                a = a[np.argsort(a[:, 0])]
                ax_w.plot(year(a[:, 0]), a[:, 1], mk, ms=2.6, mfc='white', mec=c, mew=.65, zorder=2)
        ax_w.tick_params(pad=1, length=2); ax_w.yaxis.set_major_locator(plt.MaxNLocator(3, prune='both'))
        if i == 0 and not stack:
            ax_lc.set_title('ZTF', loc='left', fontsize=8.5, pad=2)
            ax_lc.text(.03, .06, 'g', color=green, transform=ax_lc.transAxes, fontsize=8, va='bottom', fontweight='bold')
            ax_lc.text(.12, .06, 'r', color=pink, transform=ax_lc.transAxes, fontsize=8, va='bottom', fontweight='bold')
            ax_w.set_title('WISE', loc='left', fontsize=8.5, pad=2)
        if i == 0 and stack:
            ax_lc.text(.03, .06, 'g', color=green, transform=ax_lc.transAxes, fontsize=10, va='bottom', fontweight='bold')
            ax_lc.text(.12, .06, 'r', color=pink, transform=ax_lc.transAxes, fontsize=10, va='bottom', fontweight='bold')
            ax_w.text(.03, .94, 'W1', color=purple, transform=ax_w.transAxes, fontsize=10, va='top')
            ax_w.text(.18, .94, 'W2', color=w2c, transform=ax_w.transAxes, fontsize=10, va='top')
            ax_sp.set_title('Spectra' if not stack else '', loc='left', fontsize=8.5, pad=2)
        if stack:
            pass                          # year tick labels on every block, no axis title
        elif i == n-1:
            ax_lc.set_xlabel('Year', labelpad=1); ax_w.set_xlabel('Year', labelpad=1)
        else:
            ax_lc.set_xticklabels([]); ax_w.set_xticklabels([])
        # spectra: optical archival + NGPS + SPHEREx on one log-wavelength axis
        eps = spectra(tid, z)
        for e in eps:
            ax_sp.plot(e['wave_um'], e['mjy'], color=e['color'], lw=.9 if e['ngps'] else .7, zorder=3 if e['ngps'] else 2)
            for ax in (ax_lc, ax_w): ax.axvline(float(year(e['mjd'])), color=e['color'], lw=.8, ls='--', alpha=.85, zorder=.5)
        sp = spherex(tid)
        fit_json = HERE/f'inputs/spherex_fits/{tid}.json'
        if sp is not None:
            ax_sp.errorbar(sp.wavelength_um, sp.flux_mjy, yerr=sp.flux_err_mjy, fmt='o', ms=1.6, lw=.3, color=teal, mec='white', mew=.2, alpha=.9, zorder=2.5, rasterized=True)
            if fit_json.exists():
                fit = json.loads(fit_json.read_text()); curve = pd.read_csv(HERE/f'inputs/spherex_fits/{tid}_curve.csv')
                ax_sp.plot(curve.wavelength_observed_um, curve.total_mjy, color='#242424', lw=1.0, ls=(0, (4, 2.5)), zorder=2.7)
                fit_label = 'Mean continuum'
            mjds = np.sort(sp.mjds.unique())
            for m0 in [mjds.min(), mjds.max()]:
                for ax in (ax_lc, ax_w): ax.axvline(float(year(m0)), color=teal, lw=.8, ls=':', zorder=.5)
        else:
            ax_sp.text(2.2, .5, 'SPHEREx: extraction in progress', color=teal, fontsize=7.5, ha='center', transform=ax_sp.get_xaxis_transform(), va='center')
        ax_sp.set(xscale='log', yscale='log', xlim=(.34, 5.1))
        allf = np.concatenate([e['mjy'][np.isfinite(e['mjy'])] for e in eps]+([sp.flux_mjy.values] if sp is not None else []))
        lo, hi = np.nanpercentile(allf, [1, 99.7]); ax_sp.set_ylim(lo/1.6, hi*2.2)
        ax_sp.xaxis.set_major_locator(FixedLocator([.4, .6, 1, 2, 3, 5])); ax_sp.xaxis.set_major_formatter(FuncFormatter(lambda v, p: f'{v:g}'))
        ax_sp.xaxis.set_minor_formatter(NullFormatter())
        ax_sp.yaxis.set_major_locator(FixedLocator([.1, .5, 1, 5, 10, 50])); ax_sp.yaxis.set_major_formatter(FuncFormatter(lambda v, p: f'{v:g}'))
        ax_sp.yaxis.set_minor_formatter(NullFormatter()); ax_sp.tick_params(pad=1, length=2, which='both')
        ax_sp.tick_params(axis='y', which='minor', length=0)
        ax_sp.axvline(.98, color='#c8ccd0', lw=.5)
        if i == n-1: ax_sp.set_xlabel('Observed wavelength (µm)', labelpad=1)
        ax_sp.set_ylabel('mJy', labelpad=1)
        handles = [Line2D([], [], color=e['color'], lw=1.2, label=e['label']) for e in eps]
        if sp is not None:
            handles.append(Line2D([], [], marker='o', color=teal, ls='', ms=2.5, label='SPHEREx 2025–26'))
            if fit_json.exists(): handles.append(Line2D([], [], color='#242424', lw=1.0, ls=(0, (4, 2.5)), label=fit_label))
        ax_sp.legend(handles=handles, loc='lower right', bbox_to_anchor=(1.0, .0), fontsize=9.5, frameon=True, facecolor='white',
                     edgecolor='none', framealpha=.85, ncol=3, handlelength=1.1, handletextpad=.3, columnspacing=.7, borderaxespad=.2, borderpad=.2)
        mark_lines(ax_sp, z)
        ax_lc.set_ylabel('ZTF mag' if stack else 'mag', labelpad=1, fontsize=10); ax_w.set_ylabel('WISE mJy' if stack else 'mJy', labelpad=1, fontsize=10)
        inventory.append(dict(id=tid,optical_epochs=[e['label'] for e in eps],
                              spherex_points=len(sp) if sp is not None else 0,
                              expected_lines=sorted(SELECTED_LINES)))
    out = HERE/('fig1_targets.pdf' if stack else 'older/fig1_row_preview.pdf')
    fig.savefig(out, dpi=240, bbox_inches='tight', pad_inches=.02)
    if stack:
        (HERE/'inputs/fig1_targets_provenance.json').write_text(json.dumps(dict(
            generator='make_fig1_five.py stack',targets=inventory,
            presentation='Three examples with smaller cutouts, enlarged axes and legends, six selected line identifications, and unchanged epoch selection and flux normalization.'),indent=2)+'\n')
    print(out.name, 'written')


if __name__ == '__main__':
    import sys
    main(sys.argv[1] if len(sys.argv) > 1 else 'stack')

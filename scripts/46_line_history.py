"""Emission-line flux and continuum-slope history of one target from the cached spectra.

For every night with a spectrum, fit a local linear continuum plus one Gaussian to the chosen
line in the observed frame, and a straight line to the continuum in an observed-frame window
with the emission line masked. Noise is estimated per spectrum from the residual scatter, so
the error bars are internal fit uncertainties, not a full spectrophotometric error budget.
Absolute fluxes depend on the fibre aperture of each survey (SDSS legacy 3", BOSS/eBOSS and
SDSS-V 2", DESI 1.5") and on each survey's flux calibration; equivalent width and the
normalised slope are the aperture-robust quantities.

Usage: python scripts/46_line_history.py P1793 [--line MgII] [--slope 4000 5000]
Writes observing/sep23/figures/<name>_<line>_history.{png,csv} and a fit-inspection sheet.
"""
import argparse
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter, MaxNLocator
from scipy.optimize import curve_fit
from scipy.stats import pearsonr

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'data/reselection_2026-09-20'
DEST = ROOT/'observing/sep23/figures'
C_KMS = 299792.458
ZTF_G_PIVOT = 4770.  # Angstrom; ZTF g effective wavelength
ZTF_R_PIVOT = 6400.  # Angstrom; ZTF r effective wavelength; spectroscopic continuum taken as the median in 6350-6450 A
# rest wavelength, fit window, continuum side windows (rest Angstrom)
LINES = {
    'MgII': dict(rest=2798.75, window=(2650, 2960), sides=[(2650, 2720), (2900, 2960)], label='Mg II 2798'),
    'Hbeta': dict(rest=4862.68, window=(4700, 5130), sides=[(4700, 4760), (5080, 5130)], label='Hβ 4861'),
    'Halpha': dict(rest=6564.61, window=(6350, 6800), sides=[(6350, 6450), (6720, 6800)], label='Hα 6563'),
}
APERTURE = {'legacy': 'SDSS 3" fibre', 'eboss': 'eBOSS 2" fibre', 'AQMES-Medium': 'SDSS-V 2" fibre', 'bhm_aqmes': 'SDSS-V 2" fibre', 'dark': 'DESI 1.5" fibre', 'bright': 'DESI 1.5" fibre'}


def load_epochs(name):
    """One spectrum per night: per-night SDSS reductions and DESI coadds, with SDSS-V field quality."""
    public = json.loads((ROOT/'data/spectra_dl'/f'{name}.json').read_text())
    internal_path = OUT/'sdssv_spectra'/f'{name}.json'
    internal = json.loads(internal_path.read_text()) if internal_path.exists() else []
    quality = {int(round(r['mjd'])): r.get('fieldquality') for r in internal}
    from spectral_utils import selected_records
    epochs = []
    for r in selected_records(public+internal):
        if r.get('mjd') is None:
            continue
        key = int(np.floor(r['mjd']))
        meta = r.get('meta', {})
        program = r.get('program') or meta.get('programname') or ''
        epochs.append(dict(mjd=float(r['mjd']), night=key, source=r.get('source'), program=program, aperture=APERTURE.get(program, 'unknown'),
                           fieldquality=quality.get(key), proprietary=bool(r.get('proprietary')), zwarning=meta.get('zwarning'),
                           pipeline_class=meta.get('class'), pipeline_z=meta.get('z'), sn_median=r.get('sn_median_all') or meta.get('sn_median_all'),
                           wave=np.asarray(r['wave'], float), flux=np.asarray(r['flux'], float)))
    return epochs


def model(w, a, b, amp, mu, sigma):
    return a+b*(w-mu)+amp*np.exp(-0.5*((w-mu)/sigma)**2)


def fit_line(w, f, z, line):
    """Single-Gaussian diagnostic, not a broad/narrow decomposition or state classification."""
    L = LINES[line]
    lo, hi = [v*(1+z) for v in L['window']]
    sel = (w >= lo) & (w <= hi) & np.isfinite(f)
    if line == 'Hbeta':
        for a, b in [(4945, 4973), (4995, 5022)]:
            sel &= ~((w >= a*(1+z)) & (w <= b*(1+z)))
    x, y = w[sel], f[sel]
    if sel.sum() < 15:
        return None
    mu0 = L['rest']*(1+z)
    side = np.zeros_like(x, bool)
    for a, b in L['sides']:
        side |= (x >= a*(1+z)) & (x <= b*(1+z))
    if any(((x >= a*(1+z)) & (x <= b*(1+z))).sum() < 3 for a, b in L['sides']):
        return None
    cont = np.polyfit(x[side]-mu0, y[side], 1)
    resid = y-np.polyval(cont, x-mu0)
    amp0 = max(resid[np.abs(x-mu0) < 40*(1+z)].max(), 1e-3)
    p0 = [cont[1], cont[0], amp0, mu0, 25*(1+z)]
    bounds = ([-np.inf, -np.inf, 0, mu0-25*(1+z), 4*(1+z)], [np.inf, np.inf, np.inf, mu0+25*(1+z), 70*(1+z)])
    try:
        p, _ = curve_fit(model, x, y, p0=p0, bounds=bounds, maxfev=20000)
        noise = 1.4826*np.median(np.abs(y-model(x, *p)))
        noise = max(noise, 1e-6)
        p, cov = curve_fit(model, x, y, p0=p, sigma=np.full_like(y, noise), absolute_sigma=True, bounds=bounds, maxfev=20000)
    except (RuntimeError, ValueError):
        return None
    err = np.sqrt(np.diag(cov))
    a, b, amp, mu, sigma = p
    flux = amp*sigma*np.sqrt(2*np.pi)
    flux_gradient = np.array([0., 0., sigma*np.sqrt(2*np.pi), 0., amp*np.sqrt(2*np.pi)])
    flux_err = np.sqrt(max(0., flux_gradient @ cov @ flux_gradient))
    width_gradient = np.array([0., 0., 0., -2.3548*sigma/mu**2*C_KMS, 2.3548/mu*C_KMS])
    width_err = np.sqrt(max(0., width_gradient @ cov @ width_gradient))
    cont_at_line = a
    return dict(flux=flux, flux_err=flux_err, ew_rest=flux/cont_at_line/(1+z) if cont_at_line > 0 else np.nan,
                fwhm_kms=2.3548*sigma/mu*C_KMS, fwhm_err_kms=width_err, centre_kms=(mu/mu0-1)*C_KMS,
                cont_at_line=cont_at_line, noise=noise, chi2_red=float(np.sum(((y-model(x, *p))/noise)**2)/(len(y)-5)),
                params=p, x=x, y=y, mask_lo=mu-3*sigma, mask_hi=mu+3*sigma)


def fit_slope(w, f, lo, hi, masks, noise):
    sel = (w >= lo) & (w <= hi) & np.isfinite(f)
    for a, b in masks:
        sel &= ~((w >= a) & (w <= b))
    x, y = w[sel], f[sel]
    if sel.sum() < 10:
        return None
    mid = 0.5*(lo+hi)
    p, cov = np.polyfit(x-mid, y, 1, cov=True)
    err = np.sqrt(np.diag(cov))
    ok = y > 0
    alpha, alpha_err = (np.nan, np.nan)
    if ok.sum() > 10:
        q, qcov = np.polyfit(np.log10(x[ok]), np.log10(y[ok]), 1, cov=True)
        alpha, alpha_err = q[0], np.sqrt(qcov[0, 0])
    return dict(slope_per_1000A=p[0]*1000, slope_err_per_1000A=err[0]*1000, level_mid=p[1], level_mid_err=err[1],
                norm_slope_per_1000A=p[0]*1000/p[1] if p[1] > 0 else np.nan, alpha_lambda=alpha, alpha_lambda_err=alpha_err, n_pix=int(sel.sum()), mid=mid)


def year(mjd):
    return 2000+(mjd-51544.5)/365.25


def load_ztf(name):
    """ZTF g and r light curves [[mjd, mag, err], ...]: full-cadence DR24 cache first, page payload as fallback."""
    path = ROOT/'data/ztf_cache/review_dr24_20260920'/f'{name}.csv'
    if path.exists():
        d = pd.read_csv(path)
        cols = {c.lower(): c for c in d.columns}
        tcol = next((cols[c] for c in ['mjd', 'hjd', 'mjd_obs', 'obsmjd'] if c in cols), None)
        mcol = next((cols[c] for c in ['mag', 'magpsf'] if c in cols), None)
        ecol = next((cols[c] for c in ['magerr', 'sigmag', 'sigmapsf', 'mag_err'] if c in cols), None)
        bcol = next((cols[c] for c in ['filtercode', 'band', 'filter', 'fid'] if c in cols), None)
        if tcol and mcol and bcol:
            out = {}
            for band, sub in d.groupby(bcol):
                key = str(band).lower().replace('zg', 'g').replace('zr', 'r').replace('zi', 'i')
                key = {'1': 'g', '2': 'r', '3': 'i'}.get(key, key)
                mjd = sub[tcol].to_numpy(float)
                mjd = mjd-2400000.5 if mjd.min() > 2400000 else mjd
                err = sub[ecol].to_numpy(float) if ecol else np.full(len(sub), np.nan)
                out[key] = np.column_stack([mjd, sub[mcol].to_numpy(float), err])
            if out:
                return out
    text = (OUT/'candidate_review_local.html').read_text()
    data = json.loads(re.search(r'<script id="candidate-data" type="application/json">(.*?)</script>', text, re.S).group(1))
    for target in data['targets']:
        if target['name'] == name:
            return {band: np.asarray(v, float) for band, v in (target.get('ztf') or {}).items() if v}
    return {}


def mag_to_flam(mag, pivot):
    """AB magnitude to F_lambda in 1e-17 erg s^-1 cm^-2 A^-1 at the pivot wavelength."""
    fnu = 10**(-0.4*(np.asarray(mag, float)+48.6))
    return fnu*2.99792458e18/pivot**2/1e-17


def matched_ztf(ztf, mjd, band='g', window=20.):
    if band not in ztf:
        return np.nan, 0
    v = ztf[band]
    sel = np.abs(v[:, 0]-mjd) <= window
    return (float(np.median(v[sel, 1])) if sel.any() else np.nan), int(sel.sum())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('name')
    ap.add_argument('--line', default='MgII', choices=list(LINES))
    ap.add_argument('--slope', nargs=2, type=float, default=[4000, 5000], metavar=('LO', 'HI'), help='observed-frame slope window in Angstrom')
    args = ap.parse_args()
    DEST.mkdir(parents=True, exist_ok=True)
    targets = pd.read_csv(OUT/'compact_review_objects.csv').set_index('name')
    z = float(targets.loc[args.name].z)
    epochs = load_epochs(args.name)
    ztf = load_ztf(args.name)
    rows = []
    fits = []
    for e in epochs:
        ok = np.isfinite(e['flux']) & (e['wave'] > 3650)
        w, f = e['wave'][ok], e['flux'][ok]
        line = fit_line(w, f, z, args.line)
        masks = [(line['mask_lo'], line['mask_hi'])] if line else [(LINES[args.line]['rest']*(1+z)-150, LINES[args.line]['rest']*(1+z)+150)]
        slope = fit_slope(w, f, args.slope[0], args.slope[1], masks, line['noise'] if line else None)
        row = dict(name=args.name, night_mjd=e['night'], mjd=e['mjd'], year=year(e['mjd']), source=e['source'], program=e['program'], aperture=e['aperture'],
                   fieldquality=e['fieldquality'], proprietary=e['proprietary'], pipeline_class=e['pipeline_class'], pipeline_z=e['pipeline_z'], sn_median=e['sn_median'],
                   suspect=(e['fieldquality'] == 'bad') or (e['pipeline_class'] not in (None, 'QSO')))
        if line:
            row.update({k: line[k] for k in ['flux', 'flux_err', 'ew_rest', 'fwhm_kms', 'fwhm_err_kms', 'centre_kms', 'cont_at_line', 'noise', 'chi2_red']})
        if slope:
            row.update(slope)
            row['cont_ztfg_pivot'] = slope['level_mid']+slope['slope_per_1000A']*(ZTF_G_PIVOT-slope['mid'])/1000
        gmag, n_g = matched_ztf(ztf, e['mjd'], 'g')
        rmag, n_r = matched_ztf(ztf, e['mjd'], 'r')
        red = (w >= ZTF_R_PIVOT-50) & (w <= ZTF_R_PIVOT+50)
        row.update(ztf_g_mag=gmag, ztf_g_n=n_g, ztf_g_flam=mag_to_flam(gmag, ZTF_G_PIVOT) if np.isfinite(gmag) else np.nan,
                   ztf_r_mag=rmag, ztf_r_n=n_r, ztf_r_flam=mag_to_flam(rmag, ZTF_R_PIVOT) if np.isfinite(rmag) else np.nan,
                   cont_ztfr_pivot=float(np.median(f[red])) if red.sum() > 3 else np.nan)
        rows.append(row)
        fits.append((e, line))
    table = pd.DataFrame(rows)
    stem = DEST/f'{args.name}_{args.line}_history'
    table.drop(columns=[c for c in ['params', 'x', 'y'] if c in table]).to_csv(stem.with_suffix('.csv'), index=False)

    # History figure: broken time axis so the dense recent years stay readable.
    good = table[~table.suspect]
    bad = table[table.suspect]
    gaps = np.diff(np.r_[table.year.min(), table.year.values, table.year.max()])
    segments = []
    seg_start = table.year.min()
    ys = sorted(table.year.values)
    for a, b in zip(ys, ys[1:]):
        if b-a > 2.5:
            segments.append((seg_start-0.4, a+0.4))
            seg_start = b
    segments.append((seg_start-0.4, table.year.max()+0.4))
    widths = [max(0.6, b-a) for a, b in segments]
    panels = [('flux', 'flux_err', f"{LINES[args.line]['label']} flux\n(10⁻¹⁷ erg s⁻¹ cm⁻²)"),
              ('ew_rest', None, 'rest EW (Å)'),
              ('slope_per_1000A', 'slope_err_per_1000A', f'continuum slope {int(args.slope[0])}–{int(args.slope[1])} Å\n(10⁻¹⁷ erg s⁻¹ cm⁻² Å⁻¹ per 1000 Å)'),
              ('cont_ztfg_pivot', 'level_mid_err', f'continuum at {int(ZTF_G_PIVOT)} Å\n(10⁻¹⁷ erg s⁻¹ cm⁻² Å⁻¹) · grey: ZTF g')]
    fig, axes = plt.subplots(len(panels), len(segments), figsize=(12, 11), sharey='row', gridspec_kw=dict(width_ratios=widths, wspace=0.06, hspace=0.12))
    axes = np.atleast_2d(axes)
    markers = {'legacy': 'o', 'eboss': 's', 'AQMES-Medium': 'D', 'bhm_aqmes': 'D', 'dark': '^', 'bright': '^'}
    colors = {'SDSS 3" fibre': '#1f77b4', 'eBOSS 2" fibre': '#ff7f0e', 'SDSS-V 2" fibre': '#2ca02c', 'DESI 1.5" fibre': '#9467bd', 'unknown': 'grey'}
    for i, (col, ecol, label) in enumerate(panels):
        for j, (a, b) in enumerate(segments):
            ax = axes[i, j]
            for frame, filled in [(good, True), (bad, False)]:
                for ap_label, sub in frame.groupby('aperture'):
                    if col not in sub:
                        continue
                    ax.errorbar(sub.year, sub[col], yerr=sub[ecol] if ecol and ecol in sub else None, fmt=markers.get(sub.program.iloc[0], 'o'),
                                mfc=colors.get(ap_label, 'grey') if filled else 'none', mec=colors.get(ap_label, 'grey'), ecolor=colors.get(ap_label, 'grey'),
                                ms=6, capsize=2, lw=1, label=(ap_label if (i == 0 and j == len(segments)-1 and filled) else None))
            if col == 'cont_ztfg_pivot' and 'g' in ztf:
                g = ztf['g']
                ax.plot(year(g[:, 0]), mag_to_flam(g[:, 1], ZTF_G_PIVOT), '.', color='0.6', ms=3, zorder=0, label=('ZTF g photometry' if j == len(segments)-1 else None))
            ax.set_xlim(a, b)
            ax.xaxis.set_major_locator(MaxNLocator(nbins=3 if b-a < 3 else 6))
            ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f'{v:.1f}' if (b-a) < 3 else f'{v:.0f}'))
            ax.grid(alpha=.25)
            if j == 0:
                ax.set_ylabel(label, fontsize=9)
            else:
                ax.spines['left'].set_visible(False)
                ax.tick_params(labelleft=False, left=False)
            if j < len(segments)-1:
                ax.spines['right'].set_visible(False)
            if i < len(panels)-1:
                ax.tick_params(labelbottom=False)
    handles, labels = axes[0, -1].get_legend_handles_labels()
    h2, l2 = axes[-1, -1].get_legend_handles_labels()
    if handles:
        axes[0, -1].legend(handles+h2, labels+l2, fontsize=8, loc='best')
    axes[-1, len(segments)//2].set_xlabel('year')
    fig.suptitle(f"{args.name}  z={z:.3f} · {LINES[args.line]['label']} Gaussian fit and observed-frame continuum slope · open symbols: bad field quality or non-QSO pipeline class", fontsize=10)
    fig.savefig(stem.with_suffix('.png'), dpi=140, bbox_inches='tight')
    plt.close(fig)

    # Relations: does the variability have order?
    fig, axs = plt.subplots(1, 4, figsize=(20, 4.6))
    def scatter(ax, xcol, ycol, xerr=None, yerr=None):
        for frame, filled in [(good, True), (bad, False)]:
            for ap_label, sub in frame.groupby('aperture'):
                sub = sub.dropna(subset=[xcol, ycol])
                if sub.empty:
                    continue
                ax.errorbar(sub[xcol], sub[ycol], xerr=sub[xerr] if xerr else None, yerr=sub[yerr] if yerr else None, fmt=markers.get(sub.program.iloc[0], 'o'),
                            mfc=colors.get(ap_label, 'grey') if filled else 'none', mec=colors.get(ap_label, 'grey'), ecolor=colors.get(ap_label, 'grey'), ms=6, capsize=2, lw=1)
        sub = good.dropna(subset=[xcol, ycol])
        if len(sub) > 3:
            r, pval = pearsonr(sub[xcol], sub[ycol])
            ax.set_title(f'good epochs: Pearson r = {r:.2f} (p = {pval:.1e})', fontsize=9)
        return sub
    sub = scatter(axs[0], 'level_mid', 'flux', 'level_mid_err', 'flux_err')
    if len(sub) > 3:
        q = np.polyfit(np.log10(sub.level_mid), np.log10(sub.flux), 1)
        xx = np.linspace(sub.level_mid.min()*0.8, sub.level_mid.max()*1.1, 50)
        axs[0].plot(xx, 10**np.polyval(q, np.log10(xx)), 'k--', lw=1, label=f'log–log slope {q[0]:.2f}')
        axs[0].legend(fontsize=8)
    axs[0].set_xlabel(f'continuum at {int(np.mean(args.slope))} Å (10⁻¹⁷ erg s⁻¹ cm⁻² Å⁻¹)'); axs[0].set_ylabel(f"{LINES[args.line]['label']} flux (10⁻¹⁷ erg s⁻¹ cm⁻²)")
    scatter(axs[1], 'level_mid', 'norm_slope_per_1000A', 'level_mid_err', None)
    axs[1].set_xlabel(f'continuum at {int(np.mean(args.slope))} Å'); axs[1].set_ylabel(f'normalised slope {int(args.slope[0])}–{int(args.slope[1])} Å (fraction per 1000 Å)')
    sub = scatter(axs[2], 'ztf_g_flam', 'cont_ztfg_pivot', None, 'level_mid_err')
    for frame, filled in [(good, True), (bad, False)]:
        s2 = frame.dropna(subset=['ztf_r_flam', 'cont_ztfr_pivot'])
        if not s2.empty:
            axs[2].plot(s2.ztf_r_flam, s2.cont_ztfr_pivot, 'x' if filled else '+', color='#d62728', ms=7, mew=1.3, ls='none', label=('ZTF r vs 6400 Å (good)' if filled else 'ZTF r vs 6400 Å (flagged)'))
    lim = table[['ztf_g_flam', 'cont_ztfg_pivot', 'ztf_r_flam', 'cont_ztfr_pivot']].stack().dropna()
    if len(lim):
        lo, hi = lim.min()*0.9, lim.max()*1.1
        axs[2].plot([lo, hi], [lo, hi], 'k:', lw=1, label='1:1')
        axs[2].set_xlim(lo, hi); axs[2].set_ylim(lo, hi); axs[2].legend(fontsize=8)
    axs[2].set_xlabel('ZTF within ±20 d, as F_λ at the band pivot (g 4770 Å; r 6400 Å)'); axs[2].set_ylabel('spectroscopic continuum at the same wavelength')
    # Aperture-free colour test: ZTF g-r from same-night pairs against g brightness.
    pairs = []
    if 'g' in ztf and 'r' in ztf:
        g, r = ztf['g'], ztf['r']
        for night in np.unique(np.floor(g[:, 0])):
            gs = g[np.floor(g[:, 0]) == night]
            rs = r[np.abs(r[:, 0]-night-0.5) < 1.0]
            if len(gs) and len(rs):
                pairs.append((night+0.5, np.median(gs[:, 1]), np.median(rs[:, 1]), len(gs), len(rs)))
    pairs = pd.DataFrame(pairs, columns=['mjd', 'g', 'r', 'n_g', 'n_r'])
    if len(pairs) > 5:
        pairs['gr'] = pairs.g-pairs.r
        sc = axs[3].scatter(pairs.g, pairs.gr, c=year(pairs.mjd), cmap='viridis', s=14)
        cb = fig.colorbar(sc, ax=axs[3], pad=0.01)
        cb.set_label('year', fontsize=8)
        q = np.polyfit(pairs.g, pairs.gr, 1)
        rr, pv = pearsonr(pairs.g, pairs.gr)
        xx = np.linspace(pairs.g.min(), pairs.g.max(), 20)
        axs[3].plot(xx, np.polyval(q, xx), 'k--', lw=1, label=f'slope {q[0]:.2f} mag/mag')
        axs[3].set_title(f'ZTF same-night pairs (n={len(pairs)}): Pearson r = {rr:.2f} (p = {pv:.1e})', fontsize=9)
        axs[3].invert_xaxis()
        axs[3].set_xlabel('ZTF g (mag) · brighter to the right'); axs[3].set_ylabel('ZTF g − r (mag) · bluer downward')
        axs[3].legend(fontsize=8)
        pairs.to_csv(stem.with_name(stem.name+'_ztf_pairs.csv'), index=False)
        print(f'ZTF colour test: {len(pairs)} same-night g,r pairs; g-r vs g slope {q[0]:.3f} mag/mag, Pearson r {rr:.2f} (p {pv:.1e}); g range {pairs.g.min():.2f}-{pairs.g.max():.2f}')
    else:
        axs[3].text(.5, .5, 'no same-night ZTF g,r pairs', ha='center', transform=axs[3].transAxes)
    for ax in axs:
        ax.grid(alpha=.25)
    fig.suptitle(f'{args.name}: line–continuum, colour–continuum, spectroscopic vs photometric continuum, and the photometric colour–brightness test (open = flagged epochs)', fontsize=10)
    fig.tight_layout()
    fig.savefig(stem.with_name(stem.name+'_relations.png'), dpi=140)
    plt.close(fig)

    # Inspection sheet: every fit over its data.
    n = len(fits)
    ncol = 5
    nrow = int(np.ceil(n/ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.2*ncol, 2.4*nrow), sharex=True)
    for ax, (e, line) in zip(axes.ravel(), fits):
        if line is None:
            ax.set_title(f"{e['night']} no fit", fontsize=8)
            continue
        ax.plot(line['x'], line['y'], color='k', lw=.7)
        xx = np.linspace(line['x'].min(), line['x'].max(), 400)
        ax.plot(xx, model(xx, *line['params']), color='crimson', lw=1)
        ax.plot(xx, line['params'][0]+line['params'][1]*(xx-line['params'][3]), color='grey', ls='--', lw=.8)
        flag = ' ✗' if (e['fieldquality'] == 'bad' or e['pipeline_class'] not in (None, 'QSO')) else ''
        ax.set_title(f"MJD {e['night']} {e['program']}{flag}\nF={line['flux']:.0f}±{line['flux_err']:.0f} FWHM={line['fwhm_kms']:.0f}", fontsize=7)
        ax.tick_params(labelsize=6)
    for ax in axes.ravel()[n:]:
        ax.axis('off')
    fig.suptitle(f"{args.name}: {LINES[args.line]['label']} fits (observed Å); dashed = local continuum", fontsize=10)
    fig.tight_layout()
    fig.savefig(stem.with_name(stem.name+'_fits.png'), dpi=120)
    plt.close(fig)
    show = ['night_mjd', 'year', 'program', 'fieldquality', 'pipeline_class', 'flux', 'flux_err', 'ew_rest', 'fwhm_kms', 'norm_slope_per_1000A', 'level_mid', 'cont_ztfg_pivot', 'ztf_g_flam', 'cont_ztfr_pivot', 'ztf_r_flam', 'chi2_red']
    print(table[[c for c in show if c in table]].round(2).to_string(index=False))
    print('wrote', stem.with_suffix('.png'), stem.with_suffix('.csv'), stem.with_name(stem.name+'_relations.png'), stem.with_name(stem.name+'_fits.png'))


if __name__ == '__main__':
    main()

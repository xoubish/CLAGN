"""
13g_optical_change_calibration.py  --  the CLAS+ test on our calibration set: optical change since the archival spectrum measured as
ZTF DR24 median magnitude minus the synthetic ZTF magnitude of the DR16 spectrum (13f), versus the SDSS-V outcomes (continuum
dimming, significant Hbeta change).  Compared with the PSF-based change and with the W1 quantities.  Also validates the synthetic
photometry against the SDSS pipeline spectrophotometry and against ZTF for photometrically quiet objects.
Writes data/optical_change_pool.csv (every pool object: d_g, d_r, ZTF amplitudes) for the selection step.
"""
import os, numpy as np, pandas as pd, warnings; warnings.filterwarnings('ignore')
from sklearn.metrics import roc_auc_score
from scipy.stats import spearmanr
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); DATA = os.path.join(HERE, 'data'); pd.set_option('display.width', 250)
syn = pd.read_csv(os.path.join(DATA, 'synthetic_ztf_dr16.csv')).set_index('name'); zo = pd.read_csv(os.path.join(DATA, 'ztf_objects_pool.csv')).set_index('name')
pool = pd.read_csv(os.path.join(DATA, 'parent_pool_scored.csv')); pool = pool[pool.projected.fillna(False)].copy(); pool['name'] = 'P' + pool.poolid.astype(str); pool = pool.set_index('name')
x = pool[['ra', 'dec', 'z', 'psfmag_g', 'psfmag_r', 'plate', 'mjd']].join(syn[['syn_g', 'syn_r', 'snr_g', 'snr_r', 'pipe_g', 'pipe_r', 'mjd_spec']], how='left').join(zo.drop(columns=['ra', 'dec']), how='left')
# sanity: the object table carries 0 or absurd medians for a few entries and tiny ngoodobs for others -> treat as missing
for b in ('g', 'r'):
    bad = (x[f'ztf_{b}_medianmag'] < 10) | (x[f'ztf_{b}_medianmag'] > 23) | (x[f'ztf_{b}_ngoodobs'] < 20)
    for c in ['medianmag', 'minmag', 'maxmag', 'refmag', 'magrms', 'medianabsdev', 'amp']: x.loc[bad, f'ztf_{b}_{c}'] = np.nan
    print(f'ZTF {b}: {int(bad.sum())} entries dropped by the sanity filter')
print(f'pool objects: {len(x)}; with synthetic r {x.syn_r.notna().sum()}, with ZTF r statistics {x.ztf_r_medianmag.notna().sum()}, both {(x.syn_r.notna() & x.ztf_r_medianmag.notna()).sum()}')
print(f'synthetic vs SDSS pipeline spectrophotometry: median syn_r - pipe_r {(x.syn_r - x.pipe_r).median():+.3f}, scatter {1.4826*np.nanmedian(np.abs((x.syn_r - x.pipe_r) - (x.syn_r - x.pipe_r).median())):.3f} mag; syn_r - psfmag_r median {(x.syn_r - x.psfmag_r).median():+.3f}, scatter {1.4826*np.nanmedian(np.abs((x.syn_r - x.psfmag_r) - (x.syn_r - x.psfmag_r).median())):.3f}')
for b in ('g', 'r'):
    x[f'd_{b}'] = x[f'ztf_{b}_medianmag'] - x[f'syn_{b}']              # + = fainter now (ZTF 2018-25 median) than at the archival spectrum
    x[f'd_{b}_psf'] = x[f'ztf_{b}_medianmag'] - x[f'psfmag_{b}']
    x[f'd_{b}_now'] = x[f'ztf_{b}_medianmag'] - x[f'syn_{b}']         # placeholder for a latest-year value when full light curves exist
# archival-spectrum sanity: a spectrum whose synthetic magnitude is > 1.5 mag off the SDSS PSF magnitude of the same era, or with band S/N < 3,
# is a bad or mis-calibrated spectrum (e.g. P4651: syn 22.4 vs psf 18.7, S/N 1.3) and would fake a huge brightening -> no optical change from it
badspec = ((x.syn_r - x.psfmag_r).abs() > 1.5) | ((x.syn_g - x.psfmag_g).abs() > 1.5) | (x.snr_r < 3) | (x.snr_g < 3)
print(f'archival spectra rejected for the optical-change measurement: {int(badspec.sum())} (synthetic vs PSF > 1.5 mag or band S/N < 3)')
x.loc[badspec, ['d_g', 'd_r']] = np.nan; x['badspec'] = badspec
# consistency: g and r must agree to 0.8 mag (real AGN changes are mildly bluer-when-brighter); otherwise treat the optical change as missing
both = x.d_g.notna() & x.d_r.notna(); incons = both & ((x.d_g - x.d_r).abs() > 0.8)
x['d_opt'] = np.where(both & ~incons, 0.5 * (x.d_g + x.d_r), np.where(x.d_g.notna() & ~both, x.d_g, np.where(x.d_r.notna() & ~both, x.d_r, np.nan)))
x['d_opt_flag'] = np.where(incons, 'g/r inconsistent', np.where(both, 'g+r', np.where(x.d_g.notna(), 'g only', np.where(x.d_r.notna(), 'r only', 'none'))))
print(f'optical change d_opt: {x.d_opt.notna().sum()} objects ({int(incons.sum())} g/r-inconsistent set to missing)')
x.to_csv(os.path.join(DATA, 'optical_change_pool.csv'))
c = pd.read_csv(os.path.join(DATA, 'sdssv_calibration_pool.csv')).set_index('name'); ls = pd.read_csv(os.path.join(DATA, 'sdssv_calibration_lines_summary.csv'), index_col=0)
v = c[['dr', 'dimmed_any', 'bright07', 'dW1_full', 'ampW1', 'M', 'in_region_clagn', 'sv_mjd']].join(x[['d_g', 'd_r', 'd_r_psf', 'ztf_r_amp', 'ztf_g_amp', 'ztf_r_magrms', 'ztf_r_medianmag', 'syn_r', 'psfmag_r', 'plate']], how='inner').join(ls[['hb_dim', 'hb_dim3', 'hb_bright']], how='left')
v = v[v.d_r.notna()].copy(); v['dimmed'] = v.dimmed_any.astype(bool); v['bright'] = v.bright07.astype(bool); v['hb_dim'] = v.hb_dim.fillna(False).astype(bool); v['hb_bright'] = v.hb_bright.fillna(False).astype(bool)
print(f'\ncalibration objects with d_r: {len(v)}')
print(f'ZTF-minus-synthetic change vs SDSS-V spectrophotometric change: Spearman {spearmanr(v.d_r, v.dr).correlation:.2f} (PSF-based: {spearmanr(v.d_r_psf, v.dr).correlation:.2f}); median |d_r - dr| {np.median(np.abs(v.d_r - v.dr)):.2f} mag')
q = v[v.dr.abs() < 0.15]; print(f'quiet objects (|dr| < 0.15): median d_r {q.d_r.median():+.3f}, scatter {1.4826*np.median(np.abs(q.d_r - q.d_r.median())):.3f} mag  (PSF-based scatter {1.4826*np.median(np.abs(q.d_r_psf - q.d_r_psf.median())):.3f})')
print('\nAUC:')
for nm, p in [('ZTF - synthetic r (d_r)', v.d_r), ('ZTF - synthetic g (d_g)', v.d_g), ('ZTF - PSF r', v.d_r_psf), ('ZTF r amplitude 2018-25', v.ztf_r_amp), ('ZTF r rms', v.ztf_r_magrms), ('W1 fade 2010-20', v.dW1_full), ('W1 amplitude', v.ampW1), ('M', v.M)]:
    ok = p.notna(); print(f'   {nm:28s} N={ok.sum():4d}  continuum dimmed {roc_auc_score(v.dimmed[ok], p[ok]):.3f} | Hb fell>2x {roc_auc_score(v.hb_dim[ok], p[ok]):.3f} | Hb fell>3x {roc_auc_score(v.hb_dim3[ok].fillna(False).astype(bool), p[ok]):.3f} | brightened {roc_auc_score(v.bright[ok], -p[ok] if nm.startswith(("ZTF -", "W1 fade")) else p[ok]):.3f} | Hb rose>2x {roc_auc_score(v.hb_bright[ok], -p[ok] if nm.startswith(("ZTF -", "W1 fade")) else p[ok]):.3f}')
def tab(col, bins, labels):
    b = pd.cut(v[col], bins, labels=labels); t = v.groupby(b, observed=False)[['dimmed', 'hb_dim', 'hb_dim3', 'bright', 'hb_bright']].mean().mul(100).round(1); t.insert(0, 'N', v.groupby(b, observed=False).size()); return t
print('\noutcome (%) vs ZTF - synthetic r (mag, + = fainter now):'); print(tab('d_r', [-9, -1, -0.5, -0.2, 0.2, 0.5, 1, 9], ['brighter>1', '-1..-0.5', '-0.5..-0.2', 'flat', '0.2-0.5', '0.5-1', 'fainter>1']).to_string())
print('\noutcome (%) vs ZTF - synthetic g:'); print(tab('d_g', [-9, -1, -0.5, -0.2, 0.2, 0.5, 1, 9], ['brighter>1', '-1..-0.5', '-0.5..-0.2', 'flat', '0.2-0.5', '0.5-1', 'fainter>1']).to_string())
print('\njoint triggers:')
for lab, m in [('d_r > 0.5', v.d_r > 0.5), ('d_r > 0.5 & W1 fade > 0.2', (v.d_r > 0.5) & (v.dW1_full > 0.2)), ('d_r > 0.5 & W1 fade <= 0.2', (v.d_r > 0.5) & (v.dW1_full <= 0.2)), ('d_r > 0.3 & W1 amp > 0.4', (v.d_r > 0.3) & (v.ampW1 > 0.4)), ('d_r < -0.5', v.d_r < -0.5), ('d_r < -0.5 & W1 brighten > 0.2', (v.d_r < -0.5) & (v.dW1_full < -0.2)), ('|d_r| < 0.2 & |W1| < 0.1 (control)', (v.d_r.abs() < 0.2) & (v.dW1_full.abs() < 0.1))]:
    s = v[m]; print(f'   {lab:36s} N={len(s):4d}  cont dimmed {100*s.dimmed.mean():5.1f}%  Hb fell>2x {100*s.hb_dim.mean():5.1f}%  >3x {100*s.hb_dim3.fillna(False).astype(bool).mean():5.1f}%  cont brightened {100*s.bright.mean():5.1f}%  Hb rose>2x {100*s.hb_bright.mean():5.1f}%')
d = x[~x.index.isin(pd.read_csv(os.path.join(DATA, 'sdssv_internal_epochs.csv')).name)]
print(f'\ndiscovery pool (no SDSS-V epoch) with d_r: {d.d_r.notna().sum()}; fainter by > 0.5: {(d.d_r > 0.5).sum()}, > 1: {(d.d_r > 1).sum()}; brighter by > 0.5: {(d.d_r < -0.5).sum()}')

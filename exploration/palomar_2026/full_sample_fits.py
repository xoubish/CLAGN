"""Full-sample spectral decomposition (H-beta, H-alpha, Mg II) and A/B/C/D selection experiment.

Scales the Palomar pilot decomposition (pilot_fits.py) to all 99 observed objects
and 365 retained spectra without modifying any delivered reductions or pilot outputs.
All outputs are written to exploration/palomar_2026/outputs/full_sample/.
"""
from pathlib import Path
import hashlib
import json
import time
import warnings
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u
from scipy.ndimage import gaussian_filter1d
from scipy.optimize import least_squares, nnls
from scipy.special import ndtr

import analysis as A
from paper_protocol import OUT as PILOT_OUT
import pilot_fits as PF

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUT = HERE / 'outputs/full_sample'
C = 299792.458

TELLURIC_BANDS = [(5569, 5585), (6293, 6308), (6860, 6960), (7580, 7705), (9300, 9650)]
_TEMPLATE_CACHE = {}


def get_cached_templates(domain, variant):
    """Load and smooth BC03 and Fe II templates once per (domain, variant)."""
    key = (domain, variant)
    if key in _TEMPLATE_CACHE:
        return _TEMPLATE_CACHE[key]
    cfg = PF.VARIANTS[variant]
    folder = PILOT_OUT / 'templates'
    if domain == 'Hb':
        logw = np.arange(np.log(3900), np.log(5800), 30 / C)
        ref_w = 5100.0
    elif domain == 'Ha':
        logw = np.arange(np.log(6000), np.log(7100), 30 / C)
        ref_w = 6450.0
    else:
        logw = np.arange(np.log(2450), np.log(3200), 30 / C)
        ref_w = 2800.0
    grid = np.exp(logw)
    hosts = []
    for name in ['13_ssp_5Gyr_z02.spec.gz', cfg['young']]:
        if ('raw', name) not in _TEMPLATE_CACHE:
            _TEMPLATE_CACHE[('raw', name)] = np.loadtxt(folder / name)
        t = _TEMPLATE_CACHE[('raw', name)]
        f = np.interp(grid, t[:, 0], t[:, 1])
        f = gaussian_filter1d(f, 150 / 30)
        f /= np.maximum(np.interp(ref_w, grid, f), 1e-30)
        hosts.append(f)

    if domain in ('Hb', 'Ha'):
        if ('raw', 'fe_optical.txt') not in _TEMPLATE_CACHE:
            _TEMPLATE_CACHE[('raw', 'fe_optical.txt')] = np.loadtxt(folder / 'fe_optical.txt')
        fe = _TEMPLATE_CACHE[('raw', 'fe_optical.txt')]
        ff = np.interp(grid, 10 ** fe[:, 0], fe[:, 1], left=0, right=0)
        ff = gaussian_filter1d(ff, cfg['fe_sigma'] / 30)
        if domain == 'Hb':
            ff /= np.mean(ff[(grid > 4434) & (grid < 4684)])
        else:
            norm_fe = np.mean(ff[(grid > 6200) & (grid < 6500)])
            if norm_fe > 0:
                ff /= norm_fe
        _TEMPLATE_CACHE[key] = (grid, hosts, ff)
    else:
        uv_fe_blue = np.zeros_like(grid)
        for center, amp, width_a in [
            (2599.0, 1.15, 22.0), (2618.0, 1.10, 20.0), (2632.0, 0.95, 18.0),
            (2666.0, 0.55, 22.0), (2700.0, 0.70, 24.0), (2728.0, 0.90, 18.0),
            (2749.0, 1.25, 16.0), (2769.0, 0.65, 14.0)
        ]:
            uv_fe_blue += amp * np.exp(-0.5 * ((grid - center) / width_a) ** 2)
        uv_fe_red = np.zeros_like(grid)
        for center, amp, width_a in [
            (2836.0, 0.50, 15.0), (2858.0, 0.75, 18.0), (2880.0, 0.85, 20.0),
            (2925.0, 1.10, 24.0), (2950.0, 1.25, 22.0), (2985.0, 1.35, 24.0),
            (3020.0, 0.90, 25.0)
        ]:
            uv_fe_red += amp * np.exp(-0.5 * ((grid - center) / width_a) ** 2)
        uv_fe_blue = gaussian_filter1d(uv_fe_blue, cfg['fe_sigma'] / 30)
        uv_fe_red = gaussian_filter1d(uv_fe_red, cfg['fe_sigma'] / 30)
        uv_fe_blue /= np.maximum(np.mean(uv_fe_blue[(grid > 2600) & (grid < 2750)]), 1e-20)
        uv_fe_red /= np.maximum(np.mean(uv_fe_red[(grid > 2850) & (grid < 3020)]), 1e-20)
        _TEMPLATE_CACHE[key] = (grid, hosts, uv_fe_blue, uv_fe_red)
    return _TEMPLATE_CACHE[key]


def _fast_hb_basis_init(self, wave, lsf, variant='baseline'):
    self.w = np.asarray(wave)
    self.lsf = np.asarray(lsf)
    self.edges = np.r_[
        self.w[0] - (self.w[1] - self.w[0]) / 2,
        (self.w[1:] + self.w[:-1]) / 2,
        self.w[-1] + (self.w[-1] - self.w[-2]) / 2,
    ]
    grid, hosts, ff = get_cached_templates('Hb', variant)
    self.host = [np.interp(self.w, grid, h) for h in hosts]
    self.fe = np.interp(self.w, grid, ff)


PF.Basis.__init__ = _fast_hb_basis_init


def load_all_spectra():
    """Load all 365 retained epochs across the 99 objects, recovering native DESI FITS and SDSS wdisp."""
    OUT.mkdir(parents=True, exist_ok=True)
    objects, visits, _, _ = A.census()
    spectra, _ = A.load_spectra(objects, visits)
    obj_map = objects.set_index('name')
    desi_dir = ROOT / 'data/reselection_2026-09-20/desi_direct'

    audit = []
    for s in spectra:
        status = 'original native spectrum'
        is_native = s['error_kind'] in ('native diagonal', 'native diagonal; correlated pixels')
        if s['source'] == 'DESI':
            obj = obj_map.loc[s['name']]
            matches = []
            for p in sorted(desi_dir.glob(s['name'] + '_*.fits')):
                with fits.open(p, memmap=False) as h:
                    hdr = h[0].header
                    fm = h['FIBERMAP'].data
                    sep = SkyCoord(obj.ra * u.deg, obj.dec * u.deg).separation(
                        SkyCoord(float(fm['TARGET_RA'][0]) * u.deg, float(fm['TARGET_DEC'][0]) * u.deg)
                    ).arcsec
                    if sep < 2.0 and abs(hdr.get('MINMJD', -1) - s['min_mjd']) < 0.005 and abs(hdr.get('MAXMJD', -1) - s['max_mjd']) < 0.005:
                        matches.append(p)
            if len(matches) >= 1:
                p = matches[0]
                arrays = []
                with fits.open(p, memmap=False) as h:
                    for band, lo, hi in [('B', 0, 5750), ('R', 5750, 7550), ('Z', 7550, np.inf)]:
                        w = h[band + '_WAVELENGTH'].data.astype(float)
                        f = h[band + '_FLUX'].data[0].astype(float)
                        iv = h[band + '_IVAR'].data[0].astype(float)
                        good = (iv > 0) & (h[band + '_MASK'].data[0] == 0)
                        err = np.where(good, 1.0 / np.sqrt(np.maximum(iv, 1e-300)), np.nan)
                        f[~good] = np.nan
                        r = np.maximum(h[band + '_RESOLUTION'].data[0].astype(float), 0)
                        off = np.arange(r.shape[0]) - r.shape[0] // 2
                        norm = np.maximum(r.sum(axis=0), 1e-20)
                        mean = (r * off[:, None]).sum(axis=0) / norm
                        sig = np.sqrt(np.maximum((r * off[:, None] ** 2).sum(axis=0) / norm - mean ** 2, 0)) * np.median(np.diff(w))
                        m = (w >= lo) & (w < hi)
                        arrays.append((w[m], f[m], err[m], sig[m], np.full(m.sum(), band)))
                s.update(
                    wave=np.concatenate([a[0] for a in arrays]),
                    flux=np.concatenate([a[1] for a in arrays]),
                    error=np.concatenate([a[2] for a in arrays]),
                    lsf_sigma=np.concatenate([a[3] for a in arrays]),
                    channel=np.concatenate([a[4] for a in arrays]),
                    path=str(p.relative_to(ROOT)),
                    error_kind='native diagonal',
                    alternate=None,
                    uncorrected=None,
                )
                s['quality'] = np.zeros(len(s['wave']), int)
                status = 'DESI native recovered: exact interval and coordinate match'
                s['resolution_status'] = 'DESI resolution-matrix second-moment approximation'
                is_native = True
            else:
                # Estimate empirical local continuum RMS on the 6-A binned SPARCL cache
                # so the epoch can still be inspected, while flagging is_native_error=False.
                w = np.asarray(s['wave'], float)
                f = np.asarray(s['flux'], float)
                finite = np.isfinite(f)
                if finite.sum() > 20:
                    smooth = pd.Series(f).rolling(9, center=True, min_periods=3).median().to_numpy()
                    diff = np.abs(f - smooth)
                    local_mad = pd.Series(diff).rolling(25, center=True, min_periods=5).median().bfill().ffill().to_numpy()
                    global_mad = np.nanmedian(diff) if np.isfinite(diff).any() else 1.0
                    err = 1.4826 * np.where(np.isfinite(local_mad) & (local_mad > 0), local_mad, global_mad)
                    err = np.maximum(err, 0.02 * np.abs(np.where(np.isfinite(f), f, 1.0)))
                    s['error'] = np.where(finite, err, np.nan)
                    s['error_kind'] = 'empirical continuum MAD (6A binned DESI)'
                    status = 'DESI 6A binned cache: empirical continuum MAD error (non-native)'
                else:
                    status = 'DESI native not available and insufficient points'
                is_native = False

        if s['source'] == 'SDSS' and str(s['path']).endswith('.fits'):
            with fits.open(ROOT / s['path'], memmap=False) as h:
                d = h[1].data
                if 'wdisp' in [c.lower() for c in d.names]:
                    s['lsf_sigma'] = np.asarray(d['wdisp'], float) * np.median(np.diff(d['loglam'])) * np.log(10) * s['wave']
                    s['resolution_status'] = 'SDSS native wdisp'
        if 'lsf_sigma' not in s:
            resolution = 1650 if s['source'] == 'NGPS' else (1000 if '6A binned' in s['error_kind'] else 2000)
            s['lsf_sigma'] = s['wave'] / resolution / 2.35482
            s['resolution_status'] = f'assumed constant R={resolution}'
        s['is_native_error'] = bool(is_native)
        audit.append({
            k: s[k] for k in ['name', 'epoch_id', 'source', 'night', 'mjd', 'z', 'path', 'error_kind', 'is_native_error', 'resolution_status']
        } | dict(status=status))

    pd.DataFrame(audit).to_csv(OUT / 'full_native_data_audit.csv', index=False)
    return spectra, objects


def profile_kinematics(x, broad_profile, rest_center):
    """Compute FWHM (km/s), second-moment sigma (km/s), centroid v50 (km/s), and asymmetry from a fitted broad profile."""
    total = np.trapezoid(broad_profile, x)
    if not np.isfinite(total) or total <= 1e-12 or np.max(broad_profile) <= 0:
        return dict(fwhm_kms=np.nan, sigma_kms=np.nan, v50_kms=np.nan, asymmetry=np.nan)
    v_grid = C * (x - rest_center) / rest_center
    # Fine velocity grid for accurate quantiles & FWHM
    vf = np.linspace(v_grid[0], v_grid[-1], 2001)
    pf = np.maximum(np.interp(vf, v_grid, broad_profile), 0.0)
    pmax = np.max(pf)
    above = vf[pf >= 0.5 * pmax]
    fwhm = float(above[-1] - above[0]) if len(above) >= 2 else np.nan
    cdf = np.cumsum(pf)
    cdf /= cdf[-1]
    v10, v50, v90 = np.interp([0.10, 0.50, 0.90], cdf, vf)
    mean_v = float(np.trapezoid(vf * pf, vf) / np.trapezoid(pf, vf))
    var_v = float(np.trapezoid(((vf - mean_v) ** 2) * pf, vf) / np.trapezoid(pf, vf))
    sigma_kms = float(np.sqrt(max(var_v, 0.0)))
    asym = float((v90 + v10 - 2.0 * v50) / max(v90 - v10, 1.0))
    return dict(fwhm_kms=fwhm, sigma_kms=sigma_kms, v50_kms=float(v50), asymmetry=asym)


def augment_hb_metrics(x, mat, coef, p, z):
    """Compute supplementary H-beta metrics (Fe II 4570, R_FeII, core-windowed EW, broad line kinematics) from a Basis solution."""
    parts = mat * coef[None, :]
    continuum = parts[:, :3].sum(axis=1)
    broad = parts[:, 7:9].sum(axis=1)
    fe_prof = parts[:, 3]
    fe_mask = (x >= 4434.0) & (x <= 4684.0)
    feii_4570_flux = float(np.trapezoid(fe_prof[fe_mask], x[fe_mask]) * (1.0 + z)) if fe_mask.sum() > 5 else np.nan
    hb_broad_flux = 100.0 * (coef[7] + coef[8]) * (1.0 + z)
    r_feii = float(feii_4570_flux / hb_broad_flux) if hb_broad_flux > 0 and np.isfinite(feii_4570_flux) else np.nan
    core_mask = (x >= 4750.0) & (x <= 4975.0)
    if core_mask.sum() > 5:
        xc = x[core_mask]; bc = broad[core_mask]; cc = continuum[core_mask]
        hb_flux_core_win = float(np.trapezoid(bc, xc) * (1.0 + z))
        hb_ew_core_win = float(np.trapezoid(np.divide(bc, cc, out=np.zeros_like(bc), where=cc > 0), xc))
    else:
        hb_flux_core_win = np.nan
        hb_ew_core_win = np.nan
    kin = profile_kinematics(x, broad, 4862.68)
    return dict(
        feii4570_flux=feii_4570_flux,
        r_feii=r_feii,
        hb_narrow_flux=float(10.0 * coef[4] * (1.0 + z)),
        hb_flux_core_win=hb_flux_core_win,
        hb_ew_core_win=hb_ew_core_win,
        hb_broad_fwhm_kms=kin['fwhm_kms'],
        hb_broad_sigma_kms=kin['sigma_kms'],
        hb_broad_v50_kms=kin['v50_kms'],
        hb_broad_asymmetry=kin['asymmetry'],
    )


def fit_hb_epoch(s, variant='baseline', save_dir=None, tag=None, cached_npz=None):
    """Fit or reload H-beta decomposition for a single epoch and variant."""
    x, y, e, lsf, good = PF.prepare(s)
    row = {k: s[k] for k in ['name', 'epoch_id', 'source', 'night', 'mjd', 'z', 'error_kind', 'is_native_error', 'resolution_status']}
    row.update(line_complex='Hb', variant=variant, n_pixels=int(good.sum()))
    if good.sum() < 80:
        return row | dict(status='insufficient error pixels')
    coverage = A.integration_weights(x, good, [(4750, 4950)])[1]
    ocoverage = A.integration_weights(x, good, [(4985, 5030)])[1]

    loaded_from_cache = False
    if cached_npz is not None and cached_npz.exists():
        d = np.load(cached_npz, allow_pickle=True)
        opt_x = d['parameters']
        coef = d['coefficients']
        basis = PF.Basis(x, lsf, variant)
        mat = basis.matrix(opt_x)
        res = (mat[good] @ coef - y[good]) / e[good]
        opt_success = True
        nfev = 0
        loaded_from_cache = True
    else:
        opt, mat, coef, res = PF.fit_arrays(x, y, e, lsf, good, variant)
        opt_x = opt.x
        opt_success = bool(opt.success)
        nfev = int(opt.nfev)

    near = (opt_x - PF.LOW < 0.015 * (PF.HIGH - PF.LOW)) | (PF.HIGH - opt_x < 0.015 * (PF.HIGH - PF.LOW))
    dof = max(1, int(good.sum()) - len(opt_x) - int((coef > 1e-8).sum()))
    row.update(
        status='fit',
        optimizer_success=opt_success,
        nfev=nfev,
        reduced_chi2=float(np.sum(res ** 2) / dof),
        hb_window_coverage=float(coverage),
        oiii_window_coverage=float(ocoverage),
        parameter_bounds=';'.join(np.array(PF.PARAMS)[near]),
        **PF.measurements(x, mat, coef, opt_x, s['z']),
        **augment_hb_metrics(x, mat, coef, opt_x, s['z']),
    )
    row['hb_supported'] = bool(coverage >= 0.80)
    row['hb_telluric_free'] = bool(coverage >= 0.95)
    row['oiii_supported'] = bool(ocoverage >= 0.80)
    if save_dir is not None and not loaded_from_cache:
        save_dir.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(
            save_dir / f"{s['epoch_id']}__{tag or variant}.npz",
            wave_rest=x, flux=y, error=e, good=good, lsf_sigma_rest=lsf,
            parameters=opt_x, coefficients=coef, components=mat * coef[None, :], component_names=PF.NAMES
        )
    return row


# =============================================================================
# H-alpha + [N II] + [S II] + [O I] Decomposition (6200 - 6850 A rest-frame)
# =============================================================================
HA_NAMES = [
    'agn', 'host_5gyr', 'host_young', 'feii',
    'ha_narrow', 'nii_doublet', 'sii6718', 'sii6733', 'oi_doublet',
    'ha_broad1', 'ha_broad2'
]
HA_LOW = np.array([-3.5, 0.05, -0.5, 0.4, 1.0, -2.0, -3.0])
HA_HIGH = np.array([2.5, 0.5, 0.5, 4.5, 7.5, 2.0, 3.0])


class HaBasis:
    def __init__(self, wave, lsf, variant='baseline'):
        self.w = np.asarray(wave)
        self.lsf = np.asarray(lsf)
        self.edges = np.r_[
            self.w[0] - (self.w[1] - self.w[0]) / 2,
            (self.w[1:] + self.w[:-1]) / 2,
            self.w[-1] + (self.w[-1] - self.w[-2]) / 2,
        ]
        grid, hosts, ff = get_cached_templates('Ha', variant)
        self.host = [np.interp(self.w, grid, h) for h in hosts]
        self.fe = np.interp(self.w, grid, ff)

    def gaussian(self, center, sigma, v):
        mu = center * (1 + v / C)
        inst = float(np.interp(mu, self.w, self.lsf))
        sig = np.hypot(center * sigma / C, inst)
        return np.diff(ndtr((self.edges - mu) / sig)) / np.diff(self.edges)

    def matrix(self, p):
        alpha, ns, nv, s1, s2, v1, v2 = p
        ns *= 1000; nv *= 1000; s1 *= 1000; s2 *= 1000; v1 *= 1000; v2 *= 1000
        g = self.gaussian
        nii = 10 * (g(6585.27, ns, nv) + g(6549.86, ns, nv) / 2.96)
        oi = 10 * (g(6302.05, ns, nv) + g(6365.54, ns, nv) / 3.0)
        return np.column_stack([
            (self.w / 6450.0) ** alpha,
            *self.host,
            self.fe,
            10 * g(6564.61, ns, nv),
            nii,
            10 * g(6718.29, ns, nv),
            10 * g(6732.67, ns, nv),
            oi,
            100 * g(6564.61, s1, v1),
            100 * g(6564.61, s2, v2),
        ])


def prepare_ha(s, flux_key='flux'):
    w = s['wave'] / (1 + s['z'])
    m = (w >= 6200) & (w <= 6850)
    x = w[m]; y = np.asarray(s[flux_key])[m]; e = s['error'][m]
    lsf = s['lsf_sigma'][m] / (1 + s['z'])
    lsf = np.where(np.isfinite(lsf) & (lsf > 0), lsf, x / 1650 / 2.35482)
    good = np.isfinite(y) & np.isfinite(e) & (e > 0)
    observed = x * (1 + s['z'])
    for a, b in TELLURIC_BANDS:
        good &= ~((observed >= a) & (observed <= b))
    return x, y, e, lsf, good


def fit_ha_epoch(s, variant='baseline'):
    x, y, e, lsf, good = prepare_ha(s)
    row = {k: s[k] for k in ['name', 'epoch_id', 'source', 'night', 'mjd', 'z', 'error_kind', 'is_native_error', 'resolution_status']}
    row.update(line_complex='Ha', variant=variant, n_pixels=int(good.sum()))
    if good.sum() < 60:
        return row | dict(status='insufficient error pixels')
    coverage = A.integration_weights(x, good, [(6450, 6680)])[1]
    scoverage = A.integration_weights(x, good, [(6705, 6750)])[1]
    if coverage < 0.60:
        return row | dict(status='insufficient H-alpha coverage', ha_window_coverage=float(coverage))

    basis = HaBasis(x, lsf, variant)
    yg = y[good]; eg = e[good]

    def solve(p, return_parts=False):
        mat = basis.matrix(p)
        a = mat[good] / eg[:, None]
        b = yg / eg
        norms = np.maximum(np.linalg.norm(a, axis=0), 1e-12)
        coef = nnls(a / norms, b, maxiter=500)[0] / norms
        residual = (mat[good] @ coef - yg) / eg
        return (mat, coef, residual) if return_parts else residual

    seeds = [[-1.5, 0.16, 0.0, 1.2, 2.8, 0.0, 0.0], [-0.5, 0.25, 0.0, 0.8, 4.0, -0.2, 0.2]]
    attempts = [
        least_squares(solve, p, bounds=(HA_LOW, HA_HIGH), max_nfev=100, ftol=2e-6, xtol=2e-6, gtol=2e-6, x_scale='jac')
        for p in seeds
    ]
    opt = min(attempts, key=lambda r: np.sum(r.fun ** 2))
    mat, coef, res = solve(opt.x, True)
    parts = mat * coef[None, :]
    continuum = parts[:, :3].sum(axis=1)
    broad = parts[:, 9:11].sum(axis=1)
    ew = float(np.trapezoid(np.divide(broad, continuum, out=np.zeros_like(broad), where=continuum > 0), x))
    ha_broad = float(100.0 * (coef[9] + coef[10]) * (1.0 + s['z']))
    ha_narrow = float(10.0 * coef[4] * (1.0 + s['z']))
    nii_flux = float(10.0 * coef[5] * (1.0 + s['z']))
    sii_flux = float(10.0 * (coef[6] + coef[7]) * (1.0 + s['z']))
    agn6450 = float(coef[0])
    total6450 = float(np.interp(6450.0, x, continuum))
    kin = profile_kinematics(x, broad, 6564.61)
    dof = max(1, int(good.sum()) - len(opt.x) - int((coef > 1e-8).sum()))
    row.update(
        status='fit',
        optimizer_success=bool(opt.success),
        nfev=int(opt.nfev),
        reduced_chi2=float(np.sum(res ** 2) / dof),
        ha_window_coverage=float(coverage),
        sii_window_coverage=float(scoverage),
        ha_supported=bool(coverage >= 0.90),
        sii_supported=bool(scoverage >= 0.90),
        ha_flux=ha_broad,
        ha_narrow_flux=ha_narrow,
        nii6585_flux=nii_flux,
        sii_flux=sii_flux,
        ha_sii=ha_broad / sii_flux if sii_flux > 0 else np.nan,
        ha_nii=ha_broad / nii_flux if nii_flux > 0 else np.nan,
        ha_ew_total_cont=ew,
        agn6450=agn6450,
        total6450=total6450,
        host_fraction6450=1.0 - agn6450 / total6450 if total6450 > 0 else np.nan,
        ha_broad_fwhm_kms=kin['fwhm_kms'],
        ha_broad_sigma_kms=kin['sigma_kms'],
        ha_broad_v50_kms=kin['v50_kms'],
        ha_broad_asymmetry=kin['asymmetry'],
    )
    return row


# =============================================================================
# Mg II + UV Fe II Decomposition (2600 - 3050 A rest-frame)
# =============================================================================
MGII_LOW = np.array([-3.5, 0.08, -0.6, 0.5, 1.2, -2.0, -3.0])
MGII_HIGH = np.array([2.5, 0.6, 0.6, 5.0, 8.0, 2.0, 3.0])


class MgIIBasis:
    def __init__(self, wave, lsf, variant='baseline'):
        self.w = np.asarray(wave)
        self.lsf = np.asarray(lsf)
        self.edges = np.r_[
            self.w[0] - (self.w[1] - self.w[0]) / 2,
            (self.w[1:] + self.w[:-1]) / 2,
            self.w[-1] + (self.w[-1] - self.w[-2]) / 2,
        ]
        grid, hosts, uv_fe_blue, uv_fe_red = get_cached_templates('MgII', variant)
        self.host = [np.interp(self.w, grid, h) for h in hosts]
        self.fe_blue = np.interp(self.w, grid, uv_fe_blue)
        self.fe_red = np.interp(self.w, grid, uv_fe_red)

    def gaussian(self, center, sigma, v):
        mu = center * (1 + v / C)
        inst = float(np.interp(mu, self.w, self.lsf))
        sig = np.hypot(center * sigma / C, inst)
        return np.diff(ndtr((self.edges - mu) / sig)) / np.diff(self.edges)

    def matrix(self, p):
        alpha, ns, nv, s1, s2, v1, v2 = p
        ns *= 1000; nv *= 1000; s1 *= 1000; s2 *= 1000; v1 *= 1000; v2 *= 1000
        g = self.gaussian
        narrow = 10 * (0.55 * g(2796.35, ns, nv) + 0.45 * g(2803.53, ns, nv))
        broad1 = 100 * (0.5 * g(2796.35, s1, v1) + 0.5 * g(2803.53, s1, v1))
        broad2 = 100 * (0.5 * g(2796.35, s2, v2) + 0.5 * g(2803.53, s2, v2))
        return np.column_stack([
            (self.w / 2800.0) ** alpha,
            *self.host,
            self.fe_blue,
            self.fe_red,
            narrow,
            broad1,
            broad2,
        ])


def prepare_mgii(s, flux_key='flux'):
    w = s['wave'] / (1 + s['z'])
    m = (w >= 2600) & (w <= 3050)
    x = w[m]; y = np.asarray(s[flux_key])[m]; e = s['error'][m]
    lsf = s['lsf_sigma'][m] / (1 + s['z'])
    lsf = np.where(np.isfinite(lsf) & (lsf > 0), lsf, x / 1650 / 2.35482)
    good = np.isfinite(y) & np.isfinite(e) & (e > 0)
    observed = x * (1 + s['z'])
    for a, b in TELLURIC_BANDS:
        good &= ~((observed >= a) & (observed <= b))
    return x, y, e, lsf, good


def fit_mgii_epoch(s, variant='baseline'):
    x, y, e, lsf, good = prepare_mgii(s)
    row = {k: s[k] for k in ['name', 'epoch_id', 'source', 'night', 'mjd', 'z', 'error_kind', 'is_native_error', 'resolution_status']}
    row.update(line_complex='MgII', variant=variant, n_pixels=int(good.sum()))
    if good.sum() < 50:
        return row | dict(status='insufficient error pixels')
    coverage = A.integration_weights(x, good, [(2720, 2880)])[1]
    if coverage < 0.65:
        return row | dict(status='insufficient Mg II coverage', mgii_window_coverage=float(coverage))

    basis = MgIIBasis(x, lsf, variant)
    yg = y[good]; eg = e[good]

    def solve(p, return_parts=False):
        mat = basis.matrix(p)
        a = mat[good] / eg[:, None]
        b = yg / eg
        norms = np.maximum(np.linalg.norm(a, axis=0), 1e-12)
        coef = nnls(a / norms, b, maxiter=500)[0] / norms
        residual = (mat[good] @ coef - yg) / eg
        return (mat, coef, residual) if return_parts else residual

    seeds = [[-1.5, 0.20, 0.0, 1.4, 3.2, 0.0, 0.0], [-0.8, 0.30, 0.0, 0.9, 4.5, -0.2, 0.2]]
    attempts = [
        least_squares(solve, p, bounds=(MGII_LOW, MGII_HIGH), max_nfev=100, ftol=2e-6, xtol=2e-6, gtol=2e-6, x_scale='jac')
        for p in seeds
    ]
    opt = min(attempts, key=lambda r: np.sum(r.fun ** 2))
    mat, coef, res = solve(opt.x, True)
    parts = mat * coef[None, :]
    continuum = parts[:, :3].sum(axis=1)
    broad = parts[:, 6:8].sum(axis=1)
    ew = float(np.trapezoid(np.divide(broad, continuum, out=np.zeros_like(broad), where=continuum > 0), x))
    mgii_broad = float(100.0 * (coef[6] + coef[7]) * (1.0 + s['z']))
    mgii_narrow = float(10.0 * coef[5] * (1.0 + s['z']))
    fe_uv = float(np.trapezoid(parts[:, 3] + parts[:, 4], x) * (1.0 + s['z']))
    agn2800 = float(coef[0])
    total2800 = float(np.interp(2800.0, x, continuum))
    kin = profile_kinematics(x, broad, 2799.94)
    dof = max(1, int(good.sum()) - len(opt.x) - int((coef > 1e-8).sum()))
    row.update(
        status='fit',
        optimizer_success=bool(opt.success),
        nfev=int(opt.nfev),
        reduced_chi2=float(np.sum(res ** 2) / dof),
        mgii_window_coverage=float(coverage),
        mgii_supported=bool(coverage >= 0.92),
        mgii_flux=mgii_broad,
        mgii_narrow_flux=mgii_narrow,
        feii_uv_flux=fe_uv,
        mgii_ew_total_cont=ew,
        agn2800=agn2800,
        total2800=total2800,
        mgii_broad_fwhm_kms=kin['fwhm_kms'],
        mgii_broad_sigma_kms=kin['sigma_kms'],
        mgii_broad_v50_kms=kin['v50_kms'],
        mgii_broad_asymmetry=kin['asymmetry'],
    )
    return row


def build_comparisons(df_fits, epochs_df, line_prefix, keys, support_col):
    """Build multi-epoch comparison ratios and model sensitivity envelopes for a given line complex."""
    rows = []
    good = df_fits[df_fits.status.eq('fit')].copy()
    for (name, variant), g in good.groupby(['name', 'variant']):
        old = g[g.source.ne('NGPS')].sort_values('mjd')
        old_native = old[old.is_native_error]
        new = g[g.source.eq('NGPS')].sort_values('mjd')
        ep_obj = epochs_df[epochs_df.name.eq(name)]
        for n in new.itertuples(index=False):
            before_all = ep_obj[ep_obj.source.ne('NGPS') & ep_obj.max_mjd.lt(n.mjd)].sort_values('max_mjd')
            refs = []
            if len(before_all):
                latest_id = before_all.iloc[-1].epoch_id
                f_lat = old[old.epoch_id.eq(latest_id)]
                if len(f_lat):
                    refs.append(('latest_archive', f_lat.iloc[0]))
                if len(before_all) > 1:
                    earliest_id = before_all.iloc[0].epoch_id
                    f_ear = old[old.epoch_id.eq(earliest_id)]
                    if len(f_ear):
                        refs.append(('earliest_archive', f_ear.iloc[0]))
            # Also guarantee a strict native-error archival reference comparison
            before_nat = old_native[old_native.mjd < n.mjd].sort_values('mjd')
            if len(before_nat):
                refs.append(('latest_native_archive', before_nat.iloc[-1]))
            repeat = new[new.mjd < n.mjd]
            if len(repeat):
                refs.append(('previous_ngps', repeat.iloc[-1]))

            for role, o in refs:
                r = dict(
                    name=name,
                    line_complex=line_prefix,
                    variant=variant,
                    comparison=role,
                    epoch_id=n.epoch_id,
                    reference_epoch=o.epoch_id,
                    night=n.night,
                    old_source=o.source,
                    both_native_error=bool(n.is_native_error and o.is_native_error),
                    baseline_rest_years=float((n.mjd - o.mjd) / 365.25 / (1.0 + n.z)),
                    both_line_supported=bool(getattr(n, support_col) and o[support_col]),
                    max_reduced_chi2=float(max(n.reduced_chi2, o.reduced_chi2)),
                    both_optimizers_success=bool(n.optimizer_success and o.optimizer_success),
                )
                if 'oiii_supported' in g.columns:
                    r['both_oiii_supported'] = bool(n.oiii_supported and o['oiii_supported'])
                for k in keys:
                    nv = getattr(n, k, np.nan)
                    ov = o.get(k, np.nan)
                    r[f'old_{k}'] = float(ov) if np.isfinite(ov) else np.nan
                    r[f'new_{k}'] = float(nv) if np.isfinite(nv) else np.nan
                    r[f'{k}_ratio'] = float(nv / ov) if (np.isfinite(nv) and np.isfinite(ov) and ov > 0) else np.nan
                    r[f'{k}_delta'] = float(nv - ov) if (np.isfinite(nv) and np.isfinite(ov)) else np.nan
                for k in [f'{line_prefix.lower()}_broad_fwhm_kms', f'{line_prefix.lower()}_broad_v50_kms']:
                    if k in g.columns:
                        nv = getattr(n, k, np.nan)
                        ov = o.get(k, np.nan)
                        r[f'old_{k}'] = float(ov) if np.isfinite(ov) else np.nan
                        r[f'new_{k}'] = float(nv) if np.isfinite(nv) else np.nan
                        r[f'{k}_delta'] = float(nv - ov) if (np.isfinite(nv) and np.isfinite(ov)) else np.nan
                rows.append(r)

    pairs = pd.DataFrame(rows)
    if pairs.empty:
        return pairs, pd.DataFrame()

    env = []
    for (name, epoch, reference, role), g in pairs.groupby(['name', 'epoch_id', 'reference_epoch', 'comparison']):
        base_rows = g[g.variant.eq('baseline')]
        b = base_rows.iloc[0] if len(base_rows) else g.iloc[0]
        r = dict(
            name=name,
            line_complex=line_prefix,
            epoch_id=epoch,
            reference_epoch=reference,
            comparison=role,
            night=b['night'],
            old_source=b['old_source'],
            both_native_error=bool(g.both_native_error.all()),
            baseline_rest_years=float(b['baseline_rest_years']),
            n_variants=int(len(g)),
            both_line_supported=bool(g.both_line_supported.all()),
            max_reduced_chi2=float(g.max_reduced_chi2.max()),
            both_optimizers_success=bool(g.both_optimizers_success.all()),
        )
        if 'both_oiii_supported' in g.columns:
            r['both_oiii_supported'] = bool(g.both_oiii_supported.all())
        for k in keys:
            r[f'old_{k}_baseline'] = float(b[f'old_{k}'])
            r[f'new_{k}_baseline'] = float(b[f'new_{k}'])
            r[f'{k}_ratio_baseline'] = float(b[f'{k}_ratio'])
            r[f'{k}_ratio_min'] = float(g[f'{k}_ratio'].min())
            r[f'{k}_ratio_max'] = float(g[f'{k}_ratio'].max())
            r[f'{k}_delta_baseline'] = float(b[f'{k}_delta'])
        env.append(r)
    return pairs, pd.DataFrame(env)


def run_selection_experiment(hb_env, mgii_env, analysis_sample):
    """Evaluate the A/B/C/D manifold selection experiment on the full fitted sample."""
    # Take one primary comparison per object: latest_native_archive (latest NGPS visit if observed twice)
    sub_hb = hb_env[
        hb_env.comparison.eq('latest_native_archive') & hb_env.n_variants.eq(3)
    ].sort_values('epoch_id').groupby('name', as_index=False).last()

    merged = analysis_sample.merge(sub_hb, on='name', how='left', suffixes=('', '_hb'))

    # Also merge Mg II primary comparison for high-z objects
    if not mgii_env.empty:
        sub_mg = mgii_env[
            mgii_env.comparison.eq('latest_native_archive') & mgii_env.both_line_supported & mgii_env.n_variants.eq(3)
        ].sort_values('epoch_id').groupby('name', as_index=False).last()
        mg_cols = ['name', 'epoch_id', 'reference_epoch', 'old_mgii_ew_total_cont_baseline', 'new_mgii_ew_total_cont_baseline',
                   'mgii_ew_total_cont_ratio_baseline', 'mgii_ew_total_cont_ratio_min', 'mgii_ew_total_cont_ratio_max']
        merged = merged.merge(sub_mg[mg_cols], on='name', how='left', suffixes=('', '_mgii'))

    # Envelope-robust change indicators for H-beta
    # Use core-windowed EW ratio (4750-4975 A) alongside total EW ratio so ultra-wide broad2 boundary hits cannot create false positives
    for factor, label in [(1.5, '1p5x'), (2.0, '2p0x')]:
        ew_up = (merged['hb_ew_total_cont_ratio_min'] >= factor) & (merged['hb_ew_core_win_ratio_min'] >= 1.35)
        ew_dn = (merged['hb_ew_total_cont_ratio_max'] <= (1.0 / factor)) & (merged['hb_ew_core_win_ratio_max'] <= (1.0 / 1.35))
        oiii_up = merged['hb_oiii_ratio_min'] >= (factor if factor == 2.0 else 1.40)
        oiii_dn = merged['hb_oiii_ratio_max'] <= (1.0 / (factor if factor == 2.0 else 1.40))
        merged[f'hb_ew_robust_{label}'] = ew_up | ew_dn
        merged[f'hb_oiii_robust_{label}'] = oiii_up | oiii_dn
        merged[f'hb_both_robust_{label}'] = (ew_up & oiii_up) | (ew_dn & oiii_dn)

    merged['log10_hb_ew_ratio'] = np.log10(merged['hb_ew_core_win_ratio_baseline'].where(merged['hb_ew_core_win_ratio_baseline'] > 0))
    merged['abs_log10_hb_ew_ratio'] = merged['log10_hb_ew_ratio'].abs()
    merged['log10_hb_oiii_ratio'] = np.log10(merged['hb_oiii_ratio_baseline'].where(merged['hb_oiii_ratio_baseline'] > 0))
    merged['abs_log10_hb_oiii_ratio'] = merged['log10_hb_oiii_ratio'].abs()

    # Directional alignment with pre-run manifold neighborhood direction (>0 strengthening vs <0 weakening)
    valid_dir = merged['manifold_direction'].notna() & (merged['manifold_direction'] != 0) & merged['log10_hb_ew_ratio'].notna()
    merged['direction_aligned_ew'] = np.where(
        valid_dir,
        np.sign(merged['log10_hb_ew_ratio']) == np.sign(merged['manifold_direction']),
        np.nan
    )
    merged['direction_aligned_oiii'] = np.where(
        valid_dir & merged['log10_hb_oiii_ratio'].notna(),
        np.sign(merged['log10_hb_oiii_ratio']) == np.sign(merged['manifold_direction']),
        np.nan
    )
    merged.to_csv(OUT / 'full_sample_object_summary.csv', index=False)

    # Group-level summary table across A, B, C, D, and outside pre-run review
    hb_valid = merged[merged['hb_ew_total_cont_ratio_baseline'].notna()].copy()
    group_rows = []
    for grp, g in hb_valid.groupby('experiment_group'):
        g_clean = g[g.known_state_status.ne('catalog-confirmed CLAGN')]
        g_strict = g_clean[g_clean.both_line_supported.eq(True)]
        for subset_label, sub in [
            ('all_hb_fitted', g),
            ('excluding_catalog_confirmed_clagn', g_clean),
            ('clean_coverage_gt_0p80_excl_confirmed', g_strict),
        ]:
            if len(sub) == 0:
                continue
            dir_sub = sub[sub['manifold_direction'].notna() & (sub['manifold_direction'] != 0) & sub['log10_hb_ew_ratio'].notna()]
            group_rows.append(dict(
                experiment_group=grp,
                subset=subset_label,
                n_objects=int(len(sub)),
                median_abs_log10_hb_ew_ratio=float(sub['abs_log10_hb_ew_ratio'].median()),
                iqr_abs_log10_hb_ew_ratio=float(sub['abs_log10_hb_ew_ratio'].quantile(0.75) - sub['abs_log10_hb_ew_ratio'].quantile(0.25)),
                mean_abs_log10_hb_ew_ratio=float(sub['abs_log10_hb_ew_ratio'].mean()),
                median_abs_log10_hb_oiii_ratio=float(sub['abs_log10_hb_oiii_ratio'].median()),
                n_hb_ew_robust_1p5x=int(sub['hb_ew_robust_1p5x'].sum()),
                frac_hb_ew_robust_1p5x=float(sub['hb_ew_robust_1p5x'].mean()),
                n_hb_both_robust_1p5x=int(sub['hb_both_robust_1p5x'].sum()),
                frac_hb_both_robust_1p5x=float(sub['hb_both_robust_1p5x'].mean()),
                n_hb_both_robust_2p0x=int(sub['hb_both_robust_2p0x'].sum()),
                frac_hb_both_robust_2p0x=float(sub['hb_both_robust_2p0x'].mean()),
                n_with_manifold_direction=int(len(dir_sub)),
                n_direction_aligned_ew=int((np.sign(dir_sub['log10_hb_ew_ratio']) == np.sign(dir_sub['manifold_direction'])).sum()) if len(dir_sub) else 0,
                frac_direction_aligned_ew=float((np.sign(dir_sub['log10_hb_ew_ratio']) == np.sign(dir_sub['manifold_direction'])).mean()) if len(dir_sub) else np.nan,
            ))
    df_groups = pd.DataFrame(group_rows)
    df_groups.to_csv(OUT / 'abcd_group_experiment_summary.csv', index=False)

    # Bootstrap & permutation hypothesis tests (10,000 resamples, fixed seed)
    rng = np.random.default_rng(20261104)
    clean_abcd = hb_valid[
        hb_valid.experiment_group.str.slice(0, 2).isin(['A:', 'B:', 'C:', 'D:'])
        & hb_valid.known_state_status.ne('catalog-confirmed CLAGN')
    ].copy()

    def bootstrap_diff(arr_x, arr_y, stat_fn=np.median, n_boot=10000):
        arr_x = np.asarray(arr_x, float)[np.isfinite(arr_x)]
        arr_y = np.asarray(arr_y, float)[np.isfinite(arr_y)]
        if len(arr_x) < 2 or len(arr_y) < 2:
            return dict(obs_diff=np.nan, ci_2p5=np.nan, ci_97p5=np.nan, perm_p_two_sided=np.nan)
        obs = float(stat_fn(arr_x) - stat_fn(arr_y))
        bx = stat_fn(rng.choice(arr_x, size=(n_boot, len(arr_x)), replace=True), axis=1)
        by = stat_fn(rng.choice(arr_y, size=(n_boot, len(arr_y)), replace=True), axis=1)
        diffs = bx - by
        pooled = np.concatenate([arr_x, arr_y])
        nx = len(arr_x)
        perm_diffs = []
        for _ in range(n_boot):
            p = rng.permutation(pooled)
            perm_diffs.append(stat_fn(p[:nx]) - stat_fn(p[nx:]))
        perm_diffs = np.asarray(perm_diffs)
        p_val = float(np.mean(np.abs(perm_diffs) >= abs(obs)))
        return dict(
            n_x=int(len(arr_x)), n_y=int(len(arr_y)),
            stat_x=float(stat_fn(arr_x)), stat_y=float(stat_fn(arr_y)),
            obs_diff=obs,
            ci_2p5=float(np.quantile(diffs, 0.025)),
            ci_97p5=float(np.quantile(diffs, 0.975)),
            perm_p_two_sided=p_val,
        )

    grpA = clean_abcd[clean_abcd.experiment_group.str.startswith('A:')]
    grpB = clean_abcd[clean_abcd.experiment_group.str.startswith('B:')]
    grpD = clean_abcd[clean_abcd.experiment_group.str.startswith('D:')]
    grpAB = clean_abcd[clean_abcd.experiment_group.str.slice(0, 2).isin(['A:', 'B:'])]

    in_clagn = clean_abcd[clean_abcd.in_region_clagn.eq(True)]
    out_clagn = clean_abcd[clean_abcd.in_region_clagn.eq(False)]

    pos_dir = clean_abcd[clean_abcd.manifold_direction > 0]
    neg_dir = clean_abcd[clean_abcd.manifold_direction < 0]

    forecast_sub = clean_abcd[clean_abcd.temporal_forecast_candidate.eq(True)]
    f_pos = forecast_sub[forecast_sub.manifold_direction > 0]
    f_neg = forecast_sub[forecast_sub.manifold_direction < 0]

    # Spearman rank correlation with W1 amplitude and manifold features
    from scipy.stats import spearmanr
    corrs = {}
    for feat in ['w1_p95_p05_mag', 'w1_endpoint_delta_mag', 'manifold_cl_neighbor_fraction', 'manifold_on_neighbor_fraction', 'manifold_direction']:
        m = clean_abcd[feat].notna() & clean_abcd['abs_log10_hb_ew_ratio'].notna()
        if m.sum() >= 5:
            rho, pval = spearmanr(clean_abcd.loc[m, feat], clean_abcd.loc[m, 'abs_log10_hb_ew_ratio'])
            corrs[f'{feat}__vs__abs_log10_hb_ew_ratio'] = dict(n=int(m.sum()), spearman_rho=float(rho), p_value=float(pval))
        m2 = clean_abcd[feat].notna() & clean_abcd['log10_hb_ew_ratio'].notna()
        if m2.sum() >= 5:
            rho2, pval2 = spearmanr(clean_abcd.loc[m2, feat], clean_abcd.loc[m2, 'log10_hb_ew_ratio'])
            corrs[f'{feat}__vs__signed_log10_hb_ew_ratio'] = dict(n=int(m2.sum()), spearman_rho=float(rho2), p_value=float(pval2))

    tests = dict(
        A_vs_B_median_abs_log10_hb_ew=bootstrap_diff(grpA['abs_log10_hb_ew_ratio'], grpB['abs_log10_hb_ew_ratio'], np.median),
        A_vs_D_median_abs_log10_hb_ew=bootstrap_diff(grpA['abs_log10_hb_ew_ratio'], grpD['abs_log10_hb_ew_ratio'], np.median),
        AB_vs_D_median_abs_log10_hb_ew=bootstrap_diff(grpAB['abs_log10_hb_ew_ratio'], grpD['abs_log10_hb_ew_ratio'], np.median),
        A_vs_B_frac_robust_1p5x=bootstrap_diff(grpA['hb_both_robust_1p5x'].astype(float), grpB['hb_both_robust_1p5x'].astype(float), np.mean),
        A_vs_D_frac_robust_1p5x=bootstrap_diff(grpA['hb_both_robust_1p5x'].astype(float), grpD['hb_both_robust_1p5x'].astype(float), np.mean),
        in_clagn_region_vs_outside_median_abs_log10_hb_ew=bootstrap_diff(in_clagn['abs_log10_hb_ew_ratio'], out_clagn['abs_log10_hb_ew_ratio'], np.median),
        manifold_direction_pos_vs_neg_signed_log10_hb_ew_all=bootstrap_diff(pos_dir['log10_hb_ew_ratio'], neg_dir['log10_hb_ew_ratio'], np.median),
        manifold_direction_pos_vs_neg_signed_log10_hb_ew_forecast_only=bootstrap_diff(f_pos['log10_hb_ew_ratio'], f_neg['log10_hb_ew_ratio'], np.median),
        spearman_correlations=corrs,
    )
    (OUT / 'selection_experiment_tests.json').write_text(json.dumps(tests, indent=2) + '\n')
    return merged, df_groups, tests


def main():
    t0 = time.monotonic()
    spectra, objects = load_all_spectra()
    epochs_df = pd.read_csv(HERE / 'outputs/epochs.csv')
    analysis_sample = pd.read_csv(PILOT_OUT / 'analysis_sample.csv')

    print(f"Loaded {len(spectra)} total spectral epochs across {len(objects)} objects.")

    # 1. H-beta decomposition across all epochs (all 3 variants + alternate standards)
    hb_rows = []
    hb_alt_rows = []
    pilot_fits_dir = PILOT_OUT / 'fits'
    full_fits_dir = OUT / 'fits_hb'
    for i, s in enumerate(spectra, 1):
        for variant in PF.VARIANTS:
            c_full = full_fits_dir / f"{s['epoch_id']}__{variant}.npz"
            c_pilot = pilot_fits_dir / f"{s['epoch_id']}__{variant}.npz"
            cached = c_full if c_full.exists() else (c_pilot if c_pilot.exists() else None)
            r = fit_hb_epoch(s, variant=variant, save_dir=full_fits_dir, cached_npz=cached)
            hb_rows.append(r)
        # Also fit NGPS alternate standard star for calibration envelope
        alt = s.get('alternate')
        if s['source'] == 'NGPS' and alt is not None and np.isfinite(alt).any():
            a = dict(s)
            ratio = np.divide(alt, s['flux'], out=np.full(len(alt), np.nan), where=s['flux'] != 0)
            a['flux'] = alt.copy()
            a['error'] = s['error'] * np.abs(ratio)
            for variant in PF.VARIANTS:
                c_full_alt = full_fits_dir / f"{s['epoch_id']}__{variant}_alternate.npz"
                c_pilot_alt = pilot_fits_dir / f"{s['epoch_id']}__{variant}_alternate.npz"
                cached_alt = c_full_alt if c_full_alt.exists() else (c_pilot_alt if c_pilot_alt.exists() else None)
                ra = fit_hb_epoch(a, variant=variant, save_dir=full_fits_dir, tag=variant + '_alternate',
                                  cached_npz=cached_alt)
                hb_alt_rows.append(ra)
        if i % 80 == 0 or i == len(spectra):
            n_fit = sum(1 for r in hb_rows if r['status'] == 'fit' and r['variant'] == 'baseline')
            print(f"[H-beta] Processed {i}/{len(spectra)} epochs ({n_fit} baseline fits so far, {time.monotonic()-t0:.1f}s)", flush=True)

    df_hb = pd.DataFrame(hb_rows).fillna({'parameter_bounds': ''})
    df_hb.to_csv(OUT / 'hb_fit_results.csv', index=False)
    df_hb_alt = pd.DataFrame(hb_alt_rows).fillna({'parameter_bounds': ''})
    df_hb_alt.to_csv(OUT / 'hb_alternate_standard_fits.csv', index=False)

    hb_keys = ['hb_ew_total_cont', 'hb_ew_core_win', 'hb_oiii', 'hb_flux', 'hb_flux_core_win', 'oiii5008_flux', 'agn5100', 'total5100', 'host_fraction5100', 'r_feii']
    hb_pairs, hb_env = build_comparisons(df_hb, epochs_df, 'Hb', hb_keys, 'hb_supported')

    # Add calibration stress envelope to hb_env using hb_alt_rows
    if not df_hb_alt.empty and not hb_pairs.empty:
        base_idx = df_hb[df_hb.status.eq('fit')].set_index(['epoch_id', 'variant'])
        alt_idx = df_hb_alt[df_hb_alt.status.eq('fit')].set_index(['epoch_id', 'variant'])
        cal_min, cal_max = {}, {}
        for (name, epoch, ref, comp), g in hb_pairs.groupby(['name', 'epoch_id', 'reference_epoch', 'comparison']):
            for k in ['hb_ew_total_cont', 'hb_ew_core_win', 'hb_oiii', 'hb_flux', 'agn5100']:
                vals = []
                for row in g.itertuples(index=False):
                    base_ratio = getattr(row, f'{k}_ratio')
                    if not np.isfinite(base_ratio):
                        continue
                    vals.append(base_ratio)
                    nk = (row.epoch_id, row.variant)
                    ok = (row.reference_epoch, row.variant)
                    n_mult = 1.0
                    if nk in alt_idx.index and nk in base_idx.index and base_idx.loc[nk, k] > 0:
                        n_mult = float(alt_idx.loc[nk, k] / base_idx.loc[nk, k])
                    o_mult = 1.0
                    if ok in alt_idx.index and ok in base_idx.index and base_idx.loc[ok, k] > 0:
                        o_mult = float(alt_idx.loc[ok, k] / base_idx.loc[ok, k])
                    for n_use in [1.0, n_mult]:
                        for o_use in [1.0, o_mult]:
                            if o_use > 0:
                                vals.append(base_ratio * n_use / o_use)
                if vals:
                    cal_min[(epoch, ref, comp, k)] = float(np.nanmin(vals))
                    cal_max[(epoch, ref, comp, k)] = float(np.nanmax(vals))
        for k in ['hb_ew_total_cont', 'hb_ew_core_win', 'hb_oiii', 'hb_flux', 'agn5100']:
            hb_env[f'{k}_cal_model_min'] = [
                cal_min.get((r.epoch_id, r.reference_epoch, r.comparison, k), getattr(r, f'{k}_ratio_min'))
                for r in hb_env.itertuples(index=False)
            ]
            hb_env[f'{k}_cal_model_max'] = [
                cal_max.get((r.epoch_id, r.reference_epoch, r.comparison, k), getattr(r, f'{k}_ratio_max'))
                for r in hb_env.itertuples(index=False)
            ]

    hb_pairs.to_csv(OUT / 'hb_fit_comparisons.csv', index=False)
    hb_env.to_csv(OUT / 'hb_model_sensitivity_envelope.csv', index=False)

    # 2. H-alpha decomposition across all covered epochs (reuse cached CSV if present)
    if (OUT / 'ha_fit_results.csv').exists():
        df_ha = pd.read_csv(OUT / 'ha_fit_results.csv')
    else:
        ha_rows = []
        for i, s in enumerate(spectra, 1):
            for variant in PF.VARIANTS:
                ha_rows.append(fit_ha_epoch(s, variant=variant))
        df_ha = pd.DataFrame(ha_rows)
        df_ha.to_csv(OUT / 'ha_fit_results.csv', index=False)
    ha_keys = ['ha_ew_total_cont', 'ha_sii', 'ha_nii', 'ha_flux', 'sii_flux', 'agn6450', 'total6450']
    ha_pairs, ha_env = build_comparisons(df_ha, epochs_df, 'Ha', ha_keys, 'ha_supported')
    ha_pairs.to_csv(OUT / 'ha_fit_comparisons.csv', index=False)
    ha_env.to_csv(OUT / 'ha_model_sensitivity_envelope.csv', index=False)

    # 3. Mg II decomposition across all covered epochs (reuse cached CSV if present)
    if (OUT / 'mgii_fit_results.csv').exists():
        df_mgii = pd.read_csv(OUT / 'mgii_fit_results.csv')
    else:
        mgii_rows = []
        for i, s in enumerate(spectra, 1):
            for variant in PF.VARIANTS:
                mgii_rows.append(fit_mgii_epoch(s, variant=variant))
        df_mgii = pd.DataFrame(mgii_rows)
        df_mgii.to_csv(OUT / 'mgii_fit_results.csv', index=False)
    mgii_keys = ['mgii_ew_total_cont', 'mgii_flux', 'feii_uv_flux', 'agn2800', 'total2800']
    mgii_pairs, mgii_env = build_comparisons(df_mgii, epochs_df, 'MgII', mgii_keys, 'mgii_supported')
    mgii_pairs.to_csv(OUT / 'mgii_fit_comparisons.csv', index=False)
    mgii_env.to_csv(OUT / 'mgii_model_sensitivity_envelope.csv', index=False)

    # 4. Run the full-sample A/B/C/D selection experiment
    merged, df_groups, tests = run_selection_experiment(hb_env, mgii_env, analysis_sample)
    print("\n=== Full-Sample A/B/C/D Group Summary (H-beta fitted) ===")
    print(df_groups[
        ['experiment_group', 'subset', 'n_objects', 'median_abs_log10_hb_ew_ratio', 'n_hb_ew_robust_1p5x',
         'n_hb_both_robust_1p5x', 'n_hb_both_robust_2p0x', 'n_with_manifold_direction', 'frac_direction_aligned_ew']
    ].to_string(index=False))
    print(f"\nCompleted full-sample fits in {time.monotonic()-t0:.1f}s.")


if __name__ == '__main__':
    main()

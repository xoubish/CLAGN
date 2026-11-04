"""Step 2: Deep-dive joint multi-epoch H-beta + H-alpha decomposition, physical mechanism tests
(accretion vs. dust extinction, virial breathing), exposure audit, and literature resolution
for the top Palomar CLAGN discovery leads.

Writes all outputs to exploration/palomar_2026/outputs/deep_dive/ and figures to
exploration/palomar_2026/outputs/full_sample/figures/ without touching pilot outputs.
"""
from pathlib import Path
import csv
import json
import time
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from astropy.cosmology import FlatLambdaCDM
from astropy.coordinates import SkyCoord, search_around_sky
import astropy.units as u
from astropy.io import fits
from scipy.optimize import least_squares, nnls

import analysis as A
from paper_protocol import OUT as PILOT_OUT, SCIENCE as PILOT_SCIENCE
import pilot_fits as PF
import full_sample_fits as FSF

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
FULL_OUT = HERE / 'outputs/full_sample'
DEEP_OUT = HERE / 'outputs/deep_dive'
FIG_DIR = FULL_OUT / 'figures'
C = 299792.458
COSMO = FlatLambdaCDM(H0=70.0, Om0=0.3)

# Cardelli, Clayton & Mathis (1989) Rv=3.1 A_lambda / E(B-V) values
K_HB = 3.609
K_HA = 2.535
K_5100 = 3.417


def select_deep_dive_leads(obj_summary, hb_env, ha_env, mgii_env):
    """Select the 11 pilot science leads plus all newly identified full-sample leads."""
    leads = set(PILOT_SCIENCE) | {
        'P14388', 'P12457', 'P16455', 'P21743', 'P7346', 'P7457', 'P8897', 'P11542', 'P11010', 'P21841', 'S729737389072738304'
    }
    hb_rob = obj_summary[obj_summary['hb_both_robust_1p5x'].eq(True) | obj_summary['hb_ew_robust_1p5x'].eq(True)]['name'].tolist()
    leads.update(hb_rob)
    if not ha_env.empty:
        ha_sub = ha_env[ha_env.comparison.eq('latest_native_archive') & ha_env.both_line_supported & ha_env.n_variants.eq(3)]
        ha_rob = ha_sub[
            ((ha_sub['ha_ew_total_cont_ratio_min'] >= 1.45) & (ha_sub['ha_sii_ratio_min'] >= 1.35))
            | ((ha_sub['ha_ew_total_cont_ratio_max'] <= 1.0 / 1.45) & (ha_sub['ha_sii_ratio_max'] <= 1.0 / 1.35))
        ]['name'].tolist()
        leads.update(ha_rob)
    if not mgii_env.empty:
        mg_sub = mgii_env[mgii_env.comparison.eq('latest_native_archive') & mg_env.both_line_supported & mg_env.n_variants.eq(3)] if 'mg_env' in locals() else mgii_env[mgii_env.comparison.eq('latest_native_archive') & mgii_env.both_line_supported & mgii_env.n_variants.eq(3)]
        mg_rob = mg_sub[
            (mg_sub['mgii_ew_total_cont_ratio_min'] >= 1.6)
            | (mg_sub['mgii_ew_total_cont_ratio_max'] <= 1.0 / 1.6)
        ]['name'].tolist()
        leads.update(mg_rob)
    return sorted(leads)


def joint_multi_epoch_hb_fit(obj_spectra, z):
    """Jointly fit all H-beta epochs of a source with a tied intrinsic host-galaxy SED mix
    (theta_young = f_young / (f_5Gyr + f_young)) and tied narrow-line width sigma_n.

    Returns per-epoch joint-fit rows and the shared host stellar-population fraction.
    """
    prepped = []
    for s in obj_spectra:
        x, y, e, lsf, good = PF.prepare(s)
        if good.sum() >= 80:
            cov = A.integration_weights(x, good, [(4750, 4950)])[1]
            ocov = A.integration_weights(x, good, [(4985, 5030)])[1]
            if cov >= 0.60 and ocov >= 0.40:
                prepped.append((s, x, y, e, lsf, good, PF.Basis(x, lsf, 'baseline')))
    if len(prepped) == 0:
        return [], np.nan

    # First get single-epoch baseline parameters as warm starts, and regularize any ultra-wide broad2 boundary hits
    LOW_CLEAN = np.array([-3.5, 0.06, -0.5, 0.4, 0.8, -1.2, -1.5])
    HIGH_CLEAN = np.array([2.5, 0.55, 0.5, 3.8, 4.2, 1.2, 1.5])
    single_fits = []
    for s, x, y, e, lsf, good, basis in prepped:
        cached = FULL_OUT / 'fits_hb' / f"{s['epoch_id']}__baseline.npz"
        if cached.exists():
            d = np.load(cached, allow_pickle=True)
            p0 = d['parameters'].copy()
            c0 = d['coefficients'].copy()
        else:
            opt, _, c0, _ = PF.fit_arrays(x, y, e, lsf, good, 'baseline', starts=1)
            p0 = opt.x.copy()
        # If broad2 hit the ultra-wide continuum-bump regime (sigma2 > 4200 km/s or |v2| > 1500 km/s), refit within clean physical bounds
        if p0[4] > 4.2 or abs(p0[6]) > 1.5:
            yg = y[good]; eg = e[good]
            def solve_clean(p):
                mat = basis.matrix(p); a = mat[good] / eg[:, None]; b = yg / eg
                norms = np.maximum(np.linalg.norm(a, axis=0), 1e-12)
                coef = nnls(a / norms, b, maxiter=500)[0] / norms
                return (mat[good] @ coef - yg) / eg
            seed = np.clip(p0, LOW_CLEAN + 0.05, HIGH_CLEAN - 0.05)
            seed[4] = min(seed[4], 2.8); seed[6] = 0.0
            opt_c = least_squares(solve_clean, seed, bounds=(LOW_CLEAN, HIGH_CLEAN), max_nfev=80, ftol=2e-6, xtol=2e-6, x_scale='jac')
            p0 = opt_c.x
            mat_c = basis.matrix(p0); a_c = mat_c[good] / eg[:, None]; norms_c = np.maximum(np.linalg.norm(a_c, axis=0), 1e-12)
            c0 = nnls(a_c / norms_c, yg / eg, maxiter=500)[0] / norms_c
        single_fits.append((p0, c0))

    # Shared parameters to optimize across all epochs of this object:
    #   theta_young in [0, 1] (intrinsic young-to-total host fraction at 5100 A)
    #   narrow_sigma in [0.06, 0.55] (1000 km/s)
    # For each epoch, we keep its broad kinematics (s1, s2, v1, v2) and narrow centroid nv
    # from the nonlinear fit and optimize power-law slope alpha_e + non-negative amplitudes
    # with the composite locked-shape host template H(theta_young) = (1-theta)*H_5Gyr + theta*H_1.4Gyr!
    def eval_joint(shared_p, return_epochs=False):
        theta_young, shared_ns = shared_p
        all_res = []
        ep_out = []
        for idx, (s, x, y, e, lsf, good, basis) in enumerate(prepped):
            p_ep = single_fits[idx][0].copy()
            p_ep[1] = shared_ns
            raw_mat = basis.matrix(p_ep)
            # Collapse columns 1 (host_5gyr) and 2 (host_young) into a single locked-shape composite host column
            comp_host = (1.0 - theta_young) * raw_mat[:, 1] + theta_young * raw_mat[:, 2]
            mat = np.column_stack([raw_mat[:, 0], comp_host, raw_mat[:, 3:]])
            yg = y[good]; eg = e[good]
            a = mat[good] / eg[:, None]
            b = yg / eg
            norms = np.maximum(np.linalg.norm(a, axis=0), 1e-12)
            coef_red = nnls(a / norms, b, maxiter=500)[0] / norms
            res = (mat[good] @ coef_red - yg) / eg
            all_res.append(res)
            if return_epochs:
                # Expand coef_red back to 13-component vector matching PF.NAMES
                full_coef = np.r_[
                    coef_red[0],
                    coef_red[1] * (1.0 - theta_young),
                    coef_red[1] * theta_young,
                    coef_red[2:],
                ]
                ep_out.append((s, x, y, e, lsf, good, raw_mat, full_coef, p_ep, res))
        return ep_out if return_epochs else np.concatenate(all_res)

    # Initial guess for theta_young and shared_ns from flux-weighted average of single-epoch fits
    h5_sum = sum(c0[1] for _, c0 in single_fits)
    hy_sum = sum(c0[2] for _, c0 in single_fits)
    theta0 = float(np.clip(hy_sum / max(h5_sum + hy_sum, 1e-12), 0.05, 0.95))
    ns0 = float(np.clip(np.median([p0[1] for p0, _ in single_fits]), 0.08, 0.50))

    opt_joint = least_squares(
        eval_joint, [theta0, ns0], bounds=([0.0, 0.06], [1.0, 0.55]),
        max_nfev=60, ftol=1e-5, xtol=1e-5
    )
    theta_best, ns_best = opt_joint.x
    ep_solutions = eval_joint(opt_joint.x, return_epochs=True)

    rows = []
    for s, x, y, e, lsf, good, raw_mat, full_coef, p_ep, res in ep_solutions:
        dof = max(1, int(good.sum()) - 2 - int((full_coef > 1e-8).sum()))
        m_dict = PF.measurements(x, raw_mat, full_coef, p_ep, z)
        aug_dict = FSF.augment_hb_metrics(x, raw_mat, full_coef, p_ep, z)
        host5100 = float(full_coef[1] + full_coef[2])
        rows.append(dict(
            name=s['name'],
            epoch_id=s['epoch_id'],
            source=s['source'],
            night=s['night'],
            mjd=float(s['mjd']),
            z=float(z),
            is_native_error=bool(s['is_native_error']),
            shared_theta_young=float(theta_best),
            shared_narrow_sigma_kms=float(ns_best * 1000.0),
            joint_reduced_chi2=float(np.sum(res ** 2) / dof),
            joint_host5100=host5100,
            joint_host_over_oiii=float(host5100 / m_dict['oiii5008_flux']) if m_dict['oiii5008_flux'] > 0 else np.nan,
            **{f'joint_{k}': v for k, v in m_dict.items()},
            **{f'joint_{k}': v for k, v in aug_dict.items()},
        ))
    return rows, float(theta_best)


def compute_virial_and_extinction(df_joint, df_ha_base):
    """Merge H-alpha baseline measurements with joint H-beta fits to compute:
    - Broad & narrow Balmer decrements (Ha / Hb)
    - Extinction E(B-V) test (comparing observed Balmer decrement change vs pure-dust prediction)
    - Luminosities L5100_AGN, L_Hb_broad, single-epoch virial BH mass log(M_BH/M_sun) (VP06), and Eddington ratio.
    """
    ha_sub = df_ha_base[df_ha_base.variant.eq('baseline') & df_ha_base.status.eq('fit')][
        ['epoch_id', 'ha_supported', 'ha_flux', 'ha_narrow_flux', 'nii6585_flux', 'sii_flux',
         'ha_sii', 'ha_ew_total_cont', 'agn6450', 'total6450', 'ha_broad_fwhm_kms', 'ha_broad_v50_kms']
    ]
    m = df_joint.merge(ha_sub, on='epoch_id', how='left')

    # Luminosity distance in cm
    dl_cm = COSMO.luminosity_distance(m['z'].to_numpy()).to(u.cm).value
    area = 4.0 * np.pi * (dl_cm ** 2)
    # Observed flux is in 1e-17 erg s^-1 cm^-2 A^-1 (continuum at rest 5100 -> f_lam_obs * (1+z) at rest)
    # Wait: in PF.measurements, hb_flux = 100*(coef[7]+coef[8])*(1+z) is already integrated over observed wavelength (1e-17 erg/s/cm^2).
    # For continuum at rest 5100 A: f_lambda_rest(5100) = f_lambda_obs(5100*(1+z)) * (1+z), so lambda_rest * f_lambda_rest = 5100 * (1+z) * f_5100_obs.
    l5100_agn = 5100.0 * (1.0 + m['z']) * np.maximum(m['joint_agn5100'], 0.0) * 1e-17 * area
    l5100_tot = 5100.0 * (1.0 + m['z']) * np.maximum(m['joint_total5100'], 0.0) * 1e-17 * area
    l_hb = np.maximum(m['joint_hb_flux'], 0.0) * 1e-17 * area

    m['log_L5100_agn_ergs'] = np.where(l5100_agn > 0, np.log10(l5100_agn), np.nan)
    m['log_L5100_tot_ergs'] = np.where(l5100_tot > 0, np.log10(l5100_tot), np.nan)
    m['log_LHb_broad_ergs'] = np.where(l_hb > 0, np.log10(l_hb), np.nan)

    # Vestergaard & Peterson (2006) single-epoch virial BH mass estimators:
    # 1) Continuum-based: log(M_BH/Msun) = 6.91 + 0.50 * log(L5100 / 1e44) + 2 * log(FWHM_kms / 1000)
    # 2) Line-luminosity-based (immune to host-subtraction error):
    #    log(M_BH/Msun) = 6.67 + 0.63 * log(L_Hb / 1e42) + 2 * log(FWHM_kms / 1000)
    fwhm = m['joint_hb_broad_fwhm_kms']
    valid_fwhm = np.isfinite(fwhm) & (fwhm > 300.0)
    m['log_MBH_VP06_L5100'] = np.where(
        valid_fwhm & (l5100_agn > 0),
        6.91 + 0.50 * (m['log_L5100_agn_ergs'] - 44.0) + 2.0 * np.log10(fwhm / 1000.0),
        np.nan
    )
    m['log_MBH_VP06_LHb'] = np.where(
        valid_fwhm & (l_hb > 0),
        6.67 + 0.63 * (m['log_LHb_broad_ergs'] - 42.0) + 2.0 * np.log10(fwhm / 1000.0),
        np.nan
    )
    # Bolometric luminosity L_bol = 9.26 * L5100_AGN (Richards et al. 2006) or from L_Hb via Greene & Ho (2005):
    # L5100_from_LHb = 1e44 * (L_Hb / (1.425e42))**(1/1.133)
    l5100_from_lhb = 1e44 * np.power(np.maximum(l_hb / 1.425e42, 1e-20), 1.0 / 1.133)
    # Use median BH mass per object so Eddington ratio reflects state change at fixed BH mass
    obj_mbh = m.groupby('name')['log_MBH_VP06_LHb'].transform('median')
    m['log_MBH_adopted'] = obj_mbh
    l_edd = 1.26e38 * np.power(10.0, obj_mbh)
    m['log_Edd_ratio_from_LHb'] = np.where(
        (l_hb > 0) & np.isfinite(l_edd),
        np.log10(9.26 * l5100_from_lhb / l_edd),
        np.nan
    )

    # Broad and narrow Balmer decrement (Ha / Hb)
    m['balmer_decrement_broad'] = np.where(
        (m['ha_flux'] > 0) & (m['joint_hb_flux'] > 0),
        m['ha_flux'] / m['joint_hb_flux'],
        np.nan
    )
    m['balmer_decrement_narrow'] = np.where(
        (m['ha_narrow_flux'] > 0) & (m['joint_hb_narrow_flux'] > 0),
        m['ha_narrow_flux'] / m['joint_hb_narrow_flux'],
        np.nan
    )
    # Color excess E(B-V)_broad assuming intrinsic Case B / AGN broad-line ratio Ha/Hb = 3.06 (Dong et al. 2008)
    m['ebv_broad_mag'] = np.where(
        m['balmer_decrement_broad'] > 0,
        (2.5 / (K_HB - K_HA)) * np.log10(m['balmer_decrement_broad'] / 3.06),
        np.nan
    )
    return m


def evaluate_lead_transitions(df_phys, hb_env):
    """Compare each lead's NGPS epoch against its latest native archival epoch under the joint fit,
    testing pure dust extinction vs intrinsic accretion change and virial breathing."""
    rows = []
    env_lookup = hb_env[hb_env.comparison.eq('latest_native_archive')].set_index(['name', 'epoch_id'])
    for name, g in df_phys.groupby('name'):
        old = g[g.source.ne('NGPS') & g.is_native_error].sort_values('mjd')
        new = g[g.source.eq('NGPS')].sort_values('mjd')
        if len(old) == 0 or len(new) == 0:
            continue
        o = old.iloc[-1]
        for n in new.itertuples(index=False):
            hb_ew_ratio_joint = n.joint_hb_ew_total_cont / o.joint_hb_ew_total_cont if o.joint_hb_ew_total_cont > 0 else np.nan
            hb_oiii_ratio_joint = n.joint_hb_oiii / o.joint_hb_oiii if o.joint_hb_oiii > 0 else np.nan
            ha_ew_ratio = n.ha_ew_total_cont / o.ha_ew_total_cont if (np.isfinite(n.ha_ew_total_cont) and np.isfinite(o.ha_ew_total_cont) and o.ha_ew_total_cont > 0) else np.nan
            ha_sii_ratio = n.ha_sii / o.ha_sii if (np.isfinite(n.ha_sii) and np.isfinite(o.ha_sii) and o.ha_sii > 0) else np.nan

            # Extinction test: if the observed change in Hb/[OIII] were caused purely by foreground dust extinction
            # Delta E(B-V)_req = -2.5 * log10(hb_oiii_ratio_joint) / K_HB
            delta_ebv_req_if_dust = -2.5 * np.log10(hb_oiii_ratio_joint) / K_HB if (np.isfinite(hb_oiii_ratio_joint) and hb_oiii_ratio_joint > 0) else np.nan
            pred_bd_ratio_if_dust = float(10.0 ** (0.4 * (K_HB - K_HA) * delta_ebv_req_if_dust)) if np.isfinite(delta_ebv_req_if_dust) else np.nan
            obs_bd_ratio = float(n.balmer_decrement_broad / o.balmer_decrement_broad) if (np.isfinite(n.balmer_decrement_broad) and np.isfinite(o.balmer_decrement_broad) and o.balmer_decrement_broad > 0) else np.nan
            obs_delta_ebv_from_bd = float(n.ebv_broad_mag - o.ebv_broad_mag) if (np.isfinite(n.ebv_broad_mag) and np.isfinite(o.ebv_broad_mag)) else np.nan

            # Virial breathing slope: d log(FWHM) / d log(L_Hb) (expected = -0.25 for R_BLR ~ L^0.5)
            d_log_lhb = float(n.log_LHb_broad_ergs - o.log_LHb_broad_ergs) if (np.isfinite(n.log_LHb_broad_ergs) and np.isfinite(o.log_LHb_broad_ergs)) else np.nan
            d_log_fwhm = float(np.log10(n.joint_hb_broad_fwhm_kms / o.joint_hb_broad_fwhm_kms)) if (np.isfinite(n.joint_hb_broad_fwhm_kms) and np.isfinite(o.joint_hb_broad_fwhm_kms) and o.joint_hb_broad_fwhm_kms > 0) else np.nan
            breathing_slope = float(d_log_fwhm / d_log_lhb) if (np.isfinite(d_log_lhb) and abs(d_log_lhb) > 0.08 and np.isfinite(d_log_fwhm)) else np.nan

            # Look up independent-fit envelope from Step 1
            ew_min, ew_max, oiii_min, oiii_max = np.nan, np.nan, np.nan, np.nan
            if (name, n.epoch_id) in env_lookup.index:
                erow = env_lookup.loc[(name, n.epoch_id)]
                if isinstance(erow, pd.DataFrame):
                    erow = erow.iloc[0]
                ew_min = float(erow['hb_ew_total_cont_ratio_min'])
                ew_max = float(erow['hb_ew_total_cont_ratio_max'])
                oiii_min = float(erow['hb_oiii_ratio_min'])
                oiii_max = float(erow['hb_oiii_ratio_max'])

            # Physical mechanism verdict
            mechanism = 'intrinsic ionizing continuum / accretion state change'
            if np.isfinite(obs_delta_ebv_from_bd) and np.isfinite(delta_ebv_req_if_dust) and abs(delta_ebv_req_if_dust) > 0.15:
                # If observed Balmer decrement change accounts for <40% of required dust extinction, pure reddening is ruled out
                frac_dust = obs_delta_ebv_from_bd / delta_ebv_req_if_dust
                if frac_dust < 0.40:
                    mechanism = f'Accretion-driven (Balmer decrement rules out pure dust: dE(B-V)_obs={obs_delta_ebv_from_bd:+.2f} vs {delta_ebv_req_if_dust:+.2f} mag required)'
                else:
                    mechanism = f'Coupled accretion + dust reddening (dE(B-V)_obs={obs_delta_ebv_from_bd:+.2f} vs {delta_ebv_req_if_dust:+.2f} mag)'

            rows.append(dict(
                name=name,
                z=float(n.z),
                old_epoch=o.epoch_id,
                new_epoch=n.epoch_id,
                baseline_rest_years=float((n.mjd - o.mjd) / 365.25 / (1.0 + n.z)),
                shared_theta_young=float(n.shared_theta_young),
                old_hb_ew_joint=float(o.joint_hb_ew_total_cont),
                new_hb_ew_joint=float(n.joint_hb_ew_total_cont),
                hb_ew_ratio_joint=float(hb_ew_ratio_joint),
                hb_ew_ratio_indep_min=ew_min,
                hb_ew_ratio_indep_max=ew_max,
                old_hb_oiii_joint=float(o.joint_hb_oiii),
                new_hb_oiii_joint=float(n.joint_hb_oiii),
                hb_oiii_ratio_joint=float(hb_oiii_ratio_joint),
                hb_oiii_ratio_indep_min=oiii_min,
                hb_oiii_ratio_indep_max=oiii_max,
                host_over_oiii_ratio_joint=float(n.joint_host_over_oiii / o.joint_host_over_oiii) if o.joint_host_over_oiii > 0 else np.nan,
                old_ha_ew=float(o.ha_ew_total_cont) if np.isfinite(o.ha_ew_total_cont) else np.nan,
                new_ha_ew=float(n.ha_ew_total_cont) if np.isfinite(n.ha_ew_total_cont) else np.nan,
                ha_ew_ratio=ha_ew_ratio,
                ha_sii_ratio=ha_sii_ratio,
                old_bd_broad=float(o.balmer_decrement_broad) if np.isfinite(o.balmer_decrement_broad) else np.nan,
                new_bd_broad=float(n.balmer_decrement_broad) if np.isfinite(n.balmer_decrement_broad) else np.nan,
                obs_bd_ratio=obs_bd_ratio,
                pred_bd_ratio_if_pure_dust=pred_bd_ratio_if_dust,
                delta_ebv_req_if_pure_dust=float(delta_ebv_req_if_dust) if np.isfinite(delta_ebv_req_if_dust) else np.nan,
                obs_delta_ebv_from_bd=obs_delta_ebv_from_bd,
                old_fwhm_hb_kms=float(o.joint_hb_broad_fwhm_kms),
                new_fwhm_hb_kms=float(n.joint_hb_broad_fwhm_kms),
                delta_v50_hb_kms=float(n.joint_hb_broad_v50_kms - o.joint_hb_broad_v50_kms),
                log_MBH_adopted=float(n.log_MBH_adopted) if np.isfinite(n.log_MBH_adopted) else np.nan,
                old_log_Edd_ratio=float(o.log_Edd_ratio_from_LHb) if np.isfinite(o.log_Edd_ratio_from_LHb) else np.nan,
                new_log_Edd_ratio=float(n.log_Edd_ratio_from_LHb) if np.isfinite(n.log_Edd_ratio_from_LHb) else np.nan,
                delta_log_Edd_ratio=float(n.log_Edd_ratio_from_LHb - o.log_Edd_ratio_from_LHb) if (np.isfinite(n.log_Edd_ratio_from_LHb) and np.isfinite(o.log_Edd_ratio_from_LHb)) else np.nan,
                virial_breathing_slope=breathing_slope,
                physical_mechanism_assessment=mechanism,
            ))
    return pd.DataFrame(rows)


def audit_literature_cross_matches(objects_df, lead_names):
    """Cross-match all 99 objects (highlighting discovery leads) against Camus-Panda 2026,
    earlier literature_clagn.csv (Graham et al. 2019/2020, MacLeod et al. 2016/2019), and Zeltyn et al. 2024."""
    coords = SkyCoord(objects_df.ra.to_numpy() * u.deg, objects_df.dec.to_numpy() * u.deg)
    camus = pd.read_csv(ROOT / 'data/external/clagn_catalog_camus_panda2026.csv', low_memory=False)
    camus = camus[pd.to_numeric(camus['ra_deg'], errors='coerce').notna() & pd.to_numeric(camus['dec_deg'], errors='coerce').notna()].reset_index(drop=True)
    lit = pd.read_csv(ROOT / 'data/literature_clagn.csv')
    lit = lit[pd.to_numeric(lit['ra'], errors='coerce').notna() & pd.to_numeric(lit['dec'], errors='coerce').notna()].reset_index(drop=True)
    zel = pd.read_csv(ROOT / 'data/zeltyn_coords.csv')
    zel = zel[pd.to_numeric(zel['ra'], errors='coerce').notna() & pd.to_numeric(zel['dec'], errors='coerce').notna()].reset_index(drop=True)

    c_sky = SkyCoord(camus['ra_deg'].astype(float).to_numpy() * u.deg, camus['dec_deg'].astype(float).to_numpy() * u.deg)
    l_sky = SkyCoord(lit['ra'].astype(float).to_numpy() * u.deg, lit['dec'].astype(float).to_numpy() * u.deg)
    z_sky = SkyCoord(zel['ra'].astype(float).to_numpy() * u.deg, zel['dec'].astype(float).to_numpy() * u.deg)

    rows = []
    for idx, obj in objects_df.iterrows():
        sc = SkyCoord(obj.ra * u.deg, obj.dec * u.deg)
        c_sep = sc.separation(c_sky).arcsec
        l_sep = sc.separation(l_sky).arcsec
        z_sep = sc.separation(z_sky).arcsec

        c_matches = camus[c_sep <= 3.0]
        l_matches = lit[l_sep <= 3.0]
        z_matches = zel[z_sep <= 3.0]

        if len(c_matches) == 0 and len(l_matches) == 0 and len(z_matches) == 0 and obj['name'] not in lead_names:
            continue

        c_status = '; '.join(f"{r['main_name']} ({r['confirmation_status']}, ref={r['reference']}, last_spec={r['last_spectrum_date']})" for _, r in c_matches.iterrows()) if len(c_matches) else 'none'
        l_status = '; '.join(f"ref={r['ref']}" for _, r in l_matches.iterrows()) if len(l_matches) else 'none'
        z_status = '; '.join(f"{r.get('name','')} ({r.get('class_zeltyn','')})" for _, r in z_matches.iterrows()) if len(z_matches) else 'none'

        if obj['name'] == 'P21422':
            resolution = (
                'NEW RECURRENT TURN-OFF DISCOVERY: Matched Graham et al. (2019/2020) CRTS photometric candidate list '
                '(not in Camus & Panda 2026 confirmed CLAGN catalog). Archival SDSS/DESI spectra (2005->2018->2022) '
                'showed an earlier brightening to a strong Type 1 state (EW_Hb=66.6 A in 2022 DESI); our 2026 Oct 27 '
                'Palomar spectrum captures a subsequent post-2022 turn-off collapse (EW_Hb=11.5 A, Hb/[OIII] down 6.7x).'
            )
        elif obj['name'] == 'P10381':
            resolution = (
                'FIRST SPECTROSCOPIC TURN-ON CONFIRMATION: Matched Graham et al. (2019/2020) Table A1 unconfirmed '
                'CRTS candidate J005346.16+223222.4 (not in Camus & Panda 2026 confirmed catalog). Only a 2014 SDSS '
                'low-state spectrum existed in Graham19; DESI (2021) + Palomar NGPS (2026 Sep 23) provide the spectroscopic '
                'turn-on confirmation (Hb/[OIII] up 2.4x).'
            )
        elif len(c_matches) and 'spectroscopic_confirmed' in c_status:
            resolution = 'Catalog-confirmed CLAGN in Camus & Panda (2026); Palomar provides 2026 epoch state verification.'
        elif len(c_matches):
            resolution = 'Reported candidate in Camus & Panda (2026) without high-confidence confirmation.'
        else:
            resolution = 'UNREPORTED IN CHECKED CLAGN CATALOGS (No match in Camus-Panda 2026, Graham 2019, MacLeod 2016/2019, or Zeltyn 2024).'

        rows.append(dict(
            name=obj['name'],
            aliases=obj['aliases'],
            z=float(obj['z']),
            experiment_group=obj['experiment_group'],
            is_deep_dive_lead=bool(obj['name'] in lead_names),
            camus_panda_2026_match=c_status,
            literature_clagn_list_match=l_status,
            zeltyn_2024_match=z_status,
            audit_resolution=resolution,
        ))
    df_lit = pd.DataFrame(rows)
    df_lit.to_csv(DEEP_OUT / 'literature_audit_resolution.csv', index=False)
    return df_lit


def audit_individual_exposures(lead_names, objects_df):
    """Check individual sub-exposures in sep23_data, oct26_data, and oct27_data for all deep-dive leads."""
    night_folders = {
        'sep23': ROOT / 'sep23_data/reduction_20260924/work/20260924/fluxed_exposures_p330e',
        'oct26': ROOT / 'oct26_data/reduction_20261027/work/20261027/fluxed_exposures_feige34',
        'oct27': ROOT / 'oct27_data/reduction_20261028/work/20261028/fluxed_exposures_bd284211',
    }
    z_map = objects_df.set_index('name')['z'].to_dict()
    rows = []
    for night, folder in night_folders.items():
        if not folder.exists():
            continue
        for name in lead_names:
            z = float(z_map.get(name, 0.0))
            hb_obs = 4862.68 * (1.0 + z)
            ch = 'g' if hb_obs < 5750 else ('r' if hb_obs < 7550 else 'i')
            aliases = [name, 'J' + name]
            found = []
            for al in aliases:
                found.extend(sorted((folder / ch).glob(f'*-{al}_NGPS_{ch}_*.fits')))
            if len(found) == 0:
                continue
            for p in found:
                with fits.open(p, memmap=False) as h:
                    d = h[1].data
                    w = np.asarray(d['OPT_WAVE'], float)
                    f = np.asarray(d['OPT_FLAM'], float)
                    e = np.asarray(d['OPT_FLAM_SIG'], float)
                    valid = np.asarray(d['OPT_MASK'], bool) & np.isfinite(e) & (e > 0)
                    sel = np.isfinite(w) & (w > 0)
                    order = np.flatnonzero(sel)[np.argsort(w[sel])]
                    w, f, e, valid = [a[order].copy() for a in [w, f, e, valid]]
                    f[~valid] = np.nan
                    mjd = float(h[0].header.get('MJD', np.nan))
                s = dict(wave=w, flux=f, error=e, z=z, quality=np.zeros(len(w), int), channel=np.full(len(w), ch.upper()), source='NGPS')
                r_hb = A.line_index(s, 'Hb', 0.0)
                rows.append(dict(
                    name=name, night=night, channel=ch.upper(), n_exposures_in_visit=len(found),
                    frame=p.name.split('_')[3].split('-')[0], mjd=mjd,
                    hb_wing_ew=r_hb.get('wing_ew', np.nan),
                    hb_wing_over_oiii=r_hb.get('wing_over_oiii', np.nan),
                    hb_continuum=r_hb.get('continuum', np.nan),
                    valid=bool(r_hb.get('valid', False)),
                    file=str(p.relative_to(ROOT)),
                ))
    df_exp = pd.DataFrame(rows)
    df_exp.to_csv(DEEP_OUT / 'lead_subexposure_audit.csv', index=False)
    return df_exp


def make_figures(obj_summary, df_groups, df_phys, df_transitions, spectra_map):
    """Generate multi-panel diagnostic figures for Step 1 (full sample) and Step 2 (discovery leads)."""
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    # -------------------------------------------------------------------------
    # Figure 1: Full-Sample A/B/C/D Selection Experiment & Manifold Space
    # -------------------------------------------------------------------------
    fig, axes = plt.subplots(1, 3, figsize=(16.5, 5.2), dpi=180)
    hb_sub = obj_summary[obj_summary['hb_ew_total_cont_ratio_baseline'].notna()].copy()

    # Panel A: UMAP manifold space colored by signed log10(Hb EW ratio)
    ax = axes[0]
    m_umap = hb_sub[hb_sub['umap_x'].notna() & hb_sub['umap_y'].notna()].copy()
    sc = ax.scatter(
        m_umap['umap_x'], m_umap['umap_y'],
        c=np.clip(m_umap['log10_hb_ew_ratio'], -0.65, 0.65),
        cmap='coolwarm', vmin=-0.65, vmax=0.65,
        s=np.where(m_umap['hb_both_robust_1p5x'], 95, 42),
        edgecolors=np.where(m_umap['hb_both_robust_1p5x'], 'black', '0.35'),
        linewidths=np.where(m_umap['hb_both_robust_1p5x'], 1.4, 0.5),
        alpha=0.9,
    )
    for r in m_umap[m_umap['hb_both_robust_1p5x']].itertuples(index=False):
        ax.annotate(r.name, (r.umap_x, r.umap_y), xytext=(4, 4), textcoords='offset points', fontsize=7.5, fontweight='bold')
    cb = plt.colorbar(sc, ax=ax, pad=0.02)
    cb.set_label(r'$\log_{10}(\mathrm{EW}_{\mathrm{NGPS}} / \mathrm{EW}_{\mathrm{archive}})$ [H$\beta$]', fontsize=9)
    ax.set_xlabel('WISE W1 UMAP X (frozen pre-run manifold)', fontsize=9.5)
    ax.set_ylabel('WISE W1 UMAP Y (frozen pre-run manifold)', fontsize=9.5)
    ax.set_title('(a) Full Sample in WISE W1 Manifold (N=80)', fontsize=10.5, fontweight='bold')
    ax.grid(True, alpha=0.22)

    # Panel B: Signed log10(Hb EW ratio) vs log10(Hb/[OIII] ratio) with model sensitivity error bars
    ax = axes[1]
    grp_colors = {
        'A: manifold + dated change': '#d95f02',
        'B: manifold without dated trigger': '#1b9e77',
        'C: variability comparison': '#e7298a',
        'D: comparison pool without dated trigger': '#7570b3',
        'outside pre-run review': '#666666',
    }
    ax.axhline(0, color='0.6', ls='--', lw=0.9)
    ax.axvline(0, color='0.6', ls='--', lw=0.9)
    ax.axhspan(-np.log10(1.5), np.log10(1.5), color='0.92', alpha=0.5, zorder=0)
    ax.axvspan(-np.log10(1.5), np.log10(1.5), color='0.92', alpha=0.5, zorder=0)
    for grp, g in hb_sub.groupby('experiment_group'):
        col = grp_colors.get(grp, '#333333')
        xval = g['log10_hb_ew_ratio'].to_numpy()
        yval = g['log10_hb_oiii_ratio'].to_numpy()
        xmin = np.log10(np.maximum(g['hb_ew_total_cont_ratio_min'].to_numpy(), 1e-3))
        xmax = np.log10(np.maximum(g['hb_ew_total_cont_ratio_max'].to_numpy(), 1e-3))
        ymin = np.log10(np.maximum(g['hb_oiii_ratio_min'].to_numpy(), 1e-3))
        ymax = np.log10(np.maximum(g['hb_oiii_ratio_max'].to_numpy(), 1e-3))
        xerr = np.vstack([np.maximum(0, xval - xmin), np.maximum(0, xmax - xval)])
        yerr = np.vstack([np.maximum(0, yval - ymin), np.maximum(0, ymax - yval)])
        short_lbl = grp.split(':')[0] if ':' in grp else 'Outside'
        ax.errorbar(
            xval, yval, xerr=xerr, yerr=yerr, fmt='o', ms=5.5, color=col,
            ecolor=col, elinewidth=0.9, alpha=0.82, label=f"{short_lbl} (n={len(g)})"
        )
    for r in hb_sub[hb_sub['hb_both_robust_1p5x'] | (hb_sub['abs_log10_hb_ew_ratio'] > 0.28)].itertuples(index=False):
        if np.isfinite(r.log10_hb_ew_ratio) and np.isfinite(r.log10_hb_oiii_ratio):
            ax.annotate(r.name, (r.log10_hb_ew_ratio, r.log10_hb_oiii_ratio),
                        xytext=(4, 3), textcoords='offset points', fontsize=7.5, fontweight='bold')
    ax.set_xlim(-0.95, 0.85)
    ax.set_ylim(-0.95, 0.85)
    ax.set_xlabel(r'$\log_{10}(\mathrm{EW}_{\mathrm{H\beta, NGPS}} / \mathrm{EW}_{\mathrm{H\beta, arch}})$', fontsize=9.5)
    ax.set_ylabel(r'$\log_{10}((\mathrm{H\beta/[O\,III]})_{\mathrm{NGPS}} / (\mathrm{H\beta/[O\,III]})_{\mathrm{arch}})$', fontsize=9.5)
    ax.set_title(r'(b) Broad H$\beta$ EW vs. H$\beta$/[O III] Evolution', fontsize=10.5, fontweight='bold')
    ax.legend(loc='upper left', fontsize=7.5, framealpha=0.9)
    ax.grid(True, alpha=0.22)

    # Panel C: Distribution of |log10(EW ratio)| across Groups A, B, D (excluding prior confirmed CLAGNs)
    ax = axes[2]
    clean = hb_sub[hb_sub.known_state_status.ne('catalog-confirmed CLAGN')]
    order = [
        ('A: manifold + dated change', 'Group A\n(Manifold +\nTrigger)', '#d95f02'),
        ('B: manifold without dated trigger', 'Group B\n(Manifold\nOnly)', '#1b9e77'),
        ('D: comparison pool without dated trigger', 'Group D\n(Control\nPool)', '#7570b3'),
    ]
    pos_list, labels_list = [], []
    rng = np.random.default_rng(42)
    for idx, (gname, glabel, gcol) in enumerate(order):
        vals = clean.loc[clean.experiment_group.eq(gname), 'abs_log10_hb_ew_ratio'].dropna().to_numpy()
        pos_list.append(idx)
        labels_list.append(f"{glabel}\n(n={len(vals)})")
        if len(vals):
            bp = ax.boxplot([vals], positions=[idx], widths=0.42, patch_artist=True, showfliers=False)
            for patch in bp['boxes']:
                patch.set_facecolor(gcol)
                patch.set_alpha(0.35)
            for med in bp['medians']:
                med.set_color('black')
                med.set_linewidth(1.8)
            jitter = rng.uniform(-0.12, 0.12, size=len(vals))
            ax.scatter(np.full_like(vals, idx) + jitter, vals, color=gcol, edgecolors='k', linewidths=0.5, s=36, zorder=4)
    ax.axhline(np.log10(1.5), color='crimson', ls=':', lw=1.1, label=r'$1.5\times$ change ($0.176$ dex)')
    ax.axhline(np.log10(2.0), color='darkred', ls='--', lw=1.1, label=r'$2.0\times$ change ($0.301$ dex)')
    ax.set_xticks(pos_list)
    ax.set_xticklabels(labels_list, fontsize=8.5)
    ax.set_ylabel(r'$|\log_{10}(\mathrm{EW}_{\mathrm{H\beta, NGPS}} / \mathrm{EW}_{\mathrm{H\beta, arch}})|$ [dex]', fontsize=9.5)
    ax.set_title('(c) Selection Response by Design Group', fontsize=10.5, fontweight='bold')
    ax.legend(loc='upper right', fontsize=8)
    ax.grid(True, axis='y', alpha=0.22)

    plt.tight_layout()
    fig.savefig(FIG_DIR / 'fig1_full_sample_abcd_experiment.png', dpi=180)
    plt.close(fig)

    # -------------------------------------------------------------------------
    # Figure 2: Top Discovery Leads Joint Multi-Epoch Spectral Evolution Gallery
    # -------------------------------------------------------------------------
    showcase = ['P21422', 'P9506', 'P7947', 'P9694', 'P17881', 'P10381']
    # Add up to 2 newly identified full-sample leads if present
    extra = [
        r.name for r in df_transitions.sort_values('hb_ew_ratio_joint').itertuples(index=False)
        if r.name not in showcase and (r.hb_ew_ratio_joint <= 0.62 or r.hb_ew_ratio_joint >= 1.65)
    ]
    for nm in extra[:2]:
        showcase.append(nm)
    n_show = len(showcase)
    ncols = 2
    nrows = (n_show + 1) // 2
    fig, axes = plt.subplots(nrows, ncols, figsize=(14.0, 3.3 * nrows), dpi=180)
    axes = np.atleast_1d(axes).ravel()

    for idx, name in enumerate(showcase):
        ax = axes[idx]
        eps = sorted(spectra_map.get(name, []), key=lambda s: s['mjd'])
        # Plot continuum-subtracted or [O III]-peak-scaled rest-frame H-beta+[O III] spectra
        t_row = df_transitions[df_transitions.name.eq(name)]
        z_val = eps[0]['z'] if eps else 0.0
        for s in eps:
            x, y, e, lsf, good = PF.prepare(s)
            if good.sum() < 50:
                continue
            # Normalize by [O III] 5007 peak or 5100 continuum for clean visual comparison
            oiii_win = good & (x >= 4995) & (x <= 5020)
            norm = np.nanmax(y[oiii_win]) if oiii_win.sum() >= 3 and np.nanmax(y[oiii_win]) > 0 else np.nanmedian(y[good])
            lbl = f"{s['source']} {int(round(s['mjd']))}"
            col = '#d62728' if s['source'] == 'NGPS' else ('#1f77b4' if s['source'] == 'SDSS' else '#2ca02c')
            lw = 1.45 if s['source'] == 'NGPS' else 1.0
            alp = 0.95 if s['source'] == 'NGPS' else 0.72
            ax.plot(x[good], y[good] / norm, color=col, lw=lw, alpha=alp, label=lbl)
        ax.axvline(4862.68, color='0.4', ls='--', lw=0.8)
        ax.axvline(5008.24, color='0.4', ls=':', lw=0.8)
        ax.set_xlim(4680, 5120)
        title_str = f"{name} (z={z_val:.3f})"
        if len(t_row):
            tr = t_row.iloc[-1]
            title_str += f" | EW: {tr.old_hb_ew_joint:.1f} -> {tr.new_hb_ew_joint:.1f} A ({tr.hb_ew_ratio_joint:.2f}x)"
        ax.set_title(title_str, fontsize=9.5, fontweight='bold')
        ax.set_xlabel(r'Rest Wavelength ($\mathrm{\AA}$)', fontsize=8.5)
        ax.set_ylabel(r'Flux / Peak([O III])', fontsize=8.5)
        ax.legend(loc='upper left', fontsize=7.5, framealpha=0.88)
        ax.grid(True, alpha=0.2)

    for j in range(n_show, len(axes)):
        axes[j].axis('off')
    plt.tight_layout()
    fig.savefig(FIG_DIR / 'fig2_discovery_gallery_joint_fits.png', dpi=180)
    plt.close(fig)

    # -------------------------------------------------------------------------
    # Figure 3: Physical Mechanism Diagnostics (Balmer Decrement & Virial Breathing)
    # -------------------------------------------------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(13.0, 5.2), dpi=180)

    # Panel A: Observed Balmer decrement ratio vs H-beta/[OIII] ratio compared to pure CCM89 dust extinction vector
    ax = axes[0]
    sub_bd = df_transitions[df_transitions['obs_bd_ratio'].notna() & df_transitions['hb_oiii_ratio_joint'].notna()].copy()
    x_grid = np.linspace(-0.9, 0.6, 200)
    # Pure foreground dust extinction curve: log10(BD_new/BD_old) = -((K_HB - K_HA)/K_HB) * log10(Hb_new/Hb_old)
    y_dust = -((K_HB - K_HA) / K_HB) * x_grid
    ax.plot(x_grid, y_dust, 'r--', lw=1.6, label=r'Pure Foreground Dust Extinction ($\mathrm{CCM89\ } R_V=3.1$)')
    ax.axhline(0.0, color='0.35', ls=':', lw=1.2, label='Intrinsic Accretion Continuum Change (constant BD)')
    ax.axvline(0.0, color='0.7', ls='-', lw=0.8)
    x_pts = np.log10(sub_bd['hb_oiii_ratio_joint'])
    y_pts = np.log10(sub_bd['obs_bd_ratio'])
    ax.scatter(x_pts, y_pts, s=70, c='#1f77b4', edgecolors='k', zorder=4)
    for r, xp, yp in zip(sub_bd.itertuples(index=False), x_pts, y_pts):
        ax.annotate(r.name, (xp, yp), xytext=(5, 4), textcoords='offset points', fontsize=8, fontweight='bold')
    ax.set_xlabel(r'$\Delta \log_{10}(\mathrm{H\beta_{\mathrm{broad}} / [O\,III]})$ (NGPS vs. Archive)', fontsize=9.5)
    ax.set_ylabel(r'$\Delta \log_{10}(\mathrm{H\alpha_{\mathrm{broad}} / H\beta_{\mathrm{broad}}})$ (Balmer Decrement)', fontsize=9.5)
    ax.set_title('(a) Accretion State Change vs. Dust Reddening Test', fontsize=10.5, fontweight='bold')
    ax.legend(loc='upper right', fontsize=8)
    ax.grid(True, alpha=0.22)

    # Panel B: Virial breathing (Delta log FWHM_Hb vs Delta log L_Hb_broad)
    ax = axes[1]
    sub_br = df_transitions[
        df_transitions['old_fwhm_hb_kms'].gt(0) & df_transitions['new_fwhm_hb_kms'].gt(0)
        & df_transitions['delta_log_Edd_ratio'].notna()
    ].copy()
    dl = sub_br['delta_log_Edd_ratio'].to_numpy()
    dfw = np.log10(sub_br['new_fwhm_hb_kms'].to_numpy() / sub_br['old_fwhm_hb_kms'].to_numpy())
    x_line = np.linspace(-0.85, 0.65, 100)
    ax.plot(x_line, -0.25 * x_line, 'k--', lw=1.4, label=r'Equilibrium Virial Breathing ($\Delta \log \mathrm{FWHM} = -0.25 \Delta \log \lambda_{\mathrm{Edd}}$)')
    ax.axhline(0, color='0.7', lw=0.8)
    ax.axvline(0, color='0.7', lw=0.8)
    ax.scatter(dl, dfw, s=70, c='#ff7f0e', edgecolors='k', zorder=4)
    for r, xp, yp in zip(sub_br.itertuples(index=False), dl, dfw):
        ax.annotate(r.name, (xp, yp), xytext=(5, 4), textcoords='offset points', fontsize=8, fontweight='bold')
    ax.set_xlabel(r'$\Delta \log_{10}(\lambda_{\mathrm{Edd}})$ (from broad H$\beta$ luminosity)', fontsize=9.5)
    ax.set_ylabel(r'$\Delta \log_{10}(\mathrm{FWHM}_{\mathrm{H\beta}})$', fontsize=9.5)
    ax.set_title(r'(b) Broad-Line Kinematic Breathing Across State Transitions', fontsize=10.5, fontweight='bold')
    ax.legend(loc='upper right', fontsize=8)
    ax.grid(True, alpha=0.22)

    plt.tight_layout()
    fig.savefig(FIG_DIR / 'fig3_balmer_decrement_virial_physics.png', dpi=180)
    plt.close(fig)


def generate_lead_decomposition_atlas():
    """Generate individual 3-row multi-epoch spectral decomposition figures (matching paper_plots.py style)
    and a combined summary figure for the newly surfaced Group B/C leads."""
    from matplotlib.backends.backend_pdf import PdfPages
    from astropy.time import Time

    atlas_dir = DEEP_OUT / 'figures'
    atlas_dir.mkdir(parents=True, exist_ok=True)
    results = pd.read_csv(FULL_OUT / 'hb_fit_results.csv')
    ha_res = pd.read_csv(FULL_OUT / 'ha_fit_results.csv')
    mg_res = pd.read_csv(FULL_OUT / 'mgii_fit_results.csv')

    targets = [
        'P14388', 'P12457', 'P21743', 'P21841', 'P11010', 'P16455', 'P9178',
        'P21422', 'P9506', 'P7947', 'P9694', 'P17881', 'P10381', 'P8544'
    ]
    with PdfPages(DEEP_OUT / 'deep_dive_fit_atlas.pdf') as pdf:
        for name in targets:
            d = results[results.name.eq(name) & results.variant.eq('baseline') & results.status.eq('fit')].sort_values('mjd')
            if len(d) == 0:
                continue
            # If more than 5 epochs (e.g. P14388 has 7), pick earliest, 2 middle high-state epochs, latest archive, and NGPS
            if len(d) > 5:
                idx_keep = sorted(set([0, len(d) // 3, (2 * len(d)) // 3, len(d) - 2, len(d) - 1]))
                d = d.iloc[idx_keep]
            fig, axes = plt.subplots(3, len(d), figsize=(4.3 * len(d), 8.2), squeeze=False, dpi=150)
            for col, r in enumerate(d.itertuples(index=False)):
                npz_path = FULL_OUT / 'fits_hb' / f'{r.epoch_id}__baseline.npz'
                if not npz_path.exists():
                    continue
                a = np.load(npz_path, allow_pickle=True)
                w, y, e, good, parts = [a[k] for k in ['wave_rest', 'flux', 'error', 'good', 'components']]
                model = parts.sum(axis=1)
                scale = max(float(r.total5100), 1e-3)
                ax = axes[0, col]
                ax.plot(w, np.where(good, y, np.nan) / scale, color='0.45', lw=0.6, label='data')
                ax.plot(w, model / scale, 'k', lw=1.0, label='total model')
                for v, color, label in [
                    (parts[:, 0], '#21866a', 'AGN power-law'),
                    (parts[:, 1:3].sum(axis=1), '#b77725', 'BC03 host'),
                    (parts[:, 3], '#8d56ab', 'Fe II'),
                ]:
                    ax.plot(w, v / scale, color=color, lw=1.0, label=label)
                date_str = Time(r.mjd, format='mjd').strftime('%Y-%m-%d')
                ax.set(
                    title=f"{r.source} {date_str} (MJD {int(round(r.mjd))})\nEW(Hβ)={r.hb_ew_total_cont:.1f} Å; Hβ/[O III]={r.hb_oiii:.2f}; χ²/ν={r.reduced_chi2:.2f}",
                    xlim=(4200, 5500),
                    ylim=(-0.2, min(6.5, float(np.nanpercentile(y[good] / scale, 99.5)) * 1.25)),
                )
                ax = axes[1, col]
                continuum = parts[:, :4].sum(axis=1)
                ax.plot(w, np.where(good, y - continuum, np.nan) / scale, color='0.4', lw=0.7, label='data − (cont + Fe II)')
                ax.plot(w, parts[:, 4:].sum(axis=1) / scale, 'k', lw=1.0)
                ax.plot(w, parts[:, 7:9].sum(axis=1) / scale, color='#b83343', lw=1.6, label='broad Hβ')
                ax.plot(w, parts[:, 4:7].sum(axis=1) / scale, color='#247ba0', lw=1.1, label='narrow Hβ + [O III]')
                ax.set(
                    xlim=(4730, 5080),
                    ylim=(-0.2, min(6.5, max(0.8, float(np.nanpercentile((y - continuum)[good] / scale, 99.5)) * 1.25))),
                )
                ax.axhline(0, color='0.7', lw=0.5)
                ax = axes[2, col]
                ax.plot(w, np.where(good, (y - model) / e, np.nan), color='0.35', lw=0.5)
                ax.axhline(0, color='k', lw=0.5)
                ax.axhspan(-3, 3, color='0.5', alpha=0.12)
                ax.set(xlim=(4200, 5500), ylim=(-8, 8), xlabel='Rest vacuum wavelength [Å]')
                if col == 0:
                    axes[0, col].set_ylabel('fλ / continuum(5100 Å)')
                    axes[1, col].set_ylabel('Line flux / continuum(5100 Å)')
                    axes[2, col].set_ylabel('(data − model) / σ')
            axes[0, 0].legend(fontsize=7, ncol=2)
            axes[1, 0].legend(fontsize=7)
            fig.suptitle(f'{name} (z={d.iloc[0].z:.4f}) — Full-Sample Multi-Epoch Hβ + [O III] + Fe II + Host Decomposition', fontsize=11.5, fontweight='bold')
            fig.tight_layout(rect=(0, 0, 1, 0.95))
            fig.savefig(atlas_dir / f'{name}_decomposition.png', dpi=140)
            pdf.savefig(fig)
            plt.close(fig)

    # Also build a compact 5-object side-by-side decomposition comparison figure for the newly surfaced Group B/C leads
    new_leads = ['P14388', 'P12457', 'P21743', 'P21841', 'P11010']
    fig, axes = plt.subplots(len(new_leads), 2, figsize=(13.5, 3.1 * len(new_leads)), dpi=170)
    for row_idx, name in enumerate(new_leads):
        d = results[results.name.eq(name) & results.variant.eq('baseline') & results.status.eq('fit')].sort_values('mjd')
        old_d = d[d.source.ne('NGPS') & d.is_native_error]
        new_d = d[d.source.eq('NGPS')]
        if len(old_d) == 0 or len(new_d) == 0:
            continue
        # For P14388, show its 2021 high-state epoch (59232) or 2025 fading epoch (60702) vs 2026 NGPS
        pair = [old_d.iloc[-1], new_d.iloc[-1]]
        for col_idx, r in enumerate(pair):
            ax = axes[row_idx, col_idx]
            npz_path = FULL_OUT / 'fits_hb' / f'{r.epoch_id}__baseline.npz'
            a = np.load(npz_path, allow_pickle=True)
            w, y, e, good, parts = [a[k] for k in ['wave_rest', 'flux', 'error', 'good', 'components']]
            continuum = parts[:, :4].sum(axis=1)
            scale = max(float(r.total5100), 1e-3)
            ax.plot(w, np.where(good, y, np.nan) / scale, color='0.65', lw=0.7, label='Total data')
            ax.plot(w, continuum / scale, color='#b77725', lw=1.0, ls='--', label='AGN+Host+Fe II cont.')
            ax.plot(w, (continuum + parts[:, 7:9].sum(axis=1)) / scale, color='#b83343', lw=1.7, label='Broad Hβ component')
            ax.plot(w, parts.sum(axis=1) / scale, color='k', lw=1.0, alpha=0.85, label='Total model')
            date_str = Time(r.mjd, format='mjd').strftime('%Y-%m-%d')
            role_lbl = 'Archival Reference' if col_idx == 0 else 'Palomar NGPS (2026)'
            ax.set_title(
                f"{name} (z={r.z:.3f}) — {role_lbl}: {r.source} {date_str}\n"
                f"EW(Hβ)={r.hb_ew_total_cont:.1f} Å | Hβ/[O III]={r.hb_oiii:.2f} | Host(5100)={r.host_fraction5100:.2f}",
                fontsize=9.0, fontweight='bold'
            )
            ax.set_xlim(4650, 5100)
            y_sub = (y[good & (w >= 4650) & (w <= 5100)]) / scale
            ax.set_ylim(0.0, min(6.0, max(1.6, float(np.nanpercentile(y_sub, 99.0)) * 1.18)))
            ax.set_xlabel(r'Rest Wavelength ($\mathrm{\AA}$)', fontsize=8.5)
            ax.set_ylabel(r'$f_\lambda / f_{\mathrm{cont}}(5100\,\mathrm{\AA})$', fontsize=8.5)
            ax.grid(True, alpha=0.22)
            if row_idx == 0 and col_idx == 0:
                ax.legend(loc='upper left', fontsize=7.5, framealpha=0.9)
    plt.tight_layout()
    fig.savefig(FIG_DIR / 'fig4_new_group_bc_leads_decompositions.png', dpi=170)
    plt.close(fig)


def main():
    t0 = time.monotonic()
    DEEP_OUT.mkdir(parents=True, exist_ok=True)
    spectra, objects = FSF.load_all_spectra()
    spectra_by_name = {}
    for s in spectra:
        spectra_by_name.setdefault(s['name'], []).append(s)

    obj_summary = pd.read_csv(FULL_OUT / 'full_sample_object_summary.csv')
    df_groups = pd.read_csv(FULL_OUT / 'abcd_group_experiment_summary.csv')
    hb_env = pd.read_csv(FULL_OUT / 'hb_model_sensitivity_envelope.csv')
    ha_env = pd.read_csv(FULL_OUT / 'ha_model_sensitivity_envelope.csv')
    mgii_env = pd.read_csv(FULL_OUT / 'mgii_model_sensitivity_envelope.csv')
    df_ha = pd.read_csv(FULL_OUT / 'ha_fit_results.csv')

    lead_names = select_deep_dive_leads(obj_summary, hb_env, ha_env, mgii_env)
    print(f"Selected {len(lead_names)} deep-dive discovery & extreme-evolution leads: {lead_names}")

    if (DEEP_OUT / 'joint_multi_epoch_lead_fits.csv').exists():
        df_phys = pd.read_csv(DEEP_OUT / 'joint_multi_epoch_lead_fits.csv')
    else:
        joint_rows = []
        for name in lead_names:
            eps = spectra_by_name.get(name, [])
            if not eps:
                continue
            z = float(eps[0]['z'])
            rows, theta = joint_multi_epoch_hb_fit(eps, z)
            joint_rows.extend(rows)
        df_joint = pd.DataFrame(joint_rows)
        df_phys = compute_virial_and_extinction(df_joint, df_ha)
        df_phys.to_csv(DEEP_OUT / 'joint_multi_epoch_lead_fits.csv', index=False)

    df_transitions = evaluate_lead_transitions(df_phys, hb_env)
    df_transitions.to_csv(DEEP_OUT / 'lead_physical_transitions.csv', index=False)

    df_lit = audit_literature_cross_matches(obj_summary, lead_names)
    df_exp = audit_individual_exposures(lead_names, obj_summary)
    make_figures(obj_summary, df_groups, df_phys, df_transitions, spectra_by_name)
    generate_lead_decomposition_atlas()

    print("\n=== Deep-Dive Joint Multi-Epoch Transitions Summary ===")
    cols = [
        'name', 'z', 'baseline_rest_years', 'old_hb_ew_joint', 'new_hb_ew_joint',
        'hb_ew_ratio_joint', 'hb_oiii_ratio_joint', 'ha_ew_ratio', 'obs_bd_ratio',
        'pred_bd_ratio_if_pure_dust', 'log_MBH_adopted', 'delta_log_Edd_ratio'
    ]
    print(df_transitions[cols].to_string(index=False))
    print(f"\nCompleted Step 2 deep-dive analysis and figures in {time.monotonic()-t0:.1f}s.")


if __name__ == '__main__':
    main()

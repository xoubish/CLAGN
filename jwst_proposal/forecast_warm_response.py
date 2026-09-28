"""Forecast the 2028 warm-dust flux of each target from its own recorded history.

For every target the W1 history (AllWISE 2010, NEOWISE 2014-2024) is taken as the
nuclear heating history, normalised to the 2010 AllWISE epoch. The optically thin
shell model of make_fig2_diagnostics.py is then driven by that history under three
hypotheses and evaluated at the AllWISE W3 epoch (2010) and at the MIRI epoch:

  instant   fixed dust distribution, no light-travel delays (no memory)
  delayed   fixed dust distribution with light-travel delays (memory)
  evolving  delayed, and the inner boundary follows L^1/2 as illumination arrives

The observable is log10 of the nuclear 12-um (observed W3) flux at the MIRI epoch
relative to 2010. AllWISE 2010 anchors each object, so no equilibrium covering
factor is needed; the hypotheses differ only through the history. Assumptions:
constant heating before the first W1 point and after the last one; W1 host light
is not removed (amplitudes are lower bounds); the model unit R_in,Si/c equals
6.0 K-band lags (the model's silicate/graphite inner-radius ratio) with the
K-band lag from the sample table (Koshida et al. 2014 relation); rest-frame delays
are stretched by (1+z) on the observed time axis. A constant host
contribution to W1 is removed before driving the model: the 2010 host fraction is
estimated from the AllWISE W1-W2 colour by mixing an AGN colour (W1-W2 = 1.0 Vega)
with a stellar host colour (0.05), clipped to [0, 0.9]; the light-crossing scale is
also run at 0.5x and 2x the baseline as a sensitivity range.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from astropy.cosmology import Planck18

import make_fig2_diagnostics as shell

HERE = Path(__file__).resolve().parent
INPUTS = HERE/'inputs'
MIRI_EPOCH = 2028.0
W1_ZERO_JY, W2_ZERO_JY = 309.540, 171.787
AGN_W1W2, HOST_W1W2 = 1.0, .05          # Vega colours for the two-component host estimate
SCALE_FACTORS = {'baseline': 1., 'short': .5, 'long': 2.}
SIGMA_DECOMP_DEX = .03                   # MRS nuclear 12-um flux vs AllWISE W3: decomposition + random calibration
SIGMA_HOST_FRACTION = .15                # uncertainty of the W1 host fraction (colour-based estimate)
DRIVER_GAMMA = None                      # d ln F(W1) / d ln L of the model's fixed hot dust; set in main()
N_DELAY = 48
R_GRID = np.geomspace(.1, 100., 1500)
DLOGR = np.log(R_GRID[1]/R_GRID[0])
RIN_SI_OVER_RIN_GRA = None  # filled from the model


def decimal_year(mjd):
    return 2000. + (np.asarray(mjd, float) - 51544.5)/365.25


def koshida_lag_years(r_mag, z):
    """Sample-table recipe: log10 tau_K [days] = -2.11 - 0.2 M_V, M_V ~ r - DM(z)."""
    dm = Planck18.distmod(z).value
    return 10**(-2.11 - .2*(r_mag - dm))/365.25


def host_fraction(aw):
    """Host share of the 2010 W1 flux from the AllWISE W1-W2 colour (two-component mix)."""
    observed = (W2_ZERO_JY/W1_ZERO_JY)*10**(.4*(aw['w1mpro']-aw['w2mpro']))
    agn = (W2_ZERO_JY/W1_ZERO_JY)*10**(.4*AGN_W1W2)
    host = (W2_ZERO_JY/W1_ZERO_JY)*10**(.4*HOST_W1W2)
    f_agn = (observed-host)/(agn-host)
    return float(np.clip(1-f_agn, 0., .9))


def history(sample_row, binned, allwise, host=0.):
    """Return (years, L/L2010) from AllWISE + NEOWISE W1 with a constant host removed."""
    aw = allwise.loc[allwise['id'].eq(sample_row['id'])].iloc[0]
    t0 = float(decimal_year(aw['w1mjdmean']))
    f0 = W1_ZERO_JY*1e3*10**(-.4*aw['w1mpro'])
    neo = binned.loc[binned['id'].eq(sample_row['id']) & binned['band'].eq('W1')].sort_values('mjd')
    t = np.r_[t0, decimal_year(neo['mjd'].values)]
    f = np.r_[f0, W1_ZERO_JY*1e3*10**(-.4*neo['mag'].values)]
    nuclear = np.clip((f/f0-host)/(1-host), .05, None)
    if DRIVER_GAMMA:
        nuclear = nuclear**(1./DRIVER_GAMMA)   # invert the hot-dust response to the heating luminosity
    return t, nuclear, t0


def model_flux(times, hist_t, hist_l, rest_wave_um, rin_si_years, hypothesis, z=0.):
    """Nuclear F_nu (model units) at rest_wave_um for the given observed-frame epochs.

    rin_si_years is the rest-frame light-crossing time of the initial silicate inner
    radius; delays are multiplied by (1+z) to place them on the observed time axis.
    """
    q = {g['name']: float(np.interp(rest_wave_um, shell.WAVE, g['q'])) for g in shell.GRAINS}
    out = np.zeros(len(times))
    if hypothesis == 'instant':
        delays = np.zeros((len(R_GRID), 1))
    else:
        k = (np.arange(N_DELAY)+.5)/N_DELAY
        delays = 2*R_GRID[:, None]*k[None, :]*rin_si_years*(1+z)  # observed-frame years
    weight_r = R_GRID*DLOGR
    for i, t in enumerate(times):
        lum = np.interp(t - delays, hist_t, hist_l)                  # constant outside the record
        total = 0.
        for g in shell.GRAINS:
            power = shell.HEATING*g['q_heat']/R_GRID[:, None]**2*lum
            temp = np.exp(np.interp(np.log(power), g['log_cool'], g['log_t']))
            planck = shell.planck_nu(np.array([rest_wave_um]), temp)[..., 0]
            if hypothesis == 'evolving':
                mask = R_GRID[:, None] >= g['inner_radius']*np.sqrt(lum)
            else:
                mask = (R_GRID >= g['inner_radius'])[:, None]
            contrib = np.mean(planck*mask, axis=1)*weight_r*g['mass_fraction']/g['rho']*q[g['name']]
            total += contrib.sum()
        out[i] = total
    return out


def hot_response_exponent():
    """d ln F(rest 3.0 um) / d ln L for the fixed distribution in equilibrium."""
    f = [model_flux(np.array([0.]), np.array([-1., 1.]), np.array([r, r]), 3.0, 1., 'instant')[0] for r in (.5, 2.)]
    return float(np.log(f[1]/f[0])/np.log(4.))


def main():
    global RIN_SI_OVER_RIN_GRA, DRIVER_GAMMA
    grains = {g['name']: g for g in shell.GRAINS}
    RIN_SI_OVER_RIN_GRA = grains['Sil_21.gz']['inner_radius']/grains['Gra_21.gz']['inner_radius']
    DRIVER_GAMMA = hot_response_exponent()
    print(f'hot-dust response exponent gamma = {DRIVER_GAMMA:.3f}; heating history = (nuclear W1)^(1/gamma)')
    sample = pd.read_csv(INPUTS/'jwst_sample_cycle6.csv')
    binned = pd.read_csv(HERE/'review/history_audit/binned_measurements.csv')
    allwise = pd.read_csv(INPUTS/'sample_allwise_psd.csv')
    results = []
    series = {}
    for _, row in sample.iterrows():
        tau_k = row['tau_dust_yr']
        if not np.isfinite(tau_k):
            tau_k = koshida_lag_years(row['r_mag'], row['z'])
        rin_si_years = RIN_SI_OVER_RIN_GRA*tau_k
        aw = allwise.loc[allwise['id'].eq(row['id'])].iloc[0]
        host = host_fraction(aw)
        t, l, t0 = history(row, binned, allwise, host)
        rest = 12./(1+row['z'])
        epochs = np.array([t0, MIRI_EPOCH])
        # Host-fraction uncertainty propagated through the no-memory prediction.
        inst = []
        for h in [max(host-SIGMA_HOST_FRACTION, 0.), min(host+SIGMA_HOST_FRACTION, .9)]:
            th, lh, _ = history(row, binned, allwise, h)
            v = model_flux(epochs, th, lh, rest, rin_si_years, 'instant', row['z'])
            inst.append(np.log10(v[1]/v[0]))
        sigma_host = float(abs(inst[1]-inst[0])/2)
        sigma_w3 = float(np.sqrt((.4*aw['w3sigmpro'])**2 + SIGMA_DECOMP_DEX**2 + sigma_host**2))
        entry = dict(id=row['id'], target=row['target'], z=row['z'], family=row['family'],
                     tau_k_years=float(tau_k), rin_si_years=float(rin_si_years),
                     w1_host_fraction_2010=host, w1w2_2010=float(aw['w1mpro']-aw['w2mpro']),
                     w3_mag_2010=float(aw['w3mpro']), w3_sigma_mag=float(aw['w3sigmpro']), sigma_dex=sigma_w3,
                     sigma_host_dex=sigma_host,
                     anchor_year=float(t0), last_w1_year=float(t[-1]),
                     nuclear_ratio_last_over_2010=float(l[-1]), nuclear_min_ratio=float(l.min()),
                     nuclear_max_ratio=float(l.max()), rest_wave_um=float(rest))
        for name, factor in SCALE_FACTORS.items():
            pred = {h: model_flux(epochs, t, l, rest, factor*rin_si_years, h, row['z'])
                    for h in ['instant', 'delayed', 'evolving']}
            suffix = '' if name == 'baseline' else f'_{name}'
            for h, v in pred.items():
                entry[f'log_ratio_{h}{suffix}'] = float(np.log10(v[1]/v[0]))
            entry[f'delayed_minus_instant{suffix}'] = entry[f'log_ratio_delayed{suffix}']-entry[f'log_ratio_instant{suffix}']
            entry[f'evolving_minus_delayed{suffix}'] = entry[f'log_ratio_evolving{suffix}']-entry[f'log_ratio_delayed{suffix}']
            if name == 'baseline':
                baseline = pred
        pred = baseline
        # Sign-corrected signals: positive when the warm flux lags the nuclear change.
        sign = -np.sign(entry['log_ratio_instant']) if entry['log_ratio_instant'] != 0 else 1.
        entry['sign'] = float(sign)
        entry['memory_signal'] = sign*entry['delayed_minus_instant']
        entry['evolution_signal'] = sign*entry['evolving_minus_delayed']
        results.append(entry)
        grid = np.arange(2009.5, MIRI_EPOCH+.01, .25)
        series[row['id']] = dict(years=grid.tolist(), history_t=t.tolist(), history_l=l.tolist(),
                                 **{h: (model_flux(grid, t, l, rest, rin_si_years, h, row['z'])/pred[h][0]).tolist()
                                    for h in ['instant', 'delayed', 'evolving']})
        print(f"{row['id']} z={row['z']:.3f} tauK={tau_k:.2f} R_Si/c={rin_si_years:.1f} yr host={host:.2f} nuc last/2010={l[-1]:.2f} "
              f"log(2028/2010): instant {entry['log_ratio_instant']:+.3f} delayed {entry['log_ratio_delayed']:+.3f} "
              f"evolving {entry['log_ratio_evolving']:+.3f}  |D-I| {abs(entry['delayed_minus_instant']):.3f} |E-D| {abs(entry['evolving_minus_delayed']):.3f}")
    table = pd.DataFrame(results)
    table.to_csv(INPUTS/'warm_response_forecast.csv', index=False)
    (INPUTS/'warm_response_forecast_series.json').write_text(json.dumps(dict(
        model=__doc__, rin_si_over_rin_gra=RIN_SI_OVER_RIN_GRA, miri_epoch=MIRI_EPOCH,
        n_delay=N_DELAY, r_grid=[float(R_GRID[0]), float(R_GRID[-1]), len(R_GRID)], series=series)))
    err = np.sqrt((table['sigma_dex']**2).sum())/len(table)
    summary = {}
    for name in SCALE_FACTORS:
        suffix = '' if name == 'baseline' else f'_{name}'
        mem = (table['sign']*table[f'delayed_minus_instant{suffix}']).mean()
        evo = (table['sign']*table[f'evolving_minus_delayed{suffix}']).mean()
        summary[name] = dict(mean_memory_signal_dex=float(mem), mean_evolution_signal_dex=float(evo),
                             mean_error_dex=float(err), memory_sigma=float(mem/err), evolution_sigma=float(evo/err),
                             n_memory_over_3sigma=int((table[f'delayed_minus_instant{suffix}'].abs() > 3*table['sigma_dex']).sum()),
                             n_evolution_over_3sigma=int((table[f'evolving_minus_delayed{suffix}'].abs() > 3*table['sigma_dex']).sum()))
        print(f"{name:8s}: mean memory {mem:+.3f} dex ({mem/err:.1f} sigma), mean evolution {evo:+.3f} ({evo/err:.1f} sigma); "
              f"per-object >3sigma: memory {summary[name]['n_memory_over_3sigma']}, evolution {summary[name]['n_evolution_over_3sigma']}; mean error {err:.3f}")
    (INPUTS/'warm_response_forecast_summary.json').write_text(json.dumps(dict(
        summary=summary, host_estimate=f'AllWISE W1-W2 two-component mix, AGN {AGN_W1W2}, host {HOST_W1W2} (Vega), clipped to [0,0.9]',
        sigma_per_object='sqrt((0.4*w3sigmpro)^2 + 0.03^2 + sigma_host^2) dex, sigma_host from +-0.15 in the W1 host fraction through the no-memory prediction; common-mode calibration cancels in the sign-weighted mean',
        scale_factors=SCALE_FACTORS, rin_si_over_rin_gra=RIN_SI_OVER_RIN_GRA, driver_gamma=DRIVER_GAMMA,
        n_pair_over_3sigma_measurement_only=int((table[['log_ratio_instant','log_ratio_delayed','log_ratio_evolving']].max(axis=1)
            -table[['log_ratio_instant','log_ratio_delayed','log_ratio_evolving']].min(axis=1)
            > 3*np.sqrt(table['sigma_dex']**2-table['sigma_host_dex']**2)).sum()),
        n_pair_over_3sigma_total=int((table[['log_ratio_instant','log_ratio_delayed','log_ratio_evolving']].max(axis=1)
            -table[['log_ratio_instant','log_ratio_delayed','log_ratio_evolving']].min(axis=1) > 3*table['sigma_dex']).sum())), indent=2))


if __name__ == '__main__':
    main()

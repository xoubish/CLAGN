"""Scenario grid for the warm-response forecast: which assumptions move the signal.

Each scenario changes one ingredient of forecast_warm_response.py (or all of them):
  host+     W1 host fraction at its upper bound (+0.15, capped at 0.9)
  radius x3 warm-dust light-crossing scale three times the baseline
  driver    invert W1 with the model's own fixed-dust hot response, F_W1 ~ L^gamma
  trend     continue the last three observed years of the nuclear W1 trend to the MIRI
            epoch (log-linear, capped at a factor 2), instead of holding it constant
Outputs inputs/warm_response_scenarios.csv with ensemble statistics per scenario.
"""
import json
import numpy as np
import pandas as pd

import forecast_warm_response as F
import make_fig2_diagnostics as shell

HERE, INPUTS = F.HERE, F.INPUTS


def hot_response_exponent():
    """d ln F(3.4 um) / d ln L for the fixed distribution in equilibrium (rest 3.0 um)."""
    lo, hi = .5, 2.
    f = [F.model_flux(np.array([0.]), np.array([-1., 1.]), np.array([r, r]), 3.0, 1., 'instant')[0] for r in (lo, hi)]
    return float(np.log(f[1]/f[0])/np.log(hi/lo))


def run(scenario, gamma):
    sample = pd.read_csv(INPUTS/'jwst_sample_cycle6.csv')
    binned = pd.read_csv(HERE/'review/history_audit/binned_measurements.csv')
    allwise = pd.read_csv(INPUTS/'sample_allwise_psd.csv')
    base = pd.read_csv(INPUTS/'warm_response_forecast.csv').set_index('id')
    rows = []
    for _, row in sample.iterrows():
        tau = row['tau_dust_yr'] if np.isfinite(row['tau_dust_yr']) else F.koshida_lag_years(row['r_mag'], row['z'])
        aw = allwise.loc[allwise['id'].eq(row['id'])].iloc[0]
        host = F.host_fraction(aw)
        if 'host+' in scenario:
            host = min(host+F.SIGMA_HOST_FRACTION, .9)
        t, l, t0 = F.history(row, binned, allwise, host)
        # ('driver' inversion is now part of the baseline via F.DRIVER_GAMMA)
        if 'trend' in scenario:
            recent = t >= t[-1]-3.
            if recent.sum() >= 3:
                slope = np.polyfit(t[recent], np.log10(l[recent]), 1)[0]
                extrap = np.clip(10**(slope*(F.MIRI_EPOCH-t[-1])), .5, 2.)*l[-1]
                t = np.r_[t, F.MIRI_EPOCH]; l = np.r_[l, extrap]
        scale = F.RIN_SI_OVER_RIN_GRA*tau*(3. if 'radius' in scenario else 1.)
        rest = 12./(1+row['z']); epochs = np.array([t0, F.MIRI_EPOCH])
        pred = {h: np.log10(np.divide(*F.model_flux(epochs, t, l, rest, scale, h, row['z'])[::-1]))
                for h in ['instant', 'delayed', 'evolving']}
        sign = -np.sign(pred['instant']) if pred['instant'] != 0 else 1.
        rows.append(dict(id=row['id'], scenario=scenario, sigma=base.loc[row['id'], 'sigma_dex'],
                         lag=sign*(pred['delayed']-pred['instant']), evolution=sign*(pred['evolving']-pred['delayed']),
                         spread=max(pred.values())-min(pred.values()), **{f'log_{h}': v for h, v in pred.items()}))
    return pd.DataFrame(rows)


def main():
    grains = {g['name']: g for g in shell.GRAINS}
    F.RIN_SI_OVER_RIN_GRA = grains['Sil_21.gz']['inner_radius']/grains['Gra_21.gz']['inner_radius']
    gamma = hot_response_exponent()
    F.DRIVER_GAMMA = gamma   # baseline now includes the inversion; the driver scenario is then identical
    print(f'model hot-dust response exponent gamma = {gamma:.2f} (L ~ W1^(1/gamma))')
    scenarios = ['baseline', 'host+', 'radius x3', 'trend', 'host+ radius x3 trend']
    tables, summary = [], []
    for sc in scenarios:
        tb = run(sc, gamma); tables.append(tb)
        err = np.sqrt((tb['sigma']**2).sum())/len(tb)
        summary.append(dict(scenario=sc, mean_lag=tb['lag'].mean(), lag_sigma=tb['lag'].mean()/err,
                            mean_evolution=tb['evolution'].mean(), evolution_sigma=tb['evolution'].mean()/err,
                            n_lag_gt_3sigma=int((tb['lag'].abs() > 3*tb['sigma']).sum()),
                            n_evolution_gt_3sigma=int((tb['evolution'].abs() > 3*tb['sigma']).sum()),
                            n_spread_gt_3sigma=int((tb['spread'] > 3*tb['sigma']).sum()),
                            median_spread=tb['spread'].median()))
        s = summary[-1]
        print(f"{sc:32s} lag {s['mean_lag']:+.3f} ({s['lag_sigma']:.1f}σ)  evolution {s['mean_evolution']:+.3f} ({s['evolution_sigma']:.1f}σ)  "
              f"per-object >3σ: lag {s['n_lag_gt_3sigma']:2d}  evolution {s['n_evolution_gt_3sigma']:2d}  any-pair spread {s['n_spread_gt_3sigma']:2d}  median spread {s['median_spread']:.3f}")
    pd.concat(tables).to_csv(INPUTS/'warm_response_scenarios.csv', index=False)
    (INPUTS/'warm_response_scenarios_summary.json').write_text(json.dumps(dict(gamma=gamma, summary=summary), indent=2))


if __name__ == '__main__':
    main()

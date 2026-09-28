"""Run the two-epoch warm-response forecast over the eligible candidate pool.

Same model and assumptions as forecast_warm_response.py, applied to every pool object
with z<0.3, NEOWISE W1 amplitude >= 0.37 mag, >= 15 visits and an AllWISE match.
The K-band lag uses the Koshida relation with the SDSS r magnitude. Outputs
inputs/pool_forecast.csv with the three predictions, the spread and the measurement
error per object, so a separation-optimised sample can be compared with the chosen 24.
"""
import sys
import numpy as np
import pandas as pd

import forecast_warm_response as F
import make_fig2_diagnostics as shell

S = '/private/tmp/claude-502/-Users-shemmati-Desktop-CLAGN/3bcb6145-cee1-4816-b434-a34bd415e899/scratchpad'


def main(pool_csv=f'{S}/pool_eligible.csv', allwise_csv=f'{S}/pool_allwise.csv', visits_csv='../data/neowise_visits_poolall.csv'):
    grains = {g['name']: g for g in shell.GRAINS}
    F.RIN_SI_OVER_RIN_GRA = grains['Sil_21.gz']['inner_radius']/grains['Gra_21.gz']['inner_radius']
    F.DRIVER_GAMMA = F.hot_response_exponent()
    pool = pd.read_csv(pool_csv)
    allwise = pd.read_csv(allwise_csv).drop_duplicates('name').set_index('name')
    visits = pd.read_csv(F.HERE/visits_csv)
    visits = visits[visits['name'].isin(pool['name'])]
    rows = []
    for i, r in pool.iterrows():
        if r['name'] not in allwise.index:
            continue
        aw = allwise.loc[r['name']]
        if not np.isfinite([aw['w1mpro'], aw['w2mpro'], aw['w3mpro'], aw['w3sigmpro']]).all():
            continue
        v = visits[visits['name'].eq(r['name'])].sort_values('mjd')
        if len(v) < 15:
            continue
        host = F.host_fraction(aw)
        t0 = float(F.decimal_year(aw['w1mjdmean']))
        f0 = F.W1_ZERO_JY*1e3*10**(-.4*aw['w1mpro'])
        t = np.r_[t0, F.decimal_year(v['mjd'].values)]
        f = np.r_[f0, F.W1_ZERO_JY*1e3*10**(-.4*v['w1'].values)]
        nuclear = np.clip((f/f0-host)/(1-host), .05, None)**(1./F.DRIVER_GAMMA)
        tau = F.koshida_lag_years(r['psfmag_r']-r.get('extinction_r', 0.), r['z'])
        rin = F.RIN_SI_OVER_RIN_GRA*tau
        rest = 12./(1+r['z']); epochs = np.array([t0, F.MIRI_EPOCH])
        pred = {h: np.log10(np.divide(*F.model_flux(epochs, t, nuclear, rest, rin, h, r['z'])[::-1]))
                for h in ['instant', 'delayed', 'evolving']}
        # host-fraction term through the no-memory prediction
        inst = []
        for hh in [max(host-.15, 0.), min(host+.15, .9)]:
            nh = np.clip((f/f0-hh)/(1-hh), .05, None)**(1./F.DRIVER_GAMMA)
            vv = F.model_flux(epochs, t, nh, rest, rin, 'instant', r['z']); inst.append(np.log10(vv[1]/vv[0]))
        sig_host = abs(inst[1]-inst[0])/2
        sig_meas = float(np.hypot(.4*aw['w3sigmpro'], .03))
        vals = np.array(list(pred.values()))
        rows.append(dict(name=r['name'], ra=r['ra'], dec=r['dec'], z=r['z'], r_mag=r['psfmag_r'], tau_k=tau,
                         w1_amp_mag=r['w1_amp_visits'], w3_mjy=1e3*31.674*10**(-.4*aw['w3mpro']), w3_sigma_mag=aw['w3sigmpro'],
                         host=host, nuclear_last=float(nuclear[-1]), sigma_meas=sig_meas, sigma_host=sig_host,
                         sigma_total=float(np.hypot(sig_meas, sig_host)), spread=float(vals.max()-vals.min()),
                         **{f'log_{h}': float(val) for h, val in pred.items()}))
        if len(rows) % 200 == 0:
            print(len(rows), 'done', flush=True)
    out = pd.DataFrame(rows)
    out['signif_meas'] = out['spread']/out['sigma_meas']; out['signif_total'] = out['spread']/out['sigma_total']
    out.to_csv(F.INPUTS/'pool_forecast.csv', index=False)
    print('pool forecast rows', len(out))


if __name__ == '__main__':
    main(*sys.argv[1:])

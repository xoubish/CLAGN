"""Plan September promotions from the frozen backup pool at continuum S/N=5.

Writes a private, reviewable plan; does not render pages or change telescope CSVs.
The archived eight-target packet fixes the eligible promotion pool reproducibly.
"""
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import hashlib
import importlib
import json
import math
import re

import astropy.units as u
import numpy as np
import pandas as pd
from scipy.optimize import Bounds, LinearConstraint, milp
from scipy.sparse import coo_matrix

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'data/reselection_2026-09-20'
DEST = ROOT / 'observing/sep23'
SEED = ROOT / 'archive/2026-09-20/before_snr5_packet/sep23/packet.json'
CACHE = OUT / 'sep23_snr5_etc'
TZ = 'America/Los_Angeles'
SETTINGS = dict(goal_snr=5, exposures=2, sky_V=18., seeing_zenith_500nm=1.3,
                binspect=2, binspat=2, slit_arcsec=1., minimum_each_seconds=60,
                round_each_seconds=30, maximum_each_seconds=900,
                overhead_minutes=10, slot_round_minutes=5,
                normalization='Archival local Hbeta continuum; no brightness forecast',
                extraction='Central slice only; optimal point-source extraction')


def exposure(target):
    model = importlib.import_module('22_september_etc')
    ref = model.reference(target)
    if ref is None:
        return target['name'], {}
    mjd, mag, private, flux = ref
    signature = hashlib.sha256(json.dumps([SETTINGS, target['z'], ref], sort_keys=True).encode()).hexdigest()
    path = CACHE / f"{target['name']}.json"
    if path.exists():
        cached = json.loads(path.read_text())
        if cached['signature'] == signature:
            return target['name'], cached['plans']
    for ch, rn, scale in zip(model.CFG.channels, [2.8, 7.8, 3.7, 4.6], [.193, .189, .186, .186]):
        model.CFG.readnoise[ch] = rn*u.count/u.pix
        model.CFG.platescale[ch] = scale*u.arcsec/u.pix
    wave = 486.27*(1+target['z'])
    channels = [ch for ch in ['R', 'I', 'G', 'U']
                if model.CFG.channelRange[ch][0].to_value(u.nm) <= wave-4
                and model.CFG.channelRange[ch][1].to_value(u.nm) >= wave+4]
    if not channels:
        return target['name'], {}
    channel = channels[0]
    plans = {}
    for tier, x in [('preferred', 1.5), ('extended', 1.8)]:
        cmd = [channel, str(wave-4), str(wave+4), 'SNR', str(5/np.sqrt(2)),
               '-slit', 'SET', '1', '-binspect', '2', '-binspat', '2',
               '-seeing', '1.3', '500', '-airmass', str(x), '-skymag', '18',
               '-mag', str(mag), '-magsystem', 'AB', '-magfilter', 'match', '-noslicer']
        args = model.ETC.parser.parse_args(cmd)
        model.ETC.check_inputs_add_units(args)
        solved = float(model.ETC.main(args, quiet=True)['exptime'].to_value(u.s))
        each = max(60, math.ceil(solved/30)*30)
        if each > 900:
            continue
        def forward(magnitude):
            trial = cmd.copy()
            trial[3:5] = ['EXPTIME', str(each)]
            trial[trial.index('-mag')+1] = str(magnitude)
            args = model.ETC.parser.parse_args(trial)
            model.ETC.check_inputs_add_units(args)
            return float(model.ETC.main(args, quiet=True)['SNR'].value)*np.sqrt(2)
        achieved = forward(mag)
        assert achieved >= 4.99, (target['name'], tier, achieved)
        duration = 5*math.ceil((2*each/60+10)/5)
        plans[tier] = dict(exposures=2, seconds_each=each, integration_minutes=2*each/60,
                           visit_minutes=duration, airmass=x, goal_snr=5,
                           predicted_snr=achieved, predicted_snr_fainter0p5=forward(mag+.5),
                           reference_mjd=mjd, reference_private=bool(private), continuum_AB=mag,
                           channel=channel, low_nm=wave-4, high_nm=wave+4,
                           sky_V=18., seeing_zenith=1.3, overhead_minutes=10, status='scenario')
    path.write_text(json.dumps(dict(signature=signature, plans=plans), indent=2))
    return target['name'], plans


def main():
    CACHE.mkdir(parents=True, exist_ok=True)
    old = json.loads(SEED.read_text())
    original = {v['name'] for v in old['primaries']}
    former_backups = {v['name'] for v in old['backups']}
    targets = pd.read_csv(OUT/'compact_review_objects.csv').set_index('name', drop=False)
    science = pd.read_csv(OUT/'three_night_review/science_and_sensitivity.csv').set_index('name')
    opts = json.loads((OUT/'airmass_options.json').read_text())
    public = json.loads(re.search(r'<script id="candidate-data" type="application/json">(.*?)</script>',
                                  (ROOT/'docs/index.html').read_text(), re.S).group(1))
    public_names = {t['name'] for t in public['targets']}
    september = [n for n in targets.index if any(w['night']=='sep23' for w in opts['windows'].get(n, []))]
    with ProcessPoolExecutor(max_workers=4) as executor:
        plans = dict(executor.map(exposure, targets.loc[september].to_dict('records')))
    start = pd.Timestamp('2026-09-23 20:16', tz=TZ)
    end = pd.Timestamp('2026-09-24 00:08', tz=TZ)
    minutes = int((end-start).total_seconds()/60)
    admitted = []; excluded = []
    for name in sorted(original | former_backups):
        t = targets.loc[name]; s = science.loc[name]; why = []
        if name not in original and t.pool_role != 'manifold':
            why.append('Comparison reserve; prioritize manifold candidates for promotions')
        if name not in public_names or any(p['reference_private'] for p in plans.get(name, {}).values()):
            why.append('Retained as a local backup; current shared primary packet requires public identity and exposure reference')
        if t.field_status != 'clear':
            why.append('Field requires further review')
        if pd.notna(s.ztf_r_latest180_mag) and s.ztf_r_latest180_mag >= 19:
            why.append('Latest cached r does not meet brightness criterion')
        if not plans.get(name):
            why.append('No acceptable two-exposure ETC plan')
        if why:
            excluded.append(dict(name=name, reasons=why))
        else:
            admitted.append(name)
    assert original <= set(admitted), 'An original primary needs explicit review'
    choices = []
    for name in admitted:
        t = targets.loc[name]; s = science.loc[name]
        window = next(w for w in opts['windows'][name] if w['night']=='sep23')
        # This is a scheduling utility, not a calibrated changing-look probability.
        score = 20+float(s.review_order_score)+5*bool(t.balmer_pair_in_range)+5*float(t.manifold_cl_neighbor_fraction)
        for tier in ['preferred', 'extended']:
            plan = plans[name].get(tier)
            if plan is None:
                continue
            duration = plan['visit_minutes']
            ranges = [(pd.Timestamp(r['start_utc'], tz='UTC'), pd.Timestamp(r['end_utc'], tz='UTC'))
                      for r in window['tier_ranges'][tier]]
            for offset in range(minutes-duration+1):
                a = start+pd.Timedelta(minutes=offset); b = a+pd.Timedelta(minutes=duration)
                if not any(lo<=a and b<=hi for lo, hi in ranges):
                    continue
                value = score-(2 if tier=='extended' else 0)-offset*.00001
                choices.append(dict(name=name, tier=tier, offset=offset, duration=duration,
                                    start_pdt=a.strftime('%Y-%m-%d %H:%M'), utility=value))
    rr=[]; cc=[]
    name_index = {n:minutes+i for i,n in enumerate(admitted)}
    for j, c in enumerate(choices):
        rr.extend(range(c['offset'], c['offset']+c['duration'])); cc.extend([j]*c['duration'])
        rr.append(name_index[c['name']]); cc.append(j)
    matrix = coo_matrix((np.ones(len(rr)), (rr,cc)), shape=(minutes+len(admitted),len(choices))).tocsc()
    lower = np.zeros(matrix.shape[0]); upper = np.ones(matrix.shape[0])
    for name in original:
        lower[name_index[name]]=1
    result = milp(-np.array([c['utility'] for c in choices]), integrality=np.ones(len(choices)),
                  bounds=Bounds(0,1), constraints=LinearConstraint(matrix, lower, upper),
                  options=dict(time_limit=45, mip_rel_gap=.002))
    assert result.x is not None, result.message
    selected = [choices[i] for i in np.flatnonzero(result.x>.5)]
    selected.sort(key=lambda c:c['offset'])
    assert original <= {c['name'] for c in selected}
    assert len({c['name'] for c in selected}) == len(selected)
    for a,b in zip(selected,selected[1:]):
        assert a['offset']+a['duration']<=b['offset']
    output = dict(settings=SETTINGS, plans=plans, sequence=selected,
                  original_primaries=sorted(original), former_backups=sorted(former_backups),
                  promoted=[c['name'] for c in selected if c['name'] not in original],
                  admission_notes=excluded,
                  eligible_not_scheduled=sorted(set(admitted)-{c['name'] for c in selected}),
                  reserved=dict(start_pdt='2026-09-24 00:08', end_pdt='2026-09-24 00:28', minutes=20,
                                purpose='Deeper exposures or delays; final standard starts at 00:28'),
                  solver=dict(message=result.message, relative_gap=float(result.mip_gap)))
    (DEST/'snr5_plan.json').write_text(json.dumps(output,indent=2))
    pd.DataFrame(selected).to_csv(DEST/'snr5_scheduling_choices.csv', index=False)
    print(pd.DataFrame(selected).to_string(index=False),flush=True)
    print('Promoted:',output['promoted'],'Solver:',output['solver'],flush=True)


if __name__=='__main__':
    main()

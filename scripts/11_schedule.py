"""
11_schedule.py  --  time-aware allocation: turn the priority ranking into an observing sequence for each night.

04_score_tiers.py fills each night's time budget by priority per hour but does not know *when* a target is up.  For a
half night in September the RA 15.5-17.5h targets all set within the first two hours, and on the October nights the
RA 7-11h targets only rise after midnight, so a budget-based list leaves half the night idle and half the list unobservable.

This script walks each night in 5-minute steps (Palomar, astronomical twilight to the window end) and at every free moment
starts the best candidate observable now:  score = priority_night / (t_exp + overhead)  x  2 if it is one of 04's picks
(they carry the proposal's tier floors)  x  a deadline factor up to 4 when the target sets within 2.5 h  x  0.6 if it will be
much better placed later  x  an airmass factor.  Tier caps (QUOTA in 04) are respected.  Candidates are every scored object
observable that night with moon separation >= 30 deg (12,956 on Oct 26), so the idle hours get filled with the best of the pool.
A CALSPEC spectrophotometric standard opens and closes the night (the Quicklook DRP only builds sensitivity functions from
HST CALSPEC stars: data/standards_spectrophotometric.csv).

Outputs
  data/schedule_<night>.csv         one row per block (UT/local start, end, kind, name, airmass, moon sep, plan, tier, rank)
  finders/<night>/schedule_<night>.txt   printable sequence for the control room
  data/targets_<night>.csv          updated: rank = order in the sequence (1..N), unscheduled 04 picks -> rank 0 (backup),
                                    fillers from the pool appended (origin = 'schedule-fill'), sched_start = local start time
Usage: /opt/anaconda3/bin/python 11_schedule.py [sep23 oct26 oct27]        (run after 04, before 06/06b/07 -- see rescore.sh)
"""
import os, sys, importlib
import numpy as np
import pandas as pd
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, get_body

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); DATA = os.path.join(HERE, 'data')
sys.path.insert(0, HERE)
obs = importlib.import_module('05_observability'); sc = importlib.import_module('04_score_tiers')
PALOMAR, NIGHTS, AIRMASS_MAX = obs.PALOMAR, obs.NIGHTS, obs.AIRMASS_MAX
STEP_MIN = 5
OVERHEAD_MIN = getattr(sc, 'OVERHEAD_MIN', 5)      # slew + acquisition (< 2 min) + readout
QUOTA = getattr(sc, 'QUOTA', {'T3': 3, 'T4': 4})
# ---- selection v2 (14_select_targets.py): CLAGN_SEL=v2 switches the candidate table, the target files, the tier vocabulary and the caps
V2 = os.environ.get('CLAGN_SEL', 'v1') == 'v2'
SUFFIX = '_v2' if V2 else ''
MASTER = 'candidates_v2.csv' if V2 else 'master_list_scored.csv'
TIERS = ['D1', 'D2', 'K', 'B', 'C'] if V2 else ['T1', 'T2', 'T3', 'T4']
EXEMPT_Z = 'K' if V2 else 'T3'
if V2:   # per-night caps from the v2 time fractions (D2 20 %, K 10 %, B 8 %, C 7 %) at ~13 min per target
    QUOTA = {'D2': {'sep23': 5, 'oct26': 11, 'oct27': 11}, 'K': {'sep23': 2, 'oct26': 5, 'oct27': 5}, 'B': {'sep23': 2, 'oct26': 4, 'oct27': 4}}
STD_MIN, STD_AIRMASS, STD_VMAX = 10, 1.6, 13.6     # standard: acquisition + 2 short exposures; bright enough for the full moon
PICK04_BOOST = 2.0
BASE_OVH_MIN, SLEW_MIN0, SLEW_DEG_PER_MIN = 3.5, 0.5, 20.0   # v2: acquisition+readout, plus a slew of 0.5 min + 1 min per 20 deg (dome + settle)
DEADLINE_MIN = 150.0                               # setting within this many minutes -> urgency grows to 4x
COLS04 = ['rank', 'night', 'tier', 'name', 'ra', 'dec', 'z', 'r_mag', 't_exp_min', 'exp_plan', 'prio_per_hour', 'priority_night', 'priority', 'M', 'P', 'trend',
          'years_since_last_spec', 'n_spec', 'last_class', 'clagn_score', 'zeltyn_density_ratio', 'in_region_zeltyn',
          'M_combined', 'fracflux_w1', 'n_ps1_8as', 'blend_flag', 'blend_kind', 'lines_in_ngps', 'spec_dir', 'reversal_candidate', 'notes', 'why', 'why_night']


def local(t):
    return pd.Timestamp(t.utc.datetime, tz='UTC').tz_convert('US/Pacific').strftime('%H:%M')


def candidates(night, m, picked):
    """Everything observable that night, scored like 04 (priority x moon weight, exposure model), minus earlier nights' picks."""
    hrs, sep = m[f'hrs_{night}'], m[f'moonsep_{night}']
    elig = (hrs > 0) & (sep >= sc.MIN_MOONSEP) & ((m.z <= sc.ZMAX) | (m.tier == EXEMPT_Z)) & m.tier.isin(TIERS) \
           & (m.priority > 0) & ~m.name.isin(picked)
    c = m[elig].copy()
    c['priority_night'] = (c.priority * sc.moon_weight(sep[elig].values)).round(3)
    c['t_exp_min'] = sc.exposure_minutes(c.r_mag.values, c.z.values, sep[elig].values).round(0)
    nexp = np.maximum(2, np.ceil(c.t_exp_min / 10.0)).astype(int)
    c['exp_plan'] = [f'{n} x {t/n:.0f} min' for n, t in zip(nexp, c.t_exp_min)]
    c['prio_per_hour'] = (60.0 * c.priority_night / (c.t_exp_min + OVERHEAD_MIN)).round(3)
    return c


def schedule_night(night, m, picked):
    date_str, part = NIGHTS[night]
    t0, t1, dusk, dawn = obs.night_window(date_str, part)
    n = int((t1 - t0).to(u.min).value // STEP_MIN) + 1
    times = t0 + np.arange(n) * STEP_MIN * u.min
    moon = get_body('moon', times, PALOMAR.location)

    t04 = pd.read_csv(os.path.join(DATA, f'targets_{night}{SUFFIX}.csv'))
    if 'origin' in t04:
        t04 = t04[~t04.origin.isin(['schedule-fill', 'pick-other-night'])]   # re-runnable: drop last run's appended rows
    if 'rank04' not in t04:                                          # first run: remember 04's own ranking so reruns see the same picks
        t04['rank04'] = t04['rank']
    t04['rank04'] = t04.rank04.fillna(0).astype(int)
    pick04 = set(t04[t04.rank04 > 0].name) - picked          # a pick already scheduled on an earlier night must not be placed again
    # picks of the other nights that are still unscheduled compete here with the same boost (Oct 26 leftovers -> Oct 27 morning)
    for other in NIGHTS:
        po = os.path.join(DATA, f'targets_{other}{SUFFIX}.csv')
        if other != night and os.path.exists(po):
            o = pd.read_csv(po, usecols=lambda col: col in ('name', 'rank', 'rank04', 'origin'))
            rk = o['rank04'] if 'rank04' in o else o['rank']
            pick04 |= set(o.name[(rk.fillna(0) > 0) & (o.get('origin', pd.Series('04', index=o.index)).fillna('04') == '04')]) - picked
    c = candidates(night, m, picked - pick04)
    c['pick04'] = c.name.isin(pick04)
    # 04's picks keep their own (identical formula) exposure numbers; make sure none is lost to the eligibility filter
    lost = set(t04[t04.rank04 > 0].name) - set(c.name)
    if lost:
        print(f'   {night}: {len(lost)} of this night\'s 04 picks not eligible here (moon/z/priority): {sorted(lost)[:5]}')
    c = c.reset_index(drop=True)
    coords = SkyCoord(c.ra.values * u.deg, c.dec.values * u.deg)
    aa = PALOMAR.altaz(times[np.newaxis, :], coords[:, np.newaxis])
    secz = np.where(aa.alt.deg > 5, aa.secz.value, np.inf)
    moonsep_grid = obs.apparent_moon_separation(coords[:, np.newaxis], times).deg
    ok = (secz < AIRMASS_MAX) & (moonsep_grid >= sc.MIN_MOONSEP)
    invalid_prefix = np.c_[np.zeros(len(c), dtype=int), np.cumsum(~ok, axis=1)]
    texp = c.t_exp_min.values; nblk = np.ceil((texp + OVERHEAD_MIN) / STEP_MIN).astype(int)
    last_ok = np.array([np.max(np.where(ok[i])[0]) if ok[i].any() else -1 for i in range(len(c))])
    best_x = np.array([np.min(secz[i]) for i in range(len(c))])
    # 04's picks carry the proposal's tier floors: boost them, the scarce T2/T3/T4 picks more, so they survive the time crunch
    if V2:   # v2: the selection's picks are boosted per stratum so the small blind (B) and control (C) strata survive the time crunch; their caps bound them
        BOOST_V2 = {'D1': PICK04_BOOST, 'D2': PICK04_BOOST, 'K': 2 * PICK04_BOOST, 'B': 5 * PICK04_BOOST, 'C': 5 * PICK04_BOOST}
        boost = np.where(c.pick04.values, np.array([BOOST_V2.get(str(t), PICK04_BOOST) for t in c.tier.values]), 1.0)
    else:
        boost = np.where(c.pick04.values, np.where(c.tier.values.astype(str) == 'T1', PICK04_BOOST, 2 * PICK04_BOOST), 1.0)
    base = c.priority_night.values / (texp + OVERHEAD_MIN) * boost

    std = pd.read_csv(os.path.join(DATA, 'standards_spectrophotometric.csv'))
    scd = SkyCoord(std.ra.values * u.deg, std.dec.values * u.deg)
    ssecz = np.where(PALOMAR.altaz(times[np.newaxis, :], scd[:, np.newaxis]).alt.deg > 5, PALOMAR.altaz(times[np.newaxis, :], scd[:, np.newaxis]).secz.value, np.inf)

    rows = []; done = np.zeros(len(c), bool); counts = {k: 0 for k in QUOTA}
    def block(kind, name, ra, dec, k, nb, minutes, X, msep, plan, label='', tier='', prio=np.nan, z=np.nan, rmag=np.nan, origin='', slew=np.nan):
        rows.append(dict(night=night, kind=kind, name=name, ra=ra, dec=dec, slew_deg=(round(float(slew), 1) if np.isfinite(slew) else np.nan), start_ut=times[k].utc.iso[:16], start_local=local(times[k]),
                         end_local=local(times[min(k + nb, n - 1)]), minutes=int(minutes), airmass=(round(float(X), 2) if np.isfinite(X) else np.nan), moonsep=(round(float(msep)) if np.isfinite(msep) else np.nan),
                         plan=plan, label=label, tier=tier, priority=prio, z=z, r_mag=rmag, origin=origin))
    cur = [None]                                                     # current telescope pointing (SkyCoord), for the slew term
    def add_std(k, label):
        cand = np.where((ssecz[:, k] < STD_AIRMASS) & (std.V.fillna(15).values < STD_VMAX))[0]
        if not len(cand):
            cand = np.where(ssecz[:, k] < 2.0)[0]
        j = cand[np.argmin(ssecz[cand, k] + 0.05 * std.V.fillna(15).values[cand])]
        nb = int(np.ceil(STD_MIN / STEP_MIN))
        block('standard', std.name[j], std.ra[j], std.dec[j], k, nb, STD_MIN, ssecz[j, k], obs.apparent_moon_separation(scd[j],times[k]).deg,
              f'CALSPEC {std.calspec[j]} · {std.sptype[j]} V={std.V[j]:.1f} · 2 x 30-120 s', label, rmag=std.V[j], slew=(cur[0].separation(scd[j]).deg if cur[0] is not None else np.nan))
        cur[0] = scd[j]
        return k + nb

    k = add_std(0, 'start of night')
    k_end = n - 1 - int(np.ceil(STD_MIN / STEP_MIN))                  # keep the last STD_MIN for the closing standard
    while k < k_end:
        if V2:   # slew-aware: the time a candidate costs now depends on where the telescope is
            sep_now = cur[0].separation(coords).deg if cur[0] is not None else np.zeros(len(c))
            tot_now = texp + BASE_OVH_MIN + SLEW_MIN0 + sep_now / SLEW_DEG_PER_MIN
            nblk_now = np.ceil(tot_now / STEP_MIN).astype(int); base_now = c.priority_night.values * boost / tot_now
        else:
            sep_now = np.full(len(c), np.nan); tot_now = texp + OVERHEAD_MIN; nblk_now = nblk; base_now = base
        end_idx = np.minimum(k+nblk_now, n-1)
        complete_block_ok = (invalid_prefix[np.arange(len(c)),end_idx+1]-invalid_prefix[:,k]) == 0
        fits = complete_block_ok & ~done & (k + nblk_now <= k_end)
        for tier, cap in QUOTA.items():
            cap = cap.get(night, 99) if isinstance(cap, dict) else cap
            if counts.get(tier, 0) >= cap:
                fits &= (c.tier.values != tier)
        if not fits.any():
            block('gap', '', np.nan, np.nan, k, 1, STEP_MIN, np.nan, np.nan, ''); k += 1; continue
        ttd = (last_ok - k) * STEP_MIN
        urgency = 1.0 + 3.0 * np.clip(1.0 - ttd / DEADLINE_MIN, 0.0, 1.0)
        patience = np.where((secz[:, k] - best_x > 0.3) & (ttd > DEADLINE_MIN), 0.6, 1.0)
        airm = np.clip(1.3 - 0.3 * secz[:, k], 0.5, 1.0)
        score = np.where(fits, base_now * urgency * patience * airm, -1.0)
        i = int(np.argmax(score)); nb = nblk_now[i]; r = c.iloc[i]
        block('primary' if r.pick04 else 'filler', r['name'], r.ra, r.dec, k, nb, tot_now[i], secz[i, k], moonsep_grid[i,k],
              r.exp_plan, '', r.tier, float(r.priority_night), r.z, r.r_mag, '04' if r.pick04 else 'schedule-fill', slew=sep_now[i])
        done[i] = True; counts[r.tier] = counts.get(r.tier, 0) + 1; k += nb; cur[0] = coords[i]
    add_std(k_end, 'end of night')

    S = pd.DataFrame(rows)
    sched = S[S.kind.isin(['primary', 'filler'])].reset_index(drop=True)
    idle = int(S[S.kind == 'gap'].minutes.sum())
    unplaced = sorted(set(t04[t04.rank04 > 0].name) - set(sched.name))
    print(f'{night}: {local(t0)}-{local(t1)} local ({(t1 - t0).to(u.hour).value:.1f} h) | {int(sched.name.isin(t04[t04.rank04 > 0].name).sum())} of {int((t04.rank04 > 0).sum())} own 04-picks placed, {int(sched.origin.eq("04").sum() - sched.name.isin(t04[t04.rank04 > 0].name).sum())} picks from other nights, '
          f'{int(sched.origin.eq("schedule-fill").sum())} fillers from the pool, {idle} min idle | standards {", ".join(S[S.kind == "standard"].name)} | '
          f'tiers {sched.tier.value_counts().to_dict()}' + (f' | slews total {S.slew_deg.sum():.0f} deg, > 30 deg: {int((S.slew_deg > 30).sum())}' if V2 else ''))
    if unplaced:
        print(f'   unplaced 04-picks (now backups): {", ".join(unplaced)}')
    S.to_csv(os.path.join(DATA, f'schedule_{night}{SUFFIX}.csv'), index=False)

    # ---- update targets_<night>.csv: rank = sequence order; unplaced 04 picks -> backups; fillers appended with 04's columns
    order = {nm: i + 1 for i, nm in enumerate(sched.name)}
    start = dict(zip(sched.name, sched.start_local))
    t04 = t04.copy(); t04['origin'] = '04'
    t04['rank'] = t04.name.map(order).fillna(0).astype(int)
    t04['sched_start'] = t04.name.map(start).fillna('')
    # rows to append: pool fillers, and picks of other nights scheduled here (both absent from this night's 04 file)
    extra = sched[~sched.name.isin(t04.name)]
    if len(extra):
        f = c[c.name.isin(extra.name)].copy(); f['night'] = night
        f['origin'] = f.name.map(dict(zip(extra.name, extra.origin))).replace({'04': 'pick-other-night'})
        f['rank'] = f.name.map(order).astype(int); f['sched_start'] = f.name.map(start); f['rank04'] = 0
        keep = [x for x in (list(t04.columns) if V2 else COLS04) + [f'hrs_{night}', f'minX_{night}', f'moonsep_{night}', 'origin', 'sched_start'] if x in f.columns]
        keep = list(dict.fromkeys(keep))
        t04 = pd.concat([t04, f[keep]], ignore_index=True)
    t04 = t04.sort_values(['rank', 'prio_per_hour'], ascending=[True, False], key=lambda s: s.replace(0, 10**6) if s.name == 'rank' else s)
    t04.to_csv(os.path.join(DATA, f'targets_{night}{SUFFIX}.csv'), index=False)

    d = os.path.join(HERE, 'finders', night + SUFFIX); os.makedirs(d, exist_ok=True)
    with open(os.path.join(d, f'schedule_{night}.txt'), 'w') as fh:
        fh.write(f'# NGPS observing sequence {night}: {date_str} ({part} night), local Pacific times; astronomical twilight {local(dusk)} -> {local(dawn)}\n')
        fh.write(f'# window {local(t0)} -> {local(t1)}; airmass < {AIRMASS_MAX}; order = priority per hour with setting targets first; '
                 f'{idle} min unfilled\n')
        fh.write(f'{"start":>5} {"end":>5} {"kind":8} {"name":20} {"RA":>10} {"Dec":>9} {"X":>4} {"moon":>4} {"tier":4} {"#":>3} {"min":>3}  plan\n')
        seq = 0
        for r in S.itertuples():
            if r.kind == 'gap':
                continue
            seq += 1 if r.kind != 'standard' else 0
            cc = SkyCoord(r.ra * u.deg, r.dec * u.deg)
            fh.write(f'{r.start_local:>5} {r.end_local:>5} {r.kind:8} {str(r.name):20} {cc.ra.to_string(u.hour, sep=":", precision=0, pad=True):>10} '
                     f'{cc.dec.to_string(u.deg, sep=":", precision=0, alwayssign=True, pad=True):>9} {r.airmass:4.2f} {r.moonsep:4.0f} {str(r.tier):4} '
                     f'{(seq if r.kind != "standard" else 0):3d} {r.minutes:3d}  {r.plan} {r.label}\n')
    return set(sched.name)


if __name__ == '__main__':
    m = pd.read_csv(os.path.join(DATA, MASTER), low_memory=False)
    picked = set()
    for night in (sys.argv[1:] or list(NIGHTS)):
        picked |= schedule_night(night, m, picked)

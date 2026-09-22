"""Plan the September 23 visits with one instrument setting and a geometry-optimised order.

Every science target uses the same slit, binning and 300 s sub-exposures. The S/N screen runs
the official NGPS ETC with seeing scaled to the target airmass and a moonlit-sky model
evaluated at each candidate start time, so "lowest airmass" and "distance to the Moon" enter
through one physical quantity: the predicted continuum S/N per Angstrom near observed
H-beta. Primaries come from observing/sep23/user_selection.json when it exists (the user's
choice from the decision board, with optional per-target exposure counts); otherwise the
primaries of the archived packet are protected. Remaining time is filled from the admitted
pool by the existing review score. Writes a private plan; does not render pages or change
telescope CSVs.
"""
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import hashlib
import math
from types import SimpleNamespace
import importlib
import json
import re

import astropy.units as u
import numpy as np
import pandas as pd
from astropy.coordinates import AltAz, SkyCoord, get_body
from astropy.time import Time
from astropy.utils import iers
from astroplan import moon_illumination
from scipy.optimize import Bounds, LinearConstraint, milp
from scipy.sparse import coo_matrix

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'data/reselection_2026-09-20'
DEST = ROOT / 'observing/sep23'
SELECTION = DEST / 'user_selection.json'
PROTECTED_SEED = ROOT / 'archive/2026-09-20/before_slit15_packet/sep23/packet.json'
ORIGINAL_SEED = ROOT / 'archive/2026-09-20/before_snr5_packet/sep23/packet.json'
CACHE = OUT / 'sep23_slit15_etc'
TZ = 'America/Los_Angeles'
iers.conf.auto_download = False
iers.conf.auto_max_age = None
OBS = importlib.import_module('05_observability')

SETTINGS = dict(
    model_version='zenith-seeing-v2', overhead_convention='Six minutes includes normal two-exposure readouts; add 0.6 minute per extra exposure and round up.',
    mode='fixed_exposure', goal='combined continuum S/N >= 5 per Angstrom near observed H-beta',
    goal_snr=5, goal_unit='per Angstrom', exposures=2, seconds_each=300, readout_minutes=0.6,
    slit_arcsec=1.5, binspat=2, binspect=3, slitangle='PA',
    seeing_zenith_500nm=1.3, seeing_airmass_power=0.6,
    sky_model='Krisciunas & Schaefer 1991 moonlight model; k_V=0.17; dark zenith V=21.5; evaluated at mid-visit',
    moon_min_deg=40, airmass_preferred=1.5, airmass_max=1.8,
    overhead_minutes=6, visit_minutes=16, grid_minutes=4, reserve_minutes=8, end_reserve_minutes=8,
    science_start_pdt='2026-09-23 20:16', science_end_pdt='2026-09-24 00:28',
    geometry_weight=15.0,
    normalization='Archival local H-beta continuum; no brightness forecast',
    extraction='Single slit; optimal point-source extraction; official NGPS ETC with May 2026 read noise and plate scales')
TIERS = [('preferred', SETTINGS['airmass_preferred']), ('extended', SETTINGS['airmass_max'])]


def visit_minutes(nexp):
    """Conservative whole-minute visits under the documented inclusive overhead allowance.

    The normal two-exposure visit already includes readout in its six minutes.
    Each additional exposure adds integration plus a readout. This is a planning
    allowance, not a measured guarantee of acquisition/slew time.
    """
    if nexp < 1:
        raise ValueError('Exposure count must be positive')
    return math.ceil(nexp * SETTINGS['seconds_each'] / 60 + max(0, nexp - 2) * SETTINGS['readout_minutes'] + SETTINGS['overhead_minutes'])


def sky_v(alpha_deg, rho_deg, moon_alt_deg, target_alt_deg, k=0.17, vdark=21.5):
    """Johnson V sky brightness (mag/arcsec^2) with moonlight, Krisciunas & Schaefer (1991)."""
    def airmass(zenith_deg):
        return (1 - 0.96 * np.sin(np.radians(zenith_deg))**2)**-0.5
    def nanolamberts(v):
        return 34.08 * np.exp(20.7233 - 0.92104 * v)
    def magnitude(b):
        return (20.7233 - np.log(b / 34.08)) / 0.92104
    z_moon = 90 - moon_alt_deg
    z = 90 - target_alt_deg
    dark = nanolamberts(vdark) * airmass(z) * 10**(-0.4 * k * (airmass(z) - 1))
    if moon_alt_deg <= 0:
        return float(magnitude(dark))
    m = -12.73 + 0.026 * abs(alpha_deg) + 4e-9 * alpha_deg**4
    istar = 10**(-0.4 * (m + 16.57))
    f = 10**5.36 * (1.06 + np.cos(np.radians(rho_deg))**2) + 10**(6.15 - rho_deg / 40.)
    moon = f * istar * 10**(-0.4 * k * airmass(z_moon)) * (1 - 10**(-0.4 * k * airmass(z)))
    return float(magnitude(moon + dark))


def stamp(offset_minutes):
    return pd.Timestamp(SETTINGS['science_start_pdt'], tz=TZ) + pd.Timedelta(minutes=int(offset_minutes))


def label(t):
    return t.strftime('%Y-%m-%d %H:%M')


def night_geometry(targets):
    """Per-minute altitude, airmass and Moon geometry for every target across the science block."""
    start = pd.Timestamp(SETTINGS['science_start_pdt'], tz=TZ)
    end = pd.Timestamp(SETTINGS['science_end_pdt'], tz=TZ)
    minutes = int((end - start).total_seconds() / 60)
    index = pd.DatetimeIndex([start + pd.Timedelta(minutes=m) for m in range(minutes + 1)])
    times = Time(index.tz_convert('UTC').tz_localize(None).to_pydatetime())
    frame = AltAz(obstime=times, location=OBS.PALOMAR.location, pressure=0 * u.hPa)
    coords = SkyCoord(targets.ra.to_numpy() * u.deg, targets.dec.to_numpy() * u.deg)
    aa = coords[:, None].transform_to(frame)
    moon = get_body('moon', times, OBS.PALOMAR.location).transform_to(frame)
    illumination = np.asarray(moon_illumination(times), float)
    return dict(minutes=minutes, alt=aa.alt.deg, airmass=aa.secz.value, moon_sep=aa.separation(moon).deg,
                moon_alt=moon.alt.deg, phase_angle=np.degrees(np.arccos(np.clip(2 * illumination - 1, -1, 1))),
                illumination=illumination)


def feasible_runs(mask, visit):
    """Start minutes at which a full visit stays inside the mask, grouped into runs."""
    ok = np.array([mask[m:m + visit + 1].all() for m in range(len(mask) - visit)])
    runs = []
    edges = np.diff(np.r_[False, ok, False].astype(int))
    for a, b in zip(np.where(edges == 1)[0], np.where(edges == -1)[0]):
        runs.append(dict(start=label(stamp(a)), end=label(stamp(b - 1 + visit)), latest_visit_start=label(stamp(b - 1))))
    return runs


def candidate_slots(geo, i, visit):
    """Feasible visits of the given length for target i on the planning grid, with mid-visit sky and seeing."""
    alt = geo['alt'][i]
    x = geo['airmass'][i]
    sep = geo['moon_sep'][i]
    good = (alt > 0) & (x <= SETTINGS['airmass_max']) & (sep >= SETTINGS['moon_min_deg'])
    windows = {tier: feasible_runs(good & (x <= ceiling), visit) for tier, ceiling in TIERS}
    slots = []
    for offset in range(0, geo['minutes'] - visit - SETTINGS['end_reserve_minutes'] + 1, SETTINGS['grid_minutes']):
        span = slice(offset, offset + visit + 1)
        if not good[span].all():
            continue
        tier = 'preferred' if x[span].max() <= SETTINGS['airmass_preferred'] else 'extended'
        mid = offset + visit // 2
        x_mean = float(x[span].mean())
        slots.append(dict(offset=offset, duration=visit, start_pdt=label(stamp(offset)), end_pdt=label(stamp(offset + visit)), tier=tier,
                          airmass_start=float(x[offset]), airmass_end=float(x[offset + visit]),
                          airmass_mean=x_mean, airmass_max=float(x[span].max()),
                          moon_deg=float(sep[span].min()), moon_alt_deg=float(geo['moon_alt'][mid]),
                          sky_V=sky_v(geo['phase_angle'][mid], float(sep[mid]), float(geo['moon_alt'][mid]), float(alt[mid])),
                          seeing_arcsec=float(SETTINGS['seeing_zenith_500nm'] * x_mean**SETTINGS['seeing_airmass_power'])))
    return slots, windows


def evaluate(job):
    """ETC S/N for every candidate slot of one target; cached by settings, reference, exposures and geometry."""
    name, z, ref, nexp, slots = job
    signature = hashlib.sha256(json.dumps([SETTINGS, z, ref, nexp, slots], sort_keys=True).encode()).hexdigest()
    path = CACHE / f'{name}.json'
    if path.exists():
        cached = json.loads(path.read_text())
        if cached['signature'] == signature:
            return name, cached['plans'], cached['slots']
    model = importlib.import_module('22_september_etc')
    for ch, rn, scale in zip(model.CFG.channels, [2.8, 7.8, 3.7, 4.6], [.193, .189, .186, .186]):
        model.CFG.readnoise[ch] = rn * u.count / u.pix
        model.CFG.platescale[ch] = scale * u.arcsec / u.pix
    mjd, mag, private, flux = ref
    wave = 486.27 * (1 + z)
    channels = [ch for ch in ['R', 'I', 'G', 'U']
                if model.CFG.channelRange[ch][0].to_value(u.nm) <= wave - 4
                and model.CFG.channelRange[ch][1].to_value(u.nm) >= wave + 4]
    if not channels:
        path.write_text(json.dumps(dict(signature=signature, plans={}, slots={}), indent=2))
        return name, {}, {}
    channel = channels[0]
    bin_angstrom = float((model.CFG.dLambda[channel] * SETTINGS['binspect']).to_value(u.AA))
    table = {}
    for s in slots:
        cmd = [channel, str(wave - 4), str(wave + 4), 'EXPTIME', str(SETTINGS['seconds_each']),
               '-slit', 'SET', str(SETTINGS['slit_arcsec']), '-binspect', str(SETTINGS['binspect']),
               '-binspat', str(SETTINGS['binspat']), '-seeing', str(SETTINGS['seeing_zenith_500nm']), '500',
               '-airmass', f"{s['airmass_mean']:.3f}", '-skymag', f"{s['sky_V']:.2f}", '-mag', f'{mag:.4f}',
               '-magsystem', 'AB', '-magfilter', 'match', '-noslicer']
        args = model.ETC.parser.parse_args(cmd)
        model.ETC.check_inputs_add_units(args)
        per_bin = float(model.ETC.main(args, quiet=True)['SNR'].value) * np.sqrt(nexp)
        per_angstrom = per_bin / np.sqrt(bin_angstrom)
        table[s['start_pdt']] = s | dict(exposures=nexp, snr_per_bin=per_bin, snr_per_angstrom=per_angstrom,
                                         meets_goal=bool(per_angstrom >= SETTINGS['goal_snr']))
    plans = {}
    for tier, ceiling in TIERS:
        tier_slots = [v for v in table.values() if v['tier'] == tier]
        if not tier_slots:
            continue
        good = [v['snr_per_angstrom'] for v in tier_slots if v['meets_goal']]
        every = [v['snr_per_angstrom'] for v in tier_slots]
        plans[tier] = dict(exposures=nexp, seconds_each=SETTINGS['seconds_each'],
                           integration_minutes=nexp * SETTINGS['seconds_each'] / 60,
                           visit_minutes=visit_minutes(nexp), overhead_minutes=SETTINGS['overhead_minutes'],
                           airmass=ceiling, goal_snr=SETTINGS['goal_snr'], goal_unit=SETTINGS['goal_unit'],
                           predicted_snr=min(good) if good else min(every), predicted_snr_best=max(every),
                           slots_meeting_goal=len(good), slots=len(tier_slots), below_floor=not good,
                           reference_mjd=mjd, reference_private=bool(private), continuum_AB=mag,
                           channel=channel, low_nm=wave - 4, high_nm=wave + 4, bin_angstrom=bin_angstrom,
                           seeing_zenith=SETTINGS['seeing_zenith_500nm'], status='scenario')
    path.write_text(json.dumps(dict(signature=signature, plans=plans, slots=table), indent=2))
    return name, plans, table


def main():
    CACHE.mkdir(parents=True, exist_ok=True)
    previous = [v['name'] for v in json.loads(PROTECTED_SEED.read_text())['primaries']]
    original = [v['name'] for v in json.loads(ORIGINAL_SEED.read_text())['primaries']]
    selection = json.loads(SELECTION.read_text()) if SELECTION.exists() else None
    forced = [p['name'] for p in selection['primaries']] if selection else previous
    locked = (selection or {}).get('preserve_sequence', [])
    retained = set(forced) | {v['name'] for v in locked}
    nexp_override = {k: int(v) for k, v in (selection or {}).get('nexp', {}).items()}
    targets = pd.read_csv(OUT / 'compact_review_objects.csv').set_index('name', drop=False)
    science = pd.read_csv(OUT / 'three_night_review/science_and_sensitivity.csv').set_index('name')
    opts = json.loads((OUT / 'airmass_options.json').read_text())
    public = json.loads(re.search(r'<script id="candidate-data" type="application/json">(.*?)</script>',
                                  (ROOT / 'docs/index.html').read_text(), re.S).group(1))
    public_names = {t['name'] for t in public['targets']}
    september = [n for n in targets.index if any(w['night'] == 'sep23' for w in opts['windows'].get(n, []))]
    missing = set(forced) - set(september)
    assert not missing, f'Chosen primaries without a September window: {sorted(missing)}'
    geo = night_geometry(targets.loc[september])
    model = importlib.import_module('22_september_etc')
    jobs = []
    windows = {}
    notes = {}
    exposures = {}
    for i, name in enumerate(september):
        ref = model.reference(targets.loc[name].to_dict())
        if ref is None:
            notes[name] = ['No accepted continuum reference']
            continue
        nexp = nexp_override.get(name, SETTINGS['exposures'])
        exposures[name] = nexp
        slots, windows[name] = candidate_slots(geo, i, visit_minutes(nexp))
        if not slots:
            notes[name] = [f"No full {visit_minutes(nexp)}-minute visit at X<={SETTINGS['airmass_max']} and Moon>={SETTINGS['moon_min_deg']} deg"]
            continue
        mjd, mag, private, flux = ref
        jobs.append((name, float(targets.loc[name].z), (float(mjd), float(mag), bool(private), float(flux)), nexp, slots))
    with ProcessPoolExecutor(max_workers=4) as executor:
        results = list(executor.map(evaluate, jobs))
    plans = {name: p for name, p, _ in results}
    slots = {name: s for name, _, s in results}

    admitted = []
    excluded = []
    waived = {}
    for name in september:
        t = targets.loc[name]
        s = science.loc[name]
        why = list(notes.get(name, []))
        if t.pool_role != 'manifold':
            why.append('Comparison reserve; prioritize manifold candidates for promotions')
        if name not in public_names or any(p['reference_private'] for p in plans.get(name, {}).values()):
            why.append('Retained as a local backup; the shared primary packet requires public identity and exposure reference')
        if t.field_status != 'clear':
            why.append('Field requires further review')
        if pd.notna(s.ztf_r_latest180_mag) and s.ztf_r_latest180_mag >= 19:
            why.append('Latest cached r does not meet brightness criterion')
        if name in slots and not any(v['meets_goal'] for v in slots[name].values()):
            why.append(f"No slot reaches continuum S/N {SETTINGS['goal_snr']} per Angstrom at {exposures.get(name, SETTINGS['exposures'])}x{SETTINGS['seconds_each']} s")
        if name in retained:
            assert name in slots, f'Chosen primary {name} has no feasible visit: {why}'
            if why:
                waived[name] = why
            admitted.append(name)
        elif why:
            excluded.append(dict(name=name, reasons=why))
        else:
            admitted.append(name)

    minutes = geo['minutes']
    budget = minutes - SETTINGS['reserve_minutes']
    choices = []
    for name in admitted:
        t = targets.loc[name]
        s = science.loc[name]
        # Scheduling utility: the review score decides which fill targets; slot quality (S/N relative to
        # the target's own best slot) decides where. Chosen primaries are forced in regardless of value.
        value = 20 + float(s.review_order_score) + 5 * bool(t.balmer_pair_in_range) + 5 * float(t.manifold_cl_neighbor_fraction)
        usable = list(slots[name].values()) if name in retained else [v for v in slots[name].values() if v['meets_goal']]
        best = max(v['snr_per_angstrom'] for v in usable)
        # Faint chosen targets (best slot below the floor) get triple weight on slot quality: their S/N is the scarce resource.
        weight = SETTINGS['geometry_weight'] * (3.0 if (name in forced and best < SETTINGS['goal_snr']) else 1.0)
        for v in usable:
            quality = v['snr_per_angstrom'] / best
            choices.append(dict(name=name, tier=v['tier'], offset=v['offset'], duration=v['duration'], start_pdt=v['start_pdt'],
                                end_pdt=v['end_pdt'], exposures=v['exposures'], forced=name in forced,
                                utility=value + weight * quality - 1e-5 * v['offset'],
                                science_value=value, slot_quality=quality,
                                **{k: v[k] for k in ['snr_per_angstrom', 'snr_per_bin', 'meets_goal', 'airmass_mean', 'airmass_max',
                                                     'moon_deg', 'sky_V', 'seeing_arcsec']}))
    rows = []
    cols = []
    vals = []
    name_index = {n: minutes + i for i, n in enumerate(admitted)}
    budget_row = minutes + len(admitted)
    for j, c in enumerate(choices):
        rows.extend(range(c['offset'], c['offset'] + c['duration']))
        cols.extend([j] * c['duration'])
        vals.extend([1.0] * c['duration'])
        rows.extend([name_index[c['name']], budget_row])
        cols.extend([j, j])
        vals.extend([1.0, float(c['duration'])])
    matrix = coo_matrix((np.array(vals), (rows, cols)), shape=(budget_row + 1, len(choices))).tocsc()
    lower = np.zeros(matrix.shape[0])
    upper = np.ones(matrix.shape[0])
    upper[budget_row] = budget
    for name in forced:
        lower[name_index[name]] = 1
    if locked:
        # Fail closed rather than silently replace a retained primary or change its start.
        selected = []
        for v in locked:
            matching = [c for c in choices if c['name'] == v['name'] and c['start_pdt'] == v['start_pdt']]
            assert len(matching) == 1, f"Retained primary needs a revised schedule: {v}"
            selected.append(matching[0])
        selected.sort(key=lambda c: c['offset'])
        result = SimpleNamespace(message='Validated preserved sequence; no membership or start-time changes', mip_gap=0.)
    else:
        result = milp(-np.array([c['utility'] for c in choices]), integrality=np.ones(len(choices)),
                      bounds=Bounds(0, 1), constraints=LinearConstraint(matrix, lower, upper),
                      options=dict(time_limit=90, mip_rel_gap=.001))
        assert result.x is not None, f'No feasible sequence for the chosen primaries: {result.message}'
        selected = sorted((choices[i] for i in np.flatnonzero(result.x > .5)), key=lambda c: c['offset'])
    names = [c['name'] for c in selected]
    assert set(forced) <= set(names) and len(set(names)) == len(names)
    assert sum(c['duration'] for c in selected) <= budget
    for a, b in zip(selected, selected[1:]):
        assert a['offset'] + a['duration'] <= b['offset']
    gaps = []
    cursor = 0
    for c in selected + [dict(offset=minutes, duration=0)]:
        if c['offset'] > cursor:
            gaps.append(dict(start_pdt=label(stamp(cursor)), end_pdt=label(stamp(c['offset'])), minutes=c['offset'] - cursor))
        cursor = c['offset'] + c['duration']
    for rank, c in enumerate(selected, 1):
        c['rank'] = rank

    output = dict(settings=SETTINGS, plans=plans, slots=slots, windows=windows, sequence=selected,
                  protected=sorted(forced), original_primaries=sorted(original), previous_packet_primaries=sorted(previous),
                  user_selection=selection, exposures={n: exposures.get(n, SETTINGS['exposures']) for n in names},
                  promoted=[n for n in names if n not in forced], demoted=[n for n in previous if n not in names],
                  waived_rules={n: waived[n] for n in names if n in waived},
                  pool=dict(september_eligible=len(september), with_slots=len(slots), admitted=len(admitted),
                            choices=len(choices), budget_minutes=int(budget), used_minutes=int(sum(c['duration'] for c in selected))),
                  admission_notes=excluded,
                  eligible_not_scheduled=sorted(set(admitted) - set(names)),
                  reserved=dict(gaps=gaps, minutes=int(sum(g['minutes'] for g in gaps)), end_reserve_minutes=SETTINGS['end_reserve_minutes'],
                                purpose='Unscheduled time inside the science block: a 16-minute gap can take one backup visit if on schedule; the final gap before the closing standard absorbs delays; not a target row'),
                  solver=dict(message=result.message, relative_gap=float(result.mip_gap)))
    (DEST / 'snr5_plan.json').write_text(json.dumps(output, indent=2))
    frame = pd.DataFrame(selected)
    frame.to_csv(DEST / 'snr5_scheduling_choices.csv', index=False)
    table = pd.DataFrame([dict(name=n, admitted=n in admitted, **v) for n, tab in slots.items() for v in tab.values()])
    table.to_csv(DEST / 'snr5_slot_table.csv', index=False)
    show = ['rank', 'name', 'start_pdt', 'duration', 'exposures', 'tier', 'airmass_mean', 'moon_deg', 'sky_V', 'snr_per_angstrom', 'meets_goal', 'forced', 'slot_quality']
    print(frame[show].round(2).to_string(index=False), flush=True)
    print('Forced:', len(forced), 'Fill:', output['promoted'], 'Demoted from the previous packet:', output['demoted'], flush=True)
    print('Waived rules:', json.dumps(output['waived_rules'], indent=1), flush=True)
    print('Reserve minutes:', output['reserved']['minutes'], gaps, flush=True)
    print('Pool:', output['pool'], 'Solver:', output['solver'], flush=True)


if __name__ == '__main__':
    main()

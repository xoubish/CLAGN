"""Plan the September 23 visits with one instrument setting and a geometry-optimised order.

Every science target uses the same slit, binning and 2x300 s exposures. The S/N screen runs
the official NGPS ETC with seeing scaled to the target airmass and a moonlit-sky model
evaluated at each candidate start time, so "lowest airmass" and "distance to the Moon" enter
through one physical quantity: the predicted continuum S/N per Angstrom near observed
H-beta. The eleven primaries of the previous packet are protected; further targets from the
reviewed September pool fill the remaining slots. Writes a private plan; does not render
pages or change telescope CSVs.
"""
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import hashlib
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
PROTECTED_SEED = ROOT / 'archive/2026-09-20/before_slit15_packet/sep23/packet.json'
ORIGINAL_SEED = ROOT / 'archive/2026-09-20/before_snr5_packet/sep23/packet.json'
CACHE = OUT / 'sep23_slit15_etc'
TZ = 'America/Los_Angeles'
iers.conf.auto_download = False
iers.conf.auto_max_age = None
OBS = importlib.import_module('05_observability')

SETTINGS = dict(
    mode='fixed_exposure', goal='combined continuum S/N >= 5 per Angstrom near observed H-beta',
    goal_snr=5, goal_unit='per Angstrom', exposures=2, seconds_each=300,
    slit_arcsec=1.5, binspat=2, binspect=3, slitangle='PA',
    seeing_zenith_500nm=1.3, seeing_airmass_power=0.6,
    sky_model='Krisciunas & Schaefer 1991 moonlight model; k_V=0.17; dark zenith V=21.5; evaluated at mid-visit',
    moon_min_deg=40, airmass_preferred=1.5, airmass_max=1.8,
    overhead_minutes=6, visit_minutes=16, grid_minutes=4, reserve_minutes=24, end_reserve_minutes=12,
    science_start_pdt='2026-09-23 20:16', science_end_pdt='2026-09-24 00:28',
    geometry_weight=15.0,
    normalization='Archival local H-beta continuum; no brightness forecast',
    extraction='Single slit; optimal point-source extraction; official NGPS ETC with May 2026 read noise and plate scales')
TIERS = [('preferred', SETTINGS['airmass_preferred']), ('extended', SETTINGS['airmass_max'])]


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


def candidate_slots(geo, i):
    """Feasible 16-minute visits for target i on the planning grid, with mid-visit sky and seeing."""
    visit = SETTINGS['visit_minutes']
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
        slots.append(dict(offset=offset, start_pdt=label(stamp(offset)), end_pdt=label(stamp(offset + visit)), tier=tier,
                          airmass_start=float(x[offset]), airmass_end=float(x[offset + visit]),
                          airmass_mean=x_mean, airmass_max=float(x[span].max()),
                          moon_deg=float(sep[span].min()), moon_alt_deg=float(geo['moon_alt'][mid]),
                          sky_V=sky_v(geo['phase_angle'][mid], float(sep[mid]), float(geo['moon_alt'][mid]), float(alt[mid])),
                          seeing_arcsec=float(SETTINGS['seeing_zenith_500nm'] * x_mean**SETTINGS['seeing_airmass_power'])))
    return slots, windows


def evaluate(job):
    """ETC S/N for every candidate slot of one target; cached by settings, reference and geometry."""
    name, z, ref, slots = job
    signature = hashlib.sha256(json.dumps([SETTINGS, z, ref, slots], sort_keys=True).encode()).hexdigest()
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
               '-binspat', str(SETTINGS['binspat']), '-seeing', f"{s['seeing_arcsec']:.3f}", '500',
               '-airmass', f"{s['airmass_mean']:.3f}", '-skymag', f"{s['sky_V']:.2f}", '-mag', f'{mag:.4f}',
               '-magsystem', 'AB', '-magfilter', 'match', '-noslicer']
        args = model.ETC.parser.parse_args(cmd)
        model.ETC.check_inputs_add_units(args)
        per_bin = float(model.ETC.main(args, quiet=True)['SNR'].value) * np.sqrt(SETTINGS['exposures'])
        per_angstrom = per_bin / np.sqrt(bin_angstrom)
        table[s['start_pdt']] = s | dict(snr_per_bin=per_bin, snr_per_angstrom=per_angstrom,
                                         meets_goal=bool(per_angstrom >= SETTINGS['goal_snr']))
    plans = {}
    for tier, ceiling in TIERS:
        good = [v['snr_per_angstrom'] for v in table.values() if v['tier'] == tier and v['meets_goal']]
        if not good:
            continue
        plans[tier] = dict(exposures=SETTINGS['exposures'], seconds_each=SETTINGS['seconds_each'],
                           integration_minutes=SETTINGS['exposures'] * SETTINGS['seconds_each'] / 60,
                           visit_minutes=SETTINGS['visit_minutes'], overhead_minutes=SETTINGS['overhead_minutes'],
                           airmass=ceiling, goal_snr=SETTINGS['goal_snr'], goal_unit=SETTINGS['goal_unit'],
                           predicted_snr=min(good), predicted_snr_best=max(good), slots_meeting_goal=len(good),
                           reference_mjd=mjd, reference_private=bool(private), continuum_AB=mag,
                           channel=channel, low_nm=wave - 4, high_nm=wave + 4, bin_angstrom=bin_angstrom,
                           seeing_zenith=SETTINGS['seeing_zenith_500nm'], status='scenario')
    path.write_text(json.dumps(dict(signature=signature, plans=plans, slots=table), indent=2))
    return name, plans, table


def main():
    CACHE.mkdir(parents=True, exist_ok=True)
    protected = [v['name'] for v in json.loads(PROTECTED_SEED.read_text())['primaries']]
    original = [v['name'] for v in json.loads(ORIGINAL_SEED.read_text())['primaries']]
    targets = pd.read_csv(OUT / 'compact_review_objects.csv').set_index('name', drop=False)
    science = pd.read_csv(OUT / 'three_night_review/science_and_sensitivity.csv').set_index('name')
    opts = json.loads((OUT / 'airmass_options.json').read_text())
    public = json.loads(re.search(r'<script id="candidate-data" type="application/json">(.*?)</script>',
                                  (ROOT / 'docs/index.html').read_text(), re.S).group(1))
    public_names = {t['name'] for t in public['targets']}
    september = [n for n in targets.index if any(w['night'] == 'sep23' for w in opts['windows'].get(n, []))]
    assert set(protected) <= set(september), 'A protected primary is missing from the September pool'
    geo = night_geometry(targets.loc[september])
    model = importlib.import_module('22_september_etc')
    jobs = []
    windows = {}
    notes = {}
    for i, name in enumerate(september):
        ref = model.reference(targets.loc[name].to_dict())
        if ref is None:
            notes[name] = ['No accepted continuum reference']
            continue
        slots, windows[name] = candidate_slots(geo, i)
        if not slots:
            notes[name] = [f"No full {SETTINGS['visit_minutes']}-minute visit at X<={SETTINGS['airmass_max']} and Moon>={SETTINGS['moon_min_deg']} deg"]
            continue
        mjd, mag, private, flux = ref
        jobs.append((name, float(targets.loc[name].z), (float(mjd), float(mag), bool(private), float(flux)), slots))
    with ProcessPoolExecutor(max_workers=4) as executor:
        results = list(executor.map(evaluate, jobs))
    plans = {name: p for name, p, _ in results}
    slots = {name: s for name, _, s in results}

    admitted = []
    excluded = []
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
        if name in slots and not plans.get(name):
            why.append(f"No slot reaches continuum S/N {SETTINGS['goal_snr']} per Angstrom at {SETTINGS['exposures']}x{SETTINGS['seconds_each']} s")
        if why:
            excluded.append(dict(name=name, protected=name in protected, reasons=why))
        else:
            admitted.append(name)
    missing = set(protected) - set(admitted)
    assert not missing, f'Protected primaries need explicit review: {sorted(missing)}'

    minutes = geo['minutes']
    visit = SETTINGS['visit_minutes']
    max_visits = (minutes - SETTINGS['reserve_minutes']) // visit
    choices = []
    for name in admitted:
        t = targets.loc[name]
        s = science.loc[name]
        # Scheduling utility: the existing review score decides which targets; the slot quality
        # (S/N relative to the target's own best slot) decides where. Not a changing-look probability.
        value = 20 + float(s.review_order_score) + 5 * bool(t.balmer_pair_in_range) + 5 * float(t.manifold_cl_neighbor_fraction)
        good = [v for v in slots[name].values() if v['meets_goal']]
        best = max(v['snr_per_angstrom'] for v in good)
        for v in good:
            quality = v['snr_per_angstrom'] / best
            choices.append(dict(name=name, tier=v['tier'], offset=v['offset'], duration=visit, start_pdt=v['start_pdt'],
                                end_pdt=v['end_pdt'], utility=value + SETTINGS['geometry_weight'] * quality - 1e-5 * v['offset'],
                                science_value=value, slot_quality=quality,
                                **{k: v[k] for k in ['snr_per_angstrom', 'snr_per_bin', 'airmass_mean', 'airmass_max',
                                                     'moon_deg', 'sky_V', 'seeing_arcsec']}))
    rows = []
    cols = []
    name_index = {n: minutes + i for i, n in enumerate(admitted)}
    count_row = minutes + len(admitted)
    for j, c in enumerate(choices):
        rows.extend(range(c['offset'], c['offset'] + c['duration']))
        cols.extend([j] * c['duration'])
        rows.extend([name_index[c['name']], count_row])
        cols.extend([j, j])
    matrix = coo_matrix((np.ones(len(rows)), (rows, cols)), shape=(count_row + 1, len(choices))).tocsc()
    lower = np.zeros(matrix.shape[0])
    upper = np.ones(matrix.shape[0])
    upper[count_row] = max_visits
    for name in protected:
        lower[name_index[name]] = 1
    result = milp(-np.array([c['utility'] for c in choices]), integrality=np.ones(len(choices)),
                  bounds=Bounds(0, 1), constraints=LinearConstraint(matrix, lower, upper),
                  options=dict(time_limit=60, mip_rel_gap=.001))
    assert result.x is not None, result.message
    selected = sorted((choices[i] for i in np.flatnonzero(result.x > .5)), key=lambda c: c['offset'])
    names = [c['name'] for c in selected]
    assert set(protected) <= set(names) and len(set(names)) == len(names) and len(names) <= max_visits
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
                  protected=sorted(protected), original_primaries=sorted(original),
                  promoted=[n for n in names if n not in protected],
                  pool=dict(september_eligible=len(september), with_slots=len(slots), admitted=len(admitted),
                            choices=len(choices), max_visits=int(max_visits)),
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
    show = ['rank', 'name', 'start_pdt', 'tier', 'airmass_mean', 'moon_deg', 'sky_V', 'seeing_arcsec', 'snr_per_angstrom', 'slot_quality', 'science_value']
    print(frame[show].round(2).to_string(index=False), flush=True)
    print('Protected:', len(protected), 'Added:', output['promoted'], 'Reserve minutes:', output['reserved']['minutes'], gaps, flush=True)
    print('Pool:', output['pool'], 'Solver:', output['solver'], flush=True)


if __name__ == '__main__':
    main()

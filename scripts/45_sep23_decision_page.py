"""Decision board for September 23: every candidate with a feasible visit, in observability order.

Builds a local (collaboration) page listing all Sep 23 targets that fit at least one full
16-minute visit at X<=1.8 with the Moon >=40 deg, numbered P01.. in order of their earliest
feasible visit start. Each card carries what a human needs to choose a sequence: airmass and
Moon separation through the half-night, observable minutes, magnitudes, predicted S/N per
Angstrom for the fixed setting in every candidate slot, manifold position and neighbour
fractions, the 40" field, light curves and every available spectrum. Reads the plan written
by 41_september_snr5_plan.py; writes nothing that the telescope reads.
"""
from pathlib import Path
import base64, importlib, json, re
import numpy as np
import pandas as pd
import astropy.units as u
from astropy.coordinates import AltAz, SkyCoord, get_body
from astropy.time import Time
from astropy.utils import iers

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'data/reselection_2026-09-20'
DEST = ROOT/'observing/sep23'
TZ = 'America/Los_Angeles'
iers.conf.auto_download = False
iers.conf.auto_max_age = None
OBS = importlib.import_module('05_observability')
STEP = 4  # minutes between plotted geometry samples


def payload(path):
    return json.loads(re.search(r'<script id="candidate-data" type="application/json">(.*?)</script>', path.read_text(), re.S).group(1))


def minutes_since(start, text):
    return int((pd.Timestamp(text, tz=TZ)-start).total_seconds()//60)


def runs_to_minutes(start, runs, end_key='end'):
    return [[minutes_since(start, r['start']), minutes_since(start, r[end_key])] for r in runs]


def main():
    plan = json.loads((DEST/'snr5_plan.json').read_text())
    S = plan['settings']
    opts = json.loads((OUT/'airmass_options.json').read_text())
    targets = pd.read_csv(OUT/'compact_review_objects.csv').set_index('name', drop=False)
    science = pd.read_csv(OUT/'three_night_review/science_and_sensitivity.csv').set_index('name')
    local = payload(OUT/'candidate_review_local.html')
    public_names = {t['name'] for t in payload(ROOT/'docs/index.html')['targets']}
    by_name = {t['name']: t for t in local['targets']}
    names = [n for n in plan['slots'] if n in by_name]
    assert len(names) == len(plan['slots']), 'Every planned target must be in the local explorer payload'
    notes = {e['name']: e['reasons'] for e in plan['admission_notes']}
    sequence = {c['name']: c for c in plan['sequence']}

    # Geometry across the whole first half-night at one-minute resolution.
    t0, t1, _, _ = OBS.night_window('2026-09-23', 'first')
    start = pd.Timestamp(t0.utc.iso, tz='UTC').tz_convert(TZ).floor('min')
    end = pd.Timestamp(t1.utc.iso, tz='UTC').tz_convert(TZ).ceil('min')
    minutes = int((end-start).total_seconds()//60)
    index = pd.DatetimeIndex([start+pd.Timedelta(minutes=m) for m in range(minutes+1)])
    times = Time(index.tz_convert('UTC').tz_localize(None).to_pydatetime())
    frame = AltAz(obstime=times, location=OBS.PALOMAR.location, pressure=0*u.hPa)
    coords = SkyCoord(targets.loc[names].ra.to_numpy()*u.deg, targets.loc[names].dec.to_numpy()*u.deg)
    aa = coords[:, None].transform_to(frame)
    moon = get_body('moon', times, OBS.PALOMAR.location).transform_to(frame)
    alt = aa.alt.deg
    airmass = aa.secz.value
    sep = aa.separation(moon).deg
    science_block = [minutes_since(start, S['science_start_pdt']), minutes_since(start, S['science_end_pdt'])]

    rows = []
    for i, name in enumerate(names):
        t = targets.loc[name]
        s = science.loc[name]
        w = plan['windows'][name]
        slots = plan['slots'][name]
        fallback = next((v for v in opts['windows'][name] if v['night'] == 'sep23'), {}).get('tier_ranges', {}).get('fallback', [])
        good = (alt[i] > 0) & (sep[i] >= S['moon_min_deg'])
        x = airmass[i]
        visible = good & (x <= 2.0)
        xmin_index = int(np.argmin(np.where(visible, x, np.inf))) if visible.any() else None
        slot_list = sorted(slots.values(), key=lambda v: v['offset'])
        ok = [v for v in slot_list if v['meets_goal']]
        best = max(slot_list, key=lambda v: v['snr_per_angstrom'])
        first_extended = min((minutes_since(start, r['start']) for r in w['extended']), default=10**6)
        first_preferred = min((minutes_since(start, r['start']) for r in w['preferred']), default=10**6)
        seq = sequence.get(name)
        hbeta_nm = 486.27*(1+float(t.z))
        decision = dict(
            role=t.pool_role, in_sequence=(seq['rank'] if seq else None), sequence_start=(seq['start_pdt'] if seq else None),
            protected=name in plan['protected'], original=name in plan['original_primaries'], added=name in plan['promoted'],
            public_identity=name in public_names, private_reference=bool(s.reference_private) if pd.notna(s.reference_private) else False,
            r_latest_ge19=bool(pd.notna(s.ztf_r_latest180_mag) and s.ztf_r_latest180_mag >= 19),
            admitted=name not in notes, admission_reasons=notes.get(name, []),
            windows=dict(preferred=w['preferred'], extended=w['extended'], fallback=[dict(start=r['start'], end=r['end']) for r in fallback]),
            windows_min=dict(preferred=runs_to_minutes(start, w['preferred']), extended=runs_to_minutes(start, w['extended'])),
            first_start_pdt=(min(w['extended'], key=lambda r: r['start'])['start'] if w['extended'] else None),
            first_preferred_pdt=(min(w['preferred'], key=lambda r: r['start'])['start'] if w['preferred'] else None),
            minutes_x15=int((good & (x <= 1.5)).sum()), minutes_x18=int((good & (x <= 1.8)).sum()), minutes_x20=int(visible.sum()),
            xmin=(float(x[xmin_index]) if xmin_index is not None else None),
            xmin_pdt=(index[xmin_index].strftime('%H:%M') if xmin_index is not None else None),
            moon_min=float(sep[i][visible].min()) if visible.any() else None, moon_max=float(sep[i][visible].max()) if visible.any() else None,
            r_archival=float(t.r_planning), r_latest=(float(s.ztf_r_latest180_mag) if pd.notna(s.ztf_r_latest180_mag) else None),
            g_latest=(float(s.ztf_g_latest180_mag) if pd.notna(s.ztf_g_latest180_mag) else None),
            continuum_AB=float(s.continuum_AB), z=float(t.z), hbeta_nm=hbeta_nm, channel=best.get('channel', plan['plans'].get(name, {}).get('preferred', plan['plans'].get(name, {}).get('extended', {})).get('channel')),
            halpha_in_range=bool(t.balmer_pair_in_range), hbeta_in_ngps=bool(t.Hb_in_NGPS) if pd.notna(t.Hb_in_NGPS) else None,
            snr_best=float(best['snr_per_angstrom']), snr_best_pdt=best['start_pdt'][11:], snr_best_airmass=float(best['airmass_mean']), snr_best_sky=float(best['sky_V']),
            snr_min=float(min(v['snr_per_angstrom'] for v in slot_list)), n_slots=len(slot_list), n_slots_ok=len(ok),
            cl_fraction=float(t.manifold_cl_neighbor_fraction), on_fraction=float(t.manifold_on_neighbor_fraction), off_fraction=float(t.manifold_off_neighbor_fraction),
            manifold_direction=(str(t.manifold_direction) if pd.notna(t.manifold_direction) else ''), manifold_interest=(str(t.manifold_interest) if pd.notna(t.manifold_interest) else ''),
            review_score=float(s.review_order_score), science_question=str(s.science_question), state_classification=str(s.state_classification),
            known_state_status=(str(s.known_state_status) if pd.notna(s.known_state_status) else ''),
            optical_trigger=bool(s.post_spectrum_optical_trigger) if pd.notna(s.post_spectrum_optical_trigger) else False,
            trigger_direction=(str(s.trigger_direction) if pd.notna(s.trigger_direction) else ''), ir_flag=bool(s.post_spectrum_ir_flag) if pd.notna(s.post_spectrum_ir_flag) else False,
            ztf_r_change=(float(s.ztf_r_change_after_reference) if pd.notna(s.ztf_r_change_after_reference) else None),
            neowise_change=(float(s.neowise_change_after_reference_mag) if pd.notna(s.neowise_change_after_reference_mag) else None),
            reference_date=str(s.reference_date), n_spec=int(t.sdssv_n_any_epochs) if pd.notna(t.sdssv_n_any_epochs) else None,
            field_status=str(t.field_status), field_notes=str(t.field_notes) if pd.notna(t.field_notes) else '',
            curve=[[m, (round(float(x[m]), 3) if alt[i][m] > 0 and x[m] < 3.5 else None), round(float(sep[i][m]), 1)] for m in range(0, minutes+1, STEP)],
            slots=[[minutes_since(start, v['start_pdt']), round(v['snr_per_angstrom'], 2), 1 if v['tier'] == 'preferred' else 2, bool(v['meets_goal'])] for v in slot_list],
            order_key=[first_extended, first_preferred, float(t.ra)])
        rows.append((name, decision))
    rows.sort(key=lambda r: r[1]['order_key'])

    page_targets = []
    table = []
    for rank, (name, d) in enumerate(rows, 1):
        d['rank'] = rank
        d['code'] = f'P{rank:02}'
        t = json.loads(json.dumps(by_name[name]))
        for key in ['nights', 'exposure_plans']:
            t.pop(key, None)
        screen = t.get('neighbour_screen') or {}
        t['neighbour_screen'] = dict(status=screen.get('status'), flags=screen.get('flags', []), unwise=screen.get('unwise'))
        path = ROOT/'data/cutouts'/f'{name}_sdss.jpg'
        meta = path.with_suffix('.json')
        if path.exists() and meta.exists() and json.loads(meta.read_text()).get('field_arcsec') == 40:
            t['cut40'] = 'data:image/jpeg;base64,'+base64.b64encode(path.read_bytes()).decode()
        else:
            t['cut40'] = None
        t['decision'] = d
        page_targets.append(t)
        table.append(dict(code=d['code'], name=name, role=d['role'], in_current_sequence=d['in_sequence'], protected=d['protected'],
                          first_start_x18=d['first_start_pdt'], first_start_x15=d['first_preferred_pdt'],
                          minutes_x15=d['minutes_x15'], minutes_x18=d['minutes_x18'], xmin=round(d['xmin'], 3) if d['xmin'] else None, xmin_pdt=d['xmin_pdt'],
                          moon_min=round(d['moon_min'], 1) if d['moon_min'] else None, moon_max=round(d['moon_max'], 1) if d['moon_max'] else None,
                          r_archival=round(d['r_archival'], 2), r_latest=d['r_latest'], continuum_AB=round(d['continuum_AB'], 2), z=round(d['z'], 3),
                          hbeta_nm=round(d['hbeta_nm'], 1), halpha_in_range=d['halpha_in_range'],
                          snr_best=round(d['snr_best'], 1), snr_best_pdt=d['snr_best_pdt'], snr_min=round(d['snr_min'], 1), slots_meeting_floor=d['n_slots_ok'], slots=d['n_slots'],
                          cl_fraction=round(d['cl_fraction'], 3), on_fraction=round(d['on_fraction'], 3), off_fraction=round(d['off_fraction'], 3),
                          review_score=round(d['review_score'], 1), optical_trigger=d['optical_trigger'], trigger_direction=d['trigger_direction'], ir_flag=d['ir_flag'],
                          public_identity=d['public_identity'], private_reference=d['private_reference'], r_latest_ge19=d['r_latest_ge19'],
                          admitted=d['admitted'], admission_reasons=' | '.join(d['admission_reasons']), reference_date=d['reference_date'],
                          state=d['state_classification'], field_status=d['field_status'], science_question=d['science_question']))
    pd.DataFrame(table).to_csv(DEST/'sep23_decision_table.csv', index=False)

    unreachable = [n for n in targets.index if any(v['night'] == 'sep23' for v in opts['windows'].get(n, [])) and n not in plan['slots']]
    value = dict(version=local['version'], generated=pd.Timestamp.utcnow().isoformat(), access=local['access'], nights=local['nights'], manifold=local['manifold'],
                 night=dict(start_pdt=start.strftime('%Y-%m-%d %H:%M'), end_pdt=end.strftime('%Y-%m-%d %H:%M'), minutes=minutes, science_block=science_block,
                            moon_percent=local['nights']['sep23']['moon_percent'], step=STEP),
                 packet=dict(settings=S, sequence=plan['sequence'], reserved=plan['reserved']),
                 unreachable=unreachable, targets=page_targets)
    OLD = importlib.import_module('07_make_webpage')
    charts = OLD.TEMPLATE[OLD.TEMPLATE.index('function mjdToYear'):OLD.TEMPLATE.index('/* ---------- manifold thumbnail')]
    charts += '\n'+(ROOT/'web/candidate_spectra.js').read_text()
    template = (ROOT/'web/sep23_decision_template.html').read_text()
    encoded = json.dumps(value, separators=(',', ':'), allow_nan=False).replace('<', '\\u003c')
    page = template.replace('__PAYLOAD__', encoded).replace('__CHART_FUNCTIONS__', charts)
    (DEST/'sep23_decision_local.html').write_text(page)
    print(f"{len(page_targets)} candidates; page {len(page)/1e6:.1f} MB; unreachable at X<=1.8: {unreachable}")
    print(pd.DataFrame(table)[['code', 'name', 'role', 'in_current_sequence', 'first_start_x18', 'minutes_x15', 'xmin', 'xmin_pdt', 'moon_min', 'r_archival', 'continuum_AB', 'snr_best', 'cl_fraction', 'review_score']].to_string(index=False))


if __name__ == '__main__':
    main()

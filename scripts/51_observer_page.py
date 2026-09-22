"""The single observer page, in the dark night-sheet style.

Reads the candidate payloads written by 16_candidate_webpage.py (which now only supplies data),
the September plan and packet, and renders one page per access level: the local collaboration
copy with complete spectra and the backup CSV, and the public copy for GitHub Pages. Sections:
About this run, the September 23 sequence (timeline, table, backups, CSV downloads), the
candidate pool for all three nights (visibility chart, summary table) and one card per target
with geometry, light curves, spectra, image and manifold position. Stable target numbers are
assigned once, by right ascension, and hold across nights.

Usage: python scripts/51_observer_page.py
"""
import base64, json
from pathlib import Path
import importlib
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'data/reselection_2026-09-20'
DEST_LOCAL = OUT/'observer_page_local.html'
DEST_PUBLIC = ROOT/'docs/index.html'


def cut40(name):
    path = ROOT/'data/cutouts'/f'{name}_sdss.jpg'
    meta = path.with_suffix('.json')
    if path.exists() and meta.exists() and json.loads(meta.read_text()).get('field_arcsec') == 40:
        return 'data:image/jpeg;base64,'+base64.b64encode(path.read_bytes()).decode()
    return None


def sep23_slots(public=False):
    """Per-slot predicted S/N from the September planner, as minutes from the half-night start."""
    plan_path = ROOT/'observing/sep23/snr5_plan.json'
    if not plan_path.exists():
        return {}
    plan = json.loads(plan_path.read_text())
    nights = json.loads((OUT/'candidate_payload_local.json').read_text())['nights']
    start = pd.Timestamp((nights['sep23']['start_mjd']-40587)*86400, unit='s', tz='UTC').tz_convert('America/Los_Angeles')
    out = {}
    for name, table in plan['slots'].items():
        if public and (not plan['plans'].get(name) or any(p.get('reference_private') for p in plan['plans'][name].values())):
            continue
        rows = []
        for v in table.values():
            m = (pd.Timestamp(v['start_pdt'], tz='America/Los_Angeles')-start).total_seconds()/60
            rows.append([round(m, 1), round(v['snr_per_angstrom'], 2), 1 if v['tier'] == 'preferred' else 2, bool(v['meets_goal']), v['duration'], v['exposures']])
        out[name] = sorted(rows)
    return out


PUBLIC_SCIENCE = ['ztf_g_latest180_mag', 'ztf_r_latest180_mag', 'ztf_last_date', 'state_classification', 'field_notes']
PUBLIC_SCIENCE_IF_PUBLIC_REFERENCE = ['post_spectrum_optical_trigger', 'trigger_direction', 'post_spectrum_ir_flag', 'ztf_r_change_after_reference', 'neowise_change_after_reference_mag', 'continuum_AB', 'reference_date', 'review_order_score', 'hbeta_A', 'channel', 'science_question',
                                      'snr300_x1p3_sky18p5', 'snr300_x1p5_sky18', 'snr300_x1p8_sky18']


def prepare(payload, slots, public_plans=None, public_science=None):
    targets = sorted(payload['targets'], key=lambda t: t['ra'])
    if public_plans is not None:
        # Public copy: predicted S/N per tier and science fields only where the continuum reference is public;
        # photometric triggers and ZTF brightness are public data for every target.
        for t in targets:
            plans = public_plans.get(t['name'])
            if plans:
                t['exposure_plans'] = plans
            sci = (public_science or {}).get(t['name'])
            if sci:
                t['science'] = sci
    for i, t in enumerate(targets, 1):
        t['code'] = f'{i:03d}'
        t['cut40'] = cut40(t['name'])
        t['sep23_slots'] = slots.get(t['name'], [])
        screen = t.get('neighbour_screen') or {}
        t['neighbour_screen'] = dict(status=screen.get('status'), flags=screen.get('flags', []), unwise=screen.get('unwise'))
    payload['targets'] = targets
    return payload


def main():
    OLD = importlib.import_module('07_make_webpage')
    charts = OLD.TEMPLATE[OLD.TEMPLATE.index('function mjdToYear'):OLD.TEMPLATE.index('/* ---------- manifold thumbnail')]
    charts += '\n'+(ROOT/'web/candidate_spectra.js').read_text()
    template = (ROOT/'web/observer_page_template.html').read_text()
    slots = sep23_slots()
    public_slots = sep23_slots(public=True)
    local = json.loads((OUT/'candidate_payload_local.json').read_text())
    public_plans = {}
    for t in local['targets']:
        plans = t.get('exposure_plans') or {}
        if plans and not any(v.get('reference_private') for v in plans.values()):
            public_plans[t['name']] = {k: dict(airmass=v['airmass'], exposures=v['exposures'], seconds_each=v['seconds_each'], visit_minutes=v['visit_minutes'],
                                                predicted_snr_per_angstrom=v['predicted_snr_per_angstrom'], status=v['status']) for k, v in plans.items()}
    public_science = {}
    for t in local['targets']:
        sci = t.get('science') or {}
        if not sci:
            continue
        keep = {k: sci.get(k) for k in PUBLIC_SCIENCE if k in sci}
        if not sci.get('reference_private'):
            keep.update({k: sci.get(k) for k in PUBLIC_SCIENCE_IF_PUBLIC_REFERENCE if k in sci})
        keep['reference_private'] = bool(sci.get('reference_private'))
        public_science[t['name']] = {k: (None if isinstance(v, float) and v != v else v) for k, v in keep.items()}
    for source, dest, private in [(OUT/'candidate_payload_local.json', DEST_LOCAL, True), (OUT/'candidate_payload_public.json', DEST_PUBLIC, False)]:
        payload = prepare(json.loads(source.read_text()), slots if private else public_slots, None if private else public_plans, None if private else public_science)
        encoded = json.dumps(payload, separators=(',', ':'), allow_nan=False).replace('<', '\\u003c')
        page = template.replace('__PAYLOAD__', encoded).replace('__CHART_FUNCTIONS__', charts)
        if not private:
            assert 'proprietary' not in page and 'SDSS-V internal' not in page
            assert not payload['backups'] and 'sep23_backups_ngps.csv' not in payload['files']
        dest.write_text(page)
        print(f"{dest.relative_to(ROOT)}: {len(payload['targets'])} targets, {len(payload['sep23_sequence'])} sequence entries, {len(page)/1e6:.1f} MB")


if __name__ == '__main__':
    main()

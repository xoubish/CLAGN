"""Public, spectroscopically established quasars for exact September replacement visits."""
import importlib
from pathlib import Path

import numpy as np
import pandas as pd

from spectral_utils import accepted_reference, archival_records

POLICY = dict(count_per_primary=2, moon_min_deg=40, airmass_max=1.5,
              r_limit=19, strict_limits=True, match_primary_exposures=True,
              brightness='Latest available ZTF r (last 180 days of its light curve); archival r if unavailable',
              ranking='Selected literature/Zeltyn region first; then labelled on/off neighbour fraction, Balmer coverage and brightness',
              public_spectral_reference=True)
QSO_ORIGINS = {'expanded_DR16_QSO', 'full_DR16Q_catalog'}


def catalog_evidence(target, science):
    """Do not treat a generic AGN/galaxy identity as proof of a quasar classification."""
    origins = set(str(target.get('origins', '')).split(';')) | {target.get('origin')}
    if not origins & QSO_ORIGINS:
        return None
    recent = science.get('ztf_r_latest180_mag')
    use_recent = recent is not None and np.isfinite(recent)
    mag = recent if use_recent else target.get('r_planning')
    if mag is None or not np.isfinite(mag) or not mag < POLICY['r_limit']:
        return None
    region = target.get('review_region')
    in_region = region in {'literature turn-on/off', 'Zeltyn'}
    fraction = target.get('manifold_cl_neighbor_fraction', 0)
    fraction = float(fraction) if pd.notna(fraction) else 0.
    # Lexicographic region preference cannot be outweighed by brightness or S/N.
    score = 100*fraction + 10*bool(target.get('balmer_pair_in_range', False)) + 5*(19-mag)
    return dict(quasar_basis=';'.join(sorted(origins & QSO_ORIGINS)), r_mag=float(mag),
                r_source='ZTF latest available 180-day median' if use_recent else target.get('r_source', 'archival r'),
                r_date=science.get('ztf_last_date') if use_recent else None,
                region=region, in_selected_region=in_region, manifold_fraction=fraction,
                ranking_score=score)


def candidates(targets, science, public_names, primary_names):
    eligible = []
    for target in targets.to_dict('records'):
        name = target['name']
        if name in primary_names or name not in public_names or target['field_status'] != 'clear':
            continue
        evidence = catalog_evidence(target, science.loc[name])
        if evidence is None:
            continue
        records = [r for r in archival_records(name) if not r.get('proprietary')]
        reference = accepted_reference(target, records)
        if reference is None:
            continue
        record, mag, flux = reference
        evidence['reference_mjd'] = float(record['mjd'])
        eligible.append(dict(target=target, evidence=evidence,
                             reference=(float(record['mjd']), float(mag), False, float(flux))))
    return sorted(eligible, key=lambda c: (not c['evidence']['in_selected_region'],
                                          -c['evidence']['ranking_score'], c['target']['name']))


def select(eligible, primary, geometry, indices, make_visit, cache_dir):
    """Recompute depth for the primary's exposure count and geometry over its entire visit."""
    planner = importlib.import_module('41_september_snr5_plan')
    nexp = primary['plan']['exposures']
    duration = planner.visit_minutes(nexp)
    assert duration == primary['plan']['visit_minutes']
    start = primary['start_pdt']
    chosen = []
    for candidate in eligible:
        target = candidate['target']
        name = target['name']
        slots, windows = planner.candidate_slots(geometry, indices[name], duration)
        slot = next((s for s in slots if s['start_pdt'] == start), None)
        if slot is None or not slot['airmass_max'] < 1.5 or not slot['moon_deg'] > 40:
            continue
        _, plans, evaluated = planner.evaluate(
            (name, target['z'], candidate['reference'], nexp, [slot]), cache_dir=Path(cache_dir))
        visit = make_visit(target, plans, evaluated, windows, start)
        if visit is None or not visit['airmass_max_actual'] < 1.5 or not visit['moon_min'] > 40:
            continue
        assert visit['end_utc'] == primary['end_utc']
        visit['eligibility'] = candidate['evidence']
        visit['science_question'] = 'Compare broad Balmer emission with the public archival spectrum; selected using the W1 manifold prior'
        chosen.append((visit, target))
        if len(chosen) == POLICY['count_per_primary']:
            break
    return chosen

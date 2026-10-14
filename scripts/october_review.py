"""Evidence and explicit decisions for the October 14 observing review."""
from pathlib import Path
import pandas as pd
from astropy.time import Time

ROOT = Path(__file__).resolve().parents[1]
REDUCTION = ROOT / 'sep23_data/reduction_20260924'
REASSESSMENT = ROOT / 'jwst_proposal/review/sample_reassessment'
REVIEW_DATE = '2026-10-14'


def september_performance():
    summary = pd.read_csv(REDUCTION / 'products_p330e/target_summary.csv')
    raw = pd.read_csv(REDUCTION / 'audit/raw_frames.csv')
    raw = raw[raw.name.isin(summary.NAME)].sort_values('mjd')
    visits = raw.groupby('name', sort=False).agg(
        first_mjd=('mjd', 'min'), last_mjd=('mjd', 'max'),
        integration_seconds=('shutter_seconds', 'sum'), exposures=('mjd', 'size'))
    visits['start_to_next_start_minutes'] = (visits.first_mjd.shift(-1) - visits.first_mjd) * 1440
    visits['nonintegration_minutes_to_next_target'] = visits.start_to_next_start_minutes - visits.integration_seconds / 60
    metrics = dict(targets=len(summary), exposures=int(summary.NEXP.sum()),
        integration_minutes=float(summary.EXPTIME_S.sum()/60),
        elapsed_science_minutes=float((raw.mjd.iloc[-1]-raw.mjd.iloc[0])*1440+raw.shutter_seconds.iloc[-1]/60),
        median_overhead_minutes=float(visits.nonintegration_minutes_to_next_target.median()),
        p90_overhead_minutes=float(visits.nonintegration_minutes_to_next_target.quantile(.9)),
        snr_median_per_native_pixel={c:float(summary['SNR_'+c].median()) for c in 'UGRI'})
    return visits.reset_index(), metrics


def apply_review(targets):
    """Keep proposal membership; use the later history audit to refine questions."""
    audit = pd.read_csv(REASSESSMENT / 'current_target_assessment.csv').set_index('name')
    bins = pd.read_csv(ROOT / 'jwst_proposal/review/history_audit/binned_measurements.csv')
    for name in targets.index[targets.is_jwst]:
        r = targets.loc[name]
        if name in audit.index:
            a = audit.loc[name]
            targets.loc[name, 'science_question'] = (
                f"JWST {r.jwst_id}: measure current Balmer and continuum state for the dust-history experiment. "
                + str(a.assessment) + ' ' + str(a.required_check))
        # Use dated survey photometry instead of integrated catalogue galaxy r.
        b = bins[(bins.internal_id == name) & (bins.band == 'r')].sort_values('mjd')
        if len(b):
            end = b.mjd.max()
            targets.loc[name, 'r_adopted'] = float(b[b.mjd >= end-180].mag.median())
            targets.loc[name, 'brightness_source'] = 'JWST history audit: cached r 90-day bins, last 180 days'
            targets.loc[name, 'brightness_epoch'] = Time(end, format='mjd').isot[:10]
    targets['exposures'] = 2
    targets['visit_minutes'] = 14
    targets['overhead_minutes'] = 4
    targets['protected_followup'] = False
    targets.loc['P8544', ['eligible_science', 'protected_followup']] = True
    targets.loc['P8544', ['exposures', 'visit_minutes', 'overhead_minutes']] = [3, 20, 5]
    targets.loc['P8544', 'priority_reason'] = 'September profile follow-up; match September exposure depth'
    targets.loc['P8544', 'science_question'] = (
        'Repeat the September broad-Hbeta weakening candidate at 3x300 s; '
        'test persistence with matched continuum, host and FeII treatment. '
        'Preliminary profile change, not a confirmed changing-look transition.')
    targets.loc['P8544', 'science_score'] = 200
    return targets


def write_review(dest, targets, queue, windows, write_csv):
    visits, metrics = september_performance()
    write_csv(visits, dest / 'september_performance.csv')
    before = pd.read_csv(dest / 'before_review_20261014.csv').rename(columns=str.lower)
    old = before.set_index('name')
    selected = queue.set_index('name')
    decisions = []
    for name, r in targets.iterrows():
        chosen = name in selected.index
        status = ('retained' if name in old.index else 'added') if chosen else ('removed' if name in old.index else 'not selected')
        reason = r.priority_reason if chosen else (
            'September spectrum retained; no October window under adopted cuts' if r.is_jwst else
            'September spectrum already obtained; no protected October follow-up' if r.observed_sep23 else
            'Fails preparation, brightness or duplicate screen' if not r.eligible_science else
            'No complete feasible visit' if not (windows.name == name).any() else
            'Lower joint science/geometry priority within 40 targets per night')
        decisions.append(dict(name=name,decision=status,previous_night=old.loc[name,'night'] if name in old.index else '',
            scheduled_night=selected.loc[name,'night'] if chosen else '',reason=reason))
    write_csv(pd.DataFrame(decisions), dest / 'selection_decisions.csv')
    observed = targets[targets.observed_sep23].reset_index()[['name','jwst_id','protected_followup']]
    notes = {
        'P9694':'Protect Oct 27: cleanest September strengthening lead; later alerts show decline from the high state, not a monotonic rise.',
        'P10381':'Protect Oct 27: published fading history; test recovery without claiming a repeat transition.',
        'P9506':'Protect Oct 27: JWST infrared decline / optical recovery history.',
        'P7837':'Retain one JWST history-comparison repeat; no new transition established in September.',
        'P11530':'Retain September spectrum: Moon-blocked on both October nights; longer flare/fading history remains valuable.',
        'P8544':'Add Oct 27 at 3x300 s: candidate weakening; match September depth at lower airmass and larger Moon separation.',
        'P9227':'No automatic repeat: telluric-sensitive Halpha lead needs reanalysis; weaker priority than P8544.',
        'P2921':'Defer: no full October visit under normal airmass/Moon cuts; profile interpretation also needs telluric/FeII checks.',
        'P1823':'Defer: short setting window conflicts with opening calibration; September already has 4x300 s. Valuable continuum history, not a demonstrated turn-on.'}
    observed['decision'] = observed.name.map(lambda n: notes.get(n, 'Retain September spectrum; spend October time on new objects and protected follow-up.'))
    write_csv(observed, dest / 'september_followup_decisions.csv')
    replacements = pd.read_csv(REASSESSMENT / 'feasibility/ngps_visibility_summary.csv')
    replacements = replacements[~replacements.name.isin(targets.index)].copy()
    replacements['decision'] = 'Development reserve only: not in submitted 24; acquisition field and usable baseline not verified in prepared packet'
    replacements.loc[replacements.name == 'P20778', 'decision'] += '; cached r=19.37 also fails ordinary r<19 cut'
    write_csv(replacements, dest / 'development_reserves.csv')
    return metrics

"""Export the screened September working set from the local review payload.

Ordering describes observing flexibility, not discovery probability. The output
stays in the ignored research directory because spectral counts include internal
epochs. No observing sequence or exposure time is assigned.
"""
from pathlib import Path
import json
import re
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'data/reselection_2026-09-20'


def main():
    page = (OUT / 'candidate_review_local.html').read_text()
    payload = json.loads(re.search(r'<script id="candidate-data" type="application/json">(.*?)</script>', page, re.S).group(1))
    rows = []
    for target in payload['targets']:
        window = next((w for w in target['nights'] if w['night'] == 'sep23'), None)
        if window is None or target['field_status'] != 'clear':
            continue
        duration = window['longest_minutes']
        rows.append(dict(name=target['name'], pool_role=target['pool_role'],
                         ra_deg=target['ra'], dec_deg=target['dec'], z=target['z'],
                         historical_r=target['rmag'], brightness_source=target['r_source'],
                         window_tier_minutes=120 if duration >= 120 else 90 if duration >= 90 else 60,
                         longest_window_minutes=duration,
                         minutes_airmass_le1p3=window['minutes_airmass_le1p3'],
                         minimum_airmass=window['min_airmass'], minimum_moon_deg=window['moon_min'],
                         known_state_status=target['status'],
                         identified_spectral_dates=target['n_spec'],
                         loaded_spectral_dates=len((target['spec'] or {}).get('epochs', [])),
                         field_screen='No catalogue flags; visual and slit-PA review still required'))
    table = pd.DataFrame(rows).sort_values(
        ['pool_role', 'minutes_airmass_le1p3', 'longest_window_minutes', 'historical_r'],
        ascending=[True, False, False, True])
    table.to_csv(OUT / 'september_operational_shortlist.csv', index=False)
    main_pool = table[table.pool_role.eq('manifold')]
    reserves = table[table.pool_role.eq('reserve')]
    primary = main_pool[main_pool.window_tier_minutes.eq(120)]
    counts = {n: int(main_pool.longest_window_minutes.ge(n).sum()) for n in [60, 90, 120]}
    columns = ['name', 'historical_r', 'z', 'longest_window_minutes', 'minutes_airmass_le1p3', 'minimum_moon_deg', 'identified_spectral_dates']
    report = f'''# September 23 operational working set

For the **first half of the local night September 23, 2026**, {payload['nights']['sep23']['window']}.
The half-night boundary is the midpoint between astronomical dusk and dawn; confirm the observatory handover time before producing a final sequence.

Moon separation >=40 degrees and airmass <=1.5 must hold simultaneously throughout each accepted interval. Five-minute sampling accepts intervals only when both endpoints pass. Historical r<=18.5 for the main pool; r<=18 for the separate reserves. The broader projected parent was revisited, rather than restricting selection to the previous compact pool. Additional parent W1 processing is still underway, so the projected pool is not a completeness claim.

After catalogue screening, the main September pool has **{counts[120]} objects with >=120 minutes**, **{counts[90]} with >=90 minutes**, and **{counts[60]} with >=60 minutes**. These are nested counts. The public page defaults to the two-hour main set. There are **{len(reserves)} prepared reserve quasars**, outside the two candidate manifold regions, with >=120 minutes.

## Main two-hour working set

Ordered by time at airmass <=1.3, then window length and historical brightness. This is an operational order, not a science ranking.

{primary[columns].round(2).to_markdown(index=False)}

## Bright spectroscopic quasar reserves

These objects come from the public DR16 parent selected as sciencePrimary QSO with zWarning=0. Their spectra and manifold membership can support comparison observations. Being outside the candidate regions does not prove that they are nonvariable or never change state.

{reserves[columns].round(2).to_markdown(index=False)}

## Screening and remaining scientific decisions

SDSS neighbours were queried within 60 arcsec and Gaia DR3 within 120 arcsec, with available Gaia proper motions propagated to the observing epoch. Comparisons use the same photometric band. Inspection flags cover close companions, comparable or brighter neighbours near the slit, and very bright stars farther away. unWISE W1 requires a matched source, fracflux>=0.8, and zero unWISE W1 flags. These thresholds are inspection heuristics, not measured contamination limits; galaxy deblending and the planned slit position angle still need visual review. A pass does not establish uncontaminated nuclear W1 variability.

All priorities must compare dated W1/optical behaviour with the latest usable broad-line state. The manifold alone is not a transition clock. A bright reserve should enter only after the principal candidates or to answer a defined comparison question. No current on/off state, discovery probability, line S/N, exposure, or final observing sequence is assigned here. Brightness is historical, and ZTF refreshes are still in progress.

Spectral counts in this local report include available internal inventory dates; they do not imply every epoch passes broad-line quality checks. Keep this report local.
'''
    (OUT / 'SEPTEMBER_SELECTION.md').write_text(report)
    print(json.dumps(dict(main_counts=counts, prepared_reserves=len(reserves), total_screened_rows=len(table)), indent=2))


if __name__ == '__main__':
    main()

"""Full-pool NGPS target lists: every candidate observable on each night, same setting.

For each run night, one CSV in the documented NGPS target-list format containing every prepared
candidate with an observing window that night, in the order they first become observable. The
Note marks September primaries, backups, comparison reserves and plain pool members so rows
already observed can be skipped; the Comment carries the windows, predicted S/N at 2x300 s,
brightness, redshift, review score, triggers, the science question and privacy flags. These are
reserves beyond the curated primaries and backups; load rows deliberately, never the whole file.

Usage: python scripts/52_pool_csv.py
Writes observing/pool/ngps_pool_<night>.csv and observing/pool/pool_summary.csv.
"""
import importlib, json
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'data/reselection_2026-09-20'
DEST = ROOT/'observing/pool'
PACK = importlib.import_module('40_september_observing_packet')
NIGHT_LABEL = {'sep23': 'Sep 23 2026 first half', 'oct26': 'Oct 26 2026', 'oct27': 'Oct 27 2026'}


def main():
    DEST.mkdir(parents=True, exist_ok=True)
    targets = pd.read_csv(OUT/'compact_review_objects.csv').set_index('name', drop=False)
    science = pd.read_csv(OUT/'three_night_review/science_and_sensitivity.csv').set_index('name')
    opts = json.loads((OUT/'airmass_options.json').read_text())
    packet = json.loads((ROOT/'observing/sep23/packet.json').read_text())
    S = packet['settings']
    plan = json.loads((ROOT/'observing/sep23/snr5_plan.json').read_text())
    public_names = {t['name'] for t in json.loads((OUT/'candidate_payload_public.json').read_text())['targets']}
    primaries = {v['name']: v for v in packet['primaries']}
    backups = {}
    for b in packet['backups']:
        backups.setdefault(b['name'], []).append(b['replaces'])
    order = {n: i for i, n in enumerate(sorted(targets.index, key=lambda n: targets.loc[n].ra), 1)}  # stable page numbers
    summary = []
    for night in ['sep23', 'oct26', 'oct27']:
        rows = []
        for name in targets.index:
            w = next((x for x in opts['windows'].get(name, []) if x['night'] == night), None)
            if not w:
                continue
            t = targets.loc[name]
            s = science.loc[name]
            pref, ext, fb = (w['tier_ranges'].get(k, []) for k in ['preferred', 'extended', 'fallback'])
            if not (pref or ext or fb):
                continue
            first = (ext or fb)[0]
            plans = opts['plans'].get(name, {})
            snr15 = plans.get('preferred', {}).get('predicted_snr_per_angstrom')
            snr18 = plans.get('extended', {}).get('predicted_snr_per_angstrom')
            best = max([v for v in [snr15, snr18] if v is not None], default=None)
            if name in primaries and night == 'sep23':
                role = f"PRIMARY P{primaries[name]['rank']:02d}"
            elif name in backups and night == 'sep23':
                role = 'BACKUP'
            elif t.pool_role == 'reserve':
                role = 'RESERVE'
            else:
                role = 'POOL'
            code = f'{order[name]:03d}'
            span = lambda runs: '; '.join(f"{r['start']}-{r['end']}" for r in runs)
            windows = ' / '.join(f"{k} {span(runs)}" for k, runs in [('X<=1.5', pref), ('X<=1.8', ext)] if runs) or f"X<=2.0 only {span(fb)}"
            sep23_extra = ''
            if night == 'sep23' and name in plan['slots']:
                sl = max(plan['slots'][name].values(), key=lambda v: v['snr_per_angstrom'])
                sep23_extra = f"; best slot S/N {sl['snr_per_angstrom']:.1f} at {sl['start_pdt'][11:]} PDT"
            triggers = ', '.join(x for x in [('optical ' + str(s.trigger_direction)) if bool(s.post_spectrum_optical_trigger) else '', 'IR flag' if bool(s.post_spectrum_ir_flag) else ''] if x) or 'none'
            flags = '; '.join(x for x in [
                'identity not public' if name not in public_names else '',
                'private reference spectrum' if bool(s.reference_private) else '',
                'latest r >= 19' if (pd.notna(s.ztf_r_latest180_mag) and s.ztf_r_latest180_mag >= 19) else '',
                'below S/N floor at 2x300 s: consider a third exposure' if (best is not None and best < S['goal_snr']) else '',
                f"field: {t.field_status}" if t.field_status != 'clear' else ''] if x)
            fmt1 = lambda v: f'{v:.1f}' if v is not None and v == v else 'n/a'
            fmt2 = lambda v: f'{v:.2f}' if v is not None and v == v else 'n/a'
            parts = [f"{role} for {NIGHT_LABEL[night]}", f"observable {windows} PDT", f"first full visit from {first['start']}",
                     f"predicted continuum S/N per Angstrom at {S['exposures']}x{S['seconds_each']} s: {fmt1(snr15)} at X<=1.5 and {fmt1(snr18)} at X<=1.8{sep23_extra}",
                     f"role {'comparison reserve' if t.pool_role == 'reserve' else 'manifold candidate'} ({t.review_region})",
                     f"archival r {t.r_planning:.2f}", f"latest ZTF r {fmt2(s.ztf_r_latest180_mag)}", f"continuum AB {fmt2(s.continuum_AB)} at Hbeta",
                     f"z {t.z:.3f}", f"Hbeta channel {s.channel}", f"review score {fmt1(s.review_order_score)}", f"triggers {triggers}"]
            if name in backups and night == 'sep23':
                parts.append('Sep 23 backup for ' + ' '.join(backups[name]))
            comment = '; '.join(parts) + f". {s.science_question}. {t.field_notes if isinstance(t.field_notes, str) else ''} {('FLAGS: ' + flags + '.') if flags else ''} Pool row: load deliberately; skip anything already observed; same setting as the primaries."
            plan_row = dict(seconds_each=S['seconds_each'], exposures=S['exposures'], airmass=1.8 if (pref or ext) else 2.0)
            rows.append((first['start_utc'], PACK.ngps_row(dict(name=name, ra=float(t.ra), dec=float(t.dec)), plan_row, S, f"{code} {role} {first['start']}", comment)))
            summary.append(dict(night=night, code=code, name=name, role=role, first_start_pdt=first['start'], windows=windows, snr_x15=snr15, snr_x18=snr18,
                                pool_role=t.pool_role, r_archival=t.r_planning, ztf_r_latest=s.ztf_r_latest180_mag, continuum_AB=s.continuum_AB, z=t.z,
                                review_score=s.review_order_score, public_identity=name in public_names, private_reference=bool(s.reference_private), flags=flags))
        rows.sort(key=lambda x: (x[0], x[1]['name']))
        text = PACK.csv_text([r for _, r in rows])
        (DEST/f'ngps_pool_{night}.csv').write_text(text, encoding='ascii')
        print(f"{night}: {len(rows)} rows -> observing/pool/ngps_pool_{night}.csv")
    pd.DataFrame(summary).to_csv(DEST/'pool_summary.csv', index=False)
    (DEST/'README.md').write_text(
        "# Full-pool NGPS target lists\n\nOne CSV per run night with every prepared candidate that has an observing window that night, in the order they first become observable, "
        "all with the adopted setting (1.5 arcsec slit, 2x3 binning, 2x300 s, slit angle PA). The Note gives the stable page number, the role (PRIMARY Pxx, BACKUP, RESERVE, POOL) and the first start; "
        "the Comment gives the windows, predicted S/N per Angstrom at the airmass ceilings, brightness, redshift, review score, triggers, the science question and privacy flags.\n\n"
        "These are reserves beyond the curated primaries and backups in `observing/sep23/`. Load rows deliberately and skip anything already observed; never append a whole file to an automatic run. "
        "October rows are observability lists, not sequences. `pool_summary.csv` holds the same information as a table.\n")


if __name__ == '__main__':
    main()

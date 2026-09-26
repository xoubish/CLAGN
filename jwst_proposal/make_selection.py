"""Summarize the agreed Cycle 6 list without changing target membership.

The CSV is the authoritative working sample. This verifies bookkeeping, not
history classifications, event times, dust lags, or observing feasibility.
"""
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
INPUTS = HERE / 'inputs'


def main():
    path = INPUTS / 'jwst_sample_cycle6.csv'
    with path.open(newline='') as handle:
        rows = list(csv.DictReader(handle))
    families = Counter(row['family'] for row in rows)
    if len(rows) != 24 or families != {'fade': 12, 'rise': 12}:
        raise ValueError('Expected the agreed 24-target sample: 12 fade, 12 rise')
    for key in ('id', 'target', 'internal_id'):
        if len({row[key] for row in rows}) != len(rows):
            raise ValueError(f'Duplicate {key} in the working sample')
    completed = [row['id'] for row in rows if row['ngps'] == '2026-09-23']
    scheduled = [row['id'] for row in rows if row['ngps'] == '2026-10-26/27']
    if len(completed) + len(scheduled) != len(rows):
        raise ValueError('Unrecognized NGPS status; review the observing summary')
    summary = dict(
        selection='Agreed working design: 24 AGN, 12 fading and 12 rising; individual source audit pending.',
        source='inputs/jwst_sample_cycle6.csv',
        source_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        instrument='MIRI/MRS', visits_per_target=1, total=len(rows),
        family_counts=dict(families), requested_category='Medium',
        medium_charged_time_h=dict(lower_exclusive=50, upper_inclusive=130),
        charged_time_h=None, ETC_APT_validated=False,
        data_availability=dict(
            basis='User confirmation on 2026-09-26; not a file-by-file completeness audit.',
            spherex_all_targets=True, ztf_wise_manifold_all_targets=True),
        ngps=dict(completed_ids=completed, completed_count=len(completed),
                  scheduled_ids=scheduled, scheduled_count=len(scheduled),
                  scheduled_dates=['2026-10-26', '2026-10-27'],
                  note='CSV date identifies the September observing night; individual exposure UTC dates can be September 24.'),
        targets=[dict(id=row['id'], target=row['target'], internal_id=row['internal_id'],
                      family=row['family'], reversal=row['reversal'] == 'True') for row in rows],
        history_classifications_validated=False, response_phase_coverage_validated=False,
        statistical_power_validated=False,
        timing_note='Existing CSV estimates are provisional. Every warm delay implies ten times the hot delay; that assumption requires validation.',
    )
    (INPUTS / 'selection_summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    print(f'Working sample: {len(rows)} targets, 12 fading + 12 rising; '
          f'NGPS {len(completed)} completed + {len(scheduled)} scheduled.')


if __name__ == '__main__':
    main()

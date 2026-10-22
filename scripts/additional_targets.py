"""Attach explicitly requested archive reviews without changing observing queues."""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def attach_additional_targets(payload):
    targets = {t['name']: t for t in payload['targets']}
    for path in sorted((ROOT/'data/additional_targets').glob('*/target_payload.json')):
        target = json.loads(path.read_text())
        archive_path = path.with_name('archive_payload.json')
        if archive_path.exists():
            archive = json.loads(archive_path.read_text())
            review = target.setdefault('archive_review', {})
            review['summary'] = archive.pop('archive_summary') + ' ' + review.get('summary', '')
            review['note'] = archive.pop('archive_note') + ' ' + review.get('note', '')
            review.setdefault('sources', []).extend(archive.pop('archive_sources'))
            target.update(archive)
        name = target['name']
        from october_additions import DEFERRED
        if name in DEFERRED:
            target['archive_review']['note'] = DEFERRED[name] + ' ' + target['archive_review'].get('note', '')
        if name in targets:
            raise ValueError(f'Additional target duplicates an existing identity: {name}')
        # Protect against adding a known candidate under a different alias.
        import astropy.units as u
        from astropy.coordinates import SkyCoord
        existing = list(targets.values())
        if existing:
            coords = SkyCoord([t['ra'] for t in existing]*u.deg, [t['dec'] for t in existing]*u.deg)
            sep = coords.separation(SkyCoord(target['ra']*u.deg, target['dec']*u.deg)).arcsec
            if sep.min() < 1:
                raise ValueError(f'Additional target {name} overlaps {existing[sep.argmin()]["name"]}')
        targets[name] = target
    payload['targets'] = list(targets.values())
    for night in payload['nights']:
        payload['nights'][night]['count'] = sum(any(w['night']==night and
            (w.get('ranges') or w.get('tier_ranges',{}).get('fallback')) for w in t['nights']) for t in targets.values())
    return payload

"""Attach explicitly requested archive reviews without changing observing queues."""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def attach_additional_targets(payload):
    targets = {t['name']: t for t in payload['targets']}
    for path in sorted((ROOT/'data/additional_targets').glob('*/target_payload.json')):
        target = json.loads(path.read_text())
        name = target['name']
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

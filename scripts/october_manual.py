"""Supplemental target settings and scheduling constraints."""
import json
from pathlib import Path

import pandas as pd

CONFIG = Path(__file__).resolve().parents[1]/'observing/oct26_27_2026/manual_visits.json'


def visits():
    return json.loads(CONFIG.read_text()) if CONFIG.exists() else []


def filter_edges(targets, slots, edges):
    reserved = visits()
    result = []
    for edge in edges:
        i, sid = edge[:2]
        slot = slots[sid]
        start = slot['start']
        end = start + pd.Timedelta(minutes=float(targets.iloc[i].visit_minutes))
        if any(slot['night'] == r['night'] and start < pd.Timestamp(r['end_pdt'])
               and end > pd.Timestamp(r['start_pdt']) for r in reserved):
            continue
        moves = [r['relocated_visit'] for r in reserved if r.get('relocated_visit', {}).get('name') == targets.index[i]]
        if moves and not all(slot['night'] == r['night'] and start == pd.Timestamp(r['start_pdt']) for r in moves):
            continue
        result.append(edge)
    return result


def pause_note(row):
    # Supplemental targets are ordinary upload rows, not pause instructions.
    return ''


def upload_rows(regular, night, pack, regular_row):
    """Merge supplemental targets into the chronological telescope sequence."""
    combined = [(r, False) for r in regular]
    combined += [(r, True) for r in visits() if r['night'] == night]
    combined.sort(key=lambda item: item[0]['start_pdt'])
    result = []
    for order, (r, supplemental) in enumerate(combined, 1):
        note = f"{order:02d} SCI {r['start_pdt'][11:16]}"
        if supplemental:
            assert 0 < r['seconds_each'] <= r['max_exposure_seconds']
            assert r['upload_name'].isascii() and r['upload_name'].isalnum()
            settings = dict(slit_arcsec=r['slit_width_arcsec'], binspat=r['binspat'],
                            binspect=r['binspect'], slitangle=r['pa_deg'])
            row = pack.ngps_row(dict(r, name=r['upload_name']),
                dict(seconds_each=r['seconds_each'], exposures=r['exposures'], airmass=1.8),
                settings, note, f"{r['start_pdt'][11:16]} PDT")
        else:
            row = regular_row(r, note)
        result.append(row)
    assert len({r['name'] for r in result}) == len(result)
    return result


def sequence_rows(regular, night, supplemental=None):
    """A unified worksheet in the same order as the telescope upload."""
    combined = [dict(r) for r in regular]
    for r in visits() if supplemental is None else supplemental:
        if r['night'] != night:continue
        row = {key:r.get(key, '') for key in regular[0]}
        row.update(name=r['upload_name'], slit_arcsec=r['slit_width_arcsec'],
                   airmass_limit=1.8, priority_reason=r['science_priority'])
        combined.append(row)
    combined.sort(key=lambda r:r['start_pdt'])
    for order,r in enumerate(combined,1):r['order']=order
    return combined

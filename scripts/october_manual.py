"""Reserved manual acquisitions whose exposure settings are not yet executable."""
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
    return ' '.join(f"PAUSE after this visit: manual {r['name']} block {r['start_pdt'][11:16]}-{r['end_pdt'][11:16]} PDT; exposure unset; see acquisition note."
                    for r in visits() if row['night'] == r['night'] and row['end_pdt'] == r['start_pdt'])

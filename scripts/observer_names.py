"""Public display names and matching downloads, retaining stable internal links."""
import csv
import io
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def attach_names(payload):
    with (ROOT/'jwst_proposal/inputs/jwst_sample_cycle6.csv').open() as f:
        names = {r['internal_id']:r['internal_id'] if r['internal_id'].startswith('J') else 'J'+r['internal_id']
                 for r in csv.DictReader(f)}
    payload['display_names'] = names
    for r in payload['all_targets']:
        r['display_name'] = names.get(r['name'], r['name'])
    aliases = dict(names)
    pattern = re.compile(r'(?<![\w])('+'|'.join(map(re.escape, sorted(aliases,key=len,reverse=True)))+r')(?![\w])')
    for filename, text in list(payload['files'].items()):
        if not filename.endswith('.csv'):continue
        rows = list(csv.reader(io.StringIO(text)))
        if not rows:continue
        if filename.endswith('_coordinates.csv'):
            for r in rows:r[0] = aliases.get(r[0],r[0])
        else:
            header = [h.upper() for h in rows[0]]
            for row in rows[1:]:
                for j,key in enumerate(header):
                    if key in ['NAME','REPLACES_NAME','TARGET','INTERNAL_ID']:
                        row[j] = aliases.get(row[j],row[j])
                    elif key in ['COMMENT','NOTE','FLAGS','PRIORITY_REASON','SCIENCE_QUESTION']:
                        row[j] = pattern.sub(lambda m:aliases[m[0]],row[j])
        out = io.StringIO()
        csv.writer(out,lineterminator='\n').writerows(rows)
        payload['files'][filename] = out.getvalue()
        for directory in ['october','targets']:
            dest = ROOT/'docs'/directory/filename
            if dest.exists():dest.write_text(out.getvalue())

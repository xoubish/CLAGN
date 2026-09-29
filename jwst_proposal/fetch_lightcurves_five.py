"""ZTF public release + ALeRCE alerts + AllWISE multi-epoch photometry for the five NGPS targets.
Writes inputs/lightcurves/<id>_ztf_irsa.csv, <id>_alerce.csv, <id>_allwise_mep.csv (same formats as
fetch_fig1_lightcurves.py, which handles P9694 = R01)."""
import io, json
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen, Request
import pandas as pd

HERE = Path(__file__).resolve().parent; OUT = HERE/'inputs/lightcurves'; OUT.mkdir(exist_ok=True)
IDS = ['F02', 'F04', 'F06', 'F07', 'R01']
def get(url, params=None):
    if params: url = url+'?'+urlencode(params)
    return urlopen(Request(url, headers={'User-Agent': 'clagn/1.0'}), timeout=300).read().decode()

def main():
    s = pd.read_csv(HERE/'inputs/jwst_sample_cycle6.csv'); five = s[s.id.isin(IDS)]
    for _, r in five.iterrows():
        rad = 3/3600
        txt = get('https://irsa.ipac.caltech.edu/cgi-bin/ZTF/nph_light_curves', dict(POS=f'CIRCLE {r.ra} {r.dec} {rad:.5f}', BANDNAME='g,r', FORMAT='csv'))
        (OUT/f'{r.id}_ztf_irsa.csv').write_text(txt); nz = len(pd.read_csv(io.StringIO(txt))) if 'oid' in txt[:200] else 0
        objs = json.loads(get('https://api.alerce.online/ztf/v1/objects/', dict(ra=r.ra, dec=r.dec, radius=3, page_size=20)))
        frames = []
        for it in objs['items']:
            rows = json.loads(get(f"https://api.alerce.online/ztf/v1/objects/{it['oid']}/detections"))
            f = pd.DataFrame(rows); f.insert(0, 'oid', it['oid']); frames.append(f)
        na = 0
        if frames:
            a = pd.concat(frames, ignore_index=True); a.to_csv(OUT/f'{r.id}_alerce.csv', index=False); na = len(a)
        q = ('SELECT mjd,w1mpro_ep,w1sigmpro_ep,w2mpro_ep,w2sigmpro_ep,qi_fact,cc_flags,saa_sep,moon_masked FROM allwise_p3as_mep WHERE '
             f"CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{r.ra},{r.dec},{rad:.6f}))=1")
        txt = get('https://irsa.ipac.caltech.edu/TAP/sync', dict(QUERY=q, RESPONSEFORMAT='csv'))
        (OUT/f'{r.id}_allwise_mep.csv').write_text(txt); nm = len(pd.read_csv(io.StringIO(txt))) if not txt.startswith('<?xml') else -1
        print(r.id, r.internal_id, 'ztf rows', nz, 'alerce dets', na, 'allwise mep', nm, flush=True)

if __name__ == '__main__':
    main()

"""
00_literature_clagn.py  --  literature CLAGN coordinates (the Sample A Turn-on/Turn-off sources), for flagging
pool objects that are already published changing-look AGNs.  Uses the Fornax sample_selection getters
(NED / VizieR / SIMBAD); each is wrapped so one failing service does not stop the rest.

Output: data/literature_clagn.csv  (ra, dec, ref)  -- deduplicated at 2 arcsec.
"""
import os, sys, warnings
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, 'code_src'))
warnings.filterwarnings('ignore')
from astropy.coordinates import SkyCoord
import astropy.units as u
import sample_selection as ss

getters = [ss.get_lamassa_sample, ss.get_macleod16_sample, ss.get_ruan_sample, ss.get_macleod19_sample, ss.get_sheng_sample,
           ss.get_green_sample, ss.get_lyu_sample, ss.get_lopeznavas_sample, ss.get_hon_sample, ss.get_yang_sample,
           ss.get_graham_sample]
coords, labels = [], []
for g in getters:
    n0 = len(coords)
    try:
        g(coords, labels, verbose=0)
        print(f'  {g.__name__:28s} +{len(coords)-n0}')
    except Exception as e:
        print(f'  {g.__name__:28s} FAILED: {type(e).__name__}: {str(e)[:100]}')
        del coords[n0:]; del labels[n0:]

sc = SkyCoord(coords)
df = pd.DataFrame({'ra': sc.ra.deg, 'dec': sc.dec.deg, 'ref': labels})
# dedupe at 2 arcsec, keep the first reference
keep = np.ones(len(df), bool)
idx, sep, _ = sc.match_to_catalog_sky(sc, nthneighbor=2)
for i in range(len(df)):
    if keep[i] and sep[i].arcsec < 2 and idx[i] > i:
        keep[idx[i]] = False
df = df[keep].reset_index(drop=True)
out = os.path.join(HERE, 'data', 'literature_clagn.csv')
df.to_csv(out, index=False)
print(f'wrote {out}: {len(df)} unique literature CLAGNs; by ref: {df.ref.value_counts().to_dict()}')

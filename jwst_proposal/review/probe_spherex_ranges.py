from pathlib import Path
import pandas as pd, requests, time
p=Path(__file__).resolve().parent.parent/'local_data/spherex/sample'
url=pd.read_csv(p/'F01_index.csv').access_url.iloc[0]
for name,u in [('parent',url),('cutout',url+'?center=119.89561,32.362025&size=100arcsec')]:
    start=time.time()
    with requests.get(u,headers={'Range':'bytes=0-2879'},stream=True,timeout=40) as r:
        data=next(r.iter_content(2880))
        print(name,r.status_code,dict(r.headers),len(data),'seconds',time.time()-start,flush=True)

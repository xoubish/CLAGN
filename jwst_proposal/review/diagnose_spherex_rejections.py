"""Inspect a few saved rejected R12 outcomes without changing extraction cuts."""
from pathlib import Path
import json
import gzip
import sys
import numpy as np
import pandas as pd
import requests
from astropy.io import fits
from astropy.wcs import WCS

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent))
from spherex_extract import BAD_MASK, SOURCE_BIT


def main():
    folder=HERE.parent/'local_data/spherex/rejection_checks'
    folder.mkdir(parents=True,exist_ok=True)
    row=pd.read_csv(HERE.parent/'inputs/jwst_sample_cycle6.csv').set_index('id').loc['R12']
    records=[json.loads(f.read_text()) for f in sorted((HERE.parent/'local_data/spherex/sample/R12').glob('*.json'))]
    chosen=[]
    for det in [1,4,6]:
        matches=[r for r in records if r['status']=='quality_rejected' and f'D{det}_' in r['url']]
        if matches: chosen.append(matches[0])
    results=[]
    for record in chosen:
        name=record['url'].split('/')[-1].split('?')[0]
        path=folder/name
        if not path.exists():
            response=requests.get(record['url'],timeout=(12,45));response.raise_for_status()
            raw=response.content
            path.write_bytes(gzip.decompress(raw) if raw[:2]==b'\x1f\x8b' else raw)
        with fits.open(path) as h:
            im=h['IMAGE'].data;var=h['VARIANCE'].data;flags=h['FLAGS'].data
            x,y=[float(v) for v in WCS(h['IMAGE'].header).all_world2pix(row.ra,row.dec,0)]
            if not (2 <= x <= im.shape[1]-3 and 2 <= y <= im.shape[0]-3):
                result=dict(file=name,url=record['url'],detector=int(h['IMAGE'].header['DETECTOR']),
                            x=x,y=y,shape=list(im.shape),reason='source outside permitted cutout margin')
                results.append(result);print(json.dumps(result),flush=True)
                continue
            yy,xx=np.mgrid[:im.shape[0],:im.shape[1]];rr=np.hypot(xx-x,yy-y)
            good=np.isfinite(im)&np.isfinite(var)&(var>0)&((flags&BAD_MASK)==0)
            ann=(rr>2.5)&(rr<7)
            bkg=good&ann&((flags&SOURCE_BIT)==0)
            cy,cx=int(round(y)),int(round(x))
            center_flags=int(flags[cy,cx])
            result=dict(file=name,url=record['url'],detector=int(h['IMAGE'].header['DETECTOR']),
                        x=x,y=y,shape=list(im.shape),center_good=bool(good[cy,cx]),
                        center_flag_bits=[i for i in range(32) if center_flags&(1<<i)],
                        annulus_pixels=int(ann.sum()),good_annulus_pixels=int((good&ann).sum()),
                        background_pixels=int(bkg.sum()),required_background_pixels=15,
                        annulus_source_flagged=int((ann&((flags&SOURCE_BIT)!=0)).sum()),
                        center_flux_finite=bool(np.isfinite(im[cy,cx])))
            results.append(result);print(json.dumps(result),flush=True)
    report=HERE/'spherex_sample_review/rejection_checks/diagnostics.json'
    report.parent.mkdir(parents=True,exist_ok=True)
    report.write_text(json.dumps(results,indent=2)+'\n')


if __name__=='__main__':main()

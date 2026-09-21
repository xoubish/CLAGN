"""Refresh SDSS-V epochs from compact rolling-master per-field summaries.

Downloads every current northern-accessible field summary, including bad-quality
epochs for known targets; quality/position flags are retained rather than hiding
potentially interesting dates. All products are private research data.
"""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor,as_completed
import importlib,json,threading,gzip,io,time
import numpy as np
import pandas as pd
import requests
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u

ROOT=Path(__file__).resolve().parent;OUT=ROOT/'data/reselection_2026-09-20'
CACHE=OUT/'sdssv_rolling_fields';CACHE.mkdir(exist_ok=True)
ARCH=importlib.import_module('03f_sdssv_internal');LOCAL=threading.local()


def one(row):
    field,mjd=int(row['FIELD']),int(row['MJD']);path=CACHE/f'spAll-{field:06d}-{mjd}.fits.gz'
    url=f'{ARCH.REDUX}/master/spectra/daily/full/{field//1000:03d}XXX/{field:06d}/{mjd}/{path.name}'
    if not hasattr(LOCAL,'session'):
        LOCAL.session=requests.Session();LOCAL.session.auth=ARCH.auth()
    for attempt in range(3):
        try:
            if not path.exists():
                response=LOCAL.session.get(url,timeout=60);response.raise_for_status();raw=response.content
                if not raw.startswith(b'\x1f\x8b'):raise ValueError('Not a gzip spectrum summary')
                with fits.open(io.BytesIO(gzip.decompress(raw))) as h:d=ARCH.spall_table(h[1].data)
                temp=path.with_suffix('.tmp');temp.write_bytes(raw);temp.replace(path)
            else:
                with fits.open(path) as h:d=ARCH.spall_table(h[1].data)
            d['archive_version']='master';d['source']='SDSS-V rolling master 2026-09-20';d['proprietary']=True
            return d,dict(field=field,mjd=mjd,status='available',rows=len(d))
        except Exception as exc:
            if attempt==2:return pd.DataFrame(),dict(field=field,mjd=mjd,status='query failed',reason=type(exc).__name__)
            time.sleep(attempt+1)


def main():
    with fits.open(OUT/'fieldlist-master-20260920.fits') as h:
        d=h[1].data
        # An X<=1.5 target at Palomar has dec >= approximately -14.8 degrees.
        # The -17-degree field-center threshold includes the full fiber patrol
        # footprint plus margin. No RA or target-state restriction is applied.
        rows=[dict(FIELD=int(r['FIELD']),MJD=int(r['MJD'])) for r in d if r['DECCEN']>-17 and str(r['STATUS1D']).strip().lower()=='done']
    results=[];pieces=[]
    with ThreadPoolExecutor(max_workers=6) as ex:
        futures=[ex.submit(one,row) for row in rows]
        for i,f in enumerate(as_completed(futures),1):
            data,status=f.result();results.append(status)
            if len(data):pieces.append(data)
            if i%20==0 or i==len(futures):
                (OUT/'sdssv_rolling_field_status.json').write_text(json.dumps(dict(requested=len(rows),attempted=i,results=results),indent=2));print('Rolling fields',i,'/',len(rows),flush=True)
    if not pieces:return
    fresh=pd.concat(pieces,ignore_index=True);fresh=fresh[np.isfinite(fresh.ra)&np.isfinite(fresh.dec)&fresh.dec.between(-90,90)&fresh.catalogid.gt(0)].copy()
    fresh['key']=np.where(fresh.sdss_id>0,fresh.sdss_id,-fresh.catalogid)
    fresh['metadata_quality_ok']=fresh.zwarning.eq(0)&fresh.sn_median_all.ge(5)&fresh.fieldquality.eq('good')
    aliases=pd.read_parquet(OUT/'parent_aliases.parquet');mapping=aliases.drop_duplicates('name').set_index('name').canonical_name
    fresh['canonical_name']=('V'+fresh.key.astype(str)).map(mapping)
    parent=pd.read_parquet(OUT/'parent_with_manifold_descriptors.parquet')
    coords=SkyCoord(parent.ra.to_numpy()*u.deg,parent.dec.to_numpy()*u.deg)
    missing=fresh.canonical_name.isna()
    if missing.any():
        idx,sep,_=SkyCoord(fresh.loc[missing,'ra'].to_numpy()*u.deg,fresh.loc[missing,'dec'].to_numpy()*u.deg).match_to_catalog_sky(coords)
        indices=fresh.index[missing];valid=sep.arcsec<=2
        fresh.loc[indices[valid],'canonical_name']=parent.name.to_numpy()[idx[valid]]
    fresh.to_parquet(OUT/'sdssv_rolling_field_epochs.parquet',index=False)
    # New identities remain explicit; they are never described as searched on
    # the manifold until their W1 data and projection have actually been made.
    unmatched=fresh[fresh.canonical_name.isna()&fresh.metadata_quality_ok&fresh.cls.eq('QSO')&fresh.z.between(.001,1.14)].copy()
    unmatched=unmatched.sort_values('mjd').drop_duplicates('key',keep='last')
    unmatched['synthetic_r']=22.5-2.5*np.log10(unmatched.spectroflux_r.where(unmatched.spectroflux_r>0))
    unmatched.to_csv(OUT/'new_rolling_agn_identities_pending.csv',index=False)
    matched=fresh[fresh.canonical_name.notna()].copy()
    old=pd.read_parquet(OUT/'sdssv_epochs_matched.parquet');old['archive_version']=old.get('archive_version',old.run2d_row)
    combined=pd.concat([old,matched],ignore_index=True)
    combined=combined.drop_duplicates(['canonical_name','field','mjd','catalogid','archive_version'],keep='last')
    positions=parent.set_index('name').loc[combined.canonical_name]
    offsets=SkyCoord(combined.ra.to_numpy()*u.deg,combined.dec.to_numpy()*u.deg).separation(SkyCoord(positions.ra.to_numpy()*u.deg,positions.dec.to_numpy()*u.deg)).arcsec
    combined['position_offset_arcsec']=offsets;combined['metadata_quality_ok']&=offsets<=1
    combined.to_parquet(OUT/'sdssv_epochs_matched.parquet',index=False)
    print('Matched rolling epochs',len(matched),'new AGN identities pending W1',len(unmatched),'combined daily dates',combined[['canonical_name','mjd']].drop_duplicates().shape[0],flush=True)


if __name__=='__main__':main()

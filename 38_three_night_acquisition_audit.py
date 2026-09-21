"""Report actual per-target acquisition coverage, including empty/failed queries."""
from pathlib import Path
from datetime import datetime,timezone
import argparse,importlib,json
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parent;OUT=ROOT/'data/reselection_2026-09-20'


def read(path,default):
    return json.loads(path.read_text()) if path.exists() else default


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--require-baselines',action='store_true');args=parser.parse_args()
    targets=pd.read_csv(OUT/'compact_review_objects.csv')
    web=importlib.import_module('16_candidate_webpage');old=importlib.import_module('07_make_webpage')
    exposure=importlib.import_module('22_september_etc')
    wise=old.load_wise(str(ROOT/'data/wise_cache/pool'))
    zcoords=pd.read_csv(ROOT/'data/zeltyn_coords.csv')
    wise.update(old.load_wise(str(ROOT/'data/wise_cache/zeltyn'),{i+1:n for i,n in enumerate(zcoords.name)}))
    expansion=pd.read_parquet(OUT/'three_night_w1.parquet')
    counts=expansion.groupby('name').size().to_dict()
    visits={}
    for path in sorted((ROOT/'data').glob('neowise_visits_*.csv'))+sorted(OUT.glob('neowise_visits_prepared*.csv')):
        try:
            df=pd.read_csv(path)
            for n,g in df.groupby('name'):visits[n]=int(len(g))
        except pd.errors.EmptyDataError:continue
    desi_queries={}
    for path in (OUT/'three_night_desi').glob('*.json'):
        d=read(path,{})
        if isinstance(d,dict) and 'targets' in d:
            for n in d['targets']:
                if d.get('status')=='available' or n not in desi_queries:desi_queries[n]=d.get('status')
    inventory_path=ROOT/'data/spectra_epochs_three_night_desi.csv'
    inventory=pd.read_csv(inventory_path,dtype={'targetid':str}) if inventory_path.exists() else pd.DataFrame(columns=['name'])
    rows=[]
    for target in targets.itertuples():
        n=target.name
        public=read(ROOT/'data/spectra_dl'/f'{n}.json',[])
        private=read(OUT/'sdssv_spectra'/f'{n}.json',[])
        spectra=web.spectrum_payload(public+private)
        epochs=(spectra or {}).get('epochs',[])
        usable=[]
        for r in public+private:
            if r.get('coadd') and r.get('source')!='DESI':continue
            meta=r.get('meta',{});warning=meta.get('zwarning');sn=r.get('sn_median_all',meta.get('sn_median_all'))
            if r.get('metadata_quality_ok') is False or (warning is not None and warning!=0):continue
            if sn is not None and (not np.isfinite(sn) or sn<5):continue
            if sum(v is not None and np.isfinite(v) for v in r.get('flux',[]))>=20:usable.append(r)
        query=read(OUT/'three_night_archive_queries'/f'{n}_dr20_sdss.json',{})
        qfile=OUT/'three_night_archive_queries'/f'{n}_dr20_sdss.csv'
        q=pd.read_csv(qfile) if qfile.exists() else pd.DataFrame(columns=['sas_url'])
        cached={r.get('url') for r in public if r.get('cache_identity_version')==2}
        missing=[url for url in q.sas_url if isinstance(url,str) and url.startswith('https://') and url not in cached]
        ztf=read(ROOT/'data/ztf_cache/review_dr24_20260920'/f'{n}.json',{})
        image=read(ROOT/'data/cutouts'/f'{n}_sdss_wide.json',{})
        desi=inventory[inventory.name.eq(n)]
        expected={(str(r.targetid),str(r.survey),str(r.program)) for r in desi.itertuples()}
        loaded={(str(r.get('meta',{}).get('specid')),str(r.get('survey')),str(r.get('program'))) for r in public if r.get('source')=='DESI' and r.get('date_verified')}
        rows.append(dict(name=n,prepared_for_nights=target.prepared_for_nights,pool_role=target.pool_role,
            sdss_dr20_query=query.get('status','not queried'),sdss_dr20_identified_files=len(q),sdss_dr20_missing_files=len(missing),
            desi_query=desi_queries.get(n,'not queried'),desi_identified_coadds=len(expected),desi_missing_coadds=len(expected-loaded),
            sdssv_cached_dates=len({int(r['mjd']) for r in private}),plotted_traces=len(epochs),plotted_dates=len({r['epoch_day'] for r in epochs}),
            accepted_prior_files=len(usable),hbeta_continuum_baseline=exposure.reference(target._asdict()) is not None,w1_samples=counts.get(n,len(wise.get(n,{}).get('W1',[]))),neowise_visits=visits.get(n,0),
            ztf_status=ztf.get('status','not queried'),ztf_points=ztf.get('n_kept',0),ztf_last_mjd=ztf.get('last_mjd'),
            sdss_cutout_status=image.get('status','not retrieved')))
    audit=pd.DataFrame(rows);audit.to_csv(OUT/'prepared_acquisition_audit.csv',index=False)
    stats=dict(updated_utc=datetime.now(timezone.utc).isoformat(),targets=len(audit),targets_with_accepted_prior=int(audit.accepted_prior_files.gt(0).sum()),targets_with_hbeta_continuum_baseline=int(audit.hbeta_continuum_baseline.sum()),
        targets_with_w1=int(audit.w1_samples.gt(0).sum()),targets_with_neowise=int(audit.neowise_visits.gt(0).sum()),
        sdss_dr20_queries=audit.sdss_dr20_query.value_counts().to_dict(),sdss_dr20_missing_files=int(audit.sdss_dr20_missing_files.sum()),
        desi_queries=audit.desi_query.value_counts().to_dict(),desi_identified_coadds=int(audit.desi_identified_coadds.sum()),desi_missing_coadds=int(audit.desi_missing_coadds.sum()),
        sdssv_targets=int(audit.sdssv_cached_dates.gt(0).sum()),sdssv_dates=int(audit.sdssv_cached_dates.sum()),
        plotted_traces=int(audit.plotted_traces.sum()),ztf_status=audit.ztf_status.value_counts().to_dict(),
        sdss_cutouts=audit.sdss_cutout_status.value_counts().to_dict())
    (OUT/'prepared_acquisition_audit.json').write_text(json.dumps(stats,indent=2))
    print(json.dumps(stats,indent=2),flush=True)
    if args.require_baselines:
        assert audit.accepted_prior_files.gt(0).all(),'Some prepared objects lack an accepted loaded baseline; inspect audit before rebuilding the page.'
        assert audit.hbeta_continuum_baseline.all(),'Some prepared objects lack the Hbeta continuum baseline; reprepare the list after fetching files.'


if __name__=='__main__':main()

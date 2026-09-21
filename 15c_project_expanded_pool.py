"""Resume unWISE W1 acquisition and projection of newly feasible parent objects.

Downloads each sky partition once for a fixed queue, then applies the saved paper
model with its existing GP/normalization routine. --rmax is an acquisition batch,
not removal from the parent. Unknown-magnitude objects stay explicitly pending.
"""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor,as_completed
import argparse,json,os,pickle,sys,time,warnings
import numpy as np
import pandas as pd
import astropy.units as u
from astropy.coordinates import SkyCoord,search_around_sky
import hpgeom
import pyarrow as pa
import pyarrow.dataset as ds
import pyarrow.fs as fs

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'data/reselection_2026-09-20'
sys.path.insert(0,str(ROOT/'code_src'))
os.environ.setdefault('NUMBA_CACHE_DIR',str(OUT/'numba_cache'))
os.environ.setdefault('MPLCONFIGDIR',str(OUT/'matplotlib_cache'))
os.environ.setdefault('NUMBA_NUM_THREADS','8')


def fetch(rmax,workers,targets_path=None,batch_name=None,reuse_batch=None):
    todo=pd.read_csv(targets_path or OUT/'projection_todo.csv')
    todo=todo[todo.r_planning.le(rmax)].copy().reset_index(drop=True)
    batch=OUT/(batch_name or f'wise_r{str(rmax).replace(".","p")}')
    batch.mkdir(exist_ok=True)
    (batch/'status.json').unlink(missing_ok=True)
    if (batch/'targets.csv').exists():
        prior=pd.read_csv(batch/'targets.csv')
        if set(prior.name)!=set(todo.name):raise ValueError('Queue changed; use a new batch path to preserve cache provenance')
        todo=prior
    else:todo.to_csv(batch/'targets.csv',index=False)
    print('New objects in acquisition batch:',len(todo),flush=True)
    entries=[]
    for i,r in todo.iterrows():
        pixels=hpgeom.query_circle(32,r.ra,r.dec,1/3600,nest=True,lonlat=True,inclusive=True)
        entries.extend((int(p),i) for p in pixels)
    loc=pd.DataFrame(entries,columns=['pixel','idx'])
    s3=fs.S3FileSystem(region='us-west-2',anonymous=True,request_timeout=60,connect_timeout=15,
                       retry_strategy=fs.AwsStandardS3RetryStrategy(max_attempts=2))
    base='nasa-irsa-wise/unwise/neo7/catalogs/time_domain/healpix_k5/unwise-neo7-time_domain-healpix_k5.parquet'
    print('Sky partitions:',loc.pixel.nunique(),flush=True)
    if reuse_batch:
        previous=pd.read_csv(OUT/reuse_batch/'targets.csv').set_index('name')
        assert set(todo.name)<=set(previous.index),'A partition cache cannot cover targets absent from its original queue.'
        positions=previous.loc[todo.name,['ra','dec']].to_numpy()
        assert np.allclose(positions,todo[['ra','dec']].to_numpy(),atol=1e-9,rtol=0),'Cache coordinates changed.'
    # Exact linear unit conversion used in the original pipeline (W1 Vega->AB).
    factor=10**(-.4*(22.5+2.699-23.9))/1000
    def one(pixel,ind):
        path=batch/f'pixel_{pixel:05d}.parquet'
        if path.exists():
            try:return pixel,len(pd.read_parquet(path))
            except Exception:pass  # interrupted writes are fetched again
        sub=todo.iloc[ind]
        if reuse_batch:
            previous=OUT/reuse_batch/path.name
            if previous.exists():
                hit=pd.read_parquet(previous);hit=hit[hit.name.isin(sub.name)]
                hit.to_parquet(path,index=False);return pixel,len(hit)
        filt=(ds.field('primary')==1)&(ds.field('band')==1)
        err=None
        for attempt in range(3):
            try:
                # Inspect only this partition instead of downloading the
                # 140-MB all-sky metadata file before any useful work begins.
                dataset=ds.dataset(base+f'/healpix_k0={pixel//1024}/healpix_k5={pixel}',filesystem=s3,format='parquet')
                t=dataset.to_table(filter=filt,columns=['ra','dec','flux','dflux','MJDMEAN'],use_threads=False).to_pandas()
                coords=SkyCoord(t.ra.to_numpy()*u.deg,t.dec.to_numpy()*u.deg)
                targets=SkyCoord(sub.ra.to_numpy()*u.deg,sub.dec.to_numpy()*u.deg)
                ii,jj,sep,_=search_around_sky(targets,coords,1*u.arcsec)
                hit=t.iloc[jj].copy();hit['name']=sub.name.to_numpy()[ii];hit['match_arcsec']=sep.arcsec
                hit=hit[(hit.flux>0)&(hit.dflux>0)&np.isfinite(hit.flux)&np.isfinite(hit.dflux)]
                hit['flux']=hit.flux*factor;hit['err']=hit.dflux*factor
                hit=hit.rename(columns={'MJDMEAN':'time'})
                tmp=path.with_suffix('.tmp')
                hit[['name','time','flux','err','match_arcsec']].to_parquet(tmp,index=False)
                tmp.replace(path)
                return pixel,len(hit)
            except Exception as e:err=e;time.sleep(1+attempt)
        raise RuntimeError(f'pixel {pixel}: {err}')
    failed=[];start=time.time();done=0
    with ThreadPoolExecutor(max_workers=workers) as ex:
        tasks={ex.submit(one,int(p),g.idx.to_numpy()):int(p) for p,g in loc.groupby('pixel')}
        for f in as_completed(tasks):
            done+=1
            try:f.result()
            except Exception as e:failed.append({'pixel':tasks[f],'error':str(e)});print('FAILED',str(e),flush=True)
            if done%20==0 or done==len(tasks):print(f'{done}/{len(tasks)} partitions; {time.time()-start:.0f}s; failures {len(failed)}',flush=True)
    (batch/'status.json').write_text(json.dumps({'failed':failed,'completed':len(tasks)-len(failed),'targets':len(todo)},indent=2))
    if failed:raise SystemExit('Some partitions failed; rerun to resume.')
    print('Acquisition complete',flush=True)


def project(rmax,batch_name=None):
    from AGNzoo_functions import unify_lc_gp,stat_bands,combine_bands,normalize_clipmax_objects,dtw_distance
    batch=OUT/(batch_name or f'wise_r{str(rmax).replace(".","p")}')
    status=json.loads((batch/'status.json').read_text())
    if status['failed']:raise ValueError('Complete acquisition before treating batch as surveyed')
    targets=pd.read_csv(batch/'targets.csv')
    lc=pd.concat([pd.read_parquet(f) for f in sorted(batch.glob('pixel_*.parquet'))],ignore_index=True)
    # Exact duplicate matches can occur on pixel boundaries; no temporal coaddition.
    lc=lc.drop_duplicates(['name','time','flux','err']).sort_values(['name','time'])
    history=lc.groupby('name').agg(n_w1=('time','size'),w1_first_mjd=('time','min'),w1_last_mjd=('time','max'))
    history['w1_span_days']=history.w1_last_mjd-history.w1_first_mjd
    history['gp_extrapolated_fraction']=(1-history.w1_span_days/4000).clip(lower=0)
    lc['objectid']=pd.Categorical(lc.name,categories=targets.name).codes
    lc['label']='expanded';lc['band']='W1'
    lc=lc.set_index(['objectid','label','band','time'])[['flux','err']]
    with open(ROOT/'data/umap_w1_model.pkl','rb') as f:model=pickle.load(f)
    # The saved model uses exact all-reference DTW distances. Compute the identical
    # metric in parallel; keep UMAP's neighbor selection and transform unchanged.
    enable_parallel_exact_dtw()
    pieces=[];ids=np.array(lc.index.get_level_values('objectid').unique())
    for start in range(0,len(ids),500):
        out=batch/f'projection_{start:06d}.csv'
        if out.exists():pieces.append(pd.read_csv(out));continue
        sub=lc[lc.index.get_level_values('objectid').isin(ids[start:start+500])]
        objs,dobjs,labels,keeps=unify_lc_gp(sub,['W1'],xres=160,numplots=0,low_limit_size=5)
        if not len(objs):continue
        fvar,maxarr,meanarr=stat_bands(objs,dobjs,['W1'],sigmacl=5)
        dat=normalize_clipmax_objects(combine_bands(objs,['W1']),maxarr,band=0)
        valid=np.isfinite(dat).all(axis=1)
        chosen=np.array(sub.index.get_level_values('objectid').unique())[np.asarray(keeps)][valid]
        emb=model.transform(dat[valid])
        d=pd.DataFrame({'name':targets.name.to_numpy()[chosen],'umap_x':emb[:,0],'umap_y':emb[:,1],
                        'fvar_w1':fvar[0][valid],'mean_w1_mjy':meanarr[0][valid]})
        d=d.join(history,on='name')
        d.to_csv(out,index=False);pieces.append(d)
        print('Projected',min(start+500,len(ids)),'/',len(ids),flush=True)
    result=pd.concat(pieces,ignore_index=True) if pieces else pd.DataFrame(columns=['name','umap_x','umap_y'])
    result.to_csv(batch/'projected.csv',index=False)
    missing=targets[~targets.name.isin(result.name)].copy()
    missing['projection_status']=np.where(missing.name.isin(targets.name.to_numpy()[ids]),'insufficient or invalid W1 samples','no W1 match at 1 arcsec')
    missing.to_csv(batch/'not_projected.csv',index=False)
    print('Projected',len(result),'of',len(targets),'; unprojected',len(missing),flush=True)


def enable_parallel_exact_dtw():
    import numba
    import umap.umap_ as module
    from AGNzoo_functions import dtw_distance
    original=module.pairwise_distances
    @numba.njit(parallel=True)
    def matrix(x,y):
        out=np.empty((len(x),len(y)),dtype=np.float64)
        for i in numba.prange(len(x)):
            for j in range(len(y)):
                out[i,j]=dtw_distance(x[i],y[j])
        return out
    def parallel(x,y=None,metric=None,**kwargs):
        if y is not None and getattr(metric,'__name__','')=='dtw_distance':
            return matrix(np.asarray(x),np.asarray(y))
        return original(x,y,metric=metric,**kwargs)
    module.pairwise_distances=parallel
    return original


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('stage',choices=['fetch','project']);ap.add_argument('--rmax',type=float,default=19.5);ap.add_argument('--workers',type=int,default=6)
    ap.add_argument('--targets');ap.add_argument('--batch');ap.add_argument('--reuse-batch')
    args=ap.parse_args()
    if args.stage=='fetch':
        try:fetch(args.rmax,args.workers,args.targets,args.batch,args.reuse_batch)
        except Exception as exc:
            batch=OUT/(args.batch or f'wise_r{str(args.rmax).replace(".","p")}')
            batch.mkdir(parents=True,exist_ok=True)
            (batch/'status.json').write_text(json.dumps({'failed':[{'fatal_error':str(exc)}]},indent=2))
            raise
    else:project(args.rmax,args.batch)

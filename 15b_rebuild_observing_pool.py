"""Rebuild a broad manifold-first parent and independent observing feasibility.

Stages: sdssv, parent, geometry, export. Local outputs may contain collaboration
metadata and remain in the git-ignored data/reselection_2026-09-20 directory.
No old probability, optical trigger, old Moon weight or schedule enters selection.
"""
from pathlib import Path
import argparse, importlib, json, warnings
import numpy as np
import pandas as pd
import astropy.units as u
from astropy.io import fits
from astropy.coordinates import SkyCoord, AltAz, get_body, search_around_sky
from astropy.time import Time
from astropy.utils import iers
from scipy.spatial import cKDTree
from astroplan import moon_illumination

ROOT=Path(__file__).resolve().parent
DATA=ROOT/'data'
OUT=DATA/'reselection_2026-09-20'
OBS=importlib.import_module('05_observability')
iers.conf.auto_download=False
iers.conf.auto_max_age=None


def sdssv():
    """Reopen the complete cached SDSS-V catalog, independent of the old pool."""
    f=importlib.import_module('03f_sdssv_internal')
    path=DATA/'sdssv_internal/daily/spAll-lite-v6_2_1.fits.gz'
    def select(a):
        ra,dec=f.radec(a)
        carton=np.char.strip(a['FIRSTCARTON'])
        ot=np.char.strip(a['OBJTYPE'])
        bhm=(np.char.startswith(carton,b'bhm')|np.char.startswith(carton,b'open')|(ot==b'QSO')|(ot==b'GALAXY'))&~np.char.startswith(carton,b'mwm')
        science=~(np.char.startswith(ot,b'SKY')|np.char.startswith(ot,b'SPECTROPHOTO'))
        return bhm&science&np.isfinite(ra)&np.isfinite(dec)&(abs(dec)<=90)&(a['CATALOGID']>0)&(a['Z']>.001)&(a['Z']<1.14)
    d=f.stream_spall(str(path),select,convert=f.spall_table)
    d['source']='SDSS-V internal v6_2_1';d['proprietary']=True
    d['key']=np.where(d.sdss_id>0,d.sdss_id,-d.catalogid)
    d=d.drop_duplicates(['key','field','mjd','spec_file'])
    # Metadata quality only: this is not a broad-line-state classification.
    d['metadata_quality_ok']=(d.zwarning==0)&(d.sn_median_all>=5)&d.fieldquality.eq('good')
    d.to_parquet(OUT/'sdssv_epochs.parquet',index=False)
    print('SDSS-V reopened:',len(d),'epochs;',d.key.nunique(),'object keys; latest MJD',d.mjd.max(),flush=True)


def parent():
    frames=[]
    old=pd.read_csv(DATA/'parent_pool_scored.csv')
    old['name']='P'+old.poolid.astype(str);old['origin']='old_DR16_parent'
    old['reference_mjd']=old.mjd
    cols=['name','ra','dec','z','psfmag_r','reference_mjd','origin','umap_x','umap_y','clagn_score','in_region_clagn','in_region_zeltyn','poolid']
    frames.append(old[cols])
    q=pd.read_csv(OUT/'dr16_balmer_parent.csv',dtype={'specobjid':str,'objid':str})
    q['name']='S'+q.specobjid;q['origin']='expanded_DR16_'+q['class'];q['reference_mjd']=q.mjd
    frames.append(q[['name','ra','dec','z','psfmag_r','reference_mjd','origin']])
    # Supplement the query with DR16Q objects whose primary pipeline z/class differ.
    with fits.open(DATA/'external/dr16q_prop_Oct23_2022.fits',memmap=True) as h:
        t=h[1].data;z=np.asarray(t['Z_DR16Q'],float);m=(z>.001)&(z<1.14)
        cat=pd.DataFrame({'name':['Q'+str(x).strip() for x in t['SDSS_NAME'][m]],'ra':np.asarray(t['RA'][m],float),'dec':np.asarray(t['DEC'][m],float),
                          'z':z[m],'reference_mjd':np.asarray(t['MJD'][m],float),'origin':'full_DR16Q_catalog'})
        catalog_counts={'full_DR16Q_catalog':len(t),'DR16Q_Hbeta_coverage_branch':int(m.sum())}
    frames.append(cat)
    gal=pd.read_csv(DATA/'turnon_parent_galaxies.csv').rename(columns={'psfMag_r':'psfmag_r'})
    gal['origin']='old_galaxy_parent';gal['reference_mjd']=gal.mjd
    frames.append(gal[['name','ra','dec','z','psfmag_r','reference_mjd','origin']])
    z=pd.read_csv(DATA/'zeltyn_embedding_full.csv').rename(columns={'Name':'name'})
    z['origin']='Zeltyn24'
    frames.append(z[['name','ra','dec','z','origin','umap_x','umap_y','clagn_score','in_region_clagn','in_region_zeltyn']])
    if (OUT/'sampleA_coordinates.csv').exists():
        a=pd.read_csv(DATA/'sampleA_embedding_objectid.csv').merge(pd.read_csv(OUT/'sampleA_coordinates.csv'),on='objectid',validate='one_to_one')
        a['name']='A'+a.objectid.astype(str);a['origin']='paper_SampleA'
        frames.append(a)
    sv=pd.read_parquet(OUT/'sdssv_epochs.parquet')
    # Position medians suppress occasional mispositioned-fiber epochs for matching.
    s=sv.groupby('key').agg(ra=('ra','median'),dec=('dec','median'),z=('z','median')).reset_index()
    s['name']='V'+s.key.astype(str);s['origin']='SDSSV_internal_parent'
    frames.append(s[['name','ra','dec','z','origin']])
    raw=pd.concat(frames,ignore_index=True,sort=False)
    raw=raw[np.isfinite(raw.ra)&np.isfinite(raw.dec)&(abs(raw.dec)<=90)].reset_index(drop=True)
    c=SkyCoord(raw.ra.to_numpy()*u.deg,raw.dec.to_numpy()*u.deg)
    # Merge aliases by spherical position. Preserve all input rows for provenance.
    ii,jj,_,_=search_around_sky(c,c,2*u.arcsec)
    leader=np.arange(len(raw))
    def root(i):
        while leader[i]!=i:
            leader[i]=leader[leader[i]];i=leader[i]
        return i
    for i,j in zip(ii,jj):
        if i<j:
            a,b=root(i),root(j)
            if a!=b:leader[max(a,b)]=min(a,b)
    raw['group']=[root(i) for i in range(len(raw))]
    # Prefer existing projected identities, then stable parent order; no science score.
    raw['has_projection']=raw.umap_x.notna()&raw.umap_y.notna()
    reps=raw.sort_values('has_projection',ascending=False,kind='stable').drop_duplicates('group').set_index('group')
    canonical=reps.name.to_dict();raw['canonical_name']=raw.group.map(canonical)
    raw.to_parquet(OUT/'parent_aliases.parquet',index=False)
    p=reps.copy()
    # Fill missing photometry/redshift from aliases without dropping unknown objects.
    for col in ['psfmag_r','z']:
        good=raw[col].between(8,30) if col=='psfmag_r' else ((raw[col]>0)&np.isfinite(raw[col]))
        valid=raw[col].where(good)
        repgood=p[col].between(8,30) if col=='psfmag_r' else ((p[col]>0)&np.isfinite(p[col]))
        p[col]=p[col].where(repgood,valid.groupby(raw.group).first())
    p['aliases']=raw.groupby('group').name.agg(lambda x:';'.join(dict.fromkeys(x)))
    p['origins']=raw.groupby('group').origin.agg(lambda x:';'.join(sorted(set(x))))
    p=p.rename(columns={'psfmag_r':'r_catalog'})
    p['r_planning']=p.r_catalog;p['r_source']='historical SDSS catalog';p['r_epoch_mjd']=np.nan
    # Historical survey medians are more useful than an old image, but never 'current'.
    aliasmap=raw.drop_duplicates('name').set_index('name').canonical_name
    for filename in ['ztf_objects_pool.csv','ztf_objects_turnon.csv']:
        tab=pd.read_csv(DATA/filename)
        if 'name' not in tab:continue
        magcol=next((x for x in ['ztf_r_medianmag','r_medianmag'] if x in tab),None)
        if magcol:
            tab['canonical_name']=tab.name.map(aliasmap)
            vals=tab.dropna(subset=['canonical_name',magcol]).groupby('canonical_name')[magcol].median()
            v=p.name.map(vals);ok=v.between(8,30)
            p.loc[ok,'r_planning']=v[ok];p.loc[ok,'r_source']='historical ZTF catalog median'
    # Join newer spectra as dated metadata, not an automatic exclusion or score.
    sv['canonical_name']=('V'+sv.key.astype(str)).map(aliasmap)
    sv=sv.dropna(subset=['canonical_name'])
    sv.to_parquet(OUT/'sdssv_epochs_matched.parquet',index=False)
    for suffix,sub in [('any',sv),('quality',sv[sv.metadata_quality_ok])]:
        last=sub.sort_values(['mjd','sn_median_all']).drop_duplicates('canonical_name',keep='last').set_index('canonical_name')
        for col in ['mjd','sn_median_all','class','fieldquality']:
            if col not in last and col=='class':col='cls'
            p[f'sdssv_latest_{suffix}_{col}']=p.name.map(last[col])
        p[f'sdssv_n_{suffix}_epochs']=p.name.map(sub.groupby('canonical_name').mjd.nunique()).fillna(0).astype(int)
        flux=p.name.map(last.spectroflux_r)
        mag=22.5-2.5*np.log10(flux.where(flux>0))
        p[f'sdssv_latest_{suffix}_synthetic_r']=mag
        if suffix=='quality':
            fallback=p.r_planning.isna()&mag.between(8,30)
            p.loc[fallback,'r_planning']=mag[fallback]
            p.loc[fallback,'r_source']='SDSS-V synthetic aperture r; calibration unverified'
            p.loc[fallback,'r_epoch_mjd']=p.loc[fallback,'sdssv_latest_quality_mjd']
    p['Hb_A']=4861.33*(1+p.z);p['Ha_A']=6562.8*(1+p.z)
    p['Hb_in_NGPS']=p.Hb_A.between(3050,10400)
    p['Hb_line_edge_flag']=~p.Hb_A.between(3150,10300)
    p['outside_paper_z_domain']=p.z>1
    p['manifold_status']=np.where(p.has_projection,'projected','awaiting W1/projection')
    p=p.reset_index(drop=True).drop(columns=['group'],errors='ignore')
    p.to_parquet(OUT/'parent.parquet',index=False)
    catalog_counts.update({'input_catalog_rows':len(raw),'unique_coordinate_groups':len(p),'existing_projected':int(p.has_projection.sum()),
                           'original_paper_coordinates_loaded':(OUT/'sampleA_coordinates.csv').exists(),
                           'coordinate_match_radius_arcsec':2,'new_SDSS_query_rows':len(q)})
    (OUT/'parent_summary.json').write_text(json.dumps(catalog_counts,indent=2))
    print(json.dumps(catalog_counts,indent=2),flush=True)


def longest_run(mask,dt):
    cur=np.zeros(mask.shape[0]);best=cur.copy()
    for j,w in enumerate(dt):
        cur=np.where(mask[:,j],cur+w,0);best=np.maximum(best,cur)
    return best


def geometry():
    p=pd.read_parquet(OUT/'parent.parquet')
    nights=[];pieces=[]
    for night,(date,part) in OBS.NIGHTS.items():
        t0,t1,_,_=OBS.night_window(date,part)
        duration=(t1-t0).to_value(u.min)
        edges=np.r_[np.arange(0,duration,5),duration]
        times=t0+edges*u.min;dt=np.diff(edges)
        frame=AltAz(obstime=times,location=OBS.PALOMAR.location,pressure=0*u.hPa)
        moon=get_body('moon',times,OBS.PALOMAR.location).transform_to(frame)
        nights.append(dict(night=night,start_utc=t0.isot,end_utc=t1.isot,minutes=duration,
                           moon_illumination=float(moon_illumination(t0+(t1-t0)/2)),moon_max_alt=float(moon.alt.deg.max())))
        for start in range(0,len(p),2500):
            chunk=p.iloc[start:start+2500];c=SkyCoord(chunk.ra.to_numpy()*u.deg,chunk.dec.to_numpy()*u.deg)
            aa=c[:,None].transform_to(frame);x=aa.secz.value;sep=aa.separation(moon).deg
            valid=(aa.alt.deg>0)&(x<=2)&(x>0)
            d=pd.DataFrame({'name':chunk.name.values,'night':night})
            vx=np.where(valid,x,np.inf);best=vx.argmin(axis=1);row=np.arange(len(chunk))
            d['min_airmass']=np.where(valid.any(axis=1),vx.min(axis=1),np.nan)
            d['best_airmass_utc']=[times[k].isot if good else '' for k,good in zip(best,valid.any(axis=1))]
            d['moon_at_best_airmass']=np.where(valid.any(axis=1),sep[row,best],np.nan)
            d['moon_min_visible']=np.where(valid.any(axis=1),np.where(valid,sep,np.inf).min(axis=1),np.nan)
            d['moon_max_visible']=np.where(valid.any(axis=1),np.where(valid,sep,-np.inf).max(axis=1),np.nan)
            for label,am,ms in [('visible',2,0),('moon30',2,30),('moon40',2,40),('preferred',1.5,60),('usable',1.8,40)]:
                grid=valid&(x<=am)&(sep>=ms)
                # Accept interval only if both endpoints satisfy constraints.
                intervals=grid[:,:-1]&grid[:,1:]
                d[f'{label}_minutes']=(intervals*dt).sum(axis=1)
                d[f'{label}_longest_minutes']=longest_run(intervals,dt)
                first=np.where(intervals.any(axis=1),intervals.argmax(axis=1),-1)
                last=np.where(intervals.any(axis=1),intervals.shape[1]-np.flip(intervals,axis=1).argmax(axis=1),-1)
                if label in ['moon40','usable']:
                    d[f'{label}_first_utc']=[times[k].isot if k>=0 else '' for k in first]
                    d[f'{label}_last_utc']=[times[k].isot if k>=0 else '' for k in last]
            pieces.append(d)
            if start%25000==0:print(night,start,'/',len(p),flush=True)
    g=pd.concat(pieces,ignore_index=True)
    g.to_parquet(OUT/'night_geometry.parquet',index=False)
    pd.DataFrame(nights).to_csv(OUT/'night_windows.csv',index=False)
    print('Computed',len(g),'object-night geometries',flush=True)


def export():
    p=pd.read_parquet(OUT/'parent.parquet');g=pd.read_parquet(OUT/'night_geometry.parquet')
    # Apply completed expansion batches, without hiding unprojected targets.
    for file in sorted(OUT.glob('wise_r*/projected.csv')):
        new=pd.read_csv(file).set_index('name')
        match=p.name.isin(new.index)
        for col in ['umap_x','umap_y','fvar_w1','mean_w1_mjy','n_w1','w1_first_mjd','w1_last_mjd','w1_span_days','gp_extrapolated_fraction']:
            if col not in new:continue
            p.loc[match,col]=p.loc[match,'name'].map(new[col])
        p.loc[match,'has_projection']=True
        p.loc[match,'manifold_status']='projected in expanded pool'
        missing_file=file.with_name('not_projected.csv')
        if missing_file.exists():
            missing=pd.read_csv(missing_file).set_index('name')
            incomplete=p.name.isin(missing.index)&~p.has_projection
            p.loc[incomplete,'manifold_status']=p.loc[incomplete,'name'].map(missing.projection_status)
    # Manifold descriptors from the paper's labeled training sample, not empirical
    # SDSS-V outcome probabilities. Exclude a matched training self when available.
    a=pd.read_csv(DATA/'sampleA_embedding_objectid.csv')
    tree=cKDTree(a[['umap_x','umap_y']].values)
    use=p.has_projection
    # Re-evaluate region masks in the same saved embedding for new projections.
    _,xe,ye=np.histogram2d(a.umap_x,a.umap_y,bins=10)
    def bins(x,y):
        return np.clip(np.searchsorted(xe[1:],x),0,9),np.clip(np.searchsorted(ye[1:],y),0,9)
    ax,ay=bins(a.umap_x,a.umap_y);px,py=bins(p.loc[use,'umap_x'],p.loc[use,'umap_y'])
    for col in ['in_region_clagn','in_region_zeltyn']:
        mask=np.zeros((10,10),bool);sel=a[col].fillna(False).astype(bool).to_numpy();mask[ax[sel],ay[sel]]=True
        p.loc[use,col]=mask[px,py]
    dist,idx=tree.query(p.loc[use,['umap_x','umap_y']].values,k=51)
    selfmatch=dist[:,0]<1e-6
    idx=np.where(selfmatch[:,None],idx[:,1:],idx[:,:50])
    for label,bit in [('on',16),('off',32)]:
        flags=(a.label_bits.to_numpy().astype(int)&bit)>0
        p.loc[use,f'manifold_{label}_neighbor_fraction']=flags[idx].mean(axis=1)
    either=(a.label_bits.to_numpy().astype(int)&48)>0
    p.loc[use,'manifold_cl_neighbor_fraction']=either[idx].mean(axis=1)
    p['manifold_direction']=p.manifold_on_neighbor_fraction-p.manifold_off_neighbor_fraction
    p['manifold_interest']=np.where(use,np.where(p.in_region_clagn.fillna(False)|(p.manifold_cl_neighbor_fraction>=.10),'literature neighborhood',
                              np.where(p.in_region_zeltyn.fillna(False),'Zeltyn neighborhood','other projected')),'projection pending')
    m=g.merge(p,on='name',validate='many_to_one')
    m['geometry_class']=np.select([m.usable_longest_minutes>=45,m.moon40_longest_minutes>=30,m.moon30_longest_minutes>=20],
                                 ['good window','workable window','challenging window'],default='no adequate window')
    m['brightness_class']=np.select([m.r_planning.isna(),m.r_planning<=18.5,m.r_planning<=19.5,m.r_planning<=20.5],
                                   ['unknown','bright','moderate','faint'],default='very faint')
    m['spectral_weight_status']='pending dated broad-line-state review'
    # Preserve all objects, including faint and unknown brightness. These operational
    # subsets are review queues, not an approved science selection.
    m.to_parquet(OUT/'all_object_nights.parquet',index=False)
    visible=m[m.moon30_longest_minutes>=20].copy()
    visible.to_parquet(OUT/'observable_object_nights.parquet',index=False)
    for night in OBS.NIGHTS:
        d=visible[visible.night==night].sort_values(['has_projection','manifold_cl_neighbor_fraction','usable_longest_minutes'],ascending=False)
        d.to_csv(OUT/f'observable_{night}.csv',index=False)
    # W1 acquisition can be postponed until after deterministic feasibility to avoid
    # downloading curves for never-observable sources; the parent remains unchanged.
    feasible=visible[(visible.moon40_longest_minutes>=30)&(visible.r_planning.le(20.5)|visible.r_planning.isna())]
    todo=feasible[~feasible.has_projection].drop_duplicates('name')
    todo[['name','ra','dec','z','r_planning','brightness_class','origins']].to_csv(OUT/'projection_todo.csv',index=False)
    p.to_parquet(OUT/'parent_with_manifold_descriptors.parquet',index=False)
    summary=[]
    for night,d in m.groupby('night',sort=False):
        v=d[d.moon40_longest_minutes>=30]
        summary.append({'night':night,'unique_parent':len(d),'visible_20min':int((d.visible_longest_minutes>=20).sum()),
                        'moon30_20min':int((d.moon30_longest_minutes>=20).sum()),'moon40_30min':len(v),
                        'moon40_r_le19p5':int(v.r_planning.le(19.5).sum()),'moon40_unknown_r':int(v.r_planning.isna().sum()),
                        'projected_moon40_r_le19p5':int((v.has_projection&v.r_planning.le(19.5)).sum()),
                        'literature_neighborhood_moon40_r_le19p5':int((v.manifold_interest.eq('literature neighborhood')&v.r_planning.le(19.5)).sum())})
    pd.DataFrame(summary).to_csv(OUT/'selection_funnel.csv',index=False)
    print(pd.DataFrame(summary).to_string(index=False),flush=True)
    print('New feasible objects awaiting W1 projection:',len(todo),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('stage',choices=['sdssv','parent','geometry','export'])
    args=parser.parse_args();OUT.mkdir(parents=True,exist_ok=True)
    globals()[args.stage]()

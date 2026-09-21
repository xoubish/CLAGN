"""Build the provisional manifold candidate explorer, with explicit public fields.

Public output: docs/index.html and web/clagn_night_sheet.html.
Collaboration review: git-ignored reselection directory only. Never publish it.
Input tables must describe the same snapshot; this does not run the old scheduler.
"""
from pathlib import Path
from datetime import datetime, timezone
import base64
import importlib
import json
import warnings

import numpy as np
import pandas as pd
import astropy.units as u
from astropy.coordinates import SkyCoord, AltAz, get_body
from astropy.time import Time
from astropy.utils import iers
from astroplan import moon_illumination

ROOT = Path(__file__).resolve().parent
DATA = ROOT/'data'
OUT = DATA/'reselection_2026-09-20'
OLD = importlib.import_module('07_make_webpage')
OBS = importlib.import_module('05_observability')
iers.conf.auto_download = False
iers.conf.auto_max_age = None
SELECTION=json.loads((OUT/'review_selection.json').read_text()) if (OUT/'review_selection.json').exists() else dict(moon_min_deg=60,airmass_max=1.5,minimum_window_minutes=60,version='manifold-review-2026-09-20')


def native(value):
    if isinstance(value, dict):
        return {str(k): native(v) for k,v in value.items()}
    if isinstance(value, (list, tuple, np.ndarray)):
        return [native(v) for v in value]
    if isinstance(value, np.generic):
        value = value.item()
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value


def local(t):
    return pd.Timestamp(t.utc.datetime, tz='UTC').tz_convert('America/Los_Angeles').strftime('%b %d %H:%M PDT')


def spectrum_payload(records):
    """All available daily epochs, with date-only labels and no hidden epoch cap."""
    best={}
    for record in records:
        if record.get('coadd') and record.get('source','SDSS')!='DESI':
            continue  # All-epoch stacks reuse daily data and lack a single epoch.
        wave=record.get('wave',[]);flux=record.get('flux',[])
        if len(wave)!=len(flux) or sum(v is not None and np.isfinite(v) for v in flux)<20:
            continue
        mjd=record.get('mjd')
        dated=mjd is not None and np.isfinite(mjd) and mjd>40000
        day=int(np.floor(mjd)) if dated else None
        key=(record.get('source','SDSS'),day)
        # Prefer the explicit current reduction when a cached older daily record
        # represents the same source and night; never count it as a new epoch.
        def quality(r):
            sn=r.get('sn_median_all',r.get('meta',{}).get('sn_median_all'))
            return (r.get('run2d')=='v6_2_1',float(sn) if sn is not None and np.isfinite(sn) else -1)
        if key not in best or quality(record)>quality(best[key]):best[key]=record
    epochs=[]
    for key,r in sorted(best.items(),key=lambda p:(p[0][1] is None,p[0][1] or 0,p[0][0])):
        day=key[1]
        date=Time(day,format='mjd').strftime('%Y-%m-%d') if day is not None else 'Date unavailable'
        flagged=r.get('metadata_quality_ok') is False
        epochs.append(dict(label=date,date=date,mjd=r.get('mjd'),epoch_day=day,
            wave=r['wave'],flux=r['flux'],quality_flag=flagged,
            quality_note='Metadata quality checks failed; inspect before interpreting.' if flagged else '',
            coadd=bool(r.get('coadd',False))))
    return native(dict(epochs=epochs)) if epochs else None


def reconcile_spectral_dates(target):
    """Keep inventory dates distinct from available files and active plot traces."""
    plotted={e['epoch_day'] for e in (target.get('spec') or {}).get('epochs',[]) if e['epoch_day'] is not None}
    existing={int(e['mjd']) for e in target['epochs']}
    for day in sorted(plotted-existing):
        target['epochs'].append(dict(mjd=day,date=Time(day,format='mjd').strftime('%Y-%m-%d'),src='Archive'))
    target['epochs'].sort(key=lambda e:e['mjd'])
    for e in target['epochs']:e['file_available']=int(e['mjd']) in plotted
    target['n_spec']=len(target['epochs'])
    target['n_spec_available_dates']=len(plotted)


def review_ztf(name):
    """Use the refreshed curve, with explicit query/coverage status and dates."""
    tag='review_dr24_20260920'
    path=DATA/'ztf_cache'/tag/f'{name}.json'
    info=json.loads(path.read_text()) if path.exists() else {}
    if info.get('status') in ['available','no usable photometry']:
        series=OLD.ztf_series(name,tags=(tag,)) if info['status']=='available' else {}
        meta={key:info.get(key) for key in ['status','collection','queried_utc','n_g','n_r','association_warning']}
        first,last=info.get('first_mjd'),info.get('last_mjd')
    else:
        # An empty cache in one old directory must not hide a valid later cache.
        options=[(old,OLD.ztf_series(name,tags=(old,))) for old in ['pool','v2','calib','zeltyn']]
        old,series=max(options,key=lambda item:sum(len(v) for v in item[1].values()))
        first=last=None
        if series:
            cached=pd.read_csv(DATA/'ztf_cache'/old/f'{name}.csv')
            usable=cached[cached.filtercode.isin(['zg','zr']) & np.isfinite(cached.mag) & np.isfinite(cached.mjd)]
            first,last=float(usable.mjd.min()),float(usable.mjd.max())
        meta=dict(status='refresh failed; showing earlier cache' if info.get('status')=='query failed' and series
                  else 'query failed' if info.get('status')=='query failed'
                  else 'earlier cache' if series else 'not fetched',
                  collection=None,attempted_collection=info.get('collection'),queried_utc=info.get('queried_utc'))
    meta['first_date']=Time(first,format='mjd').strftime('%Y-%m-%d') if first is not None else None
    meta['last_date']=Time(last,format='mjd').strftime('%Y-%m-%d') if last is not None else None
    return series,meta


def observing_windows(targets):
    coords = SkyCoord(targets.ra.to_numpy()*u.deg,targets.dec.to_numpy()*u.deg)
    windows = {n: [] for n in targets.name}
    meta = {}
    expected = pd.read_csv(OUT/'compact_review_object_nights.csv')
    for key,(date,part) in OBS.NIGHTS.items():
        t0,t1,_,_ = OBS.night_window(date,part)
        duration = (t1-t0).to_value(u.min)
        edges = np.r_[np.arange(0,duration,5),duration]
        times = t0+edges*u.min
        frame = AltAz(obstime=times,location=OBS.PALOMAR.location,pressure=0*u.hPa)
        altaz = coords[:,None].transform_to(frame)
        x = altaz.secz.value
        sep = altaz.separation(get_body('moon',times,OBS.PALOMAR.location).transform_to(frame)).deg
        grid = (altaz.alt.deg>0)&(x>0)&(x<=SELECTION['airmass_max'])&(sep>=SELECTION['moon_min_deg'])
        accepted = grid[:,:-1]&grid[:,1:]
        for i,name in enumerate(targets.name):
            changes = np.diff(np.r_[False,accepted[i],False].astype(int))
            runs = list(zip(np.where(changes==1)[0],np.where(changes==-1)[0]))
            longest = max((edges[b]-edges[a] for a,b in runs),default=0)
            minimum=SELECTION.get('reserve_minimum_window_minutes',30) if targets.iloc[i].get('pool_role')=='reserve' else SELECTION['minimum_window_minutes']
            if longest < minimum:
                continue
            ref = expected[(expected.name==name)&(expected.night==key)]
            assert len(ref)==1 and abs(float(ref.iloc[0].preferred_longest_minutes)-longest)<.01
            # Display every qualifying continuous window, never bridge a gap.
            usable_runs = [(a,b) for a,b in runs if edges[b]-edges[a]>=minimum]
            indices = np.concatenate([np.arange(a,b+1) for a,b in usable_runs])
            windows[name].append(dict(night=key,longest_minutes=round(longest,1),
                ranges=[dict(start=local(times[a]),end=local(times[b]),start_utc=times[a].isot,end_utc=times[b].isot,
                             minutes=round(edges[b]-edges[a],1)) for a,b in usable_runs],
                min_airmass=round(float(x[i,indices].min()),3),moon_min=round(float(sep[i,indices].min()),1),
                minutes_airmass_le1p3=round(float(((grid[i,:-1]&grid[i,1:]&(x[i,:-1]<=1.3)&(x[i,1:]<=1.3))*np.diff(edges)).sum()),1),
                moon_max=round(float(sep[i,indices].max()),1),
                curve=[[round(float(t.mjd),5),round(float(xx),3) if 0<xx<4 else None,round(float(ss),2)]
                       for t,xx,ss in zip(times,x[i],sep[i])]))
        meta[key] = dict(label={'sep23':'Sep 23','oct26':'Oct 26','oct27':'Oct 27'}[key],date=date,
            part='first half' if part=='first' else 'full night',window=f'{local(t0)} → {local(t1)}',
            start_mjd=float(t0.mjd),end_mjd=float(t1.mjd),mjd=float((t0+(t1-t0)/2).mjd),
            moon_percent=round(100*float(moon_illumination(t0+(t1-t0)/2))),
            count=sum(any(w['night']==key for w in v) for v in windows.values()))
    assert sum(map(len,windows.values()))==len(expected)
    assert all(windows.values())
    return windows,meta


def main():
    targets = pd.read_csv(OUT/'compact_review_objects.csv')
    audit = pd.read_csv(OUT/'compact_spectral_audit.csv').set_index('name')
    assert targets.name.is_unique and set(targets.name)==set(audit.index), 'Refresh spectral audit for this compact pool.'
    # Current snapshot has exclusively public DR16 identities and public photometry.
    # Future additions need the same provenance checks before website publication.
    public_names=set(targets.loc[
        targets.origin.isin(['old_DR16_parent','Zeltyn24','expanded_DR16_QSO','expanded_DR16_GALAXY','full_DR16Q_catalog','old_galaxy_parent']) &
        targets.r_source.isin(['historical SDSS catalog','historical ZTF catalog median']),'name'])
    windows,nights = observing_windows(targets)
    print('Validated observing windows for',len(targets),'candidates',flush=True)
    inventory = pd.read_csv(OUT/'compact_spectral_inventory_records.csv')
    allowed = ['parent reference spectra','spectra_epochs_pool.csv','spectra_epochs_zeltyn.csv','spectra_epochs_v2pub.csv','spectra_epochs_three_night_public.csv','spectra_epochs_three_night_public_dr19.csv','spectra_epochs_three_night_public_dr20.csv','spectra_epochs_three_night_desi.csv']
    public = inventory[inventory.inventory.isin(allowed)].copy()
    public = public.groupby(['target_name','epoch_day'],as_index=False).agg(
        source=('survey_family',lambda v:'/'.join(sorted(set(v)))))
    full_epochs = pd.read_csv(OUT/'compact_spectral_epochs.csv')
    known = pd.read_csv(OUT/'compact_known_state_matches.csv').fillna('')
    unwise_path=OUT/'neighbour_unwise_manifest.csv'
    unwise_screen=pd.read_csv(unwise_path).set_index('name').to_dict('index') if unwise_path.exists() else {}
    image_status={}
    manifest=OUT/'image_manifest.csv'
    if manifest.exists():
        image_status=pd.read_csv(manifest).set_index('name').to_dict('index')
    wise = OLD.load_wise(str(DATA/'wise_cache/pool'))
    zcoords=pd.read_csv(DATA/'zeltyn_coords.csv')
    wise.update(OLD.load_wise(str(DATA/'wise_cache/zeltyn'),{i+1:n for i,n in enumerate(zcoords.name)}))
    expansion=OUT/'three_night_w1.parquet'
    if expansion.exists():
        for name,g in pd.read_parquet(expansion).groupby('name'):
            if name in set(targets.name):
                wise[name]={'W1':g.sort_values('time')[['time','flux','err']].round(4).values.tolist()}
    neo = {}
    for tag in ['pool','poolall','zeltyn','three_night']:
        path = DATA/f'neowise_visits_{tag}.csv'
        if path.exists():
            v = pd.read_csv(path)
            for name,g in v[v.name.isin(targets.name)].groupby('name'):
                g = g.sort_values('mjd')
                flux = 309.54e3*10**(-.4*g.w1.to_numpy())
                neo[name] = [[round(a,1),round(b,4),round(c,4)] for a,b,c in
                            zip(g.mjd,flux,flux*.921*g.w1err.fillna(.05))]
    a = pd.read_csv(DATA/'sampleA_embedding_objectid.csv')
    _,xe,ye = np.histogram2d(a.umap_x,a.umap_y,bins=10)
    ix = np.clip(np.searchsorted(xe[1:],a.umap_x),0,9)
    iy = np.clip(np.searchsorted(ye[1:],a.umap_y),0,9)
    rects=[]
    for col,kind in [('in_region_clagn','literature'),('in_region_zeltyn','Zeltyn')]:
        rects.extend(dict(kind=kind,x0=xe[i],x1=xe[i+1],y0=ye[j],y1=ye[j+1])
                     for i,j in sorted(set(zip(ix[a[col]],iy[a[col]]))))
    manifold=dict(background=a[['umap_x','umap_y']].round(3).values.tolist(),rects=rects,
                  xlim=[xe[0],xe[-1]],ylim=[ye[0],ye[-1]])
    items=[];all_records={}
    for r in targets.itertuples():
        name=r.name
        coord=SkyCoord(r.ra*u.deg,r.dec*u.deg)
        ep=public[public.target_name==name].sort_values('epoch_day')
        histories=[dict(mjd=int(v.epoch_day),date=Time(v.epoch_day,format='mjd').strftime('%Y-%m-%d'),src=v.source)
                   for v in ep.itertuples()]
        spec=None
        path=DATA/'spectra_dl'/f'{name}.json'
        records=[]
        if path.exists():
            records=json.loads(path.read_text())
            spec=spectrum_payload([v for v in records if not v.get('proprietary')])
        internal=OUT/'sdssv_spectra'/f'{name}.json'
        if internal.exists():records+=json.loads(internal.read_text())
        all_records[name]=records
        cut=None
        for kind in ['sdss_wide','sdss','ps1_r','ps1_g']:
            path=DATA/'cutouts'/f'{name}_{kind}.jpg'
            if path.exists() and path.stat().st_size>100:
                info={}
                metadata=path.with_suffix('.json')
                if metadata.exists():info=json.loads(metadata.read_text())
                cut=dict(source=info.get('source',{'sdss_wide':'SDSS archival image','sdss':'SDSS archival image','ps1_r':'Pan-STARRS r','ps1_g':'Pan-STARRS g'}[kind]),
                         field_arcsec=info.get('field_arcsec'),
                         image='data:image/jpeg;base64,'+base64.b64encode(path.read_bytes()).decode())
                break
        matches=known[known.name==name]
        ztf,ztf_status=review_ztf(name)
        screen_path=OUT/'neighbour_cache'/f'{name}_screen.json'
        screen=json.loads(screen_path.read_text()) if screen_path.exists() else dict(status='not checked',flags=[],neighbours=[])
        uw=unwise_screen.get(name,dict(status='not checked'))
        screen['unwise']=uw
        complete=screen.get('sdss_query')=='available' and screen.get('gaia_query')=='available' and uw.get('status')=='matched'
        wise_flag=uw.get('status')=='matched' and (uw.get('fracflux_w1',0)<.8 or uw.get('flags_unwise_w1',0)!=0)
        if wise_flag:screen['flags'].append('unWISE deblending/quality needs review; a host contribution is also possible')
        field_status='review' if screen['flags'] else 'clear' if complete else 'pending'
        screen['neighbours']=sorted(screen['neighbours'],key=lambda n:(not bool(n.get('review_reasons')),n['sep_arcsec']))[:20]
        status=audit.loc[name,'known_state_status']
        status_key={'catalog-confirmed CLAGN':'confirmed','reported candidate':'candidate',
                    'literature match; confirmation needs audit':'unverified','no match in checked catalogs':'unmatched'}[status]
        items.append(dict(name=name,jname='J'+coord.ra.to_string(u.hour,sep='',precision=2,pad=True)+coord.dec.to_string(u.deg,sep='',precision=1,alwayssign=True,pad=True),
            ra=round(r.ra,6),dec=round(r.dec,6),sex=coord.to_string('hmsdms',sep=':',precision=1),z=r.z,
            rmag=round(r.r_planning,2),r_source=r.r_source,region=r.review_region,
            ux=round(r.umap_x,4),uy=round(r.umap_y,4),on_fraction=r.manifold_on_neighbor_fraction,off_fraction=r.manifold_off_neighbor_fraction,
            status=status_key,known=matches[['catalog','catalog_name','status','transition','reference']].to_dict('records'),
            n_spec=len(histories),epochs=histories,nights=windows[name],ztf=ztf,ztf_status=ztf_status,wise=wise.get(name,{}),neo=neo.get(name,[]),
            pool_role=getattr(r,'pool_role','manifold'),field_status=field_status,neighbour_screen=screen,
            spec=spec,cut=cut,image_status=image_status.get(name,{}).get('status','not fetched'),
            lines=[dict(name=n,angstrom=round(w*(1+r.z),1),inrange=bool(3050<=w*(1+r.z)<=10400))
                                 for n,w in [('Hβ',4861.33),('[O III]',5006.84),('Hα',6562.8)] ]))
        reconcile_spectral_dates(items[-1])
    payload=native(dict(version=SELECTION['version'],selection=SELECTION,generated=datetime.now(timezone.utc).isoformat(),
                       access='public',nights=nights,manifold=manifold,targets=items))
    # Reuse the existing calibrated display units and light-curve/spectrum renderers.
    charts=OLD.TEMPLATE[OLD.TEMPLATE.index('function mjdToYear'):OLD.TEMPLATE.index('/* ---------- manifold thumbnail')]
    charts+='\n'+(ROOT/'web/candidate_spectra.js').read_text()
    template=(ROOT/'web/candidate_review_template.html').read_text()
    def html_for(value):
        encoded=json.dumps(native(value),separators=(',',':'),allow_nan=False).replace('<','\\u003c')
        return template.replace('__CHART_FUNCTIONS__',charts).replace('__PAYLOAD__',encoded)
    # Objects selected using internal-only identities or brightness remain local.
    public_payload=json.loads(json.dumps(payload))
    public_payload['targets']=[t for t in public_payload['targets'] if t['name'] in public_names]
    for night in public_payload['nights']:
        public_payload['nights'][night]['count']=sum(any(w['night']==night for w in t['nights']) for t in public_payload['targets'])
    page=html_for(public_payload)
    assert 'proprietary' not in page and 'SDSS-V internal' not in page
    for path in [ROOT/'docs/index.html',ROOT/'web/clagn_night_sheet.html']:
        path.write_text(page)
    # Complete metadata remains in the ignored local research directory.
    private=json.loads(json.dumps(payload))
    private['access']='collaboration'
    for target in private['targets']:
        e=full_epochs[full_epochs.target_name==target['name']].sort_values('epoch_day')
        target['epochs']=[dict(mjd=int(v.epoch_day),date=Time(v.epoch_day,format='mjd').strftime('%Y-%m-%d'),src=v.surveys.replace(';','/')) for v in e.itertuples()]
        target['n_spec']=len(e)
        target['sdssv_dates']=int(audit.loc[target['name'],'n_sdssv_dates'])
        target['sdssv_quality_dates']=int(audit.loc[target['name'],'n_sdssv_metadata_quality_dates'])
        target['spec']=spectrum_payload(all_records[target['name']])
        reconcile_spectral_dates(target)
    (OUT/'candidate_review_local.html').write_text(html_for(private))
    counts=dict(objects=len(items),public_multiple=sum(t['n_spec']>=2 for t in items),
                public_spectra_plotted=sum(bool(t['spec']) for t in items),ztf_curves=sum(bool(t['ztf']) for t in items),
                wise_curves=sum(bool(t['wise']) for t in items),cutouts=sum(bool(t['cut']) for t in items),
                public_bytes=len(page.encode()),night_counts={n:v['count'] for n,v in nights.items()},
                public_spectral_traces=sum(len((t['spec'] or {}).get('epochs',[])) for t in items),
                local_spectral_traces=sum(len((t['spec'] or {}).get('epochs',[])) for t in private['targets']),
                local_targets_with_spectra=sum(bool(t['spec']) for t in private['targets']))
    (OUT/'web_build_summary.json').write_text(json.dumps(counts,indent=2))
    print(json.dumps(counts,indent=2),flush=True)


if __name__=='__main__':
    main()

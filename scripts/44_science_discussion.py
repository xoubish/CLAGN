"""Build a local science comparison, provisional matches, and diagnostic plots.

This is a discussion document, not an observing queue or a CLAGN classifier.
No private information is exported to docs/ or web/.
"""
from pathlib import Path
import base64, html, importlib, json
from functools import lru_cache
import numpy as np
import pandas as pd
from astropy.time import Time
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1];DEST=ROOT/'observing/science_review'
REVIEW=importlib.import_module('42_science_reselection')
NOW=Time('2026-09-23').mjd
VISUAL_NOTES={
    'P4675':'Broad Hβ and Hα are visible in the 2022 reference, followed by optical fading. A direct weakening/disappearance test; neighbor and slit review is essential.',
    'P3456':'Broad Balmer lines remain visible in 2021. Brightening since that spectrum would test recovery/strengthening, not establish a turn-on from a demonstrated off state. Complex field.',
    'P12144':'Broad Balmer wings are visible in both plotted epochs. Test weakening after the 2021 reference; continuum normalization and aperture differences prevent reading absolute line changes from this plot.',
    'P11113':'Broad Hβ and Hα are visible in the 2018 reference. The 2026 r-band alert brightening makes this a recovery/strengthening candidate; no intervening off state is established.',
    'P12457':'The plotted Balmer profiles are predominantly narrow; weak broad emission cannot be excluded without host subtraction and native errors. A blind turn-on test, with no established dated photometric trigger.',
    'P10596':'Broad Balmer profiles persist through the 2023 reference. Recent alert detections suggest brightening but lack a robust within-alert 2025 comparison.',
    'P11530':'Broad Hβ remains visible through December 2025. Strong changes in the normalized Hα region need aperture/calibration checks; multiple spectra do not themselves imply repeated state changes.',
    'P8086':'Broad Hβ is visible in the 2021 reference. Catalog fading is followed by tentative 2026 recovery in both alert bands. Hα falls outside the planned NGPS coverage.',
    'P18343':'The recovered February 2026 spectrum still shows broad Balmer emission. This later reference supersedes the old 2022 baseline; current evidence does not establish a subsequent change.',
    'P20498':'The recovered March 2025 spectrum shows broad Balmer emission. Two later g-band alert nights are insufficient to establish the current state.'
}


@lru_cache(maxsize=1)
def wise_visits():
    frames=[]
    for path in list((ROOT/'data').glob('neowise_visits_*.csv'))+list(REVIEW.OUT.glob('neowise_visits_prepared*.csv')):
        try:
            d=pd.read_csv(path)
            if {'name','mjd','w1'}.issubset(d):frames.append(d)
        except pd.errors.EmptyDataError:pass
    d=pd.concat(frames,ignore_index=True).drop_duplicates(['name','mjd'],keep='last')
    return dict(tuple(d.groupby('name')))


def controls(e):
    """Propose pairs with explicit calipers; spectral-state matching remains open."""
    a=e[e.experiment_group.eq('A: manifold + dated change')&e.accepted_dates.gt(0)&~e.optical_association_warning]
    c=e[e.experiment_group.eq('C: variability comparison')&e.accepted_dates.gt(0)&~e.optical_association_warning]
    candidates=[]
    for x in a.itertuples():
        if x.known_state_status!='no match in checked catalogs':continue
        for y in c.itertuples():
            if y.known_state_status!='no match in checked catalogs':continue
            nights=sorted(set(x.eligible_nights.split(','))&set(y.eligible_nights.split(',')))
            if not nights or x.manifold_cl_neighbor_fraction-y.manifold_cl_neighbor_fraction<.10-1e-7:continue
            if x.post_spectrum_optical_trigger!=y.post_spectrum_optical_trigger:continue
            if x.post_spectrum_optical_trigger:
                if x.trigger_direction!=y.trigger_direction:continue
                amp1,amp2=x.optical_change_abs,y.optical_change_abs
            else:
                amp1,amp2=abs(x.neowise_change_after_reference_mag),abs(y.neowise_change_after_reference_mag)
                if np.sign(x.neowise_change_after_reference_mag)!=np.sign(y.neowise_change_after_reference_mag):continue
            dz=abs(x.z-y.z);dm=abs(x.r_planning-y.r_planning);da=abs(amp1-amp2)
            age1,age2=(NOW-x.reference_mjd)/365.25,(NOW-y.reference_mjd)/365.25
            if min(age1,age2)<=0:continue
            ratio=max(age1,age2)/min(age1,age2)
            if dz>.15 or dm>.7 or da>.3 or ratio>2.5:continue
            cost=dz/.15+dm/.7+da/.3+np.log(ratio)/np.log(2.5)
            candidates.append(dict(manifold_target=x.name,comparison_target=y.name,eligible_nights=','.join(nights),
                                   delta_z=dz,delta_r=dm,delta_variability=da,baseline_age_ratio=ratio,
                                   manifold_neighbor_fraction=x.manifold_cl_neighbor_fraction,
                                   comparison_neighbor_fraction=y.manifold_cl_neighbor_fraction,cost=cost,
                                   status='Provisional covariate match; broad-line state, cadence, host fraction and field still require review'))
    result=pd.DataFrame(candidates,columns=['manifold_target','comparison_target','eligible_nights','delta_z','delta_r','delta_variability','baseline_age_ratio','manifold_neighbor_fraction','comparison_neighbor_fraction','cost','status'])
    result.sort_values('cost').to_csv(DEST/'possible_matched_controls.csv',index=False)
    return result


def diagnostic(t,epochs):
    name=t['name'];path=DEST/'figures'/f'{name}.png';path.parent.mkdir(exist_ok=True)
    fig,grid=plt.subplots(2,2,figsize=(14,7))
    axes=[grid[0,0],grid[1,0],grid[1,1],grid[0,1]]
    phot,info=REVIEW.WEB.review_ztf(name)
    for band,color in [('g','#178c65'),('r','#d45151')]:
        v=np.asarray(phot.get(band,[]),float)
        if len(v):axes[0].scatter(Time(v[:,0],format='mjd').decimalyear,v[:,1],s=4,alpha=.4,color=color,label='catalog '+band)
    alerts=DEST/'alerts'/f'{name}_corrected.csv'
    if alerts.exists():
        d=pd.read_csv(alerts)
        for fid,band,color in [(1,'g','#178c65'),(2,'r','#d45151')]:
            a=d[d.fid.eq(fid)&d.mjd.ge(Time('2025-01-01').mjd)]
            if len(a):axes[0].scatter(Time(a.mjd.to_numpy(),format='mjd').decimalyear,a.magpsf_corr,s=13,facecolors='none',edgecolors=color,label='alert '+band)
    accepted=epochs[epochs.name.eq(name)&epochs.accepted].sort_values('mjd')
    dates=accepted.mjd.dropna().unique()
    # Preserve the photometric date range when old reference spectra predate it.
    phot_limits=axes[0].get_xlim()
    for day in dates:axes[0].axvline(Time(day,format='mjd').decimalyear,color='.4',alpha=.15,lw=.7)
    axes[0].axvline(Time('2026-09-23').decimalyear,color='#457aa4',ls='--',lw=1)
    axes[0].invert_yaxis();axes[0].set_title('Catalog photometry + separate alert detections',fontsize=9)
    if axes[0].get_legend_handles_labels()[0]:axes[0].legend(fontsize=6,loc='best')
    axes[0].set_xlabel('Year');axes[0].set_ylabel('Magnitude (total-light estimates)')
    if axes[0].collections:axes[0].set_xlim(phot_limits[0],2026.9)
    neo=wise_visits().get(name)
    if neo is not None:
        axes[3].errorbar(Time(neo.mjd.to_numpy(),format='mjd').decimalyear,neo.w1,
                        yerr=neo.w1err if 'w1err' in neo else None,fmt='o-',ms=3,lw=.7,color='#ad751c')
        for day in dates:axes[3].axvline(Time(day,format='mjd').decimalyear,color='.4',alpha=.2,lw=.7)
        axes[3].set_xlim(min(2013.5,Time(neo.mjd.min(),format='mjd').decimalyear-.2),2026.9)
        axes[3].invert_yaxis()
    else:axes[3].text(.5,.5,'NEOWISE visit photometry not cached',ha='center',transform=axes[3].transAxes)
    axes[3].axvline(2026+265/365.25,color='#457aa4',ls='--',lw=1)
    axes[3].set_title('NEOWISE W1: spectral dates marked; no extrapolation',fontsize=9)
    axes[3].set_xlabel('Year');axes[3].set_ylabel('W1 magnitude (Vega)')
    # Six-Angstrom plotting caches are adequate for inspection, not significance.
    records=REVIEW.records(name);good=[]
    for r in records:
        m=REVIEW.measure(r,t['z'])
        if not m['accepted']:continue
        a=accepted[accepted.mjd.eq(m['mjd'])&accepted.source.eq(m['source'])&accepted.proprietary.eq(m['proprietary'])]
        if a.empty:continue
        good.append((r,m))
    good.sort(key=lambda pair:pair[1]['mjd'])
    unique={}
    for r,m in good:unique.setdefault((int(m['mjd']),m['source']),(r,m))
    good=list(unique.values())
    if len(good)>5:good=[good[j] for j in np.unique(np.linspace(0,len(good)-1,5).astype(int))]
    colors=plt.cm.viridis(np.linspace(.08,.9,max(1,len(good))))
    for (r,m),color in zip(good,colors):
        w=np.asarray(r['wave'],float)/(1+t['z']);f=np.asarray(r['flux'],float)
        for ax,lo,hi in [(axes[1],4700,5150),(axes[2],6400,6800)]:
            mask=np.isfinite(f)&(w>=lo)&(w<=hi)
            ax.plot(w[mask],f[mask]/m['continuum_flux'],lw=.85,color=color,label=m['date']+('*' if m['coadd'] else ''))
    for ax,lo,hi,lines,title in [(axes[1],4700,5150,[4862.7,5008.2],'Hβ and [O III]'),(axes[2],6400,6800,[6564.6],'Hα and neighboring lines')]:
        ax.set_xlim(lo,hi);ax.set_title(title,fontsize=9);ax.set_xlabel('Rest wavelength (Å)')
        for wave in lines:ax.axvline(wave,color='.6',ls=':',lw=.7)
        if ax.lines:
            vals=np.concatenate([line.get_ydata() for line in ax.lines if len(line.get_ydata())>3]) if any(len(line.get_ydata())>3 for line in ax.lines) else np.array([0,2])
            vals=vals[np.isfinite(vals)]
            if len(vals):ax.set_ylim(min(-.2,float(np.percentile(vals,1))),max(2,float(np.percentile(vals,98))*1.15))
    if good:axes[1].legend(fontsize=6,loc='upper left')
    axes[1].set_ylabel('Flux / local Hβ continuum')
    fig.suptitle(f"{name} | z={t['z']:.3f} | {t['experiment_group']} | reference {t.get('reference_date','unknown')}",fontsize=10)
    fig.text(.01,.005,'Shape comparison only: continua normalized separately; host, Fe II and aperture differences remain. * dated coadd. Alert sampling is incomplete; no alerts does not mean no change.',fontsize=7)
    fig.tight_layout(rect=(0,.045,1,.95));fig.savefig(path,dpi=125);plt.close(fig)
    return path


def main():
    e=pd.read_csv(DEST/'all_candidate_evidence.csv',low_memory=False)
    epochs=pd.read_csv(DEST/'spectral_epoch_audit.csv')
    alerts=pd.read_csv(DEST/'public_alert_summary.csv')
    e=e.merge(alerts,on='name',how='left')
    e['selection_group']=e.experiment_group
    known=e.known_state_status.eq('catalog-confirmed CLAGN')
    e.loc[known,'experiment_group']='E: known CLAGN follow-up'
    e.loc[known,'palomar_test']='Measure recurrence or broad-line evolution relative to the latest reliable spectrum; audit the published transition before claiming another state change.'
    e['recent_alert_clue']='No robust recent change established by this alert comparison; missing detections are not evidence of stability.'
    for band in ['g','r']:
        flag=e[f'{band}_2026_minus_2025'].abs().ge(.3)
        e.loc[flag,'recent_alert_clue']='Alert-selected 2026 versus 2025 measurements differ by at least 0.3 mag in '+band+'; check reference calibration and sampling before treating this as a trigger.'
    pairs=controls(e)
    # Deliberately small science discussion sets; no claim of schedule optimality.
    september=['P3456','P4675','P12144','P10825','P9595','P10640','V116651104','P9227','P9506','P11113','P10596','P12457','S1319703171177670656','P11530','P12889']
    october=['P21422','P21743','P16145','P15439','P7346','P8086','P16590','P18343','P20498','P22470','P23155']
    # These are explicit review choices, never constraints on a later schedule.
    selected=list(dict.fromkeys(september+october+e.loc[e.old_primary,'name'].tolist()+pairs.manifold_target.head(8).tolist()+pairs.comparison_target.head(8).tolist()))
    q=e.set_index('name',drop=False).loc[[n for n in selected if n in set(e.name)]].copy()
    q['discussion_role']=np.where(q.name.isin(september),'September science comparison',np.where(q.name.isin(october),'October science comparison','Original target or control comparison'))
    q['visual_review_note']=q.name.map(VISUAL_NOTES).fillna('Not visually classified in this pass; inspect profiles and native errors.')
    q.to_csv(DEST/'proposed_science_discussion.csv',index=False)
    for i,t in enumerate(q.to_dict('records'),1):
        diagnostic(t,epochs);print('Diagnostic',i,'/',len(q),t['name'],flush=True)
    show=['name','discussion_role','experiment_group','reference_date','accepted_dates','optical_timing_test','trigger_direction','optical_change_abs','neowise_change_after_reference_mag','last_alert_date','g_2026_minus_2025','r_2026_minus_2025','field_status','z','r_planning','last_reference_private','known_state_status']
    (DEST/'science_discussion_table.md').write_text(q[show].to_markdown(index=False)+'\n')
    def esc(v):return html.escape(str(v))
    def fmt(v):return 'unknown' if pd.isna(v) else f'{v:.2f}' if isinstance(v,(float,np.floating)) else str(v)
    cards=[]
    for t in q.to_dict('records'):
        cut=ROOT/'data/cutouts'/f"{t['name']}_sdss.jpg"
        im='<img class="cut" src="data:image/jpeg;base64,'+base64.b64encode(cut.read_bytes()).decode()+'" alt="Cached SDSS field">' if cut.exists() else '<p>SDSS cutout not yet cached.</p>'
        facts=f"Reference {t['reference_date']}; {t['accepted_dates']} accepted spectral dates/interval endpoints. Optical screen: {t['optical_timing_test']}. Δmag={fmt(t.get('optical_change_abs'))}; latest public alert {fmt(t.get('last_alert_date'))}."
        cards.append(f'<article data-night="{esc(t["eligible_nights"])}" data-name="{esc(t["name"])}"><h2>{esc(t["name"])} <small>{esc(t["discussion_role"])}</small></h2><p>{esc(t["experiment_group"])} · {esc(t["known_state_status"])}</p><p>{esc(facts)}</p><p><b>Test:</b> {esc(t["palomar_test"])}</p><p><b>Recent alerts:</b> {esc(t["recent_alert_clue"])}</p><p><b>Field:</b> {esc(t["field_status"])} — {esc(t["field_notes"])}</p><p><b>Qualitative profile review:</b> {esc(t["visual_review_note"])}</p>{im}<img class="diagnostic" loading="lazy" src="figures/{esc(t["name"])}.png" alt="Light curves and archival Balmer spectral comparison for {esc(t["name"])}"></article>')
    page='''<!doctype html><html><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>CLAGN science review</title><style>body{background:#101722;color:#dce5ef;font:16px/1.5 system-ui;max-width:1500px;margin:auto;padding:24px}a{color:#9acaff}h1{font-size:30px}small{font-size:14px;color:#a7b9ca}article{border:1px solid #354658;background:#182230;padding:20px;margin:24px 0;border-radius:8px}.diagnostic{width:100%;background:white}.cut{width:180px;display:block;margin-bottom:12px}button,input{font:inherit;padding:8px;background:#223449;color:white;border:1px solid #567}.notice{border-left:4px solid #caa56b;padding:12px}nav{position:sticky;top:0;background:#101722;padding:10px;display:flex;gap:10px;flex-wrap:wrap}</style></head><body><h1>Science review — independent of the observing queue</h1><p class="notice">Local collaboration review. These are candidate tests, not confirmed changing-look classifications or scheduled primaries. Date-aware spectra include usable short coadds; reductions from the same night are not independent transitions. No object is protected because it appeared in an earlier packet.</p><p><a href="README.md">Scientific assessment and experiment protocol</a> · <a href="proposed_science_discussion.csv">Discussion table</a> · <a href="possible_matched_controls.csv">Provisional control matches</a> · <a href="all_candidate_evidence.csv">Full evidence inventory</a></p><nav><button onclick="filter('')">All</button><button onclick="filter('sep23')">September 23 eligible</button><button onclick="filter('oct26')">October 26 eligible</button><button onclick="filter('oct27')">October 27 eligible</button><input id="search" placeholder="Find target" oninput="filter(current)"></nav>'''+''.join(cards)+'''<script>let current='';function filter(n){current=n;const s=document.getElementById('search').value.toLowerCase();document.querySelectorAll('article').forEach(a=>a.hidden=!(a.dataset.night.includes(n)&&a.dataset.name.toLowerCase().includes(s)))}</script></body></html>'''
    (DEST/'index.html').write_text(page)
    print('Control pairs passing preliminary calipers',len(pairs),flush=True)


if __name__=='__main__':main()

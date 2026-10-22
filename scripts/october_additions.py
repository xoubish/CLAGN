"""Explicit October 22 additions, separate from the historical October 14 audit."""
from pathlib import Path
import json
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'data/additional_targets'
DEST=ROOT/'observing/oct26_27_2026'
DATE='2026-10-22'
DEFERRED = {
    'HNGC-0814': 'Deferred by user: public SDSS-V spectrum from 2024-12-08 already available; no established variability motivating a repeat. Non-variability is not proven.',
    'HNGC-1377': 'Deferred by user: archival optical spectra already available, including 6dFGS, nuclear/integrated spectra and DBSP; no established variability motivating a repeat. Non-variability is not proven.',
}


def apply_additions(targets):
    targets['requested_addition']=False
    targets['requested_night']=''
    targets['airmass_limit']=np.where(targets.is_jwst,2.,1.8)
    path=BASE/'helou_review/feasibility.csv'
    if not path.exists():raise FileNotFoundError('Run 64_helou_targets.py before rebuilding the October 22 queue')
    rows=pd.read_csv(path).to_dict('records')
    h=json.loads((BASE/'H5MUSES-327/summary.json').read_text())
    rows.append(dict(name='H5MUSES-327',ra=h['ra'],dec=h['dec'],z=np.nan,
        central_r=h['ztf']['r']['latest180_mag'],exposures=2,airmass_limit=1.8,
        viable_nights='oct26',snr_scenario=np.nan))
    for r in rows:
        name=r['name']; payload=json.loads((BASE/name/'target_payload.json').read_text())
        viable=isinstance(r['viable_nights'],str) and bool(r['viable_nights'])
        if viable and name!='H5MUSES-327':
            assert r['central_r']<19 and r['snr_scenario']>=5, (name,r)
        nexp=int(r['exposures']);overhead=6 if nexp==4 else 4
        row=dict(ra=r['ra'],dec=r['dec'],z=r['z'],r_planning=r['central_r'],r_adopted=r['central_r'],
            brightness_source=payload['r_source'],brightness_epoch='archival; not an October 2026 measurement',
            coordinate_source='Spitzer 5MUSES' if name=='H5MUSES-327' else 'SIMBAD ICRS nuclear position; see cached identity review',
            is_jwst=False,jwst_id='',jwst_family='',jwst_target='',observed_sep23=False,
            eligible_science=viable,requested_addition=viable,requested_night='oct26' if name=='H5MUSES-327' else '',
            airmass_limit=r['airmass_limit'],pool_role='requested',field_status='review',
            field_notes=payload['science']['field_notes'],science_question=payload['science']['science_question'],
            science_score=250.,priority_reason='Explicit October 22 requested optical observation',
            exposures=nexp,visit_minutes=nexp*5+overhead,overhead_minutes=overhead,protected_followup=False)
        if name=='HIRAS-F23290-1252':
            row['priority_reason']='Requested obscured nucleus: 4x300 s at X<=1.5 for marginal blue continuum'
        if name in DEFERRED:
            row.update(eligible_science=False, requested_addition=False, priority_reason=DEFERRED[name])
        targets.loc[name,list(row)]=list(row.values())
    return targets


def snapshot():
    path=DEST/'before_additions_20261022.csv'
    if not path.exists():path.write_bytes((DEST/'all_80_targets.csv').read_bytes())


def previous_names():
    return set(pd.read_csv(DEST/'before_additions_20261022.csv').NAME)


def write_audit(targets,queue,write_csv):
    old=pd.read_csv(DEST/'before_additions_20261022.csv').rename(columns=str.lower).set_index('name')
    new=queue.set_index('name')
    rows=[]
    for name in sorted(set(old.index)|set(new.index)|set(DEFERRED)):
        status='retained' if name in old.index and name in new.index else 'added' if name in new.index else 'returned to candidate pool'
        if name in DEFERRED:status='deferred by user; retained in candidate catalogue'
        rows.append(dict(name=name,decision=status,previous_night=old.loc[name,'night'] if name in old.index else '',
            scheduled_night=new.loc[name,'night'] if name in new.index else '',
            previous_start=old.loc[name,'start_pdt'] if name in old.index else '',start_pdt=new.loc[name,'start_pdt'] if name in new.index else '',
            deferral_reason=DEFERRED.get(name,'')))
    write_csv(pd.DataFrame(rows),DEST/'addition_decisions_20261022.csv')
    selected=queue[queue.name.str.startswith('H')]
    lines=['# October 22: requested H targets','',
        'All five proposal identities remain in the candidate catalogue with H prefixes. Two are scheduled: HIC-1953 and HIRAS-F23290-1252. HUGC-03097 has no complete Moon-safe visit on either night. HNGC-0814 and HNGC-1377 are deferred at the user\'s request because optical spectra are already available and no variability motivating a repeat is established; this is not proof of non-variability. H5MUSES-327 is also scheduled on October 26, giving three H targets in total. The two nights remain at 40 targets each, including all previously protected observations.','',
        selected[['name','night','start_pdt','exposures','r_adopted','airmass_max_actual','moon_min']].to_markdown(index=False,floatfmt='.2f'),'',
        'IRAS F23290 is resolved as F23290-1252 by its IRAS 25um flux and matching SIMBAD redshift. UGC 03097 uses SIMBAD z=0.012014; the proposal value 0.01917 conflicts with both SIMBAD and LASr (0.01197).','',
        'Brightness for the galaxies is measured from calibrated Legacy Survey images in a central 1.5 x 3 arcsec rectangle, taking the faintest of sampled position angles. It includes host starlight. This is not the whole-galaxy magnitude or the obscured-core flux. No dereddening/aperture correction is applied. Broadband-normalized point-source ETC values are scenarios; image seeing, centering, resolved structure, absorption features, bright-sky spectrum and extinction remain limitations.','',
        'HIRAS-F23290-1252 needs 4x300 s at X<=1.5 (26 minutes including 6 minutes overhead), rather than the normal 2x300 s. P8544 remains 3x300 s on October 27. Check the first exposure and use available slack for inadequate counts; optical non-detection of the buried core does not reject its activity.','',
        'H5MUSES-327 has r~15.76 and small optical scatter, but Gaia parallax/proper motion strongly favors a foreground stellar optical counterpart. Its observation tests identity; it is not an established variable AGN.','',
        '[Full before/after audit](addition_decisions_20261022.csv). Lower-ranked ordinary targets displaced to keep 40/night remain in the catalogue. The original October 14 review remains separately recorded.','',
        'Evidence: data/additional_targets/helou_review/{proposal.txt,targets.json,iras_prefix.txt,feasibility.csv}; per-target aperture_photometry.json, target_payload.json, visibility.json and public catalogue responses.']
    (DEST/'H_target_additions.md').write_text('\n'.join(lines)+'\n')

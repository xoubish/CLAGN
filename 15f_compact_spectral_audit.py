"""Count dated optical spectral coverage and catalog CLAGN labels, locally.

Counts are lower bounds from the available inventories, not a complete archive
search or a new broad-line state classification. All outputs remain git-ignored.
"""
from pathlib import Path
from datetime import datetime, timezone
import json
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord, search_around_sky
import astropy.units as u

ROOT = Path(__file__).resolve().parent
DATA = ROOT / 'data'
OUT = DATA / 'reselection_2026-09-20'


def main():
    targets = pd.read_csv(OUT / 'compact_review_objects.csv')
    names = set(targets.name)
    aliases = pd.read_parquet(OUT / 'parent_aliases.parquet')
    amap = aliases.drop_duplicates('name').set_index('name').canonical_name
    refs = aliases[aliases.canonical_name.isin(names) & aliases.reference_mjd.between(40000, 70000)].copy()
    refs = refs.rename(columns={'canonical_name': 'target_name', 'reference_mjd': 'mjd'})
    refs['survey_family'] = 'SDSS'; refs['inventory'] = 'parent reference spectra'
    pieces = [refs[['target_name','mjd','survey_family','inventory']]]
    sv = pd.read_parquet(OUT / 'sdssv_epochs_matched.parquet')
    sv = sv[sv.canonical_name.isin(names)].copy()
    sv['target_name'] = sv.canonical_name
    sv['survey_family'] = 'SDSS'; sv['inventory'] = 'SDSS-V internal tagged/rolling daily'
    pieces.append(sv[['target_name','mjd','survey_family','inventory','metadata_quality_ok']])
    coverage = set()
    for path in sorted(DATA.glob('spectra_epochs_*.csv')):
        e = pd.read_csv(path, low_memory=False)
        e['target_name'] = e.name.str.strip().map(amap)
        e = e[e.target_name.isin(names)].copy()
        if e.empty:
            continue
        if path.name in ['spectra_epochs_pool.csv','spectra_epochs_zeltyn.csv','spectra_epochs_v2pub.csv','spectra_epochs_three_night_public.csv','spectra_epochs_three_night_public_dr19.csv','spectra_epochs_three_night_public_dr20.csv']:
            coverage.update(e.target_name)
        optical = e.source.eq('DESI')
        if 'instrument' in e:
            optical |= e.source.eq('SDSS') & e.instrument.str.lower().isin(['sdss','boss'])
        else:
            # These local products were built from optical BOSS/SDSS catalogs.
            optical |= e.source.eq('SDSS')
        e = e[optical & pd.to_numeric(e.mjd, errors='coerce').between(40000,70000)].copy()
        if 'coadd' in e:
            # Epoch and all-epoch coadds can reuse daily observations.
            e = e[~e.coadd.fillna('').astype(str).str.lower().isin(['epoch','allepoch'])]
        if 'ra' in e and 'dec' in e and len(e):
            pos = targets.set_index('name').loc[e.target_name]
            sep = SkyCoord(e.ra.to_numpy()*u.deg,e.dec.to_numpy()*u.deg).separation(
                SkyCoord(pos.ra.to_numpy()*u.deg,pos.dec.to_numpy()*u.deg)).arcsec
            e = e[sep <= 2].copy()
        e['survey_family'] = e.source; e['inventory'] = path.name
        cols = ['target_name','mjd','survey_family','inventory']
        cols += [x for x in ['min_mjd','max_mjd','coadd_numnight'] if x in e]
        pieces.append(e[cols])
    records = pd.concat(pieces,ignore_index=True)
    records['epoch_day'] = np.floor(records.mjd).astype(int)
    records.to_csv(OUT / 'compact_spectral_inventory_records.csv',index=False)
    # Same observing date is one epoch even when present in several reductions.
    # DESI mean dates summarize coadds; never count their component nights here.
    epochs = records.groupby(['target_name','epoch_day'],as_index=False).agg(
        surveys=('survey_family',lambda x:';'.join(sorted(set(x)))),
        inventories=('inventory',lambda x:';'.join(sorted(set(x)))))
    epochs.to_csv(OUT / 'compact_spectral_epochs.csv',index=False)
    result = targets[['name','ra','dec','r_planning','review_region','eligible_nights']].copy()
    result['n_optical_epoch_dates_min'] = result.name.map(epochs.groupby('target_name').size()).fillna(0).astype(int)
    for label, sub in [('sdssv',sv), ('sdssv_metadata_quality',sv[sv.metadata_quality_ok])]:
        result['n_'+label+'_dates'] = result.name.map(sub.groupby('target_name').mjd.nunique()).fillna(0).astype(int)
    result['has_cached_public_inventory'] = result.name.isin(coverage)
    result['has_desi_coadd'] = result.name.isin(records.loc[records.survey_family.eq('DESI'),'target_name'])
    result['first_epoch_day'] = result.name.map(epochs.groupby('target_name').epoch_day.min())
    result['latest_epoch_day'] = result.name.map(epochs.groupby('target_name').epoch_day.max())
    result['baseline_years'] = (result.latest_epoch_day-result.first_epoch_day)/365.25

    coords = SkyCoord(targets.ra.to_numpy()*u.deg,targets.dec.to_numpy()*u.deg)
    matches = []
    catalogs = [
        ('Camus-Panda local compilation', DATA/'external/clagn_catalog_camus_panda2026.csv', 'ra_deg','dec_deg'),
        ('earlier literature list', DATA/'literature_clagn.csv','ra','dec'),
        ('Zeltyn 2024', DATA/'zeltyn_coords.csv','ra','dec')]
    for source, path, ra, dec in catalogs:
        cat = pd.read_csv(path)
        cat = cat[pd.to_numeric(cat[ra],errors='coerce').notna() & pd.to_numeric(cat[dec],errors='coerce').notna()].reset_index(drop=True)
        ii,jj,sep,_ = search_around_sky(coords,SkyCoord(cat[ra].to_numpy()*u.deg,cat[dec].to_numpy()*u.deg),2*u.arcsec)
        for i,j,s in zip(ii,jj,sep.arcsec):
            row = cat.iloc[j]
            matches.append(dict(name=targets.iloc[i]['name'],catalog=source,sep_arcsec=s,
                catalog_name=row.get('main_name',row.get('name','')),
                status=row.get('confirmation_status',row.get('class_zeltyn','literature-list match; confirmation not audited')),
                transition=row.get('transition_type',''),reference=row.get('reference',row.get('ref','Zeltyn 2024'))))
    matches = pd.DataFrame(matches)
    matches.to_csv(OUT/'compact_known_state_matches.csv',index=False)
    confirmed = set(matches.loc[matches.status.isin(['spectroscopic_confirmed','CL-AGN']),'name'])
    candidate = set(matches.loc[matches.status.str.contains('candidate',case=False,na=False),'name'])-confirmed
    other = set(matches.name)-confirmed-candidate
    result['known_state_status'] = np.select([result.name.isin(confirmed),result.name.isin(candidate),result.name.isin(other)],
        ['catalog-confirmed CLAGN','reported candidate','literature match; confirmation needs audit'],default='no match in checked catalogs')
    result.to_csv(OUT/'compact_spectral_audit.csv',index=False)
    stats = dict(updated_utc=datetime.now(timezone.utc).isoformat(),n_targets=len(result),catalog_confirmed_clagn=len(confirmed),reported_candidates=len(candidate),
        additional_literature_matches=len(other),no_known_catalog_match=int((~result.name.isin(set(matches.name))).sum()),
        multiple_optical_epoch_dates_min=int(result.n_optical_epoch_dates_min.ge(2).sum()),
        one_optical_epoch_date_identified=int(result.n_optical_epoch_dates_min.eq(1).sum()),
        sdssv_any=int(result.n_sdssv_dates.ge(1).sum()),sdssv_multiple=int(result.n_sdssv_dates.ge(2).sum()),
        sdssv_metadata_quality_any=int(result.n_sdssv_metadata_quality_dates.ge(1).sum()),
        sdssv_metadata_quality_multiple=int(result.n_sdssv_metadata_quality_dates.ge(2).sum()),
        has_desi_coadd=int(result.has_desi_coadd.sum()),has_cached_public_inventory=int(result.has_cached_public_inventory.sum()))
    (OUT/'compact_spectral_audit_summary.json').write_text(json.dumps(stats,indent=2))
    description = f'''# Spectral history of the compact pool

Snapshot: {stats['updated_utc']}; {len(result)} distinct targets. This audit pertains
to this saved compact pool, not to subsequently added projections.

## Previous changing-look reports

The checked local compilations match {len(confirmed)} catalog-confirmed CLAGNs,
{len(candidate)} additional reported candidates, and {len(other)} additional
literature-list objects whose confirmation has not been audited. These are
positional matches within 2 arcsec, not new classifications from spectral fits.
No match is not proof that a source is unpublished. Membership in a manifold
candidate region does not itself establish a previous or current state change.

{matches[matches.name.isin(confirmed)][['name','catalog_name','transition','reference']].drop_duplicates().to_markdown(index=False)}

## Spectral coverage

- At least {stats['multiple_optical_epoch_dates_min']} objects have optical spectra
  at two or more distinct dates identified in the available inventories.
- Only one date is currently identified for {stats['one_optical_epoch_date_identified']};
  this does not mean that no additional spectrum exists.
- {stats['sdssv_any']} have SDSS-V daily spectra; {stats['sdssv_multiple']} have more
  than one SDSS-V date. Of these, {stats['sdssv_metadata_quality_any']} have at least
  one and {stats['sdssv_metadata_quality_multiple']} have at least two dates passing
  warning=0, median S/N>=5, and good field quality. These are metadata checks, not
  broad-line S/N or state validation.
- {stats['has_desi_coadd']} have a DESI coadd in the cached inventory.
- {int((result.n_optical_epoch_dates_min.ge(2) & result.known_state_status.eq('no match in checked catalogs')).sum())}
  repeat-spectrum targets have no match in the checked changing-look catalogs.

Counts combine parent reference dates, matched internal SDSS-V tagged/rolling daily rows,
and cached public optical inventories. Same-day entries are merged across
reductions and surveys; all-epoch and SDSS-V epoch coadds are excluded to avoid
double-counting the daily data. APOGEE and MaNGA entries are excluded. Each DESI
coadd contributes its mean date once; its component nights are not counted as
separately usable spectra. Dates are inventory detections, not downloaded or
visually verified usable spectra.

Only {stats['has_cached_public_inventory']} targets have a matched cached public
archive inventory. This remains a lower-bound count; consult the acquisition
manifests for empty searches, unavailable files, and remaining archive failures.
The rolling daily metadata have been reconciled where field summaries are
accessible. Literature reports may refer to spectra absent from this cache.
Multiple spectra alone do not establish a broad-line transition. No scientific
weights were assigned by this audit.

Files: `compact_spectral_audit.csv` (one row per source),
`compact_known_state_matches.csv` (all catalog matches),
`compact_spectral_epochs.csv` (deduplicated dates), and
`compact_spectral_inventory_records.csv` (provenance before date deduplication).
'''
    (OUT/'COMPACT_SPECTRAL_AUDIT.md').write_text(description)
    print(json.dumps(stats,indent=2))
    print(result.groupby('known_state_status').n_optical_epoch_dates_min.agg(['count','min','max']).to_string())


if __name__ == '__main__':
    main()

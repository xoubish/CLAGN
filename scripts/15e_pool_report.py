"""Write an auditable phase-one summary and plots, including incomplete coverage."""
from pathlib import Path
import json,os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'data/reselection_2026-09-20'


def compact_review(p, g):
    """Operational review subset; neither a yield model nor a final target list."""
    config=OUT/'review_selection.json'
    if config.exists() and json.loads(config.read_text()).get('selection_phase')=='prepared':
        return  # A descriptive parent report must not replace prepared lists.
    if config.exists() and json.loads(config.read_text()).get('moon_min_deg')==40:
        import importlib
        importlib.import_module('23_operational_review').build_review(p)
        return
    eligible = p[p.has_projection & p.r_planning.le(18.5) &
                 (p.in_region_clagn.eq(True) | p.in_region_zeltyn.eq(True))].copy()
    windows = g[g.preferred_longest_minutes.ge(60)].merge(eligible, on='name')
    windows['review_region'] = np.select(
        [windows.in_region_clagn.eq(True) & windows.in_region_zeltyn.eq(True), windows.in_region_clagn.eq(True)],
        ['both', 'literature turn-on/off'], default='Zeltyn')
    # Sorting is deterministic, not a science ranking. Preserve every passing night.
    windows = windows.sort_values(['night', 'ra', 'name'])
    windows.to_csv(OUT/'compact_review_object_nights.csv', index=False)
    nights = windows.groupby('name').night.agg(lambda x: ','.join(sorted(set(x))))
    objects = windows.sort_values(['name', 'preferred_longest_minutes', 'night'],
                                 ascending=[True, False, True]).drop_duplicates('name').copy()
    objects = objects.rename(columns={'night': 'best_night_by_window_length'})
    objects['eligible_nights'] = objects.name.map(nights)
    objects = objects.sort_values(['ra', 'name'])
    objects.to_csv(OUT/'compact_review_objects.csv', index=False)
    counts = windows.groupby('night').name.nunique().reindex(['sep23','oct26','oct27'], fill_value=0)
    report = f'''# Compact working pool

{len(objects):,} distinct objects, with {len(windows):,} qualifying object-night windows.
This is a provisional first spectral-review pass, not the final observing list.

Selection is the union of the existing literature turn-on/off and Zeltyn manifold
regions, historical r <= 18.5, and at least 60 continuous minutes with airmass <= 1.5
and apparent Moon separation >= 60 degrees on at least one allocated night.
Both geometric conditions must hold simultaneously at every tested interval endpoint.
The brighter limit and stronger geometry cut reduce exposure costs under the bright Moon;
they do not establish broad-line S/N or replace a line-specific exposure calculation.

{counts.rename('objects').to_frame().to_markdown()}

The two regions remain labeled separately. Their selection and validation histories
differ; membership does not imply equal discovery yield, a calibrated probability,
or a prediction of the current broad-line state. No arbitrary neighbor-fraction
threshold or top-N science ranking was introduced to force this size.

`compact_review_objects.csv` has one row per source, all qualifying nights, and
geometry from its longest qualifying window. `compact_review_object_nights.csv`
retains every qualifying night. Neither file assigns an observing time.

Coverage remains incomplete until the expanded W1 batch finishes. The same cuts
are reapplied when the report refreshes, so newly projected sources can enter and
the count can grow. Missing projections are not classified as uninteresting.
The full parent and broader queues remain available for scientifically justified
exceptions; this working cut is not a deletion from the search population.

Next inspect actual light curves and dated spectra, especially the latest SDSS-V
broad-line state, to establish what a Palomar spectrum would newly test. Check
current brightness, redshift/classification, line coverage, and attainable limits
before assigning science weights. All products remain local.
'''
    (OUT/'COMPACT_REVIEW.md').write_text(report)
    print('Compact first-pass review:', len(objects), 'distinct objects;', counts.to_dict(), flush=True)


def main():
    p=pd.read_parquet(OUT/'parent_with_manifold_descriptors.parquet')
    g=pd.read_parquet(OUT/'night_geometry.parquet')
    compact_review(p, g)
    funnel=pd.read_csv(OUT/'selection_funnel.csv')
    aliases=pd.read_parquet(OUT/'parent_aliases.parquet',columns=['name','canonical_name'])
    amap=aliases.drop_duplicates('name').set_index('name').canonical_name
    old=[]
    for n in ['sep23','oct26','oct27']:
        t=pd.read_csv(ROOT/f'data/targets_{n}.csv');old.extend(t.loc[t['rank']>0,'name'])
    old=set(pd.Series(old).map(amap).dropna())
    v=g[g.moon40_longest_minutes>=30].merge(p,on='name')
    v['was_in_old_primary_list']=v.name.isin(old)
    view=v[v.has_projection & v.r_planning.le(19.5)].copy()
    # Raw manifold support is an ordering aid, never a forecast or success probability.
    view=view.sort_values(['manifold_cl_neighbor_fraction','usable_longest_minutes'],ascending=False)
    cols=['name','night','ra','dec','z','r_planning','r_source','r_epoch_mjd','manifold_interest',
          'manifold_on_neighbor_fraction','manifold_off_neighbor_fraction','manifold_cl_neighbor_fraction',
          'moon40_longest_minutes','usable_longest_minutes','preferred_longest_minutes','moon_min_visible',
          'moon_max_visible','min_airmass','sdssv_n_quality_epochs','sdssv_latest_quality_mjd',
          'was_in_old_primary_list','origins','aliases']
    view[cols].to_csv(OUT/'manifold_observable_review_queue.csv',index=False)
    fig,axs=plt.subplots(1,3,figsize=(15,5),sharex=True,sharey=True,layout='constrained')
    a=pd.read_csv(ROOT/'data/sampleA_embedding_objectid.csv')
    for ax,n in zip(axs,['sep23','oct26','oct27']):
        ax.scatter(a.umap_x,a.umap_y,c='lightgrey',s=5,alpha=.5)
        s=view[view.night==n]
        sc=ax.scatter(s.umap_x,s.umap_y,c=s.manifold_direction,s=7,cmap='coolwarm',vmin=-.5,vmax=.5,alpha=.65)
        prev=s[s.was_in_old_primary_list];ax.scatter(prev.umap_x,prev.umap_y,facecolors='none',edgecolors='black',s=40,lw=.7)
        ax.set_title(f'{n}\n{len(s):,} projected, observable targets; r≤19.5',fontsize=11);ax.set_xlabel('W1 UMAP 1')
    axs[0].set_ylabel('W1 UMAP 2')
    fig.colorbar(sc,ax=axs.ravel().tolist(),label='Turn-on fraction − turn-off fraction',shrink=.8)
    fig.suptitle('Reopened manifold pool; black circles mark previous primary targets\nColors describe labeled neighbors, not probabilities of the current state',fontsize=12)
    fig.savefig(OUT/'reopened_manifold.png',dpi=150,bbox_inches='tight');plt.close(fig)
    status=[]
    for batch in sorted(OUT.glob('wise_r*')):
        if not batch.is_dir():continue
        nfiles=len(list(batch.glob('pixel_*.parquet')))
        target_count=len(pd.read_csv(batch/'targets.csv')) if (batch/'targets.csv').exists() else 0
        projected=len(pd.read_csv(batch/'projected.csv')) if (batch/'projected.csv').exists() else 0
        status.append(f'- {batch.name}: {target_count:,} queued objects; {nfiles:,} downloaded sky partitions; {projected:,} completed projections.')
    txt=f'''# Rebuilt parent and observing feasibility

This is the manifold-first parent reconstruction requested on 2026-09-20. **Scientific weights and a final observing schedule have not been assigned.** The existing published schedules remain historical products with known erroneous Moon distances; these files are the replacement analysis inputs.

## Parent and scope

- {len(p):,} coordinate groups after merging 1,127,067 input rows at 2 arcsec. Close pairs and questionable catalog identifications still require review.
- Public query: all RA and declination, QSO or AGN-subclass GALAXY, 0.001<z<1.14; no brightness, variability, Moon, or redshift-quality veto. This is the Hbeta-accessible follow-up branch, not a claim to include every AGN or every redshift.
- Added the corresponding branch of the local 750,414-object DR16Q catalog, the entire cached galaxy parent, Zeltyn sources, and objects from the complete cached SDSS-V catalog, including previously unmatched objects.
- The SDSS-V input contains 917,117 distinct cached epochs for 573,582 object keys, through MJD 61235. Tagged v6_2_1 is the newest available tagged reduction. The rolling master catalog is newer (server modification 2026-09-20) and must be reconciled before final spectral-recency weighting.
- The paper's original object-ID/coordinate table is absent, so its training light curves cannot all be added as pointing targets yet. The saved paper manifold is used as the reference representation. Their absence is not interpreted as non-observability.
- {int(p.has_projection.sum()):,} objects currently have a usable saved or newly computed W1 projection. Others remain explicit pending cases, not low-manifold-score objects.

## Observability and brightness

Apparent Moon and target coordinates are evaluated in the same topocentric frame. Every accepted five-minute interval passes the constraints at both endpoints; the shorter final interval uses its actual duration. Visibility, 30/40/60-degree Moon cuts, and X<2/1.8/1.5 alternatives are retained separately. Local cached Earth-orientation data give approximately arcsecond-level limitations, negligible for these selection thresholds.

Brightness is an explicitly labeled historical SDSS or ZTF measurement, or an SDSS-V aperture synthetic magnitude when imaging is missing. The latter has not been independently spectrophotometrically validated. Unknown and faint objects remain in the parent. The r≤20.5 W1 acquisition batch is a practical first sensitivity tier, not a deletion of fainter/unknown objects. Actual exposure costs await line-specific ETC evaluation and updated photometry.

{funnel.to_markdown(index=False)}

These columns are successive operational counts; nightly lists overlap. No row count is an expected CLAGN yield. The literature-neighborhood label means a paper-defined region or at least 10% CL-labeled neighbors among 50 reference objects; other projected regions remain available. Neighbor fractions are descriptive, not calibrated spectral-state probabilities. Objects above z=1 are flagged as outside the paper's nominal redshift domain, and Hbeta near the instrument edge is flagged separately.

## W1 expansion status

{chr(10).join(status)}

New light curves use the same unWISE 2010–2020 catalog, original GP/normalization routine, and saved model as the existing pool. Exact DTW distances are computed in parallel; a regression comparison with the original transform gave identical embeddings on the verification batch. The model is not retrained to create a different map. Sampling span and extrapolation fraction are recorded for new projections.

## Files for the next decision

- `manifold_observable_review_queue.csv`: every currently projected object with r≤19.5 and a ≥30-minute Moon>40°, X<2 window, across all manifold regions. It is sorted by raw manifold support for inspection, not by finalized science return.
- `observable_sep23.csv`, `observable_oct26.csv`, `observable_oct27.csv`: broader queues, including faint and unknown-brightness objects and challenging windows.
- `parent.parquet`, `parent_aliases.parquet`, `night_geometry.parquet`: retained parent, aliases/provenance, and all tested object-night combinations.
- `sdssv_epochs_matched.parquet`: dated SDSS-V inventory. `metadata_quality_ok` checks warning, S/N and field-quality fields; it does **not** classify broad-line state or prove cross-epoch flux calibration.
- `projection_todo.csv`: feasible objects still missing W1 projections, including unknown brightness.
- `recovered_from_old_moon_*.csv`: previously excluded objects rescued by the corrected Moon geometry.
- `reopened_manifold.png`: the reopened feasible pool on the fixed W1 map.

All products here remain local and git-ignored because some contain collaboration metadata. No public website was published or modified.

## Spectral weighting comes next

For objects supported by the manifold, assess the latest usable broad-line state and its date, the timing of W1/NEOWISE changes relative to that spectrum, optical evidence after the latest spectrum, and what Palomar could newly distinguish. Preserve separate outcomes for a first transition, persistence, recurrence, and a test of an uncertain prediction. Do not automatically remove recent SDSS-V targets or reward a long gap. Establish usable broad-line baselines and attainable line limits before converting these ingredients into observing weights.
'''
    (OUT/'POOL_REBUILD.md').write_text(txt)
    print('Report written;',len(view),'projected bright object-night rows;',view.name.nunique(),'distinct objects',flush=True)


if __name__=='__main__':main()

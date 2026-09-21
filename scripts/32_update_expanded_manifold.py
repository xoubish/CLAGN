"""Apply completed W1 batches without recomputing obsolete sky cuts.

The parent remains broad. Internal galaxy/stellar targets without AGN evidence
are retained for audit but do not enter the observing candidate selection.
"""
from pathlib import Path
import argparse,importlib,json
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'data/reselection_2026-09-20'


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--select',action='store_true');args=ap.parse_args()
    p=pd.read_parquet(OUT/'parent_with_manifold_descriptors.parquet')
    for file in sorted(OUT.glob('wise_r*/projected.csv')):
        new=pd.read_csv(file).set_index('name');match=p.name.isin(new.index)
        for col in ['umap_x','umap_y','fvar_w1','mean_w1_mjy','n_w1','w1_first_mjd','w1_last_mjd','w1_span_days','gp_extrapolated_fraction']:
            if col in new:p.loc[match,col]=p.loc[match,'name'].map(new[col])
        p.loc[match,'has_projection']=True;p.loc[match,'manifold_status']='projected in expanded pool'
    a=pd.read_csv(ROOT/'data/sampleA_embedding_objectid.csv');use=p.has_projection
    _,xe,ye=np.histogram2d(a.umap_x,a.umap_y,bins=10)
    def bins(x,y):return np.clip(np.searchsorted(xe[1:],x),0,9),np.clip(np.searchsorted(ye[1:],y),0,9)
    ax,ay=bins(a.umap_x,a.umap_y);px,py=bins(p.loc[use,'umap_x'],p.loc[use,'umap_y'])
    for col in ['in_region_clagn','in_region_zeltyn']:
        mask=np.zeros((10,10),bool);sel=a[col].fillna(False).astype(bool).to_numpy();mask[ax[sel],ay[sel]]=True;p.loc[use,col]=mask[px,py]
    distances,idx=cKDTree(a[['umap_x','umap_y']].to_numpy()).query(p.loc[use,['umap_x','umap_y']].to_numpy(),k=51)
    idx=np.where((distances[:,0]<1e-6)[:,None],idx[:,1:],idx[:,:50]);bits=a.label_bits.to_numpy().astype(int)
    for label,bit in [('on',16),('off',32),('cl',48)]:p.loc[use,f'manifold_{label}_neighbor_fraction']=((bits&bit)>0)[idx].mean(axis=1)
    p['manifold_direction']=p.manifold_on_neighbor_fraction-p.manifold_off_neighbor_fraction
    p['agn_identity_status']='Public spectroscopic QSO / AGN parent'
    private=p.origin.eq('SDSSV_internal_parent')
    sv=pd.read_parquet(OUT/'sdssv_epochs_matched.parquet')
    good=sv[sv.metadata_quality_ok].copy()
    positive=good.cls.eq('QSO')|good.subclass.fillna('').str.contains('AGN',case=False)
    candidates=good.cls.eq('GALAXY')&good.firstcarton.fillna('').str.contains(r'bhm_(aqmes|gua|rm|spiders_agn|csc)',regex=True)
    confirmed=set(good.loc[positive,'canonical_name']);candidate=set(good.loc[candidates,'canonical_name'])-confirmed
    p.loc[private,'agn_identity_status']='AGN identity not established; retain outside observing selection'
    p.loc[private&p.name.isin(candidate),'agn_identity_status']='Galaxy spectrum in AGN-targeted carton; identity needs review'
    p.loc[private&p.name.isin(confirmed),'agn_identity_status']='Internal spectrum classified QSO/AGN with accepted metadata'
    # Candidate AGN cartons alone do not satisfy the agreed AGN identity gate.
    p['agn_identity_eligible']=~private|p.name.isin(confirmed)
    p.to_parquet(OUT/'parent_with_manifold_descriptors.parquet',index=False)
    g=pd.read_parquet(OUT/'three_night_geometry_r19_30min.parquet')
    observable=p[p.name.isin(g.name)]
    candidates=observable[observable.has_projection&observable.agn_identity_eligible&(observable.in_region_clagn.eq(True)|observable.in_region_zeltyn.eq(True))]
    counts=g[g.name.isin(candidates.name)].groupby('night').size().to_dict()
    summary=dict(bright_observable_parent=len(observable),projected=int(observable.has_projection.sum()),pending_projection=int((~observable.has_projection).sum()),
                 candidate_region_objects=len(candidates),night_counts=counts,origins=candidates.origin.value_counts().to_dict(),
                 agn_identity_pending=int((~observable.agn_identity_eligible).sum()))
    (OUT/'three_night_expansion_status.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2),flush=True)
    if args.select:importlib.import_module('23_operational_review').build_review(p,g)


if __name__=='__main__':main()

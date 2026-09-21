"""Add newly identified observable rolling AGN without rebuilding old aliases.

All new identities and spectral-aperture photometry remain in the private
research directory. Existing source coordinates and public identities survive.
"""
from pathlib import Path
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'data/reselection_2026-09-20'


def main():
    additions=pd.read_csv(OUT/'new_rolling_bright_observable_agn_pending.csv')
    parent=pd.read_parquet(OUT/'parent_with_manifold_descriptors.parquet')
    additions=additions[~additions.name.isin(parent.name)].copy()
    if len(additions):
        new=additions[['name','ra','dec','z','r_planning']].copy()
        new['reference_mjd']=additions.mjd
        new['origin']='SDSSV_internal_parent';new['origins']=new.origin
        new['aliases']=new.name;new['r_catalog']=np.nan
        new['r_source']='SDSS-V synthetic aperture r; calibration unverified'
        new['r_epoch_mjd']=additions.mjd
        for col in ['has_projection','in_region_clagn','in_region_zeltyn']:new[col]=False
        new['Hb_A']=4861.33*(1+new.z);new['Ha_A']=6562.8*(1+new.z)
        new['Hb_in_NGPS']=new.Hb_A.between(3050,10400)
        new['Hb_line_edge_flag']=~new.Hb_A.between(3150,10300)
        new['outside_paper_z_domain']=new.z>1
        new['manifold_status']='awaiting W1/projection'
        new['agn_identity_status']='Internal spectrum classified QSO/AGN with accepted metadata'
        new['agn_identity_eligible']=True
        parent=pd.concat([parent,new],ignore_index=True)
        aliases=pd.read_parquet(OUT/'parent_aliases.parquet')
        alias=new.copy();alias['canonical_name']=alias.name;alias['psfmag_r']=np.nan
        alias['group']=np.arange(len(alias))+int(aliases.group.max())+1
        pd.concat([aliases,alias.reindex(columns=aliases.columns)],ignore_index=True).to_parquet(OUT/'parent_aliases.parquet',index=False)
    mapping=pd.read_csv(OUT/'new_rolling_bright_observable_agn_pending.csv').set_index('key').name
    fresh=pd.read_parquet(OUT/'sdssv_rolling_field_epochs.parquet')
    fresh=fresh[fresh.key.isin(mapping.index)].copy();fresh['canonical_name']=fresh.key.map(mapping)
    # These rows were excluded earlier because the identity did not yet exist
    # in the parent; the rolling table already carries the metadata quality cut.
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    positions=parent.set_index('name').loc[fresh.canonical_name]
    fresh['position_offset_arcsec']=SkyCoord(fresh.ra.to_numpy()*u.deg,fresh.dec.to_numpy()*u.deg).separation(SkyCoord(positions.ra.to_numpy()*u.deg,positions.dec.to_numpy()*u.deg)).arcsec
    fresh['metadata_quality_ok']&=fresh.position_offset_arcsec<=1
    old=pd.read_parquet(OUT/'sdssv_epochs_matched.parquet')
    epochs=pd.concat([old,fresh],ignore_index=True).drop_duplicates(['canonical_name','field','mjd','catalogid','archive_version'],keep='last')
    epochs.to_parquet(OUT/'sdssv_epochs_matched.parquet',index=False)
    for label,rows in [('any',epochs),('quality',epochs[epochs.metadata_quality_ok])]:
        last=rows.sort_values(['mjd','sn_median_all']).drop_duplicates('canonical_name',keep='last').set_index('canonical_name')
        for col in ['mjd','sn_median_all','cls','fieldquality']:
            parent[f'sdssv_latest_{label}_{col}']=parent.name.map(last[col])
        parent[f'sdssv_n_{label}_epochs']=parent.name.map(rows.groupby('canonical_name').mjd.nunique()).fillna(0).astype(int)
        flux=parent.name.map(last.spectroflux_r)
        parent[f'sdssv_latest_{label}_synthetic_r']=22.5-2.5*np.log10(flux.where(flux>0))
    parent.to_parquet(OUT/'parent_with_manifold_descriptors.parquet',index=False)
    g=pd.read_parquet(OUT/'three_night_geometry_r19_30min.parquet')
    new_g=pd.read_csv(OUT/'new_rolling_agn_geometry.csv')
    pd.concat([g,new_g],ignore_index=True).drop_duplicates(['name','night'],keep='last').to_parquet(OUT/'three_night_geometry_r19_30min.parquet',index=False)
    bright=pd.read_parquet(OUT/'three_night_bright_parent.parquet')
    bright=pd.concat([bright,parent[parent.name.isin(mapping.values)]],ignore_index=True).drop_duplicates('name',keep='last')
    bright.to_parquet(OUT/'three_night_bright_parent.parquet',index=False)
    print('Added parent identities',len(additions),'rolling identities represented',len(mapping),flush=True)


if __name__=='__main__':main()

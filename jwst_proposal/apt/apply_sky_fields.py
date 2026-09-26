"""Apply visually reviewed sky positions without changing science exposures."""
import csv
import io
import json
from pathlib import Path
import shutil
import xml.etree.ElementTree as E
import zipfile
from astropy.coordinates import SkyCoord
from astropy import units as u

HERE=Path(__file__).resolve().parent
NS='http://www.stsci.edu/JWST/APT'
MRS=NS+'/Template/MiriMRS'
n={'j':NS,'m':MRS}


def main():
    selections=[json.loads(p.read_text()) for p in sorted((HERE/'sky_fields').glob('*_selection.json'))]
    assert len(selections)==24 and all(r['visual_review']=='accepted' for r in selections)
    chosen={r['id']:r for r in selections}
    path=HERE/'clagn24_miri_draft.aptx'
    backup=HERE/'work/before_sky_review';backup.mkdir(parents=True,exist_ok=True)
    for suffix in ['aptx','xml','times','timing.json']:
        src=HERE/f'clagn24_miri_draft.{suffix}'
        if not (backup/src.name).exists():shutil.copy2(src,backup/src.name)
    with zipfile.ZipFile(path) as z:
        files={name:z.read(name) for name in z.namelist()}
    name=next(k for k in files if k.endswith('.xml'))
    for _,(prefix,uri) in E.iterparse(io.BytesIO(files[name]),events=['start-ns']):
        E.register_namespace(prefix,uri)
    root=E.fromstring(files[name]);coordinates={};rows=[]
    for t in root.findall('j:Targets/j:Target',n):
        label=t.findtext('j:TargetName',namespaces=n)
        if not label.endswith('-SKY'):continue
        r=chosen[label[:-4]];s=r['selected'];c=SkyCoord(s['ra_deg']*u.deg,s['dec_deg']*u.deg)
        text=c.ra.to_string(unit=u.hourangle,sep=' ',precision=3,pad=True)+' '+c.dec.to_string(sep=' ',precision=2,pad=True,alwayssign=True)
        coordinates[label]=text
        t.find('j:EquatorialCoordinates',n).set('Value',text)
        t.find('j:Extended',n).text='YES'
        t.find('j:Comments',n).text=f"Sky field selected from AllWISE/2MASS catalogs and WISE W1/W3/W4 + 2MASS J images. Offset {s['radius_arcsec']} arcsec at PA {s['pa_deg']} deg. Nearest AllWISE/2MASS sources {s['nearest_allwise_arcsec']:.1f}/{s['nearest_2mass_arcsec']:.1f} arcsec. See sky_fields/sky_field_review.pdf."
        rows.append(dict(id=r['id'],target=r['target'],ra_deg=s['ra_deg'],dec_deg=s['dec_deg'],coordinates=text,
                         offset_arcsec=s['radius_arcsec'],pa_deg=s['pa_deg'],nearest_allwise_arcsec=s['nearest_allwise_arcsec'],
                         nearest_2mass_arcsec=s['nearest_2mass_arcsec'],moved=r['moved'],visual_review='accepted'))
    group=root.find('j:DataRequests/j:ObservationGroup',n)
    group.find('j:Label',n).text='24 AGN - MIRI/MRS'
    group.find('j:Comments',n).text='12 fading and 12 rising AGNs. Uniform A/B/C groups 27/60/27, four dithers and one integration. Dedicated sky fields selected from archival infrared catalogs and images; source/sky pairs consecutive and noninterruptible.'
    for obs in root.findall('.//j:ObservationGroup/j:Observation',n):
        if '-SKY' not in obs.findtext('j:TargetID',namespaces=n):continue
        m=obs.find('j:Template/m:MiriMRS',n)
        assert m.findtext('m:AcqTargetID',namespaces=n)=='NONE'
        m.find('m:Dithers/m:MrsDitherSpecification/m:OptimizedFor',n).text='EXTENDED SOURCE'
    # Coordinates and dithers have changed; cached scheduling must be recomputed.
    for parent in root.iter():
        for child in list(parent):
            if child.tag=='{'+NS+'}ToolData':parent.remove(child)
    E.indent(root);xml=E.tostring(root,encoding='UTF-8',xml_declaration=True)
    files[name]=xml
    with zipfile.ZipFile(path,'w',zipfile.ZIP_DEFLATED) as z:
        for key,value in files.items():z.writestr(key,value)
    (HERE/'clagn24_miri_draft.xml').write_bytes(xml)
    inv=list(csv.DictReader((HERE/'observation_inventory.csv').open()))
    for row in inv:
        if row['target'] in coordinates:row['coordinates']=coordinates[row['target']]
        row['extended']='YES' if row['role']=='sky' else 'NO'
        row['dither_optimized_for']='EXTENDED SOURCE' if row['role']=='sky' else 'POINT SOURCE'
    for dest,records in [(HERE/'observation_inventory.csv',inv),(HERE/'sky_fields/selected_sky_fields.csv',rows)]:
        with dest.open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=list(records[0]));w.writeheader();w.writerows(records)
    print('Applied 24 reviewed sky fields and extended-source background dithers. Recompute APT timing now.')


if __name__=='__main__':main()

"""Apply reviewed exposure and acquisition recipes; retain a pre-edit backup."""
import csv,io,shutil,zipfile
from pathlib import Path
import xml.etree.ElementTree as E
HERE=Path(__file__).resolve().parent
APT=HERE.parent/'apt'
NS='http://www.stsci.edu/JWST/APT';MRS=NS+'/Template/MiriMRS';n={'j':NS,'m':MRS}
recipes={r['id']:r for r in csv.DictReader((HERE/'acquisition_recipes.csv').open())}
assert len(recipes)==24
backup=APT/'work/before_technical_review';backup.mkdir(parents=True,exist_ok=True)
for suffix in ['aptx','xml','times','timing.json','pointing']:
 p=APT/f'clagn24_miri_draft.{suffix}'
 if not (backup/p.name).exists():shutil.copy2(p,backup/p.name)
path=APT/'clagn24_miri_draft.aptx'
with zipfile.ZipFile(path) as z:files={k:z.read(k) for k in z.namelist()}
name=next(k for k in files if k.endswith('.xml'))
for _,(prefix,uri) in E.iterparse(io.BytesIO(files[name]),events=['start-ns']):E.register_namespace(prefix,uri)
root=E.fromstring(files[name])
for o in root.findall('.//j:ObservationGroup/j:Observation',n):
 label=o.findtext('j:TargetID',namespaces=n).split()[-1];ident=label.removesuffix('-SKY')
 m=o.find('j:Template/m:MiriMRS',n)
 if not label.endswith('-SKY'):
  r=recipes[ident]
  for key,val in [('AcqFilter',r['filter']),('AcqReadoutPattern',r['readout']),('AcqGroups',r['groups'])]:m.find('m:'+key,n).text=val
 for e in m.findall('m:ExposureList/m:Exposure',n):
  setting=e.findtext('m:Wavelength',namespaces=n)
  groups=80 if ident=='R06' and setting=='MEDIUM(B)' else 60 if setting=='MEDIUM(B)' else 27
  if groups:
   for key in ['GroupsLong','GroupsShort']:e.find('m:'+key,n).text=str(groups)
root.find('.//j:ObservationGroup/j:Comments',n).text='24 AGN, 12 fading and 12 rising, with paired sky visits. Baseline A/B/C groups 27/60/27; R06 uses 80 groups in B; R01/P9694 uses the baseline sequence. Acquisitions use F560W or FND with flux-dependent readout settings; see review/acquisition_recipes.csv.'
for parent in root.iter():
 for child in list(parent):
  if child.tag=='{'+NS+'}ToolData':parent.remove(child)
E.indent(root);xml=E.tostring(root,encoding='UTF-8',xml_declaration=True);files[name]=xml
with zipfile.ZipFile(path,'w',zipfile.ZIP_DEFLATED) as z:
 for k,v in files.items():z.writestr(k,v)
(APT/'clagn24_miri_draft.xml').write_bytes(xml)
rows=list(csv.DictReader((APT/'observation_inventory.csv').open()))
for r in rows:
 ident=r['target'].removesuffix('-SKY')
 if ident=='R06':r['groups_B']='80'
 r['groups_A']='27'
 if r['role']=='science':
  a=recipes[ident];r['target_acquisition']=f"{a['filter']} {a['readout']} {a['groups']} groups"
with (APT/'observation_inventory.csv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
print('Updated 24 acquisitions and both science/sky exposure pairs for R01 and R06. APT timing must be recomputed.')

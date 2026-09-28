"""Check the saved APT product against recipes, prose, generator and layout."""
from pathlib import Path
import csv, hashlib, json, re, subprocess, tempfile, zipfile
import xml.etree.ElementTree as E
from datetime import datetime, timezone
from pypdf import PdfReader
HERE=Path(__file__).resolve().parent; BASE=HERE.parent; APT=BASE/'apt'
n={'j':'http://www.stsci.edu/JWST/APT','m':'http://www.stsci.edu/JWST/APT/Template/MiriMRS'}
with zipfile.ZipFile(APT/'clagn24_miri_draft.aptx') as z:
 xml=z.read(next(k for k in z.namelist() if k.endswith('.xml')))
assert xml==(APT/'clagn24_miri_draft.xml').read_bytes()
r=E.fromstring(xml); observations=r.findall('.//j:ObservationGroup/j:Observation',n)
assert len(observations)==48
inventory=list(csv.DictReader((APT/'observation_inventory.csv').open()))
recipes={a['id']:a for a in csv.DictReader((HERE/'acquisition_recipes.csv').open())}
assert len(recipes)==24 and len(inventory)==48
for obs,row in zip(observations,inventory):
 assert obs.findtext('j:Number',namespaces=n)==row['observation']
 assert obs.findtext('j:TargetID',namespaces=n).split()[-1]==row['target']
 assert obs.findtext('j:Instrument',namespaces=n)=='MIRI'
 ident=row['target'].removesuffix('-SKY');m=obs.find('j:Template/m:MiriMRS',n)
 exposures=m.findall('m:ExposureList/m:Exposure',n);assert len(exposures)==3
 for exp,key in zip(exposures,['A','B','C']):
  expected=80 if ident=='R06' and key=='B' else 60 if key=='B' else 27
  for det in ['Short','Long']:
   assert int(exp.findtext('m:Groups'+det,namespaces=n))==expected==int(row['groups_'+key])
   assert exp.findtext('m:Integrations'+det,namespaces=n)=='1'
   # APT 2026 saves the shared detector readout in ReadoutPatternLong.
   readout=exp.findtext('m:ReadoutPattern'+det,namespaces=n)
   assert readout=='FASTR1' or (det=='Short' and readout is None and exp.findtext('m:ReadoutPatternLong',namespaces=n)=='FASTR1')
 if row['role']=='science':
  a=recipes[ident]
  for key,col in [('AcqFilter','filter'),('AcqReadoutPattern','readout'),('AcqGroups','groups')]:
   assert m.findtext('m:'+key,namespaces=n)==a[col]
 else:
  assert m.findtext('m:AcqTargetID',namespaces=n)=='NONE'
  assert m.findtext('.//m:OptimizedFor',namespaces=n)=='EXTENDED SOURCE'
links=r.findall('j:LinkingRequirements/j:GroupWithinLink',n)
assert len(links)==24 and all(x.get('Sequence')==x.get('NonInterruptible')=='true' for x in links)
windows=[]
for e in r.findall('.//j:ToolValue',n):
 if (e.get('Name') or '').startswith('Visit Planner:'):
  v=E.fromstring(e.text);w=v.find('StVisitSchedulingWindows')
  pcf=w.get('StSchedulingPCF','').split()
  assert w.get('UpToDate')=='true' and any(float(x)>0 for x in pcf[1::2])
  windows.append({'visit':e.get('Name'),'engine':v.get('EngineVersion'),'up_to_date':True,'has_windows':True})
assert len(windows)==48
info=r.find('j:ProposalInformation',n)
request=float(info.findtext('j:ChargedTime',namespaces=n));assert request==67.0
assert info.findtext('j:ProposalSize',namespaces=n)=='MEDIUM'
abstract=(BASE/'apt_title_abstract.txt').read_text().split('Abstract:\n')[1].strip()
assert abstract==info.findtext('j:Abstract',namespaces=n).strip()
title=(BASE/'apt_title_abstract.txt').read_text().split('Title:\n',1)[1].split('\n\nAbstract:',1)[0].strip()
assert title==info.findtext('j:Title',namespaces=n).strip()
timing=json.loads((APT/'clagn24_miri_draft.timing.json').read_text())
assert timing['number_of_visits']==48 and timing['charged_time']==67.0
assert timing['charged_time']<=request<timing['charged_time']+.1
with tempfile.TemporaryDirectory() as tmp:
 subprocess.run(['python3',str(APT/'build_draft.py'),'--output-dir',tmp],check=True,capture_output=True)
 generated=list(csv.DictReader((Path(tmp)/'observation_inventory.csv').open()))
 assert generated==inventory
pdf=PdfReader(BASE/'proposal.pdf');pages=[p.extract_text() for p in pdf.pages];text='\n'.join(pages)
assert len(pages)==7 and 'Description of Observations' in pages[4] and 'Supplemental Information' in pages[5]
assert '67.0' in pages[4] and 'Generative AI disclosure' in pages[6] and 'template for JWST Cycle 6' in pages[6]
assert not re.search(r'NIRSpec|@|October|provisional|TODO',text,re.I)
# Bibliographic author names are permitted; screen identifying names in proposal prose.
assert 'References' in text
assert not re.search(r'Hemmati|Shoubaneh',text.split('References',1)[0],re.I)
assert not any(pdf.metadata.get(k) for k in ['/Author','/Subject'])
assert not re.search(r'Overfull|LaTeX Warning', (BASE/'proposal.log').read_text())
style_hash=hashlib.sha256((BASE/'jwstproposaltemplate_v6.sty').read_bytes()).hexdigest()
assert style_hash=='d5e8b0c1133974ac9eaa6cc27213794ceec297337c69f2ec13254a10684fdd08'
result=dict(checked_at_utc=datetime.now(timezone.utc).isoformat(),passed=True,
 observations=48,science_targets=24,science_sky_pairs=24,all_groups_and_acquisitions_match=True,
 generator_inventory_matches=True,apt_xml_matches_archive=True,abstract_matches=True,
 charged_time_h=timing['charged_time'],requested_time_h=request,science_duration_h=timing['science_time'],
 photon_collection_source_plus_sky_h=sum(sum(int(a['groups_'+k]) for k in ['A','B','C'])*4*2.77504 for a in inventory)/3600,
 scheduling=windows,pages=7,core_pages=5,style_sha256=style_hash,layout_warnings=re.findall(r'Underfull \\hbox[^\n]*', (BASE/'proposal.log').read_text()),
 automated_anonymity_screen_passed=True,visual_review_pages=[5,7],
 limitations='Focused ETC model checks; nuclear flux fractions, individual histories and extended-host acquisition centering are not fully validated.')
(HERE/'consistency.json').write_text(json.dumps(result,indent=2)+'\n')
print('PASS: 48 visits, 24 pairs, recipes/inventory/generator/abstract consistent; 67.0 h charged and requested; five core pages.')

"""Choose shortest TA recipe passing compact-source flux brackets.

Lower flux: half total W3 with Fnu proportional to lambda^2.
Upper flux: full total W3 with flat Fnu. These bracket assumptions are not
measurements of nuclear morphology or a guarantee against future variability.
"""
import sys,json,copy,csv
from pathlib import Path
import numpy as np
P=Path(__file__).resolve().parent
sys.path.insert(0,str(P.parent/'etc'))
from run_local import configure
paths=configure()
from synphot.config import conf
conf.vega_file=str(paths['PYSYN_CDBS']/'calspec/alpha_lyr_stis_010.fits')
from pandeia.engine.calc_utils import build_default_calc
from pandeia.engine.perform_calculation import perform_calculation
base=json.loads((P.parent/'etc/r06_ch4b_60groups_r09/input.json').read_text())
cfg=build_default_calc('jwst','miri','target_acq');cfg['background']=base['background'];cfg.pop('background_level',None)
sample=list(csv.DictReader((P.parent/'inputs/jwst_sample_cycle6.csv').open()))
recipes=[('f560w','fast',4),('f560w','fast',6),('f560w','fast',10),('fnd','fast',10),('fnd','fastgrpavg',10),('fnd','fastgrpavg8',10)]
cache=P/'technical_runs/ta_grid';cache.mkdir(exist_ok=True)
rows=[]
for row in sample:
 chosen=None
 for filt,read,groups in recipes:
  cases=[]
  for label,flux,index in [('faint',float(row['w3_mjy'])/2,2),('bright',float(row['w3_mjy']),0)]:
   path=cache/f"{row['id']}_{filt}_{read}_{groups}_{label}.json"
   if path.exists():r=json.loads(path.read_text())
   else:
    x=copy.deepcopy(cfg);x['scene']=copy.deepcopy(base['scene'])
    x['configuration']['instrument']['filter']=filt
    x['configuration']['detector'].update(ngroup=groups,readout_pattern=read)
    x['scene'][0]['spectrum']['normalization']['norm_flux']=flux
    x['scene'][0]['spectrum']['sed']['index']=index
    result=perform_calculation(x)
    r=json.loads(json.dumps(dict(input=x,scalar=result['scalar'],warnings=result['warnings']),default=lambda v:v.item() if isinstance(v,np.generic) else v.tolist()))
    path.write_text(json.dumps(r,indent=2)+'\n')
   cases.append(r)
  faint,bright=cases
  if faint['scalar']['sn']>=50 and bright['scalar']['fraction_saturation']<=.7 and not faint['warnings'] and not bright['warnings']:
   chosen=dict(id=row['id'],target=row['target'],filter=filt.upper(),readout=read.upper(),groups=groups,lower_snr=faint['scalar']['sn'],upper_saturation=bright['scalar']['fraction_saturation'],ta_time_s=faint['scalar']['total_exposure_time'])
   break
 if chosen is None:raise RuntimeError('No valid recipe for '+row['id'])
 rows.append(chosen);print(chosen,flush=True)
with (P/'acquisition_recipes.csv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

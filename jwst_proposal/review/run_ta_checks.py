"""F560W FAST/FULL 10-group acquisition checks for assumed compact nuclei."""
import sys,json,copy
import numpy as np
from pathlib import Path
P=Path(__file__).resolve().parent
sys.path.insert(0,str(P.parent/'etc'))
from run_local import configure
paths=configure()
from synphot.config import conf
conf.vega_file=str(paths['PYSYN_CDBS']/'calspec/alpha_lyr_stis_010.fits')
from pandeia.engine.calc_utils import build_default_calc
from pandeia.engine.perform_calculation import perform_calculation
base=json.loads((P.parent/'etc/r06_ch4b_60groups_r09/input.json').read_text())
cfg=build_default_calc('jwst','miri','target_acq')
cfg['configuration']['detector']['ngroup']=10
cfg['background']=base['background'];cfg.pop('background_level',None)
results=[]
for name,flux,index,filt,groups in [('faint_red',.95,2,'f560w',10),('bright_red',48.4,2,'f560w',10),('bright_flat',48.4,0,'f560w',10),('faint_red_4',.95,2,'f560w',4),('fnd_faint_red',.95,2,'fnd',10),('fnd_mid_red',2.,2,'fnd',10),('fnd_bright_red',48.4,2,'fnd',10),('fnd_bright_flat',48.4,0,'fnd',10)]:
 x=copy.deepcopy(cfg);x['scene']=copy.deepcopy(base['scene'])
 x['configuration']['instrument']['filter']=filt
 x['configuration']['detector']['ngroup']=groups
 x['scene'][0]['spectrum']['normalization']['norm_flux']=flux
 x['scene'][0]['spectrum']['sed']['index']=index
 dest=P/'technical_runs'/('ta_'+name);dest.mkdir(exist_ok=True)
 (dest/'input.json').write_text(json.dumps(x,indent=2)+'\n')
 if (dest/'report.json').exists():r=json.loads((dest/'report.json').read_text())
 else:
  result=perform_calculation(x)
  r=json.loads(json.dumps(dict(scalar=result['scalar'],warnings=result['warnings']), default=lambda x:x.item() if isinstance(x,np.generic) else x.tolist()))
  (dest/'report.json').write_text(json.dumps(r,indent=2)+'\n')
 results.append(dict(case=name,**r));print(name,json.dumps(r),flush=True)
(P/'ta_checks.json').write_text(json.dumps(results,indent=2)+'\n')

"""Focused sensitivity and bright-envelope checks with the installed Pandeia engine."""
import copy, json, subprocess, sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
ETC=ROOT/'etc'
OUT=ROOT/'review/technical_runs'
base=json.loads((ETC/'r06_ch4b_60groups_r09/input.json').read_text())
jobs=[]
for name,ch,disp,groups,ref,radius in [('r06_rest5','ch1','medium',60,6.0555,.3),('r06_rest8','ch2','medium',60,9.6888,.45)]:
 cfg=copy.deepcopy(base)
 cfg['configuration']['instrument'].update(aperture=ch,disperser=disp)
 cfg['configuration']['detector']['ngroup']=groups
 cfg['strategy'].update(reference_wavelength=ref,aperture_size=radius,target_source=1)
 cfg['scene'][0]['spectrum']['redshift']=.2111
 jobs.append((name,cfg))
# Additional diagnostic anchors at the faint end and high-redshift endpoint.
for name,ch,disp,groups,ref,radius,z,flux,index in [
 ('r06_rest12_80','ch3','medium',80,14.5332,.6,.2111,.95,2),
 ('r06_rest18_80','ch4','medium',80,21.7998,.9,.2111,.95,2),
 ('r06_rest5_80','ch1','medium',80,6.0555,.3,.2111,.95,2),
 ('p1823_rest12_75','ch4','short',75,19.446,.85,.6205,1.3,.585),
 ('r06_rest12','ch3','medium',60,14.5332,.6,.2111,.95,2),
 ('p1823_rest12','ch4','short',27,19.446,.85,.6205,1.3,.585),
 ('p1823_rest8','ch3','short',27,12.964,.6,.6205,1.3,.585)]:
 cfg=copy.deepcopy(base)
 cfg['configuration']['instrument'].update(aperture=ch,disperser=disp)
 cfg['configuration']['detector']['ngroup']=groups
 cfg['strategy'].update(reference_wavelength=ref,aperture_size=radius,target_source=1)
 cfg['scene'][0]['spectrum']['redshift']=z
 cfg['scene'][0]['spectrum']['normalization']['norm_flux']=flux
 cfg['scene'][0]['spectrum']['sed']['index']=index
 jobs.append((name,cfg))
# Flat 48.4 mJy maximizes the short-wave flux relative to our red continuum.
# Long-wave envelope has lambda^2 shape at the largest sample W3 flux.
for ch,ref,radius in [('ch1',6.0,.3),('ch2',9.0,.45),('ch3',14.5,.6),('ch4',21.8,.9)]:
 cfg=copy.deepcopy(base)
 cfg['configuration']['instrument'].update(aperture=ch,disperser='medium')
 cfg['strategy'].update(reference_wavelength=ref,aperture_size=radius,target_source=1)
 cfg['scene'][0]['spectrum']['normalization']['norm_flux']=48.4
 cfg['scene'][0]['spectrum']['sed']['index']=0 if ch in ['ch1','ch2'] else 2
 jobs.append(('bright_envelope_'+ch+'b',cfg))
OUT.mkdir(exist_ok=True)
for name,cfg in jobs:
 dest=OUT/name
 if (dest/'report.json').exists():continue
 inp=OUT/(name+'.json');inp.write_text(json.dumps(cfg,indent=2)+'\n')
 result=subprocess.run([str(ETC/'.venv/bin/python'),str(ETC/'run_local.py'),str(inp),str(dest)],capture_output=True,text=True)
 (OUT/(name+'.log')).write_text(result.stdout+result.stderr)
 if result.returncode:raise RuntimeError(name+': '+result.stderr[-2000:])
 report=json.loads((dest/'report.json').read_text())
 print(name,json.dumps(report['scalar']),flush=True)

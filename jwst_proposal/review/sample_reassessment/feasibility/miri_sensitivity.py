"""Conditional 12-target MIRI anchor sensitivity, not an APT-approved request."""
from pathlib import Path
import sys, json, copy, os, csv
import numpy as np
OUT=Path(__file__).resolve().parent
ROOT=OUT.parents[3]
ETC=ROOT/'jwst_proposal/etc'
sys.path.insert(0,str(ETC))
from run_local import configure
paths=configure()
from synphot.config import conf
conf.vega_file=str(paths['PYSYN_CDBS']/'calspec/alpha_lyr_stis_010.fits')
from pandeia.engine.perform_calculation import perform_calculation

# Nominal intervals from the MRS documentation. Prefer maximal bin edge clearance.
BANDS=[('ch1','short',4.90,5.74,.3),('ch1','medium',5.66,6.63,.3),('ch1','long',6.53,7.65,.3),
       ('ch2','short',7.51,8.77,.45),('ch2','medium',8.67,10.13,.45),('ch2','long',10.01,11.70,.45),
       ('ch3','short',11.55,13.47,.6),('ch3','medium',13.34,15.57,.6),('ch3','long',15.41,17.98,.6),
       ('ch4','short',17.70,20.95,.9),('ch4','medium',20.69,24.48,.9),('ch4','long',24.40,27.90,.9)]
DEFAULT={'short':27,'medium':60,'long':27}

def binned_snr(result,wave):
    w,signal=result['1d']['extracted_flux'];wn,noise=result['1d']['extracted_noise']
    w=np.asarray(w);signal=np.asarray(signal);noise=np.asarray(noise)
    assert np.allclose(w,wn)
    edges=np.r_[w[0]-(w[1]-w[0])/2,(w[1:]+w[:-1])/2,w[-1]+(w[-1]-w[-2])/2]
    lo,hi=wave*(1-1/200),wave*(1+1/200)
    overlap=np.maximum(0,np.minimum(edges[1:],hi)-np.maximum(edges[:-1],lo))
    ok=(overlap>0)&np.isfinite(signal)&np.isfinite(noise)&(noise>0)
    if overlap[ok].sum()<(hi-lo)*.99:raise ValueError('Incomplete R=100 bin')
    return float(np.sum(signal[ok]*overlap[ok])/np.sqrt(np.sum((noise[ok]*overlap[ok])**2)))

def main():
    sample=list(csv.DictReader((OUT/'shortlist.csv').open()))
    base=json.loads((ETC/'r06_ch4b_60groups_r09/input.json').read_text())
    cache=OUT/'miri_cache';cache.mkdir(exist_ok=True)
    rows=[];times=[]
    for r in sample:
        r['z']=float(r['z']);r['w3_mjy']=float(r['w3_mjy'])
        adopted=DEFAULT.copy();anchors=[]
        for rest,goal in [(5.,50.),(12.,50.),(18.,30.)]:
            wave=rest*(1+r['z'])
            bands=[b for b in BANDS if b[2]<wave*(1-1/200) and wave*(1+1/200)<b[3]]
            ch,disp,lo,hi,radius=max(bands,key=lambda b:min(wave-b[2],b[3]-wave))
            cfg=copy.deepcopy(base)
            cfg['configuration']['instrument'].update(aperture=ch,disperser=disp)
            cfg['strategy'].update(reference_wavelength=wave,aperture_size=radius,target_source=1)
            spec=cfg['scene'][0]['spectrum'];spec['redshift']=r['z'];spec['normalization']['norm_flux']=.5*r['w3_mjy']
            def run(groups):
                cfg['configuration']['detector']['ngroup']=groups
                path=cache/f"{r['name']}_rest{int(rest)}_g{groups}.json"
                if path.exists():return json.loads(path.read_text())
                result=perform_calculation(cfg,dict_report=True)
                s=result['scalar'];answer=dict(input=copy.deepcopy(cfg),snr_R100=binned_snr(result,wave),
                    saturation_fraction=float(s.get('fraction_saturation',0)),warnings=result.get('warnings',{}))
                path.write_text(json.dumps(answer,indent=2)+'\n');return answer
            start=DEFAULT[disp];initial=run(start);groups=start;result=initial
            while result['snr_R100']<goal:
                # Conservative sqrt(time) step, then verify with the engine.
                groups=min(500,max(groups+1,int(np.ceil(groups*(goal/result['snr_R100'])**2))))
                result=run(groups)
                if groups==500:break
            # Find the smallest whole group count above the initial recipe that meets the goal.
            lower=start;upper=groups
            while upper-lower>1:
                mid=(lower+upper)//2;trial=run(mid)
                if trial['snr_R100']>=goal:upper=mid
                else:lower=mid
            groups=upper;result=run(groups);adopted[disp]=max(adopted[disp],groups)
            anchors.append(dict(name=r['name'],z=r['z'],rest_um=rest,observed_um=wave,
                channel=ch,setting=disp,nuclear_W3_mjy=.5*r['w3_mjy'],initial_groups=start,
                initial_snr_R100=initial['snr_R100'],goal_snr_R100=goal,
                minimum_groups=groups,verified_snr_R100=result['snr_R100'],warnings=json.dumps(result['warnings'])))
        rows.extend(anchors)
        times.append(dict(name=r['name'],**{f'groups_{k}':v for k,v in adopted.items()},
            source_plus_sky_seconds=2*4*2.77504*sum(adopted.values()),
            extra_source_plus_sky_seconds=2*4*2.77504*(sum(adopted.values())-sum(DEFAULT.values()))))
        for filename,data in [('miri_anchor_sensitivity.csv',rows),('miri_exposure_recipes.csv',times)]:
            with (OUT/filename).open('w') as handle:
                writer=csv.DictWriter(handle,fieldnames=list(data[0]));writer.writeheader();writer.writerows(data)
        print(r['name'],adopted,[(int(a['rest_um']),round(a['initial_snr_R100'],1)) for a in anchors],flush=True)

if __name__=='__main__':main()

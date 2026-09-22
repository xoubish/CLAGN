"""H-beta continuum sensitivity for the review tables; local research only.

Two sets of columns from the official NGPS ETC. The legacy 2x600 s, 1.0 arcsec, 2x2 columns
(snr_sky*) feed the frozen review_order_score and are kept unchanged. The adopted-setting
columns (snr300_*) give S/N per Angstrom for 2x300 s with the 1.5 arcsec slit, 2x3 binning and
seeing 1.3 arcsec at zenith scaled by airmass^0.6. Archival continuum normalization and explicit
sky scenarios are not a 2026 forecast.
"""
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor,as_completed
import hashlib,importlib,json,os
import numpy as np
import pandas as pd
import astropy.units as u

ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'data/reselection_2026-09-20'
CACHE=OUT/'three_night_sensitivity';CACHE.mkdir(exist_ok=True)
os.environ.setdefault('MPLCONFIGDIR',str(OUT/'matplotlib_cache'))


def one(target):
    model=importlib.import_module('22_september_etc')
    reference=model.reference(target)
    base=dict(name=target['name'],z=target['z'],r_planning=target['r_planning'],pool_role=target['pool_role'])
    if reference is None:return base|dict(status='No accepted archival continuum reference')
    mjd,mag,private,continuum=reference
    wave=4862.7*(1+target['z']);lo,hi=wave/10-4,wave/10+4
    base.update(reference_mjd=mjd,continuum_AB=mag,reference_private=private,continuum_flux_1e17=continuum,hbeta_A=wave)
    signature=hashlib.sha256(json.dumps(base,sort_keys=True).encode()+b'central-2x600-v1+adopted-2x300-1p5-2x3-v1').hexdigest()
    path=CACHE/f"{target['name']}.json"
    if path.exists():
        saved=json.loads(path.read_text())
        if saved.get('signature')==signature:return saved
    for ch,rn,scale in zip(model.CFG.channels,[2.8,7.8,3.7,4.6],[.193,.189,.186,.186]):
        model.CFG.readnoise[ch]=rn*u.count/u.pix;model.CFG.platescale[ch]=scale*u.arcsec/u.pix
    channels=[ch for ch in ['R','I','G','U'] if model.CFG.channelRange[ch][0].to_value(u.nm)<=lo and model.CFG.channelRange[ch][1].to_value(u.nm)>=hi]
    if not channels:return base|dict(status='Full H-beta window falls outside a single channel')
    def snr(ch,sky,airmass,seeing,extra):
        cmd=[ch,str(lo),str(hi),'EXPTIME','600','-slit','SET','1.0','-binspect','2','-binspat','2','-seeing',str(seeing),'500','-airmass',str(airmass),'-skymag',str(sky),'-mag',str(mag+extra),'-magsystem','AB','-magfilter','match','-noslicer']
        args=model.ETC.parser.parse_args(cmd);model.ETC.check_inputs_add_units(args)
        return float(model.ETC.main(args,quiet=True)['SNR'].value)*np.sqrt(2)
    def adopted(ch,sky,airmass,extra=0.):
        seeing=1.3*airmass**0.6
        cmd=[ch,str(lo),str(hi),'EXPTIME','300','-slit','SET','1.5','-binspect','3','-binspat','2','-seeing',f'{seeing:.3f}','500','-airmass',str(airmass),'-skymag',str(sky),'-mag',str(mag+extra),'-magsystem','AB','-magfilter','match','-noslicer']
        args=model.ETC.parser.parse_args(cmd);model.ETC.check_inputs_add_units(args)
        bin_A=float((model.CFG.dLambda[ch]*3).to_value(u.AA))
        return float(model.ETC.main(args,quiet=True)['SNR'].value)*np.sqrt(2)/np.sqrt(bin_A)
    options={ch:snr(ch,18.5,1.3,1.3,0) for ch in channels};ch=max(options,key=options.get)
    result=base|dict(status='Scenario calculation',signature=signature,channel=ch,snr_sky18p5_X1p3=options[ch],
        snr_sky18_X1p5=snr(ch,18.,1.5,1.3,0),snr_sky18_X1p5_fainter0p5=snr(ch,18.,1.5,1.3,.5),
        snr_sky18_X1p5_seeing1p8=snr(ch,18.,1.5,1.8,0),
        snr300_x1p3_sky18p5=adopted(ch,18.5,1.3),snr300_x1p5_sky18=adopted(ch,18.,1.5),snr300_x1p5_sky18_fainter0p5=adopted(ch,18.,1.5,.5),snr300_x1p8_sky18=adopted(ch,18.,1.8))
    path.write_text(json.dumps(result,indent=2,allow_nan=False));return result


def main():
    targets=pd.read_csv(OUT/'compact_review_objects.csv');rows=[]
    with ProcessPoolExecutor(max_workers=4) as ex:
        futures={ex.submit(one,r):r['name'] for r in targets.to_dict('records')}
        for i,f in enumerate(as_completed(futures),1):
            try:rows.append(f.result())
            except Exception as exc:rows.append(dict(name=futures[f],status='Calculation failed',reason=str(exc)))
            if i%40==0 or i==len(futures):
                pd.DataFrame(rows).to_csv(OUT/'three_night_2x600_sensitivity.csv',index=False);print('ETC',i,'/',len(futures),flush=True)
    (OUT/'three_night_2x600_assumptions.json').write_text(json.dumps(dict(
        legacy_columns='snr_sky*: 2x600 s, 1.0 arcsec slice, 2x2, zenith seeing passed unscaled; kept only because review_order_score is frozen on snr_sky18p5_X1p3',
        adopted_columns='snr300_*: S/N per Angstrom for 2x300 s, 1.5 arcsec slit, 2x3 binning, seeing 1.3 arcsec x airmass^0.6, sky V 18.5 at X 1.3 or 18.0 at X 1.5/1.8 (adopted 2026-09-21)',
        exposures=2,seconds_each=600,binspect=2,binspat=2,slice_arcsec=1,extraction='central slice only; optimal point-source extraction',
        quantity='Continuum S/N per 2-pixel spectral bin, averaged over +/-40 observed Angstrom around vacuum H-beta; not integrated broad-line significance.',
        normalization='Latest accepted loaded archival spectrum, rest continua 4750-4790 and 5100-5140 Angstrom; no unverified cross-band photometric rescaling.',
        sky='Official ETC sky template uniformly scaled to V=18.5 or 18 mag/arcsec^2; sensitivity cases, not a lunar sky prediction.',
        limitations='Archival aperture flux may include host light; a source may since have changed. A null broad-line result needs a line/host decomposition. Ten minutes of visit overhead is provisional.',
        reference='https://caltechopticalobservatories.github.io/NGPS/users-manual/exposure-time-calculator.html'),indent=2))


if __name__=='__main__':main()

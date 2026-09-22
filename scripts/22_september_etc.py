"""NGPS planning scenarios at the observed H-beta wavelength.

Uses the unmodified official ETC in the ignored ngps_etc cache. Times target
continuum S/N per 2-pixel spectral bin, NOT broad-line detection significance.
Sky values are explicit sensitivity scenarios, not a Moon-brightness forecast.
"""
from pathlib import Path
import sys,json,importlib,math
import numpy as np
import pandas as pd
import astropy.units as u
from astropy.time import Time

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'data/reselection_2026-09-20'
sys.path.insert(0,str(OUT/'ngps_etc'))
ETC=importlib.import_module('ETC_main')
CFG=importlib.import_module('ETC_config')

def calculate(channel,lo,hi,mag,snr,seeing,sky,airmass=1.3,binspect=2,binspat=2,noslicer=False):
    def solve(goal):
        cmd=[channel,str(lo),str(hi),'SNR',str(goal),'-slit','SET','1.0',
             '-binspect',str(binspect),'-binspat',str(binspat),'-seeing',str(seeing),'500',
             '-airmass',str(airmass),'-skymag',str(sky),'-mag',str(mag),'-magsystem','AB','-magfilter','match']
        if noslicer:cmd.append('-noslicer')
        args=ETC.parser.parse_args(cmd);ETC.check_inputs_add_units(args)
        return ETC.main(args,quiet=True)
    first=solve(snr);n=max(1,math.ceil(first['exptime'].to_value(u.s)/900))
    for _ in range(10):
        result=first if n==1 else solve(snr/np.sqrt(n))
        each=float(result['exptime'].to_value(u.s))
        if each<=900:break
        n+=1
    else:raise ValueError('Unable to split within 900 s')
    return dict(seconds=n*each,nexp=n,each_seconds=each,resolution=float(result['resolution'].value),
                bin_angstrom=float((CFG.dLambda[channel]*binspect).to_value(u.AA)))

def reference(t):
    from spectral_utils import accepted_reference
    result = accepted_reference(t)
    if result is None:
        return None
    r, mag, continuum = result
    return r['mjd'], mag, bool(r.get('proprietary')), continuum

def main():
    objects=json.loads((OUT/'sep23_evidence_objects.json').read_text())
    evidence=pd.read_csv(OUT/'sep23_photometric_evidence.csv').set_index('name')
    # May 2026 instrument specification lists channel-dependent read noise and
    # spatial scales. Keep downloaded sources unchanged; record runtime overrides.
    for ch,rn,scale in zip(CFG.channels,[2.8,7.8,3.7,4.6],[.193,.189,.186,.186]):
        CFG.readnoise[ch]=rn*u.count/u.pix
        CFG.platescale[ch]=scale*u.arcsec/u.pix
    manifest=dict(repository_sha=json.loads((OUT/'ngps_etc/github_tree.json').read_text())['sha'],
        spatial_binning=2,spectral_binning=2,slice_width_arcsec=1.,airmass=1.3,
        reference='Latest accepted loaded spectrum; total aperture continuum treated as a point source.',
        spectral_shape='Constant f_nu locally, normalized to measured continuum at H-beta.',
        magnitude_adjustment='If later ZTF exists, apply median available g/r synthetic-to-photometric offset; not a matched-aperture measurement or a 2026 brightness forecast.',
        sky='ETC Gemini template scaled to Johnson V Vega mag/arcsec^2; scenarios, not predicted lunar sky. Spectral-shape uncertainty remains.',
        overrides={'readnoise_e':dict(zip(CFG.channels,[2.8,7.8,3.7,4.6])),
                   'spatial_scale_arcsec':dict(zip(CFG.channels,[.193,.189,.186,.186]))},
        documentation='https://caltechopticalobservatories.github.io/NGPS/technical-specifications.html',
        caveats=['Point-source approximation can be optimistic for host-dominated nuclei.',
                 'Continuum SNR does not establish broad-line significance or a reliable non-detection.',
                 'Sky template does not independently scale lunar continuum and OH lines.',
                 'No throughput loss from clouds or unusual extinction included.'])
    (OUT/'sep23_etc_assumptions.json').write_text(json.dumps(manifest,indent=2))
    rows=[]
    for t in objects:
        ref=reference(t)
        if ref is None:
            rows.append(dict(name=t['name'],status='no accepted continuum reference'));continue
        mjd,mag,private,continuum=ref;e=evidence.loc[t['name']]
        offset=0.
        if pd.notna(e.ztf_last_date) and Time(e.ztf_last_date).mjd>mjd and abs(mjd-e.mjd_spec)<1:
            ds=[e.get(f'delta_{b}_latest_minus_spectrum',np.nan) for b in ['g','r']]
            ds=[v for v in ds if np.isfinite(v)]
            if ds and not bool(e.badspec):offset=float(np.median(ds))
        lam=4862.7*(1+t['z']);lo,hi=lam/10-4,lam/10+4
        channels=[ch for ch in ['R','I','G','U'] if CFG.channelRange[ch][0].to_value(u.nm)<=lo and CFG.channelRange[ch][1].to_value(u.nm)>=hi]
        if not channels:raise ValueError(f'No full H-beta window: {t["name"]}')
        ch=channels[0]
        base=dict(name=t['name'],status='scenario calculation',reference_mjd=mjd,reference_date=Time(mjd,format='mjd').strftime('%Y-%m-%d'),
                  reference_private=private,reference_continuum=continuum,continuum_AB=mag,
                  archival_photometry_offset=offset,planning_continuum_AB=mag+offset,channel=ch,lo_nm=lo,hi_nm=hi)
        for snr in [15,20]:
            for label,seeing,sky,extra in [('reference',1.3,19.,0.),('bright_sky',1.3,18.,0.),('fainter_poorer',1.8,18.,.5)]:
                result=calculate(ch,lo,hi,mag+offset+extra,snr,seeing,sky)
                rows.append(base|dict(goal_snr=snr,scenario=label,seeing_zenith=seeing,sky_V=sky,extra_mag=extra)|result)
        pd.DataFrame(rows).to_csv(OUT/'sep23_etc_scenarios.csv',index=False)
        print(t['name'],ch,round(mag+offset,2),'completed',flush=True)
    pd.DataFrame(rows).to_csv(OUT/'sep23_etc_scenarios.csv',index=False)
    # Controlled comparison at fixed wavelength and S/N per Angstrom equivalent.
    comparisons=[]
    for spat,spect in [(1,1),(2,1),(2,2),(2,3)]:
        # Hold S/N per native spectral pixel equivalent fixed, not S/N per bin.
        result=calculate('R',640,650,18.2,10*np.sqrt(spect),1.3,18.5,binspect=spect,binspat=spat)
        comparisons.append(dict(binspat=spat,binspect=spect,goal_per_bin=10*np.sqrt(spect))|result)
    pd.DataFrame(comparisons).to_csv(OUT/'sep23_binning_comparison.csv',index=False)

if __name__=='__main__':main()

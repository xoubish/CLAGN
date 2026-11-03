"""Independent R-channel exposure checks for the leading Balmer candidates.

The H-beta and continuum windows for these z~0.23--0.25 sources avoid the strong
telluric bands. These native, uncoadded checks do not apply the delivered telluric
correction and are not substituted for the delivered coadds in the fit catalog.
"""
import json
import numpy as np
import pandas as pd
from astropy.io import fits
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import analysis as A
from paper_protocol import OUT

CONFIG={
    'sep23':('sep23_data/reduction_20260924/work/20260924/fluxed_exposures_p330e','P9506'),
    'oct26':('oct26_data/reduction_20261027/work/20261027/fluxed_exposures_feige34','P7947'),
    'oct27':('oct27_data/reduction_20261028/work/20261028/fluxed_exposures_bd284211','P21422;P9506'),
}


def main():
    sample=pd.read_csv(OUT/'pilot_sample.csv').set_index('name');rows=[]
    fig,axes=plt.subplots(4,1,figsize=(11,11));i=0
    for night,(folder,names) in CONFIG.items():
        for name in names.split(';'):
            ax=axes[i];i+=1
            alias='JP9506' if name=='P9506' and night=='oct27' else name
            paths=sorted((A.ROOT/folder/'r').glob('*-'+alias+'_NGPS_r_*'))
            if len(paths)!=2:raise RuntimeError(f'Expected two native exposures for {name} {night}, got {len(paths)}')
            for path in paths:
                with fits.open(A.track(path),memmap=False) as h:
                    d=h[1].data;w=np.asarray(d['OPT_WAVE'],float)
                    f=np.asarray(d['OPT_FLAM'],float);e=np.asarray(d['OPT_FLAM_SIG'],float)
                    valid=np.asarray(d['OPT_MASK'],bool)&np.isfinite(e)&(e>0)
                    sel=np.isfinite(w)&(w>0);order=np.flatnonzero(sel)[np.argsort(w[sel])]
                    w,f,e,valid=[a[order].copy() for a in [w,f,e,valid]];f[~valid]=np.nan
                    mjd=float(h[0].header['MJD'])
                s=dict(wave=w,flux=f,error=e,z=float(sample.loc[name,'z']),quality=np.zeros(len(w),int),
                    channel=np.full(len(w),'R'),source='NGPS')
                for label,shift in A.VARIANTS.items():
                    r=A.line_index(s,'Hb',shift)
                    rows.append(dict(name=name,night=night,mjd=mjd,file=str(path.relative_to(A.ROOT)),variant=label,**r))
                r=A.line_index(s,'Hb');x=w/(1+s['z']);m=(x>4740)&(x<5145)
                ax.plot(x[m],f[m]/r['continuum'],lw=.8,label=f"frame {path.name.split('_')[3].split('-')[0]}: wing EW {r['wing_ew']:.1f} Å")
            ax.set(title=f'{name} {night}: native R exposures, independently extracted',xlim=(4740,5145),ylabel='fλ / local continuum')
            ax.legend(fontsize=8);ax.grid(alpha=.2)
    axes[-1].set_xlabel('Rest vacuum wavelength [Å]');fig.tight_layout()
    fig.savefig(OUT/'figures/independent_exposures.png',dpi=140);plt.close(fig)
    d=pd.DataFrame(rows);d.to_csv(OUT/'independent_exposure_indices.csv',index=False)
    (OUT/'exposure_inputs.json').write_text(json.dumps(A.SOURCES,indent=2)+'\n')
    print(d[d.variant.eq('nominal')][['name','night','mjd','wing_ew','wing_over_oiii','valid']].to_string(index=False))


if __name__=='__main__':main()

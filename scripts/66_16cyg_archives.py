"""Validate component-specific cached spectra and build the 16 Cyg archive card."""
from pathlib import Path
import csv
import json
import os
import shutil
os.environ.setdefault('MPLCONFIGDIR','/tmp/clagn-16cyg-matplotlib')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.io.votable import parse_single_table
from astropy.time import Time

ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/'data/additional_targets/field_2954592_505204'
CACHE=DATA/'spectra'
DEST=ROOT/'docs/october'


def main():
    for release in ['18','20']:
        result=json.loads((CACHE/f'sdss_dr{release}.json').read_text())
        assert next(t['Rows'] for t in result if t['TableName']=='Table1')==[]
    gaia={r['Source']:r for r in csv.DictReader((DATA/'gaia.csv').open())}
    assert gaia['2135550755684554112']['XPcont']=='0'
    assert gaia['2135550755684554112']['XPsamp']=='0'
    assert gaia['2135550755684554112']['RVS']=='0'
    cadc=pd.read_csv(CACHE/'cadc.csv')
    spectra={};components=[]
    for name,obs,source in [('A','1634976','2135550755683405952'),('B','1634979','2135550755683407232')]:
        with fits.open(CACHE/f'{obs}i.fits') as hdus:
            h=hdus[0].header;d=hdus[0].data.astype(float)
            assert h['OBJECT']==f'16 Cyg {name}' and h['RUNID']=='13AB06'
            assert h['INSTRUME']=='ESPaDOnS' and h['DATE-OBS']=='2013-06-29'
            assert h['COL1']=='Wavelength' and h['COL2']=='Intensity'
            assert h.comments['COL2']=='Normalized'
            wave=d[0]*10;flux=d[1];err=d[2]
            good=np.isfinite(wave)&np.isfinite(flux)&np.isfinite(err)&(err>0)
            wave,flux,err=wave[good],flux[good],err[good]
            order=np.argsort(wave);wave,flux,err=wave[order],flux[order],err[order]
            assert 3690<wave.min()<3710 and 10470<wave.max()<10500
            meta=dict(component=name,object=h['OBJECT'],observation_id=obs,
                date=h['DATE-OBS'],exposure_seconds=float(h['EXPTIME']),mode=h['INSTMODE'],
                wavelength_min=float(wave.min()),wavelength_max=float(wave.max()),
                flux_kind='Continuum-normalized; not absolute flux',
                fits_url=json.loads((CACHE/f'{obs}i.query.json').read_text())['url'])
        xp_name=f'XP_SAMPLED_Gaia_DR3_{source}.xml'
        table=parse_single_table(CACHE/xp_name).to_table()
        xp_wave=np.asarray(table['wavelength'],float)*10
        xp_flux=np.asarray(table['flux'],float)
        assert len(table)==343 and str(table['flux'].unit)=='W / (nm m2)'
        public_xp=f'16Cyg-{name}-Gaia-XP.xml';shutil.copyfile(CACHE/xp_name,DEST/public_xp)
        rows=cadc[(cadc.instrument_name=='ESPaDOnS')&(cadc.calibrationLevel==2)&
                  cadc.target_name.str.strip().str.replace('*','',regex=False).str.strip().isin([f'16 Cyg {name}',f'16Cyg{name}'])]
        meta.update(confirmed_observations=int(rows.observationID.nunique()),
            first_date=Time(rows.time_bounds_lower.min(),format='mjd').iso[:10],
            last_date=Time(rows.time_bounds_upper.max(),format='mjd').iso[:10],
            gaia_source_id=source,gaia_url='october/'+public_xp)
        components.append(meta);spectra[name]=(wave,flux,xp_wave,xp_flux)
    # Native normalized spectra have overlapping echelle orders; use equal-width
    # medians only for the overview. Empty bins remain NaN so gaps are visible.
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(3,1,figsize=(11,8.5),gridspec_kw={'height_ratios':[1,1,1]},layout='constrained')
    colors={'A':'#2176b5','B':'#c65b24'}
    for name,(wave,flux,xp_wave,xp_flux) in spectra.items():
        groups=np.floor(wave/2).astype(int);bins=np.arange(groups.min(),groups.max()+1)
        med=pd.Series(flux).groupby(groups).median().reindex(bins).to_numpy()
        axes[0].plot(bins*2+1,med,color=colors[name],lw=.7,alpha=.85,label=f'16 Cyg {name}')
        # Keep separate order segments and native samples in the line view.
        use=(wave>=6550)&(wave<=6575)
        axes[1].plot(wave[use],flux[use],color=colors[name],lw=.9,alpha=.9,label=f'16 Cyg {name}')
        axes[2].plot(xp_wave,xp_flux/1e-12,color=colors[name],lw=1.5,label=f'16 Cyg {name}')
    axes[0].set(title='CFHT / ESPaDOnS · 2013-06-29 · 2 Å median overview',ylabel='Normalized intensity',xlim=(3690,10490),ylim=(-.03,1.2))
    axes[1].set(title='Hα region · native archive wavelength scale',ylabel='Normalized intensity',xlim=(6550,6575),ylim=(0,1.12))
    axes[2].set(title='Gaia DR3 BP/RP · low-resolution mean spectra',ylabel='Flux (10⁻¹² W m⁻² nm⁻¹)',xlim=(3300,10500))
    for a in axes:a.set_xlabel('Wavelength (Å)');a.grid(alpha=.15);a.legend(loc='lower right',frameon=False,ncol=2)
    fig.suptitle('16 Cyg A and B: confirmed archival optical spectra',fontsize=15)
    fig.savefig(DEST/'16Cyg_archival_spectra.png',dpi=160)
    plt.close(fig)
    review=dict(checked_date='2026-10-26',components=components,
        c_status='No securely identified archival optical spectrum found for C in this search.',
        plot='october/16Cyg_archival_spectra.png',
        search_summary='SDSS DR18 and SDSS-V/DR20 allspec: no spectra within 60 arcsec of the system. CADC: A/B observations confirmed; no C-targeted spectrum identified. Gaia DR3: A/B BP/RP spectra downloaded; C has no published BP/RP or RVS spectrum flags. LAMOST DR5: no catalogue rows within 60 arcsec.',
        limitations='LAMOST DR9 returned an unresolved count and empty export; ELODIE timed out. These are not verified non-detections. A/C are only about 3 arcsec apart: A-targeted or unresolved system spectra are not assigned to C.',
        sources=[dict(label='CFHT archive',url='https://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/en/cfht/'),
                 dict(label='ESPaDOnS data format',url='https://www.cfht.hawaii.edu/Instruments/Spectroscopy/Espadons/CLASSICAL/'),
                 dict(label='Published A/B archival study',url='https://academic.oup.com/mnras/article/514/4/4958/6611693')])
    (CACHE/'archive_review.json').write_text(json.dumps(review,indent=2)+'\n')
    (CACHE/'README.md').write_text('''# 16 Cyg archival optical spectra

Checked October 26, 2026. Files and query metadata in this directory record the searches.

A: CFHT/ESPaDOnS 1634976i, 2013-06-29, 280 s, star-only R=80,000.
B: CFHT/ESPaDOnS 1634979i, 2013-06-29, 350 s, star-only R=80,000.
The FITS OBJECT, RUNID (13AB06), and dates match independently published A/B observations. Coordinates are approximate telescope pointings; these are identified by component names and programme, not assigned to the nearest Gaia source. No A exposure is assigned to C.

Rows 0/1/2 are wavelength (nm), normalized intensity and uncertainty, verified against FITS COL keywords and CFHT documentation. The plot converts wavelength to angstroms. The broad overview uses 2-angstrom median bins; H-alpha shows native samples. These normalized spectra cannot establish absolute flux variability. Representative exposures are displayed, not a complete multi-epoch compilation.

Gaia DR3 BP/RP sampled mean spectra were retrieved by exact source IDs for A and B. Wavelength is nm; flux is W m^-2 nm^-1. No C spectrum returned, consistent with Gaia DR3 XPcont, XPsamp and RVS flags all zero.

SDSS DR18 SpecObjAll and DR20 allspec queried a 60-arcsecond radius around the system center, covering all three components and historical proper-motion offsets: successful responses with zero rows. LAMOST DR5 V/164/dr5 also returned no rows. LAMOST DR9 search 71962 returned a sentinel total of 9223372036854775807, empty rows and an empty CSV; this search remains unresolved. ELODIE timed out, not a no-match result.

CADC spectrum-region search returned 692 rows, below the 1000-row cap. Name queries also returned unrelated V1016 Cyg and V1668 Cyg entries, which were excluded. No securely identified C-targeted spectrum was found. This is not proof that none exists in every archive.
''')
    print(json.dumps(components,indent=2))


if __name__=='__main__':main()

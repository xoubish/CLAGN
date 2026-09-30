"""One literature spectrum with a quantitative wavelength-access overlay.

No source-location or flux forecast is implied; active proposal files untouched.
"""
import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/clagn-matplotlib')
from pathlib import Path
from copy import deepcopy
import json, hashlib, subprocess
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from pypdf import PdfReader,PdfWriter,Transformation
from pdf_figure_crop import cropped_figure

BASE=Path(__file__).resolve().parents[1]
OUT=BASE/'review/fig2_three_panel/panel_a';OUT.mkdir(parents=True,exist_ok=True)
protected=[BASE/n for n in ['proposal.tex','fig2_miri.pdf']]
before={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in protected}
sample=pd.read_csv(BASE/'inputs/jwst_sample_cycle6.csv')
zmin,zmax=float(sample.z.min()),float(sample.z.max())
common_mrs=[4.9/(1+zmin),27.9/(1+zmax)]
common_spx=[.75/(1+zmin),5./(1+zmax)]

# Calibrate the published logarithmic abscissa from all five visible major
# ticks in its original embedded figure (4026 x 3355 pixels), not its fluxes.
ticks=np.array([1.5,3.,5.,10.,25.])
tick_x=np.array([3307.5,3462.5,3576.5,3731.5,3936.5])
slope,intercept=np.polyfit(np.log(ticks),tick_x,1)
residual=float(np.max(abs(slope*np.log(ticks)+intercept-tick_x)))
assert residual<.6
# Embedded-image placement from pdftocairo's page-13 SVG transformation.
image_scale=.125562;image_x=42.060844;image_y=53.1633
source=deepcopy(PdfReader(BASE/'review/miri_literature/donnan_2402.17479.pdf').pages[12])
if '/Annots' in source:del source['/Annots']
pageheight=float(source.mediabox.height)
crop=(421.,67.,550.,161.) # page points, top-down; excludes the old shared title
x0,top,x1,bottom=crop;y0,y1=pageheight-bottom,pageheight-top
source=cropped_figure(source,crop)
width,height=6.6*72,4.8*72
place=(.035,.005,.93,.93)
px,py,pw,ph=np.array(place)*[width,height,width,height]
scale=min(pw/(x1-x0),ph/(y1-y0))
tx=px+(pw-(x1-x0)*scale)/2;ty=py+(ph-(y1-y0)*scale)/2
def fx(w):
    original_x=image_x+image_scale*(slope*np.log(w)+intercept)
    return (tx+(original_x-x0)*scale)/width
def fy(pixel_y):
    original_y=pageheight-(image_y+image_scale*pixel_y)
    return (ty+(original_y-y0)*scale)/height

plt.rcParams.update({'font.family':'DejaVu Sans','pdf.fonttype':42})
fig=plt.figure(figsize=(width/72,height/72))
fig.savefig(OUT/'base.pdf');plt.close(fig)
dest=PdfReader(OUT/'base.pdf').pages[0]
dest.merge_transformed_page(source,Transformation().scale(scale).translate(tx,ty))

fig=plt.figure(figsize=(width/72,height/72))
fig.patch.set_alpha(0)
low,high=fy(625),fy(138)
# Limit the SPHEREx shade to the displayed literature plot, whose left spine
# lies at original x=3274. No assertion of complete usable coverage is made.
leftwave=float(np.exp((3277-intercept)/slope))
for lo,hi,color in [(max(leftwave,common_spx[0]),common_spx[1],'#34836b'),
                    (common_mrs[0],common_mrs[1],'#d6a132')]:
    fig.add_artist(Rectangle((fx(lo),low),fx(hi)-fx(lo),high-low,
                            transform=fig.transFigure,facecolor=color,alpha=.15,edgecolor='none'))
    fig.add_artist(Rectangle((fx(lo),high+.009),fx(hi)-fx(lo),.012,
                            transform=fig.transFigure,facecolor=color,edgecolor='none'))
fig.text((fx(leftwave)+fx(common_spx[1]))/2,high+.030,'SPHEREx',ha='center',fontsize=8,color='#28715d',weight='bold')
fig.text((fx(common_mrs[0])+fx(common_mrs[1]))/2,high+.030,'MIRI / MRS',ha='center',fontsize=8,color='#946916',weight='bold')
for wave in [9.7,18.]:
    fig.add_artist(plt.Line2D([fx(wave),fx(wave)],[low,high-.003],transform=fig.transFigure,
                            color='#946916',lw=.8,ls=':',alpha=.85))
    fig.text(fx(wave),low+.022,f'{wave:g} µm',ha='center',fontsize=8,color='#79540c',
             bbox={'facecolor':'white','alpha':.88,'edgecolor':'none','pad':1.8})
# Label band positions; this does not assert detections in the nuclear example.
pah_bands=[6.2,7.7,8.6,11.3,12.7,17.]
for wave in pah_bands:
    fig.add_artist(plt.Line2D([fx(wave),fx(wave)],[high-.105,high-.078],
                            transform=fig.transFigure,color='#82549a',lw=.85))
    fig.text(fx(wave),high-.068,f'{wave:g}',ha='left',va='bottom',rotation=55,
             fontsize=7,color='#82549a')
fig.text(fx(5.55),high-.09,'PAH',ha='right',va='center',fontsize=7.4,color='#82549a')
fig.text(fx(10.8),low+.21,'NGC 7469',fontsize=9.5,weight='bold',color='#283542',
         bbox={'facecolor':'white','alpha':.9,'edgecolor':'none','pad':1.5})
fig.text(fx(10.8),low+.17,'Donnan et al. (2024)',fontsize=8,color='#283542',
         bbox={'facecolor':'white','alpha':.9,'edgecolor':'none','pad':1.5})
fig.savefig(OUT/'overlay.pdf',transparent=True);plt.close(fig)
dest.merge_page(PdfReader(OUT/'overlay.pdf').pages[0])
writer=PdfWriter();writer.add_page(dest);writer.write(OUT/'literature_with_our_leverage.pdf')
subprocess.run(['pdftoppm','-scale-to','1600','-png','-singlefile',str(OUT/'literature_with_our_leverage.pdf'),
                str(OUT/'literature_with_our_leverage')],check=True,stderr=subprocess.DEVNULL)
assert before=={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in protected}
report=dict(n_targets=len(sample),redshift_range=[zmin,zmax],common_mrs_rest_um=common_mrs,
            pah_band_positions_rest_um=pah_bands,pah_labels_mean_detections=False,
            nominal_common_spherex_rest_um=common_spx,axis_tick_max_residual_pixels=residual,
            protected_sha256=before,active_proposal_unchanged=True,
            interpretation='Wavelength access only. No flux, feature-strength, component-fraction, or source-locus prediction.',
            source='https://arxiv.org/abs/2402.17479',figure='Figure 8, NGC 7469 spectrum',
            spherex_limit='Nominal instrument access; not a claim of complete usable sample coverage.')
(OUT/'checks.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))

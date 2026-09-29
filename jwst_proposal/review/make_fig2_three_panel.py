"""Arrange three literature examples at equal height and native aspect ratios.

Panel C is the optical and infrared portion of Lyu & Rieke Figure 15.
It is not a lag forecast or a monitoring plan for our sample.
"""
import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/clagn-matplotlib')
from pathlib import Path
from copy import deepcopy
import hashlib,json,subprocess
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from pypdf import PdfReader,PdfWriter,Transformation
from pypdf.generic import RectangleObject, ContentStream, DecodedStreamObject, NameObject

BASE=Path(__file__).resolve().parents[1];HERE=BASE/'review'
OUT=HERE/'fig2_three_panel';OUT.mkdir(exist_ok=True)
BUILD=OUT/'build';BUILD.mkdir(exist_ok=True)
protected=[BASE/'proposal.tex',BASE/'fig2_miri.pdf',
           OUT/'panel_a/literature_with_our_leverage.pdf']
before={str(p.relative_to(BASE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in protected}

def insert(dest,original,crop,place):
    src=deepcopy(original)
    if '/Annots' in src:del src['/Annots']
    sw,sh=float(src.mediabox.width),float(src.mediabox.height)
    l,t,r,b=crop;x0,y0,x1,y1=l*sw,(1-b)*sh,r*sw,(1-t)*sh
    src.cropbox=RectangleObject((x0,y0,x1,y1));src.trimbox=RectangleObject((x0,y0,x1,y1))
    dw,dh=float(dest.mediabox.width),float(dest.mediabox.height)
    x,y,w,h=place;x*=dw;y*=dh;w*=dw;h*=dh
    scale=min(w/(x1-x0),h/(y1-y0))
    dest.merge_transformed_page(src,Transformation().translate(-x0,-y0).scale(scale)
                                .translate(x+(w-(x1-x0)*scale)/2,y+(h-(y1-y0)*scale)/2))

plt.rcParams.update({'font.family':'DejaVu Sans','pdf.fonttype':42})
# Build the spatial panel without an exterior galaxy label or title.
fig=plt.figure(figsize=(5,4.4))
fig.savefig(BUILD/'spatial_layout.pdf');plt.close(fig)
spatial=PdfReader(BUILD/'spatial_layout.pdf').pages[0]
spatial_source=PdfReader(HERE/'miri_literature/gonzalez_martin_2504.01103.pdf').pages[8]
def image_region(x0,y0,x1,y1):
    return ((42.060844+.425622*x0)/595.276,
            (53.161162+.425622*y0)/782.362,
            (42.060844+.425622*x1)/595.276,
            (53.161162+.425622*y1)/782.362)
insert(spatial,spatial_source,image_region(0,29,1017,560),(0,1.79/4.4,1,2.61/4.4))
for region,x in zip([(1022,4,1193,179),(1022,181,1193,354),(1022,356,1193,531)],
                    [.12,1.80,3.48]):
    insert(spatial,spatial_source,image_region(*region),(x/5,0,1.52/5,1.60/4.4))

left=PdfReader(OUT/'panel_a/literature_with_our_leverage.pdf').pages[0]
time_reader=PdfReader(HERE/'miri_literature/alternatives_20260929/lyu_2011.07638.pdf')
paper=deepcopy(time_reader.pages[18])
# Omit the two complete long-wavelength series, including labels and error
# bars, directly from the vector drawing. Preserve any model curves behind
# their labels; white masks would erase those curves as well.
form=paper['/Resources']['/XObject']['/Im15'].get_object()
stream=ContentStream(form,time_reader)
fill=stroke=(0.,0.,0.);stack=[];filtered=[];omitted={}
def longwave(color):
    return any(max(abs(a-b) for a,b in zip(color,target))<1e-6
               for target in [(1.,.64705882,0.),(0.,.50196078,0.)])
for args,op in stream.operations:
    if op==b'q':stack.append((fill,stroke))
    elif op==b'Q':fill,stroke=stack.pop()
    elif op==b'rg':fill=tuple(map(float,args))
    elif op==b'RG':stroke=tuple(map(float,args))
    elif op==b'g':fill=(float(args[0]),)*3
    elif op==b'G':stroke=(float(args[0]),)*3
    remove=((op in [b'S',b's'] and longwave(stroke)) or
            (op in [b'f',b'f*',b'F',b'Tj',b'TJ',b'Do'] and longwave(fill)))
    if remove:
        key=op.decode();omitted[key]=omitted.get(key,0)+1
        if op in [b'S',b's',b'f',b'f*',b'F']:filtered.append(([],b'n'))
    else:filtered.append((args,op))
assert omitted.get('S')==7 and omitted.get('Do')==9,omitted
stream.operations=filtered
replacement=DecodedStreamObject()
for key,value in form.items():
    if key not in ['/Filter','/DecodeParms','/Length']:replacement[key]=value
replacement.set_data(stream.get_data())
paper['/Resources']['/XObject'][NameObject('/Im15')]=replacement
# Trim exterior whitespace; retain labels, scale bars, and original aspect ratios.
crop=(.081,.077,.558,.273)
sources=[(left,(.035,.045,.94,.955)),(spatial,(0,0,1,1)),(paper,crop)]
common_height=4.4
widths=[common_height*float(p.mediabox.width)*(c[2]-c[0])/
        (float(p.mediabox.height)*(c[3]-c[1])) for p,c in sources]
gap=.28;margin=.08;bottom=.08
total_width=sum(widths)+2*gap+2*margin;total_height=5.0
starts=[margin,margin+widths[0]+gap,margin+widths[0]+gap+widths[1]+gap]
fig=plt.figure(figsize=(total_width,total_height))
for x,w,title in zip(starts,widths,['A  Separating hot and warm dust',
                                   'B  Isolating the central emission',
                                   'C  Reading the dust’s memory']):
    fig.text((x+w/2)/total_width,.956,title,ha='center',va='center',
             fontsize=12,weight='bold',color='#283542')
fig.savefig(BUILD/'layout.pdf');plt.close(fig)
dest=PdfReader(BUILD/'layout.pdf').pages[0]
calendar_crop_top=.0925
year_header_height=common_height*(calendar_crop_top-crop[1])/(crop[3]-crop[1])
for panel,((p,c),x,w) in enumerate(zip(sources,starts,widths)):
    if panel==2:
        # Remove the original calendar labels above the optical curve and
        # reserve the same height below the infrared curve, without stretching.
        insert(dest,p,(c[0],calendar_crop_top,c[2],c[3]),
               (x/total_width,(bottom+year_header_height)/total_height,
                w/total_width,(common_height-year_header_height)/total_height))
    else:
        insert(dest,p,c,(x/total_width,bottom/total_height,w/total_width,common_height/total_height))
# Place source identifications in empty regions within the published axes.
fig=plt.figure(figsize=(total_width,total_height));fig.patch.set_alpha(0)
# Clear only the old calendar-label strip, preserving the upper y-axis text.
strip_left=(starts[2]+widths[2]*(82/612-crop[0])/(crop[2]-crop[0]))/total_width
strip_height=common_height*(81.5/792-calendar_crop_top)/(crop[3]-crop[1])/total_height
fig.add_artist(Rectangle((strip_left,(bottom+common_height)/total_height-strip_height),
                        (starts[2]+widths[2])/total_width-strip_left,strip_height+.004,
                        transform=fig.transFigure,facecolor='white',edgecolor='none'))
# Replace the parameter dump in the empty region above the N-band curves.
# Coordinates are on the original PDF page, measured from its upper left.
def time_x(x):
    return (starts[2]+widths[2]*(x/612-crop[0])/(crop[2]-crop[0]))/total_width
def time_y(y):
    return (bottom+common_height-common_height*(y/792-calendar_crop_top)/(crop[3]-crop[1]))/total_height
fig.add_artist(Rectangle((time_x(92),time_y(151)),time_x(332)-time_x(92),
                        time_y(116)-time_y(151),transform=fig.transFigure,
                        facecolor='white',edgecolor='none'))
fig.text(time_x(96),time_y(124),'Fitted delay ≈ 8 yr   •   Smoothing window ≈ 14 yr',
         fontsize=12,color='#283542',va='center')
def label(panel,x,y,name,reference):
    xx=(starts[panel]+x*widths[panel])/total_width
    yy=(bottom+y*common_height)/total_height
    box={'facecolor':'white','alpha':.92,'edgecolor':'none','pad':1.3}
    fig.text(xx,yy,name,fontsize=10,weight='bold',color='#283542',bbox=box)
    fig.text(xx,yy-.037,reference,fontsize=9,color='#283542',bbox=box)
label(1,.27,.97,'NGC 6552','González-Martín et al. (2025)')
label(2,.40,.46+year_header_height/common_height,'NGC 4151','Lyu & Rieke (2021)')
# Published tick-label centers in PDF page coordinates (pdftotext -bbox).
# Reuse their exact horizontal locations for the shared bottom calendar axis.
year_centers=[115.014363,161.859649,208.717762,255.563048,302.421161]
for year,source_x in zip([1970,1980,1990,2000,2010],year_centers):
    xx=(starts[2]+widths[2]*(source_x/612-crop[0])/(crop[2]-crop[0]))/total_width
    fig.text(xx,(bottom+year_header_height-.025)/total_height,str(year),
             ha='center',va='top',fontsize=12,fontfamily='DejaVu Serif')
fig.text((starts[2]+widths[2]*(year_centers[2]/612-crop[0])/(crop[2]-crop[0]))/total_width,
         (bottom-.035)/total_height,'Year',ha='center',va='bottom',
         fontsize=12,fontfamily='DejaVu Serif')
fig.savefig(BUILD/'identifications.pdf',transparent=True);plt.close(fig)
dest.merge_page(PdfReader(BUILD/'identifications.pdf').pages[0])
writer=PdfWriter();writer.add_page(dest);writer.write(OUT/'fig2_three_panel_draft.pdf')
subprocess.run(['pdftoppm','-scale-to','2800','-png','-singlefile',str(OUT/'fig2_three_panel_draft.pdf'),
                str(OUT/'fig2_three_panel_draft')],check=True,stderr=subprocess.DEVNULL)
after={str(p.relative_to(BASE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in protected};assert before==after
(OUT/'checks.json').write_text(json.dumps(dict(protected_files_unchanged=True,sha256=after,
    layout='Three panels in one horizontal row with equal content heights and centered titles; source aspect ratios preserved.',
    panel_height_inches=common_height,panel_widths_inches=widths,
    panel_b_source='https://arxiv.org/abs/2504.01103',panel_b_source_version='14 April 2026',
    panel_c_source='https://arxiv.org/abs/2011.07638',panel_c_figure=15,panel_c_pdf_page=19,
    panel_c_crop=crop,panel_c_time_axis='Original observed calendar years for NGC 4151, relocated below both light curves at the original horizontal positions.',
    panel_c_content='Optical B-band, N-band data and fitted reverberation model; complete 20–24 and 34–37 micron series and labels omitted for clarity.',
    panel_c_omitted_vector_operations=omitted,
    panel_c_simplified_annotation='Fitted delay ≈ 8 yr; smoothing window ≈ 14 yr. Rounded from original 2914 and 5263 days; amplitude, constant offset, and chi-square omitted.',
    interpretation='Empirical time-response precedent; no lag values or calendar windows transferred to our sample.',
    proposal_pdf_compiled=False),indent=2)+'\n')
print(OUT/'fig2_three_panel_draft.pdf')

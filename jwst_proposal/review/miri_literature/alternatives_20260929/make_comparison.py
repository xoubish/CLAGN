"""Render selected vector PDF figure regions for a literature comparison.

Source plot content is unmodified. Crops omit surrounding article text;
the source pages and full captions remain in index.html.
"""
from pathlib import Path
from io import BytesIO
from copy import deepcopy
import hashlib
import subprocess
from pypdf import PdfReader, PdfWriter, Transformation
from pypdf.generic import RectangleObject
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
names = ['proposal.tex', 'fig2_miri.pdf', 'make_fig2_miri.py']
before = {n: hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in names}
# Crop coordinates are fractional page coordinates, measured from top left.
items = [
 ('01_history', '1. Histories and dust components', 'Lyu & Rieke (2021), Fig. 10',
  HERE/'lyu_2011.07638.pdf',12,(.20,.073,.80,.429),
  'https://arxiv.org/abs/2011.07638'),
 ('02_response', '2. Which earlier illumination contributes?', 'Almeyda et al. (2020), Fig. 12',
  HERE/'almeyda_2002.12823.pdf',16,(.13,.078,.82,.548),
  'https://arxiv.org/abs/2002.12823'),
 ('03_components', '3. From a spectrum to dust components', 'Donnan et al. (2024), Fig. 8: two upper-right examples',
  HERE.parent/'donnan_2402.17479.pdf',13,(.490,.065,.925,.335),
  'https://arxiv.org/abs/2402.17479'),
 ('04_silicates', '4. Two silicate features, one diagnostic plane', 'Sirocky et al. (2008), Fig. 7',
  HERE/'sirocky_0801.4776.pdf',34,(.130,.146,.837,.531),
  'https://arxiv.org/abs/0801.4776'),
 ('05_spatial', '5. Central and surrounding emission', 'Gonzalez-Martin et al., Fig. 4 (downloaded arXiv version)',
  HERE.parent/'gonzalez_martin_2504.01103.pdf',9,(.065,.060,.935,.378),
  'https://arxiv.org/abs/2504.01103'),
 ('06_pahs', '6. PAH ratios probe the dust environment', 'Garcia-Bernete et al. (2024), GATOS V, Fig. 5',
  HERE/'gatos_2409.05686.pdf',8,(.09,.06,.91,.334),
  'https://arxiv.org/abs/2409.05686'),
]

def panel(item, width=760):
    name,title,source,path,pageno,box,url = item
    pg=deepcopy(PdfReader(path).pages[pageno-1])
    pw,ph=float(pg.mediabox.width),float(pg.mediabox.height)
    left,top,right,bottom=box
    x0,y0,x1,y1=left*pw,(1-bottom)*ph,right*pw,(1-top)*ph
    pg.cropbox=RectangleObject((x0,y0,x1,y1))
    pg.trimbox=RectangleObject((x0,y0,x1,y1))
    scale=(width-32)/(x1-x0)
    height=(y1-y0)*scale+96
    buf=BytesIO(); fig=plt.figure(figsize=(width/72,height/72))
    fig.text(16/width,1-25/height,title,fontsize=17,weight='bold')
    fig.text(16/width,1-43/height,source,fontsize=10,url=url)
    fig.text(16/width,17/height,
             'Published example; not our targets or a forecast. Full caption in the source paper.',
             fontsize=9,color='#475563')
    fig.savefig(buf,format='pdf');plt.close(fig);buf.seek(0)
    out=PdfReader(buf).pages[0]
    out.merge_transformed_page(pg,Transformation().translate(-x0,-y0).scale(scale).translate(16,34))
    return out

writer=PdfWriter(); pages=[]
for item in items:
    pg=panel(item);pages.append(pg);writer.add_page(pg)
    one=PdfWriter();one.add_page(pg);one.write(HERE/(item[0]+'.pdf'))
    subprocess.run(['pdftoppm','-scale-to','1500','-png','-singlefile',
                    str(HERE/(item[0]+'.pdf')),str(HERE/item[0])],check=True,
                   stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
writer.write(HERE/'all_six_examples.pdf')
# One vector contact sheet, suitable for zooming or printing.
sheet=PdfWriter();page=sheet.add_blank_page(width=1520,height=1870)
for i,pg in enumerate(pages):
    cellw,cellh=760,620
    scale=min((cellw-24)/float(pg.mediabox.width),(cellh-24)/float(pg.mediabox.height))
    x=(i%2)*cellw+(cellw-float(pg.mediabox.width)*scale)/2
    y=1870-(i//2+1)*cellh+(cellh-float(pg.mediabox.height)*scale)/2
    page.merge_transformed_page(pg,Transformation().scale(scale).translate(x,y))
sheet.write(HERE/'comparison_sheet.pdf')
subprocess.run(['pdftoppm','-scale-to','2000','-png','-singlefile',str(HERE/'comparison_sheet.pdf'),
                str(HERE/'comparison_sheet')],check=True,stderr=subprocess.DEVNULL)
assert before=={n: hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in names}
print('Six figure previews, six-page figure collection, and comparison sheet ready. Proposal assets unchanged.')

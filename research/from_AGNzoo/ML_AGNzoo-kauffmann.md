---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.16.1
kernelspec:
  display_name: Python 3 (ipykernel)
  language: python
  name: python3
---

# Type2 AGNs line ratios


```{code-cell} ipython3
#!pip install -r requirements.txt
import sys
import os
import re
import time

import astropy.units as u
from astropy.table import Table
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.cm as cm
from matplotlib.colors import LinearSegmentedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt

import numpy as np
import pandas as pd
sys.path.append('code_src/')
from data_structures import MultiIndexDFObject
from ML_utils import unify_lc, unify_lc_gp,unify_lc_gp_parallel, stat_bands, autopct_format, combine_bands,\
mean_fractional_variation, normalize_mean_objects, normalize_max_objects, \
normalize_clipmax_objects, shuffle_datalabel, dtw_distance, stretch_small_values_arctan, translate_bitwise_sum_to_labels, update_bitsums
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
from collections import Counter,OrderedDict

from scipy import interpolate
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RationalQuadratic, RBF
from scipy.interpolate import interp1d
from tqdm import tqdm
import random
import astropy.io.fits as fits

import umap
from sompy import * #using the SOMPY package from https://github.com/sevamoo/SOMPY

import logging

# Get the root logger
logger = logging.getLogger()
logger.setLevel(logging.ERROR)

import warnings
warnings.filterwarnings('ignore')

plt.style.use('bmh')
from matplotlib.colors import ListedColormap

# Define the hex codes as a list
colors_hex = ['#005C8A', '#1D3557', '#798777', '#B6D7A8', '#FFE88C',
              '#FF9A76', '#FF4C29', '#C92C6D', '#6B0504']

# Convert the hex codes to RGB tuples
colors_rgb = [(int(color[1:3], 16)/255, int(color[3:5], 16)/255, int(color[5:7], 16)/255) for color in colors_hex]

# Create a continuous colormap from the RGB tuples
cmap1 = LinearSegmentedColormap.from_list('Picasso', colors_rgb, N=256)

colors = ["#3F51B5", "#003153", "#0047AB", "#40826D", "#50C878", "#FFEA00", "#CC7722", "#E34234", "#E30022", "#D68A59", "#8A360F", "#826644", "#5C6BC0", "#002D62", "#0056B3", "#529987", "#66D98B", "#FFED47", "#E69C53", "#F2552C", "#EB4D55", "#E6A875", "#9F5A33", "#9C7C5B", "#2E37FE", "#76D7EA", "#007BA7", "#2E8B57", "#ADDFAD", "#FFFF31"]

color4 = ['#3182bd','#6baed6','#9ecae1','#e6550d','#fd8d3c','#fdd0a2','#31a354','#a1d99b', '#c7e9c0', '#756bb1', '#bcbddc', '#dadaeb', '#969696', '#bdbdbd','#d9d9d9']
custom_cmap = LinearSegmentedColormap.from_list("custom_theme", colors)
```

***


## 1) Loading data
Here we load a parquet file of light curves generated using the multiband_lc notebook. One can build the sample from different sources in the literature and grab the data from archives of interes.

```{code-cell} ipython3
r=0
redshifts, o3lum,o3corr, bpt1,bpt2, rml50, rmu, con, d4n,hda, vdisp = [],[],[],[],[],[],[],[],[],[],[]
with open("data/agn.dat_dr4_release.v2", 'r') as file:
    for line in file:
        parts = line.split()  # Splits the line into parts
        redshifts.append(float(parts[5]))
        o3lum.append(float(parts[6]))
        o3corr.append(float(parts[7]))
        bpt1.append(float(parts[8]))
        bpt2.append(float(parts[9]))
        rml50.append(float(parts[10]))
        rmu.append(float(parts[11]))
        con.append(float(parts[12]))
        d4n.append(float(parts[13]))
        hda.append(float(parts[14]))
        vdisp.append(float(parts[15]))
        r+=1
redshifts, o3lum,o3corr, bpt1,bpt2, rml50, rmu, con, d4n,hda, vdisp = np.array(redshifts), np.array(o3lum),np.array(o3corr), np.array(bpt1),np.array(bpt2), np.array(rml50), np.array(rmu), np.array(con), np.array(d4n),np.array(hda), np.array(vdisp)

df_lc = pd.read_parquet('data/df_lc_kauffmann.parquet')
```

```{raw-cell}
bands_inlc = ['zg','zr','zi','W1','W2']
numobjs = len(df_lc.index.get_level_values('objectid')[:].unique())
#objects,dobjects,flabels,keeps,zlist = unify_lc(df_lc, redshifts,bands_inlc,xres=160,numplots=3,low_limit_size=50) #nearest neightbor linear interpolation
#objects,dobjects,flabels,keeps,zlist = unify_lc_gp(df_lc,redshifts,bands_inlc,xres=160,numplots=5,low_limit_size=10) #Gaussian process unification
print(numobjs)
#numobjs = 10
sample_objids = df_lc.index.get_level_values('objectid').unique()[:numobjs]
df_lc_small = df_lc.loc[sample_objids]
objects,dobjects,flabels,zlist,keeps = unify_lc_gp_parallel(df_lc_small,redshifts,bands_inlc=bands_inlc,xres=180)


# calculate some basic statistics with a sigmaclipping with width 5sigma
fvar, maxarray, meanarray = stat_bands(objects,dobjects,bands_inlc,sigmacl=5)

# combine different waveband into one array
dat_notnormal = combine_bands(objects,bands_inlc)

# Normalize the combinde array by mean brightness in a waveband after clipping outliers:
datm = normalize_clipmax_objects(dat_notnormal,meanarray,band = 0)

# shuffle data incase the ML routines are sensitive to order
data,fzr,p = shuffle_datalabel(datm,flabels)
fvar_arr,maximum_arr,average_arr = fvar[:,p],maxarray[:,p],meanarray[:,p]
redshift_shuffled = zlist[p]

labc = {}  # Initialize labc to hold indices of each unique label
for index, f in enumerate(fzr):
    lab = translate_bitwise_sum_to_labels(int(f))
    for label in lab:
        if label not in labc:
            labc[label] = []  # Initialize the list for this label if it's not already in labc
        labc[label].append(index)  # Append the current index to the list of indices for this label
        
# Assuming `data` is your numpy array
nan_rows = np.any(np.isnan(data), axis=1)
clean_data = data[~nan_rows]  # Rows without NaNs

redshifts1, o3lum1,o3corr1, bpt11,bpt21, rml501, rmu1, con1, d4n1,hda1, vdisp1= redshifts[keeps], o3lum[keeps],o3corr[keeps], bpt1[keeps],bpt2[keeps], rml50[keeps], rmu[keeps], con[keeps], d4n[keeps],hda[keeps], vdisp[keeps]
redshifts2, o3lum2,o3corr2, bpt12,bpt22, rml502, rmu2, con2, d4n2,hda2, vdisp2 = redshifts1[p], o3lum1[p],o3corr1[p], bpt11[p],bpt21[p], rml501[p], rmu1[p], con1[p], d4n1[p],hda1[p], vdisp1[p]
redshifts3, o3lum3,o3corr3, bpt13,bpt23, rml503, rmu3, con3, d4n3,hda3, vdisp3 = redshifts2[~nan_rows], o3lum2[~nan_rows],o3corr2[~nan_rows], bpt12[~nan_rows],bpt22[~nan_rows], rml502[~nan_rows], rmu2[~nan_rows], con2[~nan_rows], d4n2[~nan_rows],hda2[~nan_rows], vdisp2[~nan_rows]

fvar_arr1,average_arr1 = fvar_arr[:,~nan_rows],average_arr[:,~nan_rows]
np.savez('data/kauffit_all2',data=clean_data,fvar_arr1 = fvar_arr1, average_arr1 = average_arr1,redshifts3=redshifts3,o3lum3=o3lum3,o3corr3=o3corr3,bpt13=bpt13,bpt23=bpt23,rml503=rml503,rmu3=rmu3,con3=con3,d4n3=d4n3,hda3=hda3,vdisp3=vdisp3)
```

```{code-cell} ipython3
d = np.load('data/kauffit_all.npz',allow_pickle=True)
clean_data = d['data']
fvar_arr1,average_arr1 = d['fvar_arr1'],d['average_arr1']
redshifts3, o3lum3,o3corr3, bpt13,bpt23, rml503, rmu3, con3, d4n3,hda3, vdisp3 = d['redshifts3'], d['o3lum3'],d['o3corr3'], d['bpt13'],d['bpt23'], d['rml503'], d['rmu3'], d['con3'],d['d4n3'],d['hda3'], d['vdisp3']                                                                                   
print(np.min(redshifts3),np.mean(redshifts3),len(redshifts3))
```

```{code-cell} ipython3
mapper = umap.UMAP(n_neighbors=100,min_dist=0.99,metric=dtw_distance,random_state=2).fit(clean_data)
```

```{code-cell} ipython3
from mpl_toolkits.axes_grid1 import make_axes_locatable

plt.figure(figsize=(10,6),facecolor='white')
cmap1 = 'viridis'
markersize=15

ax1 = plt.subplot(2,3,1)
ax1.set_title(r'$\rm Mean\ brightness$')
thiscolor=np.log10(np.nansum(average_arr1,axis=0))
u = (thiscolor<2) & (thiscolor>=0)
cf = ax1.scatter(mapper.embedding_[u,0],mapper.embedding_[u,1],c = thiscolor[u],s=markersize,edgecolor='k',cmap=cmap1)
plt.axis('off')
divider = make_axes_locatable(ax1)
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(cf,cax=cax)

ax0 = plt.subplot(2,3,2)
ax0.set_title(r'$\rm Mean\ Fractional\ Variation$')
thiscolor=stretch_small_values_arctan(np.nansum(fvar_arr1,axis=0)/np.count_nonzero(~np.isnan(fvar_arr1), axis=0),factor=30)
u = (thiscolor<2.0) & (thiscolor>=0.)
cf = ax0.scatter(mapper.embedding_[u,0],mapper.embedding_[u,1],c = thiscolor[u],s=markersize,edgecolor='k',cmap=cmap1)
plt.axis('off')
divider = make_axes_locatable(ax0)
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(cf,cax=cax)

ax0 = plt.subplot(2,3,3)

thiscolor=rml503
u = (thiscolor<12.5) & (thiscolor>9.5)
ax0.set_title(r'$\rm Log (M_{* 50}/M_{\odot})$',color='k')
cf = ax0.scatter(mapper.embedding_[u,0],mapper.embedding_[u,1],c = thiscolor[u],s=markersize,edgecolor='k',cmap=cmap1)
plt.axis('off')
divider = make_axes_locatable(ax0)
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(cf,cax=cax)


ax0 = plt.subplot(2,3,5)
thiscolor=o3corr3
u = (thiscolor<9) & (thiscolor>5.5)
ax0.set_title(r'$\rm Log L[OIII] (L_{\odot})$',color='k')
cf = ax0.scatter(mapper.embedding_[u,0],mapper.embedding_[u,1],c = thiscolor[u],s=markersize,edgecolor='k',cmap=cmap1)
plt.axis('off')
divider = make_axes_locatable(ax0)
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(cf,cax=cax)

ax0 = plt.subplot(2,3,6)
thiscolor=hda3
u = (thiscolor<10) & (thiscolor>-2)
ax0.set_title(r'$\rm H\delta_{A}$',color='k')
cf = ax0.scatter(mapper.embedding_[u,0],mapper.embedding_[u,1],c = thiscolor[u],s=markersize,edgecolor='k',cmap=cmap1)
plt.axis('off')
divider = make_axes_locatable(ax0)
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(cf,cax=cax)

ax0 = plt.subplot(2,3,4)

thiscolor=d4n3
u = (thiscolor<2.5) & (thiscolor>1)
ax0.set_title(r'$\rm D_{n}4000$',color='k')
cf = ax0.scatter(mapper.embedding_[u,0],mapper.embedding_[u,1],c = thiscolor[u],s=markersize,edgecolor='k',cmap=cmap1)
plt.axis('off')
divider = make_axes_locatable(ax0)
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(cf,cax=cax)


plt.tight_layout()
plt.savefig('output/kauffmann-umap.png')
```

```{code-cell} ipython3
plt.figure(figsize=(12,5))
markersize=20
mapper = umap.UMAP(n_neighbors=100,min_dist=0.99,metric='manhattan',random_state=18).fit(clean_data)

ax1 = plt.subplot(1,2,1)
ax1.set_title(r'mean brightness',size=20)
#cf = ax1.scatter(mapper.embedding_[:,0],mapper.embedding_[:,1],s=markersize,c=np.log10(np.nansum(average_arr1,axis=0)),edgecolor='k',cmap='cividis')
thiscolor=np.log10(np.nansum(average_arr1,axis=0))
u = (thiscolor<2.) & (thiscolor>=0.)
cf = ax1.scatter(mapper.embedding_[u,0],mapper.embedding_[u,1],c = thiscolor[u],s=markersize,edgecolor='k',cmap='cividis')
plt.axis('off')
divider = make_axes_locatable(ax1)
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(cf,cax=cax)

print(np.shape(fvar_arr1))
ax0 = plt.subplot(1,2,2)
ax0.set_title(r'mean fractional variation',size=20)
thiscolor=stretch_small_values_arctan(np.nansum(fvar_arr1,axis=0)/np.count_nonzero(~np.isnan(fvar_arr1), axis=0),factor=30)
u = (thiscolor<2.0) & (thiscolor>=0.)
cf = ax0.scatter(mapper.embedding_[u,0],mapper.embedding_[u,1],c = thiscolor[u],s=markersize,edgecolor='k',cmap='cividis')
plt.axis('off')
divider = make_axes_locatable(ax0)
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(cf,cax=cax)
plt.tight_layout()
#plt.savefig('umap-ztf.png')
```

```{code-cell} ipython3
from mpl_toolkits.axes_grid1 import make_axes_locatable

plt.figure(figsize=(12,5),facecolor='white')
plt.subplot(1,3,1)
thiscolor=o3lum3
u = (thiscolor<9) & (thiscolor>5)
plt.title('OIII Luminosity',color='k')
plt.scatter(mapper.embedding_[u,0],mapper.embedding_[u,1],s=8,c = thiscolor[u],edgecolor='k',cmap=cmap1)
plt.axis('off')
divider = make_axes_locatable(plt.gca())
cax = divider.append_axes("right", size="5%", pad=0.05)
cbar = plt.colorbar(cax=cax)
cbar.ax.yaxis.set_tick_params(color='k')
o = plt.setp(plt.getp(cbar.ax.axes, 'yticklabels'), color='k')

plt.subplot(1,2,2)
thiscolor=o3corr3
u = (thiscolor<9) & (thiscolor>5)
plt.title('OIII Luminosity corrected',color='k')
plt.scatter(mapper.embedding_[u,0],mapper.embedding_[u,1],s=8,c = thiscolor[u],edgecolors='k',cmap=cmap1)
plt.axis('off')
divider = make_axes_locatable(plt.gca())
cax = divider.append_axes("right", size="5%", pad=0.05)
cbar = plt.colorbar(cax=cax)
cbar.ax.yaxis.set_tick_params(color='k')
o = plt.setp(plt.getp(cbar.ax.axes, 'yticklabels'), color='k')
```

```{code-cell} ipython3
plt.figure(figsize=(12,5),facecolor='black')
plt.subplot(1,2,1)
thiscolor=bpt13
u = (thiscolor<1.5) & (thiscolor>-1)
plt.title('BPT1',color='w')
plt.scatter(mapper.embedding_[u,0],mapper.embedding_[u,1],s=8,c = thiscolor[u],edgecolors='k',cmap=cmap1)
plt.axis('off')
divider = make_axes_locatable(plt.gca())
cax = divider.append_axes("right", size="5%", pad=0.05)
cbar = plt.colorbar(cax=cax)
cbar.ax.yaxis.set_tick_params(color='white')
o = plt.setp(plt.getp(cbar.ax.axes, 'yticklabels'), color='white')

plt.subplot(1,2,2)
thiscolor=bpt23
u = (thiscolor<1.5) & (thiscolor>-1)
plt.title('BPT2',color='w')
plt.scatter(mapper.embedding_[u,0],mapper.embedding_[u,1],s=8,c = thiscolor[u],edgecolors='k',cmap=cmap1)
plt.axis('off')
divider = make_axes_locatable(plt.gca())
cax = divider.append_axes("right", size="5%", pad=0.05)
cbar = plt.colorbar(cax=cax)
cbar.ax.yaxis.set_tick_params(color='white')
o = plt.setp(plt.getp(cbar.ax.axes, 'yticklabels'), color='white')
```

```{code-cell} ipython3
plt.figure(figsize=(20,4))

ykewley = [0.61/(x - 0.47) + 1.19 if x < 0.47 else -3 for x in bpt23]
ykauffmann = [0.61/(x - 0.05) + 1.3 if x < 0.05 else -3 for x in bpt23]

usf = (bpt13<ykauffmann)&(bpt13>-2)
ucomp = (bpt13<=ykewley)&(bpt13>ykauffmann)
uagn = (bpt13>=ykewley)&(bpt13<=2)

plt.subplot(1,4,1)
x = np.linspace(-1,0.46,10)
ykewley = 0.61/(x - 0.47) + 1.19
plt.plot(x,ykewley,'r--')
x = np.linspace(-1,0.05,10)
ykauffmann = 0.61/(x - 0.05) + 1.3
plt.plot(x,ykauffmann,'r--')
plt.scatter(bpt23[usf],bpt13[usf],c ='blue',marker='o')
plt.scatter(bpt23[ucomp],bpt13[ucomp],c ='green',marker='o')
plt.scatter(bpt23[uagn],bpt13[uagn],c ='purple',marker='o')
plt.xlim([-1,1])
plt.ylim([-1,1.5])
plt.xlabel(r'$\log([\rm NII]\lambda 6584/\rm H\alpha)$',fontsize=14)
plt.ylabel(r'$\log([\rm OIII]\lambda 5007/\rm H\beta)$',fontsize=14)


plt.subplot(1,4,2)
plt.title('almost SF')
plt.scatter(mapper.embedding_[:,0],mapper.embedding_[:,1],s=8,c = 'gray',edgecolors='gray',alpha=0.5)
plt.scatter(mapper.embedding_[usf,0],mapper.embedding_[usf,1],s=8,c = 'blue',edgecolors='k')
plt.axis('off')

plt.subplot(1,4,3)
plt.title('Kauffmann')
plt.scatter(mapper.embedding_[:,0],mapper.embedding_[:,1],s=8,c = 'gray',edgecolors='gray',alpha=0.5)
plt.scatter(mapper.embedding_[ucomp,0],mapper.embedding_[ucomp,1],s=8,c = 'g',edgecolors='k')
plt.axis('off')

plt.subplot(1,4,4)
plt.title('Kewley')
plt.scatter(mapper.embedding_[:,0],mapper.embedding_[:,1],s=8,c = 'gray',edgecolors='gray',alpha=0.5)
plt.scatter(mapper.embedding_[uagn,0],mapper.embedding_[uagn,1],s=8,c = 'purple',edgecolors='k')
plt.axis('off')
```

```{code-cell} ipython3
plt.figure(figsize=(8,8))

x = np.linspace(-1,0.46,10)
ykewley = 0.61/(x - 0.47) + 1.19
plt.plot(x,ykewley,'r--')
x = np.linspace(-1,0.05,10)
ykauffmann = 0.61/(x - 0.05) + 1.3
plt.plot(x,ykauffmann,'r--')
thiscolor=d4n3
u = (thiscolor<3) & (thiscolor>0)
plt.scatter(bpt23[u],bpt13[u],c =thiscolor[u],marker='o')

plt.xlim([-1,1])
plt.ylim([-1,1.5])
plt.xlabel(r'$\log([\rm NII]\lambda 6584/\rm H\alpha)$',fontsize=14)
plt.ylabel(r'$\log([\rm OIII]\lambda 5007/\rm H\beta)$',fontsize=14)
```

```{code-cell} ipython3
plt.figure(figsize=(16,5),facecolor='black')
plt.subplot(1,3,1)
thiscolor=rml503
u = (thiscolor<13) & (thiscolor>9)
plt.title('Mass 50',color='w')
plt.scatter(mapper.embedding_[u,0],mapper.embedding_[u,1],s=8,c = thiscolor[u],edgecolors='k',cmap=cmap1)
plt.axis('off')
divider = make_axes_locatable(plt.gca())
cax = divider.append_axes("right", size="5%", pad=0.05)
cbar = plt.colorbar(cax=cax)
cbar.ax.yaxis.set_tick_params(color='white')
o = plt.setp(plt.getp(cbar.ax.axes, 'yticklabels'), color='white')

plt.subplot(1,3,2)
thiscolor=rmu3
u = (thiscolor<10) & (thiscolor>7)
plt.title('Mass Surface density',color='w')
plt.scatter(mapper.embedding_[u,0],mapper.embedding_[u,1],s=8,c = thiscolor[u],edgecolors='k',cmap=cmap1)
plt.axis('off')
divider = make_axes_locatable(plt.gca())
cax = divider.append_axes("right", size="5%", pad=0.05)
cbar = plt.colorbar(cax=cax)
cbar.ax.yaxis.set_tick_params(color='white')
o = plt.setp(plt.getp(cbar.ax.axes, 'yticklabels'), color='white')

plt.subplot(1,3,3)
thiscolor=con3
u = (thiscolor<4) & (thiscolor>1)
plt.title('Concentration',color='w')
plt.scatter(mapper.embedding_[u,0],mapper.embedding_[u,1],s=8,c = thiscolor[u],edgecolors='k',cmap=cmap1)
plt.axis('off')
divider = make_axes_locatable(plt.gca())
cax = divider.append_axes("right", size="5%", pad=0.05)
cbar = plt.colorbar(cax=cax)
cbar.ax.yaxis.set_tick_params(color='white')
o = plt.setp(plt.getp(cbar.ax.axes, 'yticklabels'), color='white')
```

```{code-cell} ipython3
plt.figure(figsize=(16,4),facecolor='black')
plt.subplot(1,3,1)
thiscolor=d4n3
u = (thiscolor<2.5) & (thiscolor>1)
plt.title('D4000?',color='w')
plt.scatter(mapper.embedding_[u,0],mapper.embedding_[u,1],s=8,c = thiscolor[u],edgecolors='k',cmap=cmap1)
plt.axis('off')
divider = make_axes_locatable(plt.gca())
cax = divider.append_axes("right", size="5%", pad=0.05)
cbar = plt.colorbar(cax=cax)
cbar.ax.yaxis.set_tick_params(color='white')
o = plt.setp(plt.getp(cbar.ax.axes, 'yticklabels'), color='white')

plt.subplot(1,3,2)
thiscolor=hda3
u = (thiscolor<10) & (thiscolor>-10)
plt.title('HDA',color='w')
plt.scatter(mapper.embedding_[u,0],mapper.embedding_[u,1],s=8,c = thiscolor[u],edgecolors='k',cmap=cmap1)
plt.axis('off')
divider = make_axes_locatable(plt.gca())
cax = divider.append_axes("right", size="5%", pad=0.05)
cbar = plt.colorbar(cax=cax)
cbar.ax.yaxis.set_tick_params(color='white')
o = plt.setp(plt.getp(cbar.ax.axes, 'yticklabels'), color='white')

plt.subplot(1,3,3)
thiscolor=vdisp3
u = (thiscolor<400) & (thiscolor>0)
plt.title('Vdispersion',color='w')
plt.scatter(mapper.embedding_[u,0],mapper.embedding_[u,1],s=8,c = thiscolor[u],edgecolors='k',cmap=cmap1)
plt.axis('off')
divider = make_axes_locatable(plt.gca())
cax = divider.append_axes("right", size="5%", pad=0.05)
cbar = plt.colorbar(cax=cax)
cbar.ax.yaxis.set_tick_params(color='white')
o = plt.setp(plt.getp(cbar.ax.axes, 'yticklabels'), color='white')
```

```{code-cell} ipython3
plt.figure(figsize=(16,4),facecolor='black')
plt.subplot(1,3,1)
thiscolor=rml503
u = (thiscolor<13) & (thiscolor>8.5)
plt.title('Mass',color='w')
plt.scatter(mapper.embedding_[u,0],mapper.embedding_[u,1],s=8,c = thiscolor[u],edgecolors='k',cmap=cmap1)
plt.axis('off')
divider = make_axes_locatable(plt.gca())
cax = divider.append_axes("right", size="5%", pad=0.05)
cbar = plt.colorbar(cax=cax)
cbar.ax.yaxis.set_tick_params(color='white')
o = plt.setp(plt.getp(cbar.ax.axes, 'yticklabels'), color='white')

plt.subplot(1,3,2)
thiscolor=d4n3
u = (thiscolor<2.5) & (thiscolor>1)
plt.title('D4000',color='w')
plt.scatter(mapper.embedding_[u,0],mapper.embedding_[u,1],s=8,c = thiscolor[u],edgecolors='k',cmap=cmap1)
plt.axis('off')
divider = make_axes_locatable(plt.gca())
cax = divider.append_axes("right", size="5%", pad=0.05)
cbar = plt.colorbar(cax=cax)
cbar.ax.yaxis.set_tick_params(color='white')
o = plt.setp(plt.getp(cbar.ax.axes, 'yticklabels'), color='white')

plt.subplot(1,3,3)
thiscolor=vdisp3
u = (thiscolor<400) & (thiscolor>0)
plt.title('Vdispersion',color='w')
plt.scatter(mapper.embedding_[u,0],mapper.embedding_[u,1],s=8,c = thiscolor[u],edgecolors='k',cmap=cmap1)
plt.axis('off')
divider = make_axes_locatable(plt.gca())
cax = divider.append_axes("right", size="5%", pad=0.05)
cbar = plt.colorbar(cax=cax)
cbar.ax.yaxis.set_tick_params(color='white')
o = plt.setp(plt.getp(cbar.ax.axes, 'yticklabels'), color='white')
```

```{code-cell} ipython3
msz0,msz1 = 35,35
sm = sompy.SOMFactory.build(clean_data, mapsize=[msz0,msz1], mapshape='planar', lattice='rect', initialization='pca')
sm.train(n_job=4, shared_memory = 'no')
```

```{code-cell} ipython3
a=sm.bmu_ind_to_xy(sm.project_data(clean_data))
x,y=np.zeros(len(a)),np.zeros(len(a))
k=0
for i in a:
    x[k]=i[0]
    y[k]=i[1]
    k+=1
BPT1_r,BPT2_r=np.zeros([msz0,msz1]),np.zeros([msz0,msz1])
for i in range(msz0):
    for j in range(msz1):
        unja=(x==i)&(y==j)
        BPT1_r[i,j]=(np.nanmedian(hda3[unja]))
        BPT2_r[i,j]=(np.nanmedian(rml503[unja]))

plt.figure(figsize=(9,5))
plt.subplot(1,2,1)
cf=plt.imshow(BPT1_r,vmin=-5,vmax=5,origin='lower',cmap='viridis')
plt.axis('off')
plt.subplot(1,2,2)
cf=plt.imshow(BPT2_r,vmin=10,vmax=11.5,origin='lower',cmap='viridis')
plt.axis('off')
plt.tight_layout()
```

```{code-cell} ipython3
plt.figure(figsize=(12,10))
markersize=200

mapper = umap.UMAP(n_neighbors=50,min_dist=0.9,metric='euclidean',random_state=20).fit(data)
ax0 = plt.subplot(2,2,1)
ax0.set_title(r'Euclidean Distance, min_d=0.9, n_neighbors=50',size=12)
for label, indices in (labc.items()):
     cf = ax0.scatter(mapper.embedding_[indices,0],mapper.embedding_[indices,1],s=80,alpha=0.5,edgecolor='gray',label=label)
plt.axis('off')

mapper = umap.UMAP(n_neighbors=50,min_dist=0.9,metric='manhattan',random_state=20).fit(data)
ax0 = plt.subplot(2,2,2)
ax0.set_title(r'Manhattan Distance, min_d=0.9, n_neighbors=50',size=12)
for label, indices in (labc.items()):
     cf = ax0.scatter(mapper.embedding_[indices,0],mapper.embedding_[indices,1],s=80,alpha=0.5,edgecolor='gray',label=label)
plt.axis('off')


mapperg = umap.UMAP(n_neighbors=50,min_dist=0.9,metric=dtw_distance,random_state=20).fit(data) #this distance takes long
ax2 = plt.subplot(2,2,3)
ax2.set_title(r'DTW Distance, min_d=0.9,n_neighbors=50',size=12)
for label, indices in (labc.items()):
     cf = ax2.scatter(mapper.embedding_[indices,0],mapper.embedding_[indices,1],s=80,alpha=0.5,edgecolor='gray',label=label)
plt.axis('off')


mapper = umap.UMAP(n_neighbors=50,min_dist=0.1,metric='manhattan',random_state=20).fit(data)
ax0 = plt.subplot(2,2,4)
ax0.set_title(r'Manhattan Distance, min_d=0.1, n_neighbors=50',size=12)
for label, indices in (labc.items()):
     cf = ax0.scatter(mapper.embedding_[indices,0],mapper.embedding_[indices,1],s=80,alpha=0.5,edgecolor='gray',label=label)
plt.legend(fontsize=12)
plt.axis('off')
```

```{code-cell} ipython3

```

```{code-cell} ipython3

```

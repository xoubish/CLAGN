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

# How do AGNs selected with different techniques compare?

In this notebook I added BOSS starforming galaxies to the Kauffmann sample with no emission line ratio information except knowing that they are in the starforming part. To see if we can separate thing on bpt using the umap trained on light curves, which seems not the case.

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

```{raw-cell}
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

bands_inlc = ['zg']
numobjs = 5000#len(df_lc.index.get_level_values('objectid')[:].unique())
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
np.savez('data/kauffw2',data=clean_data,fvar_arr1 = fvar_arr1, average_arr1 = average_arr1,redshifts3=redshifts3,o3lum3=o3lum3,o3corr3=o3corr3,bpt13=bpt13,bpt23=bpt23,rml503=rml503,rmu3=rmu3,con3=con3,d4n3=d4n3,hda3=hda3,vdisp3=vdisp3)
```

```{raw-cell}
df_lc = pd.read_parquet('data/df_lc_boss_sf.parquet')

samp = pd.read_csv('data/BOSS-SF.csv')
redshifts = samp['redshift']#[objids]

bands_inlc = ['zg']
numobjs = 5000#len(df_lc.index.get_level_values('objectid')[:].unique())
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

redshifts3 = redshift_shuffled[~nan_rows]
fvar_arr1,average_arr1 = fvar_arr[:,~nan_rows],average_arr[:,~nan_rows]
np.savez('data/bossw2',data=clean_data,fvar_arr1 = fvar_arr1, average_arr1 = average_arr1,redshifts3=redshifts3)
```

```{code-cell} ipython3
d1 = np.load('data/kauffit_all.npz', allow_pickle=True)
d2 = np.load('data/bossit_all.npz', allow_pickle=True)

bptsf = np.zeros_like(d2['redshifts3'])-0.5

clean_data1 = np.concatenate([d1['data'], d2['data']], axis=0)
fvar_arr11 = np.concatenate([d1['fvar_arr1'], d2['fvar_arr1']], axis=1)
average_arr11 = np.concatenate([d1['average_arr1'], d2['average_arr1']], axis=1)
redshifts31 = np.concatenate([d1['redshifts3'], d2['redshifts3']], axis=0)
bpt131 = np.concatenate([d1['bpt13'], bptsf], axis=0)
bpt231 = np.concatenate([d1['bpt23'], bptsf], axis=0)

u = (redshifts31>0.01) & (redshifts31<0.2)
clean_data = clean_data1[u]
fvar_arr1 = fvar_arr11[:,u]
average_arr1 = average_arr11[:,u]
redshifts3 = redshifts31[u]
bpt13 = bpt131[u]
bpt23 = bpt231[u]
```

```{code-cell} ipython3
plt.figure(figsize=(16,5))
markersize=30
mapper = umap.UMAP(n_neighbors=20,min_dist=0.99,metric='manhattan',random_state=18).fit(clean_data)

ax1 = plt.subplot(1,3,1)
ax1.set_title(r'redshift',size=20)
#cf = ax1.scatter(mapper.embedding_[:,0],mapper.embedding_[:,1],s=markersize,c=np.log10(np.nansum(average_arr1,axis=0)),edgecolor='k',cmap='cividis')
thiscolor=redshifts3
u = (thiscolor<0.4) & (thiscolor>=0.)
cf = ax1.scatter(mapper.embedding_[u,0],mapper.embedding_[u,1],c = thiscolor[u],s=markersize,edgecolor='k',cmap='cividis')
plt.axis('off')
divider = make_axes_locatable(ax1)
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(cf,cax=cax)


ax1 = plt.subplot(1,3,2)
ax1.set_title(r'mean brightness',size=20)
#cf = ax1.scatter(mapper.embedding_[:,0],mapper.embedding_[:,1],s=markersize,c=np.log10(np.nansum(average_arr1,axis=0)),edgecolor='k',cmap='cividis')
thiscolor=np.log10(np.nansum(average_arr1,axis=0))
u = (thiscolor<3.) & (thiscolor>=0.)
cf = ax1.scatter(mapper.embedding_[u,0],mapper.embedding_[u,1],c = thiscolor[u],s=markersize,edgecolor='k',cmap='cividis')
plt.axis('off')
divider = make_axes_locatable(ax1)
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(cf,cax=cax)


ax0 = plt.subplot(1,3,3)
ax0.set_title(r'mean fractional variation',size=20)
thiscolor=stretch_small_values_arctan(np.nansum(fvar_arr1,axis=0)/np.count_nonzero(~np.isnan(fvar_arr1), axis=0),factor=30)
u = (thiscolor<1.2) & (thiscolor>=0.)
cf = ax0.scatter(mapper.embedding_[u,0],mapper.embedding_[u,1],c = thiscolor[u],s=markersize,edgecolor='k',cmap='cividis')
plt.axis('off')
divider = make_axes_locatable(ax0)
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(cf,cax=cax)
plt.tight_layout()
#plt.savefig('umap-ztf.png')
```

```{code-cell} ipython3
plt.figure(figsize=(20,4))

ykewley = [0.61/(x - 0.47) + 1.19 if x < 0.47 else -3 for x in bpt23]
ykauffmann = [0.61/(x - 0.05) + 1.3 if x < 0.05 else -3 for x in bpt23]

usf = (bpt131==-0.5)&(bpt231==-0.5)#(bpt13<ykauffmann)&(bpt13>-2)
ucomp = (bpt13<=ykewley)&(bpt13>ykauffmann)
uagn = (bpt13>=ykewley)&(bpt13<=2)

plt.subplot(1,4,1)
x = np.linspace(-1,0.46,30)
ykewley = 0.61/(x - 0.47) + 1.19
plt.plot(x,ykewley,'r-.',label='Keweley 2003')
x = np.linspace(-1,0.05,30)
ykauffmann = 0.61/(x - 0.05) + 1.3
plt.plot(x,ykauffmann,'r--',label='Kauffmann 2003')

thiscolor=stretch_small_values_arctan(np.nansum(fvar_arr1,axis=0)/np.count_nonzero(~np.isnan(fvar_arr1), axis=0),factor=30)
u = (thiscolor<1.2) & (thiscolor>=0.)
plt.scatter(bpt23[u],bpt13[u],c =thiscolor[u],marker='o',alpha=0.7)
plt.colorbar()
plt.legend()
plt.xlim([-1,1])
plt.ylim([-1,1.5])
plt.xlabel(r'$\log([\rm NII]\lambda 6584/\rm H\alpha)$',fontsize=14)
plt.ylabel(r'$\log([\rm OIII]\lambda 5007/\rm H\beta)$',fontsize=14)


plt.subplot(1,4,2)
plt.title('SF')
plt.scatter(mapper.embedding_[:,0],mapper.embedding_[:,1],s=38,c = 'gray',edgecolors='gray',alpha=0.5)
plt.scatter(mapper.embedding_[usf,0],mapper.embedding_[usf,1],s=8,c = 'pink',edgecolors='k')
plt.axis('off')

plt.subplot(1,4,3)
plt.title('Kauffmann')
plt.scatter(mapper.embedding_[:,0],mapper.embedding_[:,1],s=38,c = 'gray',edgecolors='gray',alpha=0.5)
plt.scatter(mapper.embedding_[ucomp,0],mapper.embedding_[ucomp,1],s=8,c = 'pink',edgecolors='k')
plt.axis('off')

plt.subplot(1,4,4)
plt.title('Kewley')
plt.scatter(mapper.embedding_[:,0],mapper.embedding_[:,1],s=38,c = 'gray',edgecolors='gray',alpha=0.5)
plt.scatter(mapper.embedding_[uagn,0],mapper.embedding_[uagn,1],s=8,c = 'pink',edgecolors='k')
plt.axis('off')
```

```{code-cell} ipython3
msz0,msz1 = 20,20
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
    
fr=np.zeros([msz0,msz1])
thiscolor=stretch_small_values_arctan(np.nansum(fvar_arr1,axis=0)/np.count_nonzero(~np.isnan(fvar_arr1), axis=0),factor=30)
for i in range(msz0):
    for j in range(msz1):
        unja=(x==i)&(y==j)
        fr[i,j]=(np.nanmedian(thiscolor[unja]))

aagn=sm.bmu_ind_to_xy(sm.project_data(clean_data[uagn]))
xagn,yagn=np.zeros(len(aagn)),np.zeros(len(aagn))
k=0
for i in aagn:
    xagn[k]=i[0]
    yagn[k]=i[1]
    k+=1
    
acomp=sm.bmu_ind_to_xy(sm.project_data(clean_data[ucomp]))
xcomp,ycomp=np.zeros(len(acomp)),np.zeros(len(acomp))
k=0
for i in acomp:
    xcomp[k]=i[0]
    ycomp[k]=i[1]
    k+=1

asf=sm.bmu_ind_to_xy(sm.project_data(clean_data[usf]))
xsf,ysf=np.zeros(len(asf)),np.zeros(len(asf))
k=0
for i in asf:
    xsf[k]=i[0]
    ysf[k]=i[1]
    k+=1


plt.figure(figsize=(12,5))
plt.subplot(1,3,3)
plt.title('BPT AGN')
cf=plt.imshow(fr,vmin=0,vmax=1.5,origin='lower',cmap='viridis')
plt.axis('off')
plt.scatter(yagn,xagn,marker='x',c='k',alpha=0.5)
plt.subplot(1,3,2)
plt.title('BPT Composite')
cf=plt.imshow(fr,vmin=0,vmax=1.5,origin='lower',cmap='viridis')
plt.axis('off')
plt.scatter(ycomp,xcomp,marker='x',c='k',alpha=0.5)
plt.subplot(1,3,1)
plt.title('BPT SF')
cf=plt.imshow(fr,vmin=0,vmax=1.5,origin='lower',cmap='viridis')
plt.axis('off')
plt.scatter(ysf,xsf,marker='x',c='k',alpha=0.5)

plt.tight_layout()
```

```{code-cell} ipython3

```

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

By the IPAC Science Platform Team, last edit: April 12th, 2024

***


## Goal
Writing these up



## Introduction

Active Galactic Nuclei (AGNs), some of the most powerful sources in the universe, emit a broad range of electromagnetic radiation, from radio waves to gamma rays. Consequently, there is a wide variety of AGN labels depending on the identification/selection scheme and the presence or absence of certain emissions (e.g., Radio loud/quiet, Quasars, Blazars, Seiferts, Changing looks). According to the unified model, this zoo of labels we see depend on a limited number of parameters, namely the viewing angle, the accretion rate, presence or lack of jets, and perhaps the properties of the host/environment (e.g., [Padovani et al. 2017](https://arxiv.org/pdf/1707.07134.pdf)). Here, we collect archival temporal data and labels from the literature to compare how some of these different labels/selection schemes compare.

We use manifold learning and dimensionality reduction to learn the distribution of AGN lightcurves observed with different facilities. We mostly focus on UMAP ([Uniform Manifold Approximation and Projection, McInnes 2020](https://arxiv.org/pdf/1802.03426.pdf)) but also show SOM ([Self organizing Map, Kohonen 1990](https://ieeexplore.ieee.org/document/58325)) examples. The reduced 2D projections from these two unsupervised ML techniques reveal similarities and overlaps of different selection techniques and coloring the projections with various statistical physical properties (e.g., mean brightness, fractional lightcurve variation) is informative of correlations of the selections technique with physics such as AGN variability. Using different parts of the EM in training (or in building the initial higher dimensional manifold) demonstrates how much information if any is in that part of the data for each labeling scheme, for example whether with ZTF optical light curves alone, we can identify sources with variability in WISE near IR bands. These techniques also have a potential for identifying targets of a specific class or characteristic for future follow up observations.

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

import umap
from sompy import * #using the SOMPY package from https://github.com/sevamoo/SOMPY

import logging

# Get the root logger
logger = logging.getLogger()
logger.setLevel(logging.ERROR)

import warnings
warnings.filterwarnings('ignore')

#plt.style.use('bmh')

colors = ["#3F51B5", "#003153", "#0047AB", "#40826D", "#50C878", "#FFEA00", "#CC7722", "#E34234", "#E30022", "#D68A59", "#8A360F", "#826644", "#5C6BC0", "#002D62", "#0056B3", "#529987", "#66D98B", "#FFED47", "#E69C53", "#F2552C", "#EB4D55", "#E6A875", "#9F5A33", "#9C7C5B", "#2E37FE", "#76D7EA", "#007BA7", "#2E8B57", "#ADDFAD", "#FFFF31"]
colors2 = [
    "#3F51B5",  # Ultramarine Blue
    "#003153",  # Prussian Blue
    "#0047AB",  # Cobalt Blue
    "#40826D",  # Viridian Green
    "#50C878",  # Emerald Green
    "#FFEA00",  # Chrome Yellow
    "#CC7722",  # Yellow Ochre
    "#E34234",  # Vermilion
    "#E30022",  # Cadmium Red
    "#D68A59",  # Raw Sienna
    "#8A360F",  # Burnt Sienna
    "#826644",  # Raw Umber
]
colors3 = [
    "#FFD700", # Glowing Yellow, reminiscent of "Starry Night" and "Sunflowers"
    "#6495ED", # Cornflower Blue, reflective of the night sky in "Starry Night"
    "#FF4500", # OrangeRed, for the vibrant tones in "Cafe Terrace at Night"
    "#006400", # DarkGreen, echoing the cypress trees and fields
    "#8B4513", # SaddleBrown, for the earth and branches in many landscapes
    "#DAA520", # Goldenrod, similar to the tones in "Wheatfield with Cypresses"
    "#00008B", # DarkBlue, for the deep night skies
    "#008000", # Green, capturing the vitality of nature scenes
    "#BDB76B", # DarkKhaki, for the muted tones in "The Potato Eaters"
    "#800080", # Purple, reflecting the subtle touches in flowers and clothing
    "#FF6347", # Tomato, for bright accents in paintings like "Red Vineyards at Arles"
    "#4682B4", # SteelBlue, for the serene moments in "Fishing Boats on the Beach"
    "#FA8072", # Salmon, for the soft glow of sunset and sunrise scenes
    "#9ACD32", # YellowGreen, for the lively foliage in many landscapes
    "#40E0D0", # Turquoise, reminiscent of the dynamic strokes in "Irises"
    "#BA55D3", # MediumOrchid, for the playful color in "Almond Blossoms"
    "#7FFF00", # Chartreuse, vibrant and lively for the touches of light
    "#ADD8E6", # LightBlue, for the peaceful skies and distant horizons
]

color4 = ['#3182bd','#6baed6','#9ecae1','#e6550d','#fd8d3c','#fdd0a2','#31a354','#a1d99b', '#c7e9c0', '#756bb1', '#bcbddc', '#dadaeb', '#969696', '#bdbdbd','#d9d9d9']
custom_cmap = LinearSegmentedColormap.from_list("custom_theme", colors2[1:])
```

```{code-cell} ipython3
samp = pd.read_csv('data/AGNsample_26Feb24.csv')

df_lc = pd.read_parquet('data/df_lc_022624.parquet')
objids = df_lc.index.get_level_values('objectid')[:].unique()
redshifts = samp['redshift']#[objids]
df_lc
```

```{code-cell} ipython3
from plot_functions import create_figure
grouped = list(df_lc.groupby('objectid'))

for ind in range(1981,1982):
    objectid, singleobj_df = grouped[ind]
    print(samp.iloc[objectid])
    _ = create_figure(df_lc = df_lc, 
                       index = ind,  
                       save_output = False,  # should the resulting plots be saved?
                      )
```

```{code-cell} ipython3
x_ztf = np.linspace(0, 1850, 175)  # For ZTF
kernel = RationalQuadratic(length_scale=1, alpha=0.1)
colors = ['#3182bd','#6baed6','#9ecae1','#e6550d','#fd8d3c','#fdd0a2','#31a354','#a1d99b', '#c7e9c0', '#756bb1', '#bcbddc', '#dadaeb', '#969696', '#bdbdbd','#d9d9d9']


for keepindex, obj in tqdm(enumerate([objectid])):    
    singleobj = df_lc.loc[obj, :, :, :]  # Extract data for the single object
    label = singleobj.index.unique('label')  # Get the label of the object
    bands = singleobj.loc[label[0], :, :].index.get_level_values('band')[:].unique()  # Extract bands

plt.figure(figsize=(6, 6))  # Set up plot if within numplots limit
plt.subplot(2,1,1)
added_labels = {}

for l, band in enumerate(['zg','zr']):
    band_lc = singleobj.loc[label[0], band, :]  # Extract light curve data for the band
    band_lc_clean = band_lc[band_lc.index.get_level_values('time') < 65000]
    x, y, dy = np.array(band_lc_clean.index.get_level_values('time') - band_lc_clean.index.get_level_values('time')[0]), np.array(band_lc_clean.flux), np.array(band_lc_clean.err)
    # Sort data based on time
    x2, y2, dy2 = x[np.argsort(x)], y[np.argsort(x)], dy[np.argsort(x)]

    # Check if there are enough points for interpolation
    if (len(x2) > 10) and not np.isnan(y2).any():
        n = np.sum(x2 == 0)
        for b in range(1, n):
            x2[::b + 1] = x2[::b + 1] + 1 * 0.001

        # Interpolate the data
        f = interpolate.interp1d(x2, y2, kind='previous', fill_value="extrapolate")
        df = interpolate.interp1d(x2, dy2, kind='previous', fill_value="extrapolate")
        l = 'nearest interpolation'
        if l not in added_labels:
            gline, = plt.plot(x_ztf, f(x_ztf), '--', label= l,color = '#92a8d1',alpha=1)
            added_labels[l] = True            
        else:
            gline, = plt.plot(x_ztf, f(x_ztf), '--',color = "#3F51B5",alpha=1)

        gcolor=gline.get_color()
        l = 'observed data'
        if l not in added_labels:
            plt.errorbar(x2, y2, dy2, capsize=1.0, marker='.',label=l, linestyle='',alpha=0.5,color='#92a8d1')
            added_labels[l] = True
        else:
            plt.errorbar(x2, y2, dy2, capsize=1.0, marker='.', linestyle='',alpha=0.5,color="#3F51B5")

        X = x2.reshape(-1, 1)        
        x_ztf = np.linspace(0,1850,175).reshape(-1, 1) # X array for interpolation
        gp = GaussianProcessRegressor(kernel=kernel, alpha=dy2**2)
        gp.fit(X, y2)
        y_pred,sigma = gp.predict(x_ztf, return_std=True)
        l = 'Gaussian Process Reg.'
        if l not in added_labels:
            gpline, = plt.plot(x_ztf,y_pred,'-',label=l,color = '#92a8d1')
            added_labels[l] = True
        else:
            gpline, = plt.plot(x_ztf,y_pred,'-',color = "#3F51B5")  
        
        gcolor= gpline.get_color()
        plt.fill_between(x_ztf.flatten(), y_pred - 1.96 * sigma,y_pred + 1.96 * sigma, alpha=0.2, color=gcolor)

plt.grid()
plt.text(20,0.1,'ZTF g,r',size=15)

#plt.xlabel(r'$\rm time(day)$',size=15)
#plt.ylim([0,0.1])
plt.ylabel(r'$\rm Flux(mJy)$',size=15)
plt.legend(loc=1)
plt.subplot(2,1,2)
x_ztf = np.linspace(0, 4000, 175)  # For ZTF


added_labels = {}

for l, band in enumerate(['W1','W2']):
    band_lc = singleobj.loc[label[0], band, :]  # Extract light curve data for the band
    band_lc_clean = band_lc[band_lc.index.get_level_values('time') < 65000]
    x, y, dy = np.array(band_lc_clean.index.get_level_values('time') - band_lc_clean.index.get_level_values('time')[0]), np.array(band_lc_clean.flux), np.array(band_lc_clean.err)
    # Sort data based on time
    x2, y2, dy2 = x[np.argsort(x)], y[np.argsort(x)], dy[np.argsort(x)]

    # Check if there are enough points for interpolation
    if (len(x2) > 10) and not np.isnan(y2).any():
        n = np.sum(x2 == 0)
        for b in range(1, n):
            x2[::b + 1] = x2[::b + 1] + 1 * 0.001

        # Interpolate the data
        f = interpolate.interp1d(x2, y2, kind='previous', fill_value="extrapolate")
        df = interpolate.interp1d(x2, dy2, kind='previous', fill_value="extrapolate")
        l = 'nearest interpolation'
        if l not in added_labels:
            gline, = plt.plot(x_ztf, f(x_ztf), '--', label= l,color = '#30D5C8',alpha=1)
            added_labels[l] = True            
        else:
            gline, = plt.plot(x_ztf, f(x_ztf), '--',color = '#7fcdbb',alpha=1)

        gcolor=gline.get_color()
        l = 'observed data'
        if l not in added_labels:
            plt.errorbar(x2, y2, dy2, capsize=1.0, marker='.',label=l, linestyle='',alpha=0.8,color=colors[0])
            added_labels[l] = True
        else:
            plt.errorbar(x2, y2, dy2, capsize=1.0, marker='.', linestyle='',alpha=0.8,color=colors[0])

        X = x2.reshape(-1, 1)        
        x_ztf = np.linspace(0,4000,175).reshape(-1, 1) # X array for interpolation
        gp = GaussianProcessRegressor(kernel=kernel, alpha=dy2**2)
        gp.fit(X, y2)
        y_pred,sigma = gp.predict(x_ztf, return_std=True)
        l = 'Gaussian Process Reg.'
        if l not in added_labels:
            gpline, = plt.plot(x_ztf,y_pred,'-',label=l,color = '#30D5C8')
            added_labels[l] = True
        else:
            gpline, = plt.plot(x_ztf,y_pred,'-',color = '#7fcdbb')

        gcolor= gpline.get_color()
        plt.fill_between(x_ztf.flatten(), y_pred - 1.96 * sigma,y_pred + 1.96 * sigma, alpha=0.2, color=gcolor)
plt.text(20,0.16,'WISE W1,W2',size=15)
plt.grid()
#plt.xlim([-10,1880])
plt.xlabel(r'$\rm time(day)$',size=15)
plt.ylabel(r'$\rm Flux(mJy)$',size=15)

plt.legend(loc=1)
plt.tight_layout()
#plt.savefig('output/unify_lc1994.png')
```

```{code-cell} ipython3
bands_inlc = ['W1']
numobjs = len(df_lc.index.get_level_values('objectid')[:].unique())
sample_objids = df_lc.index.get_level_values('objectid').unique()[:numobjs]
#df_lc_small = df_lc.loc[sample_objids]
objects,dobjects,flabels,zlist,keeps = unify_lc_gp_parallel(df_lc_small,redshifts,bands_inlc=bands_inlc,xres=160)
objects2,dobjects2,flabels2,zlist2,keeps2 = unify_lc(df_lc_small,redshifts,bands_inlc=bands_inlc,xres=160)

fvar, maxarray, meanarray = stat_bands(objects,dobjects,bands_inlc,sigmacl=5)
dat_notnormal = combine_bands(objects,bands_inlc) 
datm = normalize_clipmax_objects(dat_notnormal,meanarray,band = -1)

# shuffle data incase the ML routines are sensitive to order
data,fzr,p = shuffle_datalabel(datm,flabels)
fvar_arr,maximum_arr,average_arr = fvar[:,p],maxarray[:,p],meanarray[:,p]
#redshift_shuffled = zlist[p]

labc = {}  # Initialize labc to hold indices of each unique label
for index, f in enumerate(fzr):
    lab = translate_bitwise_sum_to_labels(int(f))
    for label in lab:
        if label not in labc:
            labc[label] = []  # Initialize the list for this label if it's not already in labc
        labc[label].append(index)  # Append the current index to the list of indices for this label

fvar2, maxarray2, meanarray2 = stat_bands(objects2,dobjects2,bands_inlc,sigmacl=5)
dat_notnormal2 = combine_bands(objects2,bands_inlc) 
datm2 = normalize_clipmax_objects(dat_notnormal2,meanarray2,band = -1)
data2,fzr2,p2 = shuffle_datalabel(datm2,flabels2)
fvar_arr2,maximum_arr2,average_arr2 = fvar2[:,p2],maxarray2[:,p2],meanarray2[:,p2]
#redshift_shuffled = zlist[p]

labc2 = {}  # Initialize labc to hold indices of each unique label
for index, f in enumerate(fzr2):
    lab = translate_bitwise_sum_to_labels(int(f))
    for label in lab:
        if label not in labc2:
            labc2[label] = []  # Initialize the list for this label if it's not already in labc
        labc2[label].append(index)  # Append the current index to the list of indices for this label

  

#Assuming `data` is your numpy array
#nan_rows = np.any(np.isnan(data), axis=1)
#clean_data = data[~nan_rows]  # Rows without NaNs
#clean_fzr = fzr[~nan_rows]
#labc3 = {}  # Initialize labc to hold indices of each unique label
#for index, f in enumerate(clean_fzr):
#    lab = translate_bitwise_sum_to_labels(int(f))
#    for label in lab:
#        if label not in labc3:
#            labc3[label] = []  # Initialize the list for this label if it's not already in labc
#        labc3[label].append(index)  # Append the current index to the list of indices for this label

#fvar_arr1,average_arr1,redshift1 = fvar_arr[:,~nan_rows],average_arr[:,~nan_rows],redshift_shuffled[~nan_rows]
#np.savez('data/sampleA_w1',data=clean_data,fzr = clean_fzr,fvar_arr1 = fvar_arr1, average_arr1 = average_arr1,redshift_shuffled=redshift1,labc = labc3)

#d = np.load('data/sampleA_w1.npz',allow_pickle=True)
#data = d['data']
#fvar_arr,average_arr = d['fvar_arr1'],d['average_arr1']
#redshifts,fzr,labc = d['redshift_shuffled'], d['fzr'],d['labc']                                                                             
#print(np.min(redshifts),np.mean(redshifts),len(redshifts))
```

```{code-cell} ipython3
#mapper_e = umap.UMAP(n_neighbors=100,min_dist=0.99,metric='euclidean',random_state=20).fit(data)
#mapper_m = umap.UMAP(n_neighbors=100,min_dist=0.9,metric='manhattan',random_state=20).fit(data)
#mapper_d = umap.UMAP(n_neighbors=100,min_dist=0.99,metric=dtw_distance,random_state=1).fit(data)

#mapper_e2 = umap.UMAP(n_neighbors=100,min_dist=0.99,metric='euclidean',random_state=20).fit(data2)
#mapper_m2 = umap.UMAP(n_neighbors=100,min_dist=0.9,metric='manhattan',random_state=20).fit(data2)
#mapper_d2 = umap.UMAP(n_neighbors=100,min_dist=0.99,metric=dtw_distance,random_state=1).fit(data2)

qq='k'
plt.figure(figsize=(12,8))

markersize=100
cmap1 = 'viridis'

ax1 = plt.subplot(2,3,4)
thiscolor=stretch_small_values_arctan(np.nansum(fvar_arr,axis=0),factor=15)
u = (thiscolor<2.) & (thiscolor>=0)
cf = ax1.scatter(mapper_e.embedding_[u,0],mapper_e.embedding_[u,1],c = thiscolor[u],s=markersize,edgecolor=qq,cmap=cmap1)
#plt.axis('off')
ax1.tick_params(axis='both',          # Changes apply to both x and y-axis
               which='both',         # Both major and minor ticks are affected
               bottom=False,         # Ticks along the bottom edge are off
               top=False,            # Ticks along the top edge are off
               left=False,           # Ticks along the left edge are off
               right=False)          # Ticks along the right edge are off

# Optional: Hide the tick labels (if you want to keep the ticks but remove labels, comment these out)
ax1.set_xticklabels([])
ax1.set_yticklabels([])
ax1.set_ylabel(r'$\rm GP\ regression$',size=20)

ax1 = plt.subplot(2,3,5)
thiscolor=stretch_small_values_arctan(np.nansum(fvar_arr,axis=0),factor=15)
u = (thiscolor<2.) & (thiscolor>=0)
cf = ax1.scatter(mapper_m.embedding_[u,0],mapper_m.embedding_[u,1],c = thiscolor[u],s=markersize,edgecolor=qq,cmap=cmap1)
ax1.tick_params(axis='both',          # Changes apply to both x and y-axis
               which='both',         # Both major and minor ticks are affected
               bottom=False,         # Ticks along the bottom edge are off
               top=False,            # Ticks along the top edge are off
               left=False,           # Ticks along the left edge are off
               right=False)          # Ticks along the right edge are off

ax1.set_xticklabels([])
ax1.set_yticklabels([])

ax1 = plt.subplot(2,3,6)
thiscolor=stretch_small_values_arctan(np.nansum(fvar_arr,axis=0),factor=15)
u = (thiscolor<2.) & (thiscolor>=0)
cf = ax1.scatter(mapper_d.embedding_[u,0],mapper_d.embedding_[u,1],c = thiscolor[u],s=markersize,edgecolor=qq,cmap=cmap1)
ax1.tick_params(axis='both',          # Changes apply to both x and y-axis
               which='both',         # Both major and minor ticks are affected
               bottom=False,         # Ticks along the bottom edge are off
               top=False,            # Ticks along the top edge are off
               left=False,           # Ticks along the left edge are off
               right=False)          # Ticks along the right edge are off

ax1.set_xticklabels([])
ax1.set_yticklabels([])


ax1 = plt.subplot(2,3,1)
ax1.set_title(r'$\rm Euclidean$',size=20)
thiscolor=stretch_small_values_arctan(np.nansum(fvar_arr2,axis=0),factor=15)
u = (thiscolor<2.) & (thiscolor>=0)
cf = ax1.scatter(mapper_e2.embedding_[u,0],mapper_e2.embedding_[u,1],c = thiscolor[u],s=markersize,edgecolor=qq,cmap=cmap1)
ax1.tick_params(axis='both',          # Changes apply to both x and y-axis
               which='both',         # Both major and minor ticks are affected
               bottom=False,         # Ticks along the bottom edge are off
               top=False,            # Ticks along the top edge are off
               left=False,           # Ticks along the left edge are off
               right=False)          # Ticks along the right edge are off

ax1.set_xticklabels([])
ax1.set_yticklabels([])
ax1.set_ylabel(r'$\rm NN\ linear\ interpolation$',size=20)

ax1 = plt.subplot(2,3,2)
ax1.set_title(r'$\rm Manhattan$',size=20)
thiscolor=stretch_small_values_arctan(np.nansum(fvar_arr2,axis=0),factor=15)
u = (thiscolor<2.) & (thiscolor>=0)
cf = ax1.scatter(mapper_m2.embedding_[u,0],mapper_m2.embedding_[u,1],c = thiscolor[u],s=markersize,edgecolor=qq,cmap=cmap1)
ax1.tick_params(axis='both',          # Changes apply to both x and y-axis
               which='both',         # Both major and minor ticks are affected
               bottom=False,         # Ticks along the bottom edge are off
               top=False,            # Ticks along the top edge are off
               left=False,           # Ticks along the left edge are off
               right=False)          # Ticks along the right edge are off

ax1.set_xticklabels([])
ax1.set_yticklabels([])

ax1 = plt.subplot(2,3,3)
ax1.set_title(r'$\rm DTW$',size=20)
thiscolor=stretch_small_values_arctan(np.nansum(fvar_arr2,axis=0),factor=15)
u = (thiscolor<2.) & (thiscolor>=0)
cf = ax1.scatter(mapper_d2.embedding_[u,0],mapper_d2.embedding_[u,1],c = thiscolor[u],s=markersize,edgecolor=qq,cmap=cmap1)
ax1.tick_params(axis='both',          # Changes apply to both x and y-axis
               which='both',         # Both major and minor ticks are affected
               bottom=False,         # Ticks along the bottom edge are off
               top=False,            # Ticks along the top edge are off
               left=False,           # Ticks along the left edge are off
               right=False)          # Ticks along the right edge are off

ax1.set_xticklabels([])
ax1.set_yticklabels([])

plt.tight_layout()
plt.subplots_adjust(hspace=0, wspace=0)

plt.savefig('output/umap_params.png')
```

```{code-cell} ipython3
mapper = umap.UMAP(n_neighbors=100,min_dist=0.99,metric=dtw_distance,random_state=1).fit(data)

plt.figure(figsize=(12,4))
markersize=100
cmap1 = 'viridis'

ax1 = plt.subplot(1,3,1)
ax1.set_title(r'$\rm Mean\ brightness$')
thiscolor=np.log10(np.nansum(average_arr,axis=0))
u = (thiscolor<2) & (thiscolor>=-2)
cf = ax1.scatter(mapper.embedding_[u,0],mapper.embedding_[u,1],c = thiscolor[u],s=markersize,edgecolor='k',cmap=cmap1)
plt.axis('off')
divider = make_axes_locatable(ax1)
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(cf,cax=cax)


ax1 = plt.subplot(1,3,3)
ax1.set_title(r'$\rm Mean\ Fractional\ Variation$')
thiscolor=stretch_small_values_arctan(np.nansum(fvar_arr,axis=0),factor=15)
u = (thiscolor<2.) & (thiscolor>=0)
cf = ax1.scatter(mapper.embedding_[u,0],mapper.embedding_[u,1],c = thiscolor[u],s=markersize,edgecolor='k',cmap=cmap1)
plt.axis('off')
divider = make_axes_locatable(ax1)
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(cf,cax=cax)


ax1 = plt.subplot(1,3,2)
ax1.set_title(r'$\rm Redshift$')
thiscolor=redshift_shuffled
u = (thiscolor<2) & (thiscolor>=0)
cf = ax1.scatter(mapper.embedding_[u,0],mapper.embedding_[u,1],c = thiscolor[u],s=markersize,edgecolor='k',cmap=cmap1)
plt.axis('off')
divider = make_axes_locatable(ax1)
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(cf,cax=cax)

plt.tight_layout()
plt.savefig('output/umap-w1-sampleA-1.png')
```

```{code-cell} ipython3
# Calculate 2D histogram
hist, x_edges, y_edges = np.histogram2d(mapper.embedding_[:, 0], mapper.embedding_[:, 1], bins=10)
plt.figure(figsize=(15,10))
i=1
laborder = ['SDSS_QSO','WISE_Variable','Optical_Variable','Galex_Variable','SPIDER_AGN','SPIDER_AGNBL','SPIDER_QSOBL','SPIDER_BL','Turn-on','Turn-off','TDE','Fermi_Blazars']
for label in laborder:
    if label in labc:
        indices = labc[label]
        hist_per_cluster, _, _ = np.histogram2d(mapper.embedding_[indices,0], mapper.embedding_[indices,1], bins=(x_edges, y_edges))
        prob = hist_per_cluster / hist
        plt.subplot(3,4,i)
        plt.title(label)
        plt.contourf(x_edges[:-1], y_edges[:-1], prob.T, levels=20, alpha=0.8,cmap=custom_cmap)
        plt.colorbar()
        plt.axis('off')
        #cf = ax0.scatter(mapper.embedding_[indices,0],mapper.embedding_[indices,1],s=80,alpha=0.5,edgecolor='gray',label=label,c=colors[i-1])
        i+=1
        
ax2 = plt.subplot(3,4,12)
ax2.set_title('sample origin',size=20)
counts = 2
for label, indices in labc.items():
    cf = ax2.scatter(mapper.embedding_[indices,0],mapper.embedding_[indices,1],s=markersize,c = color4[counts],alpha=0.8,edgecolor='k',label=label)
    counts+=1
plt.legend(loc=4,fontsize=8)
plt.axis('off')

plt.tight_layout()
plt.savefig('output/umap-w1-sampleA-2.png')
```

```{code-cell} ipython3
bands_inlc = ['zg','zr','zi','W1','W2']
numobjs = len(df_lc.index.get_level_values('objectid')[:].unique())
sample_objids = df_lc.index.get_level_values('objectid').unique()[:numobjs]
df_lc_small = df_lc.loc[sample_objids]
objects,dobjects,flabels,zlist,keeps = unify_lc_gp_parallel(df_lc_small,redshifts,bands_inlc=bands_inlc,xres=160)

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
```

```{code-cell} ipython3
#mapper = umap.UMAP(n_neighbors=100,min_dist=0.99,metric=dtw_distance,random_state=3).fit(data)
#mapper = umap.UMAP(n_neighbors=100,min_dist=0.99,metric='manhattan',random_state=20).fit(data)

plt.figure(figsize=(12,4))
markersize=100
cmap1 = 'viridis'

ax1 = plt.subplot(1,3,1)
ax1.set_title(r'$\rm Mean\ brightness$')
thiscolor=np.log10(np.nansum(average_arr,axis=0))
u = (thiscolor<2) & (thiscolor>=-2)
cf = ax1.scatter(mapper.embedding_[u,0],mapper.embedding_[u,1],c = thiscolor[u],s=markersize,edgecolor='k',cmap=cmap1)
plt.axis('off')
divider = make_axes_locatable(ax1)
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(cf,cax=cax)


ax1 = plt.subplot(1,3,3)
ax1.set_title(r'$\rm Mean\ Fractional\ Variation$')
thiscolor=stretch_small_values_arctan(np.nansum(fvar_arr,axis=0),factor=3)
u = (thiscolor<1.5) & (thiscolor>=0)
cf = ax1.scatter(mapper.embedding_[u,0],mapper.embedding_[u,1],c = thiscolor[u],s=markersize,edgecolor='k',cmap=cmap1)
plt.axis('off')
divider = make_axes_locatable(ax1)
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(cf,cax=cax)


ax1 = plt.subplot(1,3,2)
ax1.set_title(r'$\rm Redshift$')
thiscolor=redshift_shuffled
u = (thiscolor<0.8) & (thiscolor>=0)
cf = ax1.scatter(mapper.embedding_[u,0],mapper.embedding_[u,1],c = thiscolor[u],s=markersize,edgecolor='k',cmap=cmap1)
plt.axis('off')
divider = make_axes_locatable(ax1)
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(cf,cax=cax)

plt.tight_layout()
plt.savefig('output/umap-ztfw-sampleA-1.png')
```

```{code-cell} ipython3
plt.figure(figsize=(12,8))
markersize=100

hist, x_edges, y_edges = np.histogram2d(mapper.embedding_[:, 0], mapper.embedding_[:, 1], bins=10)
plt.figure(figsize=(15,10))
i=1
laborder = ['SDSS_QSO','WISE_Variable','Optical_Variable','Galex_Variable','SPIDER_AGN','SPIDER_AGNBL','SPIDER_QSOBL','SPIDER_BL','Turn-on','Turn-off','TDE','Fermi_Blazars']
for label in laborder:
    if label in labc:
        indices = labc[label]
        hist_per_cluster, _, _ = np.histogram2d(mapper.embedding_[indices,0], mapper.embedding_[indices,1], bins=(x_edges, y_edges))
        prob = hist_per_cluster / hist
        plt.subplot(3,4,i)
        plt.title(label)
        plt.contourf(x_edges[:-1], y_edges[:-1], prob.T, levels=20, alpha=0.8,cmap=custom_cmap)
        plt.colorbar()
        plt.axis('off')
        #cf = ax0.scatter(mapper.embedding_[indices,0],mapper.embedding_[indices,1],s=80,alpha=0.5,edgecolor='gray',label=label,c=colors[i-1])
        i+=1
ax2 = plt.subplot(3,4,12)
ax2.set_title('sample origin',size=20)
counts = 1
for label, indices in labc.items():
    cf = ax2.scatter(mapper.embedding_[indices,0],mapper.embedding_[indices,1],s=markersize,c = color4[counts],alpha=0.8,edgecolor='k',label=label)
    counts+=1
plt.legend(loc=3,fontsize=8)
plt.axis('off')

plt.tight_layout()
plt.savefig('output/umap-ztfw-sampleA-2.png')
```

```{code-cell} ipython3
df_lc.index.get_level_values('band')[:].unique()
```

```{code-cell} ipython3
bands_inlc = ['G', 'BP', 'RP', 'panstarrs y', 'panstarrs i', 'panstarrs z','panstarrs r', 'panstarrs g', 'W1', 'W2', 'zg', 'zi', 'zr']
numobjs = len(df_lc.index.get_level_values('objectid')[:].unique())
sample_objids = df_lc.index.get_level_values('objectid').unique()[:numobjs]
df_lc_small = df_lc.loc[sample_objids]
objects,dobjects,flabels,zlist,keeps = unify_lc_gp_parallel(df_lc_small,redshifts,bands_inlc=bands_inlc,xres=60)

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

nan_rows = np.any(np.isnan(data), axis=1)
clean_data = data[~nan_rows]  # Rows without NaNs
fvar_arr3,average_arr3 = fvar_arr[:,~nan_rows],average_arr[:,~nan_rows]
redshifts3= redshift_shuffled[~nan_rows]

clean_fzr = fzr[~nan_rows]
labc3 = {}  # Initialize labc to hold indices of each unique label
for index, f in enumerate(clean_fzr):
    lab = translate_bitwise_sum_to_labels(int(f))
    for label in lab:
        if label not in labc3:
            labc3[label] = []  # Initialize the list for this label if it's not already in labc
        labc3[label].append(index) 
```

```{code-cell} ipython3
#mapper = umap.UMAP(n_neighbors=100,min_dist=0.99,metric=dtw_distance,random_state=3).fit(data)
mapper = umap.UMAP(n_neighbors=100,min_dist=0.99,metric='manhattan',random_state=20).fit(clean_data)

plt.figure(figsize=(12,4))
markersize=100
cmap1 = 'viridis'

ax1 = plt.subplot(1,3,1)
ax1.set_title(r'$\rm Mean\ brightness$')
thiscolor=np.log10(np.nansum(average_arr3,axis=0))
u = (thiscolor<2) & (thiscolor>=-2)
cf = ax1.scatter(mapper.embedding_[u,0],mapper.embedding_[u,1],c = thiscolor[u],s=markersize,edgecolor='k',cmap=cmap1)
plt.axis('off')
divider = make_axes_locatable(ax1)
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(cf,cax=cax)


ax1 = plt.subplot(1,3,3)
ax1.set_title(r'$\rm Mean\ Fractional\ Variation$')
thiscolor=stretch_small_values_arctan(np.nansum(fvar_arr3,axis=0),factor=3)
u = (thiscolor<1.5) & (thiscolor>=0)
cf = ax1.scatter(mapper.embedding_[u,0],mapper.embedding_[u,1],c = thiscolor[u],s=markersize,edgecolor='k',cmap=cmap1)
plt.axis('off')
divider = make_axes_locatable(ax1)
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(cf,cax=cax)


ax1 = plt.subplot(1,3,2)
ax1.set_title(r'$\rm Redshift$')
thiscolor=redshifts3
u = (thiscolor<0.8) & (thiscolor>=0)
cf = ax1.scatter(mapper.embedding_[u,0],mapper.embedding_[u,1],c = thiscolor[u],s=markersize,edgecolor='k',cmap=cmap1)
plt.axis('off')
divider = make_axes_locatable(ax1)
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(cf,cax=cax)

plt.tight_layout()
plt.savefig('output/umap-all-sampleA-1.png')
```

```{code-cell} ipython3
plt.figure(figsize=(12,8))
markersize=100

hist, x_edges, y_edges = np.histogram2d(mapper.embedding_[:, 0], mapper.embedding_[:, 1], bins=10)
plt.figure(figsize=(15,10))
i=1
laborder = ['SDSS_QSO','WISE_Variable','Optical_Variable','Galex_Variable','SPIDER_AGN','SPIDER_AGNBL','SPIDER_QSOBL','SPIDER_BL','Turn-on','Turn-off','TDE','Fermi_Blazars']
for label in laborder:
    if label in labc3:
        indices = labc3[label]
        hist_per_cluster, _, _ = np.histogram2d(mapper.embedding_[indices,0], mapper.embedding_[indices,1], bins=(x_edges, y_edges))
        prob = hist_per_cluster / hist
        plt.subplot(3,4,i)
        plt.title(label)
        plt.contourf(x_edges[:-1], y_edges[:-1], prob.T, levels=20, alpha=0.8,cmap=custom_cmap)
        plt.colorbar()
        plt.axis('off')
        #cf = ax0.scatter(mapper.embedding_[indices,0],mapper.embedding_[indices,1],s=80,alpha=0.5,edgecolor='gray',label=label,c=colors[i-1])
        i+=1
ax2 = plt.subplot(3,4,12)
ax2.set_title('sample origin',size=20)
counts = 1
for label, indices in labc3.items():
    cf = ax2.scatter(mapper.embedding_[indices,0],mapper.embedding_[indices,1],s=markersize,c = color4[counts],alpha=0.8,edgecolor='k',label=label)
    counts+=1
plt.legend(loc=3,fontsize=8)
plt.axis('off')

plt.tight_layout()
plt.savefig('output/umap-all-sampleA-2.png')
```

```{code-cell} ipython3
bands_inlc = ['zg']
numobjs = len(df_lc.index.get_level_values('objectid')[:].unique())
sample_objids = df_lc.index.get_level_values('objectid').unique()[:numobjs]
df_lc_small = df_lc.loc[sample_objids]
objects,dobjects,flabels,zlist,keeps = unify_lc_gp_parallel(df_lc_small,redshifts,bands_inlc=bands_inlc,xres=160)

# calculate some basic statistics with a sigmaclipping with width 5sigma
fvar, maxarray, meanarray = stat_bands(objects,dobjects,bands_inlc,sigmacl=5)

# combine different waveband into one array
dat_notnormal = combine_bands(objects,bands_inlc) 

# Normalize the combinde array by mean brightness in a waveband after clipping outliers:
datm = normalize_clipmax_objects(dat_notnormal,meanarray,band = -1)

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
```

```{code-cell} ipython3
nan_rows = np.any(np.isnan(data), axis=1)
clean_data = data[~nan_rows]  # Rows without NaNs
fvar_arr3,average_arr3 = fvar_arr[:,~nan_rows],average_arr[:,~nan_rows]
redshifts3= redshift_shuffled[~nan_rows]

clean_fzr = fzr[~nan_rows]
labc3 = {}  # Initialize labc to hold indices of each unique label
for index, f in enumerate(clean_fzr):
    lab = translate_bitwise_sum_to_labels(int(f))
    for label in lab:
        if label not in labc3:
            labc3[label] = []  # Initialize the list for this label if it's not already in labc
        labc3[label].append(index) 
np.savez('data/sampleA_zg',data=clean_data,fzr = clean_fzr,fvar_arr = fvar_arr3, average_arr = average_arr3,redshift_shuffled=redshifts3,labc = labc3)
```

```{code-cell} ipython3
#mapper = umap.UMAP(n_neighbors=100,min_dist=0.99,metric=dtw_distance,random_state=3).fit(data)
mapper = umap.UMAP(n_neighbors=100,min_dist=0.99,metric='manhattan',random_state=20).fit(clean_data)

plt.figure(figsize=(12,4))
markersize=100
cmap1 = 'viridis'

ax1 = plt.subplot(1,3,1)
ax1.set_title(r'$\rm Mean\ brightness$')
thiscolor=np.log10(np.nansum(average_arr3,axis=0))
u = (thiscolor<2) & (thiscolor>=-2)
cf = ax1.scatter(mapper.embedding_[u,0],mapper.embedding_[u,1],c = thiscolor[u],s=markersize,edgecolor='k',cmap=cmap1)
plt.axis('off')
divider = make_axes_locatable(ax1)
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(cf,cax=cax)


ax1 = plt.subplot(1,3,3)
ax1.set_title(r'$\rm Mean\ Fractional\ Variation$')
thiscolor=stretch_small_values_arctan(np.nansum(fvar_arr3,axis=0),factor=3)
u = (thiscolor<1.5) & (thiscolor>=0)
cf = ax1.scatter(mapper.embedding_[u,0],mapper.embedding_[u,1],c = thiscolor[u],s=markersize,edgecolor='k',cmap=cmap1)
plt.axis('off')
divider = make_axes_locatable(ax1)
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(cf,cax=cax)


ax1 = plt.subplot(1,3,2)
ax1.set_title(r'$\rm Redshift$')
thiscolor=redshifts3
u = (thiscolor<0.8) & (thiscolor>=0)
cf = ax1.scatter(mapper.embedding_[u,0],mapper.embedding_[u,1],c = thiscolor[u],s=markersize,edgecolor='k',cmap=cmap1)
plt.axis('off')
divider = make_axes_locatable(ax1)
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(cf,cax=cax)

plt.tight_layout()
plt.savefig('output/umap-ztfg-sampleA-1.png')
```

```{code-cell} ipython3
plt.figure(figsize=(12,8))
markersize=100

hist, x_edges, y_edges = np.histogram2d(mapper.embedding_[:, 0], mapper.embedding_[:, 1], bins=10)
plt.figure(figsize=(15,10))
i=1
laborder = ['SDSS_QSO','WISE_Variable','Optical_Variable','Galex_Variable','SPIDER_AGN','SPIDER_AGNBL','SPIDER_QSOBL','SPIDER_BL','Turn-on','Turn-off','TDE','Fermi_Blazars']
for label in laborder:
    if label in labc3:
        indices = labc3[label]
        hist_per_cluster, _, _ = np.histogram2d(mapper.embedding_[indices,0], mapper.embedding_[indices,1], bins=(x_edges, y_edges))
        prob = hist_per_cluster / hist
        plt.subplot(3,4,i)
        plt.title(label)
        plt.contourf(x_edges[:-1], y_edges[:-1], prob.T, levels=20, alpha=0.8,cmap=custom_cmap)
        plt.colorbar()
        plt.axis('off')
        #cf = ax0.scatter(mapper.embedding_[indices,0],mapper.embedding_[indices,1],s=80,alpha=0.5,edgecolor='gray',label=label,c=colors[i-1])
        i+=1
ax2 = plt.subplot(3,4,12)
ax2.set_title('sample origin',size=20)
counts = 1
for label, indices in labc3.items():
    cf = ax2.scatter(mapper.embedding_[indices,0],mapper.embedding_[indices,1],s=markersize,c = color4[counts],alpha=0.8,edgecolor='k',label=label)
    counts+=1
plt.legend(loc=3,fontsize=8)
plt.axis('off')

plt.tight_layout()
plt.savefig('output/umap-ztfg-sampleA-2.png')
```

```{code-cell} ipython3
msz0,msz1 = 15,15
sm = sompy.SOMFactory.build(data, mapsize=[msz0,msz1], mapshape='planar', lattice='rect', initialization='pca')
sm.train(n_job=4, shared_memory = 'no')
```

```{code-cell} ipython3
a=sm.bmu_ind_to_xy(sm.project_data(data))
x,y=np.zeros(len(a)),np.zeros(len(a))
k=0
for i in a:
    x[k]=i[0]
    y[k]=i[1]
    k+=1
med_r=np.zeros([msz0,msz1])
fvar_new = stretch_small_values_arctan(np.nansum(redshift_shuffled,axis=0),factor=1)
for i in range(msz0):
    for j in range(msz1):
        unja=(x==i)&(y==j)
        med_r[i,j]=(np.nanmedian(redshift_shuffled[unja]))


plt.figure(figsize=(18,12))
plt.subplot(3,4,1)
plt.title('SDSS',fontsize=15)
cf=plt.imshow(med_r,origin='lower',cmap='viridis')
plt.axis('off')

def get_indices_by_label(labc, label):
    return labc.get(label, [])

# Example usage
u = labc.get('SDSS_QSO',[])
dsdss = data[u,:]
asdss=sm.bmu_ind_to_xy(sm.project_data(dsdss))
x,y=np.zeros(len(asdss)),np.zeros(len(asdss))
k=0
for i in asdss:
    x[k]=i[0]
    y[k]=i[1]
    k+=1
plt.plot(y,x,'rx',alpha=0.8)

plt.subplot(3,4,2)
plt.title('WISE_Variable',fontsize=15)
cf=plt.imshow(med_r,origin='lower',cmap='viridis')
plt.axis('off')
u = labc.get('WISE_Variable',[])
dwise = data[u,:]
awise=sm.bmu_ind_to_xy(sm.project_data(dwise))
x,y=np.zeros(len(awise)),np.zeros(len(awise))
k=0
for i in awise:
    x[k]=i[0]
    y[k]=i[1]
    k+=1
plt.plot(y,x,'rx',alpha=0.8)

plt.subplot(3,4,3)
plt.title('Optical_Variable',fontsize=15)
cf=plt.imshow(med_r,origin='lower',cmap='viridis')
plt.axis('off')
u = labc.get('Optical_Variable',[])
dcic = data[u,:]
acic=sm.bmu_ind_to_xy(sm.project_data(dcic))
x,y=np.zeros(len(acic)),np.zeros(len(acic))
k=0
for i in acic:
    x[k]=i[0]
    y[k]=i[1]
    k+=1
plt.plot(y,x,'rx',alpha=0.8)


plt.subplot(3,4,4)
plt.title('Galex_Variable',fontsize=15)
cf=plt.imshow(med_r,origin='lower',cmap='viridis')
plt.axis('off')
u = labc.get('Galex_Variable',[])
dcic = data[u,:]
acic=sm.bmu_ind_to_xy(sm.project_data(dcic))
x,y=np.zeros(len(acic)),np.zeros(len(acic))
k=0
for i in acic:
    x[k]=i[0]
    y[k]=i[1]
    k+=1
plt.plot(y,x,'rx',alpha=0.8)

plt.subplot(3,4,5)
plt.title('SPIDER_AGN',fontsize=15)
cf=plt.imshow(med_r,origin='lower',cmap='viridis')
plt.axis('off')
u = labc.get('SPIDER_AGN',[])
dcic = data[u,:]
acic=sm.bmu_ind_to_xy(sm.project_data(dcic))
x,y=np.zeros(len(acic)),np.zeros(len(acic))
k=0
for i in acic:
    x[k]=i[0]
    y[k]=i[1]
    k+=1
plt.plot(y,x,'rx',alpha=0.8)

plt.subplot(3,4,6)
plt.title('SPIDER_AGNBL',fontsize=15)
cf=plt.imshow(med_r,origin='lower',cmap='viridis')
plt.axis('off')
u = labc.get('SPIDER_AGNBL',[])
dcic = data[u,:]
acic=sm.bmu_ind_to_xy(sm.project_data(dcic))
x,y=np.zeros(len(acic)),np.zeros(len(acic))
k=0
for i in acic:
    x[k]=i[0]
    y[k]=i[1]
    k+=1
plt.plot(y,x,'rx',alpha=0.8)

plt.subplot(3,4,7)
plt.title('SPIDER_QSOBL',fontsize=15)
cf=plt.imshow(med_r,origin='lower',cmap='viridis')
plt.axis('off')
u = labc.get('SPIDER_QSOBL',[])
dcic = data[u,:]
acic=sm.bmu_ind_to_xy(sm.project_data(dcic))
x,y=np.zeros(len(acic)),np.zeros(len(acic))
k=0
for i in acic:
    x[k]=i[0]
    y[k]=i[1]
    k+=1
plt.plot(y,x,'rx',alpha=0.8)

plt.subplot(3,4,8)
plt.title('SPIDER_BL',fontsize=15)
cf=plt.imshow(med_r,origin='lower',cmap='viridis')
plt.axis('off')
u = labc.get('SPIDER_BL',[])
dcic = data[u,:]
acic=sm.bmu_ind_to_xy(sm.project_data(dcic))
x,y=np.zeros(len(acic)),np.zeros(len(acic))
k=0
for i in acic:
    x[k]=i[0]
    y[k]=i[1]
    k+=1
plt.plot(y,x,'rx',alpha=0.8)

plt.subplot(3,4,9)
plt.title('Turn-off',fontsize=15)
cf=plt.imshow(med_r,origin='lower',cmap='viridis')
plt.axis('off')
u = labc.get('Turn-off',[])
dcic = data[u,:]
acic=sm.bmu_ind_to_xy(sm.project_data(dcic))
x,y=np.zeros(len(acic)),np.zeros(len(acic))
k=0
for i in acic:
    x[k]=i[0]
    y[k]=i[1]
    k+=1
plt.plot(y,x,'rx',alpha=0.8)


plt.subplot(3,4,10)
plt.title('Turn-on',fontsize=15)
cf=plt.imshow(med_r,origin='lower',cmap='viridis')
plt.axis('off')
u = labc.get('Turn-on',[])
dcic = data[u,:]
acic=sm.bmu_ind_to_xy(sm.project_data(dcic))
x,y=np.zeros(len(acic)),np.zeros(len(acic))
k=0
for i in acic:
    x[k]=i[0]
    y[k]=i[1]
    k+=1
plt.plot(y,x,'rx',alpha=0.8)

plt.subplot(3,4,11)
plt.title('TDE',fontsize=15)
cf=plt.imshow(med_r,origin='lower',cmap='viridis')
plt.axis('off')
u = labc.get('TDE',[])
dcic = data[u,:]
acic=sm.bmu_ind_to_xy(sm.project_data(dcic))
x,y=np.zeros(len(acic)),np.zeros(len(acic))
k=0
for i in acic:
    x[k]=i[0]
    y[k]=i[1]
    k+=1
plt.plot(y,x,'rx',alpha=0.8)

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
# Create the figure
fig = plt.figure(figsize=(12, 8))
gs = gridspec.GridSpec(2, 4, height_ratios=[1, 1], width_ratios=[1,1,1,1])

# Create the subplots
ax0 = fig.add_subplot(gs[0:1, 0:2])  
ax1 = fig.add_subplot(gs[0:1, 2:4])
ax2 = fig.add_subplot(gs[1:2, 0:2])  
ax3 = fig.add_subplot(gs[1:2, 2:4])

objid = df_lc.index.get_level_values('objectid')[:].unique()

seen = Counter()
for (objectid, label), singleobj in df_lc.groupby(level=["objectid", "label"]):
    bitwise_sum = int(label)
    active_labels = translate_bitwise_sum_to_labels(bitwise_sum)
    #active_labels = translate_bitwise_sum_to_labels(label[0])
    seen.update(active_labels)
#changing order of labels in dictionary only for text to be readable on the plot
key_order = ('SDSS_QSO','SPIDER_BL','SPIDER_QSOBL', 'SPIDER_AGNBL',
             'WISE_Variable','Optical_Variable','Galex_Variable','Turn-on', 'Turn-off','TDE')
new_queue = OrderedDict()
for k in key_order:
    new_queue[k] = seen[k]
    
h = ax0.pie(new_queue.values(),labels=new_queue.keys(),autopct=autopct_format(new_queue.values()), textprops={'fontsize': 12},startangle=210,  labeldistance=1.1, wedgeprops = { 'linewidth' : 3, 'edgecolor' : 'white' }, colors=color4[2:])

seen2 = Counter()
for f in fzr:
    active_labels = translate_bitwise_sum_to_labels(int(f))
    seen2.update(active_labels)

new_queue = OrderedDict()
for k in key_order:
    new_queue[k] = seen2[k]


h = ax2.pie(new_queue.values(),labels=new_queue.keys(),autopct=autopct_format(new_queue.values()), textprops={'fontsize': 12},startangle=180,  labeldistance=1.1, wedgeprops = { 'linewidth' : 3, 'edgecolor' : 'white' }, colors=color4[2:])

 
#####################################################################################################
seen = Counter()
seen = df_lc.reset_index().groupby('band').objectid.nunique().to_dict()

cadence = dict((el,[]) for el in seen.keys())
timerange = dict((el,[]) for el in seen.keys())

for (_, band), times in df_lc.reset_index().groupby(["objectid", "band"]).time:
    cadence[band].append(len(times))
    if times.max() - times.min() > 0:
        timerange[band].append(np.round(times.max() - times.min(), 1))


i=0
colorlabel =[0,0,0,1,1,2,2,2,2,2,3,3,3]

for el in cadence.keys():
    #print(el,len(cadence[el]),np.mean(cadence[el]),np.std(cadence[el]))
    #print(el,len(timerange[el]),np.mean(timerange[el]),np.std(timerange[el]))
    ax1.scatter(np.mean(cadence[el]),np.mean(timerange[el]),s=seen[el],alpha=0.7,c=color4[colorlabel[i]],label=el)
    ax1.errorbar(np.mean(cadence[el]),np.mean(timerange[el]),label=el,yerr=np.std(timerange[el]),xerr=np.std(cadence[el]),alpha=0.2,c=color4[colorlabel[i]])

    i+=1
ax1.annotate('ZTF', # text to display
             (110, 1300),        # text location
             size=12, rotation=40 )
ax1.annotate('WISE', # text to display
             (20, 3800),        # text location
             size=12, rotation=40 )
ax1.annotate('GAIA', # text to display
             (25, 600),        # text location
             size=12, rotation=40 )
ax1.annotate('Pan-STARRS', # text to display
             (22, 1600),        # text location
             size=12, rotation=40 )

ax1.set_xlabel(r'$\rm Average\ number\ of\ visits$',size=15)
ax1.set_ylabel(r'$\rm Average\ baseline\ (days)$',size=15)
ax1.set_xscale('log')

#####################################################################################################
#ax3.hist(redshift_shuffled,label='final')
samp = pd.read_csv('data/AGNsample_26Feb24.csv')

for col in range(2,13):
    u = (samp.iloc[:, col]==1)
    ax3.hist(redshifts[u],histtype='step',label=samp.columns[col])
    
ax3.set_xlabel(r'$\rm Redshifts$',size=15)
ax3.set_ylabel(r'$\rm counts$',size=15)
plt.legend()
plt.tight_layout()
#plt.savefig('output/sample.png')
```

```{code-cell} ipython3
d = np.load('GP_ZTFWISE.npz',allow_pickle=True)
dd = d['data']
dl = d['labc']
```

```{code-cell} ipython3
mapper = umap.UMAP(n_neighbors=50,min_dist=0.9,metric='euclidean',random_state=20).fit(dd)
```

```{code-cell} ipython3
# Calculate 2D histogram
hist, x_edges, y_edges = np.histogram2d(mapper.embedding_[:, 0], mapper.embedding_[:, 1], bins=12)
plt.figure(figsize=(15,12))
i=1
ax0 = plt.subplot(4,4,12)
laborder = ['SDSS_QSO','WISE_Variable','Optical_Variable','Galex_Variable','SPIDER_AGN','SPIDER_AGNBL','SPIDER_QSOBL','SPIDER_BL','Turn-on','Turn-off','TDE']
for label in laborder:
    if label in labc:
        indices = labc[label]
        hist_per_cluster, _, _ = np.histogram2d(mapper.embedding_[indices,0], mapper.embedding_[indices,1], bins=(x_edges, y_edges))
        prob = hist_per_cluster / hist
        plt.subplot(4,4,i)
        plt.title(label)
        plt.contourf(x_edges[:-1], y_edges[:-1], prob.T, levels=20, alpha=0.8,cmap=custom_cmap)
        plt.colorbar()
        plt.axis('off')
        #cf = ax0.scatter(mapper.embedding_[indices,0],mapper.embedding_[indices,1],s=80,alpha=0.5,edgecolor='gray',label=label,c=colors[i-1])
        i+=1
ax0.legend(loc=4,fontsize=7)
ax0.axis('off')
plt.tight_layout()
#plt.savefig('output/umap2-ztfwise-gp.png')
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
# Create the figure
fig = plt.figure(figsize=(10, 10))
gs = gridspec.GridSpec(3, 3, height_ratios=[2, 2, 2], width_ratios=[1.5,1.5,2])

# Create the subplots
ax0 = fig.add_subplot(gs[0:2, 0:2])  # Top left, taking two-thirds of the width
ax1 = fig.add_subplot(gs[0:2, -1])  # Top right, taking one-third of the width
bottom_ax = fig.add_subplot(gs[2, :])      # Bottom, spanning the full width

objid = df_lc.index.get_level_values('objectid')[:].unique()
seen = Counter()


for b in objid:
    singleobj = df_lc.loc[b,:,:,:]
    label = singleobj.index.unique('label')
    # Translate the bitwise sum back to active labels
    bitwise_sum = int(label[0])  # Convert to integer
    active_labels = translate_bitwise_sum_to_labels(bitwise_sum)
    #active_labels = translate_bitwise_sum_to_labels(label[0])
    seen.update(active_labels)
#changing order of labels in dictionary only for text to be readable on the plot
key_order = ('SDSS_QSO','SPIDER_BL','SPIDER_QSOBL', 'SPIDER_AGNBL',
             'WISE_Variable','Optical_Variable','Galex_Variable','Turn-on', 'Turn-off','TDE','Fermi_Blazars')
new_queue = OrderedDict()
for k in key_order:
    new_queue[k] = seen[k]


h = ax0.pie(new_queue.values(),labels=new_queue.keys(),autopct=autopct_format(new_queue.values()), textprops={'fontsize': 12},startangle=210,  labeldistance=1.1, wedgeprops = { 'linewidth' : 3, 'edgecolor' : 'white' }, colors=colors3[:])


# Plot on the second subplot in the first row
seen = Counter()
for b in objid:
    singleobj = df_lc.loc[b,:,:,:]
    label = singleobj.index.unique('label')
    bands = singleobj.loc[label[0],:,:].index.get_level_values('band')[:].unique()
    seen.update(bands)
    
#####################################################################################################

cadence = dict((el,[]) for el in seen.keys())
timerange = dict((el,[]) for el in seen.keys())

for b in objid:
    singleobj = df_lc.loc[b,:,:,:]
    label = singleobj.index.unique('label')
    bband = singleobj.index.unique('band')
    for bb in bband:
        bands = singleobj.loc[label[0],bb,:].index.get_level_values('time')[:]
        #bands.values
        #print(bb,len(bands[:]),np.round(bands[:].max()-bands[:].min(),1))
        cadence[bb].append(len(bands[:]))
        if bands[:].max()-bands[:].min()>0:
            timerange[bb].append(np.round(bands[:].max()-bands[:].min(),1))

i=0
colorlabel =[0,1,2,3,3,3,4,4,4,4,4,2,5,5,5,6,7]
for el in cadence.keys():
    #print(el,len(cadence[el]),np.mean(cadence[el]),np.std(cadence[el]))
    #print(el,len(timerange[el]),np.mean(timerange[el]),np.std(timerange[el]))
    ax1.scatter(np.mean(cadence[el]),np.mean(timerange[el]),s=len(timerange[el]),alpha=0.7,c=colors3[colorlabel[i]+1])
    ax1.errorbar(np.mean(cadence[el]),np.mean(timerange[el]),label=el,yerr=np.std(timerange[el]),xerr=np.std(cadence[el]),alpha=0.2,c=colors3[colorlabel[i]+1])
    #print(el,np.mean(cadence[el]))

    i+=1
    
ax1.annotate('ZTF', # text to display
             (120, 1300),        # text location
             size=12, rotation=40 )
ax1.annotate('WISE', # text to display
             (20, 3800),        # text location
             size=12, rotation=40 )
ax1.annotate('GAIA', # text to display
             (15, 700),        # text location
             size=12, rotation=40 )
ax1.annotate('Pan-STARRS', # text to display
             (22, 1600),        # text location
             size=12, rotation=40 )
ax1.annotate('IceCube', # text to display
             (1, 1500),        # text location
             size=12, rotation=40 )
ax1.annotate('F814W', # text to display
             (6, 3000),        # text location
             size=12, rotation=40 )
ax1.annotate('Fermi', # text to display
             (1, 30),        # text location
             size=12, rotation=40 )
#ax1.annotate('GRB', # text to display
#             (1, 10),        # text location
#             size=12, rotation=40 )

#ax1.legend()
ax1.set_xlabel(r'$\rm Average\ number\ of\ visits$',size=15)
ax1.set_ylabel(r'$\rm Average\ baseline\ (days)$',size=15)
ax1.set_xscale('log')

#####################################################################################################
seen = Counter()
for b in objid:
    singleobj = df_lc.loc[b,:,:,:]
    label = singleobj.index.unique('label')
    bands = singleobj.loc[label[0],:,:].index.get_level_values('band')[:].unique()
    seen.update(bands)
    
key_order = ('IceCube','zg', 'zr', 'zi',  'G', 'BP', 'RP', 'panstarrs g',
             'panstarrs r', 'panstarrs i', 'panstarrs z', 'panstarrs y', 'F814W','W1','W2','FERMIGTRIG','SAXGRBMGRB')

#changing order of labels in dictionary only for text to be readable on the plot
new_queue = OrderedDict()
for k in key_order:
    new_queue[k] = seen[k]


tiklabels = ['IceCube','ZTF g','ZTF r','ZTF i','GAIA G', 'GAIA BP', 'GAIA RP',
             'Pan-STARRS g','Pan-STARRS r','Pan-STARRS i','Pan-STARRS z','Pan-STARRS y','F814W','WISE 1','WISE 2','Fermi','GRB']

h = bottom_ax.bar(range(len(tiklabels)), new_queue.values(),color= "#3F51B5")
bottom_ax.set_xticks(range(len(tiklabels)),tiklabels,fontsize=15,rotation=90)
bottom_ax.set_ylabel(r'$\rm number\ of\ lightcurves$',size=15)
plt.tight_layout()

#plt.savefig('output/sample.png')
```

```{code-cell} ipython3

```

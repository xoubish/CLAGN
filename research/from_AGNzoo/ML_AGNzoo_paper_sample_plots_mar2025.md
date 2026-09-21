---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.16.7
kernelspec:
  display_name: Python 3 (ipykernel)
  language: python
  name: python3
---

# How do AGNs selected with different techniques compare?

We use manifold learning and dimensionality reduction to learn the distribution of AGN lightcurves observed with different facilities. We mostly focus on UMAP ([Uniform Manifold Approximation and Projection, McInnes 2020](https://arxiv.org/pdf/1802.03426.pdf)). The reduced 2D projections from the unsupervised ML techniques reveal similarities and overlaps of different selection techniques and coloring the projections with various statistical physical properties (e.g., mean brightness, fractional lightcurve variation) is informative of correlations of the selections technique with physics such as AGN variability. Using different parts of the EM in training (or in building the initial higher dimensional manifold) demonstrates how much information if any is in that part of the data for each labeling scheme, for example whether with ZTF optical light curves alone, we can identify sources with variability in WISE near IR bands. These techniques also have a potential for identifying targets of a specific class or characteristic for future follow up observations.

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

color4 = ['#3182bd','#6baed6','#9ecae1','#e6550d','#fd8d3c','#fdd0a2','#31a354','#a1d99b', '#c7e9c0', '#756bb1', '#bcbddc', '#dadaeb', '#969696', '#bdbdbd','#d9d9d9','b','r','g']
custom_cmap = LinearSegmentedColormap.from_list("custom_theme", color4[1:])
```

```{code-cell} ipython3
samp = pd.read_csv('data/AGNsample_March7.csv')

df_lc = pd.read_parquet('data/Dave_df_march7.parquet')
objids = df_lc.index.get_level_values('objectid')[:].unique()
redshifts = samp['redshift']#[objids]
df_lc
```

```{code-cell} ipython3
unique_labels = df_lc.index.get_level_values('label').unique()
print("Unique labels:", unique_labels)

for u in unique_labels:
    print(translate_bitwise_sum_to_labels(u))
```

```{code-cell} ipython3
from plot_functions import create_figure
grouped = list(df_lc.groupby('objectid'))

#1526
for ind in range(2118,2119):
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
    if (len(x2) > 5) and not np.isnan(y2).any():
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
plt.text(20,0.92,'ZTF g,r',size=15)

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
    if (len(x2) > 5) and not np.isnan(y2).any():
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
plt.text(20,5.6,'WISE W1,W2',size=15)
plt.grid()
#plt.xlim([-10,1880])
plt.xlabel(r'$\rm time(day)$',size=15)
plt.ylabel(r'$\rm Flux(mJy)$',size=15)

plt.legend(loc=1)
plt.tight_layout()
#plt.savefig('output/unify_lc1619.png')
```

# ZTF bands only

```{code-cell} ipython3
bands_inlc = ['zg','zr','zi']
numobjs = len(df_lc.index.get_level_values('objectid')[:].unique())
sample_objids = df_lc.index.get_level_values('objectid').unique()[:numobjs]
objects,dobjects,flabels,zlist,keeps = unify_lc_gp_parallel(df_lc,redshifts,bands_inlc=bands_inlc,xres=120)

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
```

```{code-cell} ipython3
mapper = umap.UMAP(n_neighbors=100,min_dist=0.99,metric=dtw_distance,random_state=2).fit(data)
#mapper = umap.UMAP(n_neighbors=100,min_dist=0.99,metric='manhattan',random_state=3).fit(data)

plt.figure(figsize=(12,4))
markersize=100
cmap1 = 'viridis'

ax1 = plt.subplot(1,2,1)
ax1.set_title(r'$\rm Mean\ brightness$')
thiscolor=np.log10(np.nansum(average_arr,axis=0))
u = (thiscolor<2) & (thiscolor>=-2)
cf = ax1.scatter(mapper.embedding_[u,0],mapper.embedding_[u,1],c = thiscolor[u],s=markersize,edgecolor='k',cmap=cmap1)
plt.axis('off')
divider = make_axes_locatable(ax1)
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(cf,cax=cax)


ax1 = plt.subplot(1,2,2)
ax1.set_title(r'$\rm Mean\ Fractional\ Variation$')
thiscolor=stretch_small_values_arctan(np.nansum(fvar_arr,axis=0),factor=15)
u = (thiscolor<2.) & (thiscolor>=0)
cf = ax1.scatter(mapper.embedding_[u,0],mapper.embedding_[u,1],c = thiscolor[u],s=markersize,edgecolor='k',cmap=cmap1)
plt.axis('off')
divider = make_axes_locatable(ax1)
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(cf,cax=cax)


plt.tight_layout()
#plt.savefig('output/new_ztf.png')
```

```{code-cell} ipython3
# Initialize plot grid
fig, axes = plt.subplots(3, 4, figsize=(15, 9))  # 3 rows, 5 columns
axes = axes.flatten()  # Convert to a list for easier indexing

# Subset labels
laborder = ['SDSS_QSO',
            'bigMAC_binary','Graham15','Charisi16','Ward22_wise','Ward22_ztf',
            'Chen20','Liu19','PG1302','OJ287','bigMAC_dual','Rodriguez06']

# Compute background color map for full dataset
thiscolor = stretch_small_values_arctan(np.nansum(fvar_arr, axis=0), factor=3)
u = (thiscolor < 2) & (thiscolor >= 0)  # Apply mask

# Iterate over labels and plot
for i, label in enumerate(laborder):
    if label in labc:
        indices = labc[label]
        subset_transformed = mapper.transform(data[indices])  # Transform data

        # Select subplot
        ax = axes[i]
        ax.set_title(label)

        # Background density map
        cf = ax.scatter(mapper.embedding_[u, 0], mapper.embedding_[u, 1], 
                        c=thiscolor[u], s=10, edgecolor='none', cmap=cmap1, alpha=0.6)

        # Overlay subset scatter plot
        ax.scatter(subset_transformed[:, 0], subset_transformed[:, 1], 
                   s=50, alpha=0.8, edgecolor='black', color='red', label=label)

        # Hide axes for cleaner visualization
        ax.axis('off')

        # Add colorbar
        #divider = make_axes_locatable(ax)
        #cax = divider.append_axes("right", size="5%", pad=0.05)
        #fig.colorbar(cf, cax=cax)

# Adjust layout
plt.tight_layout()
plt.savefig('newway_ZTF.png')
```

```{code-cell} ipython3
# Assuming 'mapper.embedding_' is your data and 'labc' is your dictionary of labels to indices
hist, x_edges, y_edges = np.histogram2d(mapper.embedding_[:, 0], mapper.embedding_[:, 1], bins=11)
plt.figure(figsize=(12, 6))

# Define groups of labels
group_labels = {
    'SDSS QSOs': ['SDSS_QSO'],
    'MBHB Candidates': ['bigMAC_binary','OJ287','PG1302','Charisi16', 'Chen20', 'Graham15', 'Liu19', 'Ward22_wise', 'Ward22_ztf'],
}

# Custom colormap for visual consistency
custom_cmap = 'viridis'  # Replace with your colormap of choice

# Create subplots for each group
i = 1
for group_name, labels in group_labels.items():
    group_indices = np.hstack([labc[label] for label in labels if label in labc])
    hist_per_group, _, _ = np.histogram2d(mapper.embedding_[group_indices, 0], mapper.embedding_[group_indices, 1], bins=(x_edges, y_edges))
    prob = hist_per_group / hist
    plt.subplot(2, 3, i)  # Adjust the subplot layout as needed
    plt.title(group_name)
    plt.contourf(x_edges[:-1], y_edges[:-1], prob.T, levels=15, alpha=0.8, cmap=custom_cmap)
    plt.colorbar()
    plt.axis('off')
    i += 1

plt.scatter(-13,31,marker='x',color='r',label='PG 1302-2102')
plt.scatter(-6.5,29.5,marker='x',color='r',label='OJ 287')
plt.legend(loc=4)

group_labels = {
    'Confirmed Dual-MBHs': ['bigMAC_dual', 'Rodriguez06']#
}

# Custom colormap for visual consistency
custom_cmap = 'viridis'  # Replace with your colormap of choice

# Create subplots for each group
i = 1
for group_name, labels in group_labels.items():
    group_indices = np.hstack([labc[label] for label in labels if label in labc])
    hist_per_group, _, _ = np.histogram2d(mapper.embedding_[group_indices, 0], mapper.embedding_[group_indices, 1], bins=(x_edges, y_edges))
    prob = hist_per_group / hist
    plt.subplot(2, 3, 3)  # Adjust the subplot layout as needed
    plt.title(group_name)
    plt.contourf(x_edges[:-1], y_edges[:-1], prob.T, levels=12, alpha=0.8, cmap=custom_cmap)
    plt.colorbar()
    plt.axis('off')
    i += 1

plt.scatter(1,28.5,marker='x',color='r',label='0402+379')
plt.legend(loc=4)

# Define groups of labels

plt.tight_layout()
plt.savefig('ztf_hist2d.png')
```

# Wise only

```{code-cell} ipython3
bands_inlc = ['W1']
numobjs = len(df_lc.index.get_level_values('objectid')[:].unique())
sample_objids = df_lc.index.get_level_values('objectid').unique()[:numobjs]
df_lc_small = df_lc.loc[sample_objids]
objects,dobjects,flabels,zlist,keeps = unify_lc_gp_parallel(df_lc_small,redshifts,bands_inlc=bands_inlc,xres=80)

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
mapper = umap.UMAP(n_neighbors=100,min_dist=0.99,metric=dtw_distance,random_state=7).fit(data)
#mapper = umap.UMAP(n_neighbors=7,min_dist=0.999,metric='manhattan',random_state=5).fit(data)

plt.figure(figsize=(12,4))
markersize=100
cmap1 = 'viridis'

ax1 = plt.subplot(1,2,1)
ax1.set_title(r'$\rm Mean\ brightness$')
thiscolor=np.log10(np.nansum(average_arr,axis=0))
u = (thiscolor<2) & (thiscolor>=-2)
cf = ax1.scatter(mapper.embedding_[u,0],mapper.embedding_[u,1],c = thiscolor[u],s=markersize,edgecolor='k',cmap=cmap1)
plt.axis('off')
divider = make_axes_locatable(ax1)
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(cf,cax=cax)


ax1 = plt.subplot(1,2,2)
ax1.set_title(r'$\rm Mean\ Fractional\ Variation$')
thiscolor=stretch_small_values_arctan(np.nansum(fvar_arr,axis=0),factor=15)
u = (thiscolor<2.) & (thiscolor>=0)
cf = ax1.scatter(mapper.embedding_[u,0],mapper.embedding_[u,1],c = thiscolor[u],s=markersize,edgecolor='k',cmap=cmap1)
plt.axis('off')
divider = make_axes_locatable(ax1)
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(cf,cax=cax)


plt.tight_layout()
#plt.savefig('output/umap-w1-sampleA-1.png')
```

```{code-cell} ipython3
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

# Initialize plot grid
fig, axes = plt.subplots(3, 4, figsize=(15, 9))  # 3 rows, 5 columns
axes = axes.flatten()  # Convert to a list for easier indexing

# Subset labels
laborder = ['SDSS_QSO',
            'bigMAC_binary','Graham15','Charisi16','Ward22_wise','Ward22_ztf',
            'Chen20','Liu19','PG1302','OJ287','bigMAC_dual','Rodriguez06']

# Compute background color map for full dataset
thiscolor = stretch_small_values_arctan(np.nansum(fvar_arr, axis=0), factor=3)
u = (thiscolor < 1.) & (thiscolor >= 0)  # Apply mask

# Iterate over labels and plot
for i, label in enumerate(laborder):
    if label in labc:
        indices = labc[label]
        subset_transformed = mapper.transform(data[indices])  # Transform data

        # Select subplot
        ax = axes[i]
        ax.set_title(label)

        # Background density map
        cf = ax.scatter(mapper.embedding_[u, 0], mapper.embedding_[u, 1], 
                        c=thiscolor[u], s=30, edgecolor='none', cmap=cmap1, alpha=0.6)

        # Overlay subset scatter plot
        ax.scatter(subset_transformed[:, 0], subset_transformed[:, 1], 
                   s=50, alpha=0.8, edgecolor='black', color='red', label=label)

        # Hide axes for cleaner visualization
        ax.axis('off')

        # Add colorbar
        #divider = make_axes_locatable(ax)
        #cax = divider.append_axes("right", size="5%", pad=0.05)
        #fig.colorbar(cf, cax=cax)

# Adjust layout
plt.tight_layout()
plt.savefig('newway_wise.png')
```

```{code-cell} ipython3
# Assuming 'mapper.embedding_' is your data and 'labc' is your dictionary of labels to indices
hist, x_edges, y_edges = np.histogram2d(mapper.embedding_[:, 0], mapper.embedding_[:, 1], bins=12)
plt.figure(figsize=(12, 6))

# Define groups of labels
group_labels = {
    'SDSS QSOs': ['SDSS_QSO'],
    'MBHB Candidates': ['OJ287','PG1302','Charisi16', 'Chen20', 'Graham15', 'Liu19', 'Ward22_wise', 'Ward22_ztf'],#'bigMAC_binary',
}

# Custom colormap for visual consistency
custom_cmap = 'viridis'  # Replace with your colormap of choice

# Create subplots for each group
i = 1
for group_name, labels in group_labels.items():
    group_indices = np.hstack([labc[label] for label in labels if label in labc])
    hist_per_group, _, _ = np.histogram2d(mapper.embedding_[group_indices, 0], mapper.embedding_[group_indices, 1], bins=(x_edges, y_edges))
    prob = hist_per_group / hist
    plt.subplot(2, 3, i)  # Adjust the subplot layout as needed
    plt.title(group_name)
    plt.contourf(x_edges[:-1], y_edges[:-1], prob.T, levels=15, alpha=0.8, cmap=custom_cmap)
    plt.colorbar()
    plt.axis('off')
    i += 1

plt.scatter(5.5,19,marker='+',color='r',label='PG 1302-2102')
plt.scatter(8.5,16.5,marker='x',color='r',label='OJ 287')
plt.legend(loc=4)

group_labels = {
    'Confirmed Dual-MBHs': ['bigMAC_dual', 'Rodriguez06']#
}

# Custom colormap for visual consistency
custom_cmap = 'viridis'  # Replace with your colormap of choice

# Create subplots for each group
i = 1
for group_name, labels in group_labels.items():
    group_indices = np.hstack([labc[label] for label in labels if label in labc])
    hist_per_group, _, _ = np.histogram2d(mapper.embedding_[group_indices, 0], mapper.embedding_[group_indices, 1], bins=(x_edges, y_edges))
    prob = hist_per_group / hist
    plt.subplot(2, 3, 3)  # Adjust the subplot layout as needed
    plt.title(group_name)
    plt.contourf(x_edges[:-1], y_edges[:-1], prob.T, levels=12, alpha=0.8, cmap=custom_cmap)
    plt.colorbar()
    plt.axis('off')
    i += 1

plt.scatter(15.5,23.,marker='x',color='r',label='0402+379')
plt.legend(loc=4)

# Define groups of labels

plt.tight_layout()
plt.savefig('WISE_hist2d.png')
```

```{code-cell} ipython3
# Assuming 'mapper.embedding_' is your data and 'labc' is your dictionary of labels to indices
hist, x_edges, y_edges = np.histogram2d(mapper2.embedding_[:, 0], mapper2.embedding_[:, 1], bins=12)
plt.figure(figsize=(12, 6))

# Define groups of labels
group_labels = {
    'SDSS QSOs': ['SDSS_QSO'],
    'MBHB Candidates': ['OJ287','PG1302','Charisi16', 'Chen20', 'Graham15', 'Liu19', 'Ward22_wise', 'Ward22_ztf'],
}

# Custom colormap for visual consistency
custom_cmap = 'viridis'  # Replace with your colormap of choice
thiscolor = stretch_small_values_arctan(np.nansum(fvar_arr, axis=0), factor=3)


# Create subplots for each group
i = 1
for group_name, labels in group_labels.items():
    group_indices = np.hstack([labc[label] for label in labels if label in labc])
    thcol = thiscolor[group_indices]
    subset_transformed = mapper.transform(data[group_indices])  # Transform data
    plt.subplot(2, 3, i)  # Adjust the subplot layout as needed
    plt.title(group_name)
    plt.scatter(mapper.embedding_[:, 0], mapper.embedding_[:, 1],c='grey',alpha=0.3)
    plt.scatter(subset_transformed[:, 0], subset_transformed[:, 1], c=thcol, s=30, edgecolor='none', cmap=cmap1, alpha=0.6)
    plt.colorbar()
    plt.axis('off')
    i += 1


group_indices = np.hstack([labc['PG1302']])
subset_transformed = mapper.transform(data[group_indices])  # Transform data
plt.scatter(subset_transformed[:, 0], subset_transformed[:, 1],marker='+',color='r',label='PG 1302-2102')

group_indices = np.hstack([labc['OJ287']])
subset_transformed = mapper.transform(data[group_indices])  # Transform data
plt.scatter(subset_transformed[:, 0], subset_transformed[:, 1],marker='x',color='r',label='OJ 287')
plt.legend(loc=4)

group_labels = {
    'Confirmed Dual-MBHs': ['bigMAC_dual', 'Rodriguez06']#
}
for group_name, labels in group_labels.items():
    group_indices = np.hstack([labc[label] for label in labels if label in labc])
    thcol = thiscolor[group_indices]
    subset_transformed = mapper.transform(data[group_indices])  # Transform data
    plt.subplot(2, 3, 3)  # Adjust the subplot layout as needed
    plt.title(group_name)
    plt.scatter(mapper.embedding_[:, 0], mapper.embedding_[:, 1],c='grey',alpha=0.3)
    plt.scatter(subset_transformed[:, 0], subset_transformed[:, 1], c=thcol, s=30, edgecolor='none', cmap=cmap1, alpha=0.6)
    plt.colorbar()
    plt.axis('off')

group_indices = np.hstack([labc['Rodriguez06']])
subset_transformed = mapper.transform(data[group_indices])  # Transform data
plt.scatter(subset_transformed[:, 0], subset_transformed[:, 1],marker='x',color='r',label='0402+379')
plt.legend(loc=4)

# Define groups of labels

plt.tight_layout()
plt.savefig('WISE_last.png')
```

```{code-cell} ipython3
# Assuming 'mapper.embedding_' is your data and 'labc' is your dictionary of labels to indices
hist, x_edges, y_edges = np.histogram2d(mapper2.embedding_[:, 0], mapper2.embedding_[:, 1], bins=12)
plt.figure(figsize=(10, 4))

# Define groups of labels
group_labels = {
    'MBHB Candidates': ['OJ287','PG1302','Charisi16', 'Chen20', 'Graham15', 'Liu19', 'Ward22_wise', 'Ward22_ztf'],
}

# Custom colormap for visual consistency
custom_cmap = 'viridis'  # Replace with your colormap of choice
thiscolor = stretch_small_values_arctan(np.nansum(fvar_arr, axis=0), factor=3)


# Create subplots for each group
i = 1
for group_name, labels in group_labels.items():
    group_indices = np.hstack([labc[label] for label in labels if label in labc])
    thcol = thiscolor[group_indices]
    subset_transformed = mapper.transform(data[group_indices])  # Transform data
    plt.subplot(1, 2, 1)  # Adjust the subplot layout as needed
    plt.title(group_name)
    plt.scatter(mapper.embedding_[:, 0], mapper.embedding_[:, 1],c='grey',alpha=0.3)
    plt.scatter(subset_transformed[:, 0], subset_transformed[:, 1], c=thcol, s=30, edgecolor='none', cmap=cmap1)
    plt.colorbar()
    plt.axis('off')
    i += 1


group_indices = np.hstack([labc['PG1302']])
subset_transformed = mapper.transform(data[group_indices])  # Transform data
plt.scatter(subset_transformed[:, 0], subset_transformed[:, 1],marker='+',color='r',label='PG 1302-2102')

group_indices = np.hstack([labc['OJ287']])
subset_transformed = mapper.transform(data[group_indices])  # Transform data
plt.scatter(subset_transformed[:, 0], subset_transformed[:, 1],marker='x',color='r',label='OJ 287')
plt.legend(loc=3)

group_labels = {
    'Confirmed Dual-MBHs': ['bigMAC_dual', 'Rodriguez06']#
}
for group_name, labels in group_labels.items():
    group_indices = np.hstack([labc[label] for label in labels if label in labc])
    thcol = thiscolor[group_indices]
    subset_transformed = mapper.transform(data[group_indices])  # Transform data
    plt.subplot(1, 2, 2)  # Adjust the subplot layout as needed
    plt.title(group_name)
    plt.scatter(mapper.embedding_[:, 0], mapper.embedding_[:, 1],c='grey',alpha=0.3)
    plt.scatter(subset_transformed[:, 0], subset_transformed[:, 1], c=thcol, s=30, edgecolor='none', cmap=cmap1)
    plt.colorbar()
    plt.axis('off')

group_indices = np.hstack([labc['Rodriguez06']])
subset_transformed = mapper.transform(data[group_indices])  # Transform data
plt.scatter(subset_transformed[:, 0], subset_transformed[:, 1],marker='x',color='r',label='0402+379')
plt.legend(loc=3)

# Define groups of labels

plt.tight_layout()
plt.savefig('WISE_last.png')
```

# ZTF+ WISE Manifold

```{code-cell} ipython3
bands_inlc = ['zg','zr','zi','W1','W2']
numobjs = len(df_lc.index.get_level_values('objectid')[:].unique())
sample_objids = df_lc.index.get_level_values('objectid').unique()[:numobjs]
df_lc_small = df_lc.loc[sample_objids]
objects,dobjects,flabels,zlist,keeps = unify_lc_gp_parallel(df_lc_small,redshifts,bands_inlc=bands_inlc,xres=80)

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
#mapper2 = umap.UMAP(n_neighbors=50,min_dist=0.99,metric=dtw_distance,random_state=10).fit(data)
mapper2 = umap.UMAP(n_neighbors=17,min_dist=0.99,metric='manhattan',random_state=5).fit(data)

plt.figure(figsize=(12,4))
markersize=100
cmap1 = 'viridis'

ax1 = plt.subplot(1,3,1)
ax1.set_title(r'$\rm Mean\ brightness$')
thiscolor=np.log10(np.nansum(average_arr,axis=0))
u = (thiscolor<2) & (thiscolor>=-2)
cf = ax1.scatter(mapper2.embedding_[u,0],mapper2.embedding_[u,1],c = thiscolor[u],s=markersize,edgecolor='k',cmap=cmap1)
plt.axis('off')
divider = make_axes_locatable(ax1)
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(cf,cax=cax)


ax1 = plt.subplot(1,3,3)
ax1.set_title(r'$\rm Mean\ Fractional\ Variation$')
thiscolor=stretch_small_values_arctan(np.nansum(fvar_arr,axis=0),factor=3)
u = (thiscolor<1.5) & (thiscolor>=0)
cf = ax1.scatter(mapper2.embedding_[u,0],mapper2.embedding_[u,1],c = thiscolor[u],s=markersize,edgecolor='k',cmap=cmap1)
plt.axis('off')
divider = make_axes_locatable(ax1)
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(cf,cax=cax)


ax1 = plt.subplot(1,3,2)
ax1.set_title(r'$\rm Redshift$')
thiscolor=redshift_shuffled
u = (thiscolor<2) & (thiscolor>=0)
cf = ax1.scatter(mapper2.embedding_[u,0],mapper2.embedding_[u,1],c = thiscolor[u],s=markersize,edgecolor='k',cmap=cmap1)
plt.axis('off')
divider = make_axes_locatable(ax1)
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(cf,cax=cax)

plt.tight_layout()
```

## fig, axes = plt.subplots(3, 4, figsize=(15, 9))  # 3 rows, 5 columns
axes = axes.flatten()  # Convert to a list for easier indexing

# Subset labels
laborder = ['SDSS_QSO',
            'bigMAC_binary','Graham15','Charisi16','Ward22_wise','Ward22_ztf',
            'Chen20','Liu19','PG1302','OJ287','bigMAC_dual','Rodriguez06']

# Compute background color map for full dataset
thiscolor = stretch_small_values_arctan(np.nansum(fvar_arr, axis=0), factor=3)
u = (thiscolor < 1.) & (thiscolor >= 0)  # Apply mask

# Iterate over labels and plot
for i, label in enumerate(laborder):
    if label in labc:
        indices = labc[label]
        subset_transformed = mapper2.transform(data[indices])  # Transform data

        # Select subplot
        ax = axes[i]
        ax.set_title(label)

        # Background density map
        cf = ax.scatter(mapper2.embedding_[u, 0], mapper2.embedding_[u, 1], 
                        c=thiscolor[u], s=10, edgecolor='none', cmap=cmap1, alpha=0.6)

        # Overlay subset scatter plot
        ax.scatter(subset_transformed[:, 0], subset_transformed[:, 1], 
                   s=50, alpha=0.8, edgecolor='black', color='red', label=label)

        # Hide axes for cleaner visualization
        ax.axis('off')

        # Add colorbar
        #divider = make_axes_locatable(ax)
        #cax = divider.append_axes("right", size="5%", pad=0.05)
        #fig.colorbar(cf, cax=cax)

# Adjust layout
plt.tight_layout()
plt.savefig('newway_ztfwise.png')

```{code-cell} ipython3
# Assuming 'mapper.embedding_' is your data and 'labc' is your dictionary of labels to indices
hist, x_edges, y_edges = np.histogram2d(mapper2.embedding_[:, 0], mapper2.embedding_[:, 1], bins=12)
plt.figure(figsize=(12, 6))

# Define groups of labels
group_labels = {
    'SDSS QSOs': ['SDSS_QSO'],
    'MBHB Candidates': ['OJ287','PG1302','Charisi16', 'Chen20', 'Graham15', 'Liu19', 'Ward22_wise', 'Ward22_ztf'],
}

# Custom colormap for visual consistency
custom_cmap = 'viridis'  # Replace with your colormap of choice

# Create subplots for each group
i = 1
for group_name, labels in group_labels.items():
    group_indices = np.hstack([labc[label] for label in labels if label in labc])
    hist_per_group, _, _ = np.histogram2d(mapper2.embedding_[group_indices, 0], mapper2.embedding_[group_indices, 1], bins=(x_edges, y_edges))
    prob = hist_per_group / hist
    plt.subplot(2, 3, i)  # Adjust the subplot layout as needed
    plt.title(group_name)
    plt.contourf(x_edges[:-1], y_edges[:-1], prob.T, levels=15, alpha=0.8, cmap=custom_cmap)
    plt.colorbar()
    plt.axis('off')
    i += 1

plt.scatter(6.8,3,marker='+',color='r',label='PG 1302-2102')
plt.scatter(-5,11,marker='x',color='r',label='OJ 287')
plt.legend(loc=4)

group_labels = {
    'Confirmed Dual-MBHs': ['bigMAC_dual', 'Rodriguez06']#
}

# Custom colormap for visual consistency
custom_cmap = 'viridis'  # Replace with your colormap of choice

# Create subplots for each group
i = 1
for group_name, labels in group_labels.items():
    group_indices = np.hstack([labc[label] for label in labels if label in labc])
    hist_per_group, _, _ = np.histogram2d(mapper2.embedding_[group_indices, 0], mapper2.embedding_[group_indices, 1], bins=(x_edges, y_edges))
    prob = hist_per_group / hist
    plt.subplot(2, 3, 3)  # Adjust the subplot layout as needed
    plt.title(group_name)
    plt.contourf(x_edges[:-1], y_edges[:-1], prob.T, levels=12, alpha=0.8, cmap=custom_cmap)
    plt.colorbar()
    plt.axis('off')
    i += 1

plt.scatter(0.2,-3.5,marker='x',color='r',label='0402+379')
plt.legend(loc=4)

# Define groups of labels

plt.tight_layout()
plt.savefig('ZTFWISE_hist2d.png')
```

```{code-cell} ipython3
# Assuming 'mapper.embedding_' is your data and 'labc' is your dictionary of labels to indices
hist, x_edges, y_edges = np.histogram2d(mapper2.embedding_[:, 0], mapper2.embedding_[:, 1], bins=12)
plt.figure(figsize=(12, 6))

# Define groups of labels
group_labels = {
    'SDSS QSOs': ['SDSS_QSO'],
    'MBHB Candidates': ['OJ287','PG1302','Charisi16', 'Chen20', 'Graham15', 'Liu19', 'Ward22_wise', 'Ward22_ztf'],
}

# Custom colormap for visual consistency
custom_cmap = 'viridis'  # Replace with your colormap of choice
thiscolor = stretch_small_values_arctan(np.nansum(fvar_arr, axis=0), factor=3)


# Create subplots for each group
i = 1
for group_name, labels in group_labels.items():
    group_indices = np.hstack([labc[label] for label in labels if label in labc])
    thcol = thiscolor[group_indices]

    subset_transformed = mapper2.transform(data[group_indices])  # Transform data
    plt.subplot(2, 3, i)  # Adjust the subplot layout as needed
    plt.title(group_name)
    plt.scatter(mapper2.embedding_[:, 0], mapper2.embedding_[:, 1],c='grey',alpha=0.3)
    plt.scatter(subset_transformed[:, 0], subset_transformed[:, 1], c=thcol, s=30, edgecolor='none', cmap=cmap1, alpha=0.6)
    plt.colorbar()
    plt.axis('off')
    i += 1


group_indices = np.hstack([labc['PG1302']])
subset_transformed = mapper2.transform(data[group_indices])  # Transform data
plt.scatter(subset_transformed[:, 0], subset_transformed[:, 1],marker='+',color='r',label='PG 1302-2102')

group_indices = np.hstack([labc['OJ287']])
subset_transformed = mapper2.transform(data[group_indices])  # Transform data
plt.scatter(subset_transformed[:, 0], subset_transformed[:, 1],marker='x',color='r',label='OJ 287')
plt.legend(loc=4)

group_labels = {
    'Confirmed Dual-MBHs': ['bigMAC_dual', 'Rodriguez06']#
}

# Custom colormap for visual consistency
custom_cmap = 'viridis'  # Replace with your colormap of choice

# Create subplots for each group
i = 1
for group_name, labels in group_labels.items():

    group_indices = np.hstack([labc[label] for label in labels if label in labc])
    thcol = thiscolor[group_indices]
    subset_transformed = mapper2.transform(data[group_indices])  # Transform data
    plt.subplot(2, 3, 3)  # Adjust the subplot layout as needed
    plt.title(group_name)
    plt.scatter(mapper2.embedding_[:, 0], mapper2.embedding_[:, 1],c='grey',alpha=0.3)
    plt.scatter(subset_transformed[:, 0], subset_transformed[:, 1], c=thcol, s=30, edgecolor='none', cmap=cmap1, alpha=0.6)
    plt.colorbar()
    plt.axis('off')
    i += 1

group_indices = np.hstack([labc['Rodriguez06']])
subset_transformed = mapper2.transform(data[group_indices])  # Transform data
plt.scatter(subset_transformed[:, 0], subset_transformed[:, 1],marker='x',color='r',label='0402+379')
plt.legend(loc=4)

# Define groups of labels

plt.tight_layout()
plt.savefig('ZTFWISE_hist2d.png')
```

```{code-cell} ipython3
# Assuming 'mapper.embedding_' is your data and 'labc' is your dictionary of labels to indices
hist, x_edges, y_edges = np.histogram2d(mapper2.embedding_[:, 0], mapper2.embedding_[:, 1], bins=12)
plt.figure(figsize=(10, 4))

# Define groups of labels
group_labels = {
    'MBHB Candidates': ['OJ287','PG1302','Charisi16', 'Chen20', 'Graham15', 'Liu19', 'Ward22_wise', 'Ward22_ztf'],
}

# Custom colormap for visual consistency
custom_cmap = 'viridis'  # Replace with your colormap of choice
thiscolor = stretch_small_values_arctan(np.nansum(fvar_arr, axis=0), factor=3)


# Create subplots for each group
i = 1
for group_name, labels in group_labels.items():
    group_indices = np.hstack([labc[label] for label in labels if label in labc])
    thcol = thiscolor[group_indices]

    subset_transformed = mapper2.transform(data[group_indices])  # Transform data
    plt.subplot(1, 2, 1)  # Adjust the subplot layout as needed
    plt.title(group_name)
    plt.scatter(mapper2.embedding_[:, 0], mapper2.embedding_[:, 1],c='grey',alpha=0.3)
    plt.scatter(subset_transformed[:, 0], subset_transformed[:, 1], c=thcol, s=30, edgecolor='none', cmap=cmap1)
    plt.colorbar()
    plt.axis('off')
    i += 1


group_indices = np.hstack([labc['PG1302']])
subset_transformed = mapper2.transform(data[group_indices])  # Transform data
plt.scatter(subset_transformed[:, 0], subset_transformed[:, 1],marker='+',color='r',label='PG 1302-2102')

group_indices = np.hstack([labc['OJ287']])
subset_transformed = mapper2.transform(data[group_indices])  # Transform data
plt.scatter(subset_transformed[:, 0], subset_transformed[:, 1],marker='x',color='r',label='OJ 287')
plt.legend(loc=4)

group_labels = {
    'Confirmed Dual-MBHs': ['bigMAC_dual', 'Rodriguez06']#
}

# Custom colormap for visual consistency
custom_cmap = 'viridis'  # Replace with your colormap of choice

# Create subplots for each group
i = 1
for group_name, labels in group_labels.items():

    group_indices = np.hstack([labc[label] for label in labels if label in labc])
    thcol = thiscolor[group_indices]
    subset_transformed = mapper2.transform(data[group_indices])  # Transform data
    plt.subplot(1, 2, 2)  # Adjust the subplot layout as needed
    plt.title(group_name)
    plt.scatter(mapper2.embedding_[:, 0], mapper2.embedding_[:, 1],c='grey',alpha=0.3)
    plt.scatter(subset_transformed[:, 0], subset_transformed[:, 1], c=thcol, s=30, edgecolor='none', cmap=cmap1)
    plt.colorbar()
    plt.axis('off')
    i += 1

group_indices = np.hstack([labc['Rodriguez06']])
subset_transformed = mapper2.transform(data[group_indices])  # Transform data
plt.scatter(subset_transformed[:, 0], subset_transformed[:, 1],marker='x',color='r',label='0402+379')
plt.legend(loc=4)

# Define groups of labels

plt.tight_layout()
plt.savefig('ZTFWISE_hist2d.png')
```

```{code-cell} ipython3

```

---
jupyter:
  jupytext:
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.15.2
  kernelspec:
    display_name: Python 3 (ipykernel)
    language: python
    name: python3
---

# Build AGN final catalog sample 
By Shooby, Last edit Feb 26th

```python
import sys
sys.path.append('code_src/')

import os
import time
import astropy.units as u
import pandas as pd
import numpy as np
import csv
import astropy.io.fits as fits
from astropy.table import Table
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt
from data_structures import MultiIndexDFObject
from astroquery.sdss import SDSS
from astroquery.ipac.ned import Ned
from astroquery.vizier import Vizier
from astroquery.simbad import Simbad


import logging

# Get the root logger
logger = logging.getLogger()
logger.setLevel(logging.ERROR)

zmin = 0.0
zmax = 1.0

# Initialize agnlabels
#agnlabels = ['SDSS_QSO', 'WISE_Variable','Optical_Variable','Galex_Variable',
#             'Turn-on', 'Turn-off',
#             'SPIDER', 'SPIDER_AGN','SPIDER_BL','SPIDER_QSOBL','SPIDER_AGNBL', 
#             'TDE','BOSS SF','DualAGN']

agnlabels=['SDSS_QSO','WISE_Variable','Optical_Variable','Charisi16','Chen20','Graham15','Liu19','Ward22_wise','Ward22_ztf','bigMAC','Rodriguez06']

# Create an empty pandas DataFrame
columns = ['SkyCoord', 'redshift'] + agnlabels
df = pd.DataFrame(columns=columns)
# Initialize label columns to 0
for label in agnlabels:
    df[label] = 0
```

```python
# Function to check if a coordinate is close to any existing ones
def is_close(new_coord, threshold_arcsec=1):
    global df
    if df.empty:
        return False, None
    existing_coords = SkyCoord(df['SkyCoord'].tolist())
    sep = new_coord.separation(existing_coords)
    close_idx = sep < threshold_arcsec * u.arcsec
    if close_idx.any():
        return True, df[close_idx].index[0]
    else:
        return False, None

def update_or_append_multiple(ras, decs, redshifts, labels):
    global df  # Make sure df is recognized as the global DataFrame
    new_rows = []  # Prepare a list to collect new rows
    for ra, dec, redshift, label in zip(ras, decs, redshifts, labels):
        new_coord = SkyCoord(ra, dec, frame='icrs', unit='deg')
        exists, idx = is_close(new_coord)

        if exists:
            # Update existing row
            df.at[idx, label] = 1
            df.at[idx, 'redshift'] = redshift
        else:
            # Prepare a new row as a DataFrame instead of Series for consistency with pd.concat
            new_row = pd.DataFrame([{**{'SkyCoord': new_coord, 'redshift': redshift}, **{l: 1 if l == label else 0 for l in agnlabels}}])
            new_rows.append(new_row)

    # Append all new rows at once using pd.concat if there are any
    if new_rows:
        df = pd.concat([df, *new_rows], ignore_index=True)

```

## Add SDSS QSO from DR16

```python
# have to change order otherwise all low redshift
#query = f"SELECT TOP {num} specObjID, ra, dec, z FROM SpecObj \
#WHERE ( z > 0.001 AND z < {zmax} AND class='QSO' AND zWARNING=0 ) \
#ORDER BY NEWID()"
zmax = 1.0  # Example maximum redshift
num = 10000  # Total number of samples desired
num_bins = 20  # Number of bins
num_per_bin = num // num_bins  # Samples per bin

# Calculate bin edges
bin_edges = np.linspace(0.0001, zmax, num_bins + 1)

# Placeholder for results
results = pd.DataFrame()

for i in range(num_bins):
    bin_start = bin_edges[i]
    bin_end = bin_edges[i + 1]
    
    # Construct query for this bin
    query = f"""
    SELECT TOP {num_per_bin} specObjID, ra, dec, z
    FROM SpecObj
    WHERE (z > {bin_start} AND z < {bin_end} AND class='QSO' AND zWARNING=0)
    ORDER BY NEWID()
    """

    if num>0:
        res = SDSS.query_sql(query, data_release = 16)
        for r in res:
            update_or_append_multiple([r['ra']], [r['dec']], [r['z']], ['SDSS_QSO'])
print('SDSS QSO sources added: ',len(df))
hwy = plt.hist(df['redshift'], bins=30, alpha=0.75, color='blue')  # Adjust 'z' if your column name differs

```

# Variable AGNs (WISE, COSMOS VLT/Palomar, Galex)

```python
VAGN = pd.read_csv('data/WISE_MIR_variable_AGN_with_PS1_photometry_and_SDSS_redshift.csv')
uwise = (VAGN['SDSS_redshift']>0.001)&(VAGN['SDSS_redshift']<zmax)
vagn_labels = ['WISE_Variable' for ra in VAGN['SDSS_RA'][uwise]]
update_or_append_multiple(VAGN['SDSS_RA'][uwise],VAGN['SDSS_Dec'][uwise],VAGN['SDSS_redshift'][uwise],vagn_labels)
print('WISE Variable sources: ',len(vagn_labels))
```

```python
paper = Ned.query_refcode('2019A&A...627A..33D') #optically variable AGN in cosmos
up = (paper['Redshift']>0)&(paper['Redshift']<=zmax)
paper_labels = ['Optical_Variable' for ra in paper['RA'][up]]
update_or_append_multiple(paper['RA'][up], paper['DEC'][up],paper['Redshift'][up],paper_labels)
print('COSMOS VLT optical variable sources: ',len(paper_labels))

paper = Ned.query_refcode('2020ApJ...896...10B') #Palomar Variable
up = (paper['Redshift']>0)&(paper['Redshift']<=zmax)
paper_labels = ['Optical_Variable' for ra in paper['RA'][up]]
update_or_append_multiple(paper['RA'][up], paper['DEC'][up],paper['Redshift'][up],paper_labels)
print('Palomar variable sources added: ',len(paper_labels))
```

# Adding Dual AGNs from different sources

```python
#'Charisi16','Chen20','Graham15','Liu19','Ward22_wise','Ward22_ztf','bigMAC','Rodriguez06']

from astropy.coordinates import Angle

# Load data
VAGN = pd.read_csv('data/BigMAC_maintable_DR0p9.csv')
uwise = (VAGN['z1'] > 0.001) & (VAGN['z1'] < zmax) & (VAGN['ST1 Confidence Flag'] == 1)
ra_deg = [Angle(ra, unit='hourangle').deg for ra in VAGN['RA1'][uwise]]
dec_deg = [Angle(dec, unit='deg').deg for dec in VAGN['Dec1'][uwise]]
vagn_labels = ['bigMAC' for _ in ra_deg]
update_or_append_multiple(ra_deg, dec_deg, VAGN['z1'][uwise], vagn_labels)
print('Dual AGN sources in bigMAC catalog:', len(vagn_labels))


VAGN = pd.read_csv('data/Charisi16.txt', delimiter='|', comment=None, engine='python')
VAGN.columns = VAGN.columns.str.strip()
VAGN = VAGN.applymap(lambda x: x.strip() if isinstance(x, str) else x)
uwise = (VAGN['z']>0.001)&(VAGN['z']<zmax)
vagn_labels = ['Charisi16' for ra in VAGN['ra'][uwise]]
update_or_append_multiple(VAGN['ra'][uwise],VAGN['dec'][uwise],VAGN['z'][uwise],vagn_labels)
print('Charisi ',len(vagn_labels))

VAGN = pd.read_csv('data/Chen20.txt', delim_whitespace=True)
VAGN.columns = [col.strip() for col in VAGN.columns]
uwise = (VAGN['z']>0.001)&(VAGN['z']<zmax)
ra_deg = [Angle(ra, unit='hourangle').deg for ra in VAGN['rastr'][uwise]]
dec_deg = [Angle(dec, unit='deg').deg for dec in VAGN['decstr'][uwise]]
vagn_labels = ['Chen20' for ra in VAGN['rastr'][uwise]]
update_or_append_multiple(ra_deg, dec_deg, VAGN['z'][uwise], vagn_labels)
print('Chen20',len(vagn_labels))

VAGN = pd.read_csv('data/Liu19.txt', delim_whitespace=True)
uwise = (VAGN['z']>0.001)&(VAGN['z']<zmax)
vagn_labels = ['Liu19' for ra in VAGN['ra'][uwise]]
update_or_append_multiple(VAGN['ra'][uwise],VAGN['dec'][uwise],VAGN['z'][uwise],vagn_labels)
print('Liu19',len(vagn_labels))

VAGN = pd.read_csv('data/Ward22_wise.txt', delim_whitespace=True)
uwise = (VAGN['z']>0.001)&(VAGN['z']<zmax)
ra_deg = [Angle(ra, unit='hourangle').deg for ra in VAGN['rastr'][uwise]]
dec_deg = [Angle(dec, unit='deg').deg for dec in VAGN['decstr'][uwise]]
vagn_labels = ['Ward22_wise' for ra in VAGN['rastr'][uwise]]
update_or_append_multiple(ra_deg, dec_deg, VAGN['z'][uwise], vagn_labels)
print('Ward22_wise',len(vagn_labels))

VAGN = pd.read_csv('data/Ward22_ztf.txt', delim_whitespace=True)
uwise = (VAGN['z']>0.001)&(VAGN['z']<zmax)
ra_deg = [Angle(ra, unit='hourangle').deg for ra in VAGN['rastr'][uwise]]
dec_deg = [Angle(dec, unit='deg').deg for dec in VAGN['decstr'][uwise]]
vagn_labels = ['Ward22_ztf' for ra in VAGN['rastr'][uwise]]
update_or_append_multiple(ra_deg, dec_deg, VAGN['z'][uwise], vagn_labels)
print('Ward22_ztf',len(vagn_labels))

with fits.open('data/Graham15_PeriodicSample.fits') as hdul:
    vagn_data = Table(hdul[1].data)  # Assuming data is in the first extension
VAGN = vagn_data.to_pandas()

# Convert RA (hh:mm:ss) to degrees
VAGN["RA_deg"] = (VAGN["RAh"] + VAGN["RAm"]/60 + VAGN["RAs"]/3600) * 15  

# Convert Dec (dd:mm:ss) to degrees, applying the sign
VAGN["Dec_deg"] = VAGN["DEd"] + VAGN["DEm"]/60 + VAGN["DEs"]/3600
VAGN.loc[VAGN["DE-"] == "-", "Dec_deg"] *= -1  # Apply sign for negative Dec
uwise = (VAGN['z'] > 0.001) & (VAGN['z'] < zmax)
vagn_labels = ['Graham15' for _ in VAGN.loc[uwise, "RA_deg"]]
sky_coords = SkyCoord(VAGN.loc[uwise, "RA_deg"], VAGN.loc[uwise, "Dec_deg"], frame='icrs', unit='deg')
update_or_append_multiple(sky_coords.ra.deg, sky_coords.dec.deg, VAGN.loc[uwise, "z"], vagn_labels)
print('Graham15:',len(vagn_labels))

data = """
1 MBH-CSO0402+379 AT2018iih 61.455262 38.058954 0.05 Rodriguez et al., 2006 Binary
"""
lines = data.strip().split('\n')

# Loop through each line, extracting the required information
ras,decs,redshifts=[],[],[]
for line in lines:
    parts = line.split()
    id_ = parts[0]
    ra = float(parts[3].replace('−', '-'))
    dec = float(parts[4].replace('−', '-'))
    redshift = float(parts[5])
    if (redshift >0)&(redshift<zmax): 
        ras.append(ra)
        decs.append(dec)
        redshifts.append(redshift)

ras,decs,redshifts = np.array(ras),np.array(decs),np.array(redshifts)
TDE_labels = ['Rodriguez06' for ra in ras]
update_or_append_multiple(ras, decs,redshifts,TDE_labels)
```

# Save dataframe 

```python
# Assuming `df` is your pandas DataFrame
df.to_csv('data/AGNsample_dual_big.csv', index=False)

```

# Change format to ecsv and bitwise lable

```python
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.table import Table
import astropy.units as u
import numpy as np

# Read the CSV file into a DataFrame
df = pd.read_csv('data/AGNsample_dual_big.csv')

# Define a function to parse the strings into SkyCoord objects
def parse_skycoord_string(s):
    # Split the string to extract RA and Dec values
    parts = s.split(',')
    # Extract RA and Dec parts, removing unwanted characters
    ra_str = parts[1].strip().split('(')[-1]
    dec_str = parts[2].strip().split(')')[0]
    # Convert to float and create a SkyCoord object
    return SkyCoord(ra=float(ra_str)*u.deg, dec=float(dec_str)*u.deg)

# Apply this function to each row in the 'SkyCoord' column
df['SkyCoord_obj'] = df['SkyCoord'].apply(parse_skycoord_string)

# Now, 'SkyCoord_obj' column contains SkyCoord objects
# You can access RA and Dec directly from these objects if needed
# For example, to add RA and Dec as separate columns in degrees:
df['coord.ra'] = df['SkyCoord_obj'].apply(lambda x: x.ra.degree)
df['coord.dec'] = df['SkyCoord_obj'].apply(lambda x: x.dec.degree)

df.drop('SkyCoord_obj', axis=1, inplace=True)

```

```python
# Initialize agnlabels
agnlabels=['SDSS_QSO','WISE_Variable','Optical_Variable','Charisi16','Chen20','Graham15','Liu19','Ward22_wise','Ward22_ztf','bigMAC','Rodriguez06']

# Calculate the sum of label columns for each row
# Calculate the bitwise sum
bitwise_sum = np.zeros(len(df), dtype=int)
for i, label in enumerate(agnlabels):
    bitwise_sum += df[label].values * (2 ** i)

df['label'] = bitwise_sum
df['objectid'] = df.index

```

```python
import pandas as pd
from astropy.table import Table

# Select necessary columns and create a 'coord' column
df['coord'] = df['coord.ra'].astype(str) + ', ' + df['coord.dec'].astype(str)
selected_columns_df = df[['objectid', 'coord', 'label']]

# Convert the modified DataFrame to an astropy Table
t = Table.from_pandas(selected_columns_df)

# Write the table to an ECSV file
t.write('data/sample_big.ecsv', format='ascii.ecsv', overwrite=True)

```

```python
selected_columns_df = df[['objectid', 'coord.ra', 'coord.dec', 'label']]
t = Table.from_pandas(selected_columns_df)
t.write('data/BOSS_Mar5.ecsv', format='ascii.ecsv', overwrite=True)

```

```python
def translate_bitwise_sum_to_labels(bitwise_sum):
    """
    Translate a bitwise sum back to the labels which were set to 1.

    Parameters:
    - bitwise_sum: Integer, the bitwise sum representing the combination of labels.
    - labels: List of strings, the labels corresponding to each bit position.

    Returns:
    - List of strings, the labels that are set to 1.
    """
    # Initialize agnlabels
    agnlabels=['SDSS_QSO','WISE_Variable','Optical_Variable','Charisi16','Chen20','Graham15','Liu19','Ward22_wise','Ward22_ztf','bigMAC','Rodriguez06']

    active_labels = []
    for i, label in enumerate(agnlabels):
        # Check if the ith bit is set to 1
        if bitwise_sum & (1 << i):
            active_labels.append(label)
    return active_labels

# Example usage
bitwise_sum_example = 5  # For example, if the binary representation is '101'

# Translate the bitwise sum back to active labels
active_labels = translate_bitwise_sum_to_labels(bitwise_sum_example)
print("Active Labels:", active_labels)

```

```python

```

```python
# https://ui.adsabs.harvard.edu/abs/2003MNRAS.346.1055K/abstract
# https://wwwmpa.mpa-garching.mpg.de/SDSS/DR4/Data/agn.cat_dr2_desc.txt
# idpl,mjd,ifib, ra, dec, zz, o3lum,o3corr, bpt1,bpt2, rml50, rmu, con, d4n,hda, vdisp  

#r = 0
#new_rows = []  # Prepare a list to collect new rows
#with open("data/agn.dat_dr4_release.v2", 'r') as file:
#    for line in file:
#        parts = line.split()  # Splits the line into parts
#        redshift = float(parts[5])
#        transition = parts[2]
#        ra = parts[3]
#        dec = parts[4]
#        coord = SkyCoord(ra + ' ' + dec, unit=(u.deg, u.deg))
#        new_row = pd.DataFrame([{**{'SkyCoord': coord, 'redshift': redshift}, **{l: 1 if l == label else 0 for l in agnlabels}}])
#        new_rows.append(new_row)        
#        r+=1

#df = pd.concat([df, *new_rows], ignore_index=True)

#print('Kauffmann sources:',str(r))
```

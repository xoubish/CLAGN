---
jupyter:
  jupytext:
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.16.1
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
agnlabels = ['SDSS_QSO', 'WISE_Variable','Optical_Variable','Galex_Variable',
             'Turn-on', 'Turn-off',
             'SPIDER', 'SPIDER_AGN','SPIDER_BL','SPIDER_QSOBL','SPIDER_AGNBL', 
             'TDE','BOSS SF']

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
num = 1000  # Total number of samples desired
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

# SPIDERS 
SPectroscopic IDentfication of ERosita Sources- (value added SDSS DR16 which are all detected in xray) https://www.sdss.org/dr18/bhm/programs/spiders/

```python
a = fits.getdata('data/VAC_SPIDERS_2RXS_DR16.fits',1)
#print(a.columns)
print(np.unique(a['DR16_CLASS']))
#print(np.unique(a['DR16_SUBCLASS']))
print('CLASS GALAXY DR16 SPIDER:',len(a['DR16_SUBCLASS'][a['DR16_CLASS']=='GALAXY']),' subclasses: ',np.unique(a['DR16_SUBCLASS'][a['DR16_CLASS']=='GALAXY']))
print('CLASS QSO DR16 SPIDER:',len(a['DR16_SUBCLASS'][a['DR16_CLASS']=='QSO']),' subclasses: ',np.unique(a['DR16_SUBCLASS'][a['DR16_CLASS']=='QSO']))
# all spider that can be AGN/QSO
uspider = (a['DR16_CLASS']!='STAR')&(a['DR16_SUBCLASS']!='STARFORMING')&(a['DR16_SUBCLASS']!='STARFORMING BROADLINE')&(a['DR16_SUBCLASS']!='STARBURST')&(a['DR16_SUBCLASS']!='STARBURST BROADLINE')&(a['Z_BEST']>0)&(a['Z_BEST']<zmax)

# all GALAXY and QSO AGN Broadlines
uspiderbroad = (a['DR16_CLASS']!='STAR')&(a['DR16_SUBCLASS']=='AGN BROADLINE')&(a['Z_BEST']<=zmax)&(a['Z_BEST']>=0.0)
uspideragn = (a['DR16_CLASS']!='STAR')&(a['DR16_SUBCLASS']=='AGN')&(a['Z_BEST']<=zmax)&(a['Z_BEST']>=0.0)

# Gal AGNs
ugal_agn = (a['DR16_CLASS']=='GALAXY') &(a['DR16_SUBCLASS']=='AGN')&(a['Z_BEST']<=zmax)&(a['Z_BEST']>=0.0)
ugal_agnbl = (a['DR16_CLASS']=='GALAXY') &(a['DR16_SUBCLASS']=='AGN BROADLINE')&(a['Z_BEST']<=zmax)&(a['Z_BEST']>=0.0)
uqso_agn = (a['DR16_CLASS']=='QSO') &(a['DR16_SUBCLASS']=='AGN')&(a['Z_BEST']<=zmax)&(a['Z_BEST']>=0.0)
uqso_agnbl = (a['DR16_CLASS']=='QSO') &(a['DR16_SUBCLASS']=='AGN BROADLINE')&(a['Z_BEST']<=zmax)&(a['Z_BEST']>=0.0)
#print(len(a[uspider]),len(a[uspiderbroad]),len(a[uspideragn]),len(a[ugal_agn]),len(a[ugal_agnbl]),len(a[uqso_agn]),len(a[uqso_agnbl]))

```

```python
#takes long
#uspider_labels = ['SPIDER' for ra in a['SDSS_RA'][uspider]]
#update_or_append_multiple(a['SDSS_RA'][uspider],a['SDSS_DEC'][uspider],a['Z_BEST'][uspider],uspider_labels)
#print('SPIDER QSO/AGNs added :',str(len(a['SDSS_RA'][uspider])))

uspideragn_labels = ['SPIDER_AGN' for ra in a['SDSS_RA'][uspideragn]]
update_or_append_multiple(a['SDSS_RA'][uspideragn],a['SDSS_DEC'][uspideragn],a['Z_BEST'][uspideragn],uspideragn_labels)
print('SPIDER spec AGN (no BL):',str(len(a['SDSS_RA'][uspideragn])))

uspiderbl_labels = ['SPIDER_BL' for ra in a['SDSS_RA'][uspiderbroad]]
update_or_append_multiple(a['SDSS_RA'][uspiderbroad],a['SDSS_DEC'][uspiderbroad],a['Z_BEST'][uspiderbroad],uspiderbl_labels)
print('SPIDER spec BL:',str(len(a['SDSS_RA'][uspiderbroad])))

uspideragn_labels = ['SPIDER_AGNBL' for ra in a['SDSS_RA'][ugal_agnbl]]
update_or_append_multiple(a['SDSS_RA'][ugal_agnbl],a['SDSS_DEC'][ugal_agnbl],a['Z_BEST'][ugal_agnbl],uspideragn_labels)
print('SPIDER spec GAL BL:',str(len(a['SDSS_RA'][ugal_agnbl])))

uspideragn_labels = ['SPIDER_QSOBL' for ra in a['SDSS_RA'][uqso_agnbl]]
update_or_append_multiple(a['SDSS_RA'][uqso_agnbl],a['SDSS_DEC'][uqso_agnbl],a['Z_BEST'][uqso_agnbl],uspideragn_labels)
print('SPIDER spec QSO BL:',str(len(a['SDSS_RA'][uqso_agnbl])))

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

```python
paper = Ned.query_refcode('2022ApJ...933...37W') #Galex Variable
up = (paper['Redshift']>0)&(paper['Redshift']<=zmax)
paper_labels = ['Galex_Variable' for ra in paper['RA'][up]]
update_or_append_multiple(paper['RA'][up], paper['DEC'][up],paper['Redshift'][up],paper_labels)
print('Galex Variable sources added: ',len(paper_labels))
```

# Changing Looks (dividing to on and off)

```python
# one LaMassa et al. 2015 one CL at 0.31 which turned off
paper = Ned.query_refcode('2015ApJ...800..144L') 
update_or_append_multiple([paper[0]['RA']],[paper[0]['DEC']],[paper[0]['Redshift']],['Turn-off'])
print('LaMassa 2015 added sources: 1')

#------------------------------------------------------------------------------
# McLeod et al. 2016 4/6 CL that turned on/off
paper = Table.read('https://academic.oup.com/mnras/article/457/1/389/989199', htmldict={'table_id': 5}, format='ascii.html')
for i in range(len(paper)):
    redshift = paper['z\n            .'][i]
    if (redshift>0)&(redshift<zmax):
        if (paper['Max(Δg)\n            .'][i][0] =='−'):
            coord_str = paper['Name\n            .'][i]
            test_str = coord_str[0:2]+ " "+ coord_str[2:4]+ " " + coord_str[4:9] + " " + coord_str[9:12] + " " + coord_str[12:14]+ " " + coord_str[14:]
            c = SkyCoord(test_str, unit=(u.hourangle, u.deg))
            update_or_append_multiple([c.ra.deg],[c.dec.deg],[redshift],['Turn-on'])
        else:
            coord_str = paper['Name\n            .'][i]
            test_str = coord_str[0:2]+ " "+ coord_str[2:4]+ " " + coord_str[4:9] + " " + coord_str[9:12] + " " + coord_str[12:14]+ " " + coord_str[14:]
            c = SkyCoord(test_str, unit=(u.hourangle, u.deg))
            update_or_append_multiple([c.ra.deg],[c.dec.deg],[redshift],['Turn-off'])
print('McLeod 2016 added sources:',str(len(paper)))

#------------------------------------------------------------------------------
# Ruan et al. 2016 three CL which turned off
paper = Ned.query_refcode('2016ApJ...826..188R') 
update_or_append_multiple([paper[0]['RA']],[paper[0]['DEC']],[paper[0]['Redshift']],['Turn-off'])
print('Ruan 2016 added sources: 1')

#------------------------------------------------------------------------------
#Yang et al. 2018, 27 (21 new) CL-AGNs z<0.58
with open('data/Yang2018_table5.csv', mode='r', newline='') as file:
    reader = csv.DictReader(file)
    for row in reader:
        # Convert RA and DEC to degrees using SkyCoord
        coord = SkyCoord(row['RA'], row['DEC'], unit=(u.hourangle, u.deg))
        ra_deg = coord.ra.deg
        dec_deg = coord.dec.deg  
        update_or_append_multiple([ra_deg], [dec_deg],[float(row['Redshift'])],[row['Transition']])
print('Yang 2018 added sources:',str(reader.line_num))

#------------------------------------------------------------------------------
#MacLeod 2019 has turn-on/off both spectroscopically confirmed (1/16) and candidates (17/245) to be confirmed
Vizier.ROW_LIMIT = -1
catalog_list = Vizier.find_catalogs('2019ApJ...874....8M')
catalogs = Vizier.get_catalogs(catalog_list.keys())
table2 = catalogs[0]  #more than one table
g1,g2 = table2['_tab1_6'],table2['_tab1_9']
dg = g1-g2
uon_sp = (table2['z']<zmax)&(table2['z']>0)&(dg>0)&(table2['CLQ_'] > 0) & (table2['Nsigma'] > 3)
uoff_sp = (table2['z']<zmax)&(table2['z']>0)&(dg<0)&(table2['CLQ_'] > 0) & (table2['Nsigma'] > 3)
uon_all = (table2['z']<zmax)&(table2['z']>0)&(dg>0) #candidates not spec
uoff_all = (table2['z']<zmax)&(table2['z']>0)&(dg<0) #candidates not spec
cllabels = ['Turn-on' for ra in table2['_RA'][uon_sp]]
update_or_append_multiple(table2['_RA'][uon_sp], table2['_DE'][uon_sp],table2['z'][uon_sp],cllabels)
cllabels = ['Turn-off' for ra in table2['_RA'][uoff_sp]]
update_or_append_multiple(table2['_RA'][uoff_sp], table2['_DE'][uoff_sp],table2['z'][uoff_sp],cllabels)
print('McLeod 2019 added sources:',str(len(table2['_RA'][uon_sp])+len(table2['_RA'][uoff_sp])))

#------------------------------------------------------------------------------
#Graham et al. 2020 Brighter than CL AGNs, they find 111 Changing State Quasars, 
# 48 declining and 63 increasing H\beta flux.
CSQ = Table.read('https://academic.oup.com/mnras/article/491/4/4925/5634279', htmldict={'table_id': 5}, format='ascii.html')
#get coords from "name" column for this
for i in range(len(CSQ)):
    coord_str = CSQ['Name\n            .'][i]
    test_str = coord_str[6:8]+ " "+ coord_str[8:10]+ " " + coord_str[10:14] + " " + coord_str[14:17] + " " + coord_str[17:19]+ " " + coord_str[19:]
    c = SkyCoord(test_str, unit=(u.hourangle, u.deg))

    redshift = float(CSQ[i][2])
    typee = CSQ[i][9]
    if typee[0]=='−':
        update_or_append_multiple([c.ra.deg],[c.dec.deg],[redshift],['Turn-off'])
    else:
        update_or_append_multiple([c.ra.deg],[c.dec.deg],[redshift],['Turn-on'])
print('Graham 2019 added sources:',str(len(CSQ)))

#------------------------------------------------------------------------------
#Sheng et al. 2020 confirmed 6 Turn-off AGNs (one was not robust) (IR variable to no variability) 2020ApJ...889...46S
paper = Ned.query_refcode('2020ApJ...889...46S')
sheng_CLQ = [0,2,3,4,5,6]
cllabels = ['Turn-off' for ra in paper[sheng_CLQ]['RA']]
update_or_append_multiple(paper[sheng_CLQ]['RA'],paper[sheng_CLQ]['DEC'],paper[sheng_CLQ]['Redshift'],cllabels)
print('Green 2020 added sources:',str(len(cllabels)))

#------------------------------------------------------------------------------
# Green et al. 2022, Identified 19 changing look quasars from their variation 
#in the $H\beta$ line in multi-epoch SDSS-IV spectra, four of which with significant
#increase in the broad line width and the rest with dimming or disappearing.        
Vizier.ROW_LIMIT = -1
catalog_list = Vizier.find_catalogs('J/ApJ/933/180')
catalogs = Vizier.get_catalogs(catalog_list.keys())
table2 = catalogs[0]
#go to pandas to manipulate the table
green_CSQ = table2.to_pandas()
green_CSQ = green_CSQ[green_CSQ['Notes'].str.contains("CLQ", na = False)]
#pick out the coordinates from the 'SDSS' column
coord_str = green_CSQ['SDSS'].astype('string')
test_str = coord_str.str[1:3]+ " "+ coord_str.str[3:5]+ " " + coord_str.str[5:10] + " " + coord_str.str[10:13] + " " + coord_str.str[13:15]+ " " + coord_str.str[15:]
c = SkyCoord(test_str.values.tolist() , unit=(u.hourangle, u.deg))
redshifts = np.array(green_CSQ['zspec'].astype('float'))
notes = green_CSQ['Notes']
labels = ['Turn-on' if 'TurnOn' in note else 'Turn-off' for note in notes]
update_or_append_multiple(c.ra.deg,c.dec.deg,redshifts,labels)
print('Green 2022 added sources:',str(len(labels)))

#------------------------------------------------------------------------------
#get_lyu_sample(coords, labels)  #2022ApJ...927..227L on/off not indicated

#------------------------------------------------------------------------------
#Lopez et al. 2022 confirmed 4 turn-on AGN using random forest
result_table = Simbad.query_bibobj('2022MNRAS.513L..57L')
result_table = result_table[[0,1,2,3]]  #pick the correct sources by hand
redshifts = [0.1022,0.2375,0.1083,0.0849]# by hand from paper
ln_coords = [SkyCoord(ra, dec, frame='icrs', unit='deg') for ra, dec in zip(result_table['RA'], result_table['DEC'])]
for i,r in enumerate(redshifts):
    update_or_append_multiple([ln_coords[i].ra.deg],[ln_coords[i].dec.deg],[r],['Turn-on'])

print('Lopez 2022 added sources:',str(len(redshifts)))

#------------------------------------------------------------------------------
# Hon et al. 2022 searched repeated spectra for appearence/disappearence of BELs. 
# Hence sample is not biased towards brighter QSOs. Difference between TDEs and 
# turn-on BEL is timescale, so turn-on CLAGNs should be followed up to make sure they are not TDEs.
# should be 25 turn on and 2 turn off in this
r = 0
with open("data/hon22.cat", 'r') as file:
    for line in file:
        parts = line.split()  # Splits the line into parts
        coord_str = parts[0]
        redshift = float(parts[1])
        transition = parts[2]
        ra = coord_str[1:3] + 'h' + coord_str[3:5] + 'm' + coord_str[5:7] + 's'
        dec = coord_str[8:11] + 'd' + coord_str[11:13] + 'm' + coord_str[13:] + 's'
        coord = SkyCoord(ra + ' ' + dec, unit=(u.hourangle, u.deg))
        update_or_append_multiple([coord.ra.deg],[coord.dec.deg],[redshift],[transition])
        r+=1
print('Hon 2022 added sources:',str(r))

```

# TDEs

```python
data = """
1 ZTF18acaqdaa AT2018iih 262.0163662 30.6920758 0.212 van Velzen et al. (2021) TDE-He
2 ZTF18acnbpmd AT2018jbv 197.6898587 8.5678292 0.340 Hammerstein et al. (2023) TDE-featureless
3 ZTF19aabbnzo AT2018lna 105.8276892 23.0290953 0.0914 van Velzen et al. (2021) TDE-H+He
4 ZTF19aaciohh AT2019baf 268.0005082 65.6266546 0.0890 This paper; J. Somalwar et al. (2023, in preparation) Unknown
5 ZTF17aaazdba AT2019azh 123.3206388 22.6483180 0.0222 Hinkle et al. (2021) TDE-H+He
6 ZTF19aakswrb AT2019bhf 227.3165243 16.2395720 0.121 van Velzen et al. (2021) TDE-H
7 ZTF19aaniqrr AT2019cmw 282.1644974 51.0135422 0.519 This paper; J. Wise et al. (2023, in preparation) TDE-featureless
8 ZTF19aapreis AT2019dsg 314.2623552 14.2044787 0.0512 Stein et al. (2021) TDE-H+He
9 ZTF19aarioci AT2019ehz 212.4245268 55.4911223 0.0740 van Velzen et al. (2021) TDE-H
10 ZTF19abzrhgq AT2019qiz 71.6578313 −10.2263602 0.0151 Nicholl et al. (2020) TDE-H+He
11 ZTF19acspeuw AT2019vcb 189.7348778 33.1658869 0.0890 Hammerstein et al. (2023) TDE-H+He
12 ZTF20aabqihu AT2020pj 232.8956925 33.0948917 0.0680 Hammerstein et al. (2023) TDE-H+He
13 ZTF20abfcszi AT2020mot 7.8063109 85.0088329 0.0690 Hammerstein et al. (2023) TDE-H+He
14 ZTF20abgwfek AT2020neh 230.3336852 14.0696032 0.0620 Angus et al. (2022) TDE-H+He
15 ZTF20abnorit AT2020ysg 171.3584535 27.4406021 0.277 Hammerstein et al. (2023) TDE-featureless
16 ZTF20acaazkt AT2020vdq 152.2227354 42.7167535 0.0450 This paper; J. Somalwar et al. (2023, in preparation) Unknown
17 ZTF20achpcvt AT2020vwl 232.6575481 26.9824432 0.0325 Hammerstein et al. (2021a) TDE-H+He
18 ZTF20acitpfz AT2020wey 136.3578499 61.8025699 0.0274 Arcavi et al. (2020) TDE-H+He
19 ZTF20acnznms AT2020yue 165.0013942 21.1127532 0.204 This paper TDE-H?
20 ZTF20acvezvs AT2020abri 202.3219785 19.6710235 0.178 This paper Unknown
21 ZTF20acwytxn AT2020acka 238.7581288 16.3045292 0.338 Hammerstein et al. (2021b) TDE-featureless
22 ZTF21aaaokyp AT2021axu 176.6514953 30.0854257 0.192 Hammerstein et al. (2021c) TDE-H+He
23 ZTF21aakfqwq AT2021crk 176.2789219 18.5403839 0.155 This paper TDE-H+He?
24 ZTF21aanxhjv AT2021ehb 46.9492531 40.3113468 0.0180 Yao et al. (2022a) TDE-featureless
25 ZTF21aauuybx AT2021jjm 219.8777384 −27.8584845 0.153 Yao et al. (2021d) TDE-H
26 ZTF21abaxaqq AT2021mhg 4.9287185 29.3168745 0.0730 Chu et al. (2021b) TDE-H+He
27 ZTF21abcgnqn AT2021nwa 238.4636684 55.5887978 0.0470 Yao et al. (2021b) TDE-H+He
28 ZTF21abhrchb AT2021qth 302.9121723 −21.1602187 0.0805 This paper TDE-coronal
29 ZTF21abjrysr AT2021sdu 17.8496154 50.5749060 0.0590 Chu et al. (2021c) TDE-H+He
30 ZTF21abqhkjd AT2021uqv 8.1661654 22.5489257 0.106 Yao (2021) TDE-H+He
31 ZTF21abqtckk AT2021utq 229.6212498 73.3587323 0.127 This paper TDE-H
32 ZTF21abxngcz AT2021yzv 105.2774821 40.8251799 0.286 Chu et al. (2022) TDE-featureless
33 ZTF21acafvhf AT2021yte 103.7697396 12.6341503 0.0530 Yao et al. (2021a) TDE-H+He
"""

# Split the data into lines
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
TDE_labels = ['TDE' for ra in ras]

update_or_append_multiple(ras, decs,redshifts,TDE_labels)
print('TDEs added to df: ',len(TDE_labels))
```

# Fermi (gamma ray) Blazars

```python
r = 0
with open("data/fermi.cat", 'r') as file:
    for line in file:
        parts = line.split()  # Splits the line into parts
        if (float(parts[12])>0.1) and (float(parts[12])<zmax):
            ra = parts[4]
            dec = parts[5]
            coord = SkyCoord(ra + ' ' + dec, unit=(u.hourangle, u.deg))
            update_or_append_multiple([coord.ra.deg],[coord.dec.deg],[float(parts[12])],['Fermi_Blazars'])    
            r+=1
print('Fermi Blazar added sources:',str(r))
```

# WISE R90

```python
# too big of a sample million entries and no redshift info so will not worry for now.
#r90 = fits.getdata('data/J_ApJS_234_23_r90cat.dat.gz.fits.gz',1)
#r90.columns
```

# BOSS SF in BPT low z

```python
#large! and not in line with previous, its stand alone
a = fits.getdata('data/portsmouth_emlinekin_full-DR12-boss.fits')
u = (a['REDSHIFT']>=0.00015) &(a['REDSHIFT']<=0.39457)&(a['BPT']=='Star Forming')
ra,dec,zz = a['RA'][u],a['DEC'][u],a['REDSHIFT'][u]
sflab = np.full(ra.shape, 'BOSS SF')
update_or_append_multiple(ra,dec,zz, sflab)
# Assuming `df` is your pandas DataFrame
#df.to_csv('data/BOSS-SF.csv', index=False)

```

# Save dataframe 

```python
# Assuming `df` is your pandas DataFrame
df.to_csv('data/AGNsample_26Feb24.csv', index=False)

```

# Change format to ecsv and bitwise lable

```python
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.table import Table
import astropy.units as u
import numpy as np

# Read the CSV file into a DataFrame
df = pd.read_csv('data/sample.csv')

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
agnlabels = ['SDSS_QSO', 'WISE_Variable','Optical_Variable','Galex_Variable',
             'Turn-on', 'Turn-off',
             'SPIDER', 'SPIDER_AGN','SPIDER_BL','SPIDER_QSOBL','SPIDER_AGNBL', 
             'TDE','BOSS SF']
# Calculate the sum of label columns for each row
# Calculate the bitwise sum
bitwise_sum = np.zeros(len(df), dtype=int)
for i, label in enumerate(agnlabels):
    bitwise_sum += df[label].values * (2 ** i)

df['label'] = bitwise_sum
df['objectid'] = df.index

```

```python
#sampled_df = df.sample(n=10000, random_state=42)
#selected_columns_df = sampled_df[['objectid', 'coord.ra', 'coord.dec', 'label']]

t = Table.from_pandas(df)
t.write('data/sample.ecsv', format='ascii.ecsv', overwrite=True)
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
    agnlabels = ['SDSS_QSO', 'WISE_Variable','Optical_Variable','Galex_Variable',
                 'Turn-on', 'Turn-off',
                 'SPIDER', 'SPIDER_AGN','SPIDER_BL','SPIDER_QSOBL','SPIDER_AGNBL', 
                 'TDE']
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

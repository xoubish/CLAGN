#!/bin/zsh
# Driver for the steps that follow the pool unWISE fetch (02_parent_pool.py wise):
#   project -> pick the subset worth enriching -> archival spectra + current ZTF for that subset
#   -> score/tier/allocate -> finder charts.
# Run from the CLAGN folder:  ./run_after_wise.sh
set -e
PY=/opt/anaconda3/bin/python
cd "$(dirname "$0")"

echo "== 1. project pool onto the W1 manifold"
$PY 02_parent_pool.py project 2>&1 | grep -v -i "warning\|warn(\|it/s\]"

echo "== 2. subset for enrichment (T1 candidates + observable control pool)"
$PY - <<'EOF'
import pandas as pd, numpy as np
ps = pd.read_csv('data/parent_pool_scored.csv'); po = pd.read_csv('data/parent_pool_obs.csv')
d = ps.merge(po[['poolid','hrs_any','moonsep_min']], on='poolid')
d = d[d.projected.fillna(False)]
M = np.clip(np.fmax(d.zeltyn_density_ratio/3.0, d.clagn_score/0.15), 0, 2)
obs = (d.hrs_any >= 1.5) & (d.moonsep_min >= 30)          # no magnitude cut; brightness weights the ranking
d['M'] = M
B = np.select([d.psfmag_r < 18.5, d.psfmag_r < 19.0, d.psfmag_r < 19.5], [1.0, 0.8, 0.5], 0.25)
d['preprio'] = B * M
t1 = d[(M >= 1) & obs].sort_values('preprio', ascending=False).head(700)   # cap enrichment cost (~2 s + ~5 s per object)
ctrl = d[(d.clagn_score == 0) & (d.zeltyn_density_ratio <= 0.3) & ~d.in_region_zeltyn & ~d.in_region_clagn & obs & (d.psfmag_r <= 19.0)]
ctrl = ctrl.sample(n=min(150, len(ctrl)), random_state=1)
sub = pd.concat([t1, ctrl]).drop_duplicates('poolid')
sub['name'] = 'P' + sub.poolid.astype(str)
sub[['name','poolid','ra','dec','z','psfmag_r','mjd','clagn_score','zeltyn_density_ratio','in_region_zeltyn','in_region_clagn']].to_csv('data/pool_subset_for_enrich.csv', index=False)
print(f'T1 candidates (M>=1, observable, r<=19.5): {len(t1)} | control pool sampled: {len(ctrl)} | subset total: {len(sub)}')
EOF

echo "== 3. archival spectra for the subset (SDSS DR19 allspec + DESI DR1)"
$PY 03_spectra_inventory.py data/pool_subset_for_enrich.csv pool ra dec name 2>&1 | grep -v -i "warning\|warn(" | tail -6

echo "== 4. current ZTF photometry for the subset (reference epoch = archival SDSS spectrum MJD)"
$PY 03b_ztf_now.py data/pool_subset_for_enrich.csv pool ra dec name mjd 2>&1 | grep -v -i "warning\|warn(" | tail -3

echo "== 4b. NEOWISE-R single-exposure W1 to 2024 (IRSA Gator bulk upload) for the subset and the Zeltyn sample"
$PY 03c_neowise_now.py data/pool_subset_for_enrich.csv pool ra dec name 2>&1 | grep -v -i "warning\|warn(" | tail -1
$PY 03c_neowise_now.py data/zeltyn_coords.csv zeltyn ra dec name 2>&1 | grep -v -i "warning\|warn(" | tail -1

echo "== 5. score, tier, allocate"
$PY 04_score_tiers.py 2>&1 | grep -v -i "warning\|warn(" | tee data/score_tiers.log | head -20

echo "== 6. finder charts"
for n in sep23 oct26 oct27; do rm -rf finders/$n; $PY 06_finder_charts.py data/targets_$n.csv 2>&1 | grep -v -i "warning\|warn(" | tail -1; done

echo "== 6b. image cutouts (SDSS gri JPEG, ZTF g/r reference images via IRSA IBE)"
$PY 07b_cutouts.py 2>&1 | grep -v -i "warning\|warn(" | tail -1

echo "== 6c. archival spectra (SDSS all epochs via sas_url, DESI DR1 via SPARCL) for the listed targets"
$PY 03d_fetch_spectra.py 2>&1 | grep -v -i "warning\|warn(" | tail -1

echo "== 7. night sheet web page (web/clagn_night_sheet.html; republish with the Artifact tool / share the file)"
$PY 07_make_webpage.py 2>&1 | grep -v -i "warning\|warn(" | tail -1
echo "done."

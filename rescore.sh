#!/bin/zsh
# Re-run only the fast tail of the pipeline after new enrichment data land (NEOWISE, SDSS-V internal epochs, ...):
#   04 score/tier/allocate -> 06 finder charts -> 07b cutouts -> 07 night-sheet page.
# Everything upstream (projection, spectra inventory, ZTF pull) is reused from data/.  Then republish the page.
set -e
PY=/opt/anaconda3/bin/python
cd "$(dirname "$0")"
echo "== score, tier, allocate";  $PY 04_score_tiers.py 2>&1 | grep -v -i "warning\|warn(" | tee data/score_tiers.log | sed -n 1,8p
echo "== finder charts";          for n in sep23 oct26 oct27; do rm -rf finders/$n; $PY 06_finder_charts.py data/targets_$n.csv 2>&1 | grep -v -i "warning\|warn(" | tail -1; done
echo "== cutouts";                $PY 07b_cutouts.py 2>&1 | grep -v -i "warning\|warn(" | tail -1
echo "== night sheet";            $PY 07_make_webpage.py 2>&1 | grep -v -i "warning\|warn(" | tail -1
echo "done."

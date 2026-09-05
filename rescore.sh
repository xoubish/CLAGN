#!/bin/zsh
# Re-run only the fast tail of the pipeline after new enrichment data land (NEOWISE, SDSS-V internal epochs, ...):
#   04 score/tier/allocate -> 11 observing sequence (time-aware fill) -> 06 finder charts -> 06b NGPS lists -> 07b cutouts -> 03d spectra -> 07 page.
# Everything upstream (projection, spectra inventory, ZTF pull) is reused from data/.  Then republish the page.
set -e
PY=/opt/anaconda3/bin/python
cd "$(dirname "$0")"
echo "== score, tier, allocate";  $PY 04_score_tiers.py 2>&1 | grep -v -i "warning\|warn(" | tee data/score_tiers.log | sed -n 1,8p
echo "== observing sequence";     $PY 11_schedule.py 2>&1 | grep -v -i "warning\|warn(" | tee data/schedule.log | grep -E "^(sep|oct)"
echo "== finder charts";          for n in sep23 oct26 oct27; do $PY 06_finder_charts.py data/targets_$n.csv 2>&1 | grep -v -i "warning\|warn(" | tail -1; done
echo "== NGPS target lists";     $PY 06b_ngps_targetlist.py 2>&1 | tail -3     # after the finders step, which recreates finders/<night>/
echo "== cutouts";                $PY 07b_cutouts.py 2>&1 | grep -v -i "warning\|warn(" | tail -1
echo "== archival spectra";       $PY 03d_fetch_spectra.py 2>&1 | grep -v -i "warning\|warn(" | tail -1
echo "== night sheet";            $PY 07_make_webpage.py 2>&1 | grep -v -i "warning\|warn(" | tail -1
echo "done."

#!/bin/zsh
# Selection v2 chain (2026-09-11): galaxy arm -> selection -> time-aware schedule -> NGPS lists -> epochs -> spectra -> cutouts -> finders -> page.
set -e; PY=/opt/anaconda3/bin/python; cd "$(dirname "$0")"; export CLAGN_SEL=v2
echo "== 14b turn-on arm";      $PY 14b_turnon_arm.py 2>&1 | grep -v -i "warning" | head -4
echo "== 14 selection";         $PY 14_select_targets.py 2>&1 | grep -v -i "warning" | head -5
echo "== 11 schedule";          $PY 11_schedule.py 2>&1 | grep -E "^(sep|oct)" | cut -c1-140
echo "== 06b NGPS lists";       $PY 06b_ngps_targetlist.py 2>&1 | grep "rows ->"
echo "== 14c epochs";           $PY 14c_v2_epochs.py 2>&1 | tail -1
echo "== 03d spectra";          $PY 03d_fetch_spectra.py data/targets_sep23_v2.csv data/targets_oct26_v2.csv data/targets_oct27_v2.csv 2>&1 | grep -v -i "warning" | tail -1
echo "== 07b cutouts";          $PY 07b_cutouts.py data/targets_sep23_v2.csv data/targets_oct26_v2.csv data/targets_oct27_v2.csv 2>&1 | grep -v -i "warning" | tail -1
echo "== 06 finder charts";     for n in sep23 oct26 oct27; do $PY 06_finder_charts.py data/targets_${n}_v2.csv 2>&1 | grep -v -i "warning" | tail -1; done
echo "== 07 night sheet";       $PY 07_make_webpage.py 2>&1 | grep "wrote"
echo "done."

echo "== publish (official files = v2)"; ./publish_v2.sh

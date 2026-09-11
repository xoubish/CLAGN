#!/bin/zsh
# Make the v2 selection the official one (user decision 2026-09-11): copy the _v2 lane onto the file names the observing
# scripts, the artifact and GitHub Pages use.  Run after rebuild_v2.sh.
set -e; cd "$(dirname "$0")"
for n in sep23 oct26 oct27; do
  cp data/targets_${n}_v2.csv data/targets_${n}.csv
  cp data/schedule_${n}_v2.csv data/schedule_${n}.csv
  mkdir -p finders/$n; rm -f finders/$n/*.png finders/$n/*.txt finders/$n/ngps_*.csv
  cp finders/${n}_v2/*.png finders/$n/ 2>/dev/null || true
  cp finders/${n}_v2/ngps_${n}_fixed.csv finders/${n}_v2/ngps_${n}_snr.csv finders/$n/
  cp finders/${n}_v2/targetlist_${n}_v2.txt finders/$n/targetlist_${n}.txt
  cp finders/${n}_v2/schedule_${n}.txt finders/$n/schedule_${n}.txt 2>/dev/null || true
done
cp web/clagn_night_sheet_v2.html web/clagn_night_sheet.html
cp docs/index_v2.html docs/index.html
echo "official files now = selection v2: $(for n in sep23 oct26 oct27; do echo -n "$n $(awk -F, 'NR>1 && $1>0' data/targets_${n}.csv | wc -l | tr -d ' ') primaries; "; done)"

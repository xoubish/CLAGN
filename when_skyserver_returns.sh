#!/bin/zsh
# Wait for SkyServer DR19 to answer, then complete the public multi-epoch SDSS inventory for the official targets and rebuild the page.
cd "$(dirname "$0")"; PY=/opt/anaconda3/bin/python
until [ "$(curl -s --max-time 20 -o /dev/null -w '%{http_code}' 'https://skyserver.sdss.org/dr19/SkyServerWS/SearchTools/SqlSearch?cmd=SELECT%20TOP%201%20specobjid%20FROM%20allspec&format=csv')" = "200" ]; do sleep 600; done
echo "SkyServer is back: $(date)"
$PY 14c_v2_epochs.py | tail -1
$PY 03_spectra_inventory.py data/v2_targets_positions.csv v2pub ra dec name 2>&1 | grep -v -i warning | tail -3
CLAGN_SEL=v2 $PY 03d_fetch_spectra.py data/targets_sep23_v2.csv data/targets_oct26_v2.csv data/targets_oct27_v2.csv 2>&1 | grep -v -i warning | tail -1
CLAGN_SEL=v2 $PY 07_make_webpage.py 2>&1 | grep wrote
./publish_v2.sh
echo "done: $(date)"

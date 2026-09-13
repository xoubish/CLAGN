#!/bin/zsh
# Wait for SPARCL (DESI spectra) to answer, then refetch the official targets' spectra, rebuild the page and refresh the official files.
cd "$(dirname "$0")"; PY=/opt/anaconda3/bin/python
until $PY -c "from sparcl.client import SparclClient; SparclClient()" > /dev/null 2>&1; do sleep 600; done
echo "SPARCL is back: $(date)"
CLAGN_SEL=v2 $PY 03d_fetch_spectra.py data/targets_sep23_v2.csv data/targets_oct26_v2.csv data/targets_oct27_v2.csv 2>&1 | grep -v -i warning | tail -1
CLAGN_SEL=v2 $PY 07_make_webpage.py 2>&1 | grep wrote
./publish_v2.sh
echo "done: $(date)"

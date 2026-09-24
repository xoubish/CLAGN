"""Run the documented four-channel NGPS/PypeIt reduction on staged September data.

Use sep23_data/reduction_20260924/.venv/bin/python to execute this script.
Raw originals are inventoried by 53_ngps_reduction_audit.py and never modified here.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT/'sep23_data/reduction_20260924'
NIGHT = BASE/'work/20260924'
LOGS = NIGHT/'logs'
BIN = BASE/'.venv/bin'
for directory in [BASE/'cache', BASE/'config', BASE/'cache/matplotlib']:
    directory.mkdir(parents=True, exist_ok=True)
for key, value in dict(MPLBACKEND='Agg', MPLCONFIGDIR=str(BASE/'cache/matplotlib'),
    XDG_CACHE_HOME=str(BASE/'cache'), XDG_CONFIG_HOME=str(BASE/'config'),
    NGPS_PIPELINE_CONFIG=str(BASE/'ngps_pipeline.toml'), NGPS_FINALIZE_WORKERS='4',
    OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', VECLIB_MAXIMUM_THREADS='1',
    NUMEXPR_NUM_THREADS='1').items():
    os.environ[key] = value
LOGS.mkdir(parents=True, exist_ok=True)


def run(cmd, log, cwd=NIGHT):
    print('RUN', ' '.join(map(str, cmd)), flush=True)
    with (LOGS/log).open('a') as f:
        f.write('\nCOMMAND: '+' '.join(map(str, cmd))+'\n')
        f.flush()
        result = subprocess.run(list(map(str, cmd)), cwd=cwd, stdout=f, stderr=subprocess.STDOUT)
    if result.returncode:
        raise RuntimeError(f'Exit {result.returncode}; see {LOGS/log}')


def pypeit_file(ch):
    return NIGHT/f'setup_{ch}/p200_ngps_{ch}_A/p200_ngps_{ch}_A.pypeit'


def setup(ch):
    run([BIN/'pypeit_setup','-s',f'p200_ngps_{ch}','-r',NIGHT/'raw',
         '-d',NIGHT/f'setup_{ch}','-c','A','-o'],f'driver_setup_{ch}.log')
    pf = pypeit_file(ch)
    from pypeit.inputfiles import PypeItFile
    parsed = PypeItFile.from_file(str(pf))
    report = []
    for row in parsed.data:
        report.append({k:str(row[k]) for k in parsed.data.colnames})
    (BASE/f'audit/setup_{ch}.json').write_text(json.dumps(report,indent=2))
    print(ch, 'setup:', len(parsed.data), 'frames; inspect before running', flush=True)


def reduce(ch):
    pf = pypeit_file(ch)
    run([BIN/'run_pypeit',pf.name],f'driver_run_{ch}.log',cwd=pf.parent)
    print(ch,'reduction finished',flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('stage',choices=['setup','reduce'])
    p.add_argument('--channels',default='u,g,r,i')
    args=p.parse_args()
    channels=args.channels.split(',')
    assert all(c in 'ugri' and len(c)==1 for c in channels)
    with ThreadPoolExecutor(max_workers=len(channels)) as pool:
        results=list(pool.map(setup if args.stage=='setup' else reduce,channels))


if __name__=='__main__':
    main()

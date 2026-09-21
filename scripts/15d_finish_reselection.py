"""Finish a running W1 batch automatically; never publish or assign science weights.

Run with the same --rmax used for 15c fetch. Progress and failure states are local.
"""
from pathlib import Path
import argparse,json,subprocess,sys,time,os

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'data/reselection_2026-09-20'


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--rmax',type=float,default=20.5);a=ap.parse_args()
    batch=OUT/f'wise_r{str(a.rmax).replace(".","p")}'
    status=OUT/'pipeline_status.json'
    def write(stage,**extra):
        value={'stage':stage,'updated_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),'rmax_acquisition_batch':a.rmax,
               'scientific_weights':'not assigned; await spectral-state review',**extra}
        tmp=status.with_suffix('.tmp');tmp.write_text(json.dumps(value,indent=2));tmp.replace(status)
        print(value,flush=True)
    write('waiting for W1 acquisition')
    while not (batch/'status.json').exists():time.sleep(30)
    fetched=json.loads((batch/'status.json').read_text())
    if fetched['failed']:
        write('acquisition incomplete',failures=fetched['failed']);return 1
    stages=[('projecting onto saved W1 manifold',['15c_project_expanded_pool.py','project','--rmax',str(a.rmax)]),
            ('refreshing candidate tables',['15b_rebuild_observing_pool.py','export']),
            ('writing pool report',['15e_pool_report.py'])]
    for label,args in stages:
        write(label)
        with open(OUT/(args[0]+'.log'),'w') as log:
            r=subprocess.run([sys.executable,str(ROOT/'scripts'/args[0]),*args[1:]],cwd=ROOT,stdout=log,stderr=subprocess.STDOUT)
        if r.returncode:
            write('failed',step=label,returncode=r.returncode);return r.returncode
    write('acquisition batch and pool rebuild complete; spectral weighting pending')
    return 0


if __name__=='__main__':raise SystemExit(main())

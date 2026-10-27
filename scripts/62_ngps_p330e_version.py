"""Recalibrate this night's extracted spectra using one standard (default P330E) for every flux scale.

Preserves the original mixed-standard products and atmospheric transmission.
"""
import argparse
import importlib
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import shutil

env=importlib.import_module('54_reduce_ngps_sep23')

def run_channel(ch,standard):
    env.run([env.BIN/'python',Path(__file__).resolve(),'--channel',ch,'--standard',standard],
        f'{standard.lower()}_flux_coadd_{ch}.log')

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--channel',choices=list('ugri'))
    p.add_argument('--standard',default='P330E')
    args=p.parse_args();variant=args.standard.lower()
    if args.channel:
        importlib.import_module('58_ngps_central_products').flux_and_coadd(
            args.channel,standard=args.standard,variant=variant)
        return
    with ThreadPoolExecutor(max_workers=4) as pool:list(pool.map(lambda ch:run_channel(ch,args.standard),'ugri'))
    deliver=importlib.import_module('59_ngps_deliver')
    deliver.OUT=env.BASE/f'products_{variant}'
    deliver.COADDS=env.NIGHT/f'coadded_channels_{variant}'
    deliver.ASSIGNMENTS_SUFFIX='_'+variant
    # Only the flux standard changes. The already validated atmospheric
    # transmission template is retained, and does not set the flux scale.
    deliver.main()
    for name in ['extraction_qa.csv']:
        shutil.copy2(env.BASE/'products/calibration'/name,deliver.OUT/'calibration'/name)
    shutil.copy2(env.BASE/'products/plots/extraction_checks.pdf',deliver.OUT/'plots/extraction_checks.pdf')

if __name__=='__main__':main()

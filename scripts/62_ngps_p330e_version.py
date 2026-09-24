"""Recalibrate this night's extracted spectra using P330E for every flux scale.

Preserves the original mixed-standard products and atmospheric transmission.
"""
import argparse
import importlib
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import shutil

env=importlib.import_module('54_reduce_ngps_sep23')

def run_channel(ch):
    env.run([env.BIN/'python',Path(__file__).resolve(),'--channel',ch],
        f'p330e_flux_coadd_{ch}.log')

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--channel',choices=list('ugri'))
    args=p.parse_args()
    if args.channel:
        importlib.import_module('58_ngps_central_products').flux_and_coadd(
            args.channel,standard='P330E',variant='p330e')
        return
    with ThreadPoolExecutor(max_workers=4) as pool:list(pool.map(run_channel,'ugri'))
    deliver=importlib.import_module('59_ngps_deliver')
    deliver.OUT=env.BASE/'products_p330e'
    deliver.COADDS=env.NIGHT/'coadded_channels_p330e'
    deliver.ASSIGNMENTS_SUFFIX='_p330e'
    # Only the flux standard changes. The already validated BD+28 atmospheric
    # transmission template is retained, and does not set the flux scale.
    deliver.main()
    for name in ['extraction_qa.csv']:
        shutil.copy2(env.BASE/'products/calibration'/name,deliver.OUT/'calibration'/name)
    shutil.copy2(env.BASE/'products/plots/extraction_checks.pdf',deliver.OUT/'plots/extraction_checks.pdf')

if __name__=='__main__':main()

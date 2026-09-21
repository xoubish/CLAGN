"""Resume the complete parent search after the two W1 acquisitions finish.

No commit, push, or publication. Each stage logs to the ignored research
directory. A failed stage stops the chain and records the failure explicitly.
The currently rendered review remains available while acquisition proceeds.
"""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime,timezone
import hashlib,json,os,subprocess,sys,time
import pandas as pd

ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'data/reselection_2026-09-20'
STATE=OUT/'completion_pipeline_status.json'
ENV=os.environ|{'MPLCONFIGDIR':'/private/tmp/clagn_mpl','NUMBA_NUM_THREADS':'8','NUMBA_CACHE_DIR':'/private/tmp/clagn_numba'}


def state(stage,**extra):
    STATE.write_text(json.dumps(dict(stage=stage,updated_utc=datetime.now(timezone.utc).isoformat(),**extra),indent=2))
    print(stage,flush=True)


def run(script,*args):
    label=Path(script).stem+'_'+hashlib.sha256(str(args).encode()).hexdigest()[:8]
    with (OUT/f'completion_{label}.log').open('w') as log:
        subprocess.run([sys.executable,str(ROOT/'scripts'/script),*map(str,args)],cwd=ROOT,env=ENV,stdout=log,stderr=subprocess.STDOUT,check=True)


def main():
    batches=['wise_r19_three_night','wise_r19_rolling_new']
    while True:
        counts={b:len(list((OUT/b).glob('pixel_*.parquet'))) for b in batches}
        complete=True
        for b in batches:
            f=OUT/b/'status.json'
            if not f.exists():complete=False;continue
            status=json.loads(f.read_text())
            if status.get('failed'):raise RuntimeError(f'W1 acquisition has recorded failures: {b}')
            complete&='completed' in status
        if complete:break
        state('Waiting for complete W1 acquisition',cached_partitions=counts)
        time.sleep(30)
    state('Projecting the complete W1 input cohorts')
    for b in batches:run('15c_project_expanded_pool.py','project','--rmax','19','--batch',b)
    state('Merging rolling identities and completed light curves')
    run('36_merge_rolling_parent.py')
    curves=pd.concat([pd.read_parquet(f) for b in batches for f in sorted((OUT/b).glob('pixel_*.parquet'))],ignore_index=True)
    curves.drop_duplicates(['name','time','flux','err']).to_parquet(OUT/'three_night_w1.parquet',index=False)
    run('32_update_expanded_manifold.py','--select')
    state('Screening fields and preparing time-balanced lists')
    run('24_screen_neighbours.py');run('33_prepare_three_night_lists.py')
    tasks=[('28_three_night_acquisition.py','spectra','--release','dr20','--workers','6'),
           ('28_three_night_acquisition.py','ztf','--workers','6'),('17_fetch_review_spectra.py',),
           ('35_desi_direct.py',),('18_fetch_review_images.py','--wide')]
    for attempt in range(4):
        targets=pd.read_csv(OUT/'compact_review_objects.csv')
        targets.to_csv(OUT/'completed_parent_prepared_targets.csv',index=False)
        state('Retrieving spectra, ZTF and SDSS cutouts for the new lists',targets=len(targets),baseline_validation_round=attempt+1)
        with ThreadPoolExecutor(max_workers=5) as pool:
            jobs=[pool.submit(run,*task) for task in tasks]
            for job in jobs:job.result()
        run('38_three_night_acquisition_audit.py')
        audit=pd.read_csv(OUT/'prepared_acquisition_audit.csv')
        if audit.hbeta_continuum_baseline.all():break
        state('Replacing targets whose retrieved spectra cannot support the Hbeta comparison')
        run('33_prepare_three_night_lists.py')
    else:raise RuntimeError('Baseline validation still incomplete after four preparation rounds; inspect the audit.')
    existing=set()
    files=list((ROOT/'data').glob('neowise_visits_*.csv'))+list(OUT.glob('neowise_visits_prepared*.csv'))
    for f in files:
        try:existing.update(pd.read_csv(f,usecols=['name']).name)
        except (pd.errors.EmptyDataError,ValueError):continue
    missing=targets[~targets.name.isin(existing)]
    if len(missing):
        tag='prepared_full_'+hashlib.sha256('|'.join(sorted(missing.name)).encode()).hexdigest()[:10]
        path=OUT/f'{tag}_targets.csv';missing.to_csv(path,index=False)
        state('Retrieving missing final-release NEOWISE histories',targets=len(missing))
        run('03c_neowise_now.py',path,tag,'--output-dir',OUT)
    state('Auditing observations, sensitivity and spectral histories')
    run('15f_compact_spectral_audit.py');run('30_three_night_sensitivity.py');run('31_three_night_review.py')
    run('39_airmass_options.py')
    run('38_three_night_acquisition_audit.py','--require-baselines')
    state('Rebuilding local and public-data review pages')
    config=json.loads((OUT/'review_selection.json').read_text())
    config['w1_parent_search_complete']=True
    (OUT/'review_selection.json').write_text(json.dumps(config,indent=2))
    run('16_candidate_webpage.py')
    summary=pd.read_csv(OUT/'prepared_night_summary.csv')
    audit=json.loads((OUT/'prepared_acquisition_audit.json').read_text())
    report='# Completed-parent prepared review\n\n'+datetime.now(timezone.utc).isoformat()+'\n\n'
    report+='The W1 acquisition/projection batches for the current bright, observable parent are complete. Invalid or unmatched histories remain in the projection audit. Final broad-line classification, science weights and the observing sequence still require review. No files were committed, pushed or published.\n\n'
    report+=summary.to_markdown(index=False)+'\n\n'
    report+='[Complete local page](candidate_review_local.html) · [Prepared targets](compact_review_objects.csv) · [Acquisition audit](prepared_acquisition_audit.csv) · [Science and sensitivity](three_night_review/science_and_sensitivity.csv)\n\n'
    report+='Acquisition summary (empty or failed archive results are retained explicitly):\n\n```json\n'+json.dumps(audit,indent=2)+'\n```\n'
    (OUT/'FULL_PARENT_COMPLETION.md').write_text(report)
    state('Prepared review rebuilt; inspect acquisition audit for archive limitations',complete_parent_acquisition=True,targets=len(targets),committed=False,pushed=False)


if __name__=='__main__':
    try:main()
    except Exception as exc:
        state('Stopped on a recorded failure',reason=type(exc).__name__,detail=str(exc)[:500])
        raise

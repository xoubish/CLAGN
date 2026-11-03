"""Refresh public alert photometry separately from all original science inputs."""
import importlib.util
from concurrent.futures import ThreadPoolExecutor,as_completed
import pandas as pd
from paper_protocol import ROOT,OUT


def main():
    spec=importlib.util.spec_from_file_location('alert_review',ROOT/'scripts/43_public_alert_review.py')
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    mod.CACHE=OUT/'optical_alert_refresh';mod.CACHE.mkdir(parents=True,exist_ok=True)
    sample=pd.read_csv(OUT/'pilot_sample.csv');rows=[]
    with ThreadPoolExecutor(max_workers=3) as ex:
        futures=[ex.submit(mod.one,t) for t in sample.itertuples(index=False)]
        for f in as_completed(futures):
            rows.append(f.result());print(rows[-1],flush=True)
    pd.DataFrame(rows).to_csv(OUT/'optical_refresh_summary.csv',index=False)


if __name__=='__main__':main()

"""Check the reproducible pilot bundle and write a local navigation page."""
import hashlib,importlib.metadata,json,os,subprocess,sys
from datetime import datetime,timezone
from pathlib import Path
import nbformat
import numpy as np
import pandas as pd
from pypdf import PdfReader
from paper_protocol import ROOT,HERE,OUT


def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    pilot=pd.read_csv(OUT/'pilot_sample.csv');ep=pd.read_csv(OUT/'pilot_epochs.csv')
    fits=pd.read_csv(OUT/'fit_results.csv');alt=pd.read_csv(OUT/'alternate_standard_fits.csv')
    assert len(pilot)==15 and pilot.name.is_unique
    assert len(ep)==43 and ep.epoch_id.is_unique
    assert not fits.duplicated(['epoch_id','variant']).any()
    assert len(fits)==129 and fits.status.eq('fit').sum()==126
    assert set(fits[fits.status.ne('fit')].epoch_id)=={'P10381_DESI_59523'}
    assert fits.loc[fits.status.eq('fit'),'optimizer_success'].all()
    assert len(alt)==60 and alt.status.eq('fit').all() and alt.optimizer_success.all()
    assert len(list((OUT/'fits').glob('*.npz')))==186
    # Check that every latest reference is exactly the predeclared date-selected
    # reference, not the latest fit that happened to succeed.
    pairs=pd.read_csv(OUT/'fit_comparisons.csv')
    for r in pairs[pairs.comparison.eq('latest_archive')].itertuples(index=False):
        visit=ep[ep.epoch_id.eq(r.epoch_id)].iloc[0]
        old=ep[ep.name.eq(r.name)&ep.source.ne('NGPS')&ep.max_mjd.lt(visit.min_mjd)].sort_values('max_mjd')
        assert r.reference_epoch==old.iloc[-1].epoch_id
    stress=pd.read_csv(OUT/'calibration_model_envelope.csv')
    env=pd.read_csv(OUT/'model_sensitivity_envelope.csv')
    merged=stress.merge(env,on=['name','epoch_id','reference_epoch','comparison'],suffixes=('_stress','_model'))
    assert len(merged)==len(env)
    for key in ['hb_ew_total_cont','hb_oiii']:
        assert np.all(merged[key+'_ratio_min_stress']<=merged[key+'_ratio_min_model']+1e-10)
        assert np.all(merged[key+'_ratio_max_stress']>=merged[key+'_ratio_max_model']-1e-10)
    checked_inputs=0
    for filename in ['spectrum_inputs.json','exposure_inputs.json']:
        for rel,meta in json.loads((OUT/filename).read_text()).items():
            assert digest(ROOT/rel)==meta['sha256'],rel
            checked_inputs+=1
    for rel,sha in json.loads((OUT/'protocol_snapshot.json').read_text())['inputs'].items():
        assert digest(ROOT/rel)==sha,rel
        checked_inputs+=1
    for row in json.loads((OUT/'templates/provenance.json').read_text())['files']:
        assert digest(OUT/'templates'/row['file'])==row['sha256']
    assert digest(HERE/'pilot_fits.py')==json.loads((OUT/'fit_method.json').read_text())['source_sha256']
    assert len(PdfReader(OUT/'pilot_fit_atlas.pdf').pages)==15
    for name in pilot.name:assert (OUT/'figures'/f'{name}_decomposition.png').exists()
    assert len(pd.read_csv(OUT/'independent_exposure_indices.csv'))==24
    assert len(pd.read_csv(OUT/'injection_recovery.csv'))==72
    books=[]
    for prefix in ['06','07']:
        path=next((HERE/'notebooks').glob(prefix+'_*.ipynb'))
        nb=nbformat.read(path,as_version=4);cells=[c for c in nb.cells if c.cell_type=='code']
        assert all(c.execution_count is not None for c in cells)
        assert not any(o.output_type=='error' for c in cells for o in c.outputs)
        html=HERE/'outputs/html'/f'{path.stem}.html';assert html.exists()
        books.append(dict(file=path.name,executed_cells=len(cells),sha256=digest(path),html_sha256=digest(html)))
    envvars=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MPLCONFIGDIR='/private/tmp/palomar_mpl')
    test=subprocess.run([sys.executable,'-m','unittest','discover','-s',str(HERE),'-p','test_*.py','-v'],text=True,capture_output=True,env=envvars)
    (OUT/'unit_tests.txt').write_text(test.stdout+test.stderr)
    assert test.returncode==0,test.stderr
    manifest=dict(validated_utc=datetime.now(timezone.utc).isoformat(),pilot_sources=len(pilot),selected_epochs=len(ep),
        fitted_epochs=fits[fits.status.eq('fit')].epoch_id.nunique(),primary_fits=126,alternate_fits=len(alt),
        source_atlas_pages=15,independent_exposures=8,injection_trials=72,unit_test_returncode=test.returncode,
        verified_input_hash_entries=checked_inputs,notebooks=books,
        python=sys.version,packages={p:importlib.metadata.version(p) for p in ['numpy','pandas','scipy','astropy','matplotlib','nbformat','nbclient','nbconvert','requests','pyarrow','pypdf']},
        code_sha256={p.name:digest(p) for p in HERE.glob('*.py')},
        caveat='Artifact consistency and limited recovery checks; not publication-level physical validation.')
    (OUT/'pilot_validation_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    links=''.join(f'<li><a href="../html/{Path(b["file"]).stem}.html">{b["file"][:2]}: {"Sample and timing" if b["file"].startswith("06") else "Fits, plots, and scientific results"}</a></li>' for b in books)
    sources=''.join(f'<li><a href="figures/{name}_decomposition.png">{name}</a></li>' for name in pilot.name)
    (OUT/'index.html').write_text(f'''<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Palomar paper pilot</title><style>body{{font:17px/1.6 system-ui;max-width:980px;margin:3rem auto;padding:0 1.5rem;color:#223}}a{{color:#16717b}}li{{margin:.35rem 0}}.sources{{columns:3}}img{{max-width:100%}}</style>
<h1>Palomar paper pilot</h1><p>15 sources · 43 selected epochs · 42 fitted epochs · three host/Fe II models · calibration and exposure checks.</p>
<p><strong>Priority:</strong> validate P21422, P9506 and P7947 with joint Balmer fits, then measure the full observed sample. Dust physics waits for JWST. This pilot does not establish a selection advantage or new changing-look classification.</p>
<ul>{links}<li><a href="pilot_fit_atlas.pdf">15-source fit atlas (PDF)</a></li><li><a href="../../PALOMAR_PAPER.md">Publication plan and reproduction commands</a></li><li><a href="pilot_validation_manifest.json">Validation and provenance</a></li></ul>
<h2>Individual fit pages</h2><ul class="sources">{sources}</ul><h2>Primary-standard model sensitivity</h2><p>Bars span model choices, not confidence intervals.</p><img src="figures/latest_archive_model_envelope.png" alt="Broad H beta changes across model choices">
</html>''')
    print(json.dumps({k:v for k,v in manifest.items() if k not in ['python','packages','code_sha256','notebooks']},indent=2))


if __name__=='__main__':main()

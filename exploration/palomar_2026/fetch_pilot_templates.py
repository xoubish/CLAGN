"""Download documented public data templates only; no third-party code executes."""
from pathlib import Path
import hashlib
import json
import subprocess
import requests

HERE=Path(__file__).resolve().parent
DEST=HERE/'outputs/palomar_paper/templates'
FILES=['LICENSE.md','src/pyqsofit/fe_optical.txt']+[
    'src/pyqsofit/bc03/'+n for n in ['7_ssp_11Gyr_z02.spec.gz','13_ssp_5Gyr_z02.spec.gz',
                                   '19_ssp_1.4Gyr_z02.spec.gz','28_ssp_290Myr_z02.spec.gz']]


def main():
    DEST.mkdir(parents=True,exist_ok=True)
    manifest=DEST/'provenance.json'
    prior=json.loads(manifest.read_text()) if manifest.exists() else None
    version=prior['repository_commit'] if prior else subprocess.check_output(['git','ls-remote','https://github.com/legolason/PyQSOFit.git','HEAD'],text=True).split()[0]
    rows=[]
    for name in FILES:
        url=f'https://raw.githubusercontent.com/legolason/PyQSOFit/{version}/{name}'
        p=DEST/Path(name).name
        if p.exists() and prior:
            expected=next(r['sha256'] for r in prior['files'] if r['file']==p.name)
            if hashlib.sha256(p.read_bytes()).hexdigest()!=expected:
                raise RuntimeError(f'Existing template changed: {p}. Preserve it before refreshing.')
        elif p.exists():
            raise RuntimeError(f'Existing template has no provenance: {p}')
        if not p.exists():
            r=requests.get(url,timeout=45);r.raise_for_status()
            p.write_bytes(r.content)
        rows.append(dict(file=p.name,url=url,sha256=hashlib.sha256(p.read_bytes()).hexdigest(),bytes=p.stat().st_size))
        print(p.name,p.stat().st_size,flush=True)
    (DEST/'provenance.json').write_text(json.dumps(dict(repository_commit=version,files=rows,
        references=['Guo, Shen & Wang 2018 ascl:1809.008','Bruzual & Charlot 2003','Boroson & Green 1992'],
        note='Templates distributed with PyQSOFit; local pilot fitter is separate; repository license retained.'),indent=2)+'\n')


if __name__=='__main__':main()

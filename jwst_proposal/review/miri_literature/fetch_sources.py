"""Download a small, explicitly listed set of public reference PDFs.

These are literature/reference files; this does not compile the proposal.
"""
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import hashlib
import json
import requests

HERE = Path(__file__).resolve().parent
SOURCES = {
    'program_3696.pdf': 'https://www.stsci.edu/jwst/phase2-public/3696.pdf',
    'program_8245.pdf': 'https://www.stsci.edu/jwst/phase2-public/8245.pdf',
    'masterson_2503.08647.pdf': 'https://arxiv.org/pdf/2503.08647v2',
    'sanchez_saez_2607.00921.pdf': 'https://arxiv.org/pdf/2607.00921v1',
    'program_1670.pdf': 'https://www.stsci.edu/jwst/phase2-public/1670.pdf',
    'program_1328.pdf': 'https://www.stsci.edu/jwst/phase2-public/1328.pdf',
    'program_2016.pdf': 'https://www.stsci.edu/jwst/phase2-public/2016.pdf',
    'ramos_almeida_2512.02629.pdf': 'https://arxiv.org/pdf/2512.02629',
    'gonzalez_martin_2504.01103.pdf': 'https://arxiv.org/pdf/2504.01103',
    'donnan_2402.17479.pdf': 'https://arxiv.org/pdf/2402.17479',
    'armus_2209.13125.pdf': 'https://arxiv.org/pdf/2209.13125',
}


def download(item):
    name, url = item
    try:
        path = HERE / name
        if path.exists() and path.read_bytes().startswith(b'%PDF'):
            data = path.read_bytes()
            return dict(file=name, url=url, bytes=len(data),
                        sha256=hashlib.sha256(data).hexdigest(), status='cached')
        r = requests.get(url, timeout=40)
        r.raise_for_status()
        if not r.content.startswith(b'%PDF'):
            raise ValueError('Response is not a PDF')
        (HERE / name).write_bytes(r.content)
        return dict(file=name, url=url, final_url=r.url, bytes=len(r.content),
                    sha256=hashlib.sha256(r.content).hexdigest(), status='downloaded')
    except Exception as e:
        return dict(file=name, url=url, status='failed', error=str(e))


if __name__ == '__main__':
    with ThreadPoolExecutor(max_workers=4) as pool:
        records = list(pool.map(download, SOURCES.items()))
    (HERE / 'downloads.json').write_text(json.dumps(records, indent=2) + '\n')
    for r in records:
        print(r['file'], r['status'], r.get('bytes', r.get('error')))

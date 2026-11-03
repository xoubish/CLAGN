"""Execute notebooks and export standalone HTML, failing on any cell error."""
import os
from pathlib import Path
import json
import time
import sys
import nbformat
from nbclient import NotebookClient
from nbconvert import HTMLExporter

HERE = Path(__file__).resolve().parent
for key, suffix in [('MPLCONFIGDIR', 'matplotlib'), ('JUPYTER_RUNTIME_DIR', 'jupyter'), ('IPYTHONDIR', 'ipython')]:
    directory = HERE/'.cache'/suffix
    directory.mkdir(parents=True, exist_ok=True)
    os.environ[key] = str(directory)


def main():
    out = HERE/'outputs/html'
    out.mkdir(parents=True, exist_ok=True)
    result_path = HERE/'outputs/notebook_execution.json'
    results = json.loads(result_path.read_text()) if result_path.exists() and sys.argv[1:] else []
    for path in sorted((HERE/'notebooks').glob('0[0-5]_*.ipynb')):
        if sys.argv[1:] and path.name[:2] not in sys.argv[1:]:
            continue
        start = time.monotonic()
        print('Executing '+path.name, flush=True)
        book = nbformat.read(path, as_version=4)
        client = NotebookClient(book, timeout=900, kernel_name='python3', resources={'metadata': {'path': str(path.parent)}})
        client.execute()
        nbformat.write(book, path)
        exporter = HTMLExporter()
        body, _ = exporter.from_notebook_node(book)
        # Notebook artifact links are one directory up; HTML lives in outputs/html.
        body = body.replace('../outputs/', '../').replace('../../../reference/', '../../../../reference/').replace('../../../observing/', '../../../../observing/').replace('../../../jwst_proposal/', '../../../../jwst_proposal/')
        (out/f'{path.stem}.html').write_text(body)
        rows = [c for c in book.cells if c.cell_type == 'code']
        results = [r for r in results if r['notebook'] != path.name]
        results.append(dict(notebook=path.name, code_cells=len(rows),
            executed_cells=sum(c.execution_count is not None for c in rows),
            errors=sum(o.output_type=='error' for c in rows for o in c.outputs), seconds=round(time.monotonic()-start, 1)))
        result_path.write_text(json.dumps(sorted(results, key=lambda r:r['notebook']), indent=2)+'\n')
        print(results[-1], flush=True)


if __name__ == '__main__':
    main()

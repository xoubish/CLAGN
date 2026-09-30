"""Regenerate the adopted literature-based Figure 2 without compiling LaTeX."""
from pathlib import Path
import hashlib
import json
import shutil
import subprocess
import sys

BASE = Path(__file__).resolve().parent
DRAFT = BASE / 'review/fig2_three_panel'


def install():
    source = DRAFT / 'fig2_three_panel_draft.pdf'
    target = BASE / 'fig2_miri.pdf'
    shutil.copy2(source, target)
    metadata = {
        'generator': 'make_fig2_miri.py',
        'figure_type': 'Published spectral, spatial and temporal examples',
        'source_pdf': str(source.relative_to(BASE)),
        'source_sha256': hashlib.sha256(source.read_bytes()).hexdigest(),
        'caption': 'review/fig2_three_panel/caption.tex',
        'panel_a_source': 'https://arxiv.org/abs/2402.17479',
        'panel_a_figure': 8,
        'panel_a_method': 'Digitized native plot pixels, calibrated logarithmic axes, redrawn as a wide vector plot.',
        'panel_a_data': 'inputs/fig2_ngc7469_digitized.csv',
        'panel_b_source': 'https://arxiv.org/abs/2504.01103v2',
        'panel_b_publication': 'Gonzalez-Martin et al. 2025, MNRAS, 539, 2158; DOI 10.1093/mnras/staf573',
        'panel_b_version_note': 'Figure 4 from the April 2026 arXiv revision of the 2025 paper, not a different paper.',
        'panel_b_figure': 4,
        'panel_c_source': 'https://arxiv.org/abs/2011.07638',
        'panel_c_figure': 15,
        'checks': 'review/fig2_three_panel/checks.json',
        'interpretation': 'Published diagnostic precedents, not target predictions or a lag measurement from the proposed single visit.',
        'proposal_pdf_compiled': False,
        'source_embedding': 'Vector top plot and lossless cropped lower figures; source pages and off-crop text are not embedded.',
    }
    (BASE / 'inputs/fig2_miri_provenance.json').write_text(json.dumps(metadata, indent=2) + '\n')
    print(target)


def main():
    for script in ['make_fig2_spectrum.py', 'make_fig2_three_panel.py']:
        subprocess.run([sys.executable, str(BASE / 'review' / script)], check=True)
    install()


if __name__ == '__main__':
    main()

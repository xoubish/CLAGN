"""Build a local comparison of published figures; never modify proposal assets."""
from pathlib import Path
from html import escape
import hashlib
import json

HERE = Path(__file__).resolve().parent
PROPOSAL = HERE.parents[2]
protected = [PROPOSAL / name for name in ('proposal.tex', 'fig2_miri.pdf', 'make_fig2_miri.py')]
before = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in protected}

cards = [
    dict(title='1. Different wavelengths remember different parts of the history',
         source='Lyu & Rieke (2021), Figures 10 and 16', url='https://arxiv.org/abs/2011.07638',
         images=['lyu_2011.07638_p12.png', 'lyu_2011.07638_p20.png'], rank='Best physical story',
         shows='The NGC 4151 light curves are fitted with two dust components; their relative contributions change with wavelength. A second figure connects the variable spectrum to inferred spatial scales.',
         borrow='Keep a measured history above a spectrum, but use the same colors to connect hot and warm emission to their delayed responses. A compact time-weighting inset could explain why a flare, recovery, or sustained decline matters even when two endpoint fluxes agree.',
         limit='This paper has long, repeated light curves. Our single MRS visit will constrain history-conditioned spectra, not directly measure a warm-dust lag. Its measured delays and radii are not values for our sample.'),
    dict(title='2. Show the response to past illumination explicitly',
         source='Almeyda et al. (2020), Figure 12', url='https://arxiv.org/abs/2002.12823',
         images=['almeyda_2002.12823_p16.png'], rank='Companion to idea 1',
         shows='A grid of calculated response functions at 3.6, 10, and 30 microns shows how wavelength, viewing angle, and radiative transfer affect the distribution of delays.',
         borrow='Replace the large grid with two or three illustrative response kernels beside one observed history. Then mark the MIRI epoch and explain which earlier illumination contributes to the predicted spectrum.',
         limit='Kernels require physical assumptions and a luminosity-dependent timescale. A universal rule that longer wavelength always means one particular longer lag would be misleading. The 30-micron example is outside our proposed MRS coverage; use covered wavelengths in an adaptation.'),
    dict(title='3. Turn the spectrum into identifiable emission components',
         source='Donnan et al. (2024), Figure 8', url='https://arxiv.org/abs/2402.17479',
         images=['donnan_2402.17479_p13.png'], rank='Best near-term visual improvement',
         shows='Each spectrum sits beside an inferred temperature–extinction distribution. Selected cases separate contributions associated with hot AGN dust, warm AGN dust, and star formation.',
         borrow='Use one large spectrum with clearly colored hot-dust, warm-dust, and host components. Mark SPHEREx and MIRI coverage and the silicate/PAH regions. This makes the added information visible without relying on a flux-change statistic.',
         limit='These are model-dependent decompositions using joint NIRSpec/MIRI data. SPHEREx does not provide the same spectral resolution. We should not display a recovered temperature–extinction map or forecast component precision without testing our own data and model.'),
    dict(title='4. Pair the two silicate features in a physical diagnostic plane',
         source='Sirocky et al. (2008), Figure 7', url='https://arxiv.org/abs/0801.4776',
         images=['sirocky_0801.4776_p34.png'], rank='Strong idea; needs new calculations',
         shows='The 18-micron and 10-micron feature strengths are plotted together with model tracks for different dust distributions and observed ULIRGs.',
         borrow='A small S18 versus S9.7 inset could show how measuring both features constrains models beyond the broad-band continuum. A forecast point would need its joint error ellipse and the models evaluated with the same measurement definition.',
         limit='The published systems and grids are not our target sample. Our current local-continuum indices differ from the paper’s feature definitions. Do not place our mock point on these tracks or label separate regions as rising/fading or evolving/fixed without a consistent calculation.'),
    dict(title='5. Make the spatial information visibly useful',
         source='González-Martín et al., Figures 4 and 5 (downloaded arXiv version)', url='https://arxiv.org/abs/2504.01103',
         images=['gonzalez_martin_2504.01103_p9.png', 'gonzalez_martin_2504.01103_p12.png'], rank='Useful small supporting inset',
         shows='Total, point-source, and extended spectra are linked directly to their images. The maps show what the spatial decomposition removes from the central extraction.',
         borrow='A compact aperture/PSF diagram next to central and surrounding spectra would explain how we constrain host contamination. This would give the existing “host contribution” column a concrete visual counterpart.',
         limit='These nearby systems resolve finer physical scales than much of our sample. The point-source spectrum need not be pure AGN emission, and extended emission need not be purely star formation. For our targets, a spatial mock must use their redshifts, wavelength-dependent PSF, and common field of view.'),
    dict(title='6. Use PAHs as a check on the surrounding dust',
         source='García-Bernete et al. (2024), GATOS V, Figure 5', url='https://arxiv.org/abs/2409.05686',
         images=['gatos_2409.05686_p8.png'], rank='Secondary science; not the main Figure 2',
         shows='PAH band ratios from nuclear, outflow, and star-forming regions are compared with molecular-size and ionization grids.',
         borrow='Explain that PAH features characterize the dust environment and help assess contamination; they are not merely wavelengths labeled above the continuum.',
         limit='We have not forecast PAH detections or spatially resolved band-ratio precision. A whole PAH diagnostic panel would compete with the dust-memory question and require additional exposure/measurement work.'),
]

intro = ('The current Figure 2 remains in place. Four additional papers were downloaded, and '
         'previously collected MRS papers received a closer visual review. These are design precedents, '
         'not forecasts for our targets. Click an image for its full source page; original captions remain visible.')
recommendation = ('Keep the current figure as the baseline. The most promising alternative combines the '
                  'history-to-response logic of 1–2 with the readable spectral components of 3. '
                  'Retain silicate measurements and host checks as supporting constraints. '
                  'Before replacing the baseline, the alternative must make that connection clearer at proposal size '
                  'and distinguish illustrative physics from tested precision. A paired-silicate model forecast is '
                  'interesting but would need additional consistent calculations; it is not ready for insertion.')
parts = ['<!doctype html><html lang="en"><meta charset="utf-8">',
         '<meta name="viewport" content="width=device-width,initial-scale=1">',
         '<title>Figure 2: further literature ideas</title>',
         '<style>body{font:17px/1.55 system-ui,sans-serif;color:#20303b;background:#f3f5f7;margin:0} '
         'main{max-width:1250px;margin:auto;padding:30px} h1{font-size:30px} h2{font-size:22px;line-height:1.3} '
         'a{color:#145e8b} .grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(340px,1fr));gap:24px} '
         'article,.summary{background:white;border:1px solid #dce2e7;border-radius:12px;padding:23px;margin:20px 0} '
         'article p{font-size:16px} .rank{color:#876000;font-weight:650} '
         '.images{display:flex;gap:10px;align-items:start}.images a{flex:1;min-width:0} '
         'img{width:100%;max-height:440px;object-fit:contain;background:#f8f8f8} '
         'small{color:#54616b} @media(max-width:420px){main{padding:12px}.grid{display:block}} </style><main>',
         '<h1>Further ideas for Figure 2</h1><p>'+escape(intro)+'</p>',
         '<div class="summary"><b>Recommendation</b><p>'+escape(recommendation)+'</p>',
         '<a href="../../fig2_miri_preview.png">Open the current Figure 2 for comparison</a></div><div class="grid">']
md = ['# Further ideas for Figure 2', '', intro, '', '## Recommendation', '', recommendation, '']
for c in cards:
    parts += ['<article><div class="rank">'+escape(c['rank'])+'</div><h2>'+escape(c['title'])+'</h2>',
              '<p><a href="'+c['url']+'">'+escape(c['source'])+'</a></p><div class="images">']
    for img in c['images']:
        parts.append('<a href="'+img+'"><img loading="lazy" src="'+img+'" alt="Original paper page: '+escape(c['source'])+'"></a>')
    parts.append('</div>')
    md.extend(['## '+c['title'], '', '['+c['source']+']('+c['url']+')', ''])
    for key, label in [('shows','What the paper shows'),('borrow','Possible adaptation'),('limit','Limit for our proposal')]:
        parts.append('<p><b>'+label+':</b> '+escape(c[key])+'</p>')
        md.extend(['**'+label+':** '+c[key], ''])
    parts.append('</article>')
parts.append('</div><div class="summary"><b>Additional observational example</b><p>')
parts.append('<a href="armus_2209.13125_p6.png">Armus et al., Figure 2</a> compares total and nuclear MRS spectra of NGC 7469. '
             'The contrast in PAH features is especially readable. The spectra are arbitrarily scaled, so their vertical '
             'separation must not be interpreted as a measured host fraction. '
             '<a href="https://arxiv.org/abs/2209.13125">Paper</a>.</p></div>')
parts.append('<small>Literature comparison prepared 29 September 2026. Figures are original source-page previews '
             'for review, not newly measured data or proposed figure replacements. No proposal PDF was compiled.</small></main></html>')
(HERE/'index.html').write_text('\n'.join(parts))
(HERE/'README.md').write_text('\n'.join(md)+'\n')
after = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in protected}
assert before == after
for c in cards:
    for img in c['images']:
        assert (HERE/img).is_file(), img
assert (HERE/'../../fig2_miri_preview.png').is_file()
(HERE/'checks.json').write_text(json.dumps({'protected_files_unchanged':before==after,'sha256':after,
                                          'source_figures':len(cards),'proposal_pdf_compiled':False},indent=2)+'\n')
print(HERE/'index.html')

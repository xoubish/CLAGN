"""Draw the AGN-history experiment. All curves are explicitly schematic.

No measured flux, fitted delay, target prediction or significance is supplied.
The old physical dust-evolution calculation remains in make_fig2_diagnostics.py
for the exploratory P9694 audit; it is not the new proposal figure.
"""
from pathlib import Path
import csv
import hashlib
import json
import os

os.environ.setdefault('MPLCONFIGDIR', '/tmp/clagn-matplotlib')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Ellipse, FancyArrowPatch
from matplotlib.transforms import Affine2D
import numpy as np

HERE = Path(__file__).resolve().parent
BLUE, RED, GOLD = '#24658a', '#bc5039', '#bb8c35'
PURPLE, GREEN, INK = '#72518b', '#4c826e', '#283542'


def panel(ax, title):
    ax.set_title(title, loc='left', fontsize=10.5, fontweight='bold', pad=10)
    ax.spines[['top', 'right']].set_visible(False)
    ax.spines[['left', 'bottom']].set_color('#9ca4a9')
    ax.tick_params(length=3, color='#9ca4a9', labelsize=9.5)


def agn_cartoon(ax, center, bright_outer, vertical=0):
    """Identical engines and geometry; only the outer dust emissivity differs.

    This is a visual shorthand for delayed response, not an echo surface or
    spatially resolved prediction. Both inner hot components have responded.
    """
    initial_patches = len(ax.patches)
    outer = '#d58a32' if bright_outer else '#9baab6'
    for scale, alpha in [(1.16, .045), (1.08, .07), (1., .15)]:
        ax.add_patch(Ellipse((center, 0), 1.85*scale, 1.02*scale,
                            facecolor=outer, edgecolor='none', alpha=alpha))
    ax.add_patch(Ellipse((center, 0), 1.06, .57, facecolor='white', edgecolor='none'))
    # Fixed clump positions, sizes and counts for both cartoons.
    for radius, count, phase in [(.74, 22, .05), (.92, 28, .17)]:
        theta=np.linspace(0, 2*np.pi, count, endpoint=False)+phase
        for angle in theta:
            ax.add_patch(Circle((center+radius*np.cos(angle), .54*radius*np.sin(angle)),
                                .041, facecolor=outer, edgecolor='white', lw=.25,
                                alpha=.9 if bright_outer else .65, zorder=2))
    # A common hot inner rim and common central disc emphasize equal current L.
    theta=np.linspace(0,2*np.pi,22,endpoint=False)
    for angle in theta:
        ax.add_patch(Circle((center+.52*np.cos(angle), .27*np.sin(angle)),
                            .029, facecolor='#e8b051', edgecolor='none', zorder=3))
    for angle in [.45, 2.7, 4.15]:
        start=(center+.23*np.cos(angle), .12*np.sin(angle))
        end=(center+.46*np.cos(angle), .24*np.sin(angle))
        ax.add_patch(FancyArrowPatch(start,end,arrowstyle='->',mutation_scale=6,
                                    color='#b58a36',lw=.7,zorder=4))
    for width, height, color in [(.49,.23,'#fff0c5'),(.38,.16,'#e9b044'),(.26,.10,'#dc813b')]:
        ax.add_patch(Ellipse((center,0),width,height,facecolor=color,edgecolor='none',zorder=5))
    ax.add_patch(Circle((center,0),.068,facecolor=INK,edgecolor='white',lw=.45,zorder=6))

    for patch in list(ax.patches)[initial_patches:]:
        patch.set_transform(Affine2D().translate(0, vertical) + ax.transData)

def main():
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10,
                         'axes.labelsize': 10, 'axes.labelcolor': INK,
                         'text.color': INK, 'pdf.fonttype': 42, 'ps.fonttype': 42})
    sample_path = HERE/'inputs/jwst_sample_cycle6.csv'
    with sample_path.open() as handle:
        rows = list(csv.DictReader(handle))
    # Read the actual sample rather than silently assuming a redshift limit.
    zkey = next(key for key in rows[0] if key.lower() in ['z', 'redshift'])
    redshift = np.array([float(row[zkey]) for row in rows])
    common = (4.9/(1+redshift.min()), 27.9/(1+redshift.max()))
    assert common[0] < 9.7 < 18 < common[1]

    fig = plt.figure(figsize=(7.6, 2.7), facecolor='white')
    a = fig.add_axes([.015, .225, .185, .62])
    b = fig.add_axes([.295, .225, .195, .62])
    c = fig.add_axes([.570, .225, .172, .62])
    d = fig.add_axes([.822, .225, .172, .62])

    a.set_title('A  AGN history',loc='left',fontsize=10.5,fontweight='bold',pad=10)
    a.set(xlim=(-1.18,1.18),ylim=(-1.65,1.65))
    a.axis('off')
    for vertical, bright, colour, label in [(.76,False,BLUE,'Rising'),(-.88,True,RED,'Fading')]:
        agn_cartoon(a,0,bright,vertical)
        a.text(-1.05,vertical+.66,label,color=colour,fontsize=9.5,fontweight='bold')
        current=vertical+.69
        past=current+(.12 if bright else -.12)
        a.plot([.16,.43,.43,.89],[past,past,current,current],
               color=colour,lw=1.3,ls='--' if bright else '-')
        a.scatter([.89],[current],s=7,color=INK,zorder=7)
    a.text(.5,-.15,'Same current power',ha='center',fontsize=9,
           transform=a.transAxes)

    panel(b, 'B  Continuum')
    wave = np.linspace(2, common[1], 500)
    # Phenomenological unit-gain linear response to a step at t=-1:
    # R(lambda,t=0)=1+(L_past-1)*exp[-1/tau(lambda)].
    # tau increases with lambda by choice, not from a target or grain model.
    tau = .08*(wave/2)**1.6
    rise = 1 + (.35-1)*np.exp(-1/tau)
    fade = 1 + (1.8-1)*np.exp(-1/tau)
    b.axvspan(common[0], common[1], color='#eef0f2', zorder=0)
    b.axhline(1, color='#69747b', ls=':', lw=1.2)
    b.plot(wave, fade, color=RED, lw=2.2, ls='--')
    b.plot(wave, rise, color=BLUE, lw=2.2)
    b.text(6, 1.59, 'Fading', color=RED, fontsize=9.5)
    b.text(6, .49, 'Rising', color=BLUE, fontsize=9.5)
    b.text(common[1]-.5, 1.035, 'Equilibrium at\ncurrent power', ha='right',
           va='bottom', color='#465159', fontsize=9, linespacing=1.08)
    b.set(xlim=(2, common[1]), ylim=(.4, 1.85), xlabel='Rest λ (µm)',
          ylabel='Continuum / equilibrium')
    b.set_xticks([5, 10, 20]); b.set_yticks([1], ['1'])
    panel(c, 'C  Silicates')
    w = np.linspace(common[0], common[1], 600)
    continuum = .76 + .026*(w-common[0])
    features = .52*np.exp(-.5*((w-9.7)/1.05)**2) + .43*np.exp(-.5*((w-18)/1.8)**2)
    spectrum = continuum+features
    c.fill_between(w, continuum, spectrum, color=GOLD, alpha=.25)
    c.plot(w, spectrum, color=INK, lw=1.8)
    c.plot(w, continuum, color='#8b9197', lw=1.2, ls=':')
    for center in [9.7, 18]:
        c.axvline(center, color=GOLD, ls=':', lw=.8, ymax=.81)
        c.text(center, 1.74, f'{center:g}', ha='center', fontsize=9.5, color='#8c6623')
    c.set(xlim=common, ylim=(.55, 1.9), xlabel='Rest λ (µm)',
          ylabel='Flux (arbitrary)')
    c.set_xticks([5, 10, 20]); c.set_yticks([])
    panel(d, 'D  Host light')
    x = np.linspace(-3.5, 3.5, 600)
    nucleus = np.exp(-.5*(x/.43)**2)
    host = .4*np.exp(-.5*(x/1.7)**2)
    d.fill_between(x, 0, host, color=GREEN, alpha=.18)
    d.plot(x, host, color=GREEN, lw=1.6, ls='--')
    d.plot(x, nucleus, color=PURPLE, lw=1.7)
    d.plot(x, nucleus+host, color=INK, lw=1.7)
    d.annotate('Total', xy=(0, 1.4), xytext=(.65, 1.45), fontsize=9.5,
                arrowprops={'arrowstyle': '-', 'color': INK, 'lw': .7})
    d.annotate('Nucleus', xy=(-.35, .72), xytext=(-3.4, 1.10),
                color=PURPLE, fontsize=9.5, arrowprops={'arrowstyle': '-', 'color': PURPLE, 'lw': .7})
    d.text(1.0, .45, 'Host', color=GREEN, fontsize=9.5)
    d.set(xlim=(-3.5, 3.5), ylim=(0, 1.7), xlabel='Position (arb.)',
          ylabel='Brightness (arb.)')
    d.set_xticks([0], ['0']); d.set_yticks([])
    fig.savefig(HERE/'fig2_diagnostics.pdf', dpi=220)
    plt.close(fig)

    np.savetxt(HERE/'inputs/fig2_memory_response.csv', np.c_[wave,tau,rise,fade], delimiter=',',
               header='rest_um,assumed_response_time_arbitrary_units,rising_over_equilibrium,fading_over_equilibrium', comments='')
    provenance = dict(
        figure='fig2_diagnostics.pdf', generator='make_fig2_memory.py',
        scope='Conceptual AGN-memory experiment; not target spectra, measured lags, a dust radiative-transfer calculation or a population-power forecast.',
        presentation='Four panels in one horizontal row with larger labels, wider continuum panel, and an explicit equilibrium-at-current-power label on the unity baseline; cartoons stacked within panel A.',
        actual_targets_fitted=[], sample_z_range=[float(redshift.min()),float(redshift.max())],
        common_mrs_rest_um=list(common),
        panel_a=dict(current_luminosity=1, past_rising=.35, past_fading=1.8,
                     event_time=-1, observation_time=0, time_units='arbitrary',
                     cartoons='Same central-engine brightness, hot inner rim, outer dust geometry and clump count. Outer dust brightness alone differs to depict delayed illumination; not destruction/formation or a calculated equal-arrival-time surface. Miniature histories are qualitative.'),
        panel_b=dict(formula='1 + (L_past - 1) * exp(-1/tau(lambda))',
                     assumed_tau='0.08 * (lambda_rest_um / 2)**1.6 in arbitrary time units',
                     assumptions='Unit-gain exponential linear response; increasing delay with wavelength imposed for illustration. Ratios and amplitudes are not target predictions.',
                     ordinate='Continuum divided by the equilibrium prediction at the same current nuclear luminosity; unity is equilibrium.',
                     miri_shading='Common nominal rest-frame coverage for the actual sample; not a sensitivity boundary.'),
        panel_c='Arbitrary linear continuum plus Gaussian emission features at 9.7 and 18 microns; illustrates joint fitting, not a predicted feature ratio or unique geometry.',
        panel_d='Arbitrary narrow and broad Gaussian spatial components; illustrates point-source/extended decomposition, not resolved torus imaging or a guaranteed host-free spectrum.',
        limitations=['No unique accretion-history recovery.', 'No proof that fixed and evolving dust can be distinguished.',
                     'No validated signal amplitude or sample discrimination forecast.',
                     'Long-wavelength response may average over epochs before the observed light curves.'],
        references=['https://doi.org/10.3847/1538-4357/abee14', 'https://arxiv.org/abs/1909.11101',
                    'https://arxiv.org/abs/0903.2422', 'https://arxiv.org/abs/2504.01103',
                    'https://jwst-docs.stsci.edu/jwst-mid-infrared-instrument/miri-observing-modes/miri-medium-resolution-spectroscopy'],
        generator_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        sample_sha256=hashlib.sha256(sample_path.read_bytes()).hexdigest(),
        previous_figure='../archive/jwst_cleanup_20260926.zip::review/local_backups/fig2_before_agn_memory/',
        previous_spectra='inputs/fig2_diagnostics_spectra.csv (legacy fixed/evolving illustration; not plotted here)')
    (HERE/'inputs/fig2_diagnostics_provenance.json').write_text(json.dumps(provenance, indent=2)+'\n')
    print('Generated schematic AGN-memory Figure 2; common rest coverage:', common)


if __name__ == '__main__':
    main()

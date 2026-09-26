"""Reproduce the provisional sample and documented timing-review decisions.

The original pool came from the saved W1 manifold. The timing review retains
P9694 as a variable science candidate, removes the dedicated comparison arm,
and records seven additions for replicated tests.
This is not a completed selection by event age or dust-response phase.
"""
import csv
import json
from pathlib import Path

import numpy as np
from scipy.ndimage import label
from scipy.spatial import cKDTree

HERE = Path(__file__).resolve().parent
INPUTS = HERE / 'inputs'
K_VALUES = (25, 50, 75)
FLOORS = (.16, .20, .16)
MAX_DISTANCE = .75


def read_csv(path):
    with path.open(newline='') as handle:
        return list(csv.DictReader(handle))


def number(row, key):
    try:
        return float(row[key])
    except (ValueError, KeyError):
        return np.nan


def neighborhood_scores(tree, known, xy):
    distance, indices = tree.query(xy, k=max(K_VALUES))
    scores = np.column_stack([known[indices[:, :k]].mean(axis=1) for k in K_VALUES])
    return scores, distance[:, 0]


def build_selection():
    parent = read_csv(INPUTS / 'sampleA_embedding.csv')
    pool = read_csv(INPUTS / 'selection_pool.csv')
    coords = np.array([[number(r, 'umap_x'), number(r, 'umap_y')] for r in parent])
    known = np.array([r['is_known_clagn'] == 'True' for r in parent])
    tree = cKDTree(coords)
    # A fixed grid defines connected regions. Actual point scores, not grid
    # interpolation, determine whether each candidate passes the thresholds.
    gx, gy = np.meshgrid(np.linspace(0, 24, 301), np.linspace(-5, 12, 214))
    grid_scores, grid_distance = neighborhood_scores(tree, known, np.c_[gx.ravel(), gy.ravel()])
    rich = (np.all(grid_scores >= FLOORS, axis=1) & (grid_distance <= MAX_DISTANCE)).reshape(gx.shape)
    regions, _ = label(rich, structure=np.ones((3, 3)))
    supported_xy = np.c_[gx[rich], gy[rich]]
    region_tree = cKDTree(supported_xy)
    region_ids = regions[rich]
    valid = [r for r in pool if np.isfinite([number(r, 'umap_x'), number(r, 'umap_y')]).all()]
    xy = np.array([[number(r, 'umap_x'), number(r, 'umap_y')] for r in valid])
    scores, distances = neighborhood_scores(tree, known, xy)
    grid_distance, grid_index = region_tree.query(xy)
    for r, score, distance, dg, ig in zip(valid, scores, distances, grid_distance, grid_index):
        r.update({f'f{k}': float(v) for k, v in zip(K_VALUES, score)})
        r['nearest_parent_distance'] = float(distance)
        r['passes_manifold'] = bool(np.all(score >= FLOORS) and distance <= MAX_DISTANCE)
        r['component'] = int(region_ids[ig]) if r['passes_manifold'] and dg < .15 else 0
    ngps_pass = [r for r in valid if r['origin'] == 'NGPS' and r['passes_manifold'] and r['component']]
    component_ids, counts = np.unique([r['component'] for r in ngps_pass], return_counts=True)
    chosen_component = int(component_ids[np.argmax(counts)])
    # The core is the connected enriched component containing the largest
    # number of observed NGPS targets. P1823 is the stated SPHEREx exception
    # to the low-redshift, two-silicate-feature criterion.
    core = [r for r in ngps_pass if r['component'] == chosen_component
            and (0 < number(r, 'z') < .3 or r['name'] == 'P1823')]
    anchors = [r for r in valid if r['origin'] == 'catalogue'
               and r['passes_manifold'] and r['component'] == chosen_component
               and 0 < number(r, 'z') < .3
               and number(r, 'w3snr') >= 3 and number(r, 'w4snr') >= 3
               and r['wise_ccf'] == '0000' and number(r, 'wise_ext') == 0]
    # Retain all qualifying anchors in the NGPS component, rather than adding
    # an arbitrary quota from disconnected CLAGN-rich islands.
    anchors.sort(key=lambda r: (-r['f50'], r['name']))
    former_comparisons = [r for r in valid if r['origin'] == 'NGPS'
                   and 0 < number(r, 'z') < .3
                   and max(r[f'f{k}'] for k in K_VALUES) <= .08
                   and number(r, 'w3snr') >= 3 and number(r, 'w4snr') >= 3]
    # P9694 has a dated optical rise and subsequent IR brightening; it is not
    # a stable control. P9584 has no demonstrated matched-control role.
    additional_variables = [r for r in former_comparisons if r['name'] == 'P9694']
    expansion = json.loads((INPUTS / 'sample_expansion.json').read_text())
    by_name = {r['name']: r for r in valid}
    expansion_groups = []
    for addition in expansion['added']:
        row = by_name[addition['name']]
        assert 0 < number(row, 'z') < .3
        assert number(row, 'w3snr') >= 3 and number(row, 'w4snr') >= 3
        if row['origin'] == 'catalogue':
            assert row['wise_ccf'] == '0000' and number(row, 'wise_ext') == 0
        expansion_groups.append((addition['role'], [row]))
    selected = []
    for role, members in [('variable candidate', core), ('literature anchor', anchors),
                          ('variable candidate', additional_variables)] + expansion_groups:
        for r in members:
            r = dict(r)
            r.update(id=f'S{len(selected)+1:02d}', role=role,
                     cl_neighbor_fraction=r['f50'], si18_observed_um=18*(1+number(r, 'z')))
            selected.append(r)
    summary = dict(
        selection='Expanded 18-source working sample for replicated history comparisons; NGPS priority with literature anchors, no dedicated controls. Event-age/lag and class-allocation audit pending.',
        parent_rows=len(parent), labelled_clagn=int(known.sum()), parent_labelled_fraction=float(known.mean()),
        candidate_pool=len(pool), ngps_pool=sum(r['origin'] == 'NGPS' for r in pool),
        rule=dict(k_neighbors=list(K_VALUES), minimum_labelled_fractions=list(FLOORS),
                  maximum_nearest_parent_distance=MAX_DISTANCE,
                  connected_region='Component containing the most NGPS targets passing the manifold rule',
                  science_redshift='0 < z < 0.3; P1823 retained for its optical/SPHEREx history',
                  anchors='All catalogue CLAGN in the selected component with z<0.3, W3/W4 SNR>=3, AllWISE ccf=0000 and ext_flag=0',
                  historical_comparisons='NGPS targets with z<0.3, fractions<=0.08 at all three neighborhood sizes, W3/W4 SNR>=3',
                  timing_review='Retain P9694 for its 2021–2022 optical rise and brighter state through 2025; omit P9584. See timing_audit.json.',
                  expansion='Seven documented additions from sample_expansion.json. Their measured variability and archival quality determine inclusion, not manifold membership.'),
        chosen_component=chosen_component,
        ngps_core=[r['name'] for r in core], literature_anchors=[r['name'] for r in anchors],
        additional_variables=[r['name'] for r in additional_variables],
        expansion=[r['name'] for r in expansion['added']],
        all_observed_candidates=[r['name'] for r in selected if r['origin'] == 'NGPS'],
        all_literature_anchors=[r['name'] for r in selected if r['origin'] == 'catalogue'],
        comparisons=[], removed_from_request=['P9584'], total=len(selected),
        response_phase_coverage_validated=False,
        interpretation='Neighborhood label fractions record pool provenance, not response phases. Eighteen is a working design with a goal of roughly six sources per history family, not a demonstrated balanced allocation or validated power calculation.',
        ngps_audit=[dict(name=r['name'], f25=r['f25'], f50=r['f50'], f75=r['f75'],
                         component=r['component'], included=r['name'] in {s['name'] for s in selected},
                         spectral_history_note=r['spectral_history_note']) for r in valid if r['origin'] == 'NGPS'],
        grid=dict(x_min=0, x_max=24, nx=301, y_min=-5, y_max=12, ny=214),
    )
    return selected, summary, (gx, gy, rich, regions == chosen_component)


def main():
    selected, summary, _ = build_selection()
    with (INPUTS / 'jwst_sample.csv').open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(selected[0]))
        writer.writeheader()
        writer.writerows(selected)
    (INPUTS / 'selection_summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    print(f"Selection: {len(summary['all_observed_candidates'])} observed variable candidates + "
          f"{len(summary['all_literature_anchors'])} literature anchors = {len(selected)} targets; no dedicated controls")


if __name__ == '__main__':
    main()

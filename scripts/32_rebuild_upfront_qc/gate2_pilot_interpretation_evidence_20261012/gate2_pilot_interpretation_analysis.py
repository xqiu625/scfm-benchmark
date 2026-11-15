"""READ-ONLY recomputation of candidate-set structure for the Gate 2 GSM5102900_HC1
pilot (job 29882644), run directly against the real sidecar files on HPCC.

Does not re-run Scrublet, does not read outcome fields, does not touch the four
original output files. Reads:
  SCORING_MANIFEST.json   (per-prior threshold, histograms, n_called as recorded)
  PER_PRIOR_CALLS.json    (raw boolean call vectors)
  PER_PRIOR_SCORES.json   (raw float score vectors)

Computes independently, not copied from the manifest or any prior report:
  - n_called per prior, recomputed from the raw boolean vectors
  - the full 28-pair Jaccard matrix over the 8 candidate sets
  - whether every pair of candidate sets is nested (one subset of the other)
  - the number of DISTINCT candidate sets (by value, not by eyeballing a table --
    a prior round miscounted this as 3 when it is 4)
  - the global intersection and union sizes across all 8 priors
  - whether the full per-cell score ranking (argsort) is identical across all 8 priors
  - threshold/mechanical consistency: call == (score > threshold) for every cell, every prior
  - whether each row's observed and simulated histograms share identical bin edges
    (verified, not assumed -- a prior plotting script assumed this without checking)

Writes, both to --out-dir (refuses to overwrite an existing non-empty directory,
same convention as capture_preflight_GSM5102900_HC1.py):
  gate2_pilot_interpretation_recomputed.json  -- this script's own analysis
  gate2_pilot_histogram_snapshot.json         -- an immutable, self-contained copy of
    only the per-prior threshold/histogram fields SCORING_MANIFEST.json carries, bound
    to that manifest's own identity (sha256), for gate2_pilot_interpretation_plots.py
    to read instead of a live, uncommitted HPCC path.
"""
import argparse
import hashlib
import json
import os

import numpy as np


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--out-dir', required=True)
    ap.add_argument('--manifest-dir', required=True,
                     help='directory containing SCORING_MANIFEST.json, PER_PRIOR_CALLS.json, '
                          'PER_PRIOR_SCORES.json for the run to analyse')
    args = ap.parse_args()
    if os.path.exists(args.out_dir) and os.listdir(args.out_dir):
        raise SystemExit(f'REFUSING: {args.out_dir} exists and is non-empty -- '
                          'use a new output directory so no prior result is overwritten')
    os.makedirs(args.out_dir, exist_ok=True)

    manifest_path = os.path.join(args.manifest_dir, 'SCORING_MANIFEST.json')
    with open(manifest_path, 'rb') as fh:
        manifest_bytes = fh.read()
    manifest_sha256 = hashlib.sha256(manifest_bytes).hexdigest()
    manifest = json.loads(manifest_bytes)
    calls = json.load(open(os.path.join(args.manifest_dir, 'PER_PRIOR_CALLS.json')))
    scores = json.load(open(os.path.join(args.manifest_dir, 'PER_PRIOR_SCORES.json')))

    assert calls['cell_order_sha256'] == scores['cell_order_sha256'] == manifest['cell_order_sha256']
    priors = manifest['prior_grid']
    prior_keys = [str(p) for p in priors]
    n_cells = manifest['n_cells']

    row_by_rate = {r['expected_doublet_rate']: r for r in manifest['rows']}

    call_arrays = {}
    score_arrays = {}
    thresholds = {}
    histogram_snapshot_rows = {}
    for k, p in zip(prior_keys, priors):
        c = np.asarray(calls['by_prior'][k], dtype=bool)
        s = np.asarray(scores['by_prior'][k], dtype=float)
        assert c.shape[0] == n_cells and s.shape[0] == n_cells
        call_arrays[k] = c
        score_arrays[k] = s
        row = row_by_rate[p]
        thresholds[k] = row['threshold']
        # mechanical consistency: the recorded call must equal score > threshold
        # for every cell, not just in aggregate count.
        mismatch = int(np.sum(c != (s > thresholds[k])))
        if mismatch:
            raise SystemExit(f'prior {k}: {mismatch} cells disagree with call==(score>threshold)')

        obs_edges = row['observed_score_histogram']['bin_edges']
        sim_edges = row['simulated_score_histogram']['bin_edges']
        edges_shared = obs_edges == sim_edges
        histogram_snapshot_rows[k] = {
            'expected_doublet_rate': p,
            'threshold': thresholds[k],
            'n_called': row['n_called'],
            'n_cells': row['n_cells'],
            'observed_simulated_bin_edges_identical': edges_shared,
            'bin_edges': obs_edges if edges_shared else None,
            'observed_bin_edges': None if edges_shared else obs_edges,
            'simulated_bin_edges': None if edges_shared else sim_edges,
            'observed_counts': row['observed_score_histogram']['counts'],
            'simulated_counts': row['simulated_score_histogram']['counts'],
        }

    n_called_recomputed = {k: int(call_arrays[k].sum()) for k in prior_keys}
    n_called_manifest = {k: manifest['n_called'][k] for k in prior_keys}
    assert n_called_recomputed == n_called_manifest, (n_called_recomputed, n_called_manifest)

    # Distinct candidate sets, counted by value (frozenset), not by table-reading --
    # the direct cause of a prior "3 distinct, not 8" miscount (the correct count is 4).
    sets = {k: frozenset(np.flatnonzero(call_arrays[k]).tolist()) for k in prior_keys}
    distinct_sets = sorted(set(sets.values()), key=len)
    distinct_groups = [
        {'size': len(s), 'priors': sorted({k for k, v in sets.items() if v == s}, key=float)}
        for s in distinct_sets
    ]

    jaccard = {}
    nested = {}
    for i, a in enumerate(prior_keys):
        for b in prior_keys[i + 1:]:
            A, B = sets[a], sets[b]
            inter = len(A & B)
            union = len(A | B)
            j = inter / union if union else 1.0
            jaccard[f'{a}|{b}'] = j
            nested[f'{a}|{b}'] = (A <= B) or (B <= A)

    global_intersection = set.intersection(*(set(v) for v in sets.values()))
    global_union = set.union(*(set(v) for v in sets.values()))

    # argsort-identity: do all 8 score vectors rank the SAME cells in the same
    # order (stable sort, ties broken by cell index so the comparison is exact)?
    # This is a fact about rank order ONLY. It is not evidence that the score
    # VALUES share a common scale across priors -- expected_doublet_rate is an
    # input to scanpy's per-cell score formula (calculate_doublet_scores passes
    # exp_doub_rate into _nearest_neighbor_classifier), so threshold values from
    # different priors are not interchangeable cut points on one fixed axis.
    base_key = prior_keys[0]
    base_rank = np.lexsort((np.arange(n_cells), score_arrays[base_key]))
    argsort_identical_to_base = {}
    for k in prior_keys[1:]:
        rank = np.lexsort((np.arange(n_cells), score_arrays[k]))
        argsort_identical_to_base[k] = bool(np.array_equal(base_rank, rank))

    moved_out = sorted(sets['0.16'] - sets['0.2'])
    moved_in = sorted(sets['0.2'] - sets['0.16'])

    threshold_values = [thresholds[k] for k in prior_keys]
    threshold_monotonic_nondecreasing = all(
        threshold_values[i] <= threshold_values[i + 1] for i in range(len(threshold_values) - 1))

    out = {
        'source': {
            'manifest': manifest_path,
            'manifest_sha256': manifest_sha256,
            'calls_sidecar': os.path.join(args.manifest_dir, 'PER_PRIOR_CALLS.json'),
            'scores_sidecar': os.path.join(args.manifest_dir, 'PER_PRIOR_SCORES.json'),
            'cell_order_sha256': manifest['cell_order_sha256'],
            'n_cells': n_cells,
        },
        'per_prior': {
            k: {
                'threshold': thresholds[k],
                'n_called_recomputed': n_called_recomputed[k],
                'fraction_of_cells': n_called_recomputed[k] / n_cells,
            }
            for k in prior_keys
        },
        'threshold_monotonic_nondecreasing_across_prior_grid_order': threshold_monotonic_nondecreasing,
        'n_distinct_candidate_sets': len(distinct_sets),
        'distinct_candidate_set_groups': distinct_groups,
        'jaccard_pairwise_recomputed': jaccard,
        'nested_pairwise': nested,
        'all_pairs_nested': all(nested.values()),
        'global_intersection_size': len(global_intersection),
        'global_union_size': len(global_union),
        'global_intersection_cells_zero_indexed': sorted(global_intersection),
        'argsort_identical_to_prior_0.005': argsort_identical_to_base,
        'all_argsorts_identical': all(argsort_identical_to_base.values()),
        'argsort_identity_is_about_rank_order_only':
            'expected_doublet_rate feeds the per-cell score formula itself '
            '(scanpy calculate_doublet_scores -> _nearest_neighbor_classifier); '
            'argsort identity does not imply the score VALUES share one scale across priors',
        'prior_0.16_to_0.2_moved_out_count': len(moved_out),
        'prior_0.16_to_0.2_moved_in_count': len(moved_in),
        'mechanical_consistency_checked': 'call == (score > threshold) verified for every cell, every prior',
    }
    with open(os.path.join(args.out_dir, 'gate2_pilot_interpretation_recomputed.json'), 'w') as fh:
        json.dump(out, fh, indent=2, sort_keys=True)

    snapshot = {
        'schema': 'gate2-pilot-histogram-snapshot/1',
        'purpose': 'self-contained input for gate2_pilot_interpretation_plots.py; '
                   'extracted verbatim from the source manifest, bound to it by sha256',
        'source_manifest_path': manifest_path,
        'source_manifest_sha256': manifest_sha256,
        'capture': manifest['capture'],
        'cell_order_sha256': manifest['cell_order_sha256'],
        'n_cells': n_cells,
        'by_prior': histogram_snapshot_rows,
    }
    with open(os.path.join(args.out_dir, 'gate2_pilot_histogram_snapshot.json'), 'w') as fh:
        json.dump(snapshot, fh, indent=2, sort_keys=True)

    print(json.dumps({k: v for k, v in out.items() if k not in (
        'jaccard_pairwise_recomputed', 'nested_pairwise', 'global_intersection_cells_zero_indexed')},
        indent=2, sort_keys=True))
    print('jaccard_min_recomputed:', min(jaccard.values()))
    print('jaccard_max_recomputed:', max(jaccard.values()))


if __name__ == '__main__':
    main()

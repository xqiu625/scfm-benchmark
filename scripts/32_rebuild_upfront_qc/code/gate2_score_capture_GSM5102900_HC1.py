"""Single-capture Gate 2 scoring for GSM5102900_HC1. PERMIT-GATED.

    python gate2_score_capture_GSM5102900_HC1.py --out-dir <derived dir>

Runs the eight approved priors over ONE capture by calling
`gate2_stress_runner.score_capture_across_grid`, which already owns the
preprocessing, the single simulated population, the per-row record and the
frozen failure statuses. Nothing here reimplements that path.

Refuses unless a permit that is a tracked repository artifact matches this
config, these inputs and this implementation. The global
`PRODUCTION_PARAMETER_RULING` is never written; the reviewed ruling is handed to
the runner as a scoped argument only after the permit gate succeeds.

No outcome field is read or written. No cell is deleted. The output is scores,
calls and diagnostics -- a mechanical doublet candidate is not a deletion
authorisation, and no final-authorisation artifact is produced.
"""
import argparse
import gzip
import hashlib
import json
import os
import sys
import time

import gate2_scoring_config_GSM5102900_HC1 as CFG
from gate2_stress_design import RESIDUAL_RATE_GRID
from gate2_stress_runner import (Gate2RunnerError, assert_dependency_approved,
                                 dependency_identity, grid_stability)

ENTRYPOINT_VERSION = 'gate2-score-capture/GSM5102900_HC1/1'


def _write_atomic(path, payload):
    tmp = path + '.inprogress'
    with open(tmp, 'w') as fh:
        json.dump(payload, fh, indent=2, sort_keys=True, default=str)
        fh.flush()
        os.fsync(fh.fileno())
    os.replace(tmp, path)


def _load_deposit():
    """Only the three deposit files, minimised to matrix + cell identity."""
    import anndata as ad
    import numpy as np
    import scipy.io
    import scipy.sparse as sp

    CFG.assert_inputs_unchanged()
    d = CFG.DEPOSIT_DIR
    with gzip.open(os.path.join(d, f'{CFG.CAPTURE}_matrix.mtx.gz'), 'rb') as fh:
        M = scipy.io.mmread(fh).tocsr()
    with gzip.open(os.path.join(d, f'{CFG.CAPTURE}_barcodes.tsv.gz'), 'rt') as fh:
        bcs = [x.strip() for x in fh if x.strip()]
    with gzip.open(os.path.join(d, f'{CFG.CAPTURE}_features.tsv.gz'), 'rt') as fh:
        feats = [x.rstrip('\n').split('\t') for x in fh if x.strip()]
    types = {f[2] for f in feats if len(f) > 2}
    if types != {'Gene Expression'}:
        raise Gate2RunnerError(f'feature types are {types}, expected Gene Expression only')
    a = ad.AnnData(X=sp.csr_matrix(M.T))
    a.obs_names = bcs                       # cell identity only; obs carries nothing else
    a.var_names = [f[0] for f in feats]
    a.var_names_make_unique()
    if list(a.obs.columns):
        raise Gate2RunnerError(f'obs must be empty, holds {list(a.obs.columns)}')
    del np
    return a


class _GraphCapture:
    """Capture the ACTUAL neighbour graph scanpy consumes, without changing it.

    Two earlier versions of this were wrong, and both wrongs passed:

    1. `_graph_evidence()` recorded only a cell count and a row count. Neither is
       a graph fact, and naming the function that way overstated it.
    2. The replacement patched `sklearn.neighbors.NearestNeighbors.kneighbors`.
       Measured on the pinned scanpy: that captures ZERO constructions. scanpy
       reaches a `KNeighborsTransformer`, not `NearestNeighbors`, and
       `NearestNeighbors.kneighbors` IS the inherited `KNeighborsMixin.kneighbors`,
       so assigning to the subclass never shadows what the transformer resolves.
       pynndescent is also installed, so the backend is not guaranteed to be
       sklearn at all.

    The capture point is therefore the value the doublet score CONSUMES --
    `_get_indices_distances_from_sparse_matrix(knn.distances, k_adj)` in
    scanpy/preprocessing/_scrublet/core.py -- which is backend-independent. Values
    are passed through untouched.
    """

    def __init__(self):
        self.captures = []
        self._orig = None

    def __enter__(self):
        from scanpy.preprocessing._scrublet import core as _core
        self._core = _core
        self._orig = _core._get_indices_distances_from_sparse_matrix
        captures = self.captures
        orig = self._orig

        def patched(D, n_neighbors):
            idx, dist = orig(D, n_neighbors)
            captures.append({'idx': idx, 'dist': dist,
                             'n_neighbors_arg': int(n_neighbors),
                             'n_queries': int(idx.shape[0]),
                             'columns': int(idx.shape[1])})
            return idx, dist

        _core._get_indices_distances_from_sparse_matrix = patched
        return self

    def __exit__(self, *exc):
        self._core._get_indices_distances_from_sparse_matrix = self._orig
        return False

    def digest(self):
        """One record per captured graph, in call order."""
        import hashlib as _h
        import numpy as _np
        out = []
        for c in self.captures:
            idx_sha = _h.sha256(_np.ascontiguousarray(c['idx']).tobytes()).hexdigest()
            dist_sha = (_h.sha256(_np.ascontiguousarray(
                _np.round(c['dist'], 9)).tobytes()).hexdigest()
                if c['dist'] is not None else None)
            out.append({'shape': list(c['idx'].shape),
                        'n_neighbors_arg': c['n_neighbors_arg'],
                        'n_queries': c['n_queries'], 'columns': c['columns'],
                        'indices_sha256': idx_sha, 'distances_sha256': dist_sha,
                        # Recorded, not asserted: measured False on this path,
                        # contrary to the scanpy helper's own docstring.
                        'self_first_column': bool(
                            (c['idx'][:, 0] == _np.arange(c['idx'].shape[0])).all())})
        return out


def _graph_evidence(capture, n_priors, k_adj_expected, expected_shape):
    """Real graph facts, the cross-prior identity check, and the ACCEPTANCE list.

    Returns (evidence, problems). `problems` is non-empty whenever the captured
    graph fails any reviewed expectation; the caller must stop on it.

    The previous version stopped only when `identical_across_priors` was False.
    With zero captures that field is None and with one capture it is True, so a
    run that built no graph at all -- which is exactly what the old sklearn-level
    capture produced -- finalized `complete`. `k_matches_expected` was computed and
    never read.
    """
    per_call = capture.digest()
    idx_shas = {g['indices_sha256'] for g in per_call}
    dist_shas = {g['distances_sha256'] for g in per_call}
    expected_n = n_priors * CFG.GRAPH_CONSTRUCTIONS_PER_PRIOR

    problems = []
    if len(per_call) != expected_n:
        problems.append(
            f'graph constructions: observed {len(per_call)}, reviewed expectation '
            f'{expected_n} ({CFG.GRAPH_CONSTRUCTIONS_PER_PRIOR} per prior x '
            f'{n_priors} priors)')
    for i, g in enumerate(per_call):
        if g['n_neighbors_arg'] != k_adj_expected:
            problems.append(f'construction {i}: n_neighbors argument '
                            f'{g["n_neighbors_arg"]} != k_adj {k_adj_expected}')
        if tuple(g['shape']) != tuple(expected_shape):
            problems.append(f'construction {i}: shape {tuple(g["shape"])} != '
                            f'expected {tuple(expected_shape)}')
        if not g['indices_sha256']:
            problems.append(f'construction {i}: empty indices checksum')
        if not g['distances_sha256']:
            problems.append(f'construction {i}: empty distances checksum')
    if per_call and len(idx_shas) != 1:
        problems.append(f'the neighbour graph topology differed across priors: '
                        f'{len(idx_shas)} distinct index checksums')

    evidence = {
        'source': f'captured from {CFG.GRAPH_CAPTURE_POINT}; this is the graph the '
                  'doublet score consumes, independent of the neighbour backend. '
                  'Values are passed through untouched.',
        'preflight_verified': False,
        'preflight_note': 'the score-free preflight could not construct this graph; '
                          'this is a post-run pilot-acceptance observation',
        'n_graph_constructions': len(per_call),
        'expected_graph_constructions': expected_n,
        'constructions_per_prior_expected': CFG.GRAPH_CONSTRUCTIONS_PER_PRIOR,
        'per_construction': per_call,
        'k_adj_expected': k_adj_expected,
        'expected_shape': list(expected_shape),
        'columns_rule': CFG.GRAPH_COLUMNS_RULE,
        'topology_identical_across_priors': len(idx_shas) == 1 if per_call else False,
        # ENFORCED: topology (indices). RECORDED ONLY: distances. Distance
        # invariance across priors is not claimed, so it is not made a gate; this
        # field states what was observed, and nothing branches on it.
        'distances_identical_across_priors_observed':
            (len(dist_shas) == 1 if per_call else False),
        'distances_enforced': False,
        'acceptance_problems': problems,
    }
    return evidence, problems


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--out-dir', required=True)
    ap.add_argument('--permit', default=None)
    args = ap.parse_args()

    # ---- gates, before anything is loaded --------------------------------
    permit = CFG.assert_permitted(args.permit)
    ruling = CFG.scoped_ruling(permit)
    cfg = CFG.scoring_config()
    expected_dir = cfg['output']['directory']
    if os.path.abspath(args.out_dir) != os.path.abspath(expected_dir):
        raise Gate2RunnerError(
            f'--out-dir {args.out_dir} is not the derived directory {expected_dir}')
    if os.path.exists(args.out_dir) and os.listdir(args.out_dir):
        raise Gate2RunnerError(f'output reuse is forbidden: {args.out_dir} is not empty')
    os.makedirs(args.out_dir, exist_ok=True)
    dest = os.path.join(args.out_dir, 'SCORING_MANIFEST.json')
    t0 = time.time()

    import numpy as np
    import scanpy as sc
    import anndata as ad
    import scipy
    import sklearn

    identity = assert_dependency_approved(dependency_identity(sc))
    versions = {'scanpy': sc.__version__, 'anndata': ad.__version__,
                'numpy': np.__version__, 'scipy': scipy.__version__,
                'sklearn': sklearn.__version__, 'python': sys.version.split()[0]}
    drift = {k: (v, CFG.DEPENDENCY_VERSIONS[k]) for k, v in versions.items()
             if CFG.DEPENDENCY_VERSIONS.get(k) != v}
    if drift:
        raise Gate2RunnerError(f'dependency drift from the audited versions: {drift}')

    base = {
        'entrypoint': ENTRYPOINT_VERSION,
        'capture': CFG.CAPTURE, 'cohort': CFG.COHORT, 'channel_id': CFG.CHANNEL,
        'status': 'incomplete', 'stage': 'gated',
        'permit_id': permit.get('permit_id'),
        'config_sha256': CFG.config_sha256(),
        'implementation_sha256': CFG.implementation_sha256(),
        'preflight_commit': CFG.PREFLIGHT_COMMIT,
        'dependency_identity': identity, 'versions': versions,
        'production_parameter_ruling_written': False,
        'outcome_columns_read': [],
        'deletion_authorised': False,
        'started_utc': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
    }
    _write_atomic(dest, base)

    a = _load_deposit()
    base.update(stage='loaded', input_shape=list(a.shape))
    _write_atomic(dest, base)

    with _GraphCapture() as gcap:
        rows, sidecars, score_sidecars = score_capture(a, sc, ruling)
    graph, graph_problems = _graph_evidence(
        gcap, len(CFG.RESIDUAL_RATE_GRID), CFG.RESOLVED['k_adj'],
        CFG.expected_graph_shape())

    # ---- fail closed on drift the runner records -------------------------
    problems = []
    for r in rows:
        if r.get('identity_verified') is False:
            problems.append(f'cell identity mismatch at rate {r["expected_doublet_rate"]}')
        if r['n_cells'] is not None and r['n_cells'] != CFG.RESOLVED['n_cells_effective']:
            problems.append(f'dimension change: n_cells {r["n_cells"]} != '
                            f'{CFG.RESOLVED["n_cells_effective"]}')
        echo = r.get('parameters_echo_agree')
        if echo is not None and echo.get('agree') is False:
            problems.append(f'parameter echo mismatch at rate {r["expected_doublet_rate"]}')
    sims = {r['simulated_checksum'] for r in rows}
    if len(sims) != 1:
        problems.append(f'the simulated population differed across priors: {sims}')
    if problems:
        base.update(status='stopped', stage='drift_detected', problems=problems)
        _write_atomic(dest, base)
        raise SystemExit('; '.join(problems))

    # Zero, missing, extra, wrong-k, wrong-shape, empty-checksum or inconsistent
    # captures all stop the run. Truthiness on the list, never on a tri-state field
    # whose None means "nothing was captured".
    if graph_problems:
        base.update(status='stopped', stage='graph_unacceptable', graph=graph,
                    problems=graph_problems)
        _write_atomic(dest, base)
        raise SystemExit('graph acceptance failed: ' + '; '.join(graph_problems))

    stability = grid_stability(rows, sidecars)
    cell_order = [str(x) for x in a.obs_names]
    cell_order_sha = hashlib.sha256('\u0000'.join(cell_order).encode()).hexdigest()

    # Sidecars are written, validated and promoted BEFORE the manifest is
    # finalized. A previous version wrote status=complete first, so a failure on
    # the next write left a manifest asserting a completion that had not happened.
    sidecar_payloads = {
        'PER_PRIOR_CALLS.json': {
            'schema_version': 'gate2-per-prior-calls/1',
            'cell_order_sha256': cell_order_sha, 'n_cells': len(cell_order),
            'note': 'booleans in cell order; present only for examined priors',
            'by_prior': {str(k): [bool(x) for x in v] for k, v in sidecars.items()}},
        'PER_PRIOR_SCORES.json': {
            'schema_version': 'gate2-per-prior-scores/1',
            'cell_order_sha256': cell_order_sha, 'n_cells': len(cell_order),
            'note': 'observed doublet scores per cell, retained even when the '
                    'threshold failed; never inferred from calls',
            'by_prior': {str(k): [float(x) for x in v]
                         for k, v in score_sidecars.items()}},
        'CELL_ORDER.json': {
            'schema_version': 'gate2-cell-order/1',
            'cell_order_sha256': cell_order_sha, 'cells': cell_order},
    }
    outputs = _promote_sidecars(args.out_dir, sidecar_payloads, len(cell_order))

    base.update(
        status='complete', stage='scored',
        schema_version='gate2-scoring-manifest/1',
        n_rows=len(rows), prior_grid=list(RESIDUAL_RATE_GRID),
        prior_label=CFG.PRIOR_LABEL,
        rows=rows, grid_stability=stability, graph=graph,
        cell_order_sha256=cell_order_sha, n_cells=len(cell_order),
        threshold_status={str(r['expected_doublet_rate']): r['mechanical_status']
                          for r in rows},
        n_called={str(r['expected_doublet_rate']): r['n_called'] for r in rows},
        outputs=outputs,
        elapsed_seconds=round(time.time() - t0, 1))
    _write_atomic(dest, base)
    print(f'rows={len(rows)} examined='
          f'{sum(1 for r in rows if r["mechanical_status"] == "examined")} '
          f'status=complete')


def _promote_sidecars(out_dir, payloads, n_cells):
    """Write every sidecar to a temp file, validate, then promote all of them.

    What this guarantees: no sidecar appears under its final name until EVERY
    payload has been written, fsynced and re-read byte-for-byte. A validation
    failure therefore publishes nothing.

    What it does NOT guarantee: the promotion loop is a sequence of `os.replace`
    calls, not a directory-level atomic operation. A process kill between two
    replaces can leave a partial set on disk. The claim that "a failure part-way
    leaves no half-published output set" was false and has been withdrawn.

    What makes a partial set detectable is the manifest, not this function: the
    manifest is finalized `complete` only after this returns, so a partial set is
    always accompanied by a manifest that does not say `complete`. The manifest is
    the commit marker; these files are not self-certifying.
    """
    temps, outputs = {}, {}
    try:
        for name, payload in payloads.items():
            for prior, vec in payload.get('by_prior', {}).items():
                if len(vec) != n_cells:
                    raise Gate2RunnerError(
                        f'{name}: prior {prior} has {len(vec)} values for {n_cells} cells')
            tmp = os.path.join(out_dir, name + '.inprogress')
            # No `default=str`: stringifying an unexpected object would serialize
            # something other than what was computed, and pass.
            blob = json.dumps(payload, indent=2, sort_keys=True,
                              allow_nan=False).encode()
            with open(tmp, 'wb') as fh:
                fh.write(blob)
                fh.flush()
                os.fsync(fh.fileno())
            with open(tmp, 'rb') as fh:
                if fh.read() != blob:
                    raise Gate2RunnerError(f'{name}: temp file does not match its payload')
            temps[name] = tmp
            outputs[name] = {'sha256': hashlib.sha256(blob).hexdigest(),
                             'bytes': len(blob),
                             'schema_version': payload['schema_version'],
                             'n_priors': len(payload.get('by_prior', {}))}
        for name, tmp in temps.items():
            os.replace(tmp, os.path.join(out_dir, name))
    finally:
        for tmp in temps.values():
            if os.path.exists(tmp):
                os.remove(tmp)
    return outputs


def score_capture(adata, sc, ruling):
    """The production path, not a copy of it."""
    from gate2_stress_runner import score_capture_across_grid
    return score_capture_across_grid(
        CFG.COHORT, CFG.CHANNEL, adata,
        scrublet_fn=sc.pp.scrublet, sc=sc, ruling=ruling)


if __name__ == '__main__':
    main()

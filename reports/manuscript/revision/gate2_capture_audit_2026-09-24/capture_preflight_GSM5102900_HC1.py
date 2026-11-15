"""Single-capture, score-free Gate 2 preflight for GSM5102900_HC1.

    python capture_preflight_GSM5102900_HC1.py --out-dir <NEW versioned dir>

Measures feasibility on the complete 7,666-cell deposit matrix. It produces NO
doublet score, call or rate, reads NO outcome, and leaves
PRODUCTION_PARAMETER_RULING untouched.

WHY THIS VERSION REUSES INSTEAD OF REIMPLEMENTING
-------------------------------------------------
A previous draft hand-rolled the numerical path and validated something else:
it normalised the observed matrix but not the simulated one, fit a JOINT
`np.linalg.svd` over observed+simulated, and built a graph with a base `k` via
sklearn. Scanpy 1.11.3 does none of those. Measured from the pinned source
(`scanpy/preprocessing/_scrublet/pipeline.py`):

  * `mean_center`      : gene means from `_counts_obs_norm` (OBSERVED ONLY),
                         subtracted from observed AND simulated;
  * `normalize_variance`: gene variances from OBSERVED ONLY, applied to both;
  * `truncated_svd`    : `TruncatedSVD(...).fit(_counts_obs_norm)` then
                         `.transform()` each -- an OBSERVED-FIT projection.

So this version drives `scanpy.preprocessing._scrublet.core.Scrublet` itself and
stops at `set_manifold`, the last step before `calculate_doublet_scores`. The
eleven parameters are IMPORTED from `gate2_parameter_ruling`, and the seed comes
from the repository's own `capture_seed`, not invented here.

WHERE THE SCORE-FREE BOUNDARY REALLY IS
---------------------------------------
The kNN graph is built inside `calculate_doublet_scores` (it is what consumes
`k_adj = round(k * (1 + n_sim/n_obs))`). It therefore CANNOT be checksummed
without scoring. This preflight checksums the manifold -- the observed and
simulated PCA representations the graph is built from -- and records the `k_adj`
scanpy would use. The manifold is the last score-free artifact; claiming a graph
checksum here would be claiming something not computed.

INVARIANCE ACROSS THE EIGHT PRIORS IS STRUCTURAL
------------------------------------------------
`expected_doublet_rate` is consumed only by `calculate_doublet_scores`. Every
step up to the manifold is independent of it. A previous draft "proved" this with
a loop whose `rho` was unused -- which showed only that the script reruns
deterministically. This version builds the manifold ONCE, and a test asserts that
no step up to it takes the prior as input.
"""
import argparse
import hashlib
import json
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.abspath(os.path.join(HERE, '..', '..', '..', '..',
                                    'scripts', '32_rebuild_upfront_qc', 'code'))
sys.path.insert(0, CODE)

from gate2_parameter_ruling import (METHODOLOGY_PARAMETER_RULING,  # noqa: E402
                                    PARAMETER_RULING_REFERENCE,
                                    PARAMETER_RULING_VERSION,
                                    approved_for, assert_ruling_version)
from gate2_stress_design import RESIDUAL_RATE_GRID, capture_seed  # noqa: E402
from gate2_stress_runner import (APPROVED_SCANPY_VERSIONS, REQUIRED_PARAMETER_NAMES,  # noqa: E402
                                 _checksum_matrix, assert_dependency_approved,
                                 dependency_identity)

EXPECT_RULING_VERSION = 'gate2-parameter-ruling/1'
CAPTURE = 'GSM5102900_HC1'
COHORT, CHANNEL = 'qiu', 'GSM5102900'
DEPOSIT = '/bigdata/godziklab/shared/Xinru/302006/01_DATA/raw/GSE167363_qiu/extracted'
CONTAINER = ('/bigdata/godziklab/shared/Xinru/302006/01_DATA/raw/'
             'human_platelet_covid_severity_raw_counts.h5ad')
PRIOR_LABEL = 'assumed_algorithmic_stress_setting'

#: Full SHA-256 of every input, asserted BEFORE any evidence is read. Recorded by
#: verify_capture_evidence.py. A prefix is not an assertion.
EXPECTED_SHA256 = {
    f'{CAPTURE}_barcodes.tsv.gz':
        'e6dbf9f587d9a6cf3ecff31cbf4d53b09c6faf9837c82d5c1c0a26227f377a02',
    f'{CAPTURE}_features.tsv.gz':
        'c7408e82af270ca6308a8aa7ea659795a15ab76f2e3a85bd5db4cd3d53244199',
    f'{CAPTURE}_matrix.mtx.gz':
        '4ea259bfc4de02cc9799df563779fce7c9bd4c889532d3ff1178e2caf545b9f8',
}
#: The container's FULL sha256, computed read-only on 2026-09-24. A previous
#: version asserted only `f'{size}:'` with a comment promising the digest later --
#: which is path+size again, the very thing it claimed to replace.
EXPECTED_CONTAINER_BYTES = 5027039977
EXPECTED_CONTAINER_SHA256 = (
    '1a0f37902cd33f410a71428ea45b3d1573fddec393ce27f07a93074d0ed536ba')
EXPECTED_BYTES = {f'{CAPTURE}_barcodes.tsv.gz': 39604,
                  f'{CAPTURE}_features.tsv.gz': 304728,
                  f'{CAPTURE}_matrix.mtx.gz': 38785404}


def repo_provenance():
    """Commit and per-file hashes INSIDE the manifest, not only on stdout."""
    import subprocess
    repo = os.path.abspath(os.path.join(CODE, '..', '..', '..'))
    def git(*a):
        try:
            return subprocess.run(['git', '-C', repo, *a], capture_output=True,
                                  text=True, check=True).stdout.strip()
        except Exception:
            return None
    files = {
        'preflight_script': os.path.abspath(__file__),
        'wrapper': os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                'run_capture_preflight.sh'),
        'gate2_parameter_ruling': os.path.join(CODE, 'gate2_parameter_ruling.py'),
        'gate2_stress_design': os.path.join(CODE, 'gate2_stress_design.py'),
        'gate2_stress_runner': os.path.join(CODE, 'gate2_stress_runner.py'),
    }
    return {
        'repo': repo, 'commit': git('rev-parse', 'HEAD'),
        'working_tree_clean': git('status', '--porcelain') == '',
        'file_sha256': {k: (sha256(v) if os.path.exists(v) else None)
                        for k, v in files.items()},
    }


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def assert_inputs(paths):
    """Fail closed on identity BEFORE reading any evidence."""
    seen = {}
    for name, p in paths.items():
        if not os.path.exists(p):
            raise SystemExit(f'input missing: {p}')
        nbytes = os.path.getsize(p)
        if nbytes != EXPECTED_BYTES[name]:
            raise SystemExit(f'{name}: {nbytes} bytes, expected {EXPECTED_BYTES[name]}')
        got = sha256(p)
        if got != EXPECTED_SHA256[name]:
            raise SystemExit(f'{name}: sha256 {got} != expected {EXPECTED_SHA256[name]}')
        seen[name] = {'path': p, 'bytes': nbytes, 'sha256': got, 'asserted': True}
    csize = os.path.getsize(CONTAINER)
    if csize != EXPECTED_CONTAINER_BYTES:
        raise SystemExit(f'container: {csize} bytes != expected {EXPECTED_CONTAINER_BYTES}')
    cdigest = sha256(CONTAINER)
    if cdigest != EXPECTED_CONTAINER_SHA256:
        raise SystemExit(f'container: sha256 {cdigest} != expected '
                         f'{EXPECTED_CONTAINER_SHA256}')
    seen['container'] = {'path': CONTAINER, 'bytes': csize, 'sha256': cdigest,
                         'asserted': True}
    return seen


def write_atomic(path, payload):
    tmp = path + '.inprogress'
    with open(tmp, 'w') as fh:
        json.dump(payload, fh, indent=2, sort_keys=True, default=str)
        fh.flush()
        os.fsync(fh.fileno())
    os.replace(tmp, path)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--out-dir', required=True)
    args = ap.parse_args()
    if os.path.exists(args.out_dir) and os.listdir(args.out_dir):
        raise SystemExit(f'{args.out_dir} must be new or empty')
    os.makedirs(args.out_dir, exist_ok=True)
    dest = os.path.join(args.out_dir, 'CAPTURE_PREFLIGHT_MANIFEST.json')
    t0 = time.time()

    assert_ruling_version(EXPECT_RULING_VERSION)
    approved = approved_for(REQUIRED_PARAMETER_NAMES)
    provenance = repo_provenance()

    # A kill during loading, HVG or PCA must not leave nothing behind. The
    # documented lifecycle is written BEFORE the heavy work, then advanced.
    write_atomic(dest, {'preflight': 'gate2_capture_preflight/2', 'capture': CAPTURE,
                        'status': 'incomplete', 'stage': 'starting',
                        'provenance': provenance,
                        'ruling': {'version': PARAMETER_RULING_VERSION},
                        'produced_scores_calls_or_rates': False,
                        'production_parameter_ruling_written': False,
                        'started_utc': time.strftime('%Y-%m-%dT%H:%M:%SZ',
                                                     time.gmtime())})

    paths = {n: os.path.join(DEPOSIT, n) for n in EXPECTED_SHA256}
    inputs = assert_inputs(paths)

    import gzip
    import numpy as np
    import scipy.io
    import scipy.sparse as sp
    import anndata as ad
    import scanpy as sc
    import sklearn
    from scanpy.preprocessing._scrublet.core import Scrublet

    identity = assert_dependency_approved(dependency_identity(sc))

    with gzip.open(paths[f'{CAPTURE}_matrix.mtx.gz'], 'rb') as fh:
        M = scipy.io.mmread(fh).tocsr()
    with gzip.open(paths[f'{CAPTURE}_barcodes.tsv.gz'], 'rt') as fh:
        bcs = [l.strip() for l in fh if l.strip()]
    with gzip.open(paths[f'{CAPTURE}_features.tsv.gz'], 'rt') as fh:
        feats = [l.rstrip('\n').split('\t') for l in fh if l.strip()]
    a = ad.AnnData(X=sp.csr_matrix(M.T))
    a.obs_names = bcs
    a.var_names = [f[0] for f in feats]
    a.var_names_make_unique()
    shape_in = list(a.shape)

    seed = capture_seed(COHORT, CHANNEL)

    # scanpy's own scrublet preprocessing order, same as the runner's prepare_capture
    sc.pp.filter_genes(a, min_cells=3)
    sc.pp.filter_cells(a, min_genes=3)
    n_eff = int(a.n_obs)
    a.layers['raw'] = a.X.copy()
    sc.pp.normalize_total(a)
    a.layers['log1p'] = a.X.copy()
    sc.pp.log1p(a, layer='log1p')
    sc.pp.highly_variable_genes(a, layer='log1p')
    g_hvg = int(a.var['highly_variable'].sum())
    del a.layers['log1p']
    a = a[:, a.var['highly_variable']].copy()
    write_atomic(dest, {'preflight': 'gate2_capture_preflight/2', 'capture': CAPTURE,
                        'status': 'incomplete', 'stage': 'preprocessed',
                        'provenance': provenance,
                        'n_cells_effective': n_eff, 'n_usable_features_hvg': g_hvg,
                        'produced_scores_calls_or_rates': False,
                        'production_parameter_ruling_written': False})

    sim = sc.pp.scrublet_simulate_doublets(
        a, layer='raw', random_seed=seed,
        sim_doublet_ratio=approved['sim_doublet_ratio'],
        synthetic_doublet_umi_subsampling=approved['synthetic_doublet_umi_subsampling'])
    del a.layers['raw']
    if approved['log_transform']:
        sc.pp.log1p(a)
        sc.pp.log1p(sim)
    # BOTH are normalised. The previous draft normalised only the observed matrix,
    # leaving the simulated one on a different scale before stacking.
    sc.pp.normalize_total(a, target_sum=1e6)
    sc.pp.normalize_total(sim, target_sum=1e6)

    pcs = approved['n_prin_comps']
    limit = min(n_eff, g_hvg)
    feasible = pcs < limit
    n_neighbors = int(round(0.5 * np.sqrt(n_eff)))
    n_sim = int(sim.n_obs)
    k_adj = int(round(n_neighbors * (1 + n_sim / n_eff)))

    result = {
        'preflight': 'gate2_capture_preflight/2',
        'capture': CAPTURE, 'cohort': COHORT, 'channel_id': CHANNEL,
        'generated_utc': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
        'status': 'incomplete',
        'produced_scores_calls_or_rates': False,
        'outcome_columns_read': [],
        'production_parameter_ruling_written': False,
        'ruling': {'version': PARAMETER_RULING_VERSION,
                   'reference': PARAMETER_RULING_REFERENCE,
                   'imported_not_copied': True,
                   'module': 'scripts/32_rebuild_upfront_qc/code/gate2_parameter_ruling.py'},
        'provenance': provenance,
        'inputs': inputs,
        'dependency_identity': identity,
        'approved_scanpy_versions': list(APPROVED_SCANPY_VERSIONS),
        'seed': {'value': seed, 'derivation': 'capture_seed(cohort, channel_id)'},
        'input_summary': {'shape_as_loaded': shape_in, 'n_cells_effective': n_eff,
                          'n_usable_features_hvg': g_hvg, 'n_simulated': n_sim},
        'feasibility': {'rule': 'n_prin_comps < min(N_eff, G_HVG)', 'n_prin_comps': pcs,
                        'limit': limit, 'feasible': bool(feasible)},
        'parameters': dict(approved, n_neighbors=n_neighbors),
        'k_adj': {'value': k_adj,
                  'rule': 'round(n_neighbors * (1 + n_sim/n_obs)), scanpy 1.11.3'},
        'prior_grid': list(RESIDUAL_RATE_GRID), 'prior_label': PRIOR_LABEL,
        'versions': {'scanpy': sc.__version__, 'anndata': ad.__version__,
                     'numpy': np.__version__, 'scipy': scipy.__version__,
                     'sklearn': sklearn.__version__, 'python': sys.version.split()[0]},
    }
    if not feasible:
        result['stopped'] = f'PC={pcs} not < min({n_eff},{g_hvg})'
        result['status'] = 'stopped_infeasible'
        write_atomic(dest, result)
        raise SystemExit(result['stopped'])

    # ---- scanpy's OWN manifold, observed-fit, stopping before any score -----
    scrub = Scrublet(
        a.X, sim_doublet_ratio=approved['sim_doublet_ratio'],
        n_neighbors=n_neighbors, expected_doublet_rate=0.0,   # placeholder; unused here
        stdev_doublet_rate=approved['stdev_doublet_rate'], random_state=seed)
    scrub._counts_obs_norm = a.X
    scrub._counts_sim_norm = sim.X
    from scanpy.preprocessing._scrublet import pipeline as sp_pipe
    if approved['mean_center'] and approved['normalize_variance']:
        sp_pipe.zscore(scrub)
    elif approved['mean_center']:
        sp_pipe.mean_center(scrub)
    elif approved['normalize_variance']:
        sp_pipe.normalize_variance(scrub)
    sp_pipe.truncated_svd(scrub, n_prin_comps=pcs, random_state=seed) \
        if approved['mean_center'] is False else sp_pipe.pca(
            scrub, n_prin_comps=pcs, random_state=seed)

    man_obs = np.asarray(scrub.manifold_obs_)
    man_sim = np.asarray(scrub.manifold_sim_)
    result['manifold'] = {
        'built_by': 'scanpy.preprocessing._scrublet.pipeline (observed-fit)',
        'observed_shape': list(man_obs.shape), 'simulated_shape': list(man_sim.shape),
        'observed_sha256': hashlib.sha256(np.ascontiguousarray(
            np.round(man_obs, 6)).tobytes()).hexdigest(),
        'simulated_sha256': hashlib.sha256(np.ascontiguousarray(
            np.round(man_sim, 6)).tobytes()).hexdigest(),
        'all_finite': bool(np.all(np.isfinite(man_obs)) and np.all(np.isfinite(man_sim))),
        'rank_supports_n_prin_comps': bool(man_obs.shape[1] == pcs),
    }
    result['simulated_matrix_checksum'] = _checksum_matrix(sim.X)
    result['prior_independence'] = {
        'expected_doublet_rate_consumed_by': 'calculate_doublet_scores only',
        'manifold_built_once': True,
        'evidence': 'structural: no step up to the manifold takes the prior as input, '
                    'so the eight priors share one manifold by construction. A loop '
                    'over an unused rho would show only deterministic reruns.',
        'graph_checksum': None,
        'graph_checksum_reason': 'the kNN graph is built inside calculate_doublet_scores, '
                                 'which computes scores; it cannot be checksummed without '
                                 'crossing the score-free boundary',
    }
    result['status'] = 'complete'
    result['elapsed_seconds'] = round(time.time() - t0, 1)
    write_atomic(dest, result)
    print(f'N_eff={n_eff} G_HVG={g_hvg} PC={pcs} feasible={feasible} '
          f'k_adj={k_adj} manifold={man_obs.shape} status=complete')


if __name__ == '__main__':
    main()

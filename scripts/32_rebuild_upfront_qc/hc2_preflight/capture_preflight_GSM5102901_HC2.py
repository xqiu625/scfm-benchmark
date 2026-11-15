"""Single-capture, score-free Gate 2 preflight CANDIDATE for GSM5102901_HC2.

    python capture_preflight_GSM5102901_HC2.py --out-dir <NEW versioned dir>

NOT YET AUTHORISED FOR SUBMISSION. This is a candidate implementation prepared for
owner/automated review, per `review:20261012_124351` item B (5-7) and
`GATE2_NEXT_CAPTURE_DECISION_PACKET_v1.md`. It has not been run against the real
8,345-cell GSM5102901_HC2 matrix; see
`scripts/32_rebuild_upfront_qc/hc2_preflight/test_capture_preflight_hc2.py` for the
small-fixture boundary/failure-path checks that HAVE been run, on HPCC, this round.

Measures feasibility on the complete 8,345-cell deposit matrix, once run. Produces
NO doublet score, call or rate, reads NO outcome, and leaves
PRODUCTION_PARAMETER_RULING untouched.

Revised in place (2026-10-12) per `review:20261012_132816` R1/R2: the
manifold shape/finite checks were previously computed and recorded but never
gated `status`, so an invalid manifold (NaN, wrong row/column count) could still
be written as `complete`; and the wrapper's/manifest's tracked-dependency list
omitted `gate2_diagnostics.py`, which `gate2_stress_runner.py` actually imports.
Both are fixed below -- see `run_preflight`'s manifold-validity gate and
`repo_provenance`'s `files` dict.

DIRECT DERIVATIVE OF THE REVIEWED GSM5102900_HC1 PREFLIGHT
------------------------------------------------------------
This script is structurally identical to
`reports/manuscript/revision/gate2_capture_audit_2026-09-24/capture_preflight_GSM5102900_HC1.py`
(that file is untouched by this change). The only substantive differences are:

  1. Capture identity (CAPTURE/COHORT/CHANNEL) and the three deposit input hashes,
     which are GSM5102901_HC2's own, independently measured this round
     (`gate2_pilot_interpretation_evidence_20261012/gate2_hc2_candidate_evidence.json`)
     -- NOT copied from GSM5102900_HC1.
  2. The relative path from this file to `scripts/32_rebuild_upfront_qc/code` is ONE
     directory level (`../code`), not four -- this script lives in
     `scripts/32_rebuild_upfront_qc/hc2_preflight/`, not four levels down inside
     `reports/manuscript/revision/gate2_capture_audit_2026-09-24/`. A naive copy of
     the original file's HERE-relative path would silently point at a nonexistent
     directory from here; this is fixed below, not inherited.
  3. No measured value (n_cells_effective, n_usable_features_hvg, n_neighbors, k_adj,
     feasibility) is pinned in advance. All of them are computed fresh from
     GSM5102901_HC2's own data when this script actually runs, exactly as the
     original script computes them fresh for GSM5102900_HC1 -- neither script
     hardcodes these as inputs, so there is nothing here to "copy over" by mistake.

Everything else -- the Scrublet pipeline call sequence (mean_center/normalize_variance
via `zscore`, `truncated_svd`/`pca` from the OBSERVED fit, stopping at
`set_manifold`), the eleven approved parameters imported (not copied) from
`gate2_parameter_ruling`, the seed derived from `capture_seed(cohort, channel_id)`,
and the score-free boundary rationale -- is identical to the reviewed HC1 script and
is not re-explained here; see that file's own docstring.

SCORE-FREE BOUNDARY, ENFORCED
------------------------------
This script imports `scanpy.preprocessing._scrublet.core.Scrublet` and calls:
  - `Scrublet(...)` (construction only)
  - `scanpy.preprocessing._scrublet.pipeline.{zscore,mean_center,normalize_variance,
    truncated_svd,pca}` (manifold construction only)
It does NOT call, import by name, or reference `Scrublet.calculate_doublet_scores` or
`Scrublet.call_doublets` anywhere in this file. The small-fixture test suite asserts
this by monkeypatching both methods to raise if called and running the full script
against a tiny synthetic input -- see the test file for the exact assertion.
"""
import argparse
import hashlib
import json
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
# scripts/32_rebuild_upfront_qc/hc2_preflight/ -> scripts/32_rebuild_upfront_qc/code
# ONE level up, not four -- see docstring point 2.
CODE = os.path.abspath(os.path.join(HERE, '..', 'code'))
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
CAPTURE = 'GSM5102901_HC2'
COHORT, CHANNEL = 'qiu', 'GSM5102901'
DEPOSIT = '/bigdata/godziklab/shared/Xinru/302006/01_DATA/raw/GSE167363_qiu/extracted'
CONTAINER = ('/bigdata/godziklab/shared/Xinru/302006/01_DATA/raw/'
             'human_platelet_covid_severity_raw_counts.h5ad')
PRIOR_LABEL = 'assumed_algorithmic_stress_setting'

#: Full SHA-256 of every input, GSM5102901_HC2's own -- independently measured this
#: round by gate2_hc2_candidate_evidence.py, recorded in
#: gate2_pilot_interpretation_evidence_20261012/gate2_hc2_candidate_evidence.json.
#: NOT GSM5102900_HC1's hashes.
EXPECTED_SHA256 = {
    f'{CAPTURE}_barcodes.tsv.gz':
        '3113aede3fcf401d613d02ef27ab6ab3b05b568b123e21ea57cf19d1950ee644',
    f'{CAPTURE}_features.tsv.gz':
        '0774e088c47a7cd2cb99a22c7d41c50405940ec0f84682c6f1f6a7652cd23587',
    f'{CAPTURE}_matrix.mtx.gz':
        '8181ebe77603c8c607ebcf9264f18c43b8de4bf8e63bce2efb5e1768054d8fe5',
}
#: The container's FULL sha256 -- same file as GSM5102900_HC1's preflight/audit used
#: (the canonical container is shared across captures), re-asserted here rather than
#: imported, so this file's own identity contract does not depend on the HC1 file.
EXPECTED_CONTAINER_BYTES = 5027039977
EXPECTED_CONTAINER_SHA256 = (
    '1a0f37902cd33f410a71428ea45b3d1573fddec393ce27f07a93074d0ed536ba')
EXPECTED_BYTES = {f'{CAPTURE}_barcodes.tsv.gz': 42960,
                  f'{CAPTURE}_features.tsv.gz': 304728,
                  f'{CAPTURE}_matrix.mtx.gz': 39369513}


def repo_provenance():
    """Commit and per-file hashes INSIDE the manifest, not only on stdout.

    `working_tree_clean` is a whole-repo `git status --porcelain` check, a
    DIFFERENT and strictly stronger guarantee than the wrapper's per-file
    tracked/staged/unstaged check against `SCRIPT`/`WRAPPER`/`IMPORTED` only.
    A clean run here does not certify the wrapper checked every file that
    happened to be dirty elsewhere in the repo, only that none were dirty at
    all; conversely a wrapper pass does not imply this is True. Do not treat
    the two as interchangeable evidence.
    """
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
                                'run_capture_preflight_HC2.sh'),
        'gate2_parameter_ruling': os.path.join(CODE, 'gate2_parameter_ruling.py'),
        'gate2_stress_design': os.path.join(CODE, 'gate2_stress_design.py'),
        'gate2_stress_runner': os.path.join(CODE, 'gate2_stress_runner.py'),
        # gate2_stress_runner imports DIAGNOSTICS_SCHEMA_VERSION/score_distribution_summary
        # from this module at load time -- confirmed by grep against the live file,
        # not assumed -- so it decides what runs here too. Missing from this list
        # was R2 in review:20261012_132816.
        'gate2_diagnostics': os.path.join(CODE, 'gate2_diagnostics.py'),
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


def assert_inputs(paths, deposit_dir=DEPOSIT, container=CONTAINER,
                   expected_bytes=EXPECTED_BYTES, expected_sha256=EXPECTED_SHA256,
                   expected_container_bytes=EXPECTED_CONTAINER_BYTES,
                   expected_container_sha256=EXPECTED_CONTAINER_SHA256):
    """Fail closed on identity BEFORE reading any evidence.

    Takes its expectations as parameters (defaulting to the module constants) so
    the small-fixture test suite can exercise this exact function against tiny
    synthetic files and their own expected hashes, without touching the real
    8,345-cell deposit.
    """
    seen = {}
    for name, p in paths.items():
        if not os.path.exists(p):
            raise SystemExit(f'input missing: {p}')
        nbytes = os.path.getsize(p)
        if nbytes != expected_bytes[name]:
            raise SystemExit(f'{name}: {nbytes} bytes, expected {expected_bytes[name]}')
        got = sha256(p)
        if got != expected_sha256[name]:
            raise SystemExit(f'{name}: sha256 {got} != expected {expected_sha256[name]}')
        seen[name] = {'path': p, 'bytes': nbytes, 'sha256': got, 'asserted': True}
    csize = os.path.getsize(container)
    if csize != expected_container_bytes:
        raise SystemExit(f'container: {csize} bytes != expected {expected_container_bytes}')
    cdigest = sha256(container)
    if cdigest != expected_container_sha256:
        raise SystemExit(f'container: sha256 {cdigest} != expected '
                         f'{expected_container_sha256}')
    seen['container'] = {'path': container, 'bytes': csize, 'sha256': cdigest,
                         'asserted': True}
    return seen


def write_atomic(path, payload):
    tmp = path + '.inprogress'
    with open(tmp, 'w') as fh:
        json.dump(payload, fh, indent=2, sort_keys=True, default=str)
        fh.flush()
        os.fsync(fh.fileno())
    os.replace(tmp, path)


def run_preflight(*, deposit_dir, container, expected_bytes, expected_sha256,
                   expected_container_bytes, expected_container_sha256,
                   out_dir, capture=CAPTURE, cohort=COHORT, channel=CHANNEL):
    """The full preflight body, parameterized over inputs/expectations so the
    small-fixture test suite can run it end-to-end against synthetic data. The
    CLI `main()` below calls this with the real module-level constants.
    """
    if os.path.exists(out_dir) and os.listdir(out_dir):
        raise SystemExit(f'{out_dir} must be new or empty')
    os.makedirs(out_dir, exist_ok=True)
    dest = os.path.join(out_dir, 'CAPTURE_PREFLIGHT_MANIFEST.json')
    t0 = time.time()

    assert_ruling_version(EXPECT_RULING_VERSION)
    approved = approved_for(REQUIRED_PARAMETER_NAMES)
    provenance = repo_provenance()

    write_atomic(dest, {'preflight': 'gate2_capture_preflight_candidate/1', 'capture': capture,
                        'status': 'incomplete', 'stage': 'starting',
                        'provenance': provenance,
                        'ruling': {'version': PARAMETER_RULING_VERSION},
                        'produced_scores_calls_or_rates': False,
                        'production_parameter_ruling_written': False,
                        'started_utc': time.strftime('%Y-%m-%dT%H:%M:%SZ',
                                                     time.gmtime())})

    paths = {n: os.path.join(deposit_dir, n) for n in expected_sha256}
    inputs = assert_inputs(paths, deposit_dir=deposit_dir, container=container,
                            expected_bytes=expected_bytes, expected_sha256=expected_sha256,
                            expected_container_bytes=expected_container_bytes,
                            expected_container_sha256=expected_container_sha256)

    import gzip
    import numpy as np
    import scipy.io
    import scipy.sparse as sp
    import anndata as ad
    import scanpy as sc
    import sklearn
    from scanpy.preprocessing._scrublet.core import Scrublet

    identity = assert_dependency_approved(dependency_identity(sc))

    with gzip.open(paths[f'{capture}_matrix.mtx.gz'], 'rb') as fh:
        M = scipy.io.mmread(fh).tocsr()
    with gzip.open(paths[f'{capture}_barcodes.tsv.gz'], 'rt') as fh:
        bcs = [l.strip() for l in fh if l.strip()]
    with gzip.open(paths[f'{capture}_features.tsv.gz'], 'rt') as fh:
        feats = [l.rstrip('\n').split('\t') for l in fh if l.strip()]
    a = ad.AnnData(X=sp.csr_matrix(M.T))
    a.obs_names = bcs
    a.var_names = [f[0] for f in feats]
    a.var_names_make_unique()
    shape_in = list(a.shape)

    seed = capture_seed(cohort, channel)

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
    write_atomic(dest, {'preflight': 'gate2_capture_preflight_candidate/1', 'capture': capture,
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
    sc.pp.normalize_total(a, target_sum=1e6)
    sc.pp.normalize_total(sim, target_sum=1e6)

    pcs = approved['n_prin_comps']
    limit = min(n_eff, g_hvg)
    feasible = pcs < limit
    n_neighbors = int(round(0.5 * np.sqrt(n_eff)))
    n_sim = int(sim.n_obs)
    k_adj = int(round(n_neighbors * (1 + n_sim / n_eff)))

    result = {
        'preflight': 'gate2_capture_preflight_candidate/1',
        'capture': capture, 'cohort': cohort, 'channel_id': channel,
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
        'neighbor_graph': {'verified': False,
                            'reason': 'the kNN graph is built inside calculate_doublet_scores, '
                                      'which this preflight does not call; see score-free '
                                      'boundary note in this script\'s docstring'},
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
    all_finite = bool(np.all(np.isfinite(man_obs)) and np.all(np.isfinite(man_sim)))
    # Declared-shape contract only -- NOT a numerical-rank check. A shape match does
    # not prove the columns are linearly independent; introducing an actual rank
    # screen would be a new exclusion criterion requiring its own tolerance/rationale
    # (per review:20261012_132816 R1), which this fix does not do.
    expected_observed_shape = (n_eff, pcs)
    expected_simulated_shape = (n_sim, pcs)
    observed_shape_matches_expected = (man_obs.ndim == 2
                                       and tuple(man_obs.shape) == expected_observed_shape)
    simulated_shape_matches_expected = (man_sim.ndim == 2
                                        and tuple(man_sim.shape) == expected_simulated_shape)
    manifold_valid = (all_finite and observed_shape_matches_expected
                      and simulated_shape_matches_expected)
    result['manifold'] = {
        'built_by': 'scanpy.preprocessing._scrublet.pipeline (observed-fit)',
        'observed_shape': list(man_obs.shape), 'simulated_shape': list(man_sim.shape),
        'expected_observed_shape': list(expected_observed_shape),
        'expected_simulated_shape': list(expected_simulated_shape),
        'observed_sha256': hashlib.sha256(np.ascontiguousarray(
            np.round(man_obs, 6)).tobytes()).hexdigest(),
        'simulated_sha256': hashlib.sha256(np.ascontiguousarray(
            np.round(man_sim, 6)).tobytes()).hexdigest(),
        'all_finite': all_finite,
        'observed_shape_matches_expected': observed_shape_matches_expected,
        'simulated_shape_matches_expected': simulated_shape_matches_expected,
    }
    result['simulated_matrix_checksum'] = _checksum_matrix(sim.X)
    result['prior_independence'] = {
        'expected_doublet_rate_consumed_by': 'calculate_doublet_scores only',
        'manifold_built_once': True,
        'evidence': 'structural: no step up to the manifold takes the prior as input',
        'graph_checksum': None,
        'graph_checksum_reason': 'the kNN graph is built inside calculate_doublet_scores, '
                                 'which computes scores; it cannot be checksummed without '
                                 'crossing the score-free boundary',
    }
    if not manifold_valid:
        result['stopped'] = (
            f'manifold validation failed: all_finite={all_finite} '
            f'observed_shape_matches_expected={observed_shape_matches_expected} '
            f'simulated_shape_matches_expected={simulated_shape_matches_expected}')
        result['status'] = 'stopped_invalid_manifold'
        result['elapsed_seconds'] = round(time.time() - t0, 1)
        write_atomic(dest, result)
        raise SystemExit(result['stopped'])
    result['status'] = 'complete'
    result['elapsed_seconds'] = round(time.time() - t0, 1)
    write_atomic(dest, result)
    print(f'N_eff={n_eff} G_HVG={g_hvg} PC={pcs} feasible={feasible} '
          f'k_adj={k_adj} manifold={man_obs.shape} status=complete')
    return result


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--out-dir', required=True)
    args = ap.parse_args()
    run_preflight(deposit_dir=DEPOSIT, container=CONTAINER,
                  expected_bytes=EXPECTED_BYTES, expected_sha256=EXPECTED_SHA256,
                  expected_container_bytes=EXPECTED_CONTAINER_BYTES,
                  expected_container_sha256=EXPECTED_CONTAINER_SHA256,
                  out_dir=args.out_dir)


if __name__ == '__main__':
    main()

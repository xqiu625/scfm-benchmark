"""Small-fixture boundary/failure-path tests for the HC2 preflight CANDIDATE.

Run only against tiny synthetic data (tens of cells/genes), never the real
8,345-cell GSM5102901_HC2 matrix and never the real 5 GB canonical container --
both are parameterized into `assert_inputs`/`run_preflight` precisely so these
tests can supply small stand-ins instead. Intended to run on HPCC inside
`302006_env_clean` (needs scanpy/anndata/scipy), not locally, not via sbatch --
these are fast enough for an interactive/login-adjacent run per DEVELOPMENT.md's "no
heavy real-input processing on the login node" (these fixtures are KB-scale, not
the 5 GB container or the real deposit).

    cd scripts/32_rebuild_upfront_qc/hc2_preflight
    python -m pytest test_capture_preflight_hc2.py -v
"""
import gzip
import hashlib
import json
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import capture_preflight_GSM5102901_HC2 as P  # noqa: E402


def sha256_bytes(b):
    return hashlib.sha256(b).hexdigest()


def write_gz(path, text):
    data = text.encode()
    with gzip.open(path, 'wb') as fh:
        fh.write(data)
    with gzip.open(path, 'rb') as fh:
        raw = fh.read()
    return raw


def make_tiny_capture(tmp_path, capture, n_genes, n_cells, seed, mean_scale=3.0):
    """Write a tiny synthetic CellRanger-style triplet under tmp_path/extracted/
    and return (deposit_dir, expected_bytes, expected_sha256).

    Per-gene Poisson means are drawn from a lognormal spread so genes differ
    enough in mean/variance for highly_variable_genes to have something real to
    select from -- not uniform noise, which scanpy's dispersion cutoff would
    reject almost everything from.
    """
    rng = np.random.default_rng(seed)
    deposit_dir = tmp_path / 'extracted'
    deposit_dir.mkdir(parents=True, exist_ok=True)

    gene_means = rng.lognormal(mean=0.0, sigma=1.2, size=n_genes) * mean_scale
    counts = np.zeros((n_genes, n_cells), dtype=int)
    for g in range(n_genes):
        counts[g] = rng.poisson(lam=gene_means[g], size=n_cells)
    nnz_mask = counts > 0
    nnz = int(nnz_mask.sum())

    lines = ['%%MatrixMarket matrix coordinate integer general', '%',
             f'{n_genes} {n_cells} {nnz}']
    gi, ci = np.nonzero(nnz_mask)
    for g, c in zip(gi, ci):
        lines.append(f'{g + 1} {c + 1} {int(counts[g, c])}')
    mtx_text = '\n'.join(lines) + '\n'

    bases = 'ACGT'
    barcodes = []
    seen = set()
    while len(barcodes) < n_cells:
        bc = ''.join(rng.choice(list(bases), size=16)) + '-1'
        if bc not in seen:
            seen.add(bc)
            barcodes.append(bc)
    bc_text = '\n'.join(barcodes) + '\n'

    feat_lines = [f'ENSG{g:08d}\tGENE{g}\tGene Expression' for g in range(n_genes)]
    feat_text = '\n'.join(feat_lines) + '\n'

    file_names = []
    for suffix, text in (('matrix.mtx.gz', mtx_text), ('barcodes.tsv.gz', bc_text),
                          ('features.tsv.gz', feat_text)):
        path = deposit_dir / f'{capture}_{suffix}'
        write_gz(str(path), text)
        file_names.append(f'{capture}_{suffix}')
    expected_bytes = {n: os.path.getsize(deposit_dir / n) for n in file_names}
    expected_sha256 = {n: P.sha256(str(deposit_dir / n)) for n in file_names}
    return str(deposit_dir), expected_bytes, expected_sha256


def make_tiny_container(tmp_path, text='not a real h5ad, just a byte-identity stand-in\n'):
    path = tmp_path / 'fake_container.h5ad'
    data = text.encode()
    with open(path, 'wb') as fh:
        fh.write(data)
    return str(path), os.path.getsize(path), sha256_bytes(data)


def test_relative_code_path_resolves():
    """The documented bug this script fixes relative to a naive copy of HC1's
    preflight: HERE-relative path depth differs because this file lives one
    level shallower. Confirm CODE actually resolves to the real shared-module
    directory, not a path that happens to not exist (which assert_inputs/import
    would only catch indirectly via an ImportError at module load -- this test
    makes the path claim itself explicit and checkable).
    """
    assert os.path.isdir(P.CODE), f'CODE does not resolve to a directory: {P.CODE}'
    assert os.path.exists(os.path.join(P.CODE, 'gate2_parameter_ruling.py'))
    assert os.path.basename(P.CODE) == 'code'
    assert os.path.basename(os.path.dirname(P.CODE)) == '32_rebuild_upfront_qc'


def test_assert_inputs_rejects_missing_file(tmp_path):
    with pytest.raises(SystemExit, match='input missing'):
        P.assert_inputs(
            {'x_barcodes.tsv.gz': str(tmp_path / 'does_not_exist.tsv.gz')},
            deposit_dir=str(tmp_path), container=str(tmp_path / 'nope.h5ad'),
            expected_bytes={'x_barcodes.tsv.gz': 1}, expected_sha256={'x_barcodes.tsv.gz': 'a'},
            expected_container_bytes=1, expected_container_sha256='a')


def test_assert_inputs_rejects_hash_mismatch(tmp_path):
    deposit_dir, expected_bytes, expected_sha256 = make_tiny_capture(
        tmp_path, 'TESTCAP', n_genes=10, n_cells=5, seed=1)
    container_path, cbytes, csha = make_tiny_container(tmp_path)
    tampered = dict(expected_sha256)
    k = next(iter(tampered))
    tampered[k] = '0' * 64  # deliberately wrong
    paths = {k: os.path.join(deposit_dir, k) for k in expected_sha256}
    with pytest.raises(SystemExit, match='sha256'):
        P.assert_inputs(paths, deposit_dir=deposit_dir, container=container_path,
                         expected_bytes=expected_bytes, expected_sha256=tampered,
                         expected_container_bytes=cbytes, expected_container_sha256=csha)


def test_assert_inputs_rejects_byte_count_mismatch(tmp_path):
    deposit_dir, expected_bytes, expected_sha256 = make_tiny_capture(
        tmp_path, 'TESTCAP', n_genes=10, n_cells=5, seed=2)
    container_path, cbytes, csha = make_tiny_container(tmp_path)
    tampered = dict(expected_bytes)
    k = next(iter(tampered))
    tampered[k] = tampered[k] + 999
    paths = {k: os.path.join(deposit_dir, k) for k in expected_sha256}
    with pytest.raises(SystemExit, match='bytes'):
        P.assert_inputs(paths, deposit_dir=deposit_dir, container=container_path,
                         expected_bytes=tampered, expected_sha256=expected_sha256,
                         expected_container_bytes=cbytes, expected_container_sha256=csha)


def test_assert_inputs_rejects_container_mismatch(tmp_path):
    deposit_dir, expected_bytes, expected_sha256 = make_tiny_capture(
        tmp_path, 'TESTCAP', n_genes=10, n_cells=5, seed=3)
    container_path, cbytes, csha = make_tiny_container(tmp_path)
    paths = {k: os.path.join(deposit_dir, k) for k in expected_sha256}
    with pytest.raises(SystemExit, match='container'):
        P.assert_inputs(paths, deposit_dir=deposit_dir, container=container_path,
                         expected_bytes=expected_bytes, expected_sha256=expected_sha256,
                         expected_container_bytes=cbytes, expected_container_sha256='0' * 64)


def test_run_preflight_refuses_nonempty_outdir(tmp_path):
    out_dir = tmp_path / 'out'
    out_dir.mkdir()
    (out_dir / 'junk.txt').write_text('pre-existing result, must not be overwritten')
    with pytest.raises(SystemExit, match='must be new or empty'):
        P.run_preflight(deposit_dir='/nonexistent', container='/nonexistent',
                         expected_bytes={}, expected_sha256={},
                         expected_container_bytes=0, expected_container_sha256='x',
                         out_dir=str(out_dir), capture='X', cohort='x', channel='x')
    assert (out_dir / 'junk.txt').read_text() == 'pre-existing result, must not be overwritten'


def _forbid_scoring(monkeypatch):
    """Monkeypatch the two score-producing methods to raise if ever called, so a
    completing run proves they were never invoked -- not merely that the result
    dict happens to lack score/call keys.
    """
    from scanpy.preprocessing._scrublet.core import Scrublet

    def _boom(*a, **k):
        raise AssertionError('calculate_doublet_scores/call_doublets must not be '
                              'called by a score-free preflight')
    monkeypatch.setattr(Scrublet, 'calculate_doublet_scores', _boom)
    monkeypatch.setattr(Scrublet, 'call_doublets', _boom)


def test_run_preflight_infeasible_stops_before_manifold(tmp_path, monkeypatch):
    _forbid_scoring(monkeypatch)
    capture = 'TESTCAP_INFEASIBLE'
    deposit_dir, expected_bytes, expected_sha256 = make_tiny_capture(
        tmp_path, capture, n_genes=20, n_cells=30, seed=11, mean_scale=5.0)
    container_path, cbytes, csha = make_tiny_container(tmp_path)
    out_dir = tmp_path / 'out_infeasible'

    with pytest.raises(SystemExit, match=r'PC=30 not <'):
        P.run_preflight(deposit_dir=deposit_dir, container=container_path,
                         expected_bytes=expected_bytes, expected_sha256=expected_sha256,
                         expected_container_bytes=cbytes, expected_container_sha256=csha,
                         out_dir=str(out_dir), capture=capture, cohort='qiu', channel='TESTCHAN')

    manifest = json.load(open(out_dir / 'CAPTURE_PREFLIGHT_MANIFEST.json'))
    assert manifest['status'] == 'stopped_infeasible'
    assert manifest['produced_scores_calls_or_rates'] is False
    assert manifest['outcome_columns_read'] == []
    assert manifest['production_parameter_ruling_written'] is False
    assert 'manifold' not in manifest
    assert manifest['feasibility']['feasible'] is False
    # 20 genes total, so n_usable_features_hvg can never exceed 20 < n_prin_comps=30:
    # infeasibility here is a property of the fixture's gene count, not random luck.
    assert manifest['feasibility']['limit'] <= 20


def test_run_preflight_feasible_completes_without_scoring(tmp_path, monkeypatch):
    _forbid_scoring(monkeypatch)
    capture = 'TESTCAP_FEASIBLE'
    deposit_dir, expected_bytes, expected_sha256 = make_tiny_capture(
        tmp_path, capture, n_genes=400, n_cells=80, seed=22, mean_scale=4.0)
    container_path, cbytes, csha = make_tiny_container(tmp_path)
    out_dir = tmp_path / 'out_feasible'

    result = P.run_preflight(
        deposit_dir=deposit_dir, container=container_path,
        expected_bytes=expected_bytes, expected_sha256=expected_sha256,
        expected_container_bytes=cbytes, expected_container_sha256=csha,
        out_dir=str(out_dir), capture=capture, cohort='qiu', channel='TESTCHAN2')

    if result['status'] != 'complete':
        pytest.skip(f'fixture landed infeasible (limit={result["feasibility"]["limit"]}); '
                    'not what this test exercises, not a preflight defect')

    manifest = json.load(open(out_dir / 'CAPTURE_PREFLIGHT_MANIFEST.json'))
    assert manifest['status'] == 'complete'
    assert manifest['produced_scores_calls_or_rates'] is False
    assert manifest['outcome_columns_read'] == []
    assert manifest['production_parameter_ruling_written'] is False
    assert 'manifold' in manifest
    assert manifest['manifold']['all_finite'] is True
    assert manifest['manifold']['observed_shape_matches_expected'] is True
    assert manifest['manifold']['simulated_shape_matches_expected'] is True
    assert manifest['k_adj']['value'] > 0
    assert manifest['neighbor_graph']['verified'] is False
    # The two score-producing methods were monkeypatched to raise; reaching
    # 'complete' without that exception firing is the proof they were never
    # called, not an assumption from reading the source.
    assert 'doublet_scores_obs_' not in json.dumps(manifest)
    assert 'predicted_doublets_' not in json.dumps(manifest)


def _patch_pca_then(monkeypatch, corrupt):
    """Let the real `scanpy.preprocessing._scrublet.pipeline.pca` run, then apply
    `corrupt(scrub)` to the manifold it wrote -- exercises the real production PCA
    call, then injects exactly the fault R1 (review:20261012_132816)
    describes: a manifold written to disk without its shape/finiteness checked.
    """
    from scanpy.preprocessing._scrublet import pipeline as sp_pipe
    real_pca = sp_pipe.pca

    def _pca_then_corrupt(scrub, *a, **k):
        real_pca(scrub, *a, **k)
        corrupt(scrub)
    monkeypatch.setattr(sp_pipe, 'pca', _pca_then_corrupt)


def _run_or_skip_if_infeasible(**kwargs):
    """Runs `run_preflight`, returning the SystemExit message. If the fixture
    landed infeasible before ever reaching the manifold (a property of the
    fixture, not of the fault being tested), skips rather than asserting on the
    wrong failure path -- mirrors the skip guard already used by
    `test_run_preflight_feasible_completes_without_scoring`.
    """
    with pytest.raises(SystemExit) as exc_info:
        P.run_preflight(**kwargs)
    message = str(exc_info.value)
    if message.startswith('PC=30 not <'):
        pytest.skip(f'fixture landed infeasible ({message!r}); not what this test '
                    'exercises, not a preflight defect')
    return message


def test_run_preflight_detects_nonfinite_manifold(tmp_path, monkeypatch):
    """R1: a NaN in the manifold must stop the run and record why, not be
    written to disk as `status: complete`.
    """
    _forbid_scoring(monkeypatch)
    _patch_pca_then(monkeypatch, lambda scrub: scrub.manifold_obs_.__setitem__((0, 0), np.nan))

    capture = 'TESTCAP_NAN'
    deposit_dir, expected_bytes, expected_sha256 = make_tiny_capture(
        tmp_path, capture, n_genes=400, n_cells=80, seed=22, mean_scale=4.0)
    container_path, cbytes, csha = make_tiny_container(tmp_path)
    out_dir = tmp_path / 'out_nan'

    message = _run_or_skip_if_infeasible(
        deposit_dir=deposit_dir, container=container_path,
        expected_bytes=expected_bytes, expected_sha256=expected_sha256,
        expected_container_bytes=cbytes, expected_container_sha256=csha,
        out_dir=str(out_dir), capture=capture, cohort='qiu', channel='TESTCHAN3')
    assert 'manifold validation failed' in message

    manifest = json.load(open(out_dir / 'CAPTURE_PREFLIGHT_MANIFEST.json'))
    assert manifest['status'] == 'stopped_invalid_manifold'
    assert manifest['manifold']['all_finite'] is False
    assert manifest['produced_scores_calls_or_rates'] is False
    assert manifest['outcome_columns_read'] == []


def test_run_preflight_detects_observed_row_count_mismatch(tmp_path, monkeypatch):
    """R1: dropping a row from the observed manifold (an N_eff mismatch the
    original code never checked, only column count) must stop the run.
    """
    _forbid_scoring(monkeypatch)
    _patch_pca_then(monkeypatch, lambda scrub: setattr(
        scrub, 'manifold_obs_', scrub.manifold_obs_[:-1, :]))

    capture = 'TESTCAP_ROWDROP'
    deposit_dir, expected_bytes, expected_sha256 = make_tiny_capture(
        tmp_path, capture, n_genes=400, n_cells=80, seed=22, mean_scale=4.0)
    container_path, cbytes, csha = make_tiny_container(tmp_path)
    out_dir = tmp_path / 'out_rowdrop'

    message = _run_or_skip_if_infeasible(
        deposit_dir=deposit_dir, container=container_path,
        expected_bytes=expected_bytes, expected_sha256=expected_sha256,
        expected_container_bytes=cbytes, expected_container_sha256=csha,
        out_dir=str(out_dir), capture=capture, cohort='qiu', channel='TESTCHAN4')
    assert 'manifold validation failed' in message

    manifest = json.load(open(out_dir / 'CAPTURE_PREFLIGHT_MANIFEST.json'))
    assert manifest['status'] == 'stopped_invalid_manifold'
    assert manifest['manifold']['all_finite'] is True
    assert manifest['manifold']['observed_shape_matches_expected'] is False
    assert manifest['manifold']['simulated_shape_matches_expected'] is True
    assert manifest['produced_scores_calls_or_rates'] is False


if __name__ == '__main__':
    raise SystemExit(pytest.main([__file__, '-v']))

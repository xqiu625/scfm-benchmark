"""Tests for `embedding_providers.PCAEmbeddingProvider` -- the explicit PCA
CONTROL baseline, per `review:20261013_135634` section 3. Every
fixture here is a tiny, synthetic, local h5ad -- no network, no real data,
no scFM inference. Reuses `dataset_schema`/`cli`'s own input validation
(never a bare file path).
"""

import os
import sys

import numpy as np
import pytest

CODE = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "scripts",
    "34_scfm_discovery_framework",
    "code",
)
sys.path.insert(0, CODE)

import embedding_providers as E  # noqa: E402


def _base_config(tmp_path, matrix_path, **overrides):
    config = {
        "dataset_id": "test_pca_pilot",
        "species": "human",
        "tissue_or_target_cell_type": "monocyte",
        "subject_id_column": "donor_id",
        "sample_id_column": "sample_id",
        "biological_replicate_unit": "donor_id",
        "technical_batch_column": "chip_id",
        "expression_matrix_path": matrix_path,
        "counts_semantics": "raw_counts",
        "counts_slot": "X",
        "gene_id_namespace": "ensembl_gene_id",
        "cell_id_column": "cell_barcode",
        "qc_status": "provisional",
        "qc_version_ref": "none (provisional)",
    }
    config.update(overrides)
    return config


def _make_tiny_h5ad(tmp_path, *, n_vars=6, n_obs=10, name="tiny.h5ad", seed=0, low=0, high=5):
    import anndata as ad
    import pandas as pd

    cell_barcodes = [f"CELL{i}" for i in range(n_obs)]
    donor_ids = [f"D{i % 2}" for i in range(n_obs)]
    sample_ids = [f"S{i}" for i in range(n_obs)]
    chip_ids = ["C1"] * n_obs
    X = np.random.default_rng(seed).integers(low, high, size=(n_obs, n_vars)).astype("float32")
    obs = pd.DataFrame(
        {
            "cell_barcode": cell_barcodes,
            "donor_id": donor_ids,
            "sample_id": sample_ids,
            "chip_id": chip_ids,
        },
        index=[f"idx{i}" for i in range(n_obs)],
    )
    var = pd.DataFrame(index=[f"GENE{i}" for i in range(n_vars)])
    a = ad.AnnData(X=X, obs=obs, var=var)
    path = tmp_path / name
    a.write_h5ad(path)
    return str(path)


def _fake_config_path(tmp_path):
    # expression_matrix_path is always absolute in these fixtures, so this
    # path's own existence/dirname is never actually used by
    # resolve_matrix_path -- it only needs to be a string.
    return str(tmp_path / "dataset_config.json")


def test_pca_fit_transform_produces_a_readable_roundtrip_artifact(tmp_path):
    h5ad = _make_tiny_h5ad(tmp_path, n_vars=8, n_obs=12)
    config = _base_config(tmp_path, h5ad)
    provider = E.PCAEmbeddingProvider()
    result = provider.fit_transform(
        config,
        _fake_config_path(tmp_path),
        out_dir=str(tmp_path / "pca_out"),
        n_components=3,
        random_seed=0,
        svd_solver="full",
    )

    assert result.status == "complete", result.reason
    assert len(result.cell_ids) == 12
    assert len(result.feature_ids) == 8
    assert result.n_components == 3
    assert len(result.explained_variance_ratio) == 3

    manifest = E.read_embedding_result(str(tmp_path / "pca_out" / "EMBEDDING_RESULT.json"))
    assert manifest["status"] == "complete"
    arrays = E.load_embedding_matrices(manifest)
    assert arrays["embedding"].shape == (12, 3)
    assert arrays["components"].shape == (3, 8)
    assert arrays["mean"].shape == (8,)
    assert np.all(np.isfinite(arrays["embedding"]))


def test_pca_frozen_params_reproduce_bit_identical_embedding(tmp_path):
    """Per the review's own 'frozen parameters reproduce' requirement: the
    SAME dataset + SAME params run twice (different out_dirs) must produce
    numerically identical embeddings.
    """
    h5ad = _make_tiny_h5ad(tmp_path, n_vars=8, n_obs=12)
    config = _base_config(tmp_path, h5ad)
    provider = E.PCAEmbeddingProvider()
    kwargs = dict(n_components=3, random_seed=0, svd_solver="full", log1p=True)
    r1 = provider.fit_transform(
        config, _fake_config_path(tmp_path), out_dir=str(tmp_path / "run1"), **kwargs
    )
    r2 = provider.fit_transform(
        config, _fake_config_path(tmp_path), out_dir=str(tmp_path / "run2"), **kwargs
    )
    assert r1.status == r2.status == "complete"
    a1 = E.load_embedding_matrices(
        E.read_embedding_result(str(tmp_path / "run1" / "EMBEDDING_RESULT.json"))
    )
    a2 = E.load_embedding_matrices(
        E.read_embedding_result(str(tmp_path / "run2" / "EMBEDDING_RESULT.json"))
    )
    np.testing.assert_array_equal(a1["embedding"], a2["embedding"])
    assert r1.input_counts_sha256 == r2.input_counts_sha256


def test_pca_row_correspondence_survives_cell_id_reordering(tmp_path):
    """Independently reconstructs the embedding (embedding @ components +
    mean, lossless at full rank) and compares it against a SEPARATE,
    production-independent read of the raw h5ad (via `anndata.read_h5ad`,
    never the provider's own reading code) for each cell IN THE EXACT
    SHUFFLED ORDER requested -- proves a reordered cell_ids_subset is never
    silently misaligned against the matrix rows.
    """
    import anndata as ad

    h5ad = _make_tiny_h5ad(tmp_path, n_vars=6, n_obs=8)
    truth = ad.read_h5ad(h5ad)
    truth.obs_names = list(truth.obs["cell_barcode"])
    all_cell_ids = list(truth.obs_names)
    shuffled = list(reversed(all_cell_ids))
    assert shuffled != all_cell_ids

    config = _base_config(tmp_path, h5ad)
    provider = E.PCAEmbeddingProvider()
    result = provider.fit_transform(
        config,
        _fake_config_path(tmp_path),
        out_dir=str(tmp_path / "pca_out"),
        n_components=6,
        random_seed=0,
        svd_solver="full",
        cell_ids_subset=tuple(shuffled),
    )
    assert result.status == "complete", result.reason
    assert result.cell_ids == tuple(shuffled)

    manifest = E.read_embedding_result(str(tmp_path / "pca_out" / "EMBEDDING_RESULT.json"))
    arrays = E.load_embedding_matrices(manifest)
    reconstructed = arrays["embedding"] @ arrays["components"] + arrays["mean"]

    feature_ids = manifest["feature_ids"]
    ground_truth = np.asarray(truth[shuffled, feature_ids].X, dtype=np.float64)
    np.testing.assert_allclose(reconstructed, ground_truth, atol=1e-6)


def test_pca_respects_gene_axis_for_raw_x_slot(tmp_path):
    """`raw.X` must be read against `raw/var`'s gene axis, never the
    current (possibly-filtered) `var` -- mirrors `cli.py`'s own established
    raw.X-vs-var distinction.
    """
    import anndata as ad
    import pandas as pd

    n_obs = 5
    obs = pd.DataFrame(
        {
            "cell_barcode": [f"CELL{i}" for i in range(n_obs)],
            "donor_id": ["D1"] * n_obs,
            "sample_id": ["S1"] * n_obs,
            "chip_id": ["C1"] * n_obs,
        }
    )
    var_raw = pd.DataFrame(index=[f"GENE{i}" for i in range(6)])
    a = ad.AnnData(
        X=np.random.default_rng(1).integers(0, 5, size=(n_obs, 6)).astype("float32"),
        obs=obs,
        var=var_raw,
    )
    a.raw = a
    a = a[:, :2].copy()  # current var now only has GENE0, GENE1
    path = tmp_path / "withraw.h5ad"
    a.write_h5ad(path)

    config = _base_config(tmp_path, str(path), counts_slot="raw.X")
    provider = E.PCAEmbeddingProvider()
    result = provider.fit_transform(
        config,
        _fake_config_path(tmp_path),
        out_dir=str(tmp_path / "pca_out"),
        n_components=2,
        random_seed=0,
        svd_solver="full",
    )
    assert result.status == "complete", result.reason
    assert set(result.feature_ids) == {f"GENE{i}" for i in range(6)}


def test_pca_rejects_ambiguous_feature_selection(tmp_path):
    h5ad = _make_tiny_h5ad(tmp_path, n_vars=8, n_obs=10)
    config = _base_config(tmp_path, h5ad)
    provider = E.PCAEmbeddingProvider()
    result = provider.fit_transform(
        config,
        _fake_config_path(tmp_path),
        out_dir=str(tmp_path / "pca_out"),
        n_components=2,
        hvg_top_n=4,
        feature_ids_override=("GENE0", "GENE1"),
    )
    assert result.status == "failed"
    assert "at most one" in result.reason


def test_pca_rejects_too_many_components(tmp_path):
    h5ad = _make_tiny_h5ad(tmp_path, n_vars=4, n_obs=6)
    config = _base_config(tmp_path, h5ad)
    provider = E.PCAEmbeddingProvider()
    result = provider.fit_transform(
        config, _fake_config_path(tmp_path), out_dir=str(tmp_path / "pca_out"), n_components=10
    )
    assert result.status == "failed"
    assert "exceeds min" in result.reason


def test_pca_rejects_nonfinite_values_in_selected_counts(tmp_path):
    import h5py

    h5ad = _make_tiny_h5ad(tmp_path, n_vars=4, n_obs=6)
    with h5py.File(h5ad, "a") as f:
        f["X"][0, 0] = float("nan")
    config = _base_config(tmp_path, h5ad)
    provider = E.PCAEmbeddingProvider()
    result = provider.fit_transform(
        config, _fake_config_path(tmp_path), out_dir=str(tmp_path / "pca_out"), n_components=2
    )
    assert result.status == "failed"
    assert "non-finite" in result.reason


def test_pca_rejects_zero_variance_gene_when_standardize_requested(tmp_path):
    import h5py

    h5ad = _make_tiny_h5ad(tmp_path, n_vars=4, n_obs=6)
    with h5py.File(h5ad, "a") as f:
        f["X"][:, 0] = 3.0  # constant column -- zero variance
    config = _base_config(tmp_path, h5ad)
    provider = E.PCAEmbeddingProvider()
    result = provider.fit_transform(
        config,
        _fake_config_path(tmp_path),
        out_dir=str(tmp_path / "pca_out"),
        n_components=2,
        standardize=True,
    )
    assert result.status == "failed"
    assert "zero variance" in result.reason


def test_pca_hvg_top_n_selects_the_highest_variance_genes_deterministically(tmp_path):
    """Label-free feature selection: the exact top-N-by-variance gene set,
    never a disease/outcome-label-based choice (this interface has no such
    parameter at all).
    """
    import h5py

    h5ad = _make_tiny_h5ad(tmp_path, n_vars=6, n_obs=20)
    # Make GENE0's column have far higher variance than the rest, by construction.
    with h5py.File(h5ad, "a") as f:
        rng = np.random.default_rng(2)
        f["X"][:, 0] = rng.normal(loc=0, scale=50, size=20).astype("float32")
    config = _base_config(tmp_path, h5ad)
    provider = E.PCAEmbeddingProvider()
    result = provider.fit_transform(
        config,
        _fake_config_path(tmp_path),
        out_dir=str(tmp_path / "pca_out"),
        n_components=1,
        hvg_top_n=1,
        svd_solver="full",
    )
    assert result.status == "complete", result.reason
    assert result.feature_ids == ("GENE0",)


def test_pca_refuses_an_oversized_dense_read(tmp_path, monkeypatch):
    """Lowers `_MAX_DENSE_ELEMENTS` far below this tiny fixture's real size
    and confirms the explicit scope bound actually refuses it -- proves the
    bound is wired into the real read path, not merely declared.
    """
    h5ad = _make_tiny_h5ad(tmp_path, n_vars=8, n_obs=10)
    monkeypatch.setattr(E, "_MAX_DENSE_ELEMENTS", 5)
    config = _base_config(tmp_path, h5ad)
    provider = E.PCAEmbeddingProvider()
    result = provider.fit_transform(
        config, _fake_config_path(tmp_path), out_dir=str(tmp_path / "pca_out"), n_components=2
    )
    assert result.status == "failed"
    assert "element" in result.reason


def test_read_embedding_result_rejects_minimal_complete_payload_missing_everything(tmp_path):
    import json

    path = tmp_path / "bad.json"
    path.write_text(json.dumps({"schema": E.EMBEDDING_RESULT_SCHEMA, "status": "complete"}))
    with pytest.raises(SystemExit, match="requires a non-empty"):
        E.read_embedding_result(str(path))


def test_load_embedding_matrices_rejects_a_tampered_npz(tmp_path):
    """The manifest's own `embedding_sha256` must be re-checked against the
    npz file's ACTUAL current hash on read -- a tampered/stale npz file
    must be caught, not only a tampered manifest.
    """
    h5ad = _make_tiny_h5ad(tmp_path, n_vars=6, n_obs=10)
    config = _base_config(tmp_path, h5ad)
    provider = E.PCAEmbeddingProvider()
    result = provider.fit_transform(
        config, _fake_config_path(tmp_path), out_dir=str(tmp_path / "pca_out"), n_components=2
    )
    assert result.status == "complete"
    with open(result.embedding_path, "ab") as fh:
        fh.write(b"\x00")  # tamper the npz AFTER it was hashed and written
    manifest = E.read_embedding_result(str(tmp_path / "pca_out" / "EMBEDDING_RESULT.json"))
    with pytest.raises(SystemExit, match="does not match"):
        E.load_embedding_matrices(manifest)


# ===========================================================================
# F2 (review:20261013_210131) -- duplicate gene IDs must be
# refused with a clean failed record, never silently resolved by keeping
# whichever column happens to come last in the ID->position map.
# ===========================================================================


def _make_h5ad_with_gene_ids(tmp_path, X, gene_ids, *, name="dup.h5ad"):
    import anndata as ad
    import pandas as pd

    n_obs = len(X)
    obs = pd.DataFrame(
        {
            "cell_barcode": [f"CELL{i}" for i in range(n_obs)],
            "donor_id": [f"D{i % 2}" for i in range(n_obs)],
            "sample_id": [f"S{i}" for i in range(n_obs)],
            "chip_id": ["C1"] * n_obs,
        },
        index=[f"idx{i}" for i in range(n_obs)],
    )
    var = pd.DataFrame(index=list(gene_ids))
    a = ad.AnnData(X=np.asarray(X, dtype="float32"), obs=obs, var=var)
    path = tmp_path / name
    a.write_h5ad(path)
    return str(path)


def test_pca_rejects_duplicate_gene_ids_before_feature_override_selection(tmp_path):
    """Per F2: a source gene axis with `DUP, DUP, G2, G3` (two columns with
    the same ID but DIFFERENT values) must be refused outright, never
    silently collapsed onto whichever column `{g: i for i, g in
    enumerate(gene_ids)}` happened to keep last.
    """
    X = [
        [1.0, 8.0, 2.0, 3.0],
        [1.0, 8.0, 2.0, 3.0],
        [1.0, 8.0, 2.0, 3.0],
        [1.0, 8.0, 2.0, 3.0],
        [1.0, 8.0, 2.0, 3.0],
        [1.0, 8.0, 2.0, 3.0],
    ]
    h5ad = _make_h5ad_with_gene_ids(tmp_path, X, ["DUP", "DUP", "G2", "G3"])
    config = _base_config(tmp_path, h5ad)
    provider = E.PCAEmbeddingProvider()
    result = provider.fit_transform(
        config,
        _fake_config_path(tmp_path),
        out_dir=str(tmp_path / "pca_out"),
        n_components=2,
        feature_ids_override=("DUP", "G2"),
        svd_solver="full",
    )
    assert result.status == "failed"
    assert "duplicate gene ID" in result.reason
    assert "DUP" in result.reason


def test_pca_rejects_duplicate_gene_ids_even_without_feature_override(tmp_path):
    """The gene-axis validation runs BEFORE any feature-selection branch --
    including the plain "use all genes" path, which never even builds a
    `gpos` map itself.
    """
    X = [[1.0, 8.0, 2.0, 3.0]] * 6
    h5ad = _make_h5ad_with_gene_ids(tmp_path, X, ["DUP", "DUP", "G2", "G3"])
    config = _base_config(tmp_path, h5ad)
    provider = E.PCAEmbeddingProvider()
    result = provider.fit_transform(
        config,
        _fake_config_path(tmp_path),
        out_dir=str(tmp_path / "pca_out"),
        n_components=2,
        svd_solver="full",
    )
    assert result.status == "failed"
    assert "duplicate gene ID" in result.reason


def test_pca_still_accepts_ordinary_unique_gene_ids_after_the_duplicate_check(tmp_path):
    """The new F2 guard must not regress the normal, non-duplicated case --
    an explicit non-regression companion to the two tests above."""
    h5ad = _make_tiny_h5ad(tmp_path, n_vars=6, n_obs=10)
    config = _base_config(tmp_path, h5ad)
    provider = E.PCAEmbeddingProvider()
    result = provider.fit_transform(
        config,
        _fake_config_path(tmp_path),
        out_dir=str(tmp_path / "pca_out"),
        n_components=2,
        feature_ids_override=("GENE0", "GENE1"),
        svd_solver="full",
    )
    assert result.status == "complete", result.reason
    assert result.feature_ids == ("GENE0", "GENE1")


def test_mutation_kill_pca_duplicate_gene_id_guard(monkeypatch, tmp_path):
    """Disables `_validate_gene_axis_values` (as if the F2 guard did not
    exist) and confirms the duplicate-gene-ID scenario above would then
    silently pick whichever DUP column the ID->position map happens to
    keep -- reproducing the original F2 defect exactly, including the
    wrong numeric result it produced.
    """
    import anndata as ad

    monkeypatch.setattr(E, "_validate_gene_axis_values", lambda gene_ids, path: [])
    # Real row-to-row variance (unlike the all-constant fixtures above,
    # which fail before PCA ever runs in those tests) -- here PCA must
    # genuinely succeed so the SILENT wrong-column selection is visible.
    X = [
        [1.0, 8.0, 2.0, 3.0],
        [2.0, 9.0, 1.0, 4.0],
        [3.0, 7.0, 5.0, 2.0],
        [4.0, 6.0, 3.0, 6.0],
        [5.0, 5.0, 4.0, 1.0],
        [6.0, 4.0, 6.0, 5.0],
    ]
    h5ad = _make_h5ad_with_gene_ids(tmp_path, X, ["DUP", "DUP", "G2", "G3"])
    config = _base_config(tmp_path, h5ad)
    provider = E.PCAEmbeddingProvider()
    out_dir = str(tmp_path / "pca_out")
    result = provider.fit_transform(
        config,
        _fake_config_path(tmp_path),
        out_dir=out_dir,
        n_components=2,
        feature_ids_override=("DUP", "G2"),
        svd_solver="full",
    )
    assert result.status == "complete"  # with the mutant: silently "succeeds"
    manifest = E.read_embedding_result(os.path.join(out_dir, "EMBEDDING_RESULT.json"))
    arrays = E.load_embedding_matrices(manifest)
    reconstructed = arrays["embedding"] @ arrays["components"] + arrays["mean"]
    truth = ad.read_h5ad(h5ad)
    # with the mutant: `{g: i for i, g in enumerate(gene_ids)}` keeps only
    # the LAST "DUP" occurrence (var index 1), silently dropping the FIRST
    # "DUP" column (var index 0) with no indication anything was lost --
    # the F2 fix refuses this outright instead.
    np.testing.assert_allclose(reconstructed[:, 0], np.asarray(truth.X)[:, 1], atol=1e-6)
    assert np.max(np.abs(reconstructed[:, 0] - np.asarray(truth.X)[:, 0])) > 1.0


# ===========================================================================
# F3 (review:20261013_210131) -- a zero-total-variance input must
# leave a clean failed manifest, never an orphaned npz with no
# EMBEDDING_RESULT.json beside it.
# ===========================================================================


def test_pca_zero_variance_input_leaves_a_clean_failed_manifest(tmp_path):
    h5ad = _make_tiny_h5ad(tmp_path, n_vars=4, n_obs=8, low=0, high=1)  # all-zero
    config = _base_config(tmp_path, h5ad)
    provider = E.PCAEmbeddingProvider()
    out_dir = str(tmp_path / "pca_out")
    result = provider.fit_transform(
        config,
        _fake_config_path(tmp_path),
        out_dir=out_dir,
        n_components=2,
        standardize=False,
        svd_solver="full",
    )
    assert result.status == "failed"
    assert "zero total variance" in result.reason
    assert os.path.exists(os.path.join(out_dir, "EMBEDDING_RESULT.json"))
    assert not os.path.exists(os.path.join(out_dir, "embedding_matrices.npz"))
    manifest = E.read_embedding_result(os.path.join(out_dir, "EMBEDDING_RESULT.json"))
    assert manifest["status"] == "failed"


def test_mutation_kill_pca_zero_total_variance_guard(monkeypatch, tmp_path):
    """Disables `_zero_total_variance_reason` and confirms the all-zero
    input above is STILL caught (by the later finiteness safety net), but
    with a DIFFERENT, more generic reason -- proving the up-front guard is
    what actually produces the specific, actionable message, not just a
    name coincidentally present in the final reason string.
    """
    monkeypatch.setattr(E, "_zero_total_variance_reason", lambda X: None)
    h5ad = _make_tiny_h5ad(tmp_path, n_vars=4, n_obs=8, low=0, high=1)
    config = _base_config(tmp_path, h5ad)
    provider = E.PCAEmbeddingProvider()
    result = provider.fit_transform(
        config,
        _fake_config_path(tmp_path),
        out_dir=str(tmp_path / "pca_out"),
        n_components=2,
        standardize=False,
        svd_solver="full",
    )
    assert result.status == "failed"
    assert "zero total variance" not in result.reason
    assert "non-finite" in result.reason


def test_mutation_kill_pca_finiteness_check_covers_explained_variance_ratio(monkeypatch, tmp_path):
    """The ORIGINAL F3 bug: `_pca_output_finiteness_problem` narrowed to
    only check embedding/components/mean (as the pre-fix code did) lets a
    NaN `explained_variance_ratio_` reach `_write_embedding_result`, which
    writes `embedding_matrices.npz` FIRST and only then raises on the JSON
    manifest's own validation -- leaving an orphaned npz with no
    EMBEDDING_RESULT.json beside it. Bypasses the up-front zero-variance
    guard too, to isolate this specific check.
    """
    monkeypatch.setattr(E, "_zero_total_variance_reason", lambda X: None)

    def narrow_check(embedding, components, mean, explained_variance_ratio, singular_values):
        if not (
            np.all(np.isfinite(embedding))
            and np.all(np.isfinite(components))
            and np.all(np.isfinite(mean))
        ):
            return "PCA produced non-finite value(s) -- refusing to persist it"
        return None  # the ORIGINAL bug: explained_variance_ratio/singular_values never checked

    monkeypatch.setattr(E, "_pca_output_finiteness_problem", narrow_check)
    h5ad = _make_tiny_h5ad(tmp_path, n_vars=4, n_obs=8, low=0, high=1)
    config = _base_config(tmp_path, h5ad)
    provider = E.PCAEmbeddingProvider()
    out_dir = str(tmp_path / "pca_out")
    with pytest.raises(ValueError, match="refusing to write an invalid embedding result"):
        provider.fit_transform(
            config,
            _fake_config_path(tmp_path),
            out_dir=out_dir,
            n_components=2,
            standardize=False,
            svd_solver="full",
        )
    # the ORIGINAL bug symptom: npz got written, but JSON validation then
    # raised, leaving an orphaned npz with no EMBEDDING_RESULT.json.
    assert os.path.exists(os.path.join(out_dir, "embedding_matrices.npz"))
    assert not os.path.exists(os.path.join(out_dir, "EMBEDDING_RESULT.json"))


def test_pca_unexpected_backend_exception_is_recorded_then_reraised(tmp_path, monkeypatch):
    """Per F3: a genuinely unexpected backend exception (not a graceful,
    anticipated 'failed' condition) is recorded in a failed manifest FIRST,
    then re-raised -- never silently swallowed into a soft failure that
    could go unnoticed.
    """
    import sklearn.decomposition

    def boom(self, X, y=None):
        raise RuntimeError("synthetic backend crash")

    monkeypatch.setattr(sklearn.decomposition.PCA, "fit_transform", boom)
    h5ad = _make_tiny_h5ad(tmp_path, n_vars=4, n_obs=6)
    config = _base_config(tmp_path, h5ad)
    provider = E.PCAEmbeddingProvider()
    out_dir = str(tmp_path / "pca_out")
    with pytest.raises(RuntimeError, match="synthetic backend crash"):
        provider.fit_transform(
            config, _fake_config_path(tmp_path), out_dir=out_dir, n_components=2, svd_solver="full"
        )
    manifest = E.read_embedding_result(os.path.join(out_dir, "EMBEDDING_RESULT.json"))
    assert manifest["status"] == "failed"
    assert "synthetic backend crash" in manifest["reason"]


def test_mutation_kill_validate_embedding_result_skipped_on_write(monkeypatch, tmp_path):
    """Disables `validate_embedding_result` inside `_write_embedding_result`
    and confirms a deliberately invalid 'complete' result would then be
    written WITHOUT error -- proving the real validation call, not some
    other check, is what a correct write path depends on.
    """
    monkeypatch.setattr(E, "validate_embedding_result", lambda payload: [])
    bad_result = E.EmbeddingResult(status="complete", backend="x", backend_version="1")
    path = tmp_path / "out"
    os.makedirs(path)
    # With the mutant, this does NOT raise despite being structurally invalid
    # (empty cell_ids/feature_ids, no dataset_identity, n_components=0, etc.).
    written = E._write_embedding_result(
        str(path),
        bad_result,
        embedding=np.zeros((0, 1)),
        components=np.zeros((1, 0)),
        mean=np.zeros((0,)),
        singular_values=np.zeros((1,)),
    )
    assert (path / "EMBEDDING_RESULT.json").exists()
    assert written.status == "complete"

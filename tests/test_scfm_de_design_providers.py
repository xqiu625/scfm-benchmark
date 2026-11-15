"""Tests for `de_design_providers.py` -- the "B" increment (pseudobulk ->
DE input/design-matrix PREPARATION, per `review:20261014_151747`
step 2). NO statistical fitting happens here; every fixture is a tiny,
synthetic, local h5ad with HAND-COMPUTABLE sums and a hand-computable
design matrix/rank/residual-degrees-of-freedom -- no network, no real
data. `ExpressionDERankingProvider` (`gene_ranking_providers.py`) still
has no real backend; this module's own output is "DE input/design
prepared" only.
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

import de_design_providers as DD  # noqa: E402
import gene_ranking_providers as GR  # noqa: E402
import pseudobulk_providers as PB  # noqa: E402
from dataset_schema import NO_BATCH_SENTINEL  # noqa: E402


def _base_config(matrix_path, **overrides):
    config = {
        "dataset_id": "test_de_design_pilot",
        "species": "human",
        "tissue_or_target_cell_type": "monocyte",
        "subject_id_column": "donor_id",
        "sample_id_column": "sample_id",
        "biological_replicate_unit": "donor_id",
        "technical_batch_column": NO_BATCH_SENTINEL,
        "expression_matrix_path": matrix_path,
        "counts_semantics": "raw_counts",
        "counts_slot": "X",
        "gene_id_namespace": "gene_symbol",
        "cell_id_column": "cell_barcode",
        "qc_status": "provisional",
        "qc_version_ref": "none (provisional)",
    }
    config.update(overrides)
    return config


def _make_h5ad(tmp_path, *, X, cell_ids, sample_ids, donor_ids, conditions, name="de_design.h5ad"):
    import anndata as ad
    import pandas as pd

    obs = {
        "cell_barcode": cell_ids,
        "donor_id": donor_ids,
        "sample_id": sample_ids,
        "condition": conditions,
    }
    obs_df = pd.DataFrame(obs, index=[f"idx{i}" for i in range(len(cell_ids))])
    var = pd.DataFrame(index=[f"GENE{i}" for i in range(X.shape[1])])
    a = ad.AnnData(X=np.asarray(X, dtype="float32"), obs=obs_df, var=var)
    path = tmp_path / name
    a.write_h5ad(path)
    return str(path)


def _comparison(**overrides):
    base = dict(
        group_column="condition",
        reference_level="healthy",
        comparison_level="disease",
        direction_definition="disease vs. healthy, positive = higher in disease",
    )
    base.update(overrides)
    return GR.ComparisonSpec(**base)


def _sample_structure(**overrides):
    base = dict(
        donor_or_animal_id_column="donor_id",
        sample_id_column="sample_id",
        repeated_sampling_or_pairing="none_single_sample_per_donor",
        group_column="condition",
        design_type="independent_donors_unpaired",
    )
    base.update(overrides)
    return GR.SampleStructure(**base)


def _aggregate(tmp_path, out_dir_name, *, X, cell_ids, sample_ids, donor_ids, conditions):
    """Returns `(manifest, counts, pseudobulk_result_path)` -- the real
    `PSEUDOBULK_RESULT.json` path is what `DEDesignPreparationProvider.
    prepare` now takes (per `review:20261014_203915` F1), with
    `manifest`/`counts` returned alongside for tests that assert against
    the real producer's own output.
    """
    h5ad = _make_h5ad(
        tmp_path,
        X=X,
        cell_ids=cell_ids,
        sample_ids=sample_ids,
        donor_ids=donor_ids,
        conditions=conditions,
    )
    config = _base_config(h5ad)
    config_path = str(tmp_path / "dataset_config.json")
    provider = PB.PseudobulkAggregationProvider()
    result = provider.aggregate(
        config,
        config_path,
        out_dir=str(tmp_path / out_dir_name),
        cell_ids_subset=tuple(cell_ids),
        cell_type_or_state_label="monocyte",
        condition_column="condition",
    )
    assert result.status == "complete", result.reason
    path = str(tmp_path / out_dir_name / "PSEUDOBULK_RESULT.json")
    manifest = PB.read_pseudobulk_result(path)
    counts = PB.load_pseudobulk_counts(manifest)
    return manifest, counts, path


#: 4 INDEPENDENT donors, one sample each (1 cell/sample -- library size IS
#: the cell's own counts): D0/D1 healthy, D2/D3 disease.
_UNPAIRED_X = np.array([[1, 2], [3, 4], [5, 6], [7, 8]], dtype="float32")
_UNPAIRED_CELL_IDS = ["CELL0", "CELL1", "CELL2", "CELL3"]
_UNPAIRED_SAMPLE = ["S0", "S1", "S2", "S3"]
_UNPAIRED_DONOR = ["D0", "D1", "D2", "D3"]
_UNPAIRED_CONDITION = ["healthy", "healthy", "disease", "disease"]


def _unpaired_fixture(tmp_path):
    return _aggregate(
        tmp_path,
        "pb_unpaired",
        X=_UNPAIRED_X,
        cell_ids=_UNPAIRED_CELL_IDS,
        sample_ids=_UNPAIRED_SAMPLE,
        donor_ids=_UNPAIRED_DONOR,
        conditions=_UNPAIRED_CONDITION,
    )


#: 2 donors, each contributing ONE healthy + ONE disease sample (fully
#: paired): D0: S0(healthy)+S1(disease), D1: S2(healthy)+S3(disease).
_PAIRED_X = np.array([[1, 2], [3, 4], [5, 6], [7, 8]], dtype="float32")
_PAIRED_CELL_IDS = ["CELL0", "CELL1", "CELL2", "CELL3"]
_PAIRED_SAMPLE = ["S0", "S1", "S2", "S3"]
_PAIRED_DONOR = ["D0", "D0", "D1", "D1"]
_PAIRED_CONDITION = ["healthy", "disease", "healthy", "disease"]


def _paired_fixture(tmp_path):
    return _aggregate(
        tmp_path,
        "pb_paired",
        X=_PAIRED_X,
        cell_ids=_PAIRED_CELL_IDS,
        sample_ids=_PAIRED_SAMPLE,
        donor_ids=_PAIRED_DONOR,
        conditions=_PAIRED_CONDITION,
    )


# ===========================================================================
# Happy path: hand-computable design matrix, rank, residual df, contrast.
# ===========================================================================


def test_prepare_unpaired_design_hand_computable(tmp_path):
    manifest, counts, pb_path = _unpaired_fixture(tmp_path)
    comparison = _comparison()
    sample_structure = _sample_structure(design_type="independent_donors_unpaired")
    provider = DD.DEDesignPreparationProvider()
    result = provider.prepare(
        pseudobulk_result_path=pb_path,
        comparison=comparison,
        sample_structure=sample_structure,
        out_dir=str(tmp_path / "prep_unpaired"),
    )

    assert result.status == "prepared", result.reason
    # Row order: donor id ascending (each donor contributes exactly one row).
    assert result.sample_ids == ("S0", "S1", "S2", "S3")
    assert result.donor_ids == ("D0", "D1", "D2", "D3")
    assert result.group_values == ("healthy", "healthy", "disease", "disease")
    assert result.column_names == ("intercept", "condition[disease]")
    assert result.design_matrix == ((1.0, 0.0), (1.0, 0.0), (1.0, 1.0), (1.0, 1.0))
    assert result.rank == 2
    assert result.residual_degrees_of_freedom == 2  # n_samples(4) - rank(2)
    assert result.n_samples == 4
    assert result.n_distinct_donors == 4
    assert result.contrast_vector == (0.0, 1.0)
    assert "disease" in result.contrast_direction and "healthy" in result.contrast_direction
    assert result.pseudobulk_identity["aggregation_sha256"] == manifest["aggregation_sha256"]
    assert result.pseudobulk_identity["pseudobulk_manifest_path"] == os.path.abspath(pb_path)
    assert result.dataset_identity["species"] == "human"

    # Real round trip.
    path = tmp_path / "prep_unpaired" / "PREPARED_DE_DESIGN.json"
    reloaded = DD.read_prepared_de_design(str(path))
    assert reloaded["status"] == "prepared"
    prepared_counts = DD.load_prepared_de_counts(reloaded)
    # S0=[1,2], S1=[3,4], S2=[5,6], S3=[7,8] (1 cell/sample => counts == raw cell values).
    np.testing.assert_array_equal(prepared_counts, [[1.0, 2.0], [3.0, 4.0], [5.0, 6.0], [7.0, 8.0]])


def test_prepare_paired_design_hand_computable(tmp_path):
    manifest, counts, pb_path = _paired_fixture(tmp_path)
    comparison = _comparison()
    sample_structure = _sample_structure(
        design_type="paired_or_repeated_within_donor",
        repeated_sampling_or_pairing="paired:pre_post",
    )
    provider = DD.DEDesignPreparationProvider()
    result = provider.prepare(
        pseudobulk_result_path=pb_path,
        comparison=comparison,
        sample_structure=sample_structure,
        out_dir=str(tmp_path / "prep_paired"),
    )

    assert result.status == "prepared", result.reason
    # Row order: donor (D0, D1), group (healthy before disease) within each donor.
    assert result.sample_ids == ("S0", "S1", "S2", "S3")
    assert result.donor_ids == ("D0", "D0", "D1", "D1")
    assert result.group_values == ("healthy", "disease", "healthy", "disease")
    assert result.n_distinct_donors == 2
    assert len(result.column_names) == 3  # intercept, donor[D1], condition[disease]
    assert result.design_matrix == (
        (1.0, 0.0, 0.0),
        (1.0, 0.0, 1.0),
        (1.0, 1.0, 0.0),
        (1.0, 1.0, 1.0),
    )
    assert result.rank == 3
    assert result.residual_degrees_of_freedom == 1  # exactly at the >=1 boundary
    assert result.contrast_vector == (0.0, 0.0, 1.0)

    path = tmp_path / "prep_paired" / "PREPARED_DE_DESIGN.json"
    reloaded = DD.read_prepared_de_design(str(path))
    prepared_counts = DD.load_prepared_de_counts(reloaded)
    np.testing.assert_array_equal(prepared_counts, [[1.0, 2.0], [3.0, 4.0], [5.0, 6.0], [7.0, 8.0]])


def test_prepare_aligns_by_sample_id_not_by_pseudobulk_row_order(tmp_path):
    """The pseudobulk manifest's OWN `sample_rows` order need not match
    the design's row order (donor-sorted) -- `prepare` must select/
    reorder by `sample_id` LOOKUP, not by assuming position. Per
    `review:20261014_203915` F1, `prepare` no longer accepts an
    in-memory manifest a test could tamper with to produce this -- so
    this builds a REAL aggregation through the real
    aggregate -> read -> prepare path. `PseudobulkAggregationProvider`
    itself always orders `sample_rows` alphabetically BY `sample_id`
    (`pseudobulk_providers.py`: `for sid in sorted(by_sample)`), so
    `sample_id`s are chosen that sort in the OPPOSITE order from their
    donors (D0/D1/D2/D3 take sample_ids Z0/Y1/X2/W3) -- the pseudobulk's
    own row order (W3, X2, Y1, Z0) is then the exact REVERSE of the
    design's donor-sorted order (Z0, Y1, X2, W3).
    """
    X = np.array([[1, 2], [3, 4], [5, 6], [7, 8]], dtype="float32")
    cell_ids = ["CELL0", "CELL1", "CELL2", "CELL3"]
    manifest, counts, pb_path = _aggregate(
        tmp_path,
        "pb_reversed_order",
        X=X,
        cell_ids=cell_ids,
        sample_ids=["Z0", "Y1", "X2", "W3"],
        donor_ids=["D0", "D1", "D2", "D3"],
        conditions=["healthy", "healthy", "disease", "disease"],
    )
    assert [r["sample_id"] for r in manifest["sample_rows"]] == ["W3", "X2", "Y1", "Z0"]

    comparison = _comparison()
    sample_structure = _sample_structure(design_type="independent_donors_unpaired")
    provider = DD.DEDesignPreparationProvider()
    result = provider.prepare(
        pseudobulk_result_path=pb_path,
        comparison=comparison,
        sample_structure=sample_structure,
        out_dir=str(tmp_path / "prep_reversed"),
    )
    assert result.status == "prepared", result.reason
    assert result.sample_ids == ("Z0", "Y1", "X2", "W3")  # re-sorted to donor order

    path = tmp_path / "prep_reversed" / "PREPARED_DE_DESIGN.json"
    reloaded = DD.read_prepared_de_design(str(path))
    prepared_counts = DD.load_prepared_de_counts(reloaded)
    # Z0=[1,2] (CELL0), Y1=[3,4] (CELL1), X2=[5,6] (CELL2), W3=[7,8] (CELL3).
    np.testing.assert_array_equal(prepared_counts, [[1.0, 2.0], [3.0, 4.0], [5.0, 6.0], [7.0, 8.0]])


# ===========================================================================
# Explicit refusals: never silently pooled, never a partial/degenerate matrix.
# ===========================================================================


def test_prepare_rejects_declared_covariate_columns(tmp_path):
    manifest, counts, pb_path = _unpaired_fixture(tmp_path)
    comparison = _comparison()
    sample_structure = _sample_structure(
        design_type="independent_donors_unpaired", covariate_columns=("age",)
    )
    provider = DD.DEDesignPreparationProvider()
    result = provider.prepare(
        pseudobulk_result_path=pb_path,
        comparison=comparison,
        sample_structure=sample_structure,
        out_dir=str(tmp_path / "prep_covariate"),
    )
    assert result.status == "failed"
    assert "covariate_columns" in result.reason
    assert result.sample_ids == ()  # no partial/fabricated design on failure


def test_prepare_rejects_a_straddling_donor_under_unpaired_design(tmp_path):
    """A donor contributing samples to BOTH groups contradicts
    'independent_donors_unpaired' -- must be refused, never silently
    treated as two independent replicates.
    """
    X = np.array([[1, 2], [3, 4], [5, 6]], dtype="float32")
    cell_ids = ["CELL0", "CELL1", "CELL2"]
    manifest, counts, pb_path = _aggregate(
        tmp_path,
        "pb_straddle",
        X=X,
        cell_ids=cell_ids,
        sample_ids=["S0", "S1", "S2"],
        donor_ids=["D0", "D0", "D1"],
        conditions=["healthy", "disease", "disease"],
    )
    comparison = _comparison()
    sample_structure = _sample_structure(design_type="independent_donors_unpaired")
    provider = DD.DEDesignPreparationProvider()
    result = provider.prepare(
        pseudobulk_result_path=pb_path,
        comparison=comparison,
        sample_structure=sample_structure,
        out_dir=str(tmp_path / "prep_straddle"),
    )
    assert result.status == "failed"
    assert "BOTH compared groups" in result.reason


def test_prepare_rejects_incomplete_pairing_under_paired_design(tmp_path):
    """`paired_or_repeated_within_donor` supports only the FULLY-paired
    case -- a donor missing one level must be refused, never dropped.
    """
    X = np.array([[1, 2], [3, 4], [5, 6]], dtype="float32")
    cell_ids = ["CELL0", "CELL1", "CELL2"]
    manifest, counts, pb_path = _aggregate(
        tmp_path,
        "pb_unpaired_donor",
        X=X,
        cell_ids=cell_ids,
        sample_ids=["S0", "S1", "S2"],
        donor_ids=["D0", "D0", "D1"],
        conditions=["healthy", "disease", "healthy"],
    )
    comparison = _comparison()
    sample_structure = _sample_structure(
        design_type="paired_or_repeated_within_donor",
        repeated_sampling_or_pairing="paired:pre_post",
    )
    provider = DD.DEDesignPreparationProvider()
    result = provider.prepare(
        pseudobulk_result_path=pb_path,
        comparison=comparison,
        sample_structure=sample_structure,
        out_dir=str(tmp_path / "prep_unpaired_donor"),
    )
    assert result.status == "failed"
    assert "FULLY-paired" in result.reason


def test_prepare_rejects_a_repeated_sample_within_the_same_donor_group(tmp_path):
    """Two samples for the SAME donor within the SAME compared group (e.g.
    an extra visit) is unsupported this round -- must never be silently
    pooled into one row or treated as two independent donors.
    """
    X = np.array([[1, 2], [3, 4], [5, 6]], dtype="float32")
    cell_ids = ["CELL0", "CELL1", "CELL2"]
    manifest, counts, pb_path = _aggregate(
        tmp_path,
        "pb_repeat",
        X=X,
        cell_ids=cell_ids,
        sample_ids=["S0", "S1", "S2"],
        donor_ids=["D0", "D0", "D1"],
        conditions=["healthy", "healthy", "disease"],
    )
    comparison = _comparison()
    sample_structure = _sample_structure(design_type="independent_donors_unpaired")
    provider = DD.DEDesignPreparationProvider()
    result = provider.prepare(
        pseudobulk_result_path=pb_path,
        comparison=comparison,
        sample_structure=sample_structure,
        out_dir=str(tmp_path / "prep_repeat"),
    )
    assert result.status == "failed"
    assert "more than one pseudobulk sample" in result.reason


def test_prepare_rejects_a_missing_comparison_level(tmp_path):
    X = np.array([[1, 2], [3, 4]], dtype="float32")
    cell_ids = ["CELL0", "CELL1"]
    manifest, counts, pb_path = _aggregate(
        tmp_path,
        "pb_missing_group",
        X=X,
        cell_ids=cell_ids,
        sample_ids=["S0", "S1"],
        donor_ids=["D0", "D1"],
        conditions=["healthy", "healthy"],
    )
    comparison = _comparison()
    sample_structure = _sample_structure(design_type="independent_donors_unpaired")
    provider = DD.DEDesignPreparationProvider()
    result = provider.prepare(
        pseudobulk_result_path=pb_path,
        comparison=comparison,
        sample_structure=sample_structure,
        out_dir=str(tmp_path / "prep_missing_group"),
    )
    assert result.status == "failed"
    assert "zero pseudobulk samples" in result.reason


def test_prepare_rejects_a_failed_pseudobulk_aggregation(tmp_path):
    """A failed/incomplete pseudobulk aggregation must never enter DE
    design preparation at all.
    """
    h5ad = _make_h5ad(
        tmp_path,
        X=np.array([[1.0]], dtype="float32"),
        cell_ids=["CELL0"],
        sample_ids=["S0"],
        donor_ids=["D0"],
        conditions=["healthy"],
    )
    config = _base_config(h5ad)
    config_path = str(tmp_path / "dataset_config.json")
    pb_provider = PB.PseudobulkAggregationProvider()
    failed = pb_provider.aggregate(
        config,
        config_path,
        out_dir=str(tmp_path / "pb_failed"),
        cell_ids_subset=(),
        cell_type_or_state_label="monocyte",
    )  # empty subset -> fails
    assert failed.status == "failed"
    pb_path = str(tmp_path / "pb_failed" / "PSEUDOBULK_RESULT.json")

    comparison = _comparison()
    sample_structure = _sample_structure(design_type="independent_donors_unpaired")
    provider = DD.DEDesignPreparationProvider()
    result = provider.prepare(
        pseudobulk_result_path=pb_path,
        comparison=comparison,
        sample_structure=sample_structure,
        out_dir=str(tmp_path / "prep_from_failed"),
    )
    assert result.status == "failed"
    assert "never enter DE design preparation" in result.reason


def test_prepare_rejects_a_design_with_no_residual_degrees_of_freedom(tmp_path):
    """A single paired donor (2 samples, 1 healthy + 1 disease) is FULL
    RANK (2 columns, 2 samples) but leaves ZERO residual degrees of
    freedom -- full rank is never itself treated as adequate power.
    """
    X = np.array([[1, 2], [3, 4]], dtype="float32")
    cell_ids = ["CELL0", "CELL1"]
    manifest, counts, pb_path = _aggregate(
        tmp_path,
        "pb_one_donor_paired",
        X=X,
        cell_ids=cell_ids,
        sample_ids=["S0", "S1"],
        donor_ids=["D0", "D0"],
        conditions=["healthy", "disease"],
    )
    comparison = _comparison()
    sample_structure = _sample_structure(
        design_type="paired_or_repeated_within_donor",
        repeated_sampling_or_pairing="paired:pre_post",
    )
    provider = DD.DEDesignPreparationProvider()
    result = provider.prepare(
        pseudobulk_result_path=pb_path,
        comparison=comparison,
        sample_structure=sample_structure,
        out_dir=str(tmp_path / "prep_zero_df"),
    )
    assert result.status == "failed"
    assert "residual_degrees_of_freedom=0" in result.reason


def test_prepare_rejects_an_unsupported_design_type_even_if_gene_ranking_would_accept_it(
    monkeypatch,
):
    """Per the module docstring: `_SUPPORTED_DESIGN_TYPES` is declared
    independently of `gene_ranking_providers._VALID_DESIGN_TYPES` on
    purpose. Simulates a hypothetical design type gene_ranking_providers
    already accepts (so `SampleStructure.validate()` would not block it)
    but this module does not yet implement, by shrinking this module's own
    supported set and confirming the real entry point still refuses it.
    """
    monkeypatch.setattr(DD, "_SUPPORTED_DESIGN_TYPES", ("paired_or_repeated_within_donor",))
    comparison = _comparison()
    sample_structure = _sample_structure(design_type="independent_donors_unpaired")
    sample_rows = [
        {"sample_id": "S0", "donor_or_animal_id": "D0", "condition_value": "healthy"},
        {"sample_id": "S1", "donor_or_animal_id": "D1", "condition_value": "disease"},
    ]
    problems, fields = DD._prepare_de_design_core(sample_rows, comparison, sample_structure)
    assert fields == {}
    assert any("not in" in p for p in problems)
    assert any("independent_donors_unpaired" in p for p in problems)


# ===========================================================================
# F1 (`review:20261014_203915`): the public `prepare` entry point
# must actually call the real pseudobulk reader/loader, never trust a
# caller-supplied counts array independently of the manifest.
# ===========================================================================


def test_prepare_signature_has_no_caller_supplied_counts_parameter():
    """Per F1: the fix is that `prepare` can no longer be called with an
    independently-supplied `counts` array at all -- a static guard against
    the fix being silently reverted (a caller-supplied `counts` parameter
    disjoint from the manifest is exactly the surface the review's probe
    exploited: same shape, reversed/negated content, still accepted).
    """
    import inspect

    sig = inspect.signature(DD.DEDesignPreparationProvider.prepare)
    assert "counts" not in sig.parameters
    assert "pseudobulk_manifest" not in sig.parameters
    assert "pseudobulk_result_path" in sig.parameters


def test_prepare_rejects_a_tampered_on_disk_counts_npz(tmp_path):
    """Per F1's reproduction (negate the counts array, same shape): since
    `prepare` now loads counts itself via the real, hash-bound
    `pseudobulk_providers.load_pseudobulk_counts`, tampering the npz's
    actual bytes on disk is caught by that loader's own file-hash check --
    there is no longer an in-memory counts argument to tamper instead.
    """
    manifest, counts, pb_path = _unpaired_fixture(tmp_path)
    npz_path = manifest["aggregation_path"]
    with np.load(npz_path, allow_pickle=False) as npz:
        data = {k: npz[k] for k in npz.files}
    data["counts"] = -data["counts"]  # same shape, wrong content
    with open(npz_path, "wb") as fh:
        np.savez(fh, **data)  # same file path, now-mismatched hash

    comparison = _comparison()
    sample_structure = _sample_structure(design_type="independent_donors_unpaired")
    provider = DD.DEDesignPreparationProvider()
    with pytest.raises(SystemExit, match="aggregation_sha256"):
        provider.prepare(
            pseudobulk_result_path=pb_path,
            comparison=comparison,
            sample_structure=sample_structure,
            out_dir=str(tmp_path / "prep_tampered_npz"),
        )


def test_prepare_rejects_a_reordered_sample_rows_manifest(tmp_path):
    """Per F1's reproduction (reverse `sample_rows` only): reversing the
    on-disk manifest's `sample_rows` without touching the npz changes
    `identity_sha256`'s recompute, so the real `read_pseudobulk_result`
    `prepare` now calls refuses it before any design is built.
    """
    import json

    manifest, counts, pb_path = _unpaired_fixture(tmp_path)
    payload = json.loads(open(pb_path).read())
    payload["sample_rows"] = list(reversed(payload["sample_rows"]))
    with open(pb_path, "w") as fh:
        json.dump(payload, fh)

    comparison = _comparison()
    sample_structure = _sample_structure(design_type="independent_donors_unpaired")
    provider = DD.DEDesignPreparationProvider()
    with pytest.raises(SystemExit, match="failed validation"):
        provider.prepare(
            pseudobulk_result_path=pb_path,
            comparison=comparison,
            sample_structure=sample_structure,
            out_dir=str(tmp_path / "prep_tampered_rows"),
        )


def test_prepare_carries_dataset_identity_and_manifest_identity_forward(tmp_path):
    """Per F1: `PreparedDEDesign` must carry the upstream species/gene
    namespace/counts semantics forward, and the manifest's OWN path/hash
    -- never leave a downstream consumer to guess them from a file name.
    """
    manifest, counts, pb_path = _unpaired_fixture(tmp_path)
    comparison = _comparison()
    sample_structure = _sample_structure(design_type="independent_donors_unpaired")
    provider = DD.DEDesignPreparationProvider()
    result = provider.prepare(
        pseudobulk_result_path=pb_path,
        comparison=comparison,
        sample_structure=sample_structure,
        out_dir=str(tmp_path / "prep_identity"),
    )
    assert result.status == "prepared", result.reason
    assert result.dataset_identity["species"] == manifest["dataset_identity"]["species"]
    assert (
        result.dataset_identity["gene_id_namespace"]
        == manifest["dataset_identity"]["gene_id_namespace"]
    )
    assert (
        result.dataset_identity["counts_semantics"]
        == manifest["dataset_identity"]["counts_semantics"]
    )
    assert result.pseudobulk_identity["pseudobulk_manifest_sha256"] == DD.sha256_file(pb_path)
    assert result.pseudobulk_identity["pseudobulk_manifest_path"] == os.path.abspath(pb_path)

    path = tmp_path / "prep_identity" / "PREPARED_DE_DESIGN.json"
    reloaded = DD.read_prepared_de_design(str(path))
    assert reloaded["dataset_identity"]["species"] == manifest["dataset_identity"]["species"]


# ===========================================================================
# F2 (`review:20261014_203915`): a declared column name that does
# not match what the pseudobulk aggregation actually used must be refused.
# ===========================================================================


def test_prepare_rejects_a_group_column_that_does_not_match_the_aggregation(tmp_path):
    """Reproduces the review's probe: declaring `group_column` as a column
    the aggregation never used (`condition_value` was still produced from
    the real `condition` column) must be refused, not silently accepted
    just because `condition_value` happens to already be populated.
    """
    manifest, counts, pb_path = _unpaired_fixture(tmp_path)
    comparison = _comparison(group_column="unmeasured_visit")
    sample_structure = _sample_structure(design_type="independent_donors_unpaired")
    provider = DD.DEDesignPreparationProvider()
    result = provider.prepare(
        pseudobulk_result_path=pb_path,
        comparison=comparison,
        sample_structure=sample_structure,
        out_dir=str(tmp_path / "prep_bad_group_column"),
    )
    assert result.status == "failed"
    assert "do not match what the pseudobulk aggregation actually used" in result.reason
    assert "unmeasured_visit" in result.reason


def test_prepare_rejects_a_donor_column_that_does_not_match_the_aggregation(tmp_path):
    manifest, counts, pb_path = _unpaired_fixture(tmp_path)
    comparison = _comparison()
    sample_structure = _sample_structure(
        design_type="independent_donors_unpaired", donor_or_animal_id_column="not_the_real_column"
    )
    provider = DD.DEDesignPreparationProvider()
    result = provider.prepare(
        pseudobulk_result_path=pb_path,
        comparison=comparison,
        sample_structure=sample_structure,
        out_dir=str(tmp_path / "prep_bad_donor_column"),
    )
    assert result.status == "failed"
    assert "not_the_real_column" in result.reason


# ===========================================================================
# Round-trip tamper rejection, including a genuine mutation-kill of the
# npz<->manifest identity binding (same class of bug this session's A
# increment fixed in `pseudobulk_providers.load_pseudobulk_counts`).
# ===========================================================================


def test_read_prepared_de_design_rejects_a_tampered_design_matrix(tmp_path):
    import json

    manifest, counts, pb_path = _unpaired_fixture(tmp_path)
    comparison = _comparison()
    sample_structure = _sample_structure(design_type="independent_donors_unpaired")
    provider = DD.DEDesignPreparationProvider()
    out_dir = tmp_path / "prep_tamper_dm"
    provider.prepare(
        pseudobulk_result_path=pb_path,
        comparison=comparison,
        sample_structure=sample_structure,
        out_dir=str(out_dir),
    )

    path = out_dir / "PREPARED_DE_DESIGN.json"
    payload = json.loads(path.read_text())
    payload["design_matrix"][0][1] = 1.0  # tamper one cell, hash now stale
    path.write_text(json.dumps(payload))
    with pytest.raises(SystemExit, match="failed validation"):
        DD.read_prepared_de_design(str(path))


def test_read_prepared_de_design_rejects_a_tampered_rank_claim(tmp_path):
    """`rank`/`residual_degrees_of_freedom` are RECOMPUTED from
    `design_matrix` on read, never trusted as given.
    """
    import json

    manifest, counts, pb_path = _unpaired_fixture(tmp_path)
    comparison = _comparison()
    sample_structure = _sample_structure(design_type="independent_donors_unpaired")
    provider = DD.DEDesignPreparationProvider()
    out_dir = tmp_path / "prep_tamper_rank"
    provider.prepare(
        pseudobulk_result_path=pb_path,
        comparison=comparison,
        sample_structure=sample_structure,
        out_dir=str(out_dir),
    )

    path = out_dir / "PREPARED_DE_DESIGN.json"
    payload = json.loads(path.read_text())
    payload["rank"] = 1
    path.write_text(json.dumps(payload))
    with pytest.raises(SystemExit, match="failed validation"):
        DD.read_prepared_de_design(str(path))


def test_load_prepared_de_counts_rejects_a_sample_id_order_swap(tmp_path):
    """The on-disk npz's OWN stored `sample_ids` are compared DIRECTLY
    against the manifest's claimed `sample_ids` -- swapping the manifest's
    order WITHOUT touching the npz bytes changes neither the whole-file
    hash of an unrelated-but-unwritten copy nor the content hash of the
    (unchanged) array, so only this direct identity-binding check catches
    it. Mirrors `pseudobulk_providers.load_pseudobulk_counts`'s R1 fix.
    """
    import json

    manifest, counts, pb_path = _unpaired_fixture(tmp_path)
    comparison = _comparison()
    sample_structure = _sample_structure(design_type="independent_donors_unpaired")
    provider = DD.DEDesignPreparationProvider()
    out_dir = tmp_path / "prep_swap"
    provider.prepare(
        pseudobulk_result_path=pb_path,
        comparison=comparison,
        sample_structure=sample_structure,
        out_dir=str(out_dir),
    )

    path = out_dir / "PREPARED_DE_DESIGN.json"
    payload = json.loads(path.read_text())
    payload["sample_ids"] = list(reversed(payload["sample_ids"]))  # npz bytes untouched
    path.write_text(json.dumps(payload))

    with pytest.raises(SystemExit, match="do not match the manifest's sample_ids"):
        DD.load_prepared_de_counts(payload)


def test_mutation_kill_load_prepared_de_counts_identity_binding(tmp_path, monkeypatch):
    """Per `review:20261014_221155` R1, a plain `sample_ids` reorder
    is now caught by TWO independent layers: the pre-existing npz<->manifest
    direct check (`_prepared_counts_identity_binding_problems`) AND the new
    live-upstream-source recompute (`_source_binding_problems` -- reversing
    `sample_ids` alone scrambles the sample_id->donor_id/group_value pairing
    against the live pseudobulk manifest too). Demonstrating EITHER layer in
    isolation now requires disabling BOTH for this particular mutation; see
    `test_mutation_kill_load_prepared_de_counts_source_binding` below for a
    mutation `_source_binding_problems` catches that the npz-level check
    structurally cannot (it never looks at `donor_ids`/`group_values` at all).
    """
    manifest, counts, pb_path = _unpaired_fixture(tmp_path)
    comparison = _comparison()
    sample_structure = _sample_structure(design_type="independent_donors_unpaired")
    provider = DD.DEDesignPreparationProvider()
    out_dir = tmp_path / "prep_mutation"
    provider.prepare(
        pseudobulk_result_path=pb_path,
        comparison=comparison,
        sample_structure=sample_structure,
        out_dir=str(out_dir),
    )

    import json

    path = out_dir / "PREPARED_DE_DESIGN.json"
    payload = json.loads(path.read_text())
    payload["sample_ids"] = list(reversed(payload["sample_ids"]))

    with pytest.raises(SystemExit):
        DD.load_prepared_de_counts(payload)  # real guards: reject

    monkeypatch.setattr(DD, "_prepared_counts_identity_binding_problems", lambda *a, **k: [])
    with pytest.raises(SystemExit):
        DD.load_prepared_de_counts(payload)  # one mutant disabled: the OTHER guard still rejects

    monkeypatch.setattr(DD, "_source_binding_problems", lambda *a, **k: [])
    reloaded_counts = DD.load_prepared_de_counts(payload)  # BOTH mutants: silently accepted
    assert reloaded_counts.shape == (4, 2)


def test_mutation_kill_load_prepared_de_counts_source_binding(tmp_path, monkeypatch):
    """Reproduces `review:20261014_221155` R1's
    `donor_ids_replaced_order_preserved` probe case: tampering ONLY
    `donor_ids` (leaving `sample_ids`/`gene_ids` -- the only fields
    `_prepared_counts_identity_binding_problems` examines -- unchanged)
    is accepted by the npz-level check alone, but rejected by
    `_source_binding_problems`'s live-upstream-manifest lookup. Confirms
    this new guard is independently load-bearing, not merely redundant
    with the pre-existing npz check.
    """
    manifest, counts, pb_path = _unpaired_fixture(tmp_path)
    comparison = _comparison()
    sample_structure = _sample_structure(design_type="independent_donors_unpaired")
    provider = DD.DEDesignPreparationProvider()
    out_dir = tmp_path / "prep_donor_swap"
    provider.prepare(
        pseudobulk_result_path=pb_path,
        comparison=comparison,
        sample_structure=sample_structure,
        out_dir=str(out_dir),
    )

    import json

    path = out_dir / "PREPARED_DE_DESIGN.json"
    payload = json.loads(path.read_text())
    payload["donor_ids"] = ["OTHER0", "OTHER1", "OTHER2", "OTHER3"]  # sample_ids UNCHANGED

    with pytest.raises(SystemExit, match="donor_ids entry"):
        DD.load_prepared_de_counts(payload)  # real guard: rejects

    monkeypatch.setattr(DD, "_source_binding_problems", lambda *a, **k: [])
    reloaded_counts = DD.load_prepared_de_counts(payload)  # with the mutant: silently accepted
    assert reloaded_counts.shape == (4, 2)


def _build_valid_unpaired_result(tmp_path):
    """A `PreparedDEDesign` built directly (not via `prepare()`), used to
    test the public `write_prepared_de_design` entry point's OWN guards in
    isolation -- per `review:20261014_221155` R1, a caller could
    call this entry point directly, bypassing `prepare()` entirely.
    """
    import dataclasses

    manifest, counts, pb_path = _unpaired_fixture(tmp_path)
    comparison = _comparison()
    sample_structure = _sample_structure(design_type="independent_donors_unpaired")
    core_problems, fields = DD._prepare_de_design_core(
        manifest["sample_rows"], comparison, sample_structure
    )
    assert not core_problems
    design_matrix_sha256 = DD.array_content_sha256(np.array(fields["design_matrix"], dtype=float))
    result = DD.PreparedDEDesign(
        status="prepared",
        gene_ids=tuple(manifest["gene_ids"]),
        pseudobulk_identity=dict(
            pseudobulk_manifest_path=os.path.abspath(pb_path),
            pseudobulk_manifest_sha256=DD.sha256_file(pb_path),
            aggregation_path=manifest["aggregation_path"],
            aggregation_sha256=manifest["aggregation_sha256"],
            identity_sha256=manifest["identity_sha256"],
        ),
        dataset_identity={
            k: manifest["dataset_identity"][k]
            for k in ("species", "gene_id_namespace", "counts_semantics")
        },
        comparison=dataclasses.asdict(comparison),
        sample_structure=dataclasses.asdict(sample_structure),
        design_type="independent_donors_unpaired",
        backend="design_matrix_preparation_only",
        backend_version="1",
        source_code_revision_or_hash="rev",
        dependency_versions={"python": "3"},
        design_matrix_sha256=design_matrix_sha256,
        **fields,
    )
    return result, counts, pb_path, manifest


def test_write_prepared_de_design_refuses_a_counts_shape_mismatch(tmp_path):
    result, counts, pb_path, manifest = _build_valid_unpaired_result(tmp_path)
    with pytest.raises(ValueError, match="counts shape"):
        DD.write_prepared_de_design(str(tmp_path / "bad_write"), result, counts=counts[:, :1])


# ===========================================================================
# R1 (`review:20261014_221155`): a `PreparedDEDesign`'s claims
# about its upstream pseudobulk source (identity, column names, sample/
# donor/group roles, gene order, and the counts submatrix itself) are
# RE-RESOLVED against the LIVE upstream manifest, on BOTH the public write
# path and the public read/load path -- not merely checked for internal
# self-consistency. Reproduces the review's own `probe_binding.py` cases.
# ===========================================================================


def test_write_prepared_de_design_refuses_a_reversed_counts_array(tmp_path):
    """Reproduces the review's `writer_reversed_counts` case: calling the
    PUBLIC `write_prepared_de_design` entry point directly (bypassing
    `prepare()`) with a same-shape, reversed counts array must be refused
    -- not silently written under the original (correct) upstream identity.
    """
    result, counts, pb_path, manifest = _build_valid_unpaired_result(tmp_path)
    with pytest.raises(ValueError, match="does not match the submatrix freshly recomputed"):
        DD.write_prepared_de_design(
            str(tmp_path / "writer_reversed"), result, counts=counts[::-1].copy()
        )


def test_write_prepared_de_design_refuses_a_negative_counts_array(tmp_path):
    """Reproduces the review's `writer_negative_counts` case."""
    result, counts, pb_path, manifest = _build_valid_unpaired_result(tmp_path)
    with pytest.raises(ValueError, match="negative value"):
        DD.write_prepared_de_design(str(tmp_path / "writer_negative"), result, counts=-counts)


def _prepared_payload_for_tamper(tmp_path):
    import json

    result, counts, pb_path, manifest = _build_valid_unpaired_result(tmp_path)
    out_dir = tmp_path / "prep_for_tamper"
    os.makedirs(str(out_dir))
    DD.write_prepared_de_design(str(out_dir), result, counts=counts)
    path = out_dir / "PREPARED_DE_DESIGN.json"
    return path, json.loads(path.read_text())


def test_source_binding_rejects_a_species_swap(tmp_path):
    """Reproduces the review's `species_changed_human_to_mouse` case."""
    import json

    path, payload = _prepared_payload_for_tamper(tmp_path)
    payload["dataset_identity"]["species"] = "mouse"
    path.write_text(json.dumps(payload))
    reloaded = DD.read_prepared_de_design(str(path))  # internally self-consistent -- passes here
    with pytest.raises(SystemExit, match="dataset_identity\\['species'\\]"):
        DD.load_prepared_de_counts(reloaded)


def test_source_binding_rejects_a_gene_namespace_swap(tmp_path):
    """Reproduces the review's `gene_namespace_changed` case."""
    import json

    path, payload = _prepared_payload_for_tamper(tmp_path)
    payload["dataset_identity"]["gene_id_namespace"] = "ensembl_gene_id"
    path.write_text(json.dumps(payload))
    reloaded = DD.read_prepared_de_design(str(path))
    with pytest.raises(SystemExit, match="gene_id_namespace"):
        DD.load_prepared_de_counts(reloaded)


def test_source_binding_rejects_both_declared_group_columns_changed(tmp_path):
    """Reproduces the review's `both_declared_group_columns_changed` case:
    changing BOTH `comparison.group_column` and `sample_structure.
    group_column` to a column the aggregation never used -- the F2 fix
    checks this at `prepare()`-build time; this confirms it is ALSO
    checked against the LIVE upstream at read time, for an artifact that
    was tampered AFTER being legitimately prepared.
    """
    import json

    path, payload = _prepared_payload_for_tamper(tmp_path)
    payload["comparison"]["group_column"] = "unmeasured_condition"
    payload["sample_structure"]["group_column"] = "unmeasured_condition"
    path.write_text(json.dumps(payload))
    reloaded = DD.read_prepared_de_design(str(path))
    with pytest.raises(SystemExit, match="the live upstream manifest's params"):
        DD.load_prepared_de_counts(reloaded)


def test_source_binding_rejects_a_stale_upstream_hash(tmp_path):
    """Reproduces the review's `upstream_hash_changed` case."""
    import json

    path, payload = _prepared_payload_for_tamper(tmp_path)
    payload["pseudobulk_identity"]["pseudobulk_manifest_sha256"] = "0" * 64
    path.write_text(json.dumps(payload))
    reloaded = DD.read_prepared_de_design(str(path))
    with pytest.raises(SystemExit, match="pseudobulk_manifest_sha256"):
        DD.load_prepared_de_counts(reloaded)


def test_source_binding_rejects_donor_ids_replaced_order_preserved(tmp_path):
    """Reproduces the review's `donor_ids_replaced_order_preserved` case
    (same mutation as `test_mutation_kill_load_prepared_de_counts_source_
    binding` above, here asserting the real end-to-end behavior rather
    than the mutation-kill demonstration).
    """
    import json

    path, payload = _prepared_payload_for_tamper(tmp_path)
    payload["donor_ids"] = ["OTHER0", "OTHER1", "OTHER2", "OTHER3"]
    path.write_text(json.dumps(payload))
    reloaded = DD.read_prepared_de_design(str(path))
    with pytest.raises(SystemExit, match="donor_ids entry"):
        DD.load_prepared_de_counts(reloaded)


def test_source_binding_rejects_a_nonexistent_upstream_path(tmp_path):
    """Reproduces the review's `upstream_path_nonexistent` case."""
    import json

    path, payload = _prepared_payload_for_tamper(tmp_path)
    payload["pseudobulk_identity"]["pseudobulk_manifest_path"] = str(
        tmp_path / "does_not_exist.json"
    )
    path.write_text(json.dumps(payload))
    reloaded = DD.read_prepared_de_design(str(path))
    with pytest.raises(SystemExit, match="could not be resolved"):
        DD.load_prepared_de_counts(reloaded)


def test_nonneg_finite_integer_counts_problems_unit():
    """Direct unit check of the standalone numeric-contract guard (per R1):
    non-finite, negative, and non-integer-valued counts are each refused.
    """
    assert DD._nonneg_finite_integer_counts_problems(np.array([[1.0, 2.0], [3.0, 4.0]])) == []
    assert DD._nonneg_finite_integer_counts_problems(np.array([[float("nan"), 2.0]]))
    assert DD._nonneg_finite_integer_counts_problems(np.array([[-1.0, 2.0]]))
    assert DD._nonneg_finite_integer_counts_problems(np.array([[1.5, 2.0]]))


def test_nonneg_finite_integer_counts_problems_rejects_large_scale_drift(tmp_path):
    """Per `review:20261014_231259` R1/P2: `np.allclose`'s DEFAULT
    relative tolerance already accepts an off-by-one drift at a realistic
    raw-count scale (e.g. 100000 vs 100001: `rtol=1e-5` alone tolerates an
    absolute difference of ~1) -- the integer-valued check must use EXACT
    equality, not `np.allclose`.
    """
    assert DD._nonneg_finite_integer_counts_problems(np.array([[100000.0, 2.0]])) == []
    assert DD._nonneg_finite_integer_counts_problems(np.array([[100000.25, 2.0]]))  # non-integer
    # 100001 IS itself a valid integer -- this function alone cannot know it disagrees
    # with a specific upstream value; that is `_source_binding_problems`'s job (below).
    assert DD._nonneg_finite_integer_counts_problems(np.array([[100001.0, 2.0]])) == []


def test_write_prepared_de_design_refuses_a_large_integer_drift(tmp_path):
    """Reproduces `review:20261014_231259` R1/P2's exact case:
    source counts include a legitimately large integer (100000); the
    public writer, given a counts array with ONE cell incremented by 1
    (100001, same shape, still integer-valued), must be refused by
    `_source_binding_problems`'s exact-equality recompute -- `np.allclose`
    previously accepted this silently.
    """
    manifest, counts, pb_path = _aggregate(
        tmp_path,
        "pb_large_integer",
        X=np.array([[100000, 2], [100002, 4], [200000, 6], [200002, 8]]),
        cell_ids=["CELL0", "CELL1", "CELL2", "CELL3"],
        sample_ids=["S0", "S1", "S2", "S3"],
        donor_ids=["D0", "D1", "D2", "D3"],
        conditions=["healthy", "healthy", "disease", "disease"],
    )
    comparison = _comparison()
    sample_structure = _sample_structure(design_type="independent_donors_unpaired")
    provider = DD.DEDesignPreparationProvider()
    result = provider.prepare(
        pseudobulk_result_path=pb_path,
        comparison=comparison,
        sample_structure=sample_structure,
        out_dir=str(tmp_path / "prep_large_int"),
    )
    assert result.status == "prepared", result.reason
    real_counts = DD.load_prepared_de_counts(
        DD.read_prepared_de_design(str(tmp_path / "prep_large_int" / "PREPARED_DE_DESIGN.json"))
    )

    drifted = real_counts.copy()
    drifted[0, 0] += 1  # 100000 -> 100001: still a valid integer, same shape, WRONG content
    out_dir_1 = tmp_path / "writer_drift_plus_one"
    os.makedirs(
        str(out_dir_1)
    )  # so a pre-fix silent-accept would fail on a real write, not a missing dir
    with pytest.raises(ValueError, match="does not match the submatrix freshly recomputed"):
        DD.write_prepared_de_design(str(out_dir_1), result, counts=drifted)

    fractional = real_counts.copy()
    fractional[0, 0] += 0.25  # 100000.25: not even integer-valued
    out_dir_2 = tmp_path / "writer_drift_fraction"
    os.makedirs(str(out_dir_2))
    with pytest.raises(ValueError, match="non-integer-valued"):
        DD.write_prepared_de_design(str(out_dir_2), result, counts=fractional)


# ===========================================================================
# F2 (`review:20261014_203915`): semantic fields (group roles,
# contrast direction/value, donor pairing) are REBUILT from the persisted
# identity + comparison/sample_structure and compared -- not merely a
# design_matrix hash plus an estimability check, which is VACUOUS for any
# full-rank design (`c @ pinv(X) @ X == c` identically for ANY `c`).
# ===========================================================================


def _prepared_unpaired_path(tmp_path):
    manifest, counts, pb_path = _unpaired_fixture(tmp_path)
    comparison = _comparison()
    sample_structure = _sample_structure(design_type="independent_donors_unpaired")
    provider = DD.DEDesignPreparationProvider()
    out_dir = tmp_path / "prep_semantic"
    provider.prepare(
        pseudobulk_result_path=pb_path,
        comparison=comparison,
        sample_structure=sample_structure,
        out_dir=str(out_dir),
    )
    return out_dir / "PREPARED_DE_DESIGN.json"


def _tamper_and_reload(path, mutate):
    import json

    payload = json.loads(path.read_text())
    mutate(payload)
    path.write_text(json.dumps(payload))


def test_read_prepared_de_design_rejects_a_reversed_contrast(tmp_path):
    """Reproduces the review's probe: `contrast_vector` [0,1] -> [0,-1]
    still claims 'disease minus healthy' in `contrast_direction`'s text,
    but the OLD estimability-only check is vacuous for this full-rank
    design (any vector projects back onto itself) -- only rebuilding the
    expected contrast from the persisted identity/comparison catches it.
    """
    path = _prepared_unpaired_path(tmp_path)
    _tamper_and_reload(path, lambda p: p.update(contrast_vector=[0.0, -1.0]))
    with pytest.raises(SystemExit, match="failed validation"):
        DD.read_prepared_de_design(str(path))


def test_read_prepared_de_design_rejects_a_zero_contrast(tmp_path):
    path = _prepared_unpaired_path(tmp_path)
    _tamper_and_reload(path, lambda p: p.update(contrast_vector=[0.0, 0.0]))
    with pytest.raises(SystemExit, match="failed validation"):
        DD.read_prepared_de_design(str(path))


def test_read_prepared_de_design_rejects_reversed_group_values(tmp_path):
    """Reversing `group_values` alone (D0/D1 healthy, D2/D3 disease,
    donor-sorted rows) while leaving `design_matrix`/`sample_ids`/
    `donor_ids` as originally built contradicts what rebuilding the
    design from the (now-tampered) `group_values` would produce.
    """
    path = _prepared_unpaired_path(tmp_path)
    _tamper_and_reload(path, lambda p: p.update(group_values=list(reversed(p["group_values"]))))
    with pytest.raises(SystemExit, match="failed validation"):
        DD.read_prepared_de_design(str(path))


def test_read_prepared_de_design_rejects_reversed_donor_ids(tmp_path):
    path = _prepared_unpaired_path(tmp_path)
    _tamper_and_reload(path, lambda p: p.update(donor_ids=list(reversed(p["donor_ids"]))))
    with pytest.raises(SystemExit, match="failed validation"):
        DD.read_prepared_de_design(str(path))


def test_read_prepared_de_design_rejects_a_swapped_reference_and_comparison_level(tmp_path):
    """Swapping `comparison.reference_level`/`comparison_level` while
    leaving `design_matrix`/`column_names`/`contrast_vector` as originally
    built (still 'condition[disease]') must be refused: rebuilding the
    design with the swapped roles names its contrast
    'condition[healthy]' instead, which no longer matches what is stored.
    """
    path = _prepared_unpaired_path(tmp_path)

    def mutate(p):
        c = dict(p["comparison"])
        c["reference_level"], c["comparison_level"] = c["comparison_level"], c["reference_level"]
        p["comparison"] = c

    _tamper_and_reload(path, mutate)
    with pytest.raises(SystemExit, match="failed validation"):
        DD.read_prepared_de_design(str(path))


def test_mutation_kill_validate_prepared_de_design_semantic_recompute(tmp_path, monkeypatch):
    """Demonstrates the recompute-and-compare check itself is load-bearing:
    a reversed contrast is rejected by the real guard, and silently
    accepted once `_expected_design_fields_from_identity` is monkeypatched
    into a no-op recompute that just echoes the (tampered) stored fields
    back -- i.e. a guard that trusts its own input is no guard at all.
    """
    import json

    path = _prepared_unpaired_path(tmp_path)
    payload = json.loads(path.read_text())
    payload["contrast_vector"] = [0.0, -1.0]
    path.write_text(json.dumps(payload))

    with pytest.raises(SystemExit, match="failed validation"):
        DD.read_prepared_de_design(str(path))

    def _echo(sample_ids, donor_ids, group_values, comparison, sample_structure):
        return [], dict(
            sample_ids=tuple(sample_ids),
            donor_ids=tuple(donor_ids),
            group_values=tuple(group_values),
            design_matrix=tuple(tuple(r) for r in payload["design_matrix"]),
            column_names=tuple(payload["column_names"]),
            formula=payload["formula"],
            encoding_convention=payload["encoding_convention"],
            contrast_name=payload["contrast_name"],
            contrast_vector=tuple(payload["contrast_vector"]),
            contrast_direction=payload["contrast_direction"],
            rank=payload["rank"],
            residual_degrees_of_freedom=payload["residual_degrees_of_freedom"],
            n_samples=payload["n_samples"],
            n_distinct_donors=payload["n_distinct_donors"],
        )

    monkeypatch.setattr(DD, "_expected_design_fields_from_identity", _echo)
    reloaded = DD.read_prepared_de_design(str(path))  # with the mutant: silently accepted
    assert reloaded["contrast_vector"] == [0.0, -1.0]

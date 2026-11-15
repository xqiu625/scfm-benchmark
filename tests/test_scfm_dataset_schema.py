"""Dataset config schema tests for the scFM discovery framework.

No field is silently defaulted: every test here either supplies a field
explicitly or expects validate_dataset_config to name it as missing/empty.
"""

import os
import sys

CODE = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "scripts",
    "34_scfm_discovery_framework",
    "code",
)
sys.path.insert(0, CODE)

import dataset_schema as D  # noqa: E402


def _base_config(**overrides):
    config = {
        "dataset_id": "test_monocyte_pilot",
        "species": "human",
        "tissue_or_target_cell_type": "monocyte",
        "subject_id_column": "donor_id",
        "sample_id_column": "sample_id",
        "biological_replicate_unit": "donor_id",
        "technical_batch_column": "chip_id",
        "expression_matrix_path": "/does/not/need/to/exist/for/this/test.h5ad",
        "counts_semantics": "raw_counts",
        "counts_slot": "X",
        "gene_id_namespace": "ensembl_gene_id",
        "cell_id_column": "cell_barcode",
        "qc_status": "provisional",
        "qc_version_ref": "none (provisional)",
    }
    config.update(overrides)
    return config


def test_valid_config_passes():
    result = D.validate_dataset_config(_base_config())
    assert result.ok, result.issues


def test_missing_required_field_is_named():
    config = _base_config()
    del config["species"]
    result = D.validate_dataset_config(config)
    assert not result.ok
    assert any(i.field == "species" for i in result.issues)


def test_empty_required_field_is_named():
    config = _base_config(gene_id_namespace="")
    result = D.validate_dataset_config(config)
    assert not result.ok
    assert any(i.field == "gene_id_namespace" for i in result.issues)


def test_unknown_species_is_flagged_not_silently_accepted():
    config = _base_config(species="zebrafish")
    result = D.validate_dataset_config(config)
    assert not result.ok
    assert any(i.field == "species" for i in result.issues)


def test_subject_and_sample_id_columns_must_differ():
    config = _base_config(subject_id_column="donor_id", sample_id_column="donor_id")
    result = D.validate_dataset_config(config)
    assert not result.ok
    assert any(i.field == "subject_id_column" for i in result.issues)


def test_analysis_ready_requires_qc_version_ref():
    config = _base_config(qc_status="analysis-ready", qc_version_ref="")
    result = D.validate_dataset_config(config)
    assert not result.ok
    assert any("qc_version_ref" in str(i) for i in result.issues)


def test_outcome_column_is_optional_and_never_required():
    """Omitting outcome/severity entirely must not fail validation -- it is
    optional analysis metadata, never a required inference field.
    """
    config = _base_config()
    assert "outcome_or_severity_column" not in config
    result = D.validate_dataset_config(config)
    assert result.ok, result.issues


def test_mouse_species_is_known_and_accepted():
    """Species config capability must not be hardcoded to human -- the
    framework must be able to express a mouse dataset config without reading
    any real mouse data.
    """
    config = _base_config(species="mouse", tissue_or_target_cell_type="synthetic_cell_type")
    result = D.validate_dataset_config(config)
    assert result.ok, result.issues


def test_unrecognized_extra_field_is_flagged():
    config = _base_config(typo_fielld="oops")
    result = D.validate_dataset_config(config)
    assert not result.ok
    assert any(i.field == "typo_fielld" for i in result.issues)


# --- R1 (review:20261012_153316): required fields must be type-checked,
# not merely non-None/non-empty. Each of these reproduces a case the review found
# passing with `ok=True` before this fix. ---


def test_list_typed_dataset_id_is_rejected():
    config = _base_config(dataset_id=[])
    result = D.validate_dataset_config(config)
    assert not result.ok
    assert any(i.field == "dataset_id" and "string" in i.reason for i in result.issues)


def test_int_typed_subject_id_column_is_rejected():
    config = _base_config(subject_id_column=17)
    result = D.validate_dataset_config(config)
    assert not result.ok
    assert any(i.field == "subject_id_column" and "string" in i.reason for i in result.issues)


def test_bool_typed_expression_matrix_path_is_rejected():
    config = _base_config(expression_matrix_path=False)
    result = D.validate_dataset_config(config)
    assert not result.ok
    assert any(i.field == "expression_matrix_path" and "string" in i.reason for i in result.issues)


def test_all_three_type_errors_reported_together():
    """The exact combination review:20261012_153316 R1 reproduced as
    a single, previously-passing config.
    """
    config = _base_config(dataset_id=[], subject_id_column=17, expression_matrix_path=False)
    result = D.validate_dataset_config(config)
    assert not result.ok
    bad_fields = {i.field for i in result.issues}
    assert {"dataset_id", "subject_id_column", "expression_matrix_path"} <= bad_fields


def test_counts_slot_is_required():
    config = _base_config()
    del config["counts_slot"]
    result = D.validate_dataset_config(config)
    assert not result.ok
    assert any(i.field == "counts_slot" for i in result.issues)


def test_counts_slot_accepts_x_raw_x_and_named_layer():
    for slot in ("X", "raw.X", "layers/raw_counts"):
        result = D.validate_dataset_config(_base_config(counts_slot=slot))
        assert result.ok, (slot, result.issues)


def test_counts_slot_rejects_unstructured_value():
    result = D.validate_dataset_config(_base_config(counts_slot="whatever"))
    assert not result.ok
    assert any(i.field == "counts_slot" for i in result.issues)


def test_counts_slot_rejects_empty_layer_name():
    result = D.validate_dataset_config(_base_config(counts_slot="layers/"))
    assert not result.ok
    assert any(i.field == "counts_slot" for i in result.issues)

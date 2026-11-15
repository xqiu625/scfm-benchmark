"""Integration tests for the scFM discovery framework's `validate`/`plan` CLI.

Per `review:20261012_145004` item C.5 and `review:20261012_153316`
R1/R2 (reproduced defects this round fixes): legitimate input (with a REAL donor
with two samples, not cell/sample conflated with subject), a missing required
field, duplicate/missing cell IDs, a matrix with no expression slot or wrong
declared columns, type-invalid config fields, species/declared-support conflict,
outcome missing still planning, stable AND reviewable provenance (config snapshot,
real source hash, real dependency hashes), and plan not triggering execution or
network access. Invokes the REAL `cli.py` entry point via subprocess (not just the
Python functions directly), against tiny synthetic fixtures -- never real data,
never `sbatch`. Safe to run on a login node.
"""

import json
import os
import socket
import subprocess
import sys

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CODE_DIR = os.path.join(REPO_ROOT, "scripts", "34_scfm_discovery_framework", "code")
CLI = os.path.join(CODE_DIR, "cli.py")


def _base_config(tmp_path, matrix_path=None, **overrides):
    config = {
        "dataset_id": "test_monocyte_pilot",
        "species": "human",
        "tissue_or_target_cell_type": "monocyte",
        "subject_id_column": "donor_id",
        "sample_id_column": "sample_id",
        "biological_replicate_unit": "donor_id",
        "technical_batch_column": "chip_id",
        "expression_matrix_path": matrix_path or str(tmp_path / "nonexistent.h5ad"),
        "counts_semantics": "raw_counts",
        "counts_slot": "X",
        "gene_id_namespace": "ensembl_gene_id",
        "cell_id_column": "cell_barcode",
        "qc_status": "provisional",
        "qc_version_ref": "none (provisional)",
    }
    config.update(overrides)
    return config


def _write_config(tmp_path, config, name="dataset_config.json"):
    path = tmp_path / name
    path.write_text(json.dumps(config))
    return str(path)


def _make_tiny_h5ad(
    tmp_path,
    *,
    n_vars=4,
    cell_barcodes=None,
    donor_ids=None,
    sample_ids=None,
    chip_ids=None,
    name="tiny.h5ad",
    omit_var=False,
):
    """Writes a real, tiny AnnData h5ad with real donor/sample/batch obs
    columns -- not just a bare cell barcode -- so the structural check has
    something genuine to verify. Default fixture has 6 cells: donor D1
    contributes TWO samples (S1, S2) -- a real subject-with-repeat-sampling
    example, per review:20261012_153316 item 3's instruction not to
    let cell or sample stand in for an independent subject -- and donor D2
    contributes one sample (S3).
    """
    import anndata as ad
    import numpy as np
    import pandas as pd

    if cell_barcodes is None:
        cell_barcodes = [f"CELL{i}" for i in range(6)]
    n_obs = len(cell_barcodes)
    if donor_ids is None:
        donor_ids = ["D1", "D1", "D1", "D1", "D2", "D2"][:n_obs]
    if sample_ids is None:
        sample_ids = ["S1", "S1", "S2", "S2", "S3", "S3"][:n_obs]
    if chip_ids is None:
        chip_ids = ["CHIP1", "CHIP1", "CHIP2", "CHIP2", "CHIP1", "CHIP2"][:n_obs]
    assert len(donor_ids) == n_obs and len(sample_ids) == n_obs and len(chip_ids) == n_obs

    X = np.random.default_rng(0).integers(0, 5, size=(n_obs, n_vars)).astype("float32")
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
    if omit_var:
        import h5py

        with h5py.File(path, "a") as f:
            del f["var"]
    return str(path)


def run_cli(*args):
    result = subprocess.run([sys.executable, CLI, *args], capture_output=True, text=True)
    return result


# --- validate: legitimate input, with a real donor-with-two-samples example ---


def test_validate_legitimate_input_passes(tmp_path):
    h5ad = _make_tiny_h5ad(tmp_path)
    config_path = _write_config(tmp_path, _base_config(tmp_path, matrix_path=h5ad))
    result = run_cli("validate", "--dataset-config", config_path, "--check-matrix-structure")
    assert result.returncode == 0, result.stdout + result.stderr
    report = json.loads(result.stdout)
    assert report["ok"] is True
    assert report["matrix_structure_issues"] == []


# --- validate: missing required field ---


def test_validate_missing_required_field_fails_with_reason(tmp_path):
    config = _base_config(tmp_path)
    del config["species"]
    config_path = _write_config(tmp_path, config)
    result = run_cli("validate", "--dataset-config", config_path)
    assert result.returncode == 1
    report = json.loads(result.stdout)
    assert report["ok"] is False
    assert any("species" in issue for issue in report["dataset_config_issues"])


# --- R1 (review:20261012_153316): type-invalid required fields ---


def test_validate_rejects_type_invalid_fields_before_touching_any_file(tmp_path):
    config = _base_config(
        tmp_path, dataset_id=[], subject_id_column=17, expression_matrix_path=False
    )
    config_path = _write_config(tmp_path, config)
    result = run_cli("validate", "--dataset-config", config_path, "--check-matrix-structure")
    assert result.returncode == 1
    report = json.loads(result.stdout)
    assert report["ok"] is False
    bad = " ".join(report["dataset_config_issues"])
    assert "dataset_id" in bad and "subject_id_column" in bad and "expression_matrix_path" in bad
    # Must stop at config validation -- never attempt to open the (non-existent,
    # type-invalid) matrix path.
    assert report["matrix_structure_issues"] == []
    assert report["checked_items"] == []


# --- R1: declared donor/sample/batch columns that do not actually exist ---


def test_validate_detects_declared_columns_absent_from_the_file(tmp_path):
    import anndata as ad
    import numpy as np
    import pandas as pd

    # Only cell_barcode exists -- donor/sample/batch columns are NOT in obs,
    # even though the config declares them (the first version of this fixture
    # made exactly this mistake, which is why the gap went uncaught).
    obs = pd.DataFrame({"cell_barcode": ["AAAA", "BBBB", "CCCC"]})
    var = pd.DataFrame(index=["GENE0", "GENE1"])
    a = ad.AnnData(X=np.zeros((3, 2), dtype="float32"), obs=obs, var=var)
    path = tmp_path / "barcode_only.h5ad"
    a.write_h5ad(path)

    config_path = _write_config(tmp_path, _base_config(tmp_path, matrix_path=str(path)))
    result = run_cli("validate", "--dataset-config", config_path, "--check-matrix-structure")
    assert result.returncode == 1, result.stdout
    report = json.loads(result.stdout)
    issues = " ".join(report["matrix_structure_issues"])
    assert "donor_id" in issues and "sample_id" in issues and "chip_id" in issues


# --- R1: no expression slot / no gene axis at all ---


def test_validate_detects_missing_var_group(tmp_path):
    h5ad = _make_tiny_h5ad(tmp_path, omit_var=True)
    config_path = _write_config(tmp_path, _base_config(tmp_path, matrix_path=h5ad))
    result = run_cli("validate", "--dataset-config", config_path, "--check-matrix-structure")
    assert result.returncode == 1, result.stdout
    report = json.loads(result.stdout)
    assert any("var" in i for i in report["matrix_structure_issues"])


def test_validate_detects_missing_declared_counts_slot(tmp_path):
    h5ad = _make_tiny_h5ad(tmp_path)
    config_path = _write_config(
        tmp_path, _base_config(tmp_path, matrix_path=h5ad, counts_slot="layers/raw_counts")
    )
    result = run_cli("validate", "--dataset-config", config_path, "--check-matrix-structure")
    assert result.returncode == 1, result.stdout
    report = json.loads(result.stdout)
    assert any("counts_slot" in i for i in report["matrix_structure_issues"])


# --- review:20261012_163250 finding 1: raw.X must be checked against
# raw.var, never the current (possibly-filtered) var; an unresolvable index must
# fail explicitly, never be silently skipped. ---


def _make_h5ad_with_raw(
    tmp_path, *, n_genes_raw=6, n_genes_current=2, n_obs=4, name="withraw.h5ad"
):
    """A real AnnData .raw: n_genes_raw genes saved to .raw before filtering
    down to n_genes_current genes in the current (top-level) object -- the
    legitimate case the review's finding 1 describes (raw.X genuinely holds
    MORE genes than the current var, which is not a defect).
    """
    import anndata as ad
    import numpy as np
    import pandas as pd

    obs = pd.DataFrame(
        {
            "cell_barcode": [f"CELL{i}" for i in range(n_obs)],
            "donor_id": ["D1"] * n_obs,
            "sample_id": ["S1"] * n_obs,
            "chip_id": ["C1"] * n_obs,
        }
    )
    var = pd.DataFrame(index=[f"GENE{i}" for i in range(n_genes_raw)])
    a = ad.AnnData(
        X=np.random.default_rng(0).integers(0, 5, size=(n_obs, n_genes_raw)).astype("float32"),
        obs=obs,
        var=var,
    )
    a.raw = a
    a = a[:, :n_genes_current].copy()
    path = tmp_path / name
    a.write_h5ad(path)
    return str(path)


def test_validate_accepts_legitimate_raw_x_with_more_genes_than_current_var(tmp_path):
    h5ad = _make_h5ad_with_raw(tmp_path, n_genes_raw=6, n_genes_current=2)
    config_path = _write_config(
        tmp_path, _base_config(tmp_path, matrix_path=h5ad, counts_slot="raw.X")
    )
    result = run_cli("validate", "--dataset-config", config_path, "--check-matrix-structure")
    assert result.returncode == 0, result.stdout + result.stderr
    report = json.loads(result.stdout)
    assert report["matrix_structure_issues"] == []


def test_validate_detects_genuine_raw_x_var_mismatch(tmp_path):
    """A real shape mismatch -- raw.X's own row count disagreeing with obs --
    must still be caught once compared against the RIGHT axis.
    """
    import h5py

    h5ad = _make_h5ad_with_raw(tmp_path, n_genes_raw=6, n_genes_current=2, n_obs=4)
    # Corrupt raw.var's own index length so raw.X's column count (6) disagrees
    # with it (3) -- a genuine mismatch on the correct axis.
    with h5py.File(h5ad, "a") as f:
        idx_name = f["raw/var"].attrs["_index"]
        del f["raw/var"][idx_name]
        f["raw/var"].create_dataset(idx_name, data=[b"G0", b"G1", b"G2"])
    config_path = _write_config(
        tmp_path, _base_config(tmp_path, matrix_path=h5ad, counts_slot="raw.X")
    )
    result = run_cli("validate", "--dataset-config", config_path, "--check-matrix-structure")
    assert result.returncode == 1, result.stdout
    report = json.loads(result.stdout)
    assert any("raw/var" in i for i in report["matrix_structure_issues"])


def test_validate_rejects_unresolvable_var_index_instead_of_skipping(tmp_path):
    import h5py

    h5ad = _make_tiny_h5ad(tmp_path)
    with h5py.File(h5ad, "a") as f:
        del f["var"].attrs["_index"]
    config_path = _write_config(tmp_path, _base_config(tmp_path, matrix_path=h5ad))
    result = run_cli("validate", "--dataset-config", config_path, "--check-matrix-structure")
    assert result.returncode == 1, result.stdout
    report = json.loads(result.stdout)
    assert any("could not resolve" in i for i in report["matrix_structure_issues"])


def test_validate_rejects_unresolvable_obs_index_instead_of_skipping(tmp_path):
    import h5py

    h5ad = _make_tiny_h5ad(tmp_path)
    with h5py.File(h5ad, "a") as f:
        del f["obs"].attrs["_index"]
    config_path = _write_config(tmp_path, _base_config(tmp_path, matrix_path=h5ad))
    result = run_cli("validate", "--dataset-config", config_path, "--check-matrix-structure")
    assert result.returncode == 1, result.stdout
    report = json.loads(result.stdout)
    assert any("could not resolve the 'obs' index" in i for i in report["matrix_structure_issues"])


# --- review:20261012_163250 finding 2: a relative expression_matrix_path
# must resolve against the config file's own directory, not the process cwd. ---


def test_validate_relative_matrix_path_resolves_against_config_dir_not_cwd(tmp_path):
    cfg_dir = tmp_path / "cfgdir"
    cfg_dir.mkdir()
    h5ad_abs = _make_tiny_h5ad(cfg_dir, name="tiny.h5ad")
    assert os.path.dirname(h5ad_abs) == str(cfg_dir)
    config = _base_config(tmp_path, matrix_path="tiny.h5ad")
    config_path = str(cfg_dir / "config.json")
    with open(config_path, "w") as fh:
        json.dump(config, fh)

    # Run from the config's own directory AND from the repo root -- both must
    # resolve to the SAME file and report the same result.
    result_from_cfgdir = subprocess.run(
        [
            sys.executable,
            CLI,
            "validate",
            "--dataset-config",
            "config.json",
            "--check-matrix-structure",
        ],
        capture_output=True,
        text=True,
        cwd=str(cfg_dir),
    )
    result_from_root = subprocess.run(
        [
            sys.executable,
            CLI,
            "validate",
            "--dataset-config",
            config_path,
            "--check-matrix-structure",
        ],
        capture_output=True,
        text=True,
        cwd=REPO_ROOT,
    )

    assert result_from_cfgdir.returncode == 0, result_from_cfgdir.stdout + result_from_cfgdir.stderr
    assert result_from_root.returncode == 0, result_from_root.stdout + result_from_root.stderr
    report_cfgdir = json.loads(result_from_cfgdir.stdout)
    report_root = json.loads(result_from_root.stdout)
    assert (
        report_cfgdir["expression_matrix_path_resolved"]
        == report_root["expression_matrix_path_resolved"]
        == str(cfg_dir / "tiny.h5ad")
    )
    # The declared (raw) value is preserved unresolved.
    assert report_cfgdir["expression_matrix_path_declared"] == "tiny.h5ad"


def test_plan_records_both_declared_and_resolved_matrix_path(tmp_path):
    cfg_dir = tmp_path / "cfgdir"
    cfg_dir.mkdir()
    _make_tiny_h5ad(cfg_dir, name="tiny.h5ad")
    config = _base_config(tmp_path, matrix_path="tiny.h5ad")
    config_path = str(cfg_dir / "config.json")
    with open(config_path, "w") as fh:
        json.dump(config, fh)
    out_dir = tmp_path / "plan_out"
    result = subprocess.run(
        [
            sys.executable,
            CLI,
            "plan",
            "--dataset-config",
            "config.json",
            "--out",
            str(out_dir),
            "--check-matrix-structure",
        ],
        capture_output=True,
        text=True,
        cwd=str(cfg_dir),
    )
    assert result.returncode == 0, result.stdout + result.stderr
    manifest = json.loads((out_dir / "PLAN.json").read_text())
    assert (
        manifest["input_identity"]["dataset_config_snapshot"]["expression_matrix_path"]
        == "tiny.h5ad"
    )
    assert manifest["input_identity"]["expression_matrix_path_resolved"] == str(
        cfg_dir / "tiny.h5ad"
    )


# --- validate: duplicate AND empty cell IDs ---


def test_validate_detects_duplicate_cell_ids(tmp_path):
    h5ad = _make_tiny_h5ad(
        tmp_path,
        cell_barcodes=["AAAA", "BBBB", "AAAA", "CCCC"],
        donor_ids=["D1", "D1", "D1", "D1"],
        sample_ids=["S1", "S1", "S1", "S1"],
        chip_ids=["C1", "C1", "C1", "C1"],
    )
    config_path = _write_config(tmp_path, _base_config(tmp_path, matrix_path=h5ad))
    result = run_cli("validate", "--dataset-config", config_path, "--check-matrix-structure")
    assert result.returncode == 1, result.stdout + result.stderr
    report = json.loads(result.stdout)
    assert report["ok"] is False
    assert any("duplicate value" in issue for issue in report["matrix_structure_issues"])


def test_validate_detects_empty_cell_id(tmp_path):
    h5ad = _make_tiny_h5ad(
        tmp_path,
        cell_barcodes=["AAAA", "", "CCCC"],
        donor_ids=["D1", "D1", "D1"],
        sample_ids=["S1", "S1", "S1"],
        chip_ids=["C1", "C1", "C1"],
    )
    config_path = _write_config(tmp_path, _base_config(tmp_path, matrix_path=h5ad))
    result = run_cli("validate", "--dataset-config", config_path, "--check-matrix-structure")
    assert result.returncode == 1, result.stdout
    report = json.loads(result.stdout)
    assert any("missing/empty" in i for i in report["matrix_structure_issues"])


def test_validate_accepts_explicit_no_batch_sentinel(tmp_path):
    import sys as _sys

    _sys.path.insert(0, CODE_DIR)
    import dataset_schema as D

    h5ad = _make_tiny_h5ad(tmp_path)
    config_path = _write_config(
        tmp_path,
        _base_config(tmp_path, matrix_path=h5ad, technical_batch_column=D.NO_BATCH_SENTINEL),
    )
    result = run_cli("validate", "--dataset-config", config_path, "--check-matrix-structure")
    assert result.returncode == 0, result.stdout


# --- validate: species vs. declared model support ---


def test_validate_reports_model_species_evidence_without_hard_failing(tmp_path):
    """scPRINT-2's checkpoint-metadata-confirmed mouse evidence (this round's
    upgrade) must now be treated as confirmed -- not flagged as a soft issue.
    """
    h5ad = _make_tiny_h5ad(tmp_path)
    config_path = _write_config(
        tmp_path,
        _base_config(
            tmp_path,
            matrix_path=h5ad,
            species="mouse",
            tissue_or_target_cell_type="synthetic_cell_type",
        ),
    )
    result = run_cli(
        "validate",
        "--dataset-config",
        config_path,
        "--check-matrix-structure",
        "--model-registry-key",
        "scprint2_small_v2",
    )
    assert result.returncode == 0, result.stdout + result.stderr
    report = json.loads(result.stdout)
    assert report["ok"] is True
    assert report["model_report"]["species_pretraining_evidence"] == "checkpoint_metadata_confirmed"
    assert report["soft_issues"] == []


def test_validate_reports_unconfirmed_uce_species_support_as_soft_issue(tmp_path):
    h5ad = _make_tiny_h5ad(tmp_path)
    config_path = _write_config(
        tmp_path,
        _base_config(
            tmp_path,
            matrix_path=h5ad,
            species="mouse",
            tissue_or_target_cell_type="synthetic_cell_type",
        ),
    )
    result = run_cli(
        "validate",
        "--dataset-config",
        config_path,
        "--check-matrix-structure",
        "--model-registry-key",
        "uce_4layer",
    )
    assert result.returncode == 0, result.stdout
    report = json.loads(result.stdout)
    assert report["ok"] is True
    assert report["model_report"]["species_pretraining_evidence"] == "documented_unconfirmed"
    assert any("species-pretraining evidence" in s for s in report["soft_issues"])


def test_validate_reports_unknown_model_registry_key_as_hard_failure(tmp_path):
    h5ad = _make_tiny_h5ad(tmp_path)
    config_path = _write_config(tmp_path, _base_config(tmp_path, matrix_path=h5ad))
    result = run_cli(
        "validate",
        "--dataset-config",
        config_path,
        "--model-registry-key",
        "definitely_not_a_real_model",
    )
    assert result.returncode == 1
    report = json.loads(result.stdout)
    assert report["ok"] is False


# --- plan: outcome/severity missing still produces a plan ---


def test_plan_succeeds_without_outcome_column(tmp_path):
    h5ad = _make_tiny_h5ad(tmp_path)
    config = _base_config(tmp_path, matrix_path=h5ad)
    assert "outcome_or_severity_column" not in config
    config_path = _write_config(tmp_path, config)
    out_dir = tmp_path / "plan_out"
    result = run_cli("plan", "--dataset-config", config_path, "--out", str(out_dir))
    assert result.returncode == 0, result.stdout + result.stderr
    manifest = json.loads((out_dir / "PLAN.json").read_text())
    assert manifest["stage_status"] == "planned"


# --- plan: model_undecided is explicit, not a guess ---


def test_plan_without_model_key_is_explicitly_model_undecided(tmp_path):
    h5ad = _make_tiny_h5ad(tmp_path)
    config_path = _write_config(tmp_path, _base_config(tmp_path, matrix_path=h5ad))
    out_dir = tmp_path / "plan_out"
    result = run_cli("plan", "--dataset-config", config_path, "--out", str(out_dir))
    assert result.returncode == 0
    manifest = json.loads((out_dir / "PLAN.json").read_text())
    assert manifest["model_identity_or_pending"] == "model_undecided"


# --- R2 (review:20261012_153316): plan provenance must be reviewable,
# not just a stable hash of an un-recoverable config. ---


def test_plan_input_identity_contains_a_full_recoverable_config_snapshot(tmp_path):
    h5ad = _make_tiny_h5ad(tmp_path)
    config = _base_config(tmp_path, matrix_path=h5ad)
    config_path = _write_config(tmp_path, config)
    out_dir = tmp_path / "plan_out"
    result = run_cli("plan", "--dataset-config", config_path, "--out", str(out_dir))
    assert result.returncode == 0, result.stdout
    manifest = json.loads((out_dir / "PLAN.json").read_text())
    snapshot = manifest["input_identity"]["dataset_config_snapshot"]
    assert snapshot == config
    assert manifest["input_identity"]["dataset_config_path"] == os.path.abspath(config_path)


def test_plan_source_code_revision_is_not_the_literal_string_unknown(tmp_path):
    h5ad = _make_tiny_h5ad(tmp_path)
    config_path = _write_config(tmp_path, _base_config(tmp_path, matrix_path=h5ad))
    out_dir = tmp_path / "plan_out"
    result = run_cli("plan", "--dataset-config", config_path, "--out", str(out_dir))
    assert result.returncode == 0
    manifest = json.loads((out_dir / "PLAN.json").read_text())
    rev = manifest["source_code_revision_or_hash"]
    assert rev != "unknown"
    assert rev.startswith("git:") or rev.startswith("content_hash:")


def test_plan_dependencies_list_real_source_files_with_hashes(tmp_path):
    h5ad = _make_tiny_h5ad(tmp_path)
    config_path = _write_config(tmp_path, _base_config(tmp_path, matrix_path=h5ad))
    out_dir = tmp_path / "plan_out"
    result = run_cli("plan", "--dataset-config", config_path, "--out", str(out_dir))
    assert result.returncode == 0
    manifest = json.loads((out_dir / "PLAN.json").read_text())
    deps = manifest["dependencies"]
    dep_files = {d["file"] for d in deps}
    assert dep_files == {"cli.py", "dataset_schema.py", "model_registry.py", "artifact_schema.py"}
    assert all(len(d["sha256"]) == 64 for d in deps)


def test_plan_records_validation_scope_and_unchecked_state(tmp_path):
    """Without --check-matrix-structure, the plan must say so explicitly, not
    imply the matrix was checked.
    """
    h5ad = _make_tiny_h5ad(tmp_path)
    config_path = _write_config(tmp_path, _base_config(tmp_path, matrix_path=h5ad))
    out_dir = tmp_path / "plan_out"
    result = run_cli("plan", "--dataset-config", config_path, "--out", str(out_dir))
    assert result.returncode == 0
    manifest = json.loads((out_dir / "PLAN.json").read_text())
    scope = manifest["extra"]["validation_scope"]
    assert scope["matrix_structure_checked"] is False
    assert manifest["extra"]["executable"] is False


def test_plan_with_check_matrix_structure_records_model_compatibility(tmp_path):
    h5ad = _make_tiny_h5ad(tmp_path)
    config_path = _write_config(
        tmp_path,
        _base_config(
            tmp_path,
            matrix_path=h5ad,
            species="mouse",
            tissue_or_target_cell_type="synthetic_cell_type",
        ),
    )
    out_dir = tmp_path / "plan_out"
    result = run_cli(
        "plan",
        "--dataset-config",
        config_path,
        "--out",
        str(out_dir),
        "--check-matrix-structure",
        "--model-registry-key",
        "uce_4layer",
    )
    assert result.returncode == 0, result.stdout
    manifest = json.loads((out_dir / "PLAN.json").read_text())
    assert manifest["extra"]["model_compatibility"]["model_registry_key"] == "uce_4layer"
    assert manifest["extra"]["validation_scope"]["matrix_structure_checked"] is True
    assert len(manifest["extra"]["model_compatibility_open_items"]) >= 1


def test_plan_fails_closed_when_matrix_structure_checked_and_broken(tmp_path):
    import anndata as ad
    import numpy as np
    import pandas as pd

    obs = pd.DataFrame({"cell_barcode": ["AAAA", "BBBB"]})
    var = pd.DataFrame(index=["GENE0"])
    a = ad.AnnData(X=np.zeros((2, 1), dtype="float32"), obs=obs, var=var)
    path = tmp_path / "broken.h5ad"
    a.write_h5ad(path)
    config_path = _write_config(tmp_path, _base_config(tmp_path, matrix_path=str(path)))
    out_dir = tmp_path / "plan_out"
    result = run_cli(
        "plan", "--dataset-config", config_path, "--out", str(out_dir), "--check-matrix-structure"
    )
    assert result.returncode == 1
    assert not out_dir.exists() or not list(out_dir.iterdir())


# --- provenance recorded stably across repeated runs (same input -> same hash) ---


def test_plan_provenance_hash_is_stable_across_runs(tmp_path):
    h5ad = _make_tiny_h5ad(tmp_path)
    config_path = _write_config(tmp_path, _base_config(tmp_path, matrix_path=h5ad))
    out1, out2 = tmp_path / "plan_out1", tmp_path / "plan_out2"
    r1 = run_cli("plan", "--dataset-config", config_path, "--out", str(out1))
    r2 = run_cli("plan", "--dataset-config", config_path, "--out", str(out2))
    assert r1.returncode == 0 and r2.returncode == 0
    m1 = json.loads((out1 / "PLAN.json").read_text())
    m2 = json.loads((out2 / "PLAN.json").read_text())
    assert (
        m1["input_identity"]["dataset_config_sha256"]
        == m2["input_identity"]["dataset_config_sha256"]
    )
    assert m1["source_code_revision_or_hash"] == m2["source_code_revision_or_hash"]


def test_plan_refuses_to_overwrite_nonempty_out_dir(tmp_path):
    h5ad = _make_tiny_h5ad(tmp_path)
    config_path = _write_config(tmp_path, _base_config(tmp_path, matrix_path=h5ad))
    out_dir = tmp_path / "plan_out"
    out_dir.mkdir()
    (out_dir / "junk.txt").write_text("pre-existing, must not be overwritten")
    result = run_cli("plan", "--dataset-config", config_path, "--out", str(out_dir))
    assert result.returncode != 0
    assert (out_dir / "junk.txt").read_text() == "pre-existing, must not be overwritten"


def test_plan_fails_closed_on_invalid_dataset_config(tmp_path):
    config = _base_config(tmp_path)
    del config["qc_status"]
    config_path = _write_config(tmp_path, config)
    out_dir = tmp_path / "plan_out"
    result = run_cli("plan", "--dataset-config", config_path, "--out", str(out_dir))
    assert result.returncode == 1
    assert not out_dir.exists() or not list(out_dir.iterdir())


# --- plan: does not trigger execution or network access ---


def test_plan_makes_no_network_connection(tmp_path, monkeypatch):
    """Monkeypatches socket.socket.connect to raise if ever called -- proves
    plan() does not merely avoid calling a specific library, but genuinely
    never opens a network connection, during a real subprocess-free direct
    call into cmd_plan.
    """

    def _boom(self, *a, **k):
        raise AssertionError("plan() must never open a network connection")

    monkeypatch.setattr(socket.socket, "connect", _boom)

    sys.path.insert(0, CODE_DIR)
    import cli as C  # noqa: E402

    h5ad = _make_tiny_h5ad(tmp_path)
    config_path = _write_config(tmp_path, _base_config(tmp_path, matrix_path=h5ad))
    out_dir = tmp_path / "plan_out"

    class Args:
        dataset_config = config_path
        model_registry_key = None
        check_matrix_structure = False
        out = str(out_dir)

    rc = C.cmd_plan(Args())
    assert rc == 0


def test_plan_output_contains_no_sbatch_or_submission_directive(tmp_path):
    """The plan artifact itself must never contain anything that looks like a
    job-submission instruction -- it is a reviewable document, not an
    executable one.
    """
    h5ad = _make_tiny_h5ad(tmp_path)
    config_path = _write_config(tmp_path, _base_config(tmp_path, matrix_path=h5ad))
    out_dir = tmp_path / "plan_out"
    result = run_cli("plan", "--dataset-config", config_path, "--out", str(out_dir))
    assert result.returncode == 0
    text = (out_dir / "PLAN.json").read_text()
    for forbidden in ("sbatch", "srun", "SLURM_JOB_ID"):
        assert forbidden not in text


# --- review:20261012_171916 R2: a type-invalid expression_matrix_path
# must produce a structured, non-zero exit WITHOUT an unhandled traceback -- the
# dataset-config type check must run and stop validate() before path resolution
# is ever attempted. ---


def test_validate_bool_matrix_path_fails_structured_not_traceback(tmp_path):
    config = _base_config(tmp_path, expression_matrix_path=True)
    config_path = _write_config(tmp_path, config)
    result = run_cli("validate", "--dataset-config", config_path)
    assert result.returncode == 1
    assert "Traceback" not in result.stderr
    report = json.loads(result.stdout)  # must be parseable JSON, not empty
    assert report["ok"] is False
    assert any("expression_matrix_path" in issue for issue in report["dataset_config_issues"])


def test_validate_list_matrix_path_fails_structured_not_traceback(tmp_path):
    config = _base_config(tmp_path, expression_matrix_path=["a", "b"])
    config_path = _write_config(tmp_path, config)
    result = run_cli("validate", "--dataset-config", config_path)
    assert result.returncode == 1
    assert "Traceback" not in result.stderr
    report = json.loads(result.stdout)
    assert report["ok"] is False
    assert any("expression_matrix_path" in issue for issue in report["dataset_config_issues"])

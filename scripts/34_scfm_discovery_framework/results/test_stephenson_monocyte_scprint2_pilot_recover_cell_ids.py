"""Real-fixture tests for `stephenson_monocyte_scprint2_pilot_recover_cell_ids.py`.

Every fixture is a real `AnnData` written through the real `anndata` h5ad
writer/reader (not a stub dict) so the adapter is exercised against the same file
format the real pipeline produces. `apply_recovery_to_h5ad`'s happy-path test also
re-runs the EXISTING, unmodified `stephenson_monocyte_scprint2_pilot_validate.py`
against the recovered file, closing the loop end to end (DEVELOPMENT.md rule 7): the old
validator must itself report `ids_match_as_sets=True`, `hard_failures=[]` on the
output this adapter produces -- this test does not just trust the adapter's own
report.

Scope-limit test (`test_embedding_only_same_sample_swap_is_not_detected`) asserts,
rather than silently omitting, that this adapter CANNOT catch a fault that swaps two
same-sample rows' embedding values while leaving every obs column untouched -- per
`review:20261017_213117`'s explicit instruction not to overclaim
detection of that failure mode.
"""
import importlib.util
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest

THIS_DIR = Path(__file__).resolve().parent


def _load_module(name, filename):
    spec = importlib.util.spec_from_file_location(name, THIS_DIR / filename)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


recover = _load_module("recover_cell_ids_mod", "stephenson_monocyte_scprint2_pilot_recover_cell_ids.py")
validate = _load_module("validate_mod", "stephenson_monocyte_scprint2_pilot_validate.py")


N = 6
BARCODES = [f"BARCODE{i}-DONOR1" for i in range(N)]
STRIPPED = BARCODES  # cell_id_stripped already has the steph_ prefix removed
EXPECTED_IDS = [f"steph_{b}" for b in BARCODES]


def _make_expected_input_adata():
    obs = pd.DataFrame(
        {
            "cell_id_stripped": STRIPPED,
            "patient_id": [f"P{i % 2}" for i in range(N)],
            "Sample_ID": [f"S{i % 2}" for i in range(N)],
            "gate1_flag_would_delete": [False] * N,
            "decontx_estimated_contamination": [0.01 * i for i in range(N)],
        },
        index=EXPECTED_IDS,
    )
    X = np.random.default_rng(0).random((N, 4)).astype(np.float32)
    return ad.AnnData(X=X, obs=obs)


def _make_embedded_adata_uuid_corrupted(expected):
    """Mimics the REAL defect: obs_names replaced with UUID-like strings, every
    other obs column (including cell_id_stripped) untouched, embedding present."""
    embedded = expected.copy()
    embedded.obs_names = [f"uuid-{i}" for i in range(N)]
    embedded.obsm["scprint_emb"] = np.random.default_rng(1).random((N, 8)).astype(np.float32)
    return embedded


def _write(adata, path):
    adata.write_h5ad(path)
    return str(path)


def test_recover_cell_ids_happy_path():
    expected = _make_expected_input_adata()
    embedded = _make_embedded_adata_uuid_corrupted(expected)
    recovered = recover.recover_cell_ids(embedded.obs)
    assert list(recovered) == EXPECTED_IDS


def test_recover_refuses_missing_column():
    expected = _make_expected_input_adata()
    embedded = _make_embedded_adata_uuid_corrupted(expected)
    embedded.obs = embedded.obs.drop(columns=["cell_id_stripped"])
    with pytest.raises(recover.RecoveryRefused, match="no 'cell_id_stripped' column"):
        recover.recover_cell_ids(embedded.obs)


def test_recover_refuses_null_value():
    expected = _make_expected_input_adata()
    embedded = _make_embedded_adata_uuid_corrupted(expected)
    embedded.obs.loc[embedded.obs_names[0], "cell_id_stripped"] = None
    with pytest.raises(recover.RecoveryRefused, match="null value"):
        recover.recover_cell_ids(embedded.obs)


def test_recover_refuses_duplicate_stripped_id():
    expected = _make_expected_input_adata()
    embedded = _make_embedded_adata_uuid_corrupted(expected)
    embedded.obs["cell_id_stripped"] = [STRIPPED[0]] * N  # collapse all to one id
    with pytest.raises(recover.RecoveryRefused, match="duplicate recovered id"):
        recover.recover_cell_ids(embedded.obs)


def test_check_against_manifest_refuses_extra_id():
    expected = _make_expected_input_adata()
    embedded = _make_embedded_adata_uuid_corrupted(expected)
    embedded.obs.loc[embedded.obs_names[0], "cell_id_stripped"] = "NOT-IN-MANIFEST"
    recovered = pd.Index([f"steph_{v}" for v in embedded.obs["cell_id_stripped"]])
    with pytest.raises(recover.RecoveryRefused, match="extra"):
        recover.check_against_expected_manifest(recovered, embedded.obs, expected.obs)


def test_check_against_manifest_refuses_missing_id():
    expected = _make_expected_input_adata()
    embedded = _make_embedded_adata_uuid_corrupted(expected).copy()
    embedded = embedded[:-1].copy()  # drop the last cell -> incomplete recovery
    recovered = pd.Index([f"steph_{v}" for v in embedded.obs["cell_id_stripped"]])
    with pytest.raises(recover.RecoveryRefused, match="missing"):
        recover.check_against_expected_manifest(recovered, embedded.obs, expected.obs)


def test_check_against_manifest_refuses_qc_mismatch():
    """Simulates a corrupted cell_id_stripped -> steph_ mapping that happens to
    stay unique/complete but points at the WRONG manifest row's QC values -- the
    scenario the set-only check in `recover_cell_ids` cannot catch by itself."""
    expected = _make_expected_input_adata()
    embedded = _make_embedded_adata_uuid_corrupted(expected)
    # relabel so that row 0's recovered id collides semantically with row 1's
    # identity (swap which stripped id goes with which row's QC columns), while
    # keeping the overall set of recovered ids identical to the manifest.
    stripped = list(embedded.obs["cell_id_stripped"])
    stripped[0], stripped[1] = stripped[1], stripped[0]
    embedded.obs["cell_id_stripped"] = stripped
    recovered = pd.Index([f"steph_{v}" for v in stripped])
    with pytest.raises(recover.RecoveryRefused, match="diverges from the frozen prepared-input manifest"):
        recover.check_against_expected_manifest(recovered, embedded.obs, expected.obs)


def test_apply_recovery_to_h5ad_end_to_end_and_old_validator_passes(tmp_path):
    """Mutation-kill target: this test must FAIL if `apply_recovery_to_h5ad` stops
    writing a correct `obs/_index`, or if it silently skips the cross-check."""
    expected = _make_expected_input_adata()
    embedded = _make_embedded_adata_uuid_corrupted(expected)

    expected_path = _write(expected, tmp_path / "expected_input.h5ad")
    embedded_path = _write(embedded, tmp_path / "embedded_uuid.h5ad")
    out_path = str(tmp_path / "embedded_recovered.h5ad")

    # Baseline MUST be captured before the call (review:20261017_223541
    # R3): capturing it afterward cannot detect a transient modification of the
    # source made during the call.
    with open(embedded_path, "rb") as f:
        original_bytes = f.read()

    result = recover.apply_recovery_to_h5ad(embedded_path, expected_path, out_path)
    assert result["status"] == "recovered"
    assert result["n_cells"] == N
    for col, entry in result["cross_check"].items():
        assert entry.get("n_mismatching", 0) == 0, (col, entry)

    # Original FAILED artifact must be byte-for-byte untouched.
    recovered_check = ad.read_h5ad(embedded_path)
    assert list(recovered_check.obs_names) == [f"uuid-{i}" for i in range(N)]
    with open(embedded_path, "rb") as f:
        assert f.read() == original_bytes

    # Re-run the EXISTING, unmodified validator against the recovered output.
    recovered_adata = ad.read_h5ad(out_path, backed="r")
    expected_obs = ad.read_h5ad(expected_path, backed="r").obs
    report, hard_failures = validate.validate_embedding_output(recovered_adata, expected_obs)
    assert hard_failures == [], hard_failures
    assert report["ids_match_as_sets"] is True
    assert report["order_matches_expected"] is True

    # R3: a persisted provenance manifest must exist and round-trip-read cleanly,
    # binding input/output hashes, code hash, and the mapping rule -- not just the
    # terminal JSON printed to stdout.
    import hashlib
    import json as _json

    manifest_path = result["manifest_path"]
    assert manifest_path == out_path + ".manifest.json"
    with open(manifest_path) as f:
        manifest = _json.load(f)
    assert manifest["n_cells"] == N
    with open(embedded_path, "rb") as f:
        assert manifest["embedded_h5ad_sha256"] == hashlib.sha256(f.read()).hexdigest()
    with open(out_path, "rb") as f:
        assert manifest["out_h5ad_sha256"] == hashlib.sha256(f.read()).hexdigest()
    with open(THIS_DIR / "stephenson_monocyte_scprint2_pilot_recover_cell_ids.py", "rb") as f:
        assert manifest["code_sha256"] == hashlib.sha256(f.read()).hexdigest()

    # R3: no OldFormatWarning (missing encoding-type/-version) on the patched
    # obs columns -- a real anndata re-read, not a warnings-filter bypass.
    import warnings as _warnings

    with _warnings.catch_warnings(record=True) as caught:
        _warnings.simplefilter("always")
        ad.read_h5ad(out_path)
    encoding_warnings = [w for w in caught if "encoding metadata" in str(w.message)]
    assert encoding_warnings == [], [str(w.message) for w in encoding_warnings]


def test_apply_recovery_refuses_existing_destination_without_overwriting(tmp_path):
    """Mutation-kill target for R1: must FAIL if the path-conflict guard is
    removed, since `shutil.copy2` would then silently overwrite the sentinel."""
    expected = _make_expected_input_adata()
    embedded = _make_embedded_adata_uuid_corrupted(expected)

    expected_path = _write(expected, tmp_path / "expected_input.h5ad")
    embedded_path = _write(embedded, tmp_path / "embedded_uuid.h5ad")
    out_path = tmp_path / "embedded_recovered.h5ad"
    out_path.write_bytes(b"KEEP_EXISTING_RESULT")

    with pytest.raises(recover.RecoveryRefused, match="EXISTING path"):
        recover.apply_recovery_to_h5ad(embedded_path, expected_path, str(out_path))

    assert out_path.read_bytes() == b"KEEP_EXISTING_RESULT"


def test_apply_recovery_refuses_output_path_aliasing_an_input(tmp_path):
    expected = _make_expected_input_adata()
    embedded = _make_embedded_adata_uuid_corrupted(expected)

    expected_path = _write(expected, tmp_path / "expected_input.h5ad")
    embedded_path = _write(embedded, tmp_path / "embedded_uuid.h5ad")

    with pytest.raises(recover.RecoveryRefused, match="same file as input path"):
        recover.apply_recovery_to_h5ad(embedded_path, expected_path, embedded_path)


def test_apply_recovery_refuses_destination_created_after_precheck(tmp_path):
    """Mutation-kill target for the TOCTOU race review:20261018_000122's
    independent probe demonstrated: monkeypatches shutil.copy2 (no longer called in
    the fixed code path, but kept as the exact probe shape) is replaced here with a
    direct race -- another writer creates out_path AFTER _refuse_if_path_conflict's
    check has already passed, simulating a concurrent claim. The exclusive create in
    _exclusive_create_and_copy must still refuse, and the concurrent writer's content
    must survive untouched."""
    expected = _make_expected_input_adata()
    embedded = _make_embedded_adata_uuid_corrupted(expected)

    expected_path = _write(expected, tmp_path / "expected_input.h5ad")
    embedded_path = _write(embedded, tmp_path / "embedded_uuid.h5ad")
    out_path = tmp_path / "embedded_recovered.h5ad"

    real_check = recover._refuse_if_path_conflict

    def racing_check(out_h5ad_path, *other_paths):
        real_check(out_h5ad_path, *other_paths)  # passes: out_path doesn't exist yet
        Path(out_h5ad_path).write_bytes(b"ANOTHER_WRITER_RESULT")  # race: created right after

    recover._refuse_if_path_conflict = racing_check
    try:
        with pytest.raises(recover.RecoveryRefused, match="created .* between"):
            recover.apply_recovery_to_h5ad(embedded_path, expected_path, str(out_path))
    finally:
        recover._refuse_if_path_conflict = real_check

    assert out_path.read_bytes() == b"ANOTHER_WRITER_RESULT"


def test_apply_recovery_refuses_preexisting_manifest_sidecar(tmp_path):
    """Mutation-kill target: the manifest sidecar's own exclusive-create must
    refuse (and leave the sentinel untouched) even though out_h5ad_path itself
    doesn't exist yet -- review:20261018_000122's independent
    probe found the old `open(path, "w")` silently clobbered this sentinel."""
    expected = _make_expected_input_adata()
    embedded = _make_embedded_adata_uuid_corrupted(expected)

    expected_path = _write(expected, tmp_path / "expected_input.h5ad")
    embedded_path = _write(embedded, tmp_path / "embedded_uuid.h5ad")
    out_path = tmp_path / "embedded_recovered.h5ad"
    manifest_path = Path(str(out_path) + ".manifest.json")
    manifest_path.write_text("KEEP_PRIOR_MANIFEST")

    with pytest.raises(recover.RecoveryRefused, match="manifest sidecar"):
        recover.apply_recovery_to_h5ad(embedded_path, expected_path, str(out_path))

    assert manifest_path.read_text() == "KEEP_PRIOR_MANIFEST"
    assert not out_path.exists(), "the h5ad itself must not be written when its manifest path is already taken"


def test_apply_recovery_refuses_manifest_created_after_precheck(tmp_path):
    """Mutation-kill target specifically for `_exclusive_write_json` (the
    prior test above is satisfied by the earlier advisory
    `_refuse_if_manifest_exists` check alone, since that sentinel exists
    before the call even starts -- this test races a writer in AFTER that
    early check has already passed, so only the exclusive-create at write
    time can catch it)."""
    expected = _make_expected_input_adata()
    embedded = _make_embedded_adata_uuid_corrupted(expected)

    expected_path = _write(expected, tmp_path / "expected_input.h5ad")
    embedded_path = _write(embedded, tmp_path / "embedded_uuid.h5ad")
    out_path = tmp_path / "embedded_recovered.h5ad"
    manifest_path = Path(str(out_path) + ".manifest.json")

    real_check = recover._refuse_if_manifest_exists

    def racing_check(path):
        real_check(path)  # passes: manifest doesn't exist yet
        Path(path).write_text("ANOTHER_WRITER_MANIFEST")  # race: created right after

    recover._refuse_if_manifest_exists = racing_check
    try:
        with pytest.raises(recover.RecoveryRefused, match="refusing to write manifest"):
            recover.apply_recovery_to_h5ad(embedded_path, expected_path, str(out_path))
    finally:
        recover._refuse_if_manifest_exists = real_check

    assert manifest_path.read_text() == "ANOTHER_WRITER_MANIFEST"
    # The h5ad itself WAS written (the race happens after it) -- that is the
    # recovered output's own valid content, not a half-write; only its
    # manifest sidecar lost the race. Documented, not silently hidden: a
    # recovered h5ad without a manifest is visibly incomplete, not corrupt.
    assert out_path.exists()


def test_apply_recovery_refuses_missing_required_output_column(tmp_path):
    """Mutation-kill target for R2: must FAIL if `check_against_expected_manifest`
    goes back to silently skipping an absent required column instead of refusing."""
    expected = _make_expected_input_adata()
    embedded = _make_embedded_adata_uuid_corrupted(expected)
    embedded.obs = embedded.obs.drop(columns=["patient_id"])

    expected_path = _write(expected, tmp_path / "expected_input.h5ad")
    embedded_path = _write(embedded, tmp_path / "embedded_missing_col.h5ad")
    out_path = str(tmp_path / "embedded_recovered.h5ad")

    with pytest.raises(recover.RecoveryRefused, match="required column"):
        recover.apply_recovery_to_h5ad(embedded_path, expected_path, out_path)
    assert not Path(out_path).exists()


def test_check_against_manifest_refuses_duplicate_expected_manifest_id():
    """R2: a duplicate id in the expected manifest must be an explicit
    RecoveryRefused, not an uncaught ValueError from a pandas .loc broadcast."""
    expected = _make_expected_input_adata()
    embedded = _make_embedded_adata_uuid_corrupted(expected)
    dup_expected = expected[list(range(N)) + [0]].copy()  # row 0's id appears twice
    recovered = pd.Index([f"steph_{v}" for v in embedded.obs["cell_id_stripped"]])
    with pytest.raises(recover.RecoveryRefused, match="duplicate id"):
        recover.check_against_expected_manifest(recovered, embedded.obs, dup_expected.obs)


def test_recover_refuses_empty_string_stripped_id():
    expected = _make_expected_input_adata()
    embedded = _make_embedded_adata_uuid_corrupted(expected)
    embedded.obs.loc[embedded.obs_names[0], "cell_id_stripped"] = ""
    with pytest.raises(recover.RecoveryRefused, match="empty-string"):
        recover.recover_cell_ids(embedded.obs)


def test_apply_recovery_refuses_when_cross_check_mismatches(tmp_path):
    """Mutation-kill target: confirms the refusal path actually fires end to end
    through the h5ad-backed function, not only through the in-memory helper."""
    expected = _make_expected_input_adata()
    embedded = _make_embedded_adata_uuid_corrupted(expected)
    stripped = list(embedded.obs["cell_id_stripped"])
    stripped[0], stripped[1] = stripped[1], stripped[0]
    embedded.obs["cell_id_stripped"] = stripped

    expected_path = _write(expected, tmp_path / "expected_input.h5ad")
    embedded_path = _write(embedded, tmp_path / "embedded_corrupt.h5ad")
    out_path = str(tmp_path / "embedded_recovered.h5ad")

    with pytest.raises(recover.RecoveryRefused):
        recover.apply_recovery_to_h5ad(embedded_path, expected_path, out_path)


def test_embedding_only_same_sample_swap_is_not_detected(tmp_path):
    """Scope-limit test, asserted rather than silently assumed: swapping two
    same-sample rows' `scprint_emb` values while leaving every obs column
    (including `cell_id_stripped`) in place is a real failure mode this adapter
    (an obs-only recovery) CANNOT see, because it never reads `obsm` at all. The
    recovered file passes both this adapter's own checks AND the old validator,
    while the two cells' embeddings are, by construction, swapped. This failure
    mode must be caught upstream (row-order/producer-side instrumentation), not
    retroactively here -- this test exists so that claim is never silently
    overstated."""
    expected = _make_expected_input_adata()
    embedded = _make_embedded_adata_uuid_corrupted(expected)
    # same patient/sample for rows 0 and 1 by construction (i % 2 pattern with i=0,1
    # differ -- force them onto the same sample explicitly for this scenario):
    embedded.obs.loc[embedded.obs_names[0], "patient_id"] = "SAME"
    embedded.obs.loc[embedded.obs_names[1], "patient_id"] = "SAME"
    embedded.obs.loc[embedded.obs_names[0], "Sample_ID"] = "SAME"
    embedded.obs.loc[embedded.obs_names[1], "Sample_ID"] = "SAME"
    expected.obs.loc[expected.obs_names[0], "patient_id"] = "SAME"
    expected.obs.loc[expected.obs_names[1], "patient_id"] = "SAME"
    expected.obs.loc[expected.obs_names[0], "Sample_ID"] = "SAME"
    expected.obs.loc[expected.obs_names[1], "Sample_ID"] = "SAME"
    original_emb = embedded.obsm["scprint_emb"].copy()
    swapped_emb = original_emb.copy()
    swapped_emb[[0, 1]] = swapped_emb[[1, 0]]
    embedded.obsm["scprint_emb"] = swapped_emb

    expected_path = _write(expected, tmp_path / "expected_input.h5ad")
    embedded_path = _write(embedded, tmp_path / "embedded_emb_swapped.h5ad")
    out_path = str(tmp_path / "embedded_recovered.h5ad")

    result = recover.apply_recovery_to_h5ad(embedded_path, expected_path, out_path)
    assert result["status"] == "recovered"  # the adapter has no way to refuse this

    recovered_adata = ad.read_h5ad(out_path, backed="r")
    expected_obs = ad.read_h5ad(expected_path, backed="r").obs
    report, hard_failures = validate.validate_embedding_output(recovered_adata, expected_obs)
    assert hard_failures == []  # the OLD validator ALSO cannot see this -- documented, not fixed, here

    # The embedding rows for cells 0 and 1 are in fact still swapped relative to
    # their now-correct obs labels -- this is the gap left open for the
    # row-order/permutation-invariance smoke test, not closed by this file.
    recovered_emb = np.asarray(recovered_adata.obsm["scprint_emb"])
    assert np.allclose(recovered_emb[0], original_emb[1])
    assert np.allclose(recovered_emb[1], original_emb[0])

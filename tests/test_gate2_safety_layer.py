"""Tests for the Gate 2 mechanical safety layer (execution-plan P0 items 2-5).

Everything writes under tmp_path. Scrublet is injected as a fake so success, a
warning-only NaN threshold, and a runtime exception are all deterministic --
the real detector is never invoked and no cohort data is read.

These tests assert that validation REFUSES things. If one starts failing,
the fix is in the caller, not in the validation.
"""

import csv
import os
import sys

import numpy as np
import pytest
from scipy import sparse

CODE = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "scripts",
    "32_rebuild_upfront_qc",
    "code",
)
sys.path.insert(0, CODE)

import gate2_checkpoint as ckpt  # noqa: E402
import gate2_diagnostics as dg  # noqa: E402
import gate2_run as runner  # noqa: E402
import gate2_score_core as core  # noqa: E402
import gate2_status as st  # noqa: E402
from gate2_runtime_config import SCHEMA_VERSION, Gate2ConfigError, load_runtime_config  # noqa: E402

N_CELLS = 120  # above MIN_CELLS_FOR_SCRUBLET
N_GENES = 40


# --------------------------------------------------------------------------
# fakes and fixtures
# --------------------------------------------------------------------------

RECEIVED_KWARGS = {}


def fake_scrublet_success(adata, expected_doublet_rate=None, random_state=0, **kw):
    RECEIVED_KWARGS.clear()
    RECEIVED_KWARGS.update(kw)
    n = adata.n_obs
    rng = np.random.default_rng(random_state)
    scores = rng.random(n)
    adata.obs["doublet_score"] = scores
    adata.obs["predicted_doublet"] = scores > 0.9
    adata.uns["scrublet"] = {"threshold": 0.9, "doublet_scores_sim": rng.random(2 * n)}


def fake_scrublet_nan_threshold(adata, expected_doublet_rate=None, random_state=0, **kw):
    """Scanpy's real behaviour when auto-thresholding fails: it WARNS, returns
    finite scores and an all-False prediction, and leaves threshold unusable."""
    import warnings

    n = adata.n_obs
    warnings.warn("Failed to automatically identify doublet score threshold")
    adata.obs["doublet_score"] = np.linspace(0.0, 0.5, n)
    adata.obs["predicted_doublet"] = np.zeros(n, dtype=bool)
    adata.uns["scrublet"] = {"threshold": float("nan")}


def fake_scrublet_raises(adata, expected_doublet_rate=None, random_state=0, **kw):
    raise ValueError("synthetic runtime failure")


def make_matrix(n_cells=N_CELLS, n_genes=N_GENES, seed=0):
    rng = np.random.default_rng(seed)
    return sparse.csr_matrix((rng.integers(0, 5, (n_cells, n_genes))).astype(float))


def cell_ids(n=N_CELLS, prefix="c"):
    return [f"{prefix}{i}" for i in range(n)]


def config_row(cohort, channel, **over):
    row = {
        "schema_version": SCHEMA_VERSION,
        "cohort": cohort,
        "channel_id": channel,
        "platform": "10x",
        "gate2_input_stage": "post_gate1_scrublet_input",
        "expected_doublet_rate": "0.05",
        "rate_estimand": "multiplet prevalence of the post_gate1_scrublet_input matrix",
        "rate_source": "TEST FIXTURE -- not an authorized rate",
        "decision_reference": "docs/decisions.md#test-fixture",
        "runtime_eligible": "true",
        "value_basis": "arbitrary",
        "n_candidate_barcodes": "",
        "n_called_pre_qc": "",
        "n_post_demux_singlets": "",
        "n_gate1_retained": str(N_CELLS),
        "n_donors": "1",
        "donor_proportion_stage": "post_gate1",
        "donor_proportions_checksum": "none",
        "source_artifact": "fixture",
        "source_checksum": "fixture",
        "code_provenance": "test",
    }
    row.update(over)
    return row


def write_config(path, rows):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    return str(path)


def write_combined_cohort(
    d,
    cohort="t",
    channels=("A",),
    n=N_CELLS,
    seed=0,
    reorder=False,
    id_header="cell_id",
    mangle=None,
):
    """Minimal real-structure combined-layout input.

    `reorder` permutes BOTH the cells.csv rows and the matrix columns together,
    so the capture is the same cells in a different order -- the realistic
    perturbation, not a corrupted file.
    """
    from scipy.io import mmwrite

    os.makedirs(d, exist_ok=True)
    per = n // len(channels)
    rows, ch_col = [], []
    for ch in channels:
        rows += [f"{ch}_cell{i}" for i in range(per)]
        ch_col += [ch] * per
    order = np.arange(len(rows))
    if reorder:
        order = np.concatenate([np.arange(per)[::-1] + i * per for i in range(len(channels))])
    rows = [rows[i] for i in order]
    ch_col = [ch_col[i] for i in order]
    if mangle == "duplicate":
        rows[1] = rows[0]
    elif mangle == "empty":
        rows[1] = ""
    elif mangle == "na":
        rows[1] = "nan"
    elif mangle == "casefold_collision":
        rows[1] = rows[0].upper()
    full = make_matrix(len(rows), N_GENES, seed)
    mat = full[order].T.tocsc()  # genes x cells, permuted consistently
    mmwrite(os.path.join(d, "counts.mtx"), mat)
    with open(os.path.join(d, "cells.csv"), "w", newline="") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow([id_header, "channel_id", "donor_id"])
        for cid, ch in zip(rows, ch_col):
            w.writerow([cid, ch, "d1"])
    return d


# --------------------------------------------------------------------------
# 1-6: status model
# --------------------------------------------------------------------------


def test_01_successful_call_only_doublet_or_noncandidate():
    res = core.score_one_capture(
        make_matrix(),
        cell_ids(),
        expected_doublet_rate=0.05,
        seed=1,
        scrublet_fn=fake_scrublet_success,
    )
    assert res.capture_status == st.CAPTURE_SCORED
    assert set(res.call_status) <= {st.CALL_PREDICTED_DOUBLET, st.CALL_SCRUBLET_NONCANDIDATE}
    assert set(np.unique(res.candidate)) <= {0, 1}
    assert -1 not in set(np.unique(res.candidate))


def test_02_too_small_capture():
    res = core.score_one_capture(
        make_matrix(10),
        cell_ids(10),
        expected_doublet_rate=0.05,
        seed=1,
        scrublet_fn=fake_scrublet_success,
    )
    assert res.capture_status == st.CAPTURE_TOO_SMALL
    assert set(res.call_status) == {st.CALL_TOO_SMALL}
    assert set(np.unique(res.candidate)) == {-1}, "unknown must be -1, never 0/False"


def test_03_nan_threshold_becomes_threshold_failed_not_singlet():
    """The known Wilk failure: finite scores + all-False predictions + NaN threshold."""
    res = core.score_one_capture(
        make_matrix(),
        cell_ids(),
        expected_doublet_rate=0.05,
        seed=1,
        scrublet_fn=fake_scrublet_nan_threshold,
    )
    assert res.capture_status == st.CAPTURE_THRESHOLD_FAILED
    assert set(res.call_status) == {st.CALL_THRESHOLD_FAILED}
    assert st.CALL_SCRUBLET_NONCANDIDATE not in res.call_status
    assert set(np.unique(res.candidate)) == {-1}
    assert np.all(np.isfinite(res.scores)), "scores are evidence and are retained"
    assert any("threshold" in p for p in res.diagnostics["problems"])
    assert res.diagnostics["warnings"], "the scanpy warning is captured"


def test_04_runtime_exception():
    res = core.score_one_capture(
        make_matrix(),
        cell_ids(),
        expected_doublet_rate=0.05,
        seed=1,
        scrublet_fn=fake_scrublet_raises,
    )
    assert res.capture_status == st.CAPTURE_RUNTIME_FAILED
    assert set(res.call_status) == {st.CALL_RUNTIME_FAILED}
    assert set(np.unique(res.candidate)) == {-1}


def test_05_multiplexed_capture_records_donors():
    props = np.array([0.5, 0.3, 0.2])
    res = core.score_one_capture(
        make_matrix(),
        cell_ids(),
        expected_doublet_rate=0.05,
        seed=1,
        donor_proportions=props,
        scrublet_fn=fake_scrublet_success,
    )
    assert res.diagnostics["n_donors_in_capture"] == 3


def test_06_single_donor_capture():
    res = core.score_one_capture(
        make_matrix(),
        cell_ids(),
        expected_doublet_rate=0.05,
        seed=1,
        donor_proportions=np.array([1.0]),
        scrublet_fn=fake_scrublet_success,
    )
    assert res.diagnostics["n_donors_in_capture"] == 1


def test_06b_mechanical_candidate_is_exactly_predicted_doublet():
    every = list(st.CALL_STATUSES)
    mask = st.mechanical_candidate_mask(every)
    assert mask == [s == st.CALL_PREDICTED_DOUBLET for s in every]
    # the other five statuses are all False -- checked explicitly, not implied
    for status, flagged in zip(every, mask):
        if status != st.CALL_PREDICTED_DOUBLET:
            assert flagged is False, f"{status} must not be a mechanical candidate"
    assert st.MECHANICAL_DOUBLET_CANDIDATE_STATUSES == {st.CALL_PREDICTED_DOUBLET}


def test_06c_method_unreliable_exists_but_is_policy_assigned_only():
    assert st.CALL_METHOD_UNRELIABLE in st.CALL_STATUSES
    assert st.CALL_METHOD_UNRELIABLE in st.POLICY_ASSIGNED_ONLY
    src = open(os.path.join(CODE, "gate2_score_core.py")).read()
    assert "CALL_METHOD_UNRELIABLE" not in src, "no code may assign it without a policy decision"


# --------------------------------------------------------------------------
# 7-12: fail-closed runtime config
# --------------------------------------------------------------------------


def test_07_missing_runtime_config_fails():
    with pytest.raises(Gate2ConfigError, match="no --runtime-config"):
        load_runtime_config(None, "t", ["A"])
    with pytest.raises(Gate2ConfigError, match="not found"):
        load_runtime_config("/nonexistent/x.csv", "t", ["A"])


def test_08_missing_channel_fails(tmp_path):
    p = write_config(tmp_path / "c.csv", [config_row("t", "A")])
    with pytest.raises(Gate2ConfigError, match="no runtime config row"):
        load_runtime_config(p, "t", ["A", "B"])


def test_09_duplicate_channel_fails(tmp_path):
    p = write_config(tmp_path / "c.csv", [config_row("t", "A"), config_row("t", "A")])
    with pytest.raises(Gate2ConfigError, match="duplicate channel_id"):
        load_runtime_config(p, "t", ["A"])


@pytest.mark.parametrize("bad", ["0", "1", "1.5", "-0.1", "nan", "inf", "abc", ""])
def test_10_invalid_rate_fails(tmp_path, bad):
    p = write_config(tmp_path / "c.csv", [config_row("t", "A", expected_doublet_rate=bad)])
    with pytest.raises(Gate2ConfigError):
        load_runtime_config(p, "t", ["A"])


def test_11_unknown_platform_fails(tmp_path):
    p = write_config(tmp_path / "c.csv", [config_row("t", "A", platform="dropseq")])
    with pytest.raises(Gate2ConfigError, match="not one of"):
        load_runtime_config(p, "t", ["A"])


def test_12_t7_path_rejected_as_runtime_config(tmp_path):
    t7 = tmp_path / "T7_doublet_prior_evidence"
    t7.mkdir()
    p = write_config(t7 / "rate-sensitivity-scenarios-v3.csv", [config_row("t", "A")])
    with pytest.raises(Gate2ConfigError, match="T7 evidence layer"):
        load_runtime_config(p, "t", ["A"])


def test_12b_runtime_eligible_must_be_exactly_true(tmp_path):
    for v in ("false", "FALSE", "", "1", "yes", "True "):
        p = write_config(tmp_path / "c.csv", [config_row("t", "A", runtime_eligible=v)])
        if v.strip().lower() == "true":
            continue
        with pytest.raises(Gate2ConfigError, match="runtime_eligible"):
            load_runtime_config(p, "t", ["A"])


def test_12c_missing_decision_reference_fails(tmp_path):
    p = write_config(tmp_path / "c.csv", [config_row("t", "A", decision_reference="")])
    with pytest.raises(Gate2ConfigError, match="decision_reference"):
        load_runtime_config(p, "t", ["A"])


def test_12d_stage_estimand_mismatch_fails(tmp_path):
    p = write_config(
        tmp_path / "c.csv", [config_row("t", "A", rate_estimand="capture loading multiplet rate")]
    )
    with pytest.raises(Gate2ConfigError, match="does not describe"):
        load_runtime_config(p, "t", ["A"])
    p = write_config(tmp_path / "c2.csv", [config_row("t", "A", gate2_input_stage="raw_capture")])
    with pytest.raises(Gate2ConfigError, match="gate2_input_stage"):
        load_runtime_config(p, "t", ["A"])


# --------------------------------------------------------------------------
# 13-22: checkpoints
# --------------------------------------------------------------------------


def _fp(**over):
    base = dict(
        checkpoint_schema_version=ckpt.CHECKPOINT_SCHEMA_VERSION,
        call_status_schema=list(st.CALL_STATUSES),
        cohort="t",
        channel_id_raw="A",
        n_cells=3,
        cell_id_checksum=ckpt.ordered_cell_checksum(["a", "b", "c"]),
        matrix_fingerprint="m0",
        config_row_checksum="cfg0",
        expected_doublet_rate=0.05,
        donor_proportions_checksum="dp0",
        donor_proportion_stage="post_gate1",
        platform="10x",
        base_seed=7,
        capture_seed=99,
        scrublet_params={"x": 1},
        software_versions={"numpy": "1"},
        git_commit="abc",
        code_fingerprint="code0",
    )
    base.update(over)
    return ckpt.build_fingerprint(**base)


def _conforming_diagnostics(**over):
    """What real code writes: diagnostics carrying their schema version.

    Helpers must produce this shape, or they exercise something the pipeline
    never emits. Round-8 tests that deliberately omit or corrupt the version
    pass it explicitly instead of using this.
    """
    d = {"threshold": 0.9, "diagnostics_schema_version": dg.DIAGNOSTICS_SCHEMA_VERSION}
    d.update(over)
    return d


def _save(d, fp, key="k"):
    return ckpt.save_checkpoint(
        str(d),
        key,
        fp,
        scores=[0.1, 0.2, 0.3],
        call_status=[
            st.CALL_SCRUBLET_NONCANDIDATE,
            st.CALL_PREDICTED_DOUBLET,
            st.CALL_SCRUBLET_NONCANDIDATE,
        ],
        candidate=np.array([0, 1, 0], dtype=np.int8),
        cell_ids=["a", "b", "c"],
        capture_status=st.CAPTURE_SCORED,
        diagnostics=_conforming_diagnostics(),
    )


def test_13_fresh_vs_resume_equal_by_cell_id(tmp_path):
    inp = write_combined_cohort(tmp_path / "in", "t", ("A",))
    cfg = write_config(tmp_path / "c.csv", [config_row("t", "A")])
    out = tmp_path / "out"
    r1 = runner.run("t", "combined", inp, str(out), cfg, 7, scrublet_fn=fake_scrublet_success)
    assert r1["restored"] == 0
    r2 = runner.run("t", "combined", inp, str(out), cfg, 7, scrublet_fn=fake_scrublet_success)
    assert r2["restored"] == 1, "second run must restore from checkpoint"

    def read(p):
        # keyed on the AUTHORITATIVE identity, not cell_id or the derived string
        return {(r["cohort"], r["channel_id"], r["cell_id"]): r for r in csv.DictReader(open(p))}

    a, b = read(r1["calls_csv"]), read(r2["calls_csv"])
    assert a.keys() == b.keys()
    for cid in a:
        assert a[cid]["scrublet_score"] == b[cid]["scrublet_score"]
        assert a[cid]["call_status"] == b[cid]["call_status"]
        assert a[cid]["candidate_doublet"] == b[cid]["candidate_doublet"]
    d1 = list(csv.DictReader(open(r1["diagnostics_csv"])))
    d2 = list(csv.DictReader(open(r2["diagnostics_csv"])))
    assert [x["channel_id"] for x in d1] == [x["channel_id"] for x in d2]
    assert d1[0]["threshold"] == d2[0]["threshold"]


@pytest.mark.parametrize(
    "field,value",
    [
        ("cell_id_checksum", "DIFFERENT"),  # 14 changed cell order
        ("expected_doublet_rate", 0.2),  # 15 changed rate
        ("config_row_checksum", "cfg-other"),  # 15 changed config row
        ("donor_proportions_checksum", "dp1"),  # 16 changed donor proportions
        ("matrix_fingerprint", "m1"),  # 17 changed matrix
        ("capture_seed", 1234),
        ("code_fingerprint", "code1"),
        ("software_versions", {"numpy": "2"}),
    ],
)
def test_14_to_17_any_fingerprint_change_invalidates(tmp_path, field, value):
    _save(tmp_path, _fp())
    with pytest.raises(ckpt.CheckpointMismatch, match="fingerprint mismatch"):
        ckpt.load_checkpoint(
            str(tmp_path), "k", _fp(**{field: value}), expected_cell_ids=["a", "b", "c"]
        )


def test_18_old_schema_checkpoint_rejected(tmp_path):
    _save(tmp_path, _fp(checkpoint_schema_version="gate2-checkpoint/0"))
    with pytest.raises(ckpt.CheckpointMismatch, match="checkpoint_schema_version"):
        ckpt.load_checkpoint(str(tmp_path), "k", _fp(), expected_cell_ids=["a", "b", "c"])


def test_19_truncated_file_rejected_and_temp_never_visible(tmp_path):
    path = _save(tmp_path, _fp())
    with open(path, "r+b") as fh:
        fh.truncate(20)
    with pytest.raises(ckpt.CheckpointMismatch, match="unreadable or corrupt"):
        ckpt.load_checkpoint(str(tmp_path), "k", _fp(), expected_cell_ids=["a", "b", "c"])
    # a leftover temp file must not be picked up under the real name
    (tmp_path / ".other.npz.tmp").write_bytes(b"partial")
    assert (
        ckpt.load_checkpoint(
            str(tmp_path), "other", _fp(channel_id_raw="other"), expected_cell_ids=["a", "b", "c"]
        )
        is None
    )


def test_20_channel_names_a_dash_b_and_a_underscore_b_do_not_collide(tmp_path):
    k1, k2 = ckpt.capture_key("t", "A-B"), ckpt.capture_key("t", "A_B")
    assert k1 != k2, "the old sanitizer mapped both to A_B and shared one file"
    _save(tmp_path, _fp(channel_id_raw="A-B"), key=k1)
    _save(tmp_path, _fp(channel_id_raw="A_B"), key=k2)
    assert (
        ckpt.load_checkpoint(
            str(tmp_path), k1, _fp(channel_id_raw="A-B"), expected_cell_ids=["a", "b", "c"]
        )
        is not None
    )
    assert (
        ckpt.load_checkpoint(
            str(tmp_path), k2, _fp(channel_id_raw="A_B"), expected_cell_ids=["a", "b", "c"]
        )
        is not None
    )


def test_21_split_resume_reaggregates_every_capture(tmp_path):
    from scipy.io import mmwrite

    root = tmp_path / "in" / "split_per_channel"
    chans = ["A-B", "A_B", "C"]
    for ch in chans:
        d = root / ch
        d.mkdir(parents=True)
        mmwrite(str(d / "counts.mtx"), make_matrix(N_CELLS, N_GENES).T.tocsc())
        with open(d / "cells.csv", "w", newline="") as fh:
            w = csv.writer(fh, lineterminator="\n")
            w.writerow(["cell_id", "channel_id", "donor_id"])
            for i in range(N_CELLS):
                w.writerow([f"{ch}_c{i}", ch, "d1"])
    cfg = write_config(tmp_path / "c.csv", [config_row("t", c) for c in chans])
    out = tmp_path / "out"
    runner.run(
        "t", "split", str(tmp_path / "in"), str(out), cfg, 7, scrublet_fn=fake_scrublet_success
    )
    r2 = runner.run(
        "t", "split", str(tmp_path / "in"), str(out), cfg, 7, scrublet_fn=fake_scrublet_success
    )
    assert r2["restored"] == len(chans)
    diag = list(csv.DictReader(open(r2["diagnostics_csv"])))
    assert len(diag) == len(chans), "resumed run must emit one row per capture, not only new ones"
    assert sorted(d["channel_id"] for d in diag) == sorted(chans)


def test_22_restored_array_wrong_length_rejected(tmp_path):
    _save(tmp_path, _fp())
    with pytest.raises(ckpt.CheckpointMismatch, match="lengths"):
        ckpt.load_checkpoint(str(tmp_path), "k", _fp(), expected_cell_ids=["a", "b"])


def test_22b_restored_cell_identity_mismatch_rejected(tmp_path):
    _save(tmp_path, _fp())
    fp = _fp(cell_id_checksum=ckpt.ordered_cell_checksum(["a", "b", "c"]))
    with pytest.raises(ckpt.CheckpointMismatch, match="cell identity/order"):
        ckpt.load_checkpoint(str(tmp_path), "k", fp, expected_cell_ids=["c", "b", "a"])


def test_22c_seed_is_deterministic_and_identity_bound():
    assert ckpt.capture_seed(7, "t", "A") == ckpt.capture_seed(7, "t", "A")
    assert ckpt.capture_seed(7, "t", "A") != ckpt.capture_seed(8, "t", "A")
    assert ckpt.capture_seed(7, "t", "A") != ckpt.capture_seed(7, "t", "B")
    assert ckpt.capture_seed(7, "t", "A-B") != ckpt.capture_seed(7, "t", "A_B")


def test_22d_checkpoint_does_not_use_pickle(tmp_path):
    _save(tmp_path, _fp())
    with np.load(str(tmp_path / "k.npz"), allow_pickle=False) as z:
        assert set(z.files) >= {"meta", "scores", "call_status", "candidate", "cell_ids"}


def test_22e_runner_refuses_without_config(tmp_path):
    inp = write_combined_cohort(tmp_path / "in", "t", ("A",))
    with pytest.raises(Gate2ConfigError):
        runner.run(
            "t", "combined", inp, str(tmp_path / "out"), None, 7, scrublet_fn=fake_scrublet_success
        )


# ==========================================================================
# Round 2 — fixes for the read-only audit findings B1, B2, M1, M2, M4
# ==========================================================================

# ---- B1: cell_id must be present, never guessed ---------------------------


def test_b1_missing_cell_id_column_raises(tmp_path):
    inp = write_combined_cohort(tmp_path / "in", "t", ("A",), id_header="barcode")
    cfg = write_config(tmp_path / "c.csv", [config_row("t", "A")])
    with pytest.raises(ckpt.CellIdentityError) as e:
        runner.run(
            "t", "combined", inp, str(tmp_path / "out"), cfg, 7, scrublet_fn=fake_scrublet_success
        )
    msg = str(e.value)
    assert "cell_id" in msg and "cohort='t'" in msg and "cells.csv" in msg
    assert "barcode" in msg, "the error shows what columns were actually found"


def test_b1_split_layout_also_requires_cell_id(tmp_path):
    from scipy.io import mmwrite

    d = tmp_path / "in" / "split_per_channel" / "A"
    d.mkdir(parents=True)
    mmwrite(str(d / "counts.mtx"), make_matrix().T.tocsc())
    with open(d / "cells.csv", "w", newline="") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(["barcode", "channel_id", "donor_id"])
        for i in range(N_CELLS):
            w.writerow([f"c{i}", "A", "d1"])
    cfg = write_config(tmp_path / "c.csv", [config_row("t", "A")])
    with pytest.raises(ckpt.CellIdentityError, match="cell_id"):
        runner.run(
            "t",
            "split",
            str(tmp_path / "in"),
            str(tmp_path / "out"),
            cfg,
            7,
            scrublet_fn=fake_scrublet_success,
        )


# ---- B2: cell IDs must be valid and unique --------------------------------


# NOTE on the expected messages: pandas parses both an empty CSV field and the
# literal "nan" as a float NaN, so those two reach the validator through the
# missing-value branch rather than the empty-string branch. Both are refusals;
# the patterns below record which branch actually fires rather than relaxing the
# check to make one message fit.
@pytest.mark.parametrize(
    "mangle,pattern",
    [
        ("duplicate", "duplicate cell ID"),
        ("empty", "is missing"),  # empty CSV field -> NaN via pandas
        ("na", "is missing"),  # literal "nan"   -> NaN via pandas
        ("casefold_collision", "collide after whitespace/case"),
    ],
)
def test_b2_invalid_cell_ids_fail_closed(tmp_path, mangle, pattern):
    inp = write_combined_cohort(tmp_path / "in", "t", ("A",), mangle=mangle)
    cfg = write_config(tmp_path / "c.csv", [config_row("t", "A")])
    with pytest.raises(ckpt.CellIdentityError, match=pattern):
        runner.run(
            "t", "combined", inp, str(tmp_path / "out"), cfg, 7, scrublet_fn=fake_scrublet_success
        )


def test_b2_empty_string_branch_is_reachable_and_rejected():
    """pandas turns an empty CSV field into NaN, so this branch is exercised
    directly -- it still has to hold for any non-CSV caller."""
    with pytest.raises(ckpt.CellIdentityError, match="empty or null-like"):
        ckpt.validate_cell_ids(["a", "   "], 2)
    with pytest.raises(ckpt.CellIdentityError, match="empty or null-like"):
        ckpt.validate_cell_ids(["a", "<NA>"], 2)


def test_b2_length_mismatch_against_matrix_fails():
    with pytest.raises(ckpt.CellIdentityError, match="but the matrix has"):
        ckpt.validate_cell_ids(["a", "b"], 3, cohort="t", channel="A")


def test_b2_identity_is_verbatim_not_normalized():
    """Case/whitespace variants are REPORTED, never silently merged, and valid
    IDs come back byte-identical."""
    out = ckpt.validate_cell_ids(["  a  ", "B"], 2)
    assert out == ["  a  ", "B"], "no trim, no case fold"


def test_b2_none_and_nan_are_rejected():
    with pytest.raises(ckpt.CellIdentityError, match="is missing"):
        ckpt.validate_cell_ids(["a", None], 2)
    with pytest.raises(ckpt.CellIdentityError, match="is missing"):
        ckpt.validate_cell_ids(["a", float("nan")], 2)


# ---- M1: recorded Scrublet params == params actually passed ---------------


def test_m1_every_recorded_param_reaches_the_call():
    core.score_one_capture(
        make_matrix(),
        cell_ids(),
        expected_doublet_rate=0.05,
        seed=1,
        scrublet_fn=fake_scrublet_success,
    )
    assert (
        RECEIVED_KWARGS == core.SCRUBLET_PARAMS
    ), "the fingerprint records SCRUBLET_PARAMS, so every one of them must be passed"


def test_m1_params_are_real_scanpy_arguments():
    """Guards against re-introducing a parameter scanpy does not accept, which
    would make the fingerprint bind a value that changes nothing."""
    import inspect

    try:
        import scanpy as sc
    except Exception:
        pytest.skip("scanpy not installed")
    accepted = set(inspect.signature(sc.pp.scrublet).parameters)
    assert (
        set(core.SCRUBLET_PARAMS) <= accepted
    ), f"not real scanpy args: {set(core.SCRUBLET_PARAMS) - accepted}"


def test_m1_param_change_invalidates_checkpoint(tmp_path):
    _save(tmp_path, _fp(scrublet_params={"sim_doublet_ratio": 2.0}))
    with pytest.raises(ckpt.CheckpointMismatch, match="fingerprint mismatch"):
        ckpt.load_checkpoint(
            str(tmp_path),
            "k",
            _fp(scrublet_params={"sim_doublet_ratio": 3.0}),
            expected_cell_ids=["a", "b", "c"],
        )


# ---- M2: strict scores validation on restore ------------------------------


def _save_raw(
    d,
    fp,
    key="k",
    *,
    scores,
    call_status=None,
    candidate=None,
    cell_ids=("a", "b", "c"),
    capture_status=st.CAPTURE_SCORED,
    diagnostics=None,
):
    """Write a checkpoint directly, bypassing the normal path, to simulate a
    corrupt or hand-edited file."""
    import json as _json

    meta = _json.dumps(
        {
            "fingerprint": fp,
            "fingerprint_digest": ckpt._fingerprint_digest(fp),
            "capture_status": capture_status,
            "diagnostics": diagnostics if diagnostics is not None else _conforming_diagnostics(),
        },
        sort_keys=True,
        default=str,
    )
    os.makedirs(d, exist_ok=True)
    np.savez(
        os.path.join(str(d), f"{key}.npz"),
        meta=np.array(meta),
        scores=scores,
        call_status=np.asarray(call_status or [st.CALL_SCRUBLET_NONCANDIDATE] * 3, dtype="U24"),
        candidate=np.asarray(candidate if candidate is not None else [0, 0, 0], dtype=np.int8),
        cell_ids=np.asarray([str(c) for c in cell_ids], dtype="U64"),
    )


@pytest.mark.parametrize(
    "scores,pattern",
    [
        (np.array([0.1, np.nan, 0.3]), "non-finite"),
        (np.array([0.1, np.inf, 0.3]), "non-finite"),
        (np.array([[0.1, 0.2, 0.3]]), "ndim"),
        (np.array([True, False, True]), "not numeric"),
        (np.array(["a", "b", "c"]), "not numeric"),
    ],
)
def test_m2_invalid_scores_rejected_on_restore(tmp_path, scores, pattern):
    _save_raw(tmp_path, _fp(), scores=scores)
    with pytest.raises(ckpt.CheckpointMismatch, match=pattern):
        ckpt.load_checkpoint(str(tmp_path), "k", _fp(), expected_cell_ids=["a", "b", "c"])


def test_m2_scored_capture_without_finite_threshold_rejected(tmp_path):
    _save_raw(
        tmp_path,
        _fp(),
        scores=np.array([0.1, 0.2, 0.3]),
        diagnostics=_conforming_diagnostics(threshold=None),
    )
    # now raised through the unified capture contract rather than a lone check
    with pytest.raises(ckpt.CheckpointMismatch, match="requires a finite threshold"):
        ckpt.load_checkpoint(str(tmp_path), "k", _fp(), expected_cell_ids=["a", "b", "c"])


def test_m2_nan_scores_allowed_for_non_scored_capture(tmp_path):
    """too_small captures legitimately carry all-NaN scores."""
    _save_raw(
        tmp_path,
        _fp(),
        scores=np.full(3, np.nan),
        call_status=[st.CALL_TOO_SMALL] * 3,
        candidate=[-1, -1, -1],
        capture_status=st.CAPTURE_TOO_SMALL,
        diagnostics=_conforming_diagnostics(threshold=None),
    )
    got = ckpt.load_checkpoint(str(tmp_path), "k", _fp(), expected_cell_ids=["a", "b", "c"])
    assert got is not None and got["capture_status"] == st.CAPTURE_TOO_SMALL


# ---- M4: durability -------------------------------------------------------


def test_m4_failed_save_leaves_no_temp_file(tmp_path, monkeypatch):
    monkeypatch.setattr(ckpt.os, "replace", lambda *a, **k: (_ for _ in ()).throw(OSError("boom")))
    with pytest.raises(OSError):
        _save(tmp_path, _fp())
    leftovers = [f for f in os.listdir(str(tmp_path)) if f.endswith(".tmp")]
    assert leftovers == [], f"own temp must be cleaned up, found {leftovers}"
    assert not os.path.exists(os.path.join(str(tmp_path), "k.npz"))


def test_m4_foreign_temp_file_is_not_deleted(tmp_path):
    foreign = tmp_path / ".someone-else.99999.npz.tmp"
    tmp_path.mkdir(exist_ok=True)
    foreign.write_bytes(b"not ours")
    _save(tmp_path, _fp())
    assert foreign.exists(), "another writer's temp file is not ours to remove"


# ---- extended fresh vs resume --------------------------------------------


def test_ext_fresh_vs_resume_full_diagnostics_schema(tmp_path):
    inp = write_combined_cohort(tmp_path / "in", "t", ("A", "B"), n=240)
    cfg = write_config(tmp_path / "c.csv", [config_row("t", "A"), config_row("t", "B")])
    out = tmp_path / "out"
    r1 = runner.run("t", "combined", inp, str(out), cfg, 7, scrublet_fn=fake_scrublet_success)
    r2 = runner.run("t", "combined", inp, str(out), cfg, 7, scrublet_fn=fake_scrublet_success)
    assert r2["restored"] == 2

    def rows(p):
        with open(p) as fh:
            rd = csv.DictReader(fh)
            return rd.fieldnames, list(rd)

    f1, c1 = rows(r1["calls_csv"])
    f2, c2 = rows(r2["calls_csv"])
    assert f1 == f2, "column order must be identical"
    assert [r["cell_id"] for r in c1] == [r["cell_id"] for r in c2], "cell output order identical"
    for x, y in zip(c1, c2):
        assert x == y, "every field of every cell row must match"

    g1, d1 = rows(r1["diagnostics_csv"])
    g2, d2 = rows(r2["diagnostics_csv"])
    assert g1 == g2, "diagnostics schema and column order identical"
    assert [r["channel_id"] for r in d1] == [r["channel_id"] for r in d2], "capture order identical"
    ignore = {"restored_from_checkpoint"}
    for x, y in zip(d1, d2):
        assert {k: v for k, v in x.items() if k not in ignore} == {
            k: v for k, v in y.items() if k not in ignore
        }


def test_ext_reordered_input_invalidates_old_checkpoint(tmp_path):
    """(b) A real order perturbation -- same cells, permuted rows AND matrix
    columns -- must refuse the stored checkpoint rather than reuse it."""
    inp = tmp_path / "in"
    write_combined_cohort(inp, "t", ("A",))
    cfg = write_config(tmp_path / "c.csv", [config_row("t", "A")])
    out = tmp_path / "out"
    runner.run("t", "combined", str(inp), str(out), cfg, 7, scrublet_fn=fake_scrublet_success)
    write_combined_cohort(inp, "t", ("A",), reorder=True)  # same cells, new order
    with pytest.raises(ckpt.CheckpointMismatch, match="fingerprint mismatch"):
        runner.run("t", "combined", str(inp), str(out), cfg, 7, scrublet_fn=fake_scrublet_success)


def test_ext_reordered_fresh_run_scientific_result_is_order_dependent(tmp_path):
    """(c) Reported, not asserted away.

    With a clean output directory, a reordered input scores fresh. Joining by
    cell_id, the scores are NOT guaranteed identical: the detector's own
    neighbourhood construction and the injected fake are both position-
    dependent. This test RECORDS that fact rather than claiming Scrublet is
    order-insensitive -- a claim nothing here establishes.
    """
    cfg = write_config(tmp_path / "c.csv", [config_row("t", "A")])
    a_in = write_combined_cohort(tmp_path / "a", "t", ("A",))
    b_in = write_combined_cohort(tmp_path / "b", "t", ("A",), reorder=True)
    ra = runner.run(
        "t", "combined", a_in, str(tmp_path / "oa"), cfg, 7, scrublet_fn=fake_scrublet_success
    )
    rb = runner.run(
        "t", "combined", b_in, str(tmp_path / "ob"), cfg, 7, scrublet_fn=fake_scrublet_success
    )
    A = {r["cell_id"]: r for r in csv.DictReader(open(ra["calls_csv"]))}
    B = {r["cell_id"]: r for r in csv.DictReader(open(rb["calls_csv"]))}
    assert A.keys() == B.keys(), "the same cells must be present either way"
    n_diff = sum(1 for cid in A if A[cid]["scrublet_score"] != B[cid]["scrublet_score"])
    # The mechanical guarantee is identity, not numerical invariance.
    assert n_diff >= 0
    print(
        f"\n[order sensitivity] {n_diff}/{len(A)} cells differ in score after reordering; "
        "order invariance is NOT a property this layer provides or claims."
    )


# ---- candidate export -----------------------------------------------------


def test_export_candidate_is_nullable_int_never_true_false(tmp_path):
    inp = write_combined_cohort(tmp_path / "in", "t", ("A",))
    cfg = write_config(tmp_path / "c.csv", [config_row("t", "A")])
    r = runner.run(
        "t", "combined", inp, str(tmp_path / "out"), cfg, 7, scrublet_fn=fake_scrublet_nan_threshold
    )
    vals = {row["candidate_doublet"] for row in csv.DictReader(open(r["calls_csv"]))}
    assert vals == {""}, f"threshold_failed must export blank (unknown), got {vals}"
    assert "False" not in vals and "True" not in vals


def test_export_round_trip_unknown_does_not_become_false(tmp_path):
    inp = write_combined_cohort(tmp_path / "in", "t", ("A",))
    cfg = write_config(tmp_path / "c.csv", [config_row("t", "A")])
    r = runner.run(
        "t", "combined", inp, str(tmp_path / "out"), cfg, 7, scrublet_fn=fake_scrublet_nan_threshold
    )
    rows_ = list(csv.DictReader(open(r["calls_csv"])))
    assert all(x["call_status"] == st.CALL_THRESHOLD_FAILED for x in rows_)
    assert st.mechanical_candidate_mask(x["call_status"] for x in rows_) == [False] * len(rows_)
    assert (
        sum(1 for x in rows_ if x["candidate_doublet"] == "0") == 0
    ), 'unknown must not be exported as 0, which would read back as "not a doublet"'


# ==========================================================================
# Round 3 — second-opinion review: B3, B4, M3, M5, M6, M7, m1, strict bool
# ==========================================================================

# ---- B3: status combinations enforced on fresh, save AND restore ----------

CONTRADICTIONS = [
    # (capture_status, call_status, candidate, scores, threshold, why)
    (
        st.CAPTURE_SCORED,
        [st.CALL_TOO_SMALL] * 3,
        [-1, -1, -1],
        [0.1, 0.2, 0.3],
        0.9,
        "scored but every cell too_small",
    ),
    (
        st.CAPTURE_SCORED,
        [st.CALL_SCRUBLET_NONCANDIDATE] * 3,
        [-1, -1, -1],
        [0.1, 0.2, 0.3],
        0.9,
        "scored with unknown candidate",
    ),
    (
        st.CAPTURE_SCORED,
        [st.CALL_PREDICTED_DOUBLET, st.CALL_SCRUBLET_NONCANDIDATE, st.CALL_SCRUBLET_NONCANDIDATE],
        [0, 0, 0],
        [0.1, 0.2, 0.3],
        0.9,
        "predicted_doublet with candidate 0",
    ),
    (
        st.CAPTURE_SCORED,
        [st.CALL_SCRUBLET_NONCANDIDATE] * 3,
        [1, 0, 0],
        [0.1, 0.2, 0.3],
        0.9,
        "scrublet_noncandidate with candidate 1",
    ),
    (
        st.CAPTURE_THRESHOLD_FAILED,
        [st.CALL_THRESHOLD_FAILED, st.CALL_PREDICTED_DOUBLET, st.CALL_THRESHOLD_FAILED],
        [-1, 1, -1],
        [0.1, 0.2, 0.3],
        None,
        "threshold_failed retaining a predicted_doublet",
    ),
    (
        st.CAPTURE_TOO_SMALL,
        [st.CALL_SCRUBLET_NONCANDIDATE] * 3,
        [0, 0, 0],
        [np.nan] * 3,
        None,
        "too_small yielding singlets",
    ),
    (
        st.CAPTURE_RUNTIME_FAILED,
        [st.CALL_SCRUBLET_NONCANDIDATE] * 3,
        [0, 0, 0],
        [np.nan] * 3,
        None,
        "runtime_failed yielding singlets",
    ),
    (
        st.CAPTURE_TOO_SMALL,
        [st.CALL_TOO_SMALL] * 3,
        [-1, -1, -1],
        [0.1, 0.2, 0.3],
        None,
        "non-scored status carrying finite scores",
    ),
    (
        st.CAPTURE_THRESHOLD_FAILED,
        [st.CALL_THRESHOLD_FAILED] * 3,
        [-1, -1, -1],
        [0.1, 0.2, 0.3],
        0.9,
        "non-scored status carrying a finite threshold",
    ),
    (
        st.CAPTURE_THRESHOLD_FAILED,
        [st.CALL_THRESHOLD_FAILED] * 3,
        [-1, -1, -1],
        [0.1, np.nan, 0.3],
        None,
        "scores a mixture of finite and non-finite",
    ),
]


@pytest.mark.parametrize("cap,calls,cand,scores,thr,why", CONTRADICTIONS)
def test_b3_contradictory_combinations_rejected(cap, calls, cand, scores, thr, why):
    with pytest.raises(st.CaptureInvariantError):
        st.validate_capture_invariants(
            cap,
            calls,
            np.asarray(cand, dtype=np.int8),
            np.asarray(scores, dtype=float),
            threshold=thr,
            where=why,
        )


@pytest.mark.parametrize("cap,calls,cand,scores,thr,why", CONTRADICTIONS)
def test_b3_contradictions_also_rejected_on_restore(tmp_path, cap, calls, cand, scores, thr, why):
    """Same table, through the checkpoint restore path -- fresh, save and load
    must enforce one contract, not three partial ones."""
    _save_raw(
        tmp_path,
        _fp(),
        scores=np.asarray(scores, dtype=float),
        call_status=calls,
        candidate=cand,
        capture_status=cap,
        diagnostics=_conforming_diagnostics(threshold=thr),
    )
    with pytest.raises(ckpt.CheckpointMismatch, match="violates the status contract"):
        ckpt.load_checkpoint(str(tmp_path), "k", _fp(), expected_cell_ids=["a", "b", "c"])


def test_b3_save_refuses_contradiction(tmp_path):
    with pytest.raises(st.CaptureInvariantError):
        ckpt.save_checkpoint(
            str(tmp_path),
            "k",
            _fp(),
            scores=[0.1, 0.2, 0.3],
            call_status=[st.CALL_SCRUBLET_NONCANDIDATE] * 3,
            candidate=np.array([1, 1, 1], dtype=np.int8),
            cell_ids=["a", "b", "c"],
            capture_status=st.CAPTURE_SCORED,
            diagnostics=_conforming_diagnostics(),
        )


def test_b3_diagnostics_count_must_match_the_array():
    with pytest.raises(st.CaptureInvariantError, match="predicted doublets"):
        st.validate_capture_invariants(
            st.CAPTURE_SCORED,
            [st.CALL_SCRUBLET_NONCANDIDATE] * 3,
            np.zeros(3, dtype=np.int8),
            np.array([0.1, 0.2, 0.3]),
            threshold=0.9,
            diagnostics={"n_predicted_doublets": 2},
        )


def test_b3_valid_combinations_pass():
    st.validate_capture_invariants(
        st.CAPTURE_SCORED,
        [st.CALL_PREDICTED_DOUBLET, st.CALL_SCRUBLET_NONCANDIDATE],
        np.array([1, 0], dtype=np.int8),
        np.array([0.95, 0.1]),
        threshold=0.9,
        diagnostics={"n_predicted_doublets": 1},
    )
    st.validate_capture_invariants(
        st.CAPTURE_TOO_SMALL,
        [st.CALL_TOO_SMALL] * 2,
        np.array([-1, -1], dtype=np.int8),
        np.array([np.nan, np.nan]),
        threshold=None,
    )


# ---- B4: raw-count contract ----------------------------------------------


@pytest.mark.parametrize(
    "bad,pattern",
    [
        ("negative", "negative values"),
        ("fractional", "non-integer"),
        ("nan", "NaN or inf"),
        ("inf", "NaN or inf"),
    ],
)
def test_b4_non_raw_counts_rejected(bad, pattern):
    m = make_matrix().astype(float).tolil()
    if bad == "negative":
        m[0, 0] = -1.0
    elif bad == "fractional":
        m[0, 0] = 1.5  # log1p / normalized / DecontX-corrected
    elif bad == "nan":
        m[0, 0] = np.nan
    else:
        m[0, 0] = np.inf
    with pytest.raises(core.RawCountError, match=pattern):
        core.validate_raw_counts(m.tocsr())


def test_b4_integer_valued_float_matrix_is_accepted():
    """dtype is not the test -- counts are commonly stored as float."""
    core.validate_raw_counts(make_matrix().astype(float))


def test_b4_scoring_refuses_non_raw_input():
    m = make_matrix().astype(float).tolil()
    m[0, 0] = 0.37
    with pytest.raises(core.RawCountError):
        core.score_one_capture(
            m.tocsr(),
            cell_ids(),
            expected_doublet_rate=0.05,
            seed=1,
            scrublet_fn=fake_scrublet_success,
        )


# ---- strict boolean predicted_doublet -------------------------------------


def test_strict_bool_non_boolean_predictions_are_threshold_failed():
    def fake_float_preds(adata, expected_doublet_rate=None, random_state=0, **kw):
        n = adata.n_obs
        adata.obs["doublet_score"] = np.linspace(0, 1, n)
        # 0.5 and 2 would both become True under a bare .astype(bool)
        adata.obs["predicted_doublet"] = np.full(n, 0.5)
        adata.uns["scrublet"] = {"threshold": 0.9}

    res = core.score_one_capture(
        make_matrix(), cell_ids(), expected_doublet_rate=0.05, seed=1, scrublet_fn=fake_float_preds
    )
    assert res.capture_status == st.CAPTURE_THRESHOLD_FAILED
    assert any("strict boolean" in p for p in res.diagnostics["problems"])
    # the threshold Scrublet reported is kept as EVIDENCE but is not operative,
    # so a capture that failed for a non-threshold reason cannot look scored
    assert res.diagnostics["threshold"] is None
    assert res.diagnostics["threshold_observed"] == 0.9


def test_strict_bool_int_zero_one_is_accepted():
    def fake_int_preds(adata, expected_doublet_rate=None, random_state=0, **kw):
        n = adata.n_obs
        sc_ = np.linspace(0, 1, n)
        adata.obs["doublet_score"] = sc_
        adata.obs["predicted_doublet"] = (sc_ > 0.9).astype(np.int8)
        adata.uns["scrublet"] = {"threshold": 0.9, "doublet_scores_sim": sc_}

    res = core.score_one_capture(
        make_matrix(), cell_ids(), expected_doublet_rate=0.05, seed=1, scrublet_fn=fake_int_preds
    )
    assert res.capture_status == st.CAPTURE_SCORED


# ---- M5: temp files are unique per call ----------------------------------


def test_m5_temp_names_unique_within_one_process(tmp_path, monkeypatch):
    """Two saves in ONE process must not share a temp path -- PID alone could
    not distinguish threads, retries or nested calls."""
    seen = []
    real_replace = ckpt.os.replace

    def spy(src, dst):
        seen.append(str(src))
        return real_replace(src, dst)

    monkeypatch.setattr(ckpt.os, "replace", spy)
    _save(tmp_path, _fp())
    _save(tmp_path, _fp())
    assert (
        len(seen) == 2 and seen[0] != seen[1]
    ), f"same-process saves must not share a temp path: {seen}"


# ---- M7: composite identity ----------------------------------------------


def test_m7_output_carries_a_composite_identity_key(tmp_path):
    inp = write_combined_cohort(tmp_path / "in", "t", ("A", "B"), n=240)
    cfg = write_config(tmp_path / "c.csv", [config_row("t", "A"), config_row("t", "B")])
    r = runner.run(
        "t", "combined", inp, str(tmp_path / "out"), cfg, 7, scrublet_fn=fake_scrublet_success
    )
    rows_ = list(csv.DictReader(open(r["calls_csv"])))
    assert "gate2_cell_key" in rows_[0]
    keys = [x["gate2_cell_key"] for x in rows_]
    assert len(keys) == len(set(keys)), "composite key must be globally unique"
    for x in rows_[:5]:
        import json as _json

        assert x["gate2_cell_key"] == _json.dumps(
            ["t", x["channel_id"], x["cell_id"]], ensure_ascii=False, separators=(",", ":")
        )


def test_m7_same_barcode_in_two_captures_stays_distinct(tmp_path):
    """Barcodes legitimately repeat across captures; the composite key must keep
    them apart where cell_id alone would merge them."""
    from scipy.io import mmwrite

    root = tmp_path / "in" / "split_per_channel"
    for ch in ("A", "B"):
        d = root / ch
        d.mkdir(parents=True)
        mmwrite(str(d / "counts.mtx"), make_matrix().T.tocsc())
        with open(d / "cells.csv", "w", newline="") as fh:
            w = csv.writer(fh, lineterminator="\n")
            w.writerow(["cell_id", "channel_id", "donor_id"])
            for i in range(N_CELLS):
                w.writerow([f"AAACCTGAGAAACCAT-{i}", ch, "d1"])  # identical barcodes
    cfg = write_config(tmp_path / "c.csv", [config_row("t", "A"), config_row("t", "B")])
    r = runner.run(
        "t",
        "split",
        str(tmp_path / "in"),
        str(tmp_path / "out"),
        cfg,
        7,
        scrublet_fn=fake_scrublet_success,
    )
    rows_ = list(csv.DictReader(open(r["calls_csv"])))
    assert len({x["cell_id"] for x in rows_}) == N_CELLS, "cell_id alone collides, as expected"
    # proven on the AUTHORITATIVE three-column tuple, not on the derived string
    triples = {(x["cohort"], x["channel_id"], x["cell_id"]) for x in rows_}
    assert len(triples) == 2 * N_CELLS, "the (cohort, channel_id, cell_id) tuple must be unique"
    assert (
        len({x["gate2_cell_key"] for x in rows_}) == 2 * N_CELLS
    ), "the derived key must agree with the authoritative tuple"


# ---- m1: git state distinguishes the three cases --------------------------


def test_m1_git_state_distinguishes_dirty_unavailable_and_non_repo(monkeypatch):
    import subprocess

    class R:
        def __init__(self, rc, out):
            self.returncode, self.stdout = rc, out

    def fake_run(args, **kw):
        if args[:2] == ["git", "rev-parse"]:
            return R(0, "abc123\n")
        return R(0, " M some/file.py\n")

    monkeypatch.setattr(subprocess, "run", fake_run)
    assert runner._git_state() == "abc123+dirty"

    monkeypatch.setattr(
        subprocess,
        "run",
        lambda args, **kw: R(0, "abc123\n") if args[:2] == ["git", "rev-parse"] else R(0, ""),
    )
    assert runner._git_state() == "abc123"

    monkeypatch.setattr(subprocess, "run", lambda args, **kw: R(128, ""))
    assert runner._git_state() == "not-a-git-repository"

    def boom(args, **kw):
        raise FileNotFoundError()

    monkeypatch.setattr(subprocess, "run", boom)
    assert runner._git_state() == "git-unavailable"


# ---- M3: the ambiguous artifact is no longer produced ---------------------


def test_m3_materialize_no_longer_writes_prior_counts():
    src = open(os.path.join(CODE, "materialize_liu_wilk_for_gate2.py")).read()
    assert "to_csv(f'{out}/prior_counts.csv'" not in src
    assert "NO LONGER WRITTEN" in src


# ==========================================================================
# Round 4 — third-audit findings: NA crash, identity escaping, U64, bounds
# ==========================================================================

# ---- MAJOR-1: nothing escapes the validator ------------------------------


def _fake_with_preds(vals):
    def f(adata, expected_doublet_rate=None, random_state=0, **kw):
        n = adata.n_obs
        adata.obs["doublet_score"] = np.linspace(0, 1, n)
        adata.obs["predicted_doublet"] = vals(n)
        adata.uns["scrublet"] = {"threshold": 0.9, "doublet_scores_sim": np.linspace(0, 1, 2 * n)}

    return f


def test_r4_nullable_boolean_with_na_is_classified_not_crashed():
    """The exact input that made np.unique raise 'boolean value of NA is
    ambiguous' from inside the error-message construction, escaping the
    validator entirely."""
    pd = pytest.importorskip("pandas")
    res = core.score_one_capture(
        make_matrix(),
        cell_ids(),
        expected_doublet_rate=0.05,
        seed=1,
        scrublet_fn=_fake_with_preds(lambda n: pd.array([True, None] * (n // 2), dtype="boolean")),
    )
    assert res.capture_status == st.CAPTURE_THRESHOLD_FAILED
    assert set(res.call_status) == {st.CALL_THRESHOLD_FAILED}
    assert set(np.unique(res.candidate)) == {-1}
    assert res.diagnostics["problems"]


@pytest.mark.parametrize(
    "vals,label",
    [
        (lambda n: np.full(n, 0.5), "float 0.5"),
        (lambda n: (np.arange(n) % 2).astype(float), "float 0.0/1.0"),
        (lambda n: np.array(["True", "False"] * 60), "string"),
        (lambda n: np.full(n, 2, dtype=np.int64), "integer 2"),
        (lambda n: np.array([None] * n, dtype=object), "object with None"),
    ],
)
def test_r4_illegal_predictions_all_become_threshold_failed(vals, label):
    res = core.score_one_capture(
        make_matrix(),
        cell_ids(),
        expected_doublet_rate=0.05,
        seed=1,
        scrublet_fn=_fake_with_preds(vals),
    )
    assert res.capture_status == st.CAPTURE_THRESHOLD_FAILED, label
    assert set(res.call_status) == {st.CALL_THRESHOLD_FAILED}
    assert set(np.unique(res.candidate)) == {-1}
    assert res.diagnostics["threshold"] is None, "no operative threshold on a failed capture"


def test_r4_two_dimensional_predictions_refused_by_the_extractor():
    """A 2-D array cannot be assigned into adata.obs at all, so it never reaches
    the extractor through a real call -- the extractor is tested directly."""
    preds, problem = core._extract_predictions(np.zeros((10, 2), dtype=bool), 10)
    assert preds is None and "not 1-D" in problem


def test_r4_categorical_predictions_refused():
    pd = pytest.importorskip("pandas")
    res = core.score_one_capture(
        make_matrix(),
        cell_ids(),
        expected_doublet_rate=0.05,
        seed=1,
        scrublet_fn=_fake_with_preds(lambda n: pd.Categorical(["a", "b"] * (n // 2))),
    )
    assert res.capture_status == st.CAPTURE_THRESHOLD_FAILED


def test_r4_validator_exception_becomes_runtime_failed(monkeypatch):
    """Last-resort boundary: if the validator itself fails unexpectedly, the
    capture is runtime_failed with every cell unknown -- never a call."""

    def boom(adata, n_cells):
        raise RuntimeError("synthetic validator explosion")

    monkeypatch.setattr(core, "_validate_scrublet_output", boom)
    res = core.score_one_capture(
        make_matrix(),
        cell_ids(),
        expected_doublet_rate=0.05,
        seed=1,
        scrublet_fn=fake_scrublet_success,
    )
    assert res.capture_status == st.CAPTURE_RUNTIME_FAILED
    assert set(res.call_status) == {st.CALL_RUNTIME_FAILED}
    assert set(np.unique(res.candidate)) == {-1}
    assert any("validator failed" in p for p in res.diagnostics["problems"])


def test_r4_error_messages_never_touch_unvalidated_values():
    """np.unique on an unvalidated array is itself a crash path (pd.NA). Only
    CODE is checked -- the docstrings deliberately name it to explain why."""
    import ast

    src = open(os.path.join(CODE, "gate2_score_core.py")).read()
    for fn in ast.walk(ast.parse(src)):
        if not isinstance(fn, ast.FunctionDef):
            continue
        if fn.name not in ("_safe_repr_array", "_extract_predictions"):
            continue
        body = fn.body
        if body and isinstance(body[0], ast.Expr) and isinstance(body[0].value, ast.Constant):
            body = body[1:]
        code = chr(10).join(ast.unparse(n) for n in body)
        assert "unique" not in code, fn.name + " must not sort or uniq an unvalidated array"


# ---- MAJOR-2: identity ----------------------------------------------------


@pytest.mark.parametrize(
    "cohort,ch,cid",
    [
        ("t", "A|B", "c1"),
        ("t", "A", "B|c1"),
        ("t", "A,B", "c1"),
        ("t", 'A"B', "c1"),
        ("t", "A B", "c1"),
        ("t", "Aβ細胞", "c1"),
        ("t", "", "c1"),
    ],
)
def test_r4_derived_key_is_injective(cohort, ch, cid):
    """A '|' join is not injective: ('A|B','c1') and ('A','B|c1') collide. The
    canonical JSON form must keep every pair distinct."""
    import json as _json

    k = runner.derived_cell_key(cohort, ch, cid)
    assert _json.loads(k) == [cohort, ch, cid]


def test_r4_pipe_containing_fields_do_not_collide():
    a = runner.derived_cell_key("t", "A|B", "c1")
    b = runner.derived_cell_key("t", "A", "B|c1")
    assert a != b, "the unescaped join collided here; canonical JSON must not"


def test_r4_no_consumer_derives_deletion_from_the_key():
    for f in ("gate2_run.py", "gate2_status.py", "gate2_checkpoint.py", "gate2_score_core.py"):
        src = open(os.path.join(CODE, f)).read()
        for line in src.splitlines():
            if "gate2_cell_key" in line and "delet" in line.lower():
                raise AssertionError(f"{f}: deletion must not depend on the derived key: {line}")


# ---- MAJOR-3: no fixed-width truncation ----------------------------------


@pytest.mark.parametrize("length", [63, 64, 65, 256])
def test_r4_long_cell_ids_survive_save_and_restore(tmp_path, length):
    ids = [("x" * (length - 3)) + f"{i:03d}" for i in range(3)]
    assert all(len(i) == length for i in ids)
    fp = _fp(cell_id_checksum=ckpt.ordered_cell_checksum(ids))
    ckpt.save_checkpoint(
        str(tmp_path),
        "k",
        fp,
        scores=[0.1, 0.2, 0.3],
        call_status=[st.CALL_SCRUBLET_NONCANDIDATE] * 3,
        candidate=np.zeros(3, dtype=np.int8),
        cell_ids=ids,
        capture_status=st.CAPTURE_SCORED,
        diagnostics=_conforming_diagnostics(),
    )
    got = ckpt.load_checkpoint(str(tmp_path), "k", fp, expected_cell_ids=ids)
    assert got is not None and got["cell_ids"] == ids


def test_r4_unicode_cell_ids_survive(tmp_path):
    ids = ["細胞-α-1", "cellule-β-2", "ячейка-γ-3"]
    fp = _fp(cell_id_checksum=ckpt.ordered_cell_checksum(ids))
    ckpt.save_checkpoint(
        str(tmp_path),
        "k",
        fp,
        scores=[0.1, 0.2, 0.3],
        call_status=[st.CALL_SCRUBLET_NONCANDIDATE] * 3,
        candidate=np.zeros(3, dtype=np.int8),
        cell_ids=ids,
        capture_status=st.CAPTURE_SCORED,
        diagnostics=_conforming_diagnostics(),
    )
    assert ckpt.load_checkpoint(str(tmp_path), "k", fp, expected_cell_ids=ids)["cell_ids"] == ids


def test_r4_schema_version_bumped_and_v1_rejected(tmp_path):
    # Assert the durable property, not a literal: whatever the current version
    # is, it must not be /1, and /1 must be in the rejected set. Pinning the
    # literal here only created churn on every legitimate bump.
    assert ckpt.CHECKPOINT_SCHEMA_VERSION != "gate2-checkpoint/1"
    assert "gate2-checkpoint/1" in ckpt.REJECTED_SCHEMA_VERSIONS
    _save_raw(
        tmp_path,
        _fp(checkpoint_schema_version="gate2-checkpoint/1"),
        scores=np.array([0.1, 0.2, 0.3]),
    )
    with pytest.raises(ckpt.CheckpointMismatch, match="checkpoint_schema_version"):
        ckpt.load_checkpoint(str(tmp_path), "k", _fp(), expected_cell_ids=["a", "b", "c"])


# ---- MINOR-1/2: tightened invariants --------------------------------------


def test_r4_zero_cell_capture_refused():
    with pytest.raises(st.CaptureInvariantError, match="zero cells"):
        st.validate_capture_invariants(
            st.CAPTURE_SCORED, [], np.array([], dtype=np.int8), np.array([]), threshold=0.9
        )


@pytest.mark.parametrize(
    "cand,why",
    [
        (np.array(["0", "1", "0"]), "string array"),
        (np.array([0.0, 1.0, 0.0]), "float array"),
        (np.array([False, True, False]), "bool array"),
        (np.zeros((3, 1), dtype=np.int8), "2-D int8"),
        ([0, 1, 0], "python list"),
    ],
)
def test_r4_candidate_must_be_strict_int8_1d(cand, why):
    with pytest.raises(st.CaptureInvariantError):
        st.validate_capture_invariants(
            st.CAPTURE_SCORED,
            [
                st.CALL_SCRUBLET_NONCANDIDATE,
                st.CALL_PREDICTED_DOUBLET,
                st.CALL_SCRUBLET_NONCANDIDATE,
            ],
            cand,
            np.array([0.1, 0.9, 0.2]),
            threshold=0.5,
            diagnostics={"n_predicted_doublets": 1},
        )


# ---- raw-count numeric bounds ---------------------------------------------


@pytest.mark.parametrize(
    "arr,pattern",
    [
        (np.array([[True, False]]), "boolean"),
        (np.array([[1 + 0j, 2 + 0j]]), "complex"),
        (np.array([[float(2**53 + 2), 1.0]]), "above 2\\*\\*53"),
        (np.array([[1e300, 1.0]]), "above 2\\*\\*53"),
        (np.array([[-1, 2]], dtype=np.int64), "negative"),
    ],
)
def test_r4_raw_count_numeric_bounds(arr, pattern):
    with pytest.raises(core.RawCountError, match=pattern):
        core.validate_raw_counts(arr)


def test_r4_all_zero_capture_refused_as_unrunnable_not_invalid():
    from scipy import sparse as _sp

    with pytest.raises(core.RawCountError) as e:
        core.validate_raw_counts(_sp.csr_matrix((5, 5)))
    msg = str(e.value)
    assert "runnability limit" in msg, "an all-zero matrix is legal, just not scoreable"


def test_r4_value_just_below_2_53_accepted():
    core.validate_raw_counts(np.array([[float(2**53 - 2), 1.0]]))


# ---- checkpoint concurrency -----------------------------------------------


def test_r4_two_writers_same_key_do_not_share_or_delete_temps(tmp_path, monkeypatch):
    """A failing writer must delete only its own temp, even when a second
    writer is mid-flight on the SAME capture key."""
    temps = []
    real_replace = ckpt.os.replace

    def capture_then_fail(src, dst):
        temps.append(str(src))
        if len(temps) == 1:
            return real_replace(src, dst)
        raise OSError("second writer fails after the first succeeded")

    monkeypatch.setattr(ckpt.os, "replace", capture_then_fail)
    _save(tmp_path, _fp())
    with pytest.raises(OSError):
        _save(tmp_path, _fp())
    assert temps[0] != temps[1], "two writers must not share a temp path"
    assert os.path.exists(
        os.path.join(str(tmp_path), "k.npz")
    ), "the first writer's checkpoint must survive the second writer's failure"
    assert not os.path.exists(temps[1]), "the failing writer removes its own temp"


def test_r4_dir_fsync_failure_keeps_the_checkpoint(tmp_path, monkeypatch):
    """replace succeeded, so the file exists and its contents are durable.
    Raising here would tell the caller nothing was written when something was."""
    real_fsync = ckpt.os.fsync
    state = {"n": 0}

    def flaky(fd):
        state["n"] += 1
        if state["n"] >= 2:  # the directory fsync, after the file's
            raise OSError("synthetic directory fsync failure")
        return real_fsync(fd)

    monkeypatch.setattr(ckpt.os, "fsync", flaky)
    path = _save(tmp_path, _fp())
    assert os.path.exists(path), "the checkpoint stands"
    assert (
        ckpt.load_checkpoint(str(tmp_path), "k", _fp(), expected_cell_ids=["a", "b", "c"])
        is not None
    )


# ==========================================================================
# Round 5 — two-layer decision model: mechanical candidate != authorization
# ==========================================================================

GATE2_MODULES = (
    "gate2_status.py",
    "gate2_run.py",
    "gate2_score_core.py",
    "gate2_checkpoint.py",
    "gate2_runtime_config.py",
)


def test_r5_old_deletion_vocabulary_is_gone_from_active_code():
    """A breaking rename, not an alias: a caller using the old names must fail
    loudly rather than quietly receive a candidate mask it may read as
    permission."""
    assert not hasattr(st, "DELETION_ELIGIBLE")
    assert not hasattr(st, "deletion_mask")
    assert hasattr(st, "MECHANICAL_DOUBLET_CANDIDATE_STATUSES")
    assert hasattr(st, "mechanical_candidate_mask")
    for mod in GATE2_MODULES:
        code = open(os.path.join(CODE, mod)).read()
        for line in code.splitlines():
            stripped = line.strip()
            if stripped.startswith("#") or '"""' in line or "previous name" in line:
                continue  # docstrings may record the rename
            for old in ("DELETION_ELIGIBLE", "deletion_mask", "deletion_eligible"):
                assert old not in line, f"{mod}: stale name {old!r} in: {line.strip()[:80]}"


def test_r5_no_final_deletion_authorization_artifact_exists():
    """The accurate representation of the current state is ABSENCE. Not even an
    all-false column, which would later read as a completed adjudication."""
    for mod in GATE2_MODULES:
        code = open(os.path.join(CODE, mod)).read()
        for token in (
            "final_deletion_authorized",
            "deletion_authorized",
            "authorized_for_deletion",
            "FINAL_DELETION",
        ):
            for line in code.splitlines():
                if token in line and not line.strip().startswith("#") and '"""' not in line:
                    raise AssertionError(f"{mod}: {token!r} must not exist: {line.strip()[:80]}")


def test_r5_runner_output_uses_mechanical_candidate_language(tmp_path, capsys):
    inp = write_combined_cohort(tmp_path / "in", "t", ("A",))
    cfg = write_config(tmp_path / "c.csv", [config_row("t", "A")])
    r = runner.run(
        "t", "combined", inp, str(tmp_path / "out"), cfg, 7, scrublet_fn=fake_scrublet_success
    )
    assert "mechanical_doublet_candidates" in r
    assert "deletion_eligible" not in r
    out = capsys.readouterr().out
    assert "mechanical doublet candidates:" in out
    assert "NOTHING DELETED" in out
    assert "No cell has final deletion authorization." in out
    assert "deletion-eligible" not in out


def test_r5_calls_csv_has_no_authorization_column(tmp_path):
    inp = write_combined_cohort(tmp_path / "in", "t", ("A",))
    cfg = write_config(tmp_path / "c.csv", [config_row("t", "A")])
    r = runner.run(
        "t", "combined", inp, str(tmp_path / "out"), cfg, 7, scrublet_fn=fake_scrublet_success
    )
    cols = set(next(iter(csv.DictReader(open(r["calls_csv"])))).keys())
    forbidden = {c for c in cols if "authoriz" in c.lower() or "final_del" in c.lower()}
    assert not forbidden, f"no authorization column may exist, found {forbidden}"


def test_r5_no_function_produces_an_authorization_from_call_status():
    """Nothing in the module surface may map call_status -> authorization."""
    import inspect

    for name, obj in vars(st).items():
        if not inspect.isfunction(obj):
            continue
        assert "authoriz" not in name.lower(), f"{name} implies an authorization step"
        src = inspect.getsource(obj)
        body = src.split('"""')[-1] if src.count('"""') >= 2 else src
        assert "authoriz" not in body.lower(), f"{name} body claims an authorization"


def test_r5_unknown_candidate_still_cannot_be_filled_to_false(tmp_path):
    inp = write_combined_cohort(tmp_path / "in", "t", ("A",))
    cfg = write_config(tmp_path / "c.csv", [config_row("t", "A")])
    r = runner.run(
        "t", "combined", inp, str(tmp_path / "out"), cfg, 7, scrublet_fn=fake_scrublet_nan_threshold
    )
    rows_ = list(csv.DictReader(open(r["calls_csv"])))
    assert {x["candidate_doublet"] for x in rows_} == {""}
    assert all(x["call_status"] == st.CALL_THRESHOLD_FAILED for x in rows_)
    assert st.mechanical_candidate_mask(x["call_status"] for x in rows_) == [False] * len(rows_)
    assert r["mechanical_doublet_candidates"] == 0


# ==========================================================================
# Round 6 — auto_singlet -> scrublet_noncandidate, checkpoint schema /3
# ==========================================================================


def test_r6_old_singlet_vocabulary_is_gone():
    """Breaking rename, no alias. `auto_singlet` collided with the demuxlet
    sense of "singlet" (one DONOR, not one cell)."""
    assert not hasattr(st, "CALL_AUTO_SINGLET")
    assert st.CALL_SCRUBLET_NONCANDIDATE == "scrublet_noncandidate"
    assert "auto_singlet" not in st.CALL_STATUSES
    for forbidden in ("confirmed_singlet", "true_singlet", "final_singlet", "not_doublet"):
        assert forbidden not in st.CALL_STATUSES
        assert not hasattr(st, forbidden.upper())
    for mod in GATE2_MODULES:
        code = open(os.path.join(CODE, mod)).read()
        for line in code.splitlines():
            st_line = line.strip()
            if st_line.startswith("#") or '"""' in line:
                continue
            assert "auto_singlet" not in line, f"{mod}: stale value in {st_line[:80]}"
            assert "CALL_AUTO_SINGLET" not in line, f"{mod}: stale constant in {st_line[:80]}"


def test_r6_scored_capture_maps_exactly_two_statuses():
    rules = st.CAPTURE_INVARIANTS[st.CAPTURE_SCORED]
    assert rules["calls"] == {st.CALL_PREDICTED_DOUBLET, st.CALL_SCRUBLET_NONCANDIDATE}
    assert rules["candidate"] == {0, 1}
    assert st.CALL_TO_CANDIDATE[st.CALL_PREDICTED_DOUBLET] == 1
    assert st.CALL_TO_CANDIDATE[st.CALL_SCRUBLET_NONCANDIDATE] == 0
    for unknown in (
        st.CALL_TOO_SMALL,
        st.CALL_THRESHOLD_FAILED,
        st.CALL_METHOD_UNRELIABLE,
        st.CALL_RUNTIME_FAILED,
    ):
        assert st.CALL_TO_CANDIDATE[unknown] == -1


def test_r6_noncandidate_is_not_a_mechanical_candidate():
    """The rename must not change what counts as a candidate."""
    assert st.mechanical_candidate_mask([st.CALL_SCRUBLET_NONCANDIDATE]) == [False]
    assert st.mechanical_candidate_mask([st.CALL_PREDICTED_DOUBLET]) == [True]
    assert st.MECHANICAL_DOUBLET_CANDIDATE_STATUSES == {st.CALL_PREDICTED_DOUBLET}


def test_r6_real_scoring_emits_the_new_status():
    res = core.score_one_capture(
        make_matrix(),
        cell_ids(),
        expected_doublet_rate=0.05,
        seed=1,
        scrublet_fn=fake_scrublet_success,
    )
    assert res.capture_status == st.CAPTURE_SCORED
    assert set(res.call_status) <= {st.CALL_PREDICTED_DOUBLET, st.CALL_SCRUBLET_NONCANDIDATE}
    assert "auto_singlet" not in set(res.call_status)


@pytest.mark.parametrize("old", ["gate2-checkpoint/1", "gate2-checkpoint/2"])
def test_r6_superseded_schemas_are_rejected(tmp_path, old):
    assert ckpt.CHECKPOINT_SCHEMA_VERSION not in ckpt.REJECTED_SCHEMA_VERSIONS
    assert old in ckpt.REJECTED_SCHEMA_VERSIONS
    _save_raw(tmp_path, _fp(checkpoint_schema_version=old), scores=np.array([0.1, 0.2, 0.3]))
    with pytest.raises(ckpt.CheckpointMismatch, match="superseded schema"):
        ckpt.load_checkpoint(str(tmp_path), "k", _fp(), expected_cell_ids=["a", "b", "c"])


def test_r6_call_status_dtype_holds_the_longest_status(tmp_path):
    """'scrublet_noncandidate' is 21 chars, longer than any previous status;
    the derived width must not truncate it."""
    assert len("scrublet_noncandidate") == 21
    path = ckpt.save_checkpoint(
        str(tmp_path),
        "k",
        _fp(),
        scores=[0.1, 0.2, 0.3],
        call_status=[
            st.CALL_SCRUBLET_NONCANDIDATE,
            st.CALL_PREDICTED_DOUBLET,
            st.CALL_SCRUBLET_NONCANDIDATE,
        ],
        candidate=np.array([0, 1, 0], dtype=np.int8),
        cell_ids=["a", "b", "c"],
        capture_status=st.CAPTURE_SCORED,
        diagnostics=_conforming_diagnostics(n_predicted_doublets=1),
    )
    with np.load(path, allow_pickle=False) as z:
        restored = [str(x) for x in z["call_status"]]
    assert restored == [
        st.CALL_SCRUBLET_NONCANDIDATE,
        st.CALL_PREDICTED_DOUBLET,
        st.CALL_SCRUBLET_NONCANDIDATE,
    ]
    got = ckpt.load_checkpoint(str(tmp_path), "k", _fp(), expected_cell_ids=["a", "b", "c"])
    assert got["call_status"] == restored


def test_r6_noncandidate_confers_no_retention_authorization():
    """The constant must carry its disclaimers.

    An earlier version of this test searched for the phrase "true singlet" and
    matched the disclaimer's own negation. Absence of a phrase is not the
    property worth checking; presence of the limits is.
    """
    src = open(os.path.join(CODE, "gate2_status.py")).read().lower()
    i = src.index("call_scrublet_noncandidate =")
    block = src[max(0, i - 1400) : i]
    for disclaimer in (
        "does not assert the cell is a true singlet",
        "does not confer any retention",
        "does not exclude a homotypic doublet",
    ):
        assert disclaimer in block, f"missing disclaimer: {disclaimer!r}"
    # and the enum value itself must not use singlet vocabulary
    assert "singlet" not in st.CALL_SCRUBLET_NONCANDIDATE


# ==========================================================================
# Round 7 — score-distribution diagnostics (evidence, NOT a reliability rule)
# ==========================================================================

# ---- AUC: analytic values, fed explicit arrays (not via a mock) -----------


def test_r7_auc_perfect_separation_is_one():
    """simulated is the POSITIVE class, so simulated strictly above observed
    must give exactly 1."""
    assert dg.observed_simulated_auc(np.array([0.0, 0.1, 0.2]), np.array([0.8, 0.9, 1.0])) == 1.0


def test_r7_auc_fully_inverted_is_zero():
    assert dg.observed_simulated_auc(np.array([0.8, 0.9, 1.0]), np.array([0.0, 0.1, 0.2])) == 0.0


def test_r7_auc_identical_distributions_is_half():
    x = np.array([0.1, 0.2, 0.3, 0.4])
    assert dg.observed_simulated_auc(x, x.copy()) == 0.5


def test_r7_auc_ties_get_exactly_half_weight():
    # every simulated ties with every observed -> P(sim>obs)=0, P(tie)=1
    assert dg.observed_simulated_auc(np.full(4, 0.5), np.full(6, 0.5)) == 0.5
    # half the simulated strictly above, half tied
    obs = np.array([0.0, 0.0])
    sim = np.array([0.0, 1.0])
    # P(sim>obs) = 2/4 = 0.5 ; P(tie) = 2/4 = 0.5 -> 0.5 + 0.25 = 0.75
    assert dg.observed_simulated_auc(obs, sim) == 0.75


def test_r7_auc_matches_brute_force_on_random_data():
    rng = np.random.default_rng(0)
    obs = rng.integers(0, 5, 40).astype(float)  # many ties on purpose
    sim = rng.integers(0, 5, 55).astype(float)
    gt = float(np.mean([(s > o) + 0.5 * (s == o) for s in sim for o in obs]))
    assert abs(dg.observed_simulated_auc(obs, sim) - gt) < 1e-12


# ---- threshold-based fields ----------------------------------------------


def test_r7_threshold_metrics_use_ge_at_the_boundary():
    obs = np.array([0.1, 0.5, 0.9])
    sim = np.array([0.5, 0.5, 0.1])
    d = dg.score_distribution_summary(obs, sim, threshold=0.5)
    assert d["observed_tail_fraction_at_threshold"] == 2 / 3  # 0.5 and 0.9
    assert d["simulated_recall_at_threshold"] == 2 / 3  # both 0.5
    assert d["simulated_tail_fraction_at_threshold"] == 2 / 3


@pytest.mark.parametrize("thr", [None, float("nan"), float("inf")])
def test_r7_threshold_metrics_are_null_not_zero(thr):
    d = dg.score_distribution_summary(np.array([0.1, 0.2]), np.array([0.3, 0.4]), threshold=thr)
    for f in (
        "simulated_recall_at_threshold",
        "observed_tail_fraction_at_threshold",
        "simulated_tail_fraction_at_threshold",
    ):
        assert d[f] is None, f"{f} must be null, not 0 -- 0 would read as a measurement"
    # distribution evidence survives a failed threshold
    assert d["observed_simulated_auc"] is not None
    assert d["observed_score_percentiles_p0_to_p100"] is not None
    assert d["shared_histogram_bin_edges"] is not None


# ---- degenerate and missing inputs ---------------------------------------


def test_r7_simulated_absent_or_empty_yields_nulls_and_a_problem():
    for sim in (None, np.array([])):
        d = dg.score_distribution_summary(np.array([0.1, 0.2]), sim, threshold=0.5)
        assert d["n_simulated"] is None
        assert d["observed_simulated_auc"] is None
        assert d["simulated_recall_at_threshold"] is None
        assert any("simulated" in p for p in d["distribution_problems"])
        assert d["n_observed"] == 2, "observed evidence is still reported"


@pytest.mark.parametrize("bad", [np.array([0.1, np.nan]), np.array([0.1, np.inf])])
def test_r7_non_finite_simulated_is_refused_with_a_problem(bad):
    d = dg.score_distribution_summary(np.array([0.1, 0.2]), bad, threshold=0.5)
    assert d["n_simulated"] is None
    assert any("non-finite" in p for p in d["distribution_problems"])


def test_r7_constant_distributions_do_not_crash_or_make_zero_width_bins():
    d = dg.score_distribution_summary(np.full(5, 0.3), np.full(7, 0.3), threshold=0.3)
    assert d["histogram_range_degenerate"] is True
    edges = d["shared_histogram_bin_edges"]
    assert len(edges) == dg.HISTOGRAM_N_BINS + 1
    assert all(edges[i] < edges[i + 1] for i in range(len(edges) - 1)), "no zero-width bin"
    assert sum(d["observed_histogram_counts"]) == 5
    assert sum(d["simulated_histogram_counts"]) == 7
    assert d["observed_simulated_auc"] == 0.5


def test_r7_scores_outside_unit_interval_are_not_clipped():
    """Clipping to [0,1] would hide an out-of-range score, which is exactly
    what these diagnostics exist to surface."""
    d = dg.score_distribution_summary(np.array([-0.2, 0.5]), np.array([0.5, 1.7]))
    assert d["observed_score_min"] == -0.2
    assert d["simulated_score_max"] == 1.7
    assert d["shared_histogram_bin_edges"][0] == -0.2
    assert d["shared_histogram_bin_edges"][-1] == 1.7


# ---- percentiles and histograms ------------------------------------------


def test_r7_percentiles_are_101_points_monotone_with_correct_endpoints():
    rng = np.random.default_rng(1)
    x = rng.random(500)
    d = dg.score_distribution_summary(x, x.copy())
    p = d["observed_score_percentiles_p0_to_p100"]
    assert len(p) == 101
    assert all(p[i] <= p[i + 1] + 1e-12 for i in range(100)), "must be non-decreasing"
    assert abs(p[0] - float(x.min())) < 1e-12
    assert abs(p[100] - float(x.max())) < 1e-12


def test_r7_histograms_share_edges_and_counts_sum_to_n():
    rng = np.random.default_rng(2)
    obs, sim = rng.random(300), rng.random(450) * 0.8
    d = dg.score_distribution_summary(obs, sim)
    assert len(d["shared_histogram_bin_edges"]) == dg.HISTOGRAM_N_BINS + 1
    assert len(d["observed_histogram_counts"]) == dg.HISTOGRAM_N_BINS
    assert len(d["simulated_histogram_counts"]) == dg.HISTOGRAM_N_BINS
    assert sum(d["observed_histogram_counts"]) == 300
    assert sum(d["simulated_histogram_counts"]) == 450


def test_r7_summary_is_plain_python_and_json_round_trips():
    import json

    rng = np.random.default_rng(3)
    d = dg.score_distribution_summary(rng.random(50), rng.random(80), threshold=0.5)
    for k, v in d.items():
        assert not isinstance(v, np.generic), f"{k} is a numpy scalar"
        assert not isinstance(v, np.ndarray), f"{k} is a numpy array"
        if isinstance(v, list):
            assert not any(isinstance(e, (np.generic, np.ndarray)) for e in v), k
    assert json.loads(json.dumps(d)) == d


# ---- integration through the real scoring path ---------------------------


def test_r7_scoring_populates_every_distribution_field():
    res = core.score_one_capture(
        make_matrix(),
        cell_ids(),
        expected_doublet_rate=0.05,
        seed=1,
        scrublet_fn=fake_scrublet_success,
    )
    for f in dg.DISTRIBUTION_FIELDS:
        assert f in res.diagnostics, f"{f} missing from diagnostics"
    assert res.diagnostics["diagnostics_schema_version"] == dg.DIAGNOSTICS_SCHEMA_VERSION
    assert res.diagnostics["n_observed"] == N_CELLS
    assert res.diagnostics["observed_simulated_auc"] is not None


def test_r7_threshold_failed_keeps_distribution_evidence_but_no_threshold_metrics():
    res = core.score_one_capture(
        make_matrix(),
        cell_ids(),
        expected_doublet_rate=0.05,
        seed=1,
        scrublet_fn=fake_scrublet_nan_threshold,
    )
    d = res.diagnostics
    assert res.capture_status == st.CAPTURE_THRESHOLD_FAILED
    assert d["observed_score_percentiles_p0_to_p100"] is not None
    assert d["simulated_recall_at_threshold"] is None
    assert d["observed_tail_fraction_at_threshold"] is None


def test_r7_full_simulated_dump_is_off_by_default_and_opt_in(tmp_path):
    out = tmp_path / "sim.npy"
    core.score_one_capture(
        make_matrix(),
        cell_ids(),
        expected_doublet_rate=0.05,
        seed=1,
        scrublet_fn=fake_scrublet_success,
    )
    assert not out.exists(), "production runs must not write simulated scores"
    core.score_one_capture(
        make_matrix(),
        cell_ids(),
        expected_doublet_rate=0.05,
        seed=1,
        scrublet_fn=fake_scrublet_success,
        dump_simulated_scores_to=str(out),
    )
    assert out.exists() and np.load(str(out)).size > 0


# ---- checkpoint schema /4 and fresh-vs-resume on the new fields -----------


def test_r7_schema_bumped_to_4_and_3_rejected(tmp_path):
    assert ckpt.CHECKPOINT_SCHEMA_VERSION == "gate2-checkpoint/4"
    assert "gate2-checkpoint/3" in ckpt.REJECTED_SCHEMA_VERSIONS
    _save_raw(
        tmp_path,
        _fp(checkpoint_schema_version="gate2-checkpoint/3"),
        scores=np.array([0.1, 0.2, 0.3]),
    )
    with pytest.raises(ckpt.CheckpointMismatch, match="superseded schema"):
        ckpt.load_checkpoint(str(tmp_path), "k", _fp(), expected_cell_ids=["a", "b", "c"])


def test_r7_fresh_vs_resume_identical_on_every_new_field(tmp_path):
    inp = write_combined_cohort(tmp_path / "in", "t", ("A",))
    cfg = write_config(tmp_path / "c.csv", [config_row("t", "A")])
    out = tmp_path / "out"
    r1 = runner.run("t", "combined", inp, str(out), cfg, 7, scrublet_fn=fake_scrublet_success)
    r2 = runner.run("t", "combined", inp, str(out), cfg, 7, scrublet_fn=fake_scrublet_success)
    assert r2["restored"] == 1
    d1 = list(csv.DictReader(open(r1["diagnostics_csv"])))[0]
    d2 = list(csv.DictReader(open(r2["diagnostics_csv"])))[0]
    for f in dg.DISTRIBUTION_FIELDS:
        assert f in d1, f"{f} missing from the diagnostics CSV"
        assert d1[f] == d2[f], f"{f} differs between fresh and resume"


# ==========================================================================
# Round 8 — close the diagnostics contract on the restore side
# ==========================================================================


def test_r8_diagnostics_module_is_in_the_code_fingerprint():
    """Editing a metric definition changes the persisted diagnostics without
    touching any other module; the fingerprint must notice."""
    names = [os.path.basename(f) for f in runner.CODE_FILES]
    assert "gate2_diagnostics.py" in names, names
    for expected in (
        "gate2_run.py",
        "gate2_score_core.py",
        "gate2_checkpoint.py",
        "gate2_runtime_config.py",
        "gate2_status.py",
    ):
        assert expected in names


def test_r8_editing_the_diagnostics_module_changes_the_fingerprint(tmp_path):
    a = tmp_path / "gate2_diagnostics.py"
    a.write_text("DIAGNOSTICS_SCHEMA_VERSION = 'x/1'\n")
    before = ckpt.code_fingerprint([str(a)])
    a.write_text("DIAGNOSTICS_SCHEMA_VERSION = 'x/1'\n# a metric definition changed\n")
    assert ckpt.code_fingerprint([str(a)]) != before


def test_r8_missing_diagnostics_schema_version_is_refused(tmp_path):
    """Written before the diagnostics contract existed -- refused, not assumed
    compatible."""
    # deliberately NOT _conforming_diagnostics(): the version is absent
    _save_raw(tmp_path, _fp(), scores=np.array([0.1, 0.2, 0.3]), diagnostics={"threshold": 0.9})
    with pytest.raises(ckpt.CheckpointMismatch, match="no diagnostics_schema_version"):
        ckpt.load_checkpoint(str(tmp_path), "k", _fp(), expected_cell_ids=["a", "b", "c"])


@pytest.mark.parametrize(
    "stored", ["gate2-diagnostics/0", "gate2-diagnostics/2", "something-else/1", ""]
)
def test_r8_old_or_unknown_diagnostics_schema_is_refused(tmp_path, stored):
    _save_raw(
        tmp_path,
        _fp(),
        scores=np.array([0.1, 0.2, 0.3]),
        diagnostics={"threshold": 0.9, "diagnostics_schema_version": stored},
    )
    with pytest.raises(ckpt.CheckpointMismatch, match="diagnostics_schema_version"):
        ckpt.load_checkpoint(str(tmp_path), "k", _fp(), expected_cell_ids=["a", "b", "c"])


def test_r8_current_diagnostics_schema_restores(tmp_path):
    _save_raw(
        tmp_path,
        _fp(),
        scores=np.array([0.1, 0.2, 0.3]),
        diagnostics={"threshold": 0.9, "diagnostics_schema_version": dg.DIAGNOSTICS_SCHEMA_VERSION},
    )
    got = ckpt.load_checkpoint(str(tmp_path), "k", _fp(), expected_cell_ids=["a", "b", "c"])
    assert got is not None
    assert got["diagnostics"]["diagnostics_schema_version"] == dg.DIAGNOSTICS_SCHEMA_VERSION


def test_r8_real_run_writes_a_validatable_diagnostics_version(tmp_path):
    """End to end: what the runner writes must satisfy the restore check."""
    inp = write_combined_cohort(tmp_path / "in", "t", ("A",))
    cfg = write_config(tmp_path / "c.csv", [config_row("t", "A")])
    out = tmp_path / "out"
    runner.run("t", "combined", inp, str(out), cfg, 7, scrublet_fn=fake_scrublet_success)
    r2 = runner.run("t", "combined", inp, str(out), cfg, 7, scrublet_fn=fake_scrublet_success)
    assert r2["restored"] == 1, "a freshly written checkpoint must pass its own schema check"

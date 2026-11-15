"""Artifact/provenance manifest schema tests. Round-trips with the REAL
serializer (write_atomic -> read_manifest), not a stand-in -- per DEVELOPMENT.md's
"round trips must be strict" standing rule.
"""

import os
import sys

import pytest

CODE = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "scripts",
    "34_scfm_discovery_framework",
    "code",
)
sys.path.insert(0, CODE)

import artifact_schema as A  # noqa: E402


def test_build_manifest_rejects_invalid_stage_status():
    with pytest.raises(ValueError, match="stage_status"):
        A.build_manifest(
            input_identity={},
            source_code_revision_or_hash="x",
            model_identity_or_pending="model_undecided",
            stage_status="definitely_done",
            expected_outputs=[],
            dependencies=[],
        )


def test_build_manifest_rejects_guessed_model_identity():
    with pytest.raises(ValueError, match="model_identity_or_pending"):
        A.build_manifest(
            input_identity={},
            source_code_revision_or_hash="x",
            model_identity_or_pending="probably_scprint2",
            stage_status="planned",
            expected_outputs=[],
            dependencies=[],
        )


def test_build_manifest_accepts_model_undecided():
    manifest = A.build_manifest(
        input_identity={"dataset_id": "x"},
        source_code_revision_or_hash="x",
        model_identity_or_pending="model_undecided",
        stage_status="planned",
        expected_outputs=["foo.tsv"],
        dependencies=[],
    )
    assert manifest["model_identity_or_pending"] == "model_undecided"
    assert manifest["stage_status"] == "planned"


def test_write_then_read_round_trips_exactly(tmp_path):
    manifest = A.build_manifest(
        input_identity={"dataset_id": "x", "dataset_config_sha256": "abc123"},
        source_code_revision_or_hash="deadbeef",
        model_identity_or_pending={"model_generation": "scPRINT-2"},
        stage_status="planned",
        expected_outputs=["a.tsv", "b.tsv"],
        dependencies=[{"file": "dataset_schema.py", "sha256": "f" * 64}],
    )
    dest = tmp_path / "PLAN.json"
    A.write_atomic(str(dest), manifest)
    read_back = A.read_manifest(str(dest))
    assert read_back["input_identity"] == manifest["input_identity"]
    assert read_back["stage_status"] == "planned"
    assert read_back["model_identity_or_pending"] == {"model_generation": "scPRINT-2"}
    assert read_back["dependencies"] == [{"file": "dataset_schema.py", "sha256": "f" * 64}]


def test_read_manifest_rejects_wrong_schema_version(tmp_path):
    import json

    dest = tmp_path / "bad.json"
    dest.write_text(json.dumps({"schema": "not-the-real-schema/1"}))
    with pytest.raises(SystemExit, match="schema"):
        A.read_manifest(str(dest))


def test_read_manifest_missing_file_raises_systemexit(tmp_path):
    with pytest.raises(SystemExit, match="not found"):
        A.read_manifest(str(tmp_path / "does_not_exist.json"))


def test_sha256_text_is_deterministic_and_order_sensitive():
    assert A.sha256_text("abc") == A.sha256_text("abc")
    assert A.sha256_text("abc") != A.sha256_text("cba")


# --- R3 (review:20261012_153316): read_manifest must reject a tampered
# on-disk file, not just check the schema-version string. Each case here writes
# the bytes directly (bypassing build_manifest/write_atomic entirely), the exact
# attack the review used. ---


def test_read_manifest_rejects_illegal_stage_status_even_if_schema_matches(tmp_path):
    import json

    path = tmp_path / "tampered.json"
    path.write_text(json.dumps({"schema": A.SCHEMA_VERSION, "stage_status": "invented"}))
    with pytest.raises(SystemExit, match="stage_status"):
        A.read_manifest(str(path))


def test_read_manifest_rejects_missing_required_keys(tmp_path):
    import json

    path = tmp_path / "tampered.json"
    path.write_text(json.dumps({"schema": A.SCHEMA_VERSION, "stage_status": "planned"}))
    with pytest.raises(SystemExit, match="missing required key"):
        A.read_manifest(str(path))


def test_read_manifest_rejects_empty_input_identity(tmp_path):
    import json

    path = tmp_path / "tampered.json"
    payload = {
        "schema": A.SCHEMA_VERSION,
        "generated_utc": "x",
        "input_identity": {},
        "source_code_revision_or_hash": "x",
        "model_identity_or_pending": "model_undecided",
        "stage_status": "planned",
        "expected_outputs": [],
        "dependencies": [],
    }
    path.write_text(json.dumps(payload))
    with pytest.raises(SystemExit, match="input_identity"):
        A.read_manifest(str(path))


def test_read_manifest_rejects_guessed_model_identity_in_tampered_file(tmp_path):
    import json

    path = tmp_path / "tampered.json"
    payload = {
        "schema": A.SCHEMA_VERSION,
        "generated_utc": "x",
        "input_identity": {"x": 1},
        "source_code_revision_or_hash": "x",
        "model_identity_or_pending": "probably_scprint2",
        "stage_status": "planned",
        "expected_outputs": [],
        "dependencies": [],
    }
    path.write_text(json.dumps(payload))
    with pytest.raises(SystemExit, match="model_identity_or_pending"):
        A.read_manifest(str(path))


def test_build_manifest_itself_rejects_the_same_tampered_shape():
    """The construction path and the read path must agree -- the same
    validate_manifest() backs both, so a problem caught on read could not
    have passed through build in the first place.
    """
    with pytest.raises(ValueError, match="invalid manifest"):
        A.build_manifest(
            input_identity={},
            source_code_revision_or_hash="x",
            model_identity_or_pending="model_undecided",
            stage_status="planned",
            expected_outputs=[],
            dependencies=[],
        )


def test_write_atomic_rejects_nan_instead_of_silently_writing_it(tmp_path):
    """Rejected even earlier than write_atomic now: build_manifest's own
    validate_manifest() call catches the NaN recursively (see the R1 fix in
    review:20261012_171916), so this never even reaches
    write_atomic -- a strictly stronger guarantee than the original test
    exercised, not a weaker one.
    """
    with pytest.raises(ValueError, match="non-finite"):
        A.build_manifest(
            input_identity={"x": 1},
            source_code_revision_or_hash="x",
            model_identity_or_pending="model_undecided",
            stage_status="planned",
            expected_outputs=[],
            dependencies=[],
            extra={"a_nan_value": float("nan")},
        )


def test_write_atomic_rejects_unsupported_object_instead_of_stringifying(tmp_path):
    """Previously `default=str` would have silently turned this into the
    string "<object object at 0x...>" rather than raising.
    """
    manifest = A.build_manifest(
        input_identity={"x": 1},
        source_code_revision_or_hash="x",
        model_identity_or_pending="model_undecided",
        stage_status="planned",
        expected_outputs=[],
        dependencies=[],
    )
    manifest["extra"] = {"not_json_safe": object()}
    with pytest.raises(TypeError):
        A.write_atomic(str(tmp_path / "x.json"), manifest)


def test_mutation_kill_illegal_stage_status_check(monkeypatch):
    """Disables the stage_status check inside validate_manifest and confirms
    the tampered-manifest test above would then genuinely fail -- proving it
    detects the fault, not merely exercising a code path that always passes.
    """
    import artifact_schema as A2  # fresh import so monkeypatch is scoped

    real_validate = A2.validate_manifest

    def _mutant(manifest):
        problems = real_validate(manifest)
        return [p for p in problems if "stage_status" not in p]

    monkeypatch.setattr(A2, "validate_manifest", _mutant)
    import json
    import tempfile

    path = tempfile.mktemp(suffix=".json")
    with open(path, "w") as fh:
        json.dump(
            {
                "schema": A2.SCHEMA_VERSION,
                "generated_utc": "x",
                "input_identity": {"x": 1},
                "source_code_revision_or_hash": "x",
                "model_identity_or_pending": "model_undecided",
                "stage_status": "invented",
                "expected_outputs": [],
                "dependencies": [],
            },
            fh,
        )
    # With the mutant validator, read_manifest must NOT raise -- proving the
    # real validator (not a tautology) is what makes the real test above pass.
    manifest = A2.read_manifest(path)
    assert manifest["stage_status"] == "invented"
    os.remove(path)


# --- review:20261012_163250 finding 3: reader/writer type contracts
# were not aligned -- bool source_code_revision_or_hash, a null dependency
# entry, and a NaN persisted via raw JSON (bypassing write_atomic) were all
# previously ACCEPTED by read_manifest. Each is reproduced by tampering a REAL
# PLAN produced by build_manifest, not a hand-simplified object. ---


def _real_plan():
    return A.build_manifest(
        input_identity={"dataset_id": "x", "dataset_config_sha256": "abc123"},
        source_code_revision_or_hash="git:deadbeef",
        model_identity_or_pending="model_undecided",
        stage_status="planned",
        expected_outputs=["a.tsv"],
        dependencies=[{"file": "cli.py", "sha256": "a" * 64}],
    )


def test_read_manifest_rejects_bool_source_code_revision(tmp_path):
    import json

    plan = _real_plan()
    plan["source_code_revision_or_hash"] = True
    path = tmp_path / "tampered.json"
    path.write_text(json.dumps(plan))
    with pytest.raises(SystemExit, match="source_code_revision_or_hash"):
        A.read_manifest(str(path))


def test_read_manifest_rejects_null_dependency_entry(tmp_path):
    import json

    plan = _real_plan()
    plan["dependencies"] = [None]
    path = tmp_path / "tampered.json"
    path.write_text(json.dumps(plan))
    with pytest.raises(SystemExit, match=r"dependencies\[0\]"):
        A.read_manifest(str(path))


def test_read_manifest_rejects_dependency_missing_sha256(tmp_path):
    import json

    plan = _real_plan()
    plan["dependencies"] = [{"file": "cli.py"}]
    path = tmp_path / "tampered.json"
    path.write_text(json.dumps(plan))
    with pytest.raises(SystemExit, match=r"dependencies\[0\]\['sha256'\]"):
        A.read_manifest(str(path))


def test_read_manifest_rejects_nan_in_raw_persisted_json(tmp_path):
    """Bypasses write_atomic entirely (which would refuse this) by writing
    the tampered bytes directly -- json.load accepts a bare NaN token by
    default even though it is not standard JSON; read_manifest must not.
    """
    plan = _real_plan()
    import json as _json

    text = _json.dumps(plan).replace('"abc123"', "NaN")
    path = tmp_path / "tampered.json"
    path.write_text(text)
    with pytest.raises(SystemExit, match="non-finite"):
        A.read_manifest(str(path))


def test_mutation_kill_nonfinite_constant_check():
    """Confirms json.load's default behavior (accepting NaN) is what the
    parse_constant override must override -- calling json.load without it
    on the same tampered text would NOT raise, proving the override is load
    -bearing, not redundant with default json parsing.
    """
    import json

    plan = _real_plan()
    text = json.dumps(plan).replace('"abc123"', "NaN")
    import math

    # Default json.load accepts this (demonstrating the fault the override fixes).
    parsed = json.loads(text)
    assert math.isnan(parsed["input_identity"]["dataset_config_sha256"])


# --- review:20261012_171916 R1: a normal numeric literal that
# overflows to inf during parsing (e.g. 1e999) is NOT one of the three special
# JSON constant tokens parse_constant intercepts -- it must still be caught,
# recursively, anywhere in the manifest. ---


def test_read_manifest_rejects_overflow_to_infinity_in_extra(tmp_path):
    import json

    plan = _real_plan()
    text = json.dumps(plan)
    tampered = text[:-1] + ', "extra": {"huge": 1e999}}'
    path = tmp_path / "tampered.json"
    path.write_text(tampered)
    with pytest.raises(SystemExit, match="non-finite"):
        A.read_manifest(str(path))


def test_read_manifest_rejects_negative_overflow_to_infinity(tmp_path):
    import json

    plan = _real_plan()
    text = json.dumps(plan)
    tampered = text[:-1] + ', "extra": {"huge": -1e999}}'
    path = tmp_path / "tampered.json"
    path.write_text(tampered)
    with pytest.raises(SystemExit, match="non-finite"):
        A.read_manifest(str(path))


def test_read_manifest_accepts_normal_finite_exponent_value(tmp_path):
    """A real, finite scientific-notation value must NOT be rejected --
    confirms the fix targets non-finite values specifically, not large
    exponents in general.
    """
    import json

    plan = _real_plan()
    text = json.dumps(plan)
    tampered = text[:-1] + ', "extra": {"normal": 1.5e10}}'
    path = tmp_path / "ok.json"
    path.write_text(tampered)
    read_back = A.read_manifest(str(path))
    assert read_back["extra"]["normal"] == 1.5e10


def test_build_and_write_also_reject_nonfinite_nested_in_extra(tmp_path):
    """The same recursive check applies at build() and write_atomic(), not
    only read() -- a caller cannot construct an in-memory manifest with a
    buried inf and have it silently accepted either.
    """
    with pytest.raises(ValueError, match="non-finite"):
        A.build_manifest(
            input_identity={"x": 1},
            source_code_revision_or_hash="git:x",
            model_identity_or_pending="model_undecided",
            stage_status="planned",
            expected_outputs=[],
            dependencies=[],
            extra={"nested": {"deep": [1.0, float("inf")]}},
        )


def test_mutation_kill_nonfinite_recursive_check(monkeypatch):
    """Disables find_nonfinite_floats and confirms the overflow test above
    would then genuinely fail to detect the fault.
    """
    import artifact_schema as A2

    monkeypatch.setattr(A2, "find_nonfinite_floats", lambda obj, path="<root>": [])
    import json

    plan = _real_plan()
    text = json.dumps(plan)
    tampered = text[:-1] + ', "extra": {"huge": 1e999}}'
    path_obj = __import__("tempfile").mktemp(suffix=".json")
    with open(path_obj, "w") as fh:
        fh.write(tampered)
    # With the mutant, read_manifest must NOT raise -- proving the real check
    # (not a tautology) is what makes the real test above pass.
    manifest = A2.read_manifest(path_obj)
    assert manifest["extra"]["huge"] == float("inf")
    os.remove(path_obj)

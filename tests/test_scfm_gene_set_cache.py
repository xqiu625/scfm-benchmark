"""Tests for the local, versioned gene-set cache. Never touches a URL or an
online library name -- every case here is a local .gmt file on disk.
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

import gene_set_cache as G  # noqa: E402


def _write_gmt(path, sets: dict):
    lines = []
    for name, genes in sets.items():
        lines.append("\t".join([name, f"{name}_description", *genes]))
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")


def test_parse_gmt_reads_sets_correctly(tmp_path):
    gmt = tmp_path / "test.gmt"
    _write_gmt(gmt, {"SET_A": ["G1", "G2", "G3"], "SET_B": ["G4", "G5"]})
    sets = G.parse_gmt(str(gmt))
    assert sets == {"SET_A": ["G1", "G2", "G3"], "SET_B": ["G4", "G5"]}


def test_parse_gmt_missing_file_raises():
    with pytest.raises(SystemExit, match="not found"):
        G.parse_gmt("/does/not/exist.gmt")


def test_parse_gmt_rejects_malformed_line(tmp_path):
    gmt = tmp_path / "bad.gmt"
    gmt.write_text("ONLY_NAME_AND_DESC\tdesc\n")
    with pytest.raises(SystemExit, match="malformed"):
        G.parse_gmt(str(gmt))


def test_parse_gmt_rejects_duplicate_set_name(tmp_path):
    gmt = tmp_path / "dupe.gmt"
    gmt.write_text("SET_A\tdesc\tG1\tG2\nSET_A\tdesc2\tG3\tG4\n")
    with pytest.raises(SystemExit, match="duplicate gene set name"):
        G.parse_gmt(str(gmt))


def test_parse_gmt_rejects_empty_file(tmp_path):
    gmt = tmp_path / "empty.gmt"
    gmt.write_text("")
    with pytest.raises(SystemExit, match="no gene sets parsed"):
        G.parse_gmt(str(gmt))


def test_import_gene_sets_records_real_metadata_not_fabricated(tmp_path):
    gmt = tmp_path / "test.gmt"
    _write_gmt(gmt, {"SET_A": ["G1", "G2"]})
    out_dir = tmp_path / "cache"
    manifest = G.import_gene_sets(
        source_gmt_path=str(gmt),
        collection_name="TEST_COLLECTION",
        species="human",
        gene_id_namespace="gene_symbol",
        source_version="local-test-v1",
        out_dir=str(out_dir),
    )
    assert manifest["source_type"] == "local_file"
    assert manifest["collection_name"] == "TEST_COLLECTION"
    assert manifest["n_gene_sets"] == 1
    assert manifest["source_content_sha256"] == G.sha256_file(str(gmt))
    # "imported_utc" is the local import time, never a fabricated live
    # retrieval timestamp -- confirm the field name itself says so.
    assert "imported_utc" in manifest
    assert "retrieved_utc" not in manifest


def test_import_gene_sets_requires_all_identity_fields_nonempty(tmp_path):
    gmt = tmp_path / "test.gmt"
    _write_gmt(gmt, {"SET_A": ["G1"]})
    with pytest.raises(SystemExit, match="must all be non-empty"):
        G.import_gene_sets(
            source_gmt_path=str(gmt),
            collection_name="",
            species="human",
            gene_id_namespace="gene_symbol",
            source_version="v1",
            out_dir=str(tmp_path / "cache"),
        )


def test_import_gene_sets_refuses_nonempty_out_dir(tmp_path):
    gmt = tmp_path / "test.gmt"
    _write_gmt(gmt, {"SET_A": ["G1"]})
    out_dir = tmp_path / "cache"
    out_dir.mkdir()
    (out_dir / "junk.txt").write_text("pre-existing")
    with pytest.raises(SystemExit, match="must be new or empty"):
        G.import_gene_sets(
            source_gmt_path=str(gmt),
            collection_name="TEST",
            species="human",
            gene_id_namespace="gene_symbol",
            source_version="v1",
            out_dir=str(out_dir),
        )
    assert (out_dir / "junk.txt").read_text() == "pre-existing"


def test_load_cached_gene_sets_round_trips_exactly(tmp_path):
    gmt = tmp_path / "test.gmt"
    _write_gmt(gmt, {"SET_A": ["G1", "G2"], "SET_B": ["G3"]})
    out_dir = tmp_path / "cache"
    G.import_gene_sets(
        source_gmt_path=str(gmt),
        collection_name="TEST",
        species="mouse",
        gene_id_namespace="gene_symbol",
        source_version="v1",
        out_dir=str(out_dir),
    )
    gene_sets, manifest = G.load_cached_gene_sets(str(out_dir))
    assert gene_sets == {"SET_A": ["G1", "G2"], "SET_B": ["G3"]}
    assert manifest["species"] == "mouse"


def test_load_cached_gene_sets_missing_manifest_raises(tmp_path):
    with pytest.raises(SystemExit, match="no CACHE_MANIFEST.json"):
        G.load_cached_gene_sets(str(tmp_path))


def test_load_cached_gene_sets_rejects_tampered_content(tmp_path):
    """Corrupts the cached gene_sets.json AFTER import -- the manifest's own
    recorded hash must catch this, with no fallback to a live fetch. The
    tamper is syntactically valid JSON with different CONTENT (an added
    gene), not corrupted bytes -- a content-only change is the harder, more
    realistic case to catch (a bytes-level corruption would also break
    json.load on its own, which would mask whether the hash check itself
    is doing any work).
    """
    import json

    gmt = tmp_path / "test.gmt"
    _write_gmt(gmt, {"SET_A": ["G1", "G2"]})
    out_dir = tmp_path / "cache"
    manifest = G.import_gene_sets(
        source_gmt_path=str(gmt),
        collection_name="TEST",
        species="human",
        gene_id_namespace="gene_symbol",
        source_version="v1",
        out_dir=str(out_dir),
    )
    with open(manifest["cached_gene_sets_path"], "w") as fh:
        json.dump({"SET_A": ["G1", "G2", "SNEAKILY_ADDED_GENE"]}, fh)
    with pytest.raises(SystemExit, match="hash mismatch"):
        G.load_cached_gene_sets(str(out_dir))


def test_load_cached_gene_sets_rejects_manifest_missing_a_required_field(tmp_path):
    """review:20261012_220407 F3/F4: a hash match proves content
    integrity, not metadata validity -- a manifest missing a required field
    (here: 'species') must still be rejected on load even though the
    gene-sets file's own content hash is untouched and correct.
    """
    import json

    gmt = tmp_path / "test.gmt"
    _write_gmt(gmt, {"SET_A": ["G1", "G2"]})
    out_dir = tmp_path / "cache"
    G.import_gene_sets(
        source_gmt_path=str(gmt),
        collection_name="TEST",
        species="human",
        gene_id_namespace="gene_symbol",
        source_version="v1",
        out_dir=str(out_dir),
    )
    manifest_path = out_dir / "CACHE_MANIFEST.json"
    manifest = json.loads(manifest_path.read_text())
    del manifest["species"]
    manifest_path.write_text(json.dumps(manifest))
    with pytest.raises(SystemExit, match="species"):
        G.load_cached_gene_sets(str(out_dir))


def test_mutation_kill_cache_manifest_schema_validation(monkeypatch, tmp_path):
    """Disables validate_cache_manifest and confirms the test above's
    field-missing manifest would then be silently ACCEPTED (its gene-sets
    file content hash is still correct, so the hash check alone lets it
    through) -- proving the real schema validation, not the hash check, is
    what the test above actually relies on.
    """
    import json

    gmt = tmp_path / "test.gmt"
    _write_gmt(gmt, {"SET_A": ["G1", "G2"]})
    out_dir = tmp_path / "cache"
    G.import_gene_sets(
        source_gmt_path=str(gmt),
        collection_name="TEST",
        species="human",
        gene_id_namespace="gene_symbol",
        source_version="v1",
        out_dir=str(out_dir),
    )
    manifest_path = out_dir / "CACHE_MANIFEST.json"
    manifest = json.loads(manifest_path.read_text())
    del manifest["species"]
    manifest_path.write_text(json.dumps(manifest))

    monkeypatch.setattr(G, "validate_cache_manifest", lambda m: [])
    gene_sets, loaded_manifest = G.load_cached_gene_sets(str(out_dir))  # does NOT raise
    assert "species" not in loaded_manifest
    assert gene_sets == {"SET_A": ["G1", "G2"]}


def test_mutation_kill_hash_verification(monkeypatch, tmp_path):
    """Disables the hash check inside load_cached_gene_sets and confirms the
    tamper-rejection test above would then genuinely fail to detect it --
    the same content-level tamper now gets through silently.
    """
    import json

    gmt = tmp_path / "test.gmt"
    _write_gmt(gmt, {"SET_A": ["G1", "G2"]})
    out_dir = tmp_path / "cache"
    manifest = G.import_gene_sets(
        source_gmt_path=str(gmt),
        collection_name="TEST",
        species="human",
        gene_id_namespace="gene_symbol",
        source_version="v1",
        out_dir=str(out_dir),
    )
    with open(manifest["cached_gene_sets_path"], "w") as fh:
        json.dump({"SET_A": ["G1", "G2", "SNEAKILY_ADDED_GENE"]}, fh)

    real_sha256_file = G.sha256_file
    # Mutant: pretend every file hashes to whatever the manifest expects.
    monkeypatch.setattr(
        G,
        "sha256_file",
        lambda path: (
            manifest["cached_gene_sets_sha256"]
            if path == manifest["cached_gene_sets_path"]
            else real_sha256_file(path)
        ),
    )
    gene_sets, _ = G.load_cached_gene_sets(str(out_dir))
    # Did NOT raise, AND silently returned the tampered content -- proving the
    # real hash check (not this mutant) is what the real test above relies on.
    assert gene_sets["SET_A"] == ["G1", "G2", "SNEAKILY_ADDED_GENE"]

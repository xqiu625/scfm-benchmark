"""Local, versioned gene-set cache for the scFM discovery framework.

Per `review:20261012_153316` item C.4 and `review:20261012_171916`'s
next-increment instruction: imports a LOCAL `.gmt` file into a new, versioned cache
directory, recording collection name, species, gene-ID namespace, source
version/identifier, LOCAL import time (never a fabricated "retrieval" timestamp --
there is no network fetch here), and the source file's own content SHA-256.
`load_cached_gene_sets` re-verifies that hash before ever returning gene sets to a
caller, and raises -- it never silently falls back to a network fetch, because there
is no fallback path in this module at all.

Never queries a URL, never resolves an Enrichr/online library name. A `.gmt` file not
present on local disk is a `SystemExit`, not a retry-over-network.
"""
from __future__ import annotations

import hashlib
import json
import os
import time


CACHE_MANIFEST_SCHEMA = "scfm-gene-set-cache/1"

_CACHE_MANIFEST_REQUIRED_STR_FIELDS = (
    "schema", "collection_name", "species", "gene_id_namespace",
    "source_version", "source_type", "source_gmt_path_at_import",
    "source_content_sha256", "imported_utc", "cached_gene_sets_path",
    "cached_gene_sets_sha256",
)


def validate_cache_manifest(manifest: dict) -> list[str]:
    """Shared by `import_gene_sets` (write) and `load_cached_gene_sets`
    (read), per `review:20261012_220407` F3: a cache manifest that
    is missing a required field, carries a wrong-typed one, or claims an
    unrecognized schema must be rejected at BOTH boundaries, not accepted
    just because its gene-sets file's content hash happens to match -- a
    hash check proves content integrity, not metadata validity.
    """
    problems = []
    if not isinstance(manifest, dict):
        return [f"cache manifest must be a dict, got {type(manifest).__name__}"]
    if manifest.get("schema") != CACHE_MANIFEST_SCHEMA:
        problems.append(
            f"schema {manifest.get('schema')!r} != expected {CACHE_MANIFEST_SCHEMA!r}")
    for key in _CACHE_MANIFEST_REQUIRED_STR_FIELDS:
        if key not in manifest:
            problems.append(f"missing required key: {key!r}")
        elif not isinstance(manifest[key], str) or not manifest[key]:
            problems.append(
                f"{key!r} must be a non-empty string, got {manifest[key]!r}")
    if "n_gene_sets" not in manifest:
        problems.append("missing required key: 'n_gene_sets'")
    elif not isinstance(manifest["n_gene_sets"], int) or isinstance(manifest["n_gene_sets"], bool) \
            or manifest["n_gene_sets"] <= 0:
        problems.append(
            f"n_gene_sets must be a positive int, got {manifest.get('n_gene_sets')!r}")
    return problems


def sha256_file(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def parse_gmt(path: str) -> dict:
    """Parses a standard `.gmt` file: one gene set per line, tab-separated
    `name<TAB>description<TAB>gene1<TAB>gene2...`. Returns `{name: [genes]}`.
    `path` must be a local file -- this never follows a URL.
    """
    if not os.path.exists(path):
        raise SystemExit(f"gmt file not found: {path}")
    gene_sets: dict[str, list[str]] = {}
    with open(path) as fh:
        for line_no, raw_line in enumerate(fh, 1):
            line = raw_line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                raise SystemExit(
                    f"{path}:{line_no}: malformed GMT line (need name, description, "
                    f">=1 gene, got {len(parts)} field(s)): {line!r}")
            name, _description, *genes = parts
            genes = [g for g in genes if g]
            if not genes:
                raise SystemExit(f"{path}:{line_no}: gene set {name!r} has no genes")
            if name in gene_sets:
                raise SystemExit(f"{path}:{line_no}: duplicate gene set name {name!r}")
            gene_sets[name] = genes
    if not gene_sets:
        raise SystemExit(f"{path}: no gene sets parsed (empty file?)")
    return gene_sets


def import_gene_sets(*, source_gmt_path: str, collection_name: str, species: str,
                     gene_id_namespace: str, source_version: str,
                     out_dir: str) -> dict:
    """Imports a local `.gmt` file into a NEW, versioned cache directory under
    `out_dir`. Refuses to overwrite a non-empty `out_dir` -- same convention
    used throughout this project (e.g.
    `capture_preflight_GSM5102901_HC2.py`'s `run_preflight`). Returns the
    written `CACHE_MANIFEST.json` dict.
    """
    if not collection_name or not species or not gene_id_namespace or not source_version:
        raise SystemExit(
            "collection_name, species, gene_id_namespace, and source_version "
            "must all be non-empty -- none of them are inferred from the file")
    if os.path.exists(out_dir) and os.listdir(out_dir):
        raise SystemExit(f"{out_dir} must be new or empty")

    gene_sets = parse_gmt(source_gmt_path)
    source_content_sha256 = sha256_file(source_gmt_path)

    os.makedirs(out_dir, exist_ok=True)
    cached_gene_sets_filename = "gene_sets.json"
    cached_gene_sets_path = os.path.join(out_dir, cached_gene_sets_filename)
    with open(cached_gene_sets_path, "w") as fh:
        json.dump(gene_sets, fh, indent=2, sort_keys=True)
    cached_gene_sets_sha256 = sha256_file(cached_gene_sets_path)

    manifest = {
        "schema": CACHE_MANIFEST_SCHEMA,
        "collection_name": collection_name,
        "species": species,
        "gene_id_namespace": gene_id_namespace,
        "source_version": source_version,
        "source_type": "local_file",
        "source_gmt_path_at_import": os.path.abspath(source_gmt_path),
        "source_content_sha256": source_content_sha256,
        # LOCAL import time -- not a live-fetch "retrieved_utc"; there was no
        # network request to retrieve anything.
        "imported_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "n_gene_sets": len(gene_sets),
        # Stored RELATIVE to this manifest's own directory (per
        # review:20261012_220407 F4) -- resolved against the
        # manifest's location at load time, never the process's cwd at
        # import time, so a cache directory copied or loaded from a
        # different cwd still resolves correctly.
        "cached_gene_sets_path": cached_gene_sets_filename,
        "cached_gene_sets_sha256": cached_gene_sets_sha256,
    }
    problems = validate_cache_manifest(manifest)
    if problems:
        raise ValueError(f"import_gene_sets produced an invalid manifest: {problems}")
    manifest_path = os.path.join(out_dir, "CACHE_MANIFEST.json")
    tmp = manifest_path + ".inprogress"
    with open(tmp, "w") as fh:
        json.dump(manifest, fh, indent=2, sort_keys=True)
        fh.flush()
        os.fsync(fh.fileno())
    os.replace(tmp, manifest_path)
    # Returned to the caller with an ABSOLUTE path for immediate convenience
    # (e.g. printing in an example) -- the ON-DISK manifest keeps the
    # relative form above; load_cached_gene_sets re-resolves it itself.
    returned = dict(manifest)
    returned["cached_gene_sets_path"] = cached_gene_sets_path
    return returned


def load_cached_gene_sets(cache_dir: str) -> tuple[dict, dict]:
    """Reads back a previously-imported cache, VERIFYING the stored gene-sets
    file's hash against the manifest's own recorded hash before returning
    anything. A missing manifest, a missing gene-sets file, or a hash
    mismatch (corruption, hand-editing, or a swapped file) is a `SystemExit`
    -- there is no fallback to a live fetch anywhere in this function.

    Returns `(gene_sets_dict, cache_manifest_dict)`.
    """
    manifest_path = os.path.join(cache_dir, "CACHE_MANIFEST.json")
    if not os.path.exists(manifest_path):
        raise SystemExit(f"no CACHE_MANIFEST.json in {cache_dir}")
    with open(manifest_path) as fh:
        manifest = json.load(fh)
    problems = validate_cache_manifest(manifest)
    if problems:
        raise SystemExit(f"cache manifest at {manifest_path} failed validation: {problems}")
    recorded_path = manifest["cached_gene_sets_path"]
    # Resolved relative to the MANIFEST's own directory (per
    # review:20261012_220407 F4), never the caller's cwd -- an
    # already-absolute recorded path (older caches, or a caller who passed
    # one explicitly) is still honored as-is.
    gene_sets_path = (recorded_path if os.path.isabs(recorded_path)
                      else os.path.join(cache_dir, recorded_path))
    if not os.path.exists(gene_sets_path):
        raise SystemExit(f"cached gene-sets file missing: {gene_sets_path}")
    actual_hash = sha256_file(gene_sets_path)
    expected_hash = manifest.get("cached_gene_sets_sha256")
    if actual_hash != expected_hash:
        raise SystemExit(
            f"cached gene-sets file hash mismatch at {gene_sets_path}: expected "
            f"{expected_hash}, got {actual_hash} -- cache may be corrupted or "
            "tampered; refusing to use it (no fallback to a live fetch exists)")
    with open(gene_sets_path) as fh:
        gene_sets = json.load(fh)
    # The manifest returned to the CALLER carries the fully-RESOLVED absolute
    # path actually used and verified above -- useful for a result's own
    # provenance record -- while the ON-DISK file keeps the relative,
    # cwd/location-independent form.
    manifest_for_caller = dict(manifest)
    manifest_for_caller["cached_gene_sets_path"] = gene_sets_path
    return gene_sets, manifest_for_caller

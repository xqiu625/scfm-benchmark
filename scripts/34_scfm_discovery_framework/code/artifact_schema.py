"""Artifact/provenance manifest schema for the scFM discovery framework.

Per `review:20261012_145004` item C.3: input identity, source-code/config
hash, model identity or pending, stage status, expected outputs, dependencies.
`planned`/`executed`/`validated` are kept as three distinct states on purpose --
nothing in this module can produce `executed` or `validated` without a caller
explicitly asserting a real run happened (there is no code path from `plan()` to
either of those states). Mirrors the atomic-write convention already used in
`scripts/32_rebuild_upfront_qc/hc2_preflight/capture_preflight_GSM5102901_HC2.py`.
"""
from __future__ import annotations

import hashlib
import json
import os
import time

SCHEMA_VERSION = "scfm-discovery-artifact/1"

#: `planned`: a plan artifact was produced; nothing ran. `executed`: a real run
#: happened and produced raw output, not yet checked. `validated`: a real run's
#: output was independently re-checked against its own stated expectations.
#: There is deliberately no automatic path from `planned` to either later state.
STAGE_STATUSES = ("planned", "executed", "validated")

#: Top-level keys every manifest must have, independent of stage_status.
REQUIRED_MANIFEST_KEYS = ("schema", "generated_utc", "input_identity",
                          "source_code_revision_or_hash",
                          "model_identity_or_pending", "stage_status",
                          "expected_outputs", "dependencies")


def find_nonfinite_floats(obj, path: str = "<root>") -> list[str]:
    """Recursively finds any non-finite float anywhere in a nested dict/list
    structure. Catches BOTH a literal NaN/Infinity JSON token (already
    intercepted earlier by `read_manifest`'s `parse_constant`) AND a normal
    numeric literal that overflows to `inf` during parsing (e.g. `1e999`,
    which `parse_constant` never sees because it is not one of the three
    special tokens) -- per `review:20261012_171916` R1: the reader
    accepted `1e999` and produced `inf`, which `write_atomic`'s
    `allow_nan=False` would then refuse, an asymmetry between what read
    accepts and what write accepts.
    """
    import math
    problems = []
    if isinstance(obj, float):
        if not math.isfinite(obj):
            problems.append(f"{path}: non-finite float {obj!r}")
    elif isinstance(obj, dict):
        for k, v in obj.items():
            problems.extend(find_nonfinite_floats(v, f"{path}.{k}"))
    elif isinstance(obj, list):
        for i, v in enumerate(obj):
            problems.extend(find_nonfinite_floats(v, f"{path}[{i}]"))
    return problems


def validate_manifest(manifest: dict) -> list[str]:
    """Structural validation shared by build/write/read, per
    `review:20261012_153316` R3 and `review:20261012_163250`
    finding 3: a manifest missing input/source/model/output fields, carrying
    an illegal `stage_status`, or carrying a wrong-typed
    `source_code_revision_or_hash`/`dependencies` entry, must be rejected at
    every one of those boundaries, not merely at construction time (which a
    tampered on-disk file bypasses entirely). Also rejects ANY non-finite
    float anywhere in the manifest (not just the top-level literal constant
    case), per `review:20261012_171916` R1.
    """
    problems = []
    if not isinstance(manifest, dict):
        return [f"manifest must be a dict, got {type(manifest).__name__}"]
    problems.extend(find_nonfinite_floats(manifest))
    for key in REQUIRED_MANIFEST_KEYS:
        if key not in manifest:
            problems.append(f"missing required key: {key!r}")
    if "schema" in manifest and manifest["schema"] != SCHEMA_VERSION:
        problems.append(
            f"schema {manifest['schema']!r} != expected {SCHEMA_VERSION!r}")
    if "stage_status" in manifest and manifest["stage_status"] not in STAGE_STATUSES:
        problems.append(
            f"stage_status {manifest['stage_status']!r} not in {STAGE_STATUSES}")
    if "input_identity" in manifest and not isinstance(manifest["input_identity"], dict):
        problems.append("input_identity must be a dict")
    if "input_identity" in manifest and isinstance(manifest["input_identity"], dict) \
            and not manifest["input_identity"]:
        problems.append("input_identity must not be empty")
    mip = manifest.get("model_identity_or_pending")
    if "model_identity_or_pending" in manifest and mip != "model_undecided" \
            and not isinstance(mip, dict):
        problems.append(
            "model_identity_or_pending must be the literal string "
            "'model_undecided' or a dict")
    if "source_code_revision_or_hash" in manifest \
            and not isinstance(manifest["source_code_revision_or_hash"], str):
        problems.append(
            f"source_code_revision_or_hash must be a string, got "
            f"{type(manifest['source_code_revision_or_hash']).__name__}")
    if "expected_outputs" in manifest:
        if not isinstance(manifest["expected_outputs"], list):
            problems.append("expected_outputs must be a list")
        elif not all(isinstance(x, str) for x in manifest["expected_outputs"]):
            problems.append("expected_outputs entries must all be strings")
    if "dependencies" in manifest:
        deps = manifest["dependencies"]
        if not isinstance(deps, list):
            problems.append("dependencies must be a list")
        else:
            for i, dep in enumerate(deps):
                if not isinstance(dep, dict):
                    problems.append(
                        f"dependencies[{i}] must be a dict with 'file' and "
                        f"'sha256' string keys, got {type(dep).__name__}")
                    continue
                for required_key in ("file", "sha256"):
                    if not isinstance(dep.get(required_key), str):
                        problems.append(
                            f"dependencies[{i}][{required_key!r}] must be a "
                            f"string, got {type(dep.get(required_key)).__name__}")
    return problems


def _reject_nonfinite_json_constant(token: str):
    raise ValueError(
        f"non-finite JSON constant {token!r} is not accepted -- a manifest "
        "must not persist NaN/Infinity, matching write_atomic's allow_nan=False")


def sha256_file(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode()).hexdigest()


def safe_json_identity(value) -> str | None:
    """Strict identity hash: succeeds only when `value` serializes with the
    standard JSON encoder's NATIVE type support (str/int/float/bool/None/
    list/dict) -- never silently stringifies an unrecognized object via
    `default=str`. Returns `None`, never raises, when the value cannot be
    hashed this way; callers must treat `None` as "identity unavailable",
    not as an empty-string hash of nothing. Per
    `review:20261013_222450` R3; shared here (rather than defined
    separately in `gene_ranking_providers.py` and `pseudobulk_providers.py`)
    so neither evidence stream duplicates this policy.
    """
    try:
        return sha256_text(json.dumps(value, sort_keys=True))
    except TypeError:
        return None


def array_content_sha256(arr) -> str:
    """Content hash of a numpy array's actual bytes plus its shape -- binds
    an on-disk JSON record's claims ABOUT an array (row/column labels,
    derived per-row/per-column sums) to the array's REAL content, not just
    to the record's own internal self-consistency. Per
    `review:20261014_100436` F1: a hash of the npz FILE alone
    cannot catch a record whose row/column LABELS were reordered without
    touching the array bytes; `array_content_sha256` is combined with a
    recompute-from-the-loaded-array check at the actual read boundary
    (see `pseudobulk_providers.load_pseudobulk_counts`), not merely stored
    and left unread. `str(arr.shape)` is folded into the hash so two
    differently-shaped arrays with coincidentally identical bytes (e.g. an
    all-zero `(2, 3)` vs `(3, 2)` array) still hash differently.
    """
    import numpy as np
    arr = np.ascontiguousarray(arr, dtype=np.float64)
    return hashlib.sha256(str(arr.shape).encode() + b"|" + arr.tobytes()).hexdigest()


def resolve_source_revision(here: str, source_modules: tuple) -> str:
    """Prefers the real git commit of the checkout `here` lives in; falls back
    to a content hash of `source_modules` if git is unavailable -- never the
    literal string "unknown" when either is obtainable, since the running
    source is right there on disk. Factored out of `cli.py`'s original
    identically-named function so `gsea_runner.py` can bind its own source
    revision/dependency identity (per `review:20261012_220407` F2)
    without duplicating this logic or importing the `cli` module.
    """
    import subprocess
    try:
        out = subprocess.run(["git", "-C", here, "rev-parse", "HEAD"],
                             capture_output=True, text=True, timeout=10)
        if out.returncode == 0:
            commit = out.stdout.strip()
            dirty = subprocess.run(["git", "-C", here, "status", "--porcelain",
                                    "--", *source_modules],
                                   capture_output=True, text=True, timeout=10)
            suffix = "" if dirty.stdout.strip() == "" else "+dirty"
            return f"git:{commit}{suffix}"
    except Exception:
        pass
    hashes = []
    for name in source_modules:
        path = os.path.join(here, name)
        if os.path.exists(path):
            hashes.append(sha256_file(path))
        else:
            hashes.append("missing:" + name)
    return "content_hash:" + sha256_text("\n".join(hashes))


def source_module_dependencies(here: str, source_modules: tuple) -> list[dict]:
    deps = []
    for name in source_modules:
        path = os.path.join(here, name)
        deps.append({
            "file": name,
            "sha256": sha256_file(path) if os.path.exists(path) else "missing",
        })
    return deps


def build_manifest(*, input_identity: dict, source_code_revision_or_hash: str,
                    model_identity_or_pending, stage_status: str,
                    expected_outputs: list[str], dependencies: list,
                    extra: dict | None = None) -> dict:
    """`model_identity_or_pending` is either a dict (a resolved model registry
    entry, see model_registry.py) or the literal string "model_undecided" --
    never a guessed/default model.
    """
    if stage_status not in STAGE_STATUSES:
        raise ValueError(f"stage_status={stage_status!r} not in {STAGE_STATUSES}")
    if model_identity_or_pending != "model_undecided" and not isinstance(model_identity_or_pending, dict):
        raise ValueError(
            "model_identity_or_pending must be the literal string 'model_undecided' "
            "or a dict describing a resolved model registry entry")
    manifest = {
        "schema": SCHEMA_VERSION,
        "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "input_identity": input_identity,
        "source_code_revision_or_hash": source_code_revision_or_hash,
        "model_identity_or_pending": model_identity_or_pending,
        "stage_status": stage_status,
        "expected_outputs": list(expected_outputs),
        "dependencies": list(dependencies),
    }
    if extra:
        manifest["extra"] = extra
    problems = validate_manifest(manifest)
    if problems:
        raise ValueError(f"build_manifest produced an invalid manifest: {problems}")
    return manifest


def write_atomic(path: str, payload: dict) -> None:
    """No `default=` fallback: an unsupported JSON value (e.g. a NaN float, a
    set, a custom object) raises here, loudly, rather than being silently
    stringified -- DEVELOPMENT.md's "round trips must be strict" standing rule.
    Validates the payload's own structure before writing, so a caller cannot
    accidentally persist something `read_manifest` would then have to reject.
    """
    problems = validate_manifest(payload)
    if problems:
        raise ValueError(f"refusing to write an invalid manifest: {problems}")
    tmp = path + ".inprogress"
    with open(tmp, "w") as fh:
        json.dump(payload, fh, indent=2, sort_keys=True, allow_nan=False)
        fh.flush()
        os.fsync(fh.fileno())
    os.replace(tmp, path)


def read_manifest(path: str) -> dict:
    """Re-validates the FULL manifest structure on read, not just the schema
    version string -- a hand-edited or tampered on-disk file with an illegal
    `stage_status`, missing required keys, or a wrong-typed field is rejected
    here too, not only at `build_manifest()` construction time. Also rejects
    NaN/Infinity JSON constants on parse (Python's `json.load` accepts them
    by default even though they are not standard JSON and `write_atomic`
    refuses to write them) -- reader and writer agree on this, per
    `review:20261012_163250` finding 3.
    """
    if not os.path.exists(path):
        raise SystemExit(f"manifest not found: {path}")
    with open(path) as fh:
        try:
            manifest = json.load(fh, parse_constant=_reject_nonfinite_json_constant)
        except ValueError as e:
            raise SystemExit(f"manifest at {path} failed to parse: {e}")
    problems = validate_manifest(manifest)
    if problems:
        raise SystemExit(f"manifest at {path} failed validation: {problems}")
    return manifest

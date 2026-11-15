"""`validate` / `plan` entry points for the scFM discovery framework.

Per `review:20261012_145004` item C.2 and `review:20261012_153316`
R1/R2: validates config/metadata with explicit error reasons, actually checks the
declared matrix structure (not just key presence), and produces a reviewable plan
artifact with a real config snapshot and real source-code identity. Unknown model
and unverified species compatibility are reported explicitly, not silently passed.
A plan may be produced before a model is pinned (`model_undecided`), but it is
never labeled executable.

This CLI never submits a job, downloads a weight, connects to a remote database, or
reads a full real expression matrix. If an h5ad's STRUCTURE needs checking, it is
read with bounded `h5py` access to `obs`/`var`/declared-slot shapes and the one
named cell-ID column's values -- never `anndata.read_h5ad(..., backed="r")`, which
the batch-effect audit work in this project already found still loads every layer
in full, and never `X`'s actual numeric contents. A structural pass here is NOT a
numeric QC pass -- it never has been and this round does not change that.
"""
from __future__ import annotations

import argparse
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from artifact_schema import build_manifest, sha256_text, write_atomic  # noqa: E402
from artifact_schema import resolve_source_revision as _resolve_source_revision  # noqa: E402
from artifact_schema import source_module_dependencies as _source_module_dependencies  # noqa: E402
from dataset_schema import (NO_BATCH_SENTINEL, load_dataset_config,  # noqa: E402
                            validate_dataset_config)
from model_registry import REGISTRY, UNKNOWN  # noqa: E402


def resolve_matrix_path(dataset_config_path: str, config: dict) -> str | None:
    """A relative `expression_matrix_path` is resolved against the dataset
    config FILE's own directory, never the process's current working
    directory -- per `review:20261012_163250` finding 2, the same
    config run from two different cwds must resolve to the same file. An
    already-absolute path is returned unchanged. The config's own declared
    value is never mutated; this is a separate, derived value for file
    access and provenance recording.

    Raises `SystemExit` (a clean, structured failure) if `expression_matrix_path`
    is present but not a string, rather than letting `os.path.isabs` crash with
    a raw `TypeError` -- per `review:20261012_171916` R2. Normal
    callers should never actually hit this: `validate_dataset_config` already
    rejects a non-string value, and `cmd_validate`/`cmd_plan` only call this
    after confirming that check passed. This is defense in depth, not the
    primary fix (see the control-flow fix in `cmd_validate`).
    """
    raw = config.get("expression_matrix_path")
    if not raw:
        return None
    if not isinstance(raw, str):
        raise SystemExit(
            f"expression_matrix_path must be a string, got {type(raw).__name__}: {raw!r}")
    if os.path.isabs(raw):
        return raw
    config_dir = os.path.dirname(os.path.abspath(dataset_config_path))
    return os.path.join(config_dir, raw)

#: The modules a `plan`/`validate` run actually executes -- their combined content
#: hash is `source_code_revision_or_hash` when git is unavailable. Listed once here
#: so the dependency list in build_manifest() and the hash computation cannot drift
#: apart from each other.
SOURCE_MODULES = ("cli.py", "dataset_schema.py", "model_registry.py", "artifact_schema.py")


def _group_or_dataset_shape(node) -> list[int] | None:
    """Metadata-only shape of an h5py dataset OR a sparse-matrix group (which
    stores `shape` as an attr, not as `.shape`) -- never reads `data`/`indices`.
    """
    import h5py
    if isinstance(node, h5py.Dataset):
        return list(node.shape)
    if "shape" in node.attrs:
        return [int(v) for v in node.attrs["shape"]]
    return None


def resolve_counts_slot_shape(f, counts_slot: str) -> list[int] | None:
    """Resolves `"X"` / `"raw.X"` / `"layers/<name>"` to that slot's declared
    shape. Returns None if the slot does not exist in the file at all.
    """
    if counts_slot == "X":
        path = ("X",)
    elif counts_slot == "raw.X":
        path = ("raw", "X")
    elif counts_slot.startswith("layers/"):
        path = ("layers", counts_slot[len("layers/"):])
    else:
        return None
    node = f
    for part in path:
        if part not in node:
            return None
        node = node[part]
    return _group_or_dataset_shape(node)


def _index_length(f, group_path: str) -> int | None:
    """Metadata-only length of an obs/var-shaped group's own index -- shape of
    the index dataset (or its categorical `codes`), never its values.
    `group_path` may be a nested h5py path (e.g. `"raw/var"`) -- h5py resolves
    `/`-separated keys directly, so no manual traversal is needed.
    """
    if group_path not in f:
        return None
    group = f[group_path]
    index_name = group.attrs.get("_index")
    if index_name is None or index_name not in group:
        return None
    index_node = group[index_name]
    import h5py
    if isinstance(index_node, h5py.Group):
        return int(index_node["codes"].shape[0])
    return int(index_node.shape[0])


def gene_axis_group_for_slot(counts_slot: str) -> str | None:
    """Maps a counts_slot to the h5ad group that is ITS gene axis: `"X"` and
    `"layers/<name>"` are indexed by the current, possibly-filtered `var`;
    `"raw.X"` is indexed by `raw.var`, which legitimately can (and often does)
    retain more genes than the current `var` after a filtering step -- they
    must never be compared against each other
    (https://anndata.readthedocs.io/en/stable/generated/anndata.AnnData.raw.html).
    Returns None for an unrecognized slot (dataset_schema already restricts
    counts_slot's format, so this should not happen for a config that passed
    validate_dataset_config).
    """
    if counts_slot == "raw.X":
        return "raw/var"
    if counts_slot == "X" or counts_slot.startswith("layers/"):
        return "var"
    return None


def check_duplicate_cell_ids(cell_ids: list) -> list[str]:
    seen = set()
    dupes = set()
    for c in cell_ids:
        if c in seen:
            dupes.add(c)
        seen.add(c)
    return sorted(dupes)


def read_obs_column_values(path: str, column: str) -> list:
    """Reads ONE named obs column's values -- cell-ID-scale metadata, not the
    expression matrix. Handles both a plain dataset and a pandas categorical
    group (codes + categories), the two forms anndata actually writes.
    """
    import h5py
    with h5py.File(path, "r") as f:
        if "obs" not in f or column not in f["obs"]:
            raise SystemExit(f"obs column not found: {column}")
        node = f["obs"][column]
        if isinstance(node, h5py.Group):
            codes = node["codes"][()]
            categories = node["categories"][()]
            out = []
            for c in codes:
                if c < 0:
                    out.append(None)
                else:
                    v = categories[c]
                    out.append(v.decode() if isinstance(v, bytes) else v)
            return out
        values = node[()]
        return [v.decode() if isinstance(v, bytes) else v for v in values]


def validate_matrix_structure(path: str, config: dict) -> list[str]:
    """Real structural validation: the declared counts_slot must exist, be
    2-D, and its (n_obs, n_vars) must agree with obs's own index length AND
    the index length of ITS OWN gene axis (`var` for `"X"`/`"layers/*"`,
    `raw/var` for `"raw.X"` -- never compared against the wrong one: `raw.X`
    legitimately keeps more genes than the current, possibly-filtered `var`,
    per `review:20261012_163250` finding 1); every declared
    obs-identity column (subject, sample, biological-replicate, cell-ID, and
    technical-batch unless it is explicitly `NO_BATCH_SENTINEL`) must
    actually be an obs column; the declared cell-ID column's values must be
    non-empty AND unique. An index that cannot be resolved (missing/
    unresolvable `_index` attribute) is reported as an explicit failure, not
    silently skipped -- an unverifiable axis is never treated as a passed
    check.

    Returns a list of problem strings (empty = no structural problems found).
    Still never reads X's/the slot's actual numeric contents -- this checks
    shape and column presence only, not counts.
    """
    import h5py
    problems: list[str] = []
    if not os.path.exists(path):
        return [f"expression_matrix_path does not exist: {path}"]

    with h5py.File(path, "r") as f:
        n_obs_from_obs = None
        if "obs" not in f:
            problems.append("no 'obs' group in the file")
        else:
            n_obs_from_obs = _index_length(f, "obs")
            if n_obs_from_obs is None:
                problems.append(
                    "could not resolve the 'obs' index (missing or unresolvable "
                    "'_index' attribute) -- cell-count consistency cannot be verified")

        counts_slot = config.get("counts_slot")
        var_group_path = gene_axis_group_for_slot(counts_slot) if counts_slot else None
        n_vars_from_var = None
        if var_group_path is None:
            problems.append(f"counts_slot={counts_slot!r} has no known gene-axis "
                             "mapping -- cannot determine which var group to check")
        elif var_group_path not in f:
            problems.append(
                f"counts_slot={counts_slot!r} implies its gene axis lives at "
                f"{var_group_path!r}, which does not exist in the file")
        else:
            n_vars_from_var = _index_length(f, var_group_path)
            if n_vars_from_var is None:
                problems.append(
                    f"could not resolve the {var_group_path!r} index (missing or "
                    "unresolvable '_index' attribute) -- counts_slot shape "
                    "consistency cannot be verified")

        slot_shape = resolve_counts_slot_shape(f, counts_slot) if counts_slot else None
        if counts_slot and slot_shape is None:
            problems.append(f"counts_slot={counts_slot!r} does not exist in the file "
                             "-- there is no expression matrix at the declared slot")
        elif counts_slot and len(slot_shape) != 2:
            problems.append(f"counts_slot={counts_slot!r} has shape {slot_shape}, "
                             "not 2-D")
        elif counts_slot:
            n_obs_x, n_vars_x = slot_shape
            if n_obs_from_obs is not None and n_obs_x != n_obs_from_obs:
                problems.append(
                    f"counts_slot={counts_slot!r} shape[0]={n_obs_x} disagrees with "
                    f"obs index length {n_obs_from_obs}")
            if n_vars_from_var is not None and n_vars_x != n_vars_from_var:
                problems.append(
                    f"counts_slot={counts_slot!r} shape[1]={n_vars_x} disagrees with "
                    f"{var_group_path!r} index length {n_vars_from_var}")

        if "obs" in f:
            obs_keys = set(f["obs"].keys())
            identity_fields = ["subject_id_column", "sample_id_column",
                               "biological_replicate_unit", "cell_id_column"]
            for field_name in identity_fields:
                col = config.get(field_name)
                if col and col not in obs_keys:
                    problems.append(
                        f"{field_name}={col!r} not found among obs keys "
                        f"{sorted(obs_keys)}")

            batch_col = config.get("technical_batch_column")
            if batch_col and batch_col != NO_BATCH_SENTINEL and batch_col not in obs_keys:
                problems.append(
                    f"technical_batch_column={batch_col!r} not found among obs keys "
                    f"{sorted(obs_keys)} (use {NO_BATCH_SENTINEL!r} explicitly if this "
                    "dataset genuinely has no batch variation to declare)")

            cell_id_col = config.get("cell_id_column")
            if cell_id_col and cell_id_col in obs_keys:
                values = read_obs_column_values(path, cell_id_col)
                n_empty = sum(1 for v in values if v is None or v == "")
                if n_empty:
                    problems.append(
                        f"cell_id_column={cell_id_col!r} has {n_empty} missing/empty "
                        "value(s) -- every cell must have a real identity")
                dupes = check_duplicate_cell_ids(values)
                if dupes:
                    problems.append(
                        f"cell_id_column={cell_id_col!r} has {len(dupes)} duplicate "
                        f"value(s): {dupes[:5]}" + (" (truncated)" if len(dupes) > 5 else ""))

    return problems


def resolve_source_revision() -> str:
    """Thin wrapper over `artifact_schema.resolve_source_revision`, bound to
    this CLI's own HERE/SOURCE_MODULES. See that function's docstring.
    """
    return _resolve_source_revision(HERE, SOURCE_MODULES)


def source_module_dependencies() -> list[dict]:
    """Thin wrapper over `artifact_schema.source_module_dependencies`, bound
    to this CLI's own HERE/SOURCE_MODULES."""
    return _source_module_dependencies(HERE, SOURCE_MODULES)


def _model_report(config: dict, model_registry_key: str | None):
    """Shared by validate() and plan() so plan() cannot silently diverge from
    validate()'s species-compatibility check (R2's reproduced finding)."""
    if not model_registry_key:
        return None, []
    entry = REGISTRY.get(model_registry_key)
    if entry is None:
        return None, [f"model_registry_key={model_registry_key!r} not found in "
                       f"REGISTRY ({sorted(REGISTRY)})"]
    species = config.get("species")
    evidence = entry.species_pretraining_evidence.get(species, UNKNOWN)
    report = {
        "model_registry_key": model_registry_key,
        "model_generation": entry.model_generation,
        "species_requested": species,
        "species_pretraining_evidence": evidence,
    }
    soft = []
    if evidence not in ("documented_confirmed", "checkpoint_metadata_confirmed"):
        soft.append(
            f"model {entry.model_generation!r} species-pretraining evidence for "
            f"{species!r} is {evidence!r}, not confirmed -- reported explicitly, "
            "not treated as a validation failure by itself (a model MAY still be "
            "used deliberately with unconfirmed species evidence, but "
            "validate()/plan() must say so)")
    return report, soft


def _run_structural_checks(config: dict, check_matrix_structure: bool,
                           resolved_matrix_path: str | None) -> tuple[list[str], list[str]]:
    """Returns (hard_issues, checked_items) -- checked_items documents what
    was and was not actually inspected, so a plan artifact can state its own
    validation scope rather than implying everything was checked.
    `resolved_matrix_path` must already be resolved (see resolve_matrix_path)
    -- this function never re-derives it from a bare cwd-relative value.
    """
    if not check_matrix_structure:
        return [], ["dataset config fields only -- matrix structure NOT checked "
                    "(pass --check-matrix-structure to check it)"]
    if not resolved_matrix_path:
        return [], ["expression_matrix_path missing from config -- nothing to check"]
    problems = validate_matrix_structure(resolved_matrix_path, config)
    counts_slot = config.get("counts_slot")
    var_group_path = gene_axis_group_for_slot(counts_slot) if counts_slot else None
    checked = ["obs group presence and resolvable index",
               f"gene axis for counts_slot={counts_slot!r} resolved to "
               f"{var_group_path!r} and its index checked for resolvability",
               "counts_slot existence/shape vs. obs index length AND its OWN gene "
               "axis's index length (never the wrong var group)",
               "subject/sample/biological-replicate/cell-ID columns present in obs",
               "technical_batch_column present in obs (or explicit NO_BATCH_SENTINEL)",
               "cell_id_column values non-empty and unique",
               "NOTE: counts_slot's actual numeric values were never read (no "
               "numeric QC is performed by this check)"]
    return problems, checked


def cmd_validate(args) -> int:
    config = load_dataset_config(args.dataset_config)
    result = validate_dataset_config(config)

    # resolve_matrix_path() is only called once the dataset config itself is
    # known-valid (in particular, expression_matrix_path is confirmed to be a
    # string) -- calling it unconditionally let a type-invalid config (e.g.
    # expression_matrix_path=True) crash with a raw TypeError instead of a
    # structured report, per review:20261012_171916 R2. plan()
    # already had this control flow (it returns early on `not result.ok`);
    # validate() now matches it.
    resolved_matrix_path = None
    structural_issues, checked_items = ([], [])
    if result.ok:
        resolved_matrix_path = resolve_matrix_path(args.dataset_config, config)
        structural_issues, checked_items = _run_structural_checks(
            config, args.check_matrix_structure, resolved_matrix_path)

    model_report, soft_issues = _model_report(config, args.model_registry_key)
    if model_report is None and args.model_registry_key:
        # _model_report returned an error list in soft_issues' slot via the
        # "not found" path -- treat that specific case as a hard issue.
        structural_issues = structural_issues + soft_issues
        soft_issues = []

    hard_issues = [str(i) for i in result.issues] + structural_issues
    ok = len(hard_issues) == 0

    report = {
        "ok": ok,
        "dataset_config_issues": [str(i) for i in result.issues],
        "expression_matrix_path_declared": config.get("expression_matrix_path"),
        "expression_matrix_path_resolved": resolved_matrix_path,
        "matrix_structure_issues": structural_issues,
        "checked_items": checked_items,
        "model_report": model_report,
        "soft_issues": soft_issues,
    }
    print(json.dumps(report, indent=2))
    return 0 if ok else 1


def cmd_plan(args) -> int:
    config = load_dataset_config(args.dataset_config)
    result = validate_dataset_config(config)
    if not result.ok:
        print(json.dumps({"ok": False, "issues": [str(i) for i in result.issues]},
                          indent=2))
        return 1

    resolved_matrix_path = resolve_matrix_path(args.dataset_config, config)
    structural_issues, checked_items = _run_structural_checks(
        config, args.check_matrix_structure, resolved_matrix_path)
    model_report, soft_issues = _model_report(config, args.model_registry_key)
    if model_report is None and args.model_registry_key:
        print(json.dumps({"ok": False, "issues": soft_issues}, indent=2))
        return 1

    model_identity = (REGISTRY[args.model_registry_key].to_dict()
                      if args.model_registry_key else "model_undecided")

    config_text = json.dumps(config, sort_keys=True)
    manifest = build_manifest(
        input_identity={
            "dataset_id": config["dataset_id"],
            "dataset_config_path": os.path.abspath(args.dataset_config),
            "dataset_config_sha256": sha256_text(config_text),
            "dataset_config_snapshot": config,
            "expression_matrix_path_resolved": resolved_matrix_path,
        },
        source_code_revision_or_hash=resolve_source_revision(),
        model_identity_or_pending=model_identity,
        stage_status="planned",
        expected_outputs=["cells.h5ad or embeddings (not produced by plan())",
                          "gene_rankings.tsv (not produced by plan())",
                          "pathway_results.tsv (not produced by plan())"],
        dependencies=source_module_dependencies(),
        extra={
            "note": "plan() never executes inference, downloads weights, submits "
                    "a job, or connects to a network database.",
            "executable": False,
            "validation_scope": {
                "matrix_structure_checked": args.check_matrix_structure,
                "matrix_structure_issues": structural_issues,
                "checked_items": checked_items,
            },
            "model_compatibility": model_report,
            "model_compatibility_open_items": soft_issues,
        },
    )
    if structural_issues:
        print(json.dumps({"ok": False, "issues": structural_issues}, indent=2))
        return 1
    if os.path.exists(args.out) and os.listdir(args.out):
        raise SystemExit(f"{args.out} must be new or empty")
    os.makedirs(args.out, exist_ok=True)
    dest = os.path.join(args.out, "PLAN.json")
    write_atomic(dest, manifest)
    print(json.dumps(manifest, indent=2))
    print(f"plan written: {dest}")
    return 0


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    sub = ap.add_subparsers(dest="command", required=True)

    p_validate = sub.add_parser("validate")
    p_validate.add_argument("--dataset-config", required=True)
    p_validate.add_argument("--model-registry-key", default=None)
    p_validate.add_argument("--check-matrix-structure", action="store_true")
    p_validate.set_defaults(func=cmd_validate)

    p_plan = sub.add_parser("plan")
    p_plan.add_argument("--dataset-config", required=True)
    p_plan.add_argument("--model-registry-key", default=None)
    p_plan.add_argument("--check-matrix-structure", action="store_true")
    p_plan.add_argument("--out", required=True)
    p_plan.set_defaults(func=cmd_plan)

    args = ap.parse_args()
    sys.exit(args.func(args))


if __name__ == "__main__":
    main()

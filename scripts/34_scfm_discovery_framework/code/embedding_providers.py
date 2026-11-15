"""Embedding provider contracts for the scFM discovery framework.

Per `review:20261013_135634` section 3: a PCA baseline is a
CONTROL -- an explicit "what can a simple expression representation already
do on this input/downstream question" comparison -- it is never treated as
proof that a scFM is superior, and a PCA loading/dimension correlation is
never passed off as model-masking evidence (that is `gene_ranking_providers
.ModelMaskingRankingProvider`'s separate, still-not-implemented contract).

`PCAEmbeddingProvider` is a REAL, scikit-learn-backed implementation (not a
contract stub) that reuses this framework's own input validation
(`dataset_schema.validate_dataset_config`) and structural matrix checks
(`cli.resolve_matrix_path`, `cli.validate_matrix_structure`,
`cli.gene_axis_group_for_slot`) -- it is never handed a bare file path with
`counts_slot`/gene-axis/counts-semantics stripped away. Unlike `cli.py`
(which deliberately never reads a matrix's actual numeric contents -- see
that module's own docstring), this module DOES read real numeric values,
because computing a PCA embedding is exactly the point; the scope that read
is bounded to is stated explicitly below (`_MAX_DENSE_ELEMENTS`), and an
input outside that scope is REFUSED, never silently truncated.

Never modifies the input expression matrix (opened read-only throughout).
Never treats PCA as a batch-correction step. Never uses a disease/outcome
label to select features -- there is no such parameter on this interface at
all; the only feature-selection paths are "all genes in the counts_slot",
an explicit caller-supplied `feature_ids_override`, or a label-free
top-variance `hvg_top_n` selection.

Per `review:20261013_210131` F2: the counts_slot's own gene axis
(`var`/`raw/var`) is validated -- non-empty, every value a non-empty
string, and UNIQUE -- before any ID->position map, `feature_ids_override`,
or `hvg_top_n` selection is built from it. A duplicated gene ID is refused
with a clean `failed` record; it is never silently resolved by
`.var_names_make_unique()`-style suffixing, by keeping the first/last
occurrence, or by averaging/merging the duplicated columns.

Per `review:20261013_210131` F3: every array this module is about
to PERSIST (embedding, components, mean, singular values, AND explained
variance ratio) is checked for finiteness before anything is written --
not just the three checked before this round, which let a zero-total-
variance input (e.g. an all-zero counts block) reach `_write_embedding_
result` with a NaN `explained_variance_ratio_`, write `embedding_matrices
.npz`, and only then raise on the JSON manifest's own validation, leaving
an orphaned npz with no `EMBEDDING_RESULT.json` beside it. A zero-total-
variance input is now also refused explicitly and up front, with a clear
reason, rather than relying on sklearn happening to produce a NaN. A
genuinely unexpected backend exception is recorded as a `failed` result
first, then re-raised -- never silently swallowed into a soft failure that
could go unnoticed.
"""
from __future__ import annotations

import hashlib
import json
import math
import os
import platform
import sys
import time
from collections import Counter
from dataclasses import dataclass, field, replace

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from artifact_schema import find_nonfinite_floats, resolve_source_revision, sha256_file, sha256_text  # noqa: E402
from dataset_schema import validate_dataset_config  # noqa: E402
import cli as _cli  # noqa: E402

EMBEDDING_RESULT_SCHEMA = "scfm-embedding-result/1"

#: Per `review:20261013_135634`: this module's own source-identity
#: binding -- every module whose code path actually runs on every call.
SOURCE_MODULES = ("embedding_providers.py", "dataset_schema.py", "cli.py", "artifact_schema.py")

#: Explicit this-round scope bound: a dense in-memory read of a counts_slot
#: larger than this many elements (float64: this bound is ~160MB) is
#: refused outright rather than silently attempted -- per DEVELOPMENT.md's
#: standing caution against reading a large real matrix in a context that
#: can OOM. A future round that wants more must raise this deliberately
#: (and likely add a chunked/sparse-native reader), not relax a check it
#: does not also replace.
_MAX_DENSE_ELEMENTS = 20_000_000

_VALID_STATUSES = ("complete", "failed", "not_run")


def _reject_nonfinite_constant(token: str):
    raise ValueError(f"non-finite JSON constant {token!r} is not accepted")


def _dependency_versions() -> dict:
    import scipy
    import sklearn
    import h5py
    return {
        "sklearn": sklearn.__version__, "numpy": np.__version__,
        "scipy": scipy.__version__, "h5py": h5py.__version__,
        "python": platform.python_version(),
    }


def _counts_block_sha256(array: np.ndarray, cell_ids, feature_ids, counts_slot: str) -> str:
    """Binds the EXACT numeric block actually used to the cell/feature IDs
    and slot that produced it -- so a frozen re-run can tell a genuinely
    identical input block apart from one that merely has the same shape.

    Per `review:20261013_210131`'s non-blocking note (corrected
    per `review:20261013_222450` A.5): `fit_transform` calls this
    on `X` AFTER `log1p` has already been applied (when `log1p=True`) but
    BEFORE `standardize`'s mean/std scaling -- "pre-preprocessing" is
    therefore inaccurate; the exact phase is "post cell/feature selection,
    post optional log1p, pre standardize".
    """
    meta = json.dumps({
        "counts_slot": counts_slot, "cell_ids": list(cell_ids),
        "feature_ids": list(feature_ids), "shape": list(array.shape),
        "dtype": str(array.dtype),
    }, sort_keys=True)
    array_sha = hashlib.sha256(np.ascontiguousarray(array).tobytes()).hexdigest()
    return sha256_text(meta + array_sha)


def _read_var_index_values(path: str, var_group_path: str) -> list:
    """Reads the gene-ID VALUES of a counts_slot's own gene axis (`var` or
    `raw/var`, as resolved by `cli.gene_axis_group_for_slot`) -- mirrors
    `cli.read_obs_column_values`'s plain-dataset-vs-categorical-group
    handling, applied to the axis's own `_index` rather than a named obs
    column (kept as a separate, narrowly-scoped function rather than
    reusing `cli.read_obs_column_values` directly, since that function's
    contract is specifically "one named obs COLUMN", not "a group's own
    index").
    """
    import h5py
    with h5py.File(path, "r") as f:
        if var_group_path not in f:
            raise SystemExit(f"{var_group_path} not found in {path}")
        group = f[var_group_path]
        index_name = group.attrs.get("_index")
        if index_name is None or index_name not in group:
            raise SystemExit(
                f"{var_group_path} has no resolvable '_index' attribute in {path}")
        node = group[index_name]
        if isinstance(node, h5py.Group):
            codes = node["codes"][()]
            categories = node["categories"][()]
            out = []
            for c in codes:
                if c < 0:
                    raise SystemExit(
                        f"{var_group_path}'s index has an unset category code -- "
                        "every gene must have a real identity")
                v = categories[c]
                out.append(v.decode() if isinstance(v, bytes) else v)
            return out
        values = node[()]
        return [v.decode() if isinstance(v, bytes) else v for v in values]


def _validate_gene_axis_values(gene_ids: list, var_group_path: str) -> list[str]:
    """Per `review:20261013_210131` F2: validated BEFORE any
    ID->position map, `feature_ids_override`, or `hvg_top_n` selection is
    built from `gene_ids` -- a duplicated/empty/non-string gene ID must
    become a clean `failed` record here, never a silent
    `{g: i for i, g in enumerate(gene_ids)}` collapse onto the LAST
    occurrence of a repeated ID (which the prior round's `gpos` map did,
    picking a numerically-wrong column with no indication anything was
    dropped).
    """
    if not gene_ids:
        return [f"{var_group_path} has an empty gene axis -- nothing to "
               "select features from"]
    bad_idx = [i for i, g in enumerate(gene_ids) if not isinstance(g, str) or not g]
    if bad_idx:
        return [f"{var_group_path} has {len(bad_idx)} non-string/empty gene ID "
               f"value(s) at index/indices {bad_idx[:5]}"
               + (" (truncated)" if len(bad_idx) > 5 else "")
               + " -- every gene must have a real string identity"]
    counts = Counter(gene_ids)
    dupes = sorted(g for g, n in counts.items() if n > 1)
    if dupes:
        return [f"{var_group_path} has {len(dupes)} duplicate gene ID value(s): "
               f"{dupes[:5]}" + (" (truncated)" if len(dupes) > 5 else "")
               + " -- this module never silently merges, suffixes, or picks a "
               "single column for a duplicated gene ID; if gene-symbol "
               "deduplication is genuinely needed, it requires a separate, "
               "explicit, traceable mapping step, not a default taken here"]
    return []


def _zero_total_variance_reason(X: np.ndarray) -> str | None:
    """Per `review:20261013_210131` F3: an explicit, up-front
    refusal for an input with nothing for PCA to explain (e.g. an all-zero
    counts block), rather than relying on sklearn happening to produce a
    NaN `explained_variance_ratio_` that then has to be caught later.
    """
    if float(np.var(X, axis=0).sum()) == 0.0:
        return ("the selected counts block has zero total variance across "
               "all selected genes -- a PCA has nothing to explain (every "
               "value is constant); refused explicitly rather than letting "
               "the backend produce a NaN explained_variance_ratio")
    return None


def _pca_output_finiteness_problem(embedding, components, mean,
                                   explained_variance_ratio, singular_values) -> str | None:
    """Per `review:20261013_210131` F3: checks EVERY array this
    module is about to persist (not just embedding/components/mean, which
    let a NaN `explained_variance_ratio_` reach `_write_embedding_result`
    and write `embedding_matrices.npz` before the JSON manifest's own
    validation caught it, leaving an orphaned npz with no
    EMBEDDING_RESULT.json beside it).
    """
    if not (np.all(np.isfinite(embedding)) and np.all(np.isfinite(components))
           and np.all(np.isfinite(mean)) and np.all(np.isfinite(explained_variance_ratio))
           and np.all(np.isfinite(singular_values))):
        return "PCA produced non-finite value(s) -- refusing to persist it"
    return None


def _read_counts_slot_dense(path: str, counts_slot: str, row_index: list) -> np.ndarray:
    """Reads the declared counts_slot (dense `h5py.Dataset` OR an anndata
    `csr_matrix`/`csc_matrix` group) as a dense float64 array, selects
    `row_index` (cells, in the given order), and returns ALL genes -- gene
    subsetting happens afterward, in-memory (see `fit_transform`). An
    unrecognized sparse encoding, a non-2-D slot, or a block over
    `_MAX_DENSE_ELEMENTS` is refused outright, never guessed at or silently
    truncated.

    Per `review:20261013_210131`'s non-blocking note (corrected
    per `review:20261013_222450` A.5): the `_MAX_DENSE_ELEMENTS`
    check below is computed from the counts_slot's FULL `n_obs x n_vars`
    shape, and `node[()]`/`mat.todense()` reads that FULL slot into memory,
    BEFORE `row_index` ever selects a cell subset -- a `cell_ids_subset`
    narrower than the full counts_slot does NOT reduce this check or this
    read. The refusal message below no longer suggests `cell_ids_subset` as
    a way around the bound, since it is not one.
    """
    import h5py
    with h5py.File(path, "r") as f:
        if counts_slot == "X":
            node = f["X"]
        elif counts_slot == "raw.X":
            node = f["raw"]["X"]
        elif counts_slot.startswith("layers/"):
            node = f["layers"][counts_slot[len("layers/"):]]
        else:
            raise SystemExit(f"unrecognized counts_slot: {counts_slot!r}")

        if isinstance(node, h5py.Dataset):
            shape = node.shape
            if len(shape) != 2:
                raise SystemExit(f"counts_slot={counts_slot!r} has shape {shape}, not 2-D")
            n_obs, n_vars = shape
            if n_obs * n_vars > _MAX_DENSE_ELEMENTS:
                raise SystemExit(
                    f"counts_slot={counts_slot!r} has {n_obs}x{n_vars}="
                    f"{n_obs * n_vars} elements, over this round's explicit "
                    f"{_MAX_DENSE_ELEMENTS}-element dense-read bound -- this "
                    "check is against the FULL counts_slot shape; a narrower "
                    "cell_ids_subset does NOT reduce it (the full slot is read "
                    "before any row selection). Out of this round's scope (a "
                    "chunked/sparse-native reader is not implemented); never "
                    "silently truncated here")
            dense = node[()].astype(np.float64)
        else:
            encoding = node.attrs.get("encoding-type", b"")
            if isinstance(encoding, bytes):
                encoding = encoding.decode()
            if encoding not in ("csr_matrix", "csc_matrix"):
                raise SystemExit(
                    f"counts_slot={counts_slot!r} is a group with "
                    f"encoding-type={encoding!r}, which this round's reader "
                    "does not recognize (only dense datasets and csr_matrix/"
                    "csc_matrix groups are supported) -- refusing to guess its "
                    "structure")
            shape = [int(v) for v in node.attrs["shape"]]
            n_obs, n_vars = shape
            if n_obs * n_vars > _MAX_DENSE_ELEMENTS:
                raise SystemExit(
                    f"counts_slot={counts_slot!r} has {n_obs}x{n_vars}="
                    f"{n_obs * n_vars} elements, over this round's explicit "
                    f"{_MAX_DENSE_ELEMENTS}-element dense-read bound -- this "
                    "check is against the FULL counts_slot shape; a narrower "
                    "cell_ids_subset does NOT reduce it")
            import scipy.sparse as sp
            data = node["data"][()]
            indices = node["indices"][()]
            indptr = node["indptr"][()]
            cls = sp.csr_matrix if encoding == "csr_matrix" else sp.csc_matrix
            mat = cls((data, indices, indptr), shape=(n_obs, n_vars))
            dense = np.asarray(mat.todense(), dtype=np.float64)

    return dense[row_index, :]


@dataclass(frozen=True)
class EmbeddingResult:
    #: "complete": a PCA was actually fit and transformed. "failed": a
    #: structural problem (invalid config, bad slot, non-finite input,
    #: unestimable n_components, zero-variance gene under standardize=True)
    #: stopped the run. "not_run": never attempted (the caller's concern).
    status: str
    reason: str = ""
    embedding_path: str = ""
    embedding_sha256: str = ""
    cell_ids: tuple = ()
    n_components: int = 0
    #: The genes actually used as PCA INPUT features (post selection, pre
    #: any dimensionality reduction) -- never the full gene axis unless that
    #: is genuinely what was used.
    feature_ids: tuple = ()
    explained_variance_ratio: tuple = ()
    #: Binds the EXACT numeric input block actually used -- post cell/
    #: feature selection, post optional log1p, PRE standardize (see
    #: `_counts_block_sha256`'s docstring for why "pre-preprocessing" is
    #: not accurate when log1p=True). When standardize=True, note that the
    #: fit mean/std used to scale X are NOT themselves saved anywhere in
    #: this result or the persisted npz -- this hash lets a re-run confirm
    #: it started from the same raw block, but cannot be used to replay the
    #: exact standardize transform on NEW data without recomputing mean/std
    #: from this same block again.
    input_counts_sha256: str = ""
    backend: str = "unknown"
    backend_version: str = "unknown"
    #: n_components, log1p, standardize, random_seed, svd_solver,
    #: n_cells_used, n_genes_used, and (optional) hvg_top_n/
    #: feature_ids_override_provided/cell_ids_subset_provided.
    params: dict = field(default_factory=dict)
    #: dataset_id, species, gene_id_namespace, counts_slot, counts_semantics,
    #: the resolved matrix path, and the dataset config's own content hash.
    dataset_identity: dict = field(default_factory=dict)
    source_code_revision_or_hash: str = ""
    dependency_versions: dict = field(default_factory=dict)

    def to_dict(self) -> dict:
        d = dict(self.__dict__)
        d["cell_ids"] = list(self.cell_ids)
        d["feature_ids"] = list(self.feature_ids)
        d["explained_variance_ratio"] = list(self.explained_variance_ratio)
        return d


#: Required `params` keys for a 'complete' result, with their expected kind
#: ("int" excludes bool explicitly, matching gsea_runner.py's convention).
_REQUIRED_EMBEDDING_PARAMS = {
    "n_components": "int", "log1p": "bool", "standardize": "bool",
    "random_seed": "int", "svd_solver": "str", "n_cells_used": "int",
    "n_genes_used": "int",
}


def _validate_complete_embedding_params(params) -> list[str]:
    if not isinstance(params, dict) or not params:
        return ["status='complete' requires a non-empty 'params' dict"]
    problems = []
    for key, kind in _REQUIRED_EMBEDDING_PARAMS.items():
        if key not in params:
            problems.append(f"status='complete' requires params[{key!r}]")
            continue
        v = params[key]
        ok = ((kind == "bool" and isinstance(v, bool))
              or (kind == "int" and isinstance(v, int) and not isinstance(v, bool))
              or (kind == "str" and isinstance(v, str)))
        if not ok:
            problems.append(f"params[{key!r}]={v!r} must be a {kind}, got {type(v).__name__}")
    n_components = params.get("n_components")
    if isinstance(n_components, int) and not isinstance(n_components, bool) and n_components < 1:
        problems.append(f"params['n_components']={n_components!r} must be >= 1")
    return problems


def validate_embedding_result(payload: dict) -> list[str]:
    """Shared by the writer and `read_embedding_result`, matching
    `gsea_runner.validate_gsea_result`'s rigor: a tampered or merely
    incomplete on-disk file is rejected at BOTH boundaries.
    """
    problems = []
    if not isinstance(payload, dict):
        return [f"embedding result must be a dict, got {type(payload).__name__}"]
    if payload.get("schema") != EMBEDDING_RESULT_SCHEMA:
        problems.append(f"schema {payload.get('schema')!r} != expected {EMBEDDING_RESULT_SCHEMA!r}")
    status = payload.get("status")
    if status not in _VALID_STATUSES:
        problems.append(f"status {status!r} not in {_VALID_STATUSES}")
    problems.extend(find_nonfinite_floats(payload))

    if status != "complete" and not payload.get("reason"):
        problems.append(f"status={status!r} requires a non-empty 'reason'")

    if status == "complete":
        for key in ("embedding_path", "embedding_sha256", "backend", "backend_version",
                    "input_counts_sha256", "source_code_revision_or_hash"):
            v = payload.get(key)
            if not isinstance(v, str) or not v:
                problems.append(
                    f"status='complete' requires a non-empty string {key!r}, got {v!r}")

        cell_ids = payload.get("cell_ids")
        if not isinstance(cell_ids, list) or not cell_ids \
                or not all(isinstance(c, str) and c for c in cell_ids):
            problems.append(
                "status='complete' requires a non-empty list of non-empty string cell_ids")
        elif len(set(cell_ids)) != len(cell_ids):
            problems.append("cell_ids contains duplicate value(s)")

        feature_ids = payload.get("feature_ids")
        if not isinstance(feature_ids, list) or not feature_ids \
                or not all(isinstance(g, str) and g for g in feature_ids):
            problems.append(
                "status='complete' requires a non-empty list of non-empty string feature_ids")
        elif len(set(feature_ids)) != len(feature_ids):
            problems.append("feature_ids contains duplicate value(s)")

        n_components = payload.get("n_components")
        if not isinstance(n_components, int) or isinstance(n_components, bool) or n_components < 1:
            problems.append(f"n_components must be a positive int, got {n_components!r}")

        evr = payload.get("explained_variance_ratio")
        if not isinstance(evr, list) or not evr:
            problems.append(
                "status='complete' requires a non-empty 'explained_variance_ratio' list")
        elif isinstance(n_components, int) and not isinstance(n_components, bool) \
                and len(evr) != n_components:
            problems.append(
                f"explained_variance_ratio has {len(evr)} entries, expected "
                f"n_components={n_components}")
        else:
            for i, v in enumerate(evr):
                if not isinstance(v, (int, float)) or isinstance(v, bool) or not math.isfinite(v):
                    problems.append(
                        f"explained_variance_ratio[{i}]={v!r} must be a finite number")

        problems.extend(_validate_complete_embedding_params(payload.get("params")))

        dataset_identity = payload.get("dataset_identity")
        if not isinstance(dataset_identity, dict) or not dataset_identity:
            problems.append("status='complete' requires a non-empty 'dataset_identity' dict")

        dep = payload.get("dependency_versions")
        if not isinstance(dep, dict) or not dep:
            problems.append("status='complete' requires a non-empty 'dependency_versions' dict")
        elif not all(isinstance(k, str) and isinstance(v, str) for k, v in dep.items()):
            problems.append(f"dependency_versions keys and values must all be strings, got {dep!r}")
    return problems


def _write_embedding_result(out_dir: str, result: EmbeddingResult, embedding=None,
                            components=None, mean=None, singular_values=None) -> EmbeddingResult:
    if result.status == "complete":
        npz_path = os.path.join(out_dir, "embedding_matrices.npz")
        tmp = npz_path + ".inprogress"
        with open(tmp, "wb") as fh:
            np.savez(
                fh,
                embedding=np.ascontiguousarray(embedding, dtype=np.float64),
                components=np.ascontiguousarray(components, dtype=np.float64),
                mean=np.ascontiguousarray(mean, dtype=np.float64),
                singular_values=np.ascontiguousarray(singular_values, dtype=np.float64))
            fh.flush()
            os.fsync(fh.fileno())
        os.replace(tmp, npz_path)
        result = replace(result, embedding_path=npz_path, embedding_sha256=sha256_file(npz_path))

    payload = result.to_dict()
    payload["schema"] = EMBEDDING_RESULT_SCHEMA
    payload["generated_utc"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    problems = validate_embedding_result(payload)
    if problems:
        raise ValueError(f"refusing to write an invalid embedding result: {problems}")

    result_path = os.path.join(out_dir, "EMBEDDING_RESULT.json")
    tmp2 = result_path + ".inprogress"
    with open(tmp2, "w") as fh:
        json.dump(payload, fh, indent=2, sort_keys=True, allow_nan=False)
        fh.flush()
        os.fsync(fh.fileno())
    os.replace(tmp2, result_path)
    return result


def read_embedding_result(path: str) -> dict:
    if not os.path.exists(path):
        raise SystemExit(f"embedding result not found: {path}")
    with open(path) as fh:
        try:
            payload = json.load(fh, parse_constant=_reject_nonfinite_constant)
        except ValueError as e:
            raise SystemExit(f"embedding result at {path} failed to parse: {e}")
    problems = validate_embedding_result(payload)
    if problems:
        raise SystemExit(f"embedding result at {path} failed validation: {problems}")
    return payload


def load_embedding_matrices(manifest: dict) -> dict:
    """Loads the npz artifact a 'complete' manifest references, re-checks its
    file hash against `embedding_sha256` (a tampered or stale npz is caught
    here, not only at write time), and checks every array's shape against
    the manifest's own `cell_ids`/`feature_ids`/`n_components` -- the
    read-back contract this round's review explicitly asked for.
    """
    if manifest.get("status") != "complete":
        raise SystemExit(
            f"load_embedding_matrices requires status='complete', got "
            f"{manifest.get('status')!r}")
    path = manifest["embedding_path"]
    actual_sha = sha256_file(path)
    if actual_sha != manifest["embedding_sha256"]:
        raise SystemExit(
            f"{path} has file hash {actual_sha!r}, which does not match the "
            f"manifest's recorded embedding_sha256 {manifest['embedding_sha256']!r}")
    with np.load(path) as npz:
        arrays = {k: npz[k] for k in ("embedding", "components", "mean", "singular_values")}

    n_cells = len(manifest["cell_ids"])
    n_components = manifest["n_components"]
    n_features = len(manifest["feature_ids"])
    if arrays["embedding"].shape != (n_cells, n_components):
        raise SystemExit(
            f"embedding array shape {arrays['embedding'].shape} != expected "
            f"({n_cells}, {n_components})")
    if arrays["components"].shape != (n_components, n_features):
        raise SystemExit(
            f"components array shape {arrays['components'].shape} != expected "
            f"({n_components}, {n_features})")
    if arrays["mean"].shape != (n_features,):
        raise SystemExit(f"mean array shape {arrays['mean'].shape} != expected ({n_features},)")
    for name, arr in arrays.items():
        if not np.all(np.isfinite(arr)):
            raise SystemExit(f"{name!r} array in {path} contains non-finite value(s)")
    return arrays


class EmbeddingProvider:
    """Contract: input is a VALIDATED dataset config dict plus the config
    FILE's own path (so a relative `expression_matrix_path` resolves the
    same way `cli.py` already resolves it) -- never a bare matrix file path
    with `counts_slot`/gene-axis/counts-semantics stripped away.
    """

    def fit_transform(self, dataset_config: dict, dataset_config_path: str, *,
                      out_dir: str, **params) -> EmbeddingResult:
        raise NotImplementedError(
            "EmbeddingProvider.fit_transform is a contract stub -- no matrix "
            "is read here. Use PCAEmbeddingProvider for the real local "
            "implementation.")


class PCAEmbeddingProvider(EmbeddingProvider):
    """Real, scikit-learn-backed PCA baseline. See module docstring for the
    scope this round explicitly supports and explicitly refuses.
    """

    def fit_transform(self, dataset_config: dict, dataset_config_path: str, *,
                      out_dir: str, n_components: int = 10, log1p: bool = False,
                      standardize: bool = False, hvg_top_n: int | None = None,
                      feature_ids_override: tuple | None = None,
                      cell_ids_subset: tuple | None = None, random_seed: int = 0,
                      svd_solver: str = "auto") -> EmbeddingResult:
        if os.path.exists(out_dir) and os.listdir(out_dir):
            raise SystemExit(f"{out_dir} must be new or empty")
        os.makedirs(out_dir, exist_ok=True)

        import sklearn

        source_revision = resolve_source_revision(HERE, SOURCE_MODULES)
        dep_versions = _dependency_versions()
        common = dict(
            backend="sklearn.decomposition.PCA", backend_version=sklearn.__version__,
            source_code_revision_or_hash=source_revision, dependency_versions=dep_versions)

        def finish(result: EmbeddingResult, **arrays) -> EmbeddingResult:
            return _write_embedding_result(out_dir, result, **arrays)

        def fail(reason: str) -> EmbeddingResult:
            return finish(EmbeddingResult(status="failed", reason=reason, **common))

        if hvg_top_n is not None and feature_ids_override is not None:
            return fail("pass at most one of hvg_top_n / feature_ids_override, "
                       "never both -- feature selection must be unambiguous")
        if n_components < 1:
            return fail(f"n_components={n_components} must be >= 1")

        val = validate_dataset_config(dataset_config)
        if not val.ok:
            return fail(f"dataset_config failed validation: {[str(i) for i in val.issues]}")

        matrix_path = _cli.resolve_matrix_path(dataset_config_path, dataset_config)
        if not matrix_path:
            return fail("dataset_config has no expression_matrix_path -- nothing to read")

        structural_problems = _cli.validate_matrix_structure(matrix_path, dataset_config)
        if structural_problems:
            return fail(f"matrix structure invalid: {structural_problems}")

        counts_slot = dataset_config["counts_slot"]
        cell_id_col = dataset_config["cell_id_column"]
        all_cell_ids = _cli.read_obs_column_values(matrix_path, cell_id_col)

        if cell_ids_subset is None:
            selected_cell_ids = list(all_cell_ids)
            row_index = list(range(len(all_cell_ids)))
        else:
            pos = {c: i for i, c in enumerate(all_cell_ids)}
            missing = [c for c in cell_ids_subset if c not in pos]
            if missing:
                return fail(
                    f"cell_ids_subset contains {len(missing)} ID(s) not present "
                    f"in {cell_id_col!r}: {missing[:5]}")
            selected_cell_ids = list(cell_ids_subset)
            row_index = [pos[c] for c in selected_cell_ids]
        if len(set(selected_cell_ids)) != len(selected_cell_ids):
            return fail("cell_ids_subset contains duplicate value(s)")
        if not selected_cell_ids:
            return fail("no cells selected -- nothing to fit a PCA on")

        var_group_path = _cli.gene_axis_group_for_slot(counts_slot)
        # Per the GSEA chain's own R2 lesson: a structural problem this
        # module's own readers detect (oversized block, unrecognized sparse
        # encoding, non-2-D slot -- none of which validate_matrix_structure
        # above checks) must become a clean, audited 'failed' record, never
        # an uncaught SystemExit that skips writing EMBEDDING_RESULT.json at
        # all.
        try:
            gene_ids = _read_var_index_values(matrix_path, var_group_path)
        except SystemExit as e:
            return fail(str(e))

        # Per F2: validated BEFORE any ID->position map, feature_ids_override,
        # or hvg_top_n selection is built from gene_ids -- and before the
        # (potentially large) dense counts read below.
        gene_axis_problems = _validate_gene_axis_values(gene_ids, var_group_path)
        if gene_axis_problems:
            return fail("gene axis invalid: " + "; ".join(gene_axis_problems))

        try:
            raw = _read_counts_slot_dense(matrix_path, counts_slot, row_index)
        except SystemExit as e:
            return fail(str(e))

        bad = np.where(~np.isfinite(raw))
        if bad[0].size:
            return fail(
                f"counts_slot={counts_slot!r} contains {bad[0].size} non-finite "
                "value(s) in the selected cells -- refusing to fit a PCA on them")

        if feature_ids_override is not None:
            gpos = {g: i for i, g in enumerate(gene_ids)}
            missing = [g for g in feature_ids_override if g not in gpos]
            if missing:
                return fail(
                    f"feature_ids_override contains {len(missing)} gene ID(s) not "
                    f"present in the counts_slot's gene axis: {missing[:5]}")
            if len(set(feature_ids_override)) != len(feature_ids_override):
                return fail("feature_ids_override contains duplicate gene ID(s)")
            feat_idx = [gpos[g] for g in feature_ids_override]
            selected_feature_ids = list(feature_ids_override)
        elif hvg_top_n is not None:
            if hvg_top_n < 1 or hvg_top_n > len(gene_ids):
                return fail(f"hvg_top_n={hvg_top_n} out of range for {len(gene_ids)} genes")
            pre = np.log1p(raw) if log1p else raw
            variances = pre.var(axis=0)
            order = np.argsort(-variances, kind="stable")[:hvg_top_n]
            feat_idx = sorted(order.tolist())
            selected_feature_ids = [gene_ids[i] for i in feat_idx]
        else:
            feat_idx = list(range(len(gene_ids)))
            selected_feature_ids = list(gene_ids)

        X = raw[:, feat_idx]
        if log1p:
            if np.any(X < 0):
                return fail(
                    "log1p=True but the selected block contains negative "
                    "value(s) -- log1p requires non-negative counts; this is "
                    "not silently clipped or coerced")
            X = np.log1p(X)

        input_sha = _counts_block_sha256(X, selected_cell_ids, selected_feature_ids, counts_slot)

        if standardize:
            std = X.std(axis=0)
            zero_var = np.where(std == 0)[0]
            if zero_var.size:
                return fail(
                    f"{zero_var.size} selected gene(s) have zero variance across "
                    "the selected cells -- standardize=True cannot divide by "
                    "zero; exclude them explicitly (feature_ids_override) or "
                    "set standardize=False")
            X = (X - X.mean(axis=0)) / std

        n_samples, n_features = X.shape
        if n_components > min(n_samples, n_features):
            return fail(
                f"n_components={n_components} exceeds min(n_samples={n_samples}, "
                f"n_features={n_features})")

        zero_variance_reason = _zero_total_variance_reason(X)
        if zero_variance_reason is not None:
            return fail(zero_variance_reason)

        from sklearn.decomposition import PCA
        pca = PCA(n_components=n_components, svd_solver=svd_solver, random_state=random_seed)
        try:
            embedding = pca.fit_transform(X)
        except Exception as exc:
            # Per F3: an unexpected backend error is recorded as a clean
            # 'failed' result FIRST (so the audit trail has it), then
            # re-raised -- never silently swallowed into a soft failure
            # that could go unnoticed.
            fail(f"PCA backend raised an unexpected error: {type(exc).__name__}: {exc}")
            raise

        finiteness_problem = _pca_output_finiteness_problem(
            embedding, pca.components_, pca.mean_, pca.explained_variance_ratio_,
            pca.singular_values_)
        if finiteness_problem is not None:
            return fail(finiteness_problem)

        params = dict(
            n_components=n_components, log1p=log1p, standardize=standardize,
            hvg_top_n=hvg_top_n, feature_ids_override_provided=feature_ids_override is not None,
            random_seed=random_seed, svd_solver=svd_solver, n_cells_used=n_samples,
            n_genes_used=n_features, cell_ids_subset_provided=cell_ids_subset is not None)

        dataset_identity = dict(
            dataset_id=dataset_config["dataset_id"], species=dataset_config["species"],
            gene_id_namespace=dataset_config["gene_id_namespace"], counts_slot=counts_slot,
            counts_semantics=dataset_config["counts_semantics"],
            expression_matrix_path_resolved=matrix_path,
            dataset_config_sha256=sha256_text(json.dumps(dataset_config, sort_keys=True)))

        result = EmbeddingResult(
            status="complete", cell_ids=tuple(selected_cell_ids), n_components=n_components,
            feature_ids=tuple(selected_feature_ids),
            explained_variance_ratio=tuple(float(v) for v in pca.explained_variance_ratio_),
            input_counts_sha256=input_sha, params=params, dataset_identity=dataset_identity,
            **common)

        return finish(result, embedding=embedding, components=pca.components_,
                     mean=pca.mean_, singular_values=pca.singular_values_)

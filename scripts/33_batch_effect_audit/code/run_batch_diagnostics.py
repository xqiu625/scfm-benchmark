#!/usr/bin/env python3
"""Stage-33 batch-effect diagnostic CLI. v2 (2026-10-09, second review round).

Expression/PCA mode works without any FM embedding -- optional embeddings
(obsm keys, e.g. a real FM or the synthetic `X_synthFM` fixture) are
analyzed in addition when present, never required UNLESS explicitly
requested via --embeddings, in which case a missing/failed one now fails
the run (see "Status propagation" below). Analyzes each requested lineage
SEPARATELY (never pools lineages, which would make lineage identity itself
look like "batch structure"). Uses the shared, tested `scib_batch_asw`
(scripts/00_shared/code/batch_asw.py) -- not a local copy.

What this does NOT do: it does not decide whether a batch effect is
"present" in some single-number sense, does not choose a correction method,
and does not touch the real benchmark population -- the patient-balanced
sample it draws is a DIAGNOSTIC subsample, reported with its exact size and
rule, never silently substituted for the full cohort.

v2 fixes defects a 2026-10-09 review found in v1
(review:20261009_133450), confirmed independently before fixing:

- Seeding used Python's built-in `hash()`, which is per-process-randomized
  for str/tuple-of-str (PYTHONHASHSEED) -- a fixed `--seed` did NOT actually
  reproduce the same patient-balanced selection across separate process
  invocations (confirmed: `grep hash((seed` found the exact pattern).
  Replaced with `stable_seed.derive_seed()` (scripts/00_shared/code/), which
  is proven stable across PYTHONHASHSEED values by a committed test.
- UMAP coordinate CSVs carried patient_id but no cell_id -- patient_id
  repeats across cells and cannot identify a row. Fixed: cell_id is now the
  first column.
- A Harmony failure/block, or a missing/unsupported EXPLICITLY REQUESTED
  embedding, could leave the enclosing embedding `status="complete"` and the
  whole run `status="complete"` with exit 0, silently omitting the
  requested correction or comparison. Fixed: requested-but-unavailable or
  requested-but-failed work now sets that embedding's status to `failed`,
  which now DOES propagate to the top-level status and a nonzero exit (see
  `REQUIRED_WORK_SPEC` below -- explicitly distinguishes required from
  optional so a truly optional miss does not fail the whole run).
- `manifest.json`'s `code_revision` could be `"UNKNOWN (...)"` while
  `status` still read `"complete"` -- a caller had no signal that code
  identity was unresolved. Fixed: `--code-revision` can be passed explicitly
  (the deployer's job, since a loose-file scratch copy is not a git
  checkout); if omitted and git detection fails, `code_revision_verified`
  is `false` and the run now FAILS rather than silently publishing
  unresolved provenance as if it were fine.
- No no-clobber check -- a rerun into an existing --out-dir could silently
  mix with or overwrite a previous run's files, contrary to this project's
  standing "never overwrite a result" rule. Fixed: refuses if --out-dir
  already contains any file.
- `pd.DataFrame(metrics_rows)` with an empty list has no columns (not a
  stable, loadable empty schema). Fixed: explicit `METRICS_COLS`.
- A cell with outcome=="unknown" (or any caller-declared sentinel) was
  treated as a THIRD biological class by the bio-label-stratified metrics,
  rather than being excluded from biological scoring as explicitly missing.
  Fixed: `--exclude-outcome-values` (default `unknown`) filters those cells
  out of the diagnostic subsample before any outcome-stratified computation,
  with the excluded count reported, not silently absorbed.
- Representations were saved as bare `.npy` with no cell-ID binding.
  Fixed: every saved representation (uncorrected PCA/embedding AND
  Harmony-corrected) is now saved alongside its exact `cell_id` ordering in
  one `.npz`, and reloaded immediately after writing to verify identity.
- (2026-10-09, later same day) UMAP coordinates were only saved for the
  UNCORRECTED representation, making a real before/after visual comparison
  impossible. Now saved for both scopes (`save_umap_coords()`), each
  cell-ID-bound, so `make_batch_effect_plots.py` can render an actual
  before/after PNG, not just the metrics table.
- A failure during per-lineage processing before the final manifest write
  left no persisted record at all. Fixed: the manifest is now written
  (status="failed", with whatever partial content exists) in an
  exception handler that then re-raises, not only at a successful end.

v3 fixes defects a 2026-10-09 review found (review:20261009_230728),
confirmed independently before fixing:

- `--population-stage` free text could legally drift from the documented job
  spec's own wording for the same stage (confirmed: it had), and a mismatch
  only WARNED. Added `--stage-id` (scripts/00_shared/code/population_stage.py):
  a small fixed vocabulary, matched exactly against the input's own
  `uns["population_stage_id"]`, HARD FAILS on any mismatch.
- `code_revision_verified` was literally `code_revision is not None` -- ANY
  string passed "verified". Replaced with `_resolve_code_identity()`: verified
  now means a real git checkout, HEAD resolved, working tree clean, and
  matching the asserted value. A loose-file scratch deployment (no `.git`,
  used throughout this project's HPCC work) can never be "verified" this way;
  `source_file_sha256` is recorded instead as the honest substitute -- an
  immutable record of what ran, not a verified identity.
- The --out-dir non-empty check was check-then-act, not atomic -- confirmed
  as the root cause of a real race this session (an orphaned process, not
  actually dead, raced a later "clean" retry). Added
  `scripts/00_shared/code/run_ownership.py::claim_ownership()`: an OS-atomic
  exclusive create, so the loser gets FileExistsError from the syscall, never
  from a race window.
- Unknown-outcome cells were excluded from the ENTIRE diagnostic population
  (PCA, Harmony, neighbor-mixing) before sampling, coupling technical-
  integration eligibility to an unrelated biological label's availability.
  Fixed: sampling and neighbor_batch_mixing now use every lineage-matched
  cell; only batch_asw and bio_conservation (which genuinely require a real
  bio label) are restricted to `bio_known_mask`, with both cell counts
  reported separately.
- A failed corrected-UMAP did not downgrade `status`, unlike the uncorrected
  analogue -- a requested before/after plot pair could be silently missing
  its corrected half in a run still reporting "complete". Fixed: symmetric.
- The manifest dict (including `_file_sha256(args.h5ad)`) was built BEFORE
  the try/except that writes failure manifests -- a hashing failure left no
  persisted record at all. Fixed: manifest is now written as "running"
  immediately after claiming ownership, before any fallible I/O.

Usage:
    python run_batch_diagnostics.py \\
        --h5ad <path.h5ad> --out-dir <new, empty dir> \\
        --population-stage "synthetic_fixture" --stage-id synthetic_fixture \\
        [--code-revision <git-sha-or-identifier>] \\
        [--embeddings X_synthFM] [--run-harmony]
"""
from __future__ import annotations

import argparse
import json
import logging
import subprocess
import sys
import traceback
from datetime import datetime, timezone
from pathlib import Path

import anndata
import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
from sklearn.neighbors import NearestNeighbors

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "00_shared/code"))
from batch_asw import scib_batch_asw  # noqa: E402
from stable_seed import derive_seed  # noqa: E402
from population_stage import check_stage_id, KNOWN_STAGE_IDS  # noqa: E402
from run_ownership import claim_ownership, OwnershipError  # noqa: E402

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger("run_batch_diagnostics")

METRICS_COLS = ["lineage", "embedding", "scope", "metric", "stratum", "value"]

# Source files actually imported/executed by this run -- hashed into the manifest
# regardless of whether git identity can be verified (see _resolve_code_identity),
# so the EXACT content that ran is captured immutably even on a loose-file scratch
# deployment that is not a git checkout at all.
_SOURCE_FILES_RELATIVE_TO_REPO_ROOT = [
    "scripts/33_batch_effect_audit/code/run_batch_diagnostics.py",
    "scripts/00_shared/code/batch_asw.py",
    "scripts/00_shared/code/stable_seed.py",
    "scripts/00_shared/code/population_stage.py",
    "scripts/00_shared/code/run_ownership.py",
]


def _file_sha256(path: Path) -> str:
    import hashlib
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            chunk = f.read(1 << 20)
            if not chunk:
                break
            h.update(chunk)
    return h.hexdigest()


def _git_commit(repo_dir: Path) -> str | None:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=repo_dir, text=True, stderr=subprocess.DEVNULL
        ).strip()
    except Exception:
        return None


def _git_clean(repo_dir: Path) -> bool | None:
    """True if `git status --porcelain` is empty (no uncommitted changes), False if
    not, None if this is not a git checkout at all (status call itself failed)."""
    try:
        out = subprocess.check_output(
            ["git", "status", "--porcelain"], cwd=repo_dir, text=True, stderr=subprocess.DEVNULL
        )
        return len(out.strip()) == 0
    except Exception:
        return None


def _dependency_versions() -> dict:
    import platform
    versions = dict(python=platform.python_version())
    for mod_name in ("anndata", "numpy", "pandas", "scanpy", "sklearn"):
        try:
            mod = __import__(mod_name)
            versions[mod_name] = getattr(mod, "__version__", "unknown")
        except Exception:
            versions[mod_name] = "not importable"
    try:
        import harmonypy
        versions["harmonypy"] = getattr(harmonypy, "__version__", "unknown")
    except Exception:
        versions["harmonypy"] = "not importable"
    return versions


def _resolve_code_identity(repo_root: Path, asserted_revision: str | None) -> dict:
    """Distinguishes ASSERTED code identity (whatever the caller typed) from
    VERIFIED code identity (an actual, independently-checked fact).

    A 2026-10-09 review found the prior version set `code_revision_verified =
    (code_revision is not None)` -- so ANY string, including a hand-typed label like
    "e93cc6a-plus-plots" that does not correspond to any real git commit, was
    reported as "verified". Fixed: `code_revision_verified` is now True only if this
    IS a real git checkout (has `.git`), `git rev-parse HEAD` resolves, the working
    tree is clean (`git status --porcelain` empty), AND it matches the asserted
    revision if one was supplied. A loose-file scratch deployment (no `.git`, as
    used throughout this project's HPCC scratch work) can never be "verified" this
    way -- `source_file_sha256` is always recorded instead, as the honest substitute:
    not a git-verified identity, but an immutable record of exactly what ran, that
    CAN be cross-checked later against `git show <sha>:<path> | sha256sum`."""
    resolved_head = _git_commit(repo_root)
    clean = _git_clean(repo_root) if resolved_head is not None else None
    if resolved_head is not None and clean:
        verified = (asserted_revision is None) or (asserted_revision == resolved_head)
        method = "git_head_clean_tree" if verified else "asserted_revision_mismatches_git_head"
    elif resolved_head is not None and clean is False:
        verified = False
        method = "git_checkout_has_uncommitted_changes"
    else:
        verified = False
        method = "not_a_git_checkout_unverified_asserted_string"

    source_hashes = {}
    for rel in _SOURCE_FILES_RELATIVE_TO_REPO_ROOT:
        f = repo_root / rel
        source_hashes[rel] = _file_sha256(f) if f.exists() else None

    return dict(
        code_revision_asserted=asserted_revision,
        code_revision_resolved_git_head=resolved_head,
        code_revision_git_tree_clean=clean,
        code_revision_verified=verified,
        code_revision_verification_method=method,
        source_file_sha256=source_hashes,
        dependency_versions=_dependency_versions(),
    )


def patient_balanced_sample(obs: pd.DataFrame, patient_col: str, n_per_patient: int,
                             seed: int) -> np.ndarray:
    """Deterministic: one independent RandomState per patient, keyed by
    (seed, patient_id) via `stable_seed.derive_seed` (NOT Python's `hash()`,
    which is per-process-randomized -- see module docstring) -- so the
    subsample for a given patient does not depend on which other patients
    are present, on iteration order, or on the interpreter's hash seed."""
    chosen = []
    for patient_id, idx in obs.groupby(patient_col).groups.items():
        idx = np.asarray(idx)
        n_take = min(len(idx), n_per_patient)
        rng = np.random.RandomState(derive_seed(seed, "patient", patient_id))
        chosen.extend(rng.choice(idx, size=n_take, replace=False))
    return np.array(sorted(chosen, key=str))


def neighborhood_batch_mixing(X: np.ndarray, batch_labels: np.ndarray, k: int = 15) -> dict:
    """Fraction of each cell's k nearest neighbors (excluding itself) that
    carry a DIFFERENT batch label, averaged over cells. 0 = every neighbor
    shares this cell's batch (no mixing); higher = more mixing. Also reports
    the expected value under random mixing (neighbors' batch label
    independent of this cell's) given the actual batch proportions, as a
    reference point -- a raw mixing fraction means little without it."""
    n = X.shape[0]
    k_eff = min(k, n - 1)
    if k_eff < 1:
        return dict(mean_neighbor_batch_mixing=float("nan"),
                    expected_under_random_mixing=float("nan"), k_used=0,
                    reason="n_cells < 2, no neighbors possible")
    nn = NearestNeighbors(n_neighbors=k_eff + 1).fit(X)
    _, idx = nn.kneighbors(X)
    idx = idx[:, 1:]  # drop self
    own_batch = batch_labels[:, None]
    neighbor_batch = batch_labels[idx]
    frac_diff = (neighbor_batch != own_batch).mean(axis=1)
    batch_props = pd.Series(batch_labels).value_counts(normalize=True)
    expected_random = 1.0 - float((batch_props ** 2).sum())
    return dict(
        mean_neighbor_batch_mixing=float(frac_diff.mean()),
        expected_under_random_mixing=expected_random,
        k_used=k_eff,
        reason="",
    )


def bio_conservation(X: np.ndarray, bio_labels: np.ndarray, seed: int) -> dict:
    """ARI/NMI of KMeans(n_clusters=n_bio_classes) against the bio label --
    a lighter-weight bio-conservation proxy than Leiden+ARI (no neighbors-
    graph resolution to tune), explicit about that simplification. Caller
    must have already excluded any "missing outcome" sentinel (e.g.
    "unknown") from bio_labels -- this function has no policy for that and
    would otherwise treat it as a real third class."""
    classes = sorted(set(bio_labels))
    if len(classes) < 2 or X.shape[0] < len(classes):
        return dict(ari=float("nan"), nmi=float("nan"),
                     reason=f"n_classes={len(classes)}, n_cells={X.shape[0]} -- too few to cluster")
    km = KMeans(n_clusters=len(classes), random_state=seed, n_init=10).fit(X)
    return dict(
        ari=float(adjusted_rand_score(bio_labels, km.labels_)),
        nmi=float(normalized_mutual_info_score(bio_labels, km.labels_)),
        reason="",
    )


def determine_embedding_status(*, asw_skipped_groups: bool, uncorrected_umap_ok: bool,
                                 run_harmony: bool, harmony_blocked_or_failed: bool,
                                 corrected_umap_ok: bool | None) -> str:
    """Pure, directly-testable status-aggregation rule for one embedding --
    extracted from main()'s inline logic in a 2026-10-09 review so the
    corrected-UMAP-failure case (asymmetric with the uncorrected case) could
    be regression-tested deterministically rather than relying on actually
    triggering a flaky real UMAP failure. `corrected_umap_ok` is None when
    --run-harmony was not requested at all (irrelevant then)."""
    if run_harmony and harmony_blocked_or_failed:
        return "failed"  # Harmony was REQUESTED; a block/failure is not optional.
    status = "complete"
    if asw_skipped_groups or not uncorrected_umap_ok:
        status = "partial"
    if run_harmony and corrected_umap_ok is False:
        status = "partial"
    return status


def run_one_embedding(X: np.ndarray, batch_labels: np.ndarray, bio_labels: np.ndarray,
                       bio_known_mask: np.ndarray, seed: int, min_group_size: int) -> dict:
    """Two eligibility populations, kept explicitly separate (2026-10-09 review):
    `neighborhood_batch_mixing` is purely technical (batch labels only) and runs on
    EVERY cell passed in, including ones with an unknown/missing biological outcome
    -- a cell's outcome being unrecorded is not a reason to withhold it from
    technical batch-integration diagnostics. `scib_batch_asw` and
    `bio_conservation` both genuinely require a real biological label (ASW is
    stratified BY it; bio-conservation scores clustering AGAINST it), so both run
    only on `bio_known_mask` -- a cell with an unknown outcome cannot be assigned to
    a real stratification group without inventing one. The two cell counts are
    reported separately so a reader can see exactly which population backs which
    number, not infer it."""
    mixing = neighborhood_batch_mixing(X, batch_labels)
    X_known, batch_known, bio_known = X[bio_known_mask], batch_labels[bio_known_mask], bio_labels[bio_known_mask]
    asw = scib_batch_asw(X_known, batch_known, bio_known, min_group_size=min_group_size)
    bioc = bio_conservation(X_known, bio_known, seed)
    return dict(batch_asw=asw, neighborhood_mixing=mixing, bio_conservation=bioc,
                n_cells_technical=int(X.shape[0]), n_cells_bio_known=int(X_known.shape[0]))


def compute_pca(adata_sub: anndata.AnnData, n_pcs: int, seed: int) -> np.ndarray:
    a = adata_sub.copy()
    a.X = a.layers["counts"].copy() if "counts" in a.layers else a.X
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    n_pcs_eff = min(n_pcs, a.shape[0] - 1, a.shape[1] - 1)
    sc.pp.pca(a, n_comps=n_pcs_eff, random_state=seed)
    return a.obsm["X_pca"]


def compute_umap(X: np.ndarray, seed: int) -> np.ndarray | None:
    try:
        a = anndata.AnnData(X=np.zeros((X.shape[0], 1)))
        a.obsm["X_in"] = X
        n_neighbors = min(15, X.shape[0] - 1)
        if n_neighbors < 2:
            return None
        sc.pp.neighbors(a, n_neighbors=n_neighbors, use_rep="X_in", random_state=seed)
        sc.tl.umap(a, random_state=seed)
        return a.obsm["X_umap"]
    except Exception:
        log.warning("UMAP failed:\n%s", traceback.format_exc())
        return None


def run_harmony_correction(X: np.ndarray, batch_labels: np.ndarray, seed: int) -> np.ndarray:
    """Opt-in, real harmonypy integration -- jointly fit across every
    eligible cell in this diagnostic subsample. This is NOT an out-of-sample
    transform: a cohort not present when this was fit has no defined
    correction here. Callers must label results from this function
    'transductive/exploratory integration', never 'strict held-out'.

    `nclust` is passed explicitly rather than left to harmonypy's own
    default (`round(N/30)`, capped at 100): that heuristic returns 1 for any
    N < 45, and harmonypy 0.0.10 has a real bug at nclust==1 (`sigma` is
    only broadcast from a scalar to a per-cluster array when `nclust > 1`,
    so nclust==1 crashes inside `init_cluster` with
    `TypeError: 'float' object is not subscriptable` -- found by
    tests/test_evaluation_boundaries.py on a 40-cell synthetic case, not
    expected at real single-cell N, but fixed explicitly rather than only
    documented, since a genuinely small diagnostic subsample could hit it)."""
    import harmonypy

    meta = pd.DataFrame({"batch": batch_labels})
    n = X.shape[0]
    nclust = max(2, min(20, n // 10))
    ho = harmonypy.run_harmony(X, meta, ["batch"], nclust=nclust, random_state=seed, verbose=False)
    return np.asarray(ho.Z_corr).T


def save_umap_coords(out_dir: Path, name: str, umap_coords: np.ndarray, cell_ids: np.ndarray,
                       batch_labels: np.ndarray, bio_labels: np.ndarray,
                       patient_ids: np.ndarray) -> None:
    """Cell-ID-bound UMAP coordinates for both the uncorrected and (when
    --run-harmony) Harmony-corrected representation -- v2 only saved the
    uncorrected UMAP; a 2026-10-09 review asked for both so a before/after
    visual comparison is actually possible, not just the metrics table."""
    plot_df = pd.DataFrame({"cell_id": cell_ids,
                              "UMAP1": umap_coords[:, 0], "UMAP2": umap_coords[:, 1]})
    plot_df["cohort"] = batch_labels
    plot_df["outcome"] = bio_labels
    plot_df["patient_id"] = patient_ids
    plot_df.to_csv(out_dir / f"{name}_umap_coords.csv", index=False)


def save_representation(out_dir: Path, name: str, X: np.ndarray, cell_ids: np.ndarray) -> None:
    """Persist a representation bound to its exact cell-ID ordering, and
    immediately reload to verify identity/dimensions/values survive the
    round trip -- a 2026-10-09 review found bare .npy saves had no such
    binding."""
    out = out_dir / f"{name}.npz"
    np.savez(out, X=X, cell_ids=np.asarray(cell_ids, dtype=str))
    reloaded = np.load(out, allow_pickle=False)
    if reloaded["X"].shape != X.shape:
        raise ValueError(f"{out}: write/reload shape mismatch {reloaded['X'].shape} != {X.shape}")
    if not np.array_equal(reloaded["cell_ids"], np.asarray(cell_ids, dtype=str)):
        raise ValueError(f"{out}: write/reload cell_id mismatch")
    if not np.allclose(reloaded["X"], X, equal_nan=True):
        raise ValueError(f"{out}: write/reload value mismatch")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--h5ad", required=True, type=Path)
    p.add_argument("--out-dir", required=True, type=Path)
    p.add_argument("--population-stage", required=True,
                    help="Free-text label naming the true QC stage of --h5ad, e.g. "
                         "'synthetic_fixture' or 'pre-Gate-2 decontx_input PROVISIONAL'. "
                         "Purely descriptive -- checked but only WARNED on mismatch; see "
                         "--stage-id for the hard-fail machine-checkable version.")
    p.add_argument("--stage-id", required=True,
                    choices=sorted(KNOWN_STAGE_IDS),
                    help="Machine-checkable stage identifier, matched EXACTLY against the "
                         "input's own uns['population_stage_id'] -- a mismatch HARD FAILS "
                         "(unlike --population-stage above). A 2026-10-09 review found the "
                         "free-text label alone let a documented command's own wording drift "
                         "from the data's and only warn, which following the documented "
                         "command would itself trigger.")
    p.add_argument("--code-revision", default=None,
                    help="Git SHA or other identifier for the code actually running. If "
                         "omitted, this script tries `git rev-parse HEAD` in its own repo "
                         "checkout; on a loose-file scratch deployment that will fail, and "
                         "the run FAILS rather than publish unresolved provenance silently.")
    p.add_argument("--expression-source", default=None,
                    help="Expected value of the input's uns['expression_source'] (e.g. "
                         "'decontx_counts'), if the caller wants this checked. None = skip "
                         "the check (the input may not declare this field at all).")
    p.add_argument("--lineage-col", default="lineage")
    p.add_argument("--batch-col", default="cohort")
    p.add_argument("--patient-col", default="patient_id")
    p.add_argument("--outcome-col", default="outcome")
    p.add_argument("--exclude-outcome-values", nargs="*", default=["unknown"],
                    help="Outcome values treated as missing-for-biological-scoring, not a "
                         "real class -- excluded from the diagnostic subsample before any "
                         "outcome-stratified computation. Excluded count is reported.")
    p.add_argument("--lineages", nargs="*", default=None,
                    help="Lineage values to analyze; default = every value present.")
    p.add_argument("--embeddings", nargs="*", default=[],
                    help="Additional obsm keys to analyze beyond expression/PCA. These are "
                         "REQUIRED once requested -- a missing or failed one fails the run.")
    p.add_argument("--n-pcs", type=int, default=30)
    p.add_argument("--diagnostic-n-per-patient", type=int, default=20)
    p.add_argument("--min-group-size", type=int, default=10)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--run-harmony", action="store_true",
                    help="Required once set -- a blocked/failed Harmony correction fails the run.")
    args = p.parse_args()

    # 2026-10-09 review: the prior guard (check non-empty, THEN mkdir(exist_ok=True))
    # is check-then-act, not atomic -- two processes launched close together can both
    # pass the check before either writes. The non-empty check still catches leftover
    # files from a COMPLETED or crashed run (unchanged below); claim_ownership() closes
    # the race for two FRESH concurrent starts via an OS-atomic exclusive create, and
    # the ownership record (host/PID/start time/command) lets a later liveness check
    # be aimed at the correct host -- this cluster's process table is not shared
    # across its two login nodes.
    if args.out_dir.exists() and any(args.out_dir.iterdir()):
        raise FileExistsError(
            f"--out-dir {args.out_dir} already exists and is non-empty -- refusing to run "
            f"into it (this project's standing rule: never overwrite a result; use a fresh "
            f"directory)."
        )
    args.out_dir.mkdir(parents=True, exist_ok=True)
    ownership = claim_ownership(args.out_dir / ".run_owner.json", command=" ".join(sys.argv))
    log.info("Claimed run ownership: %s", ownership)

    # Manifest is constructed and written as status="running" BEFORE any fallible
    # I/O (code-identity resolution, input hashing) -- a 2026-10-09 review found the
    # prior version built this dict inline with `_file_sha256(args.h5ad)` OUTSIDE the
    # try/except that writes failure manifests, so a hashing failure (e.g. an
    # unreadable --h5ad) left no persisted record of the run having started at all.
    t0 = datetime.now(timezone.utc)
    manifest = dict(
        population_stage=args.population_stage, stage_id=args.stage_id,
        input_h5ad=str(args.h5ad), input_sha256=None,
        params={k: str(v) for k, v in vars(args).items()},
        started_at=t0.isoformat(), status="running",
        per_lineage=[], errors=[], run_ownership=ownership,
    )

    def write_manifest(status: str):
        manifest["status"] = status
        manifest["finished_at"] = datetime.now(timezone.utc).isoformat()
        with open(args.out_dir / "manifest.json", "w") as f:
            json.dump(manifest, f, indent=2, default=str)

    write_manifest("running")

    try:
        repo_root = Path(__file__).resolve().parents[3]
        identity = _resolve_code_identity(repo_root, args.code_revision)
        manifest.update(identity)
        manifest["code_revision"] = identity["code_revision_asserted"] or identity["code_revision_resolved_git_head"] or "UNRESOLVED"
        if manifest["code_revision"] == "UNRESOLVED":
            raise ValueError(
                "code identity could not be resolved at all (no --code-revision given and "
                "this is not a git checkout) -- refusing to publish a run with no code "
                "identity whatsoever. Pass --code-revision explicitly when deploying from a "
                "loose-file scratch copy rather than a git checkout."
            )
        if not identity["code_revision_verified"]:
            log.warning(
                "code_revision %r is ASSERTED, not VERIFIED (method=%s) -- proceeding (this "
                "is not fatal, unlike an entirely unresolved identity above), but "
                "source_file_sha256 in the manifest is the only immutable record of what "
                "actually ran; cross-check it against `git show <sha>:<path> | sha256sum` "
                "before treating this run's provenance as settled.",
                manifest["code_revision"], identity["code_revision_verification_method"],
            )

        manifest["input_sha256"] = _file_sha256(args.h5ad)
        write_manifest("running")

        log.info("Loading %s", args.h5ad)
        adata = anndata.read_h5ad(args.h5ad)

        if not adata.obs_names.is_unique:
            raise ValueError("Input h5ad has non-unique obs_names (cell IDs) -- refusing to "
                              "proceed; downstream cell-ID alignment would be unsound.")
        for col_arg, col_name in [("lineage_col", args.lineage_col), ("batch_col", args.batch_col),
                                    ("patient_col", args.patient_col), ("outcome_col", args.outcome_col)]:
            if col_name not in adata.obs.columns:
                raise ValueError(f"Required obs column '{col_name}' (--{col_arg.replace('_','-')}) "
                                  f"not found; available columns: {list(adata.obs.columns)}")

        if args.population_stage.strip().lower() != str(adata.uns.get("population_stage", "")).strip().lower():
            log.warning(
                "--population-stage (%r) does not match the input file's own "
                "uns['population_stage'] (%r). Proceeding (this field is descriptive only; "
                "--stage-id below is the hard-fail check), but this is a real discrepancy a "
                "caller should resolve, not silently reconcile.",
                args.population_stage, adata.uns.get("population_stage"),
            )
        check_stage_id(args.stage_id, str(adata.uns.get("population_stage_id", "")),
                        context="run_batch_diagnostics --stage-id vs input uns['population_stage_id']")
        if args.expression_source is not None and args.expression_source != adata.uns.get("expression_source"):
            raise ValueError(
                f"--expression-source {args.expression_source!r} does not match the input's "
                f"own uns['expression_source'] ({adata.uns.get('expression_source')!r}) -- "
                f"refusing to proceed with a declared-but-mismatched expression contract."
            )

        requested_lineages = args.lineages
        lineages = requested_lineages or sorted(adata.obs[args.lineage_col].unique().tolist())
        metrics_rows = []
        requested_embeddings = set(args.embeddings)

        for lineage in lineages:
            lineage_mask = (adata.obs[args.lineage_col] == lineage).values
            n_total = int(lineage_mask.sum())
            if n_total == 0:
                failed = requested_lineages is not None and lineage in requested_lineages
                manifest["per_lineage"].append(dict(
                    lineage=lineage, status="failed" if failed else "skipped",
                    reason="0 cells for this lineage" + (" (explicitly requested)" if failed else ""),
                ))
                continue
            sub = adata[lineage_mask].copy()

            # 2026-10-09 review: this used to exclude unknown-outcome cells from `sub`
            # BEFORE sampling, so they were withheld from PCA/Harmony/technical
            # neighbor-mixing too -- coupling technical-integration eligibility to an
            # unrelated biological label's availability. Fixed: the patient-balanced
            # sample is drawn from every lineage-matched cell regardless of outcome;
            # `bio_known_mask` (derived AFTER sampling, from `diag`) is the only place
            # outcome availability is consulted, and only for the two computations that
            # genuinely need a real biological label (see run_one_embedding).
            sample_idx = patient_balanced_sample(sub.obs, args.patient_col,
                                                  args.diagnostic_n_per_patient, args.seed)
            diag = sub[sample_idx].copy()
            n_excluded_cap = n_total - diag.shape[0]

            batch_labels = diag.obs[args.batch_col].values
            bio_labels = diag.obs[args.outcome_col].values
            cell_ids = diag.obs_names.to_numpy()
            bio_known_mask = ~pd.Series(bio_labels, dtype=str).isin(args.exclude_outcome_values).to_numpy()
            n_bio_known = int(bio_known_mask.sum())

            log.info("Lineage %s: %d total cells, %d patients -> %d-cell diagnostic subsample "
                      "(%d further excluded by the patient-balanced cap, NOT excluded from the "
                      "benchmark); %d of the %d diagnostic cells have a known outcome (not in "
                      "%s) and are eligible for bio-stratified scoring -- the other %d remain "
                      "eligible for technical/neighbor-mixing diagnostics, just not bio scoring",
                      lineage, n_total, sub.obs[args.patient_col].nunique(), diag.shape[0],
                      n_excluded_cap, n_bio_known, diag.shape[0], args.exclude_outcome_values,
                      diag.shape[0] - n_bio_known)

            embeddings_to_run = {"expression_pca": None}
            for key in args.embeddings:
                embeddings_to_run[key] = key

            lineage_record = dict(
                lineage=lineage, n_total_cells=n_total,
                excluded_outcome_values=list(args.exclude_outcome_values),
                n_patients=int(sub.obs[args.patient_col].nunique()),
                n_diagnostic_cells=int(diag.shape[0]), n_excluded_by_cap=int(n_excluded_cap),
                n_diagnostic_cells_bio_known=n_bio_known,
                n_diagnostic_cells_bio_unknown=int(diag.shape[0] - n_bio_known),
                diagnostic_sampling_rule=f"<= {args.diagnostic_n_per_patient} cells/patient, seeded "
                                          f"(outcome-agnostic; bio_known_mask applied only to "
                                          f"batch_asw/bio_conservation, not to sampling or to "
                                          f"neighbor_batch_mixing)",
                embeddings={},
            )

            for label, obsm_key in embeddings_to_run.items():
                is_required = (label == "expression_pca") or (label in requested_embeddings)
                try:
                    if label == "expression_pca":
                        X = compute_pca(diag, args.n_pcs, args.seed)
                    else:
                        if obsm_key not in diag.obsm:
                            lineage_record["embeddings"][label] = dict(
                                status="failed" if is_required else "unsupported",
                                reason=f"obsm key '{obsm_key}' not present"
                                       + (" (explicitly requested)" if is_required else ""))
                            continue
                        X = np.asarray(diag.obsm[obsm_key])

                    result = run_one_embedding(X, batch_labels, bio_labels, bio_known_mask,
                                                 args.seed, args.min_group_size)
                    umap_coords = compute_umap(X, args.seed)

                    corrected_result = None
                    corrected_umap_coords = None
                    harmony_blocked_or_failed = False
                    corrected_umap_ok = None
                    if args.run_harmony:
                        try:
                            X_corr = run_harmony_correction(X, batch_labels, args.seed)
                            corrected_result = run_one_embedding(
                                X_corr, batch_labels, bio_labels, bio_known_mask,
                                args.seed, args.min_group_size)
                            save_representation(args.out_dir, f"{lineage}_{label}_harmony_corrected",
                                                  X_corr, cell_ids)
                            corrected_umap_coords = compute_umap(X_corr, args.seed)
                            corrected_umap_ok = corrected_umap_coords is not None
                            if corrected_umap_coords is not None:
                                save_umap_coords(args.out_dir, f"{lineage}_{label}_harmony_corrected",
                                                   corrected_umap_coords, cell_ids, batch_labels,
                                                   bio_labels, diag.obs[args.patient_col].values)
                        except ImportError as e:
                            corrected_result = dict(status="blocked",
                                                      reason=f"harmonypy not importable: {e}")
                            harmony_blocked_or_failed = True
                        except Exception as e:
                            corrected_result = dict(status="failed", reason=str(e))
                            harmony_blocked_or_failed = True

                    status = determine_embedding_status(
                        asw_skipped_groups=bool(result["batch_asw"]["skipped_groups"]),
                        uncorrected_umap_ok=umap_coords is not None,
                        run_harmony=args.run_harmony,
                        harmony_blocked_or_failed=harmony_blocked_or_failed,
                        corrected_umap_ok=corrected_umap_ok,
                    )

                    save_representation(args.out_dir, f"{lineage}_{label}_uncorrected", X, cell_ids)
                    if umap_coords is not None:
                        save_umap_coords(args.out_dir, f"{lineage}_{label}_uncorrected",
                                           umap_coords, cell_ids, batch_labels, bio_labels,
                                           diag.obs[args.patient_col].values)

                    lineage_record["embeddings"][label] = dict(
                        status=status, uncorrected=result,
                        harmony=dict(
                            label="transductive/exploratory integration -- jointly fit on this "
                                  "diagnostic subsample, NOT an out-of-sample transform",
                            result=corrected_result,
                        ) if args.run_harmony else None,
                    )
                    for scope, res in [("uncorrected", result)] + (
                        [("harmony_corrected", corrected_result)]
                        if args.run_harmony and isinstance(corrected_result, dict)
                        and "batch_asw" in corrected_result else []
                    ):
                        for stratum, val in res["batch_asw"]["asw_batch_per_group"].items():
                            metrics_rows.append(dict(lineage=lineage, embedding=label, scope=scope,
                                                       metric="batch_asw_abs_silhouette",
                                                       stratum=stratum, value=val))
                        metrics_rows.append(dict(lineage=lineage, embedding=label, scope=scope,
                                                   metric="batch_asw_neg_scib", stratum="__all__",
                                                   value=res["batch_asw"]["batch_asw_neg_scib"]))
                        metrics_rows.append(dict(lineage=lineage, embedding=label, scope=scope,
                                                   metric="neighbor_batch_mixing", stratum="__all__",
                                                   value=res["neighborhood_mixing"]["mean_neighbor_batch_mixing"]))
                        metrics_rows.append(dict(lineage=lineage, embedding=label, scope=scope,
                                                   metric="bio_ari", stratum="__all__",
                                                   value=res["bio_conservation"]["ari"]))

                except Exception as e:
                    log.exception("Embedding %s failed for lineage %s", label, lineage)
                    lineage_record["embeddings"][label] = dict(status="failed", reason=str(e))

            manifest["per_lineage"].append(lineage_record)

        metrics_df = pd.DataFrame(metrics_rows, columns=METRICS_COLS)
        metrics_out = args.out_dir / "metrics_long.csv"
        metrics_df.to_csv(metrics_out, index=False)

        any_required_failed = any(
            emb.get("status") == "failed"
            for rec in manifest["per_lineage"] for emb in rec.get("embeddings", {}).values()
        ) or any(rec.get("status") == "failed" for rec in manifest["per_lineage"])
        any_optional_issue = any(
            emb.get("status") in ("partial", "unsupported")
            for rec in manifest["per_lineage"] for emb in rec.get("embeddings", {}).values()
        )
        if any_required_failed:
            final_status = "failed"
        elif any_optional_issue:
            final_status = "complete_with_partial_results"
        else:
            final_status = "complete"
        write_manifest(final_status)

        log.info("Wrote %s and manifest.json. status=%s", metrics_out, final_status)
        print(f"\n=== STAGE-33 DIAGNOSTIC RUN: {final_status} ===")
        print(metrics_df.to_string(index=False) if len(metrics_df) else "(no metrics rows)")

        if final_status == "failed":
            sys.exit(1)
    except Exception:
        log.exception("run_batch_diagnostics failed")
        write_manifest("failed")
        raise


if __name__ == "__main__":
    main()

# 34_scfm_discovery_framework

Model-agnostic foundation for the scFM biological-discovery framework (owner scope
pivot, 2026-10-12 — see `docs/decisions.md`). Builds the dataset/model-registry
contracts, `validate`/`plan` entry points, and artifact-provenance schema ahead of
any specific model adapter, per `review:20261012_145004` item C.

## What exists here

- `code/dataset_schema.py` — dataset config validation. Every required field is
  type-checked (string, non-empty), not just non-None; outcome/severity is
  accepted only as optional metadata; `counts_slot` (`"X"`/`"raw.X"`/
  `"layers/<name>"`) states exactly where counts live.
- `code/model_registry.py` — model/checkpoint identity, with pretraining-species
  evidence and input-species-acceptance kept as separate fields. Two entries
  (`scprint2_small_v2`, `uce_4layer`) reflect `MODEL_SELECTION_v1.md`'s
  corrected findings; neither is the framework's confirmed backbone.
  scPRINT-2's checkpoint-metadata-confirmed human+mouse coverage and its
  Ensembl-ID/HUGO-symbol input requirement were independently re-verified
  against the checkpoint file and the README FAQ, respectively.
- `code/artifact_schema.py` — provenance manifest (`planned`/`executed`/
  `validated`); nothing here can produce `executed`/`validated` without a real
  run asserting it. A shared `validate_manifest()` enforces field *types*
  (e.g. `source_code_revision_or_hash` must be a string, each `dependencies`
  entry must be `{"file": str, "sha256": str}`) at build, write, AND read —
  a hand-tampered on-disk file is rejected on read, not just on construction.
  `read_manifest` also rejects NaN/Infinity JSON constants, matching
  `write_atomic`'s `allow_nan=False`.
- `code/enrichment_providers.py` — GSEA (method) vs. GO/KEGG (annotation source)
  contracts. `AnnotationProvider` (live GO/KEGG/MSigDB fetch) is still a
  `NotImplementedError` stub this round. `GeneRanking` now requires an explicit
  `species` field (never defaulted) and `validate()` rejects non-finite
  scores, empty/non-string gene IDs, and a fully-empty ranking, without
  requiring both positive and negative values.
- `code/gene_set_cache.py` — imports a LOCAL `.gmt` file into a new, versioned
  cache directory (never a URL, never an online library name). Records
  collection/species/namespace/version and the source file's own content
  SHA-256, with a LOCAL `imported_utc` (never a fabricated live-retrieval
  timestamp — there is no network fetch here). `load_cached_gene_sets`
  re-verifies that hash before ever returning gene sets, with no fallback to
  a live fetch anywhere in the module.
- `code/gsea_runner.py` — `LocalPrerankGSEAProvider`, the real `GSEAProvider`
  implementation, backed by the installed `gseapy.prerank` (version resolved
  and recorded at call time). `gene_sets` is always passed to gseapy as an
  in-memory dict loaded from the hash-verified cache above — never a path or
  library-name string, the two shapes gseapy could read as an Enrichr
  request. Distinguishes `"complete"` / `"failed"` (invalid ranking,
  species/namespace mismatch, no gene set surviving size filtering against
  the exact ranking intersection, an identity-consistency violation, or a
  backend error) / `"no_overlap"` (zero genes shared with the whole
  collection — distinct from a completed run with no significant terms) /
  `"not_run"`. Records every parameter that could change the numeric result
  (seed, threads, permutation count, min/max set size, weight, ranking
  direction, the backend's own p-value/FDR definitions, a tied-score count)
  in the persisted result.

  Per `review:20261012_220407`'s review of the first implementation
  (F1-F5), the installed gseapy's own `base.py` was found to silently
  uppercase gene IDs for matching when the ranking isn't mostly-uppercase
  but the gene sets sampled are — an implicit identity rewrite. Every gene
  set is now pre-filtered to its EXACT, case-sensitive intersection with the
  ranking before gseapy ever sees it (with a defense-in-depth check on the
  backend's own reported leading edge), so `intersection_size` and
  `leading_edge_genes` can no longer disagree about which genes matched.
  Every run also freezes the FULL ranking (exact gene/score order, the tied
  score count, which genes were excluded from a given gene set and why) to
  `RANKING_SNAPSHOT.json`, binds the wrapper's own source revision and
  dependency versions, and a brand-new process can replay the whole run
  from those frozen artifacts alone via `code/gsea_replay.py` — never
  calling back into whatever generated the original ranking. A
  `validate_gsea_result()` shared by the writer and the reader rejects an
  incomplete "complete" record, an impossible probability field, or a
  non-finite value anywhere in the payload (not only in `results_table`).
  The gene-set cache's own on-disk path is stored relative to its manifest's
  directory, so a cache copied elsewhere or loaded from a different cwd
  still resolves. **This backend numeric chain has been run end to end on
  synthetic data; it is not yet a frozen evidence package that has been
  reviewed and accepted for a real dataset** — see "What does not exist
  here yet" below.

  Per `review:20261013_114111`'s review of the above (R1/R2/C1): a
  numerically-identical replay was not proof that it used the SAME ranking,
  cache, or source/dependency identity as the original run — a
  legitimately re-hashed but DIFFERENT ranking snapshot, a swapped cache
  with a different collection/version, and a tampered source-revision
  field were all previously accepted as "REPLAY MATCH". `gsea_replay.py`
  now compares the reconstructed ranking's semantic hash and the snapshot
  file's own hash against what the ORIGINAL result itself recorded (not
  just the snapshot's own internal self-check); the cache's full semantic
  identity (species/namespace/source_version, not only collection
  name/content hash) is bound into the result and checked before replay —
  a verbatim-copied, relocated cache still passes, a substituted one does
  not; every numeric run parameter must be present with the correct type
  (a probe found `source_code_revision_or_hash=True` and
  `params={"unrelated": 1}` both previously accepted); and a source-code or
  dependency-version drift between the original run and the replay
  environment is a strict-replay failure by default, with a separate,
  explicitly-labeled `--allow-identity-drift` flag ("REPLAY COMPARISON",
  never "REPLAY MATCH") for the deliberate cross-version-numbers-only case.
  A NaN/Inf score previously made the ranking-snapshot writer raise
  mid-write (`allow_nan=False`) BEFORE `ranking.validate()` ever ran,
  leaving a bare `.inprogress` file and no `GSEA_RESULT.json` at all —
  fixed by validating first, so an invalid ranking now gets a clean
  `status="failed"` record with no snapshot and no leftover temp file. The
  F1 pre-filter eliminates CROSS-set leakage but does not by itself resolve
  a case COLLISION within a ranking's own distinct gene IDs (e.g. `ANCHOR`
  and `Anchor` both present) — this still stays correctly fail-closed via
  the separate, defense-in-depth leading-edge check (confirmed by a direct
  regression against the real backend), and the module's own documentation
  has been corrected to state that scope accurately rather than claim the
  pre-filter alone makes it impossible.

  Per `review:20261013_135634`'s independent re-review of the
  above (R3/R4): `gsea_replay.py`'s strict replay is now explicitly
  restricted to a `status="complete"` original — a probe found that for
  `"no_overlap"` (reached before the cache was ever checked against the
  original), a cache swapped for a different collection/version still
  printed "REPLAY MATCH" because the cache-identity check only ran inside
  an `if status == "complete"` branch. Any non-`"complete"` original is now
  refused outright, before the snapshot is read or the backend is ever
  called — the original record is left exactly as it was. Separately, the
  printed "REPLAY MATCH"/"REPLAY COMPARISON" label previously depended on
  whether a drift happened to be *observed*, not on whether the caller
  explicitly asked for `--allow-identity-drift`'s cross-version-comparison
  mode — so that flag against two identical environments still printed
  "REPLAY MATCH", contradicting the CLI's own help text. The label now
  reflects the mode the caller actually requested.
- `code/cli.py` — `validate`/`plan` subcommands. Never submits a job, downloads
  a weight, or reads a full real matrix; structural h5ad checks use bounded
  `h5py` reads (never `anndata(..., backed="r")`, which the batch-effect-audit
  work in this project found still loads every layer in full). A relative
  `expression_matrix_path` resolves against the dataset CONFIG file's own
  directory, never the caller's cwd. The declared `counts_slot` is mapped to
  its OWN gene axis (`var` for `"X"`/`"layers/*"`, `raw/var` for `"raw.X"` —
  never compared against the wrong one), and an unresolvable obs/var index is
  reported as an explicit failure, never silently skipped.
- `examples/run_example.py` — a real, runnable, end-to-end demonstration: builds
  a tiny synthetic fixture (donor `D1` contributing two samples, a genuine
  biological-replicate example), then runs
  `validate --check-matrix-structure -> plan --check-matrix-structure ->
  read_manifest` and asserts the read-back object exactly equals what was
  written to disk.
- `examples/run_gsea_example.py` — a real, runnable, FULLY OFFLINE demonstration
  of the gene-set cache + GSEA chain: builds a tiny synthetic `.gmt` (one gene
  set deliberately boosted to the top of a ranking by construction, one
  unrelated — both named `TOY_*` so they are never mistaken for real GO/KEGG
  data), imports it into a local cache, runs the real `gseapy`-backed
  provider, persists `GSEA_RESULT.json`, reads it back, asserts the boosted
  set's ES AND NES both exceed the unrelated set's and that all of its genes
  are in its own leading edge, and finally replays the whole run from its
  frozen artifacts alone in a brand-new process (`code/gsea_replay.py`),
  confirming every scientific field matches exactly.
- `examples/run_pca_example.py` — a real, runnable, FULLY OFFLINE
  demonstration of `PCAEmbeddingProvider` on a tiny synthetic HUMAN dataset
  and a tiny synthetic MOUSE dataset through the SAME interface: builds a
  synthetic h5ad (`TOYGENE*`/`TOYCELL*` IDs, random counts), validates it
  via the same `dataset_schema`/`cli` checks the rest of this framework
  uses, runs the real `sklearn`-backed PCA, and reads the persisted
  manifest + npz back, asserting the round-trip shape. Explicitly labeled
  synthetic throughout — makes no claim about real study data or any
  biological finding.
- `code/gsea_replay.py` — replays a previously-run GSEA result from
  `RANKING_SNAPSHOT.json` + `GSEA_RESULT.json` alone, in a new process, never
  calling back into whatever generated the original ranking. Compares status,
  ES, NES, p-value, FDR, FWER p-value, set sizes, and leading-edge genes
  between the original and the replay. Per
  `review:20261013_114111` R1: before ever running the backend, it
  strictly checks that the reconstructed ranking and the selected gene-set
  cache's full semantic identity match what the ORIGINAL result itself
  recorded, and requires every numeric run parameter present with the
  correct type — never a numerically-coincidental match on a DIFFERENT
  analysis, cache, or environment. A relocated (copied) cache directory
  still replays; a substituted collection/version does not. A source-code
  or dependency-version drift between the original run and the replay
  environment fails a default, strict replay; `--allow-identity-drift`
  tolerates it explicitly and relabels the result "REPLAY COMPARISON"
  rather than "REPLAY MATCH" (always, per R4 above — regardless of whether
  a drift was actually observed). Per R3, only a `status="complete"`
  original can be replayed; any other status is refused outright, before
  touching the snapshot, cache, or backend.
- `code/embedding_providers.py` — `PCAEmbeddingProvider`, a REAL
  `scikit-learn`-backed PCA baseline, per `review:20261013_135634`
  section 3: an explicit CONTROL ("what can a simple expression
  representation already do on this input/downstream question"), never
  treated as proof that a scFM is superior, and never conflated with
  model-masking evidence. Reuses `dataset_schema.validate_dataset_config`
  and `cli.py`'s own structural/gene-axis checks — never given a bare file
  path. Reads the declared `counts_slot` (dense or `csr_matrix`/
  `csc_matrix`) read-only, bounded to an explicit
  `_MAX_DENSE_ELEMENTS` this round (a larger block is refused outright, not
  silently truncated). Feature selection is either "all genes", an explicit
  caller-supplied `feature_ids_override`, or a label-free top-variance
  `hvg_top_n` — there is no disease/outcome-label parameter anywhere on
  this interface. Never modifies the input matrix, never treats PCA as
  batch correction. Persists `EMBEDDING_RESULT.json` (manifest: params,
  dataset identity, source/dependency versions, an `input_counts_sha256`
  binding the EXACT post-selection numeric block used) plus
  `embedding_matrices.npz` (embedding, components, mean, singular values),
  both re-validated and hash-checked on read (`read_embedding_result` +
  `load_embedding_matrices`).

  Per `review:20261013_210131` F2: the counts_slot's own gene axis
  is now validated (non-empty, every value a non-empty string, UNIQUE)
  BEFORE any ID->position map, `feature_ids_override`, or `hvg_top_n`
  selection is built from it — a probe found a duplicated source gene ID
  silently resolved to whichever column `{g: i for i, g in
  enumerate(gene_ids)}` happened to keep last, with no indication the other
  column was dropped. A duplicate is now refused with a clean `failed`
  record; this module never auto-merges, suffixes, or averages duplicated
  gene IDs. Per F3: every array about to be PERSISTED (embedding,
  components, mean, singular values, AND explained variance ratio) is now
  checked for finiteness before anything is written — a zero-total-variance
  input (e.g. an all-zero counts block) previously reached a NaN
  `explained_variance_ratio_`, wrote `embedding_matrices.npz`, and only
  then raised on the JSON manifest's own validation, leaving an orphaned
  npz with no `EMBEDDING_RESULT.json` beside it. A zero-total-variance
  input is now refused explicitly and up front. A genuinely unexpected
  backend exception is recorded as a `failed` result first, then
  re-raised.
- `code/gene_ranking_providers.py` — the two gene-ranking EVIDENCE STREAM
  contracts, per `review:20261013_135634` section 3:
  `ExpressionDERankingProvider` ("which genes' expression changed, under an
  explicit comparison" — the intended estimator is donor/animal-aware
  pseudobulk + a shrunken LFC-style estimator, never independent-cell
  Wilcoxon or an arbitrary signed `-log10(p)` substitute by default) and
  `ModelMaskingRankingProvider` (`Delta(g) = s_full - s_masked(g)` for a
  FIXED, bound scoring function — masking a gene is masking a GENE INPUT to
  the model; this interface has no cell-selection/QC-mask parameter at
  all). Both produce the same `enrichment_providers.GeneRanking` the GSEA
  chain already consumes, bound with their own `evidence_type` and
  comparison/score/sample-structure metadata via `make_evidence_ranking`.
  This round implements TYPES/VALIDATION/SERIALIZATION ONLY — `run()`
  always returns `status="failed"` (a structural/design problem) or
  `status="not_run"` (inputs valid, no real backend yet), never a
  fabricated `"complete"` result; no main model checkpoint is pinned and no
  donor/animal grouping, masking score, or comparison has been confirmed by
  the owner.

  Per `review:20261013_210131` F1: a donor-blocked `~ donor +
  condition` estimability check (`check_design_estimability`) is no longer
  applied unconditionally to every comparison — a probe found that three
  independent healthy donors vs. three independent disease donors, an
  ordinary UNPAIRED design, was refused by a donor-blocking check that
  should never have applied to it, while a single cross-group donor out of
  several was previously enough to call a donor-blocked design "estimable"
  outright. `SampleStructure.design_type` (one of
  `"independent_donors_unpaired"` / `"paired_or_repeated_within_donor"`,
  EXPLICIT, never inferred from an `estimator=` string) now decides which
  check applies, and each branch is checked against the REAL donor/group
  input, not just trusted from the declaration. `comparison.group_column`
  disagreeing with `sample_structure.group_column` is now refused;
  declared-but-unchecked `covariate_columns` are flagged UNVERIFIED in the
  reason rather than silently passed. Per F4: `ExpressionDEResult` and
  `ModelMaskingResult` now bind a content identity of their ACTUAL inputs
  (`group_values_sha256`/`donor_values_sha256`; `perturbable_genes` + its
  hash) — a probe found two different `perturbable_genes` inputs
  serializing to an identical record.

  Per `review:20261013_222450` R1: the F1 fix above only checked
  that SOME compared-group sample existed, never that BOTH `reference_
  level` and `comparison_level` were present (a `healthy`-only input passed
  as `not_run`), and an empty-string/`None` donor or group value reached
  `sorted()` uncaught, raising a bare `TypeError`. A new `_validate_group_
  donor_values` check (length match, non-empty string values, BOTH levels
  present) now runs inside `check_design_estimability` before any donor
  grouping or `sorted()` call. Per R2: the `paired_or_repeated_within_
  donor` refusal message claimed a partially-paired design "cannot be
  estimated", which is not generally true (a `~ donor + condition` design
  matrix CAN remain full rank with some, not all, donors paired — proven
  with a full-rank 6-row/5-column counterexample in the test suite). The
  message now states this round only SUPPORTS the fully-paired case
  ("UNSUPPORTED by this round's check", never "mathematically
  inestimable"), and that a straddling/non-straddling donor is never
  silently dropped by the check. Per R3: `ExpressionDEResult`/
  `ModelMaskingResult` carried only a length+hash summary with no way to
  recover the actual input from the record, no dedicated writer/reader,
  and an identity hash built with `json.dumps(..., default=str)` (silently
  stringifies an unrecognized value). Both results now carry a round-
  trippable `input_mapping`/`perturbable_genes` field plus a recomputed-
  and-compared content hash; new `write_expression_de_result`/
  `read_expression_de_result` and `write_model_masking_result`/
  `read_model_masking_result` functions validate identically at write and
  read (mirroring `embedding_providers.write/read_embedding_result`) and
  reject a tampered mapping/hash on read; `_safe_json_identity` returns an
  explicit "identity unavailable" state rather than ever silently
  stringifying a non-JSON-native value (now `artifact_schema.safe_json_
  identity`, reused by `pseudobulk_providers.py` too).

  Per `review:20261014_100436` F3: the write/read validators were
  still asymmetric — `status='complete', ranking={}, backend='not_
  implemented'` was ACCEPTED (only `ranking is not None` was checked, and
  an empty dict is not `None`), while a GENUINE `status='failed'` record
  carrying a diagnostic `donor_value=None` (the real provider's own record
  of that exact invalid input) was REJECTED by the same strict per-row
  checks a `'complete'`/`'not_run'` record needs. Both are fixed together:
  `ranking` is now re-validated against `GeneRanking.validate()` itself;
  `backend='not_implemented'`/`backend_version='n/a'` are refused on any
  `status='complete'` record; `MaskingScoreSpec.validate()` now checks
  TYPE, not truthiness (`score_definition=True` is truthy but not a
  string); the strict per-row/per-field VALUE checks on `input_mapping`/
  `masking_score`/`comparison`/`sample_structure` now apply only when
  `status in ('complete', 'not_run')` — a `'failed'` record's raw
  (possibly invalid) input is a diagnostic snapshot, never a record
  claiming to BE validated; the key-presence/hash checks still apply to
  every status unconditionally. `n_samples` vs. `len(input_mapping)` and
  nested `sample_structure.design_type` are now re-checked on every read,
  not only at dataclass construction.

  Per `review:20261014_114357` R2/R3: `n_samples`'s own TYPE was
  never checked (`isinstance(n_samples, int) and ... != ...` had no `else`
  branch, so `n_samples='two'`/`None` silently skipped the consistency
  check — now checked unconditionally first), and the `n_samples ==
  len(input_mapping)` consistency check ran UNCONDITIONALLY, which made a
  genuinely `'failed'` length-mismatch record (`group_values`/
  `donor_values` of different lengths, so `input_mapping` is correctly
  forced empty while `n_samples` stays the raw non-zero `group_values`
  length) unwritable — the same asymmetry F3 fixed for values, recurring
  for a length relationship. `ModelMaskingResult.perturbable_genes` had the
  matching gap (its per-entry value check ran unconditionally, rejecting a
  genuine `perturbable_genes=(None,)` failed record). Fixed by new
  `input_lengths`/`input_mapping_buildable` fields (the raw lengths and
  whether they agreed, recorded ALWAYS, never dropped on failure) that the
  length-consistency check is now gated on, and by gating `perturbable_
  genes`' value check the same way `masking_score`'s fields already were.

  Per `review:20261014_151747` P2: the R3 fix above left one gap
  — `input_mapping_buildable=False` (a genuine length disagreement) was
  cross-checked against `input_mapping`'s shape, but nothing tied it to
  `status` itself. A genuinely `'failed'` length-mismatch record, written
  and read back correctly, could have ONLY its `status` field rewritten to
  `'not_run'`/`'complete'` on disk and still be ACCEPTED, since every other
  cross-check still held. Fixed: `validate_expression_de_result` now
  rejects `input_mapping_buildable=False` combined with
  `status in ('not_run', 'complete')` directly — this module defines
  `'not_run'` as "inputs PASSED validation, no backend ran yet", which an
  unbuildable input never satisfies.
- `code/de_design_providers.py` — see "What does not exist here yet"
  below for the "B" increment this file implements.
- `code/pseudobulk_providers.py` — `PseudobulkAggregationProvider`, the
  first real increment toward `ExpressionDERankingProvider`'s intended
  backend, per `review:20261013_222450` B: sample x
  (caller-selected) cell-type/state RAW count-sum aggregation ONLY — no DE
  statistic, p-value, log-fold-change, or `GeneRanking`;
  `ExpressionDERankingProvider` itself is untouched. Per Bioconductor
  OSCA's own convention, the aggregation unit is always `dataset_config
  ["sample_id_column"]`'s value (never donor/animal ID), so a donor's
  repeat visits never silently merge into one row. Requires
  `counts_semantics == "raw_counts"` and non-negative, integer-valued
  counts explicitly — never coerces by rounding or relabeling
  `counts_semantics`. The caller must explicitly select which cells go
  into the aggregation (no "all cells" default); metadata (donor,
  condition, visit/state) that conflicts across one sample_id's cells is
  refused outright, never resolved by taking the first/last value seen. A
  sample_id with zero selected cells is simply absent from the output,
  never fabricated as an observed all-zero row. Reports
  `n_pseudobulk_samples` and `n_distinct_donors` SEPARATELY (one donor can
  contribute more than one sample). Reuses `embedding_providers`'s F2
  gene-axis validation and dense counts reader, and shares the same
  strict write/read/tamper-rejection contract (`write_pseudobulk_result`/
  `read_pseudobulk_result`/`load_pseudobulk_counts`), including a
  counts-conservation check (selected-cell total == aggregated total).

  Per `review:20261014_100436` F1: the manifest's `sample_rows`/
  `gene_ids` LABELS carried no binding to the npz array's actual row/
  column CONTENT — reversing either, or editing `n_distinct_donors`/a
  row's `library_size`/`condition_value`, passed unchanged (the npz
  file-hash check proves only that the npz BYTES are unmodified, never
  that the JSON's labels for those bytes are correct). Fixed with a new
  `identity_sha256` field (a hash over a canonical snapshot of `gene_ids`/
  `gene_total_counts`/`sample_rows`/the derived counts, recomputed and
  compared on every read) plus `load_pseudobulk_counts` now ALSO
  recomputing each row's library size and each gene's total count
  directly from the loaded npz array and comparing them against the
  manifest's claims — binds the JSON labels to the array's real content,
  not just to the JSON's own self-consistency. `dataset_identity` gained
  `selected_raw_counts_sha256`/`selected_cell_metadata_sha256`/`gene_ids_
  sha256` — content hashes of the exact selected-cell raw-count block and
  per-cell metadata actually read this run (never re-reads the whole
  source file; a config hash is not a content hash of the expression
  matrix). Per F2: `cell_type_or_state_label=True` passed the entry
  truthiness check and wrote an orphan npz before the writer's own type
  check rejected it; an explicitly-requested `condition`/`visit` column
  with a missing value either silently wrote `null` or crashed `sorted()`
  with an uncaught `TypeError` on a `None`/`str` mix. Fixed by TYPE-
  checking every entry argument and refusing (never raising) a missing/
  non-string value for an explicitly-requested column before any
  grouping, plus splitting `_write_pseudobulk_result`'s validation into a
  pre-npz core check (`_validate_pseudobulk_result_core`) as a second,
  independent guard against an orphan npz.

  Per `review:20261014_114357` R1: the F1 fix above was itself
  insufficient against an EQUAL-MARGINALS swap — a fixture with two rows
  whose library sizes AND per-gene totals are identical
  (`[[1, 9], [9, 1]]`, both summing to 10) can have `sample_rows` reversed,
  `gene_ids` reversed, or a `condition_value` rewritten, while every
  recomputed row/column TOTAL (and a correspondingly self-consistently
  rehashed `identity_sha256`) still agrees — a total is a many-to-one
  projection of content, never a content identity. Fixed by saving
  `sample_ids`/`gene_ids`/a normalized per-row `sample_identity_json`
  snapshot (donor/condition/visit/cell_ids) INSIDE the npz itself (plain
  numpy string arrays, no pickle, covered by the same whole-file
  `aggregation_sha256`); `load_pseudobulk_counts` now compares the
  manifest's claimed order/metadata DIRECTLY against this independent,
  hash-protected in-npz copy (`_npz_identity_binding_problems`), not only
  against totals recomputed from `counts`. Separately, `n_pseudobulk_
  samples`/`n_distinct_donors` in `_validate_pseudobulk_result_core` had
  the same unconditional-type-check-missing shape as R2 — `read_
  pseudobulk_result` alone (no npz access) would have accepted e.g.
  `n_pseudobulk_samples="two"` unchanged; both fields are now type-checked
  unconditionally first.
- `examples/run_pseudobulk_example.py` — a standalone, runnable human+mouse
  example through the real `PseudobulkAggregationProvider` interface,
  including a failed-record read-back path (no DE statistic computed).
- `MODEL_SELECTION_v1.md` — the limited backbone shortlist (recommendation, not
  a selection).

Tests: `tests/test_scfm_*.py` (project-wide `tests/` convention), run with
`python -m pytest tests/test_scfm_*.py -v` inside `302006_env_clean` on HPCC.
Every test uses tiny synthetic fixtures — never real expression data, never
`sbatch`.

## Usage

```bash
python code/cli.py validate --dataset-config <config.json> \
    [--check-matrix-structure] [--model-registry-key scprint2_small_v2]

python code/cli.py plan --dataset-config <config.json> --out <new/empty/dir> \
    [--check-matrix-structure] [--model-registry-key scprint2_small_v2]
```

A relative `expression_matrix_path` inside the config resolves against the
config file's own directory. `counts_slot` is required and must be `"X"`,
`"raw.X"`, or `"layers/<name>"`.

To see the whole chain run end to end against a real (synthetic) config and
fixture:

```bash
python examples/run_example.py --out-dir <new/empty/dir>
```

To see the local gene-set cache + offline GSEA chain run end to end:

```bash
python examples/run_gsea_example.py --out-dir <new/empty/dir>
```

```bash
python code/gene_set_cache.py  # (not a CLI yet -- import_gene_sets()/load_cached_gene_sets()
                                #  are called directly; see run_gsea_example.py for usage)
```

To see the PCA baseline run end to end against a tiny synthetic human AND
mouse dataset through the same interface:

```bash
python examples/run_pca_example.py --out-dir <new/empty/dir>
```

To see the raw-count-sum pseudobulk aggregation run end to end against a
tiny synthetic human AND mouse dataset through the same interface,
including a genuinely invalid input producing a readable `failed` record:

```bash
python examples/run_pseudobulk_example.py --out-dir <new/empty/dir>
```

To replay a previously-run GSEA result from its frozen artifacts alone, in a
new process (never calling back into whatever generated the original
ranking):

```bash
python code/gsea_replay.py --from <original_gsea_out_dir> --out <new/empty/dir> \
    [--cache-dir <relocated_cache_dir>] [--allow-identity-drift]
```

`--cache-dir` must point at the SAME collection/version/species/namespace/
content the original run used — a copy of that exact cache directory moved
elsewhere is fine; a different collection or version is refused.
`--allow-identity-drift` is the separate, explicitly-labeled cross-version
comparison mode (reported as "REPLAY COMPARISON", never "REPLAY MATCH")
for deliberately checking whether two different environments produce the
same numbers; it is never the default.

## What does not exist here yet

No model adapter, no real inference, no Docker packaging, no live GO/KEGG
annotation fetch (`AnnotationProvider` is still a stub). No real dataset
config for the human monocyte application has been
written against real data this round, and no real (non-synthetic) gene set
or ranking has been run through the GSEA chain. `PCAEmbeddingProvider` is a
real implementation but has only been run on synthetic data this round —
no real preprocessing/dimensionality decision has been made for any real
dataset; its gene-axis duplicate check and zero-total-variance/finiteness
checks (`review:20261013_210131` F2/F3) are a bounded-input
contract, not a general-purpose large-matrix pipeline — the dense-read
bound still applies to the FULL counts_slot even when `cell_ids_subset` is
small. Neither `ExpressionDERankingProvider` nor
`ModelMaskingRankingProvider` has a real statistical/model backend yet
(see `code/gene_ranking_providers.py`'s own docstring) — both are
type/validation/serialization contracts only; `check_design_estimability`'s
`design_type` branching (F1) checks group/donor structure only, not a full
design matrix, and declared `covariate_columns` are explicitly UNVERIFIED
until a real backend accepts covariate values. `check_design_estimability`
only SUPPORTS a fully-paired design under `"paired_or_repeated_within_
donor"` (per R2 above) — a partially-paired design is refused as
unsupported by this round's check, not because it is necessarily
mathematically inestimable. `make_evidence_ranking`'s
`ranking_definition` substring check (F4) is a lint-level tripwire, not a
structured source/provenance verification. `code/pseudobulk_providers.py`
implements RAW count-sum aggregation only — no shrunken-LFC or other DE
statistic, no p-value, no `GeneRanking`; it is the first increment toward
`ExpressionDERankingProvider`'s intended backend, not that backend itself.
Its manifest's `sample_rows`/`gene_ids` order and `dataset_identity`'s
content hashes (per `review:20261014_100436` F1, directly bound
inside the npz itself per `review:20261014_114357` R1) record
what was read/computed this run and are cross-checked against the loaded
npz array's own stored identity arrays and against each other on every
read — they do NOT themselves prove that an independent re-read of the
original source file would reproduce the same content; that would be a
separate replay step this round does not implement. No DE backend (PyDESeq2, an R/DESeq2 subprocess, or a custom
shrinkage estimator) has been chosen or installed; a read-only dependency
check found the clean conda env (Python 3.10.18) has neither
`pydeseq2` nor `rpy2` installed. Per `review:20261014_151747`
step 3 and `DE_BACKEND_ENVIRONMENT_PLAN_v1.md`: PyPI's own metadata
(fetched 2026-10-14) shows `pydeseq2` raised its `requires_python` to
`>=3.11` starting at `0.5.3` -- but `0.5.2` (the immediately prior
release) still declares `requires_python>=3.10.0`, so a Python-3.10.18
install is possible WITHOUT upgrading the shared env's Python; `0.5.2` is
this plan's recommended pin, not the PyPI-latest `0.5.4`. Also per that
document: an earlier round's "R/module not found" observation was from a
single NON-interactive SSH shell and was correctly left unconfirmed either
way (`DEVELOPMENT.md` rule 3) -- this round checked an interactive LOGIN shell's
module system directly and found R IS genuinely available on this cluster
(`R/4.1.2` through `R/4.6.0`, `R/4.5.2` default), so an R/DESeq2/`rpy2`
path remains real, just not this plan's primary recommendation (PyDESeq2
avoids a second language boundary and reaches the same DESeq2 negative-
binomial-GLM-plus-apeGLM-shrinkage estimator family natively in Python).
Nothing has been installed; this is a plan only. Per
`review:20261014_203915` F3, the plan now proposes a NEW isolated
environment path (`302006_de_backend_v1`, never `302006_env_clean`) and
Python 3.11 + the current `pydeseq2==0.5.4` as the primary pin (the
isolated env's own Python removes the reason to stay on `0.5.2`), with
`0.5.2`/Python 3.10 kept as a documented fallback; a real (read-only)
`pip install --dry-run`/`conda create --dry-run` dependency resolution was
executed for both options (see the document for the important caveat on
what that dry-run does and does not prove), and `lfc_shrink`'s real
signature (`coeff` has no usable default, despite its own docstring)
was verified directly against PyDESeq2's source, not its rendered docs.
**Revised per `review:20261014_221155` R4**: re-ran the dry-run
with `--ignore-installed` (a real, full, 40-package dependency closure,
not reusing `302006_env_clean`'s already-installed packages), re-verified
`lfc_shrink`/`DeseqDataSet` directly against `v0.5.4`'s own source (not
only `v0.5.2`'s), and changed the deployment draft to an explicit
`--prefix <path>` rather than only `-n <name>`.

**New this round**: `code/de_design_providers.py` -- pseudobulk-to-DE
input/design-matrix PREPARATION only (per `review:20261014_151747`
step 2, the "B" increment; fixed per `review:20261014_203915`
F1/F2 below). `prepare()` takes a `pseudobulk_result_path` and calls
`read_pseudobulk_result`/`load_pseudobulk_counts` ITSELF -- there is no
longer a caller-supplied `counts` array independent of the manifest (F1's
fix; the review's probe exploited exactly that parameter) -- and the SAME
`ComparisonSpec`/`SampleStructure` interfaces `ExpressionDERankingProvider`
uses; supports exactly two synthetic designs
(`independent_donors_unpaired`, `paired_or_repeated_within_donor`,
matching `gene_ranking_providers._VALID_DESIGN_TYPES`); builds a real
numeric design matrix, a named target contrast, rank, and residual
degrees of freedom -- full rank alone is never treated as adequate power;
a design with fewer than 1 residual degree of freedom is refused even
when full rank. Declared covariates, repeated samples within the same
donor/group, partial pairing, and a straddling donor under the unpaired
design are all explicitly refused, never silently pooled. Counts/metadata
are aligned by `sample_id` LOOKUP, never by position. A declared
`comparison`/`sample_structure` column name that does not match what the
pseudobulk aggregation ACTUALLY used (`params`) is refused before any
design is built (F2's fix); at read time, every semantic field (`design_
matrix`/`column_names`/`formula`/`contrast_vector`/`contrast_direction`)
is REBUILT from the persisted `sample_ids`/`donor_ids`/`group_values`/
`comparison`/`sample_structure` and compared, not merely hash-checked plus
an "estimable" check -- the latter is VACUOUS for any full-rank design
(`c @ pinv(X) @ X == c` identically for ANY `c`), so a reversed or zero
contrast previously passed. `PreparedDEDesign` also now carries the
upstream `dataset_identity` (species/gene namespace/counts semantics)
forward, so a downstream consumer never has to guess them from a file
name. Output status is `'prepared'`/`'failed'` ONLY -- deliberately
distinct from `ExpressionDEResult`'s `'complete'`/`'failed'`/`'not_run'`,
so "a design matrix was prepared" can never be confused with "a DE
ranking was computed". `ExpressionDERankingProvider` itself is UNCHANGED
by this -- still no real backend, still never returns a fabricated
`'complete'` result. See `examples/run_de_design_example.py` for a
runnable human-unpaired + mouse-paired walkthrough through the real
aggregate -> prepare -> write -> read/load chain.

**R1 fix** (`review:20261014_221155`): the F1/F2 fixes above only closed the `prepare()`
entry point -- the public `write_prepared_de_design(out_dir, result, counts=...)` writer and the
public `load_prepared_de_counts(manifest)` reader could still be called directly with an
internally-self-consistent but WRONG `PreparedDEDesign`/`counts` pair (swapped species/gene
namespace/donor IDs, a stale upstream hash, a nonexistent upstream path, or same-shape reversed/
negative counts handed straight to the writer) and both silently accepted it, because neither
re-resolved the claimed upstream pseudobulk source -- only checked the payload's own internal
self-consistency. Fixed by a new `_source_binding_problems` check, shared by both the write and
the read/load paths, that re-reads the claimed `pseudobulk_manifest_path` LIVE via
`pseudobulk_providers.read_pseudobulk_result`/`load_pseudobulk_counts`, and refuses any
mismatch against that live source's own `dataset_identity`, column-name `params`, per-sample
`donor_or_animal_id`/`condition_value`, `gene_ids`, or counts content -- plus a standalone
non-negative/finite/integer-valued counts contract. All eight of the review's own `probe_binding.py`
reproduction cases (re-run VERBATIM, unmodified, against the fixed code) are now rejected.

**The offline GSEA backend's numeric chain has been run end to end and
replayed from frozen artifacts, on synthetic data only** — this is evidence
that the chain computes and reproduces correctly, not a frozen evidence
package that has been reviewed and accepted for any real dataset or
biological claim. Per the owner's current one-month-MVP framing (see
`review:20261012_172156`), the immediate next task is closing out
this offline evidence chain, not building the production adapter, a live
GO/KEGG fetch, or Docker packaging — those remain explicitly out of scope
until the owner has reviewed this round's evidence and made the backbone/
comparison-design decisions that gate them.

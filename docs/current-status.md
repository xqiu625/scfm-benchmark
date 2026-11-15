# Current Status — 302006

**Research scope updated 2026-10-12 from the owner's explicit decision.** Operational records
below retain their individual verification dates; HPCC state was not rechecked for this scope update.
This file holds *only* the current state. Update by overwriting, never by
appending. History belongs in `docs/experiment-log.md`; decisions in `docs/decisions.md`; the
plan in `docs/execution-plan.md`.

Replaces `docs/memory-bank/{MEMORY,activecontext,progress}.md`, archived 2026-09-21 under
`archive/superseded_2026-09-21/` — those files must not be read as current state.

---

## What the project is

**Owner decision, 2026-10-15:** the primary within-COVID monocyte contrast is now explicitly
**(Severe + Critical) minus (Mild + Moderate)**. Critical pooling is confirmed; older pending
wording below is superseded by this decision. The owner also requests a later **sepsis ICU vs
COVID-19 ICU** comparison. Actual ICU location at sampling must be tracked separately from
severity. **The owner has now explicitly selected ICU-only sampling for the current severe/critical
arm:** require Severe/Critical AND confirmed ICU location at that sample's collection time.
The reference arm remains Mild/Moderate. Unknown location does not meet confirmed ICU eligibility;
retain unresolved records for audit. No new analysis has run on this basis. See `docs/decisions.md`.

**A reusable scFM biological-discovery framework:** use pretrained representations to identify
cell states and disease-associated gene programs, connect gene-level evidence to GSEA and
GO/KEGG annotations, and produce traceable pathway insights. Deliver through a Docker-based
framework with explicit model and species adapters.

The owner has replaced the benchmark/mortality-ranking objective. Quality controls and expression
baselines remain necessary to validate discoveries, but are not a separate benchmark workstream.
The previous fixed five-model set and target venue are not new-framework commitments.

**Model scope (owner clarification, 2026-10-12): one primary scFM for the framework.** Do not
implement every available model or make cross-model consensus a required output. Select the
backbone using pathway-interpretation suitability, human/mouse applicability, reproducibility and
deployment feasibility. Its exact checkpoint is not yet selected. A small, targeted comparison can
support the paper's choice of backbone; superiority is a result to test, not a premise. Distinguish
backbone evidence from improvements due to the surrounding framework.

**Limited shortlist (2026-10-12, corrected twice same day — `review:20261012_145004`,
then `153316.md`):** `scripts/34_scfm_discovery_framework/MODEL_SELECTION_v1.md` recommends
**scPRINT-2 as a candidate worth a minimal validation first, UCE as the alternative/control** — a
recommendation, not a selection, and not a claim that either is already proven or working on mouse
data. Geneformer/scGPT/STATE are structurally ruled out as the cross-species backbone (no official
mouse-trained checkpoint found for any of the three); their existing human embeddings are historical
artifacts against the old cohort/QC version, not validated against the new monocyte input.
**scPRINT-2's `small-v2.ckpt` checkpoint's own saved metadata (not README prose) declares human and
mouse among 16 organisms** — independently verified this round by two different parties via a safe
byte-level scan of the checkpoint's pickle (no deserialization): this is `checkpoint_metadata_confirmed`
species-*coverage-list* evidence specifically, not a per-species training-cell count, train/val
split, or any downstream efficacy claim. UCE's species evidence remains
`documented_unconfirmed` (the review reports reading the UCE Nature paper directly and finding
human+mouse named; this has not yet been independently re-verified against the primary source).
**Order:** framework first; human **monocytes only** as the first real-data application.
Existing platelet results remain preserved. Patient identity, rather than individual cells, defines biological
replication; specific contrasts and validation designs remain to be agreed before real analysis.

---

## Where the work actually stands

The framework direction and scope are owner-confirmed. **Framework foundation implemented
2026-10-12, corrected twice the same day** (`scripts/34_scfm_discovery_framework/`):
`review:20261012_153316` fixed `validate` accepting type-invalid config fields and a
matrix with donor/sample/batch columns that did not exist in the file, `plan`'s provenance not
being independently reviewable, `read_manifest` only checking the schema-version string, and the
gene-ranking validator silently accepting NaN/empty rankings.
`review:20261012_163250` then found and fixed 3 further gaps: `raw.X` was checked against
the wrong gene axis (the current, possibly-filtered `var` instead of `raw.var` — a legitimate
`raw.X` with MORE genes than the current `var` was wrongly rejected) and an unresolvable obs/var
index was silently skipped rather than failed; a relative `expression_matrix_path` resolved against
the process's cwd instead of the config file's own directory (the same config gave different
results depending on where it was run from); and the manifest reader/writer type contracts were not
fully aligned (a boolean `source_code_revision_or_hash`, a `null` dependency entry, and a hand-raw
NaN written to disk were all accepted on read). All fixed and re-verified on HPCC — **91 tests pass**
(up from 43 → 78 → 91 across the three rounds), including out-of-band mutation-kill checks that
disable the structural-check / gene-axis / path-resolution logic in scratch mutants and confirm the
targeted tests genuinely fail. A real, runnable example (`examples/run_example.py`) now builds a
tiny synthetic fixture and runs the full `validate --check-matrix-structure -> plan
--check-matrix-structure -> read_manifest` chain end to end, asserting the read-back object exactly
equals what was written.

Two more small-scope gaps found and fixed same day (`review:20261012_171916`): a numeric
literal that overflows during JSON parsing (e.g. `1e999`) produced `inf` and was accepted on read,
an asymmetry with `write_atomic`'s `allow_nan=False` — fixed with a recursive non-finite-float check
applied at build/write/read; and `validate`'s control flow called the path resolver before confirming
the config was type-valid, crashing with a raw `TypeError` on a non-string `expression_matrix_path`
instead of a structured report — fixed by matching `plan`'s existing fail-fast-on-invalid-config
order. **98 tests pass** after those two fixes (up from 91).

**Local gene-set cache + offline GSEA implemented** (`code/gene_set_cache.py`, `code/gsea_runner.py`,
2026-10-12, same day): imports a local `.gmt` into a versioned, hash-verified cache (never a URL,
never a fabricated live-retrieval timestamp); `LocalPrerankGSEAProvider` runs the installed
`gseapy.prerank` against that cache and the existing `GeneRanking`/`GSEAProvider` contracts, with
species/namespace consistency checks (fails closed, no silent remap) and a `no_overlap` state kept
distinct from a completed run with no significant terms. A second runnable example
(`examples/run_gsea_example.py`) builds a tiny synthetic gene set with a known enrichment direction,
runs the real backend, and verifies the result — not a mocked call. **131 tests pass total.** Model
choice (scPRINT-2 vs. UCE) still does not block this work; it only gates the specific production
*adapter*. See `review:20261012_140650`/`145004.md`/`153316.md`/`163250.md`/`171916.md` in
the parent project directory for the model-scope clarification and its corrections; earlier design
handoffs apply only where consistent with this single-backbone scope.

**GSEA chain's identity/provenance gaps found and fixed** (`review:20261012_220407`,
2026-10-12, same day, F1-F5): the installed gseapy's own `base.py` was found (by reading its real
source, not assuming behaviour) to silently uppercase gene IDs for matching when a ranking isn't
mostly-uppercase but the gene sets sampled are, producing a self-contradictory result (a 1-gene
case-sensitive intersection reported alongside a 12-gene case-folded leading edge) — fixed by
pre-filtering every gene set to its exact, case-sensitive ranking intersection before gseapy ever
sees it, plus a defense-in-depth check on the backend's own reported leading edge (F1). A ranking's
hash alone verified content but could not recover it, and no replay command existed — fixed by
freezing the full ranking (exact gene/score order, tied-score count) to `RANKING_SNAPSHOT.json` on
every run, binding the wrapper's source revision and dependency versions into the result, and adding
`code/gsea_replay.py`, which replays a run from its frozen artifacts alone in a brand-new process and
compares every scientific field (F2). The GSEA result reader accepted a payload with only
`{schema, status="complete"}`, an impossible `p_value=2.0`, and a `params` field that overflowed to
`inf` — fixed with one `validate_gsea_result()` shared by the writer and reader, covering the whole
payload (F3). The gene-set cache stored its own file path relative to the import-time cwd, breaking
a cache loaded from a different directory or copied elsewhere — fixed by resolving it relative to
the manifest's own directory, with a cache-manifest schema validator shared by import and load (F4).
When every gene set was filtered out by size after a real overlap existed, the backend raised an
unhandled `LookupError` and no result was ever written — fixed by replicating the size-acceptance
check before calling the backend, recording each excluded gene set's reason, and (for any other
genuinely unexpected backend exception) writing an audited `failed` record before re-raising rather
than swallowing it (F5). **153 tests pass total** (up from 131), including out-of-band mutation-kill
checks (reverting the F4 path-resolution fix and the F5 size-filter precheck in scratch mutants) that
confirm the new regression tests genuinely fail without the fix. The offline example now also
replays its own run from frozen artifacts in a new process and confirms an exact match. **This is
the backend numeric chain running and reproducing correctly on synthetic data — it is evidence
toward, not yet, a frozen offline evidence package reviewed and accepted for any real dataset.**

**Replay-identity and failure-ordering gaps found and fixed** (`review:20261013_114111`,
2026-10-13, R1/R2/C1): an independent review of the above round found that `gsea_replay.py`'s
"REPLAY MATCH" only proved the numbers agreed, never that the replay used the SAME ranking, cache,
or source/dependency identity as the original run — a legitimately re-hashed but DIFFERENT snapshot
(same scores, different `ranking_definition`), a swapped cache with a different collection/version,
and a tampered source-code-revision field were all accepted as a match (R1). Fixed by having replay
compare the reconstructed ranking's semantic hash and the snapshot file's own hash against what the
ORIGINAL result itself recorded (not just the snapshot's internal self-check); by recording the
cache's full semantic identity (species/namespace/source_version, not only collection name/content
hash) in the result and refusing a `--cache-dir` whose cache doesn't match it (a verbatim-copied,
relocated cache still passes); by requiring every numeric run parameter present with the correct
type (`validate_gsea_result` now rejects a `params` dict missing a required key, or a boolean where
an identity string belongs — a probe found `source_code_revision_or_hash=True` and
`params={"unrelated": 1}` both previously accepted); and by treating a source-code/dependency-
version drift as a strict-replay failure by default, with a separately labeled `--allow-identity-
drift` ("REPLAY COMPARISON", never "REPLAY MATCH") for the deliberate cross-version-numbers-only
case. `SOURCE_MODULES` was also missing `artifact_schema.py`, which `gsea_runner.py` imports and
calls on every run — added. A NaN/Inf score made `write_ranking_snapshot()` raise mid-write
(`allow_nan=False`) BEFORE `ranking.validate()` ever ran, leaving a bare `.inprogress` file and no
`GSEA_RESULT.json` at all — fixed by validating first; an invalid ranking now gets a clean
`status="failed"` record with no snapshot and no leftover temp file (R2). A case COLLISION within a
ranking's own distinct gene IDs (e.g. `ANCHOR` and `Anchor` both present) can still trigger gseapy's
internal case-folding even after the F1 pre-filter, which only eliminates CROSS-set leakage — this
stays correctly fail-closed via the existing defense-in-depth leading-edge check (confirmed by a
direct regression against the real backend), but the module docstring's "nothing left to expand"
phrasing overclaimed the pre-filter's own scope and has been corrected (C1). **165 tests pass total**
(up from 153, 12 new), including two genuine out-of-band scratch mutants (reverting the new replay
ranking-identity check, and reverting the validate-before-snapshot ordering) that reproduce the
original failure modes exactly (a silent pass for the former, the literal pre-fix
`ValueError: Out of range float values are not JSON compliant: nan` for the latter) when deleted,
plus inline monkeypatch mutation-kills for the new cache-identity check. **Still the same standing
caveat as above: this closes software evidence-chain gaps on synthetic data; it is not a claim about
any real dataset or biological finding.**

**Replay-status scoping, replay-label fix, and the first PCA/gene-ranking-interface round**
(`review:20261013_135634`, 2026-10-13, R3/R4 + PCA/ranking interfaces): an independent
re-review of commit `131a794` found two further `gsea_replay.py` gaps and specified the next-phase
PCA baseline / gene-ranking interfaces. **R3**: the cache-identity check only ran for a
`status="complete"` original; a probe found a `"no_overlap"` original replayed against a cache
swapped for a different collection/version still printed "REPLAY MATCH" (both are vacuously
zero-overlap, so the numeric fields agreed even though the cache differed). Fixed by restricting
strict replay to `status="complete"` originals only — any other status is refused outright
(`SystemExit`, non-zero exit) before the snapshot is read, the cache is checked, or the backend is
ever called; the original record is left untouched. **R4**: the printed "REPLAY MATCH"/"REPLAY
COMPARISON" label depended on whether `identity_drifts` was non-empty (i.e. whether a drift was
*observed*), not on whether `--allow-identity-drift` was actually passed — so that flag against two
identical environments still printed "REPLAY MATCH", contradicting the CLI's own help text's promise
("REPLAY COMPARISON, never REPLAY MATCH, even when every number agrees"). Fixed by keying the label
on `allow_identity_drift` itself. **205 tests pass total** (up from 165: 7 new replay tests + 33 new
PCA/ranking-interface tests), including a real out-of-band HPCC scratch mutant (not merely an inline
monkeypatch) that reverted a genuinely NEW bug found and fixed this round (see below) and reproduced
the exact original crash.

This round also implements the first piece of the post-GSEA work the owner's one-month MVP calls
for: `code/embedding_providers.PCAEmbeddingProvider`, a REAL `scikit-learn`-backed PCA baseline — an
explicit CONTROL, never evidence that a scFM is superior — that reuses `dataset_schema`/`cli.py`'s
own input validation (never a bare file path), reads the declared `counts_slot` read-only (dense or
`csr_matrix`/`csc_matrix`, bounded to an explicit `_MAX_DENSE_ELEMENTS` this round), and persists a
hash-checked, round-trip-readable manifest + npz artifact. Label-free feature selection only (all
genes, an explicit caller-supplied list, or top-variance `hvg_top_n` — no disease/outcome-label
parameter exists on this interface). **A real bug was found and fixed during this round's own
testing** (not by automated review — by the implementer's own test run): a structural problem the new readers
detect (an oversized block, an unrecognized sparse encoding) raised a bare, uncaught `SystemExit`
instead of the established "every failed run leaves a clean, audited record" pattern this project's
GSEA chain already follows (the same class of defect R2 fixed for `write_ranking_snapshot`) — fixed
by catching it and converting it into a `status="failed"` `EmbeddingResult`; independently confirmed
via a genuine out-of-band HPCC scratch mutant that reverting the fix reproduces the exact original
uncaught-crash behavior. `code/gene_ranking_providers.py` implements the two gene-ranking EVIDENCE
STREAM contracts the review's section 3 specified, correcting the PRIOR handoff's own overclaim
(masking framed as a cell-level mask, the model stream framed loosely as "embedding-dimension
attribution"): `ExpressionDERankingProvider` (an explicit biological comparison; a
`check_design_estimability` guard refuses a donor-blocked design where every donor contributes to
only one compared group) and `ModelMaskingRankingProvider` (`Delta(g) = s_full - s_masked(g)` for a
FIXED, bound scoring function — masking a GENE INPUT, never a cell-selection mask; the interface has
no cell-selection parameter at all). **Both are type/validation/serialization contracts only this
round** — `run()` always returns `status="failed"` or `status="not_run"`, never a fabricated
`"complete"` result; no main model checkpoint is pinned, and no donor/animal grouping, masking score,
or comparison has been confirmed by the owner (open decision table in this round's
implementation handoff). A small, explicitly-synthetic human + mouse PCA example
(`examples/run_pca_example.py`) runs the same interface end to end for both species. **Same standing
caveat: this is software-contract and offline-baseline work on synthetic data — it is not a real
scFM/PCA comparison, not a real gene-expression or model-masking finding, and not a claim about any real dataset.**

**F1-F4 design-contract and evidence-identity fixes** (`review:20261013_210131`, 2026-10-13):
an independent review of commit `14a5037` ran a bounded synthetic probe against the new PCA/
ranking-interface round and found four gaps, none touching the already-accepted GSEA R3/R4 work.
**F1**: `check_design_estimability` applied a donor-blocked `~ donor + condition` check
unconditionally — three independent healthy donors vs. three independent disease donors (an
ordinary unpaired design) was wrongly refused, while a single cross-group donor out of several was
wrongly enough to call a full donor-blocked design "estimable." Fixed by adding an EXPLICIT,
code-recognized `SampleStructure.design_type` (`"independent_donors_unpaired"` /
`"paired_or_repeated_within_donor"`, never inferred from an `estimator=` string) that decides which
check runs, each checked against the REAL donor/group input; `comparison.group_column` disagreeing
with `sample_structure.group_column` is now refused; declared-but-unchecked `covariate_columns` are
flagged UNVERIFIED in the reason. **F2**: a source gene axis with a duplicate ID (e.g. two `"DUP"`
columns with different values) was silently resolved by `{g: i for i, g in enumerate(gene_ids)}`
onto whichever column came LAST, with no indication the other was dropped. Fixed by validating the
gene axis (non-empty, string, UNIQUE) before any ID->position map/feature-selection branch is built
— a duplicate is now a clean `failed` record, never auto-merged/suffixed/averaged. **F3**: an
all-zero (zero-total-variance) counts block produced a NaN `explained_variance_ratio_`, which
`_write_embedding_result` wrote `embedding_matrices.npz` for BEFORE the JSON manifest's own
validation raised — leaving an orphaned npz with no `EMBEDDING_RESULT.json`. Fixed by checking
EVERY array about to be persisted (not just embedding/components/mean) for finiteness before
anything is written, plus an explicit up-front zero-total-variance refusal; an unexpected backend
exception is now recorded as a `failed` result first, then re-raised. **F4**: `ModelMaskingResult`
for `perturbable_genes=("A",)` vs. `("B",)` serialized to an IDENTICAL record (the gene set was
validated but never stored). Fixed by binding each result's ACTUAL input content identity
(`perturbable_genes` + its hash; `group_values_sha256`/`donor_values_sha256`) — `GeneRanking`
itself (shared with the already-tested GSEA chain) is deliberately left unchanged, and
`make_evidence_ranking`'s substring check is now documented as a lint-level tripwire, not a
structured provenance check. **226 tests pass total** (up from 205: 21 new F1-F4 tests), with
dedicated mutation-kill tests for every fix, including a real HPCC-observed `RuntimeWarning:
invalid value encountered in divide` confirming the zero-variance mutation test genuinely drove
sklearn into the original NaN path. The standalone human+mouse PCA example was re-run end to end
and produces numerically identical output to the prior round (`explained_variance_ratio`
unchanged for both species), confirming the normal-input path is unaffected. **Same standing
caveat: software-contract fixes on synthetic data — no real dataset, no real DE/model backend, no
new scientific claim.**

**R1-R3 entry-point/evidence-envelope fixes and a first pseudobulk-aggregation increment**
(`review:20261013_222450`, 2026-10-13): an independent review of commit `0fec319`
found the F1-F4 round above only partially closed F1/F4, and specified a first real increment
for `ExpressionDERankingProvider`'s backend. **R1**: `ExpressionDERankingProvider.run`'s entry
point only checked that SOME compared-group sample existed, never that BOTH `reference_level`
and `comparison_level` were present (a `healthy`-only input passed as `not_run`); an empty-string
or `None` donor/group value reached Python's `sorted()` uncaught, raising a bare `TypeError`
instead of a structured `failed` record. Fixed by a new `_validate_group_donor_values` check
(length match, non-empty string values, BOTH levels present) that runs inside
`check_design_estimability` before any donor grouping or `sorted()` call. **R2**: the
`paired_or_repeated_within_donor` refusal message claimed a partially-paired design "cannot be
estimated", which is not generally true (demonstrated with the review's own counterexample: a
6-row/5-column `~ donor + condition` design matrix with one donor pair fully paired and two
donors contributing to only one level is full rank). Fixed by rewording the refusal as a
SUPPORT-SCOPE limit ("UNSUPPORTED by this round's check", not mathematically inestimable) and
stating explicitly that straddling/non-straddling donors are never silently dropped. **R3**:
`ExpressionDEResult`/`ModelMaskingResult` carried only a length+hash summary, with no way to
recover the actual ordered input from the record, no dedicated writer/reader, and an identity
hash built with `json.dumps(..., default=str)` (silently stringifies an unrecognized value
instead of refusing it). Fixed by adding an `input_mapping`/`perturbable_genes` round-trippable
field plus a recomputed-and-compared content hash to both results, new
`write_expression_de_result`/`read_expression_de_result` and
`write_model_masking_result`/`read_model_masking_result` functions (validated identically at
write and read, mirroring `embedding_providers.write/read_embedding_result`), and a
`_safe_json_identity` helper that returns an explicit "identity unavailable" state instead of
ever silently stringifying a non-JSON-native value. Also fixed the three outstanding
documentation-accuracy items from the prior round (`input_counts_sha256`'s actual log1p/
standardize phase, the unsaved standardize mean/std, and the dense-read bound's refusal message
no longer suggesting `cell_ids_subset` as a workaround it is not). **A new
`pseudobulk_providers.py` module** implements the first real backend increment: sample x
(caller-selected) cell-type/state RAW count-sum aggregation (`dataset_config["sample_id_column"]`
is always the aggregation unit, never donor/animal ID, so repeat visits never silently merge),
with explicit sample/donor/condition/visit identity, conflicting per-sample metadata refused,
a counts-conservation check, and the same strict write/read/tamper-rejection contract as the
other providers. It does **not** compute any DE statistic, p-value, log-fold-change, or
`GeneRanking` — `ExpressionDERankingProvider` itself is untouched. **266 tests pass total** (up
from 226: 18 new gene-ranking-provider tests, 21 new pseudobulk tests), with dedicated
mutation-kill tests for every new guard. The standalone human+mouse PCA example and the GSEA
example were both re-run end to end on HPCC and still pass. **Same standing caveat: software-
contract fixes and a synthetic-data aggregation increment — no real dataset, no real DE
statistic, no new scientific claim.**

**F1-F3 pseudobulk-identity/type-validation fixes** (`review:20261014_100436`,
2026-10-14): an independent review of commit `32b22fa` found the pseudobulk increment's identity
binding and the evidence-envelope's status contract were only partially closed. **F1**: the
`PSEUDOBULK_RESULT.json` manifest's `sample_rows`/`gene_ids` LABELS carried no binding to the
npz array's actual row/column CONTENT — reversing `sample_rows` or `gene_ids`, or editing
`n_distinct_donors`/a row's `library_size`/`condition_value`, passed both the writer's and the
reader's validation unchanged (the npz file-hash check proves only that the npz bytes are
unmodified, never that the JSON's labels for those bytes are correct). Fixed with two
independent mechanisms: a new `identity_sha256` field (a hash over a canonical snapshot of
`gene_ids`/`gene_total_counts`/`sample_rows`/the derived counts, recomputed and compared on every
read) catches any field changing without every dependent field being rewritten consistently;
`load_pseudobulk_counts` now additionally recomputes each row's library size and each gene's
total count DIRECTLY from the loaded npz array and compares them against the manifest's claims,
binding the JSON labels to the array's real content rather than only to the JSON's own
self-consistency. `dataset_identity` also gained `selected_raw_counts_sha256`/
`selected_cell_metadata_sha256`/`gene_ids_sha256` — content hashes of the exact selected-cell raw
counts block and per-cell metadata actually read this run (not a re-read of the whole source
file; a config hash alone is not a content hash of the expression matrix). **F2**: type/null
validation ran too late or not at all in three ways — `cell_type_or_state_label=True` passed the
entry truthiness check and ran to npz-write time before the writer's own type check rejected it,
leaving an orphan npz with no manifest; a sample whose selected cells all had a missing
`condition`/`visit` value (the column explicitly requested, not merely unprovided) silently wrote
`condition_value: null`; a sample mixing a missing value with a real string value crashed
`_sample_metadata_conflict_problem`'s `sorted()` call with an uncaught `TypeError`. Fixed by
validating `cell_type_or_state_label`/`cell_ids_subset`/`condition_column`/
`visit_or_state_column`'s TYPES at entry, refusing (never raising) any selected cell whose
explicitly-requested condition/visit value is missing/non-string before any grouping or
`sorted()` call, and splitting `_write_pseudobulk_result`'s validation into a pre-npz core check
(everything except `aggregation_path`/`aggregation_sha256`) as a second, independent guard
against an orphan npz. **F3**: the evidence envelope's write/read validators were asymmetric — an
illegitimate `status='complete'` record (`ranking={}`, or `backend='not_implemented'`) was
accepted, while a genuinely `'failed'` record carrying a diagnostic `None` donor_value (the REAL
provider's own record of that invalid input) was REJECTED by the same strict checks a
`'complete'`/`'not_run'` record needs. Fixed together: `ranking` is now re-validated against
`GeneRanking.validate()` itself (an empty dict is no longer "non-null"), `backend='not_implemented'`/
`backend_version='n/a'` are refused on a `'complete'` record, `MaskingScoreSpec.validate()` now
checks TYPE not truthiness (`score_definition=True` is truthy but not a string), and
`validate_expression_de_result`/`validate_model_masking_result` now gate their strict per-row/
per-field VALUE checks to `status in ('complete', 'not_run')` — a `'failed'` record's
`input_mapping`/`masking_score` is a diagnostic snapshot of the actual (possibly invalid) input,
never a record that claims to BE a validated one; `n_samples` vs. `len(input_mapping)` and nested
`sample_structure.design_type`/`comparison` field checks now run on every read, not only at
construction. The shared `_safe_json_identity` helper moved to `artifact_schema.py`
(`safe_json_identity`) so `pseudobulk_providers.py` can reuse it without a backwards layering
dependency. A new standalone `run_pseudobulk_example.py` demonstrates the real interface end to
end for both species, including the failed-record read-back path. **288 tests pass total** (up
from 266: 13 new pseudobulk tests, 9 new gene-ranking tests), with dedicated mutation-kill tests
for the new `identity_sha256`/pre-npz-validation/ranking-shape guards. A read-only DE-backend
dependency check (per the review's instruction) confirmed independently: the clean conda env is
Python 3.10.18 with neither `pydeseq2` nor `rpy2` installed, and this shell's `Rscript` was not
found on `PATH` — this does NOT establish that no R/DESeq2 path exists on the cluster, only that
this one non-interactive shell did not find one. **Same standing caveat: software-contract fixes
on synthetic data — no real dataset, no real DE/model backend, no new scientific claim.**

**R1-R3 equal-marginals/type-check/failed-record fixes** (`review:20261014_114357`,
2026-10-14): an independent review of commit `0a40d4f` reproduced its own probe and found F1/F2/F3
were only PARTIALLY closed. **R1**: the F1 npz-identity binding was itself insufficient against an
EQUAL-MARGINALS swap — a fixture with two rows whose library sizes AND per-gene totals are
identical (`[[1, 9], [9, 1]]`, both summing to 10) could have `sample_rows` reversed, `gene_ids`
reversed, or a `condition_value` rewritten, with every recomputed total (and a self-consistently
rehashed `identity_sha256`) still agreeing — a total is a many-to-one projection of content, never
a content identity. Fixed by saving `sample_ids`/`gene_ids`/a normalized per-row identity snapshot
(donor/condition/visit/cell_ids) INSIDE the npz itself (plain numpy string arrays, no pickle,
covered by the same whole-file `aggregation_sha256`); `load_pseudobulk_counts` now compares the
manifest's claimed order/metadata DIRECTLY against this independent, hash-protected in-npz copy.
**R2**: `n_samples` (expression-DE) and `n_pseudobulk_samples`/`n_distinct_donors` (pseudobulk) had
the same `isinstance(x, int) and ... != ...` shape with no `else` branch, so a non-int value
silently SKIPPED the consistency check instead of being rejected for its type — now checked
unconditionally first in all three places. **R3**: the `n_samples == len(input_mapping)` check ran
UNCONDITIONALLY, which made a genuinely `'failed'` length-mismatch record (`group_values`/
`donor_values` of different lengths, so `input_mapping` is correctly forced empty while `n_samples`
stays the raw non-zero `group_values` length) unwritable; `perturbable_genes`' per-entry value
check had the matching gap for `ModelMaskingResult`. Fixed by new `input_lengths`/`input_mapping_
buildable` fields (recorded ALWAYS, never dropped on failure) that the length-consistency check is
now gated on, and by gating `perturbable_genes`' value check the same way `masking_score`'s fields
already were. **300 tests pass total** (up from 288: 6 new pseudobulk tests, 6 new gene-ranking
tests), with the 11 new behavioral tests independently confirmed to fail against the pre-fix
(`0a40d4f`) code for the exact reason each fix addresses, plus a dedicated mutation-kill test for
the new npz-identity-binding guard. **Same standing caveat: software-contract fixes on synthetic
data — no real dataset, no real DE/model backend, no new scientific claim.**

**P2 status-contract fix, DE design-matrix-preparation increment ("B"), and DE-backend
environment plan** (`review:20261014_151747`, 2026-10-14): an independent review of
commit `3131146` confirmed R1/R2/R3 above now closed (independently re-verified the equal-
marginals rejections, the non-int type rejections, and the two genuinely-`'failed'` round-trips),
and found one new P2 gap: `input_mapping_buildable=False` was cross-checked against
`input_mapping`'s shape but never tied to `status` itself, so a genuinely `'failed'`
length-mismatch record could have ONLY its `status` rewritten to `'not_run'`/`'complete'` on disk
and still be accepted. Fixed by rejecting `input_mapping_buildable=False` combined with
`status in ('not_run', 'complete')` directly (this module defines `'not_run'` as "inputs passed
validation, no backend ran yet", which an unbuildable input never satisfies) — independently
confirmed via a genuine before/after run of the 3 new tests against the unmodified `3131146` code
(all 3 fail for the exact reason the fix addresses).

**New `code/de_design_providers.py`** (the "B" increment, gated on the above passing, per the
review's own A/B instruction): pseudobulk-to-DE input/design-matrix PREPARATION only, no
statistical fitting. Consumes an ALREADY-validated pseudobulk aggregation
(`read_pseudobulk_result` + `load_pseudobulk_counts`) and the SAME `ComparisonSpec`/
`SampleStructure` interfaces `ExpressionDERankingProvider` uses (never a duplicate pair of
types); supports exactly two synthetic designs (`independent_donors_unpaired`,
`paired_or_repeated_within_donor`); builds a real numeric design matrix, column names, formula,
encoding convention, a named target contrast (with direction "comparison minus reference"), rank,
and residual degrees of freedom. Declared covariates, a repeated sample within the same
donor/group, partial pairing, a straddling donor under the unpaired design, a missing comparison
level, a non-full-rank design, and a design with fewer than 1 residual degree of freedom (full
rank is explicitly never treated as adequate power) are all refused, never silently pooled or
accepted. Counts/metadata are aligned by `sample_id` LOOKUP, never by position (tested by
reversing the pseudobulk manifest's own row order and confirming the design still aligns
correctly). Persistence mirrors `pseudobulk_providers.py`'s own R1 fix: the selected/reordered
counts submatrix is saved to its own npz with its OWN stored `sample_ids`/`gene_ids`, covered by
both a whole-file hash and a content-only hash, so a manifest-only sample-id reordering (the same
class of bug R1 fixed) is caught directly — proven via a dedicated mutation-kill test. Output
status is `'prepared'`/`'failed'` ONLY, deliberately distinct from `ExpressionDEResult`'s
`'complete'`/`'failed'`/`'not_run'`, so "a design matrix was prepared" is never confused with "a
DE ranking was computed"; `ExpressionDERankingProvider` itself is unchanged — still no real
backend, still never returns a fabricated `'complete'` result. 17 new tests, all using
hand-computable synthetic fixtures (paired and unpaired, with hand-verified rank/residual-df/
contrast values).

**New `scripts/34_scfm_discovery_framework/DE_BACKEND_ENVIRONMENT_PLAN_v1.md`** (step 3, a plan
only — nothing installed/upgraded this round): recommends pinning `pydeseq2==0.5.2` (NOT the
PyPI-latest `0.5.4`, which requires Python >=3.11 — PyPI's own metadata, fetched 2026-10-14, shows
`0.5.2` is the last release still declaring `requires_python>=3.10.0`, matching the clean env's
Python 3.10.18). Also corrects a prior round's framing: that round's "R/module not found" was from
a single non-interactive SSH shell and was correctly left unconfirmed either way; this round
checked an interactive LOGIN shell's module system directly and found R genuinely IS available on
this cluster (`R/4.1.2` through `R/4.6.0`), so an R/DESeq2/`rpy2` path remains real, just not the
primary recommendation (PyDESeq2 avoids a second language boundary). Documents a draft
`pip install --dry-run` command (not run), an LFC-shrinkage plan (PyDESeq2's `lfc_shrink`, an
apeGLM-prior method, recording the shrunken LFC alongside the raw Wald LFC/p-value, never
overwriting it), and a positive/null synthetic verification plan (an injected-effect positive
control; a null control whose calibration is checked, never requiring a single draw's LFC to
equal exactly 0; a paired-design null that shuffles labels WITHIN each donor's pair, never across
donors). **Full test suite: 320 tests pass** (up from 300: 3 new gene-ranking tests for P2, 17 new
de-design tests — correcting an earlier "64 → 70 / +6" mislabel; the actual counts are
`67 → 70` gene-ranking and `0 → 17` de-design, `3 + 17 = 20` new, `300 + 20 = 320`). **Same
standing caveat: software-contract/design-matrix work on synthetic data — no real dataset, no
real DE statistical backend, no new scientific claim; the owner's pending decisions (real
comparison, donor/visit structure, main model/checkpoint, masking score, formal QC/validation
design, and whether/when to actually install a DE backend) remain unaffected.**

**F1/F2/F3 fixes to the "B" increment, a new DE-design standalone example, and an isolated
DE-backend environment plan** (`review:20261014_203915`, 2026-10-14): an independent
review of commit `51d6b75` found two new P1 gaps in `de_design_providers.py` and one P2 gap in
`DE_BACKEND_ENVIRONMENT_PLAN_v1.md`. **F1** — `prepare()` took an independently-supplied `counts`
ndarray alongside the manifest `dict`, checking only `status`/`shape` agreement; same-shape
reversed/negated counts, or a `sample_rows`-reordered manifest, were accepted. Fixed:
`prepare()` now takes `pseudobulk_result_path` and calls `pseudobulk_providers.
read_pseudobulk_result`/`load_pseudobulk_counts` ITSELF — there is no longer a caller-supplied
counts parameter to tamper with independently of the manifest; `PreparedDEDesign` also now
carries the manifest's own path/hash and the upstream `dataset_identity`
(species/gene namespace/counts semantics) forward. **F2** — declared `comparison`/
`sample_structure` column names were never checked against the pseudobulk aggregation's ACTUAL
`params`, and read-time validation recomputed rank/residual-df/estimability from `design_matrix`/
`contrast_vector` but never rebuilt the EXPECTED values from the persisted
`sample_ids`/`donor_ids`/`group_values`/`comparison` — and the "estimable" check is VACUOUS for
any full-rank design (`c @ pinv(X) @ X == c` identically for ANY `c`), so a reversed or zero
contrast vector, reversed group/donor roles, or a swapped reference/comparison level all
previously passed. Fixed: a column-name mismatch against upstream `params` is refused before any
design is built, and `_validate_prepared_de_design_core` now reruns `_prepare_de_design_core` on
the persisted identity and rejects ANY mismatch against the stored semantic fields. **F3** — the
environment plan's only install draft pointed at the shared `302006_env_clean`, and quoted
`lfc_shrink(coeff=None, adapt=True)` though `coeff` has no usable default (verified directly
against PyDESeq2 v0.5.2's own source). Fixed: the plan now proposes a NEW isolated environment
(`302006_de_backend_v1`, Python 3.11 + `pydeseq2==0.5.4` primary, `0.5.2`/Python 3.10 documented
fallback — the isolated env's own Python removes the reason to stay on `0.5.2`), backed by a real
executed `pip install --dry-run`/`conda create --dry-run` dependency resolution (the latter was
killed by OOM on the login node and is reported as inconclusive, not negative), and a corrected
backend-call draft (DataFrame counts/design + numeric contrast, explicit `lfc_shrink(coeff=...)`).
All seven of the review's own reproduction mutations (reversed/negated counts, reordered
sample_rows, wrong declared columns, reversed/zero/swapped contrast and group/donor roles) were
independently re-run against the fixed code and are now rejected (previously all accepted). New
standalone example `examples/run_de_design_example.py` (human unpaired + mouse fully-paired,
through the real `aggregate -> prepare -> write -> read/load` chain, including the
sample-id-vs-donor-id alignment case and the F2 column-mismatch refusal). **Full test suite: 331
tests pass** (up from 320: 11 new de-design tests). **Same standing caveat applies.**

**Monocyte severity pilot — science-decision package** (`review:20261014_204536`,
2026-10-14): the owner recorded a new decision (`docs/decisions.md`) selecting the question —
within one COVID-19 cohort, severe vs. mild/moderate, monocyte cell states/programs — `confirmed`.
Four documents were produced this round, metadata reads and literature checks only, no real
DE/model run: **`results/monocyte_pilot_feasibility_v1.tsv`** — real HPCC metadata (obs-only
h5py reads, no X-matrix load) showing the monocyte-cell-type-annotated universe in this pipeline
is currently exactly {Ren (`302004data15`, 53 patients: 16 mild/0 moderate/28 severe/9 fatal),
Stephenson (`302004data17`, 90 patients: 32 mild/28 moderate/22 severe/8 fatal), Qiu-sepsis
(excluded by disease)} — no other admitted COVID cohort has a monocyte extraction built yet, and
real monocyte-subtype labels (classical CD14 vs. CD16) exist only in each cohort's PER-cohort raw
file, not the currently-used merged object. **`MONOCYTE_PILOT_PREREG_v1.md`** — cohort (Stephenson
`proposed` primary, Ren `proposed` alternate), severity-value handling (`fatal`/`recovered`
folding), subtype scope, repeated-visit incidence, and the design-matrix/contrast direction
(severe minus mild/moderate, `confirmed`) each marked `confirmed`/`proposed`/`pending_owner`
explicitly; reuses `de_design_providers.py`'s own existing refusal (a donor fixed effect on a
one-sample-per-donor design saturates the model) as the stated reason an independent-donor
comparison must use the unpaired design type, never a donor-fixed-effect model. **`POSITIVE_
CONTROL_EVIDENCE_v1.md`** — three literature candidates independently re-verified via `WebSearch`
(2026-10-14, title/journal/DOI/core-finding level, not full-text): Schulte-Schrepping (*Cell*
2020) and Combes (*Nature* 2021) both confirm an ISG/interferon program HIGH in mild COVID,
largely absent in severe; Lee (*Sci Immunol* 2020) independently confirms the OPPOSITE — type I
IFN co-occurring with TNF/IL-1β inflammation specifically in SEVERE classical monocytes — a real,
confirmed literature conflict, not a summarization artifact; no positive-control direction frozen.
Schulte-Schrepping is cited here as literature ONLY — `docs/decisions.md`'s 2026-09-19 permanent
exclusion of that cohort as a DATA source is unchanged. **`MODEL_DECISION_v1.md`** — restates
scPRINT-2/`small-v2.ckpt`'s already-verified checkpoint/species status from `MODEL_SELECTION_v1.md`
with an explicit UCE trigger condition (activated only on a real scPRINT-2 runnability failure,
never on a comparative score), and proposes separating a masking score's TARGET (model-derived,
computed without reference to any pre-chosen gene set) from its VALIDATION (an after-the-fact
overlap check against the positive-control set), specifically to avoid the circularity the review
flagged (defining a masking score using the same set that is also the positive control). **Same
standing caveat: a feasibility/decision package, not an executed pilot — no real DE, state
discovery, or pathway run; cohort, severity mapping, subtype scope, positive-control direction,
checkpoint confirmation, and the validation split all remain the owner's to confirm.**

**R1-R4 follow-up fixes** (`review:20261014_221155`, 2026-10-14): an independent review
of commit `90c00d3` found the F1/F2 fixes above only closed `prepare()`'s entry point; the
public `write_prepared_de_design`/`load_prepared_de_counts` could still be called directly with
an internally-consistent but wrong artifact (swapped species/namespace/donor IDs, a stale
upstream hash, a nonexistent upstream path, or reversed/negative counts handed straight to the
writer) and both silently accepted it. **R1, fixed**: a new `_source_binding_problems` check,
shared by the write and read/load paths, re-resolves the claimed upstream pseudobulk manifest
LIVE and refuses any mismatch against its real `dataset_identity`/column-name `params`/per-sample
donor-group identity/`gene_ids`/counts content, plus a standalone non-negative/finite/integer
counts contract — all eight of the review's `probe_binding.py` cases, re-run VERBATIM against the
fixed code, are now rejected (10 new/updated tests; 341 tests pass, up from 331). **R2**: implementation assistant's
own prior round had already found the monocyte-canonical object's `covid_severity` was a
collapsed survived/fatal binary, but the "richer" field it substituted (joined from the
PLATELET-level canonical object) ALSO disagreed with Stephenson's OWN native field. A real,
obs-only (`h5py`, X never loaded) extraction from Stephenson's raw object
(`covid_portal_210320_with_raw.h5ad`) — after one OOM-killed login-node attempt was fixed by
masking before building any full-atlas DataFrame — found the cohort's own
`Status_on_day_collection_summary` field gives Asymptomatic=9/Mild=23/Moderate=30/Severe=13/
Critical=15 (not the earlier mild=32/moderate=28/severe=22), with a CONFIRMED-clean identity join
(128,503/128,503 cells matched, 0 conflicts) to the existing monocyte extraction, 12/90 patients
with a real D0/D28 repeated visit (severity consistent across both visits for all 12), and Site
as a real 3-way confound (Cambridge/Newcastle/Sanger). Ren's equivalent check was attempted and
found genuinely blocked (its local raw files carry no severity field at all). **R3**: gave
`MODEL_DECISION_v1.md` a concrete masking-score contract (a frozen linear severity-direction
score `s` fit on development donors, `Δ(g)=s_full-s_masked(g)` on held-out donors, PCA run
through the identical logic) and corrected "scPRINT-2 standing default" to "candidate" (matching
`docs/decisions.md`'s own wording); downloaded and froze the real MSigDB `HALLMARK_INTERFERON_
ALPHA_RESPONSE` gene set (97 genes, CC-BY-4.0, 96/97 covered in Stephenson's own gene panel) as
`gene_sets/`, and softened `POSITIVE_CONTROL_EVIDENCE_v1.md`'s "direct contradiction" framing to
an acknowledged cross-cohort disagreement, plus an explicit caution that a non-detected program is
not itself evidence of pipeline failure. **R4**: re-ran the dependency dry-run with
`--ignore-installed` (a real, full 40-package closure, not reusing `302006_env_clean`'s installed
packages), re-verified the backend call against `pydeseq2==0.5.4`'s own source directly, and
switched the deployment draft to an explicit `--prefix` path. **New feasibility table**: `results/
monocyte_pilot_feasibility_v2.tsv` supersedes v1's Stephenson row with the verified numbers and
flags Ren's as `UNVERIFIED`. **Same standing caveat: engineering fixes + metadata verification
only — no real DE, state discovery, or pathway run; cohort/severity-mapping/subtype/positive-
control/checkpoint decisions remain the owner's.**

**Incremental follow-up round** (`review:20261014_231259`, 2026-10-14): a further
independent review of `546bee0` found two local issues and one real further science finding.
**R1/P2 (fixed)**: the R1 source-binding counts check used `np.allclose`'s default relative
tolerance, which already accepts an off-by-one drift at realistic raw-count scale (100000 vs.
100001) and a non-integer value (100000.25) — fixed to EXACT equality; genuine before/after
confirmed against the unmodified `546bee0` baseline. **343 tests pass** (up from 341).
**Metadata-script fix**: `stephenson_raw_metadata_extraction.py`'s categorical decoder had a
latent bug (the newer-encoding branch didn't handle a missing code the same safe way the legacy
branch did) — fixed; confirmed NONE of the actually-read columns ever hit it, so prior real
findings are unaffected; its duplicate/unmatched/conflict/field-consistency checks are now hard
failures, not print-and-continue. **Real science finding**: Stephenson's 90 `patient_id` strings
were found to OVERCOUNT true biological donors — a direct crosswalk against the cohort's own
depositor SDRF (`Characteristics[individual]`) found 3 confirmed hidden same-donor groups (1
triple + 2 pairs) invisible to `patient_id`-string matching, 2 of which showed a REAL severity
PROGRESSION between visits. Corrected, true-donor-level (n=86) severity: Asymptomatic=9/Mild=23/
Moderate=30/Severe=13/Critical=11 (down from Critical=15). New `results/
monocyte_pilot_feasibility_v3.tsv` and `results/stephenson_donor_identity_roster_v1.tsv` (a
per-sample draft identity roster) supersede v2. `MODEL_DECISION_v1.md`'s masking-score contract
sharpened (PCA/scFM each fit their own `w`/`b` in their own space; explicit cell→donor
aggregation; one primary statistic chosen in advance; signed `Δ(g)` vs. `|Δ(g)|` distinguished;
pretrain-corpus overlap `unknown`, not asserted absent; gene-set coverage split into "final DE
universe" vs. "model-measurable" vs. the preliminary raw-panel check). New
`DE_BACKEND_EXECUTION_PACKET_v1.md` consolidates the dependency-resolution evidence/isolated-env
path/execution commands for the next phase's review. **Same standing caveat applies.**

**Round following `review:20261015_104226` and `_114223.md`, 2026-10-15.** Four
scoped workstreams: (1) the SDRF donor-identity crosswalk and roster are now generated by ONE
script so they cannot diverge again; **86 is confirmed the SDRF-matched true-donor count
outright** (both previously-"unresolved" Cambridge samples resolve unambiguously via their own
`patient_id`'s other SDRF rows) — the remaining open item is a *different*, already-tracked
question (which 10x pool/batch those 2 samples belong to, `stephenson_identity_decision_sheet_
v3.md`), not the donor count. (2) `scripts/34_scfm_discovery_framework/results/
stephenson_author_code_investigation_v1.md`: read the authors' own `COVIDPBMC` repo
(`myeloid_signature_scoring_publish.ipynb`, commit `9f51fe9a`) — their Fig.2f/Table 7 IFN
scoring uses a hardcoded 99-gene GO:0034340 list (only 24/99 overlap this project's frozen
97-gene Hallmark set — confirmed NOT the same asset); the 6-to-3 severity-collapse rule and the
real `severity`/`Days_from_onset` definitions are confirmed NOT in the public repo (a named gap,
not an unresolved search). `POSITIVE_CONTROL_EVIDENCE_v1.md` upgraded accordingly; IFN-down
remains a `CANDIDATE`, not frozen. (3) `ExpressionDERankingProvider.run_pydeseq2_from_prepared_
design(...)` implemented in `gene_ranking_providers.py` — consumes the already-validated
`PreparedDEDesign` as-is, keeps raw Wald/LFC separate from shrunken LFC, returns an honest
`not_run` on the real (confirmed) `ImportError` since `pydeseq2` is still not installed anywhere
accessible; the numerical pydeseq2 call path itself was NOT executed this round (code-reviewed
and skip-gated only). 348→(+6 own, 76 file total)/1 skipped overall. (4)
`DE_BACKEND_EXECUTION_PACKET_v1.md` made concrete (real absolute paths, hashes re-verified
unchanged via live `sha256sum`, draft SLURM table, one scoped authorization framing) —
still nothing installed.

**Separately, Cindy made two direct scientific decisions in a automated review conversation** (recorded in
`docs/decisions.md` 2026-10-15, authored by automated review, confirmed directly by the user in this
session after an attribution mix-up — see that file's correction entry): the main comparison is
now **(Severe + Critical, confirmed ICU location at sampling) vs (Mild + Moderate)**, and a
later sepsis-ICU vs. COVID-19-ICU comparison is a confirmed future direction. implementation assistant implemented
the mapping (`results/stephenson_icu_severity_grouping.py`, derived only from the real
`Status_on_day_collection` string, never from the severity grade) and the evidence report
(`results/stephenson_icu_severity_grouping_v1.tsv`/`_summary.json`, 9 passing tests incl. a
mutation-kill check): in this cohort, EVERY Severe/Critical sample already has confirmed ICU
status (no non-ICU/unknown-location exclusions triggered) — donor-level (any-visit, not yet
deduplicated by visit): comparison=25, reference=53, excluded-other=9 of 86 true donors, with
exactly 1 donor (`sdrf::C-8821`) whose visits land in BOTH arms — the concrete case the next open
decision (visit selection, §6 of `MONOCYTE_PILOT_PREREG_v1.md`) needs to resolve.
`MONOCYTE_PILOT_PREREG_v1.md` §3/§4/§6/§7 updated to lock this comparison in. **Same standing
caveat applies: no real DE run, no SLURM, no install, no push.**

**Follow-up round, `review:20261015_142903`, same day.** A mid-round attribution
mix-up (implementation assistant briefly, and wrongly, suspected a subagent had fabricated the above two
decisions) was corrected and documented in `docs/decisions.md` — the decisions are real,
written by automated review in its own concurrent conversation with Cindy. Then: **R1/R2 fixed** — a probe
found `run_pydeseq2_from_prepared_design` trusted a caller-supplied manifest/counts without
re-validation (same-shape content tamper, a fractional count silently truncated by the `int64`
cast, a flipped contrast sign, and an all-zero design matrix all reached the backend call
unrejected); fixed by re-running `de_design_providers.validate_prepared_de_design`/
`_source_binding_problems` (the SAME checks `read_prepared_de_design`/`load_prepared_de_counts`
already apply) at the consumption entry point via a lazy import, plus the numeric contrast
vector fixed from a plain `list` to `np.asarray(..., dtype=float)` (pydeseq2 v0.5.4's own
constructor, verified against its source, misreads a list as a 3-item categorical contrast and
raises `IndexError`). 6 new tests; before/after confirmed against the unmodified `bb9666d`
baseline. **R3**: the synthetic fixture's prior "null" gene was a hard-coded CONSTANT across
every sample, which distorted size-factor normalization (independently confirmed: it normalized
to log2FC≈-2.68, not ≈0) — replaced with a 7-gene design (known up/down effects + 4 stable,
non-constant, group-mean-matched background genes), independently sanity-checked by hand
(median-of-ratios: size factors 0.98-1.02, |log2FC|<0.02 for all null/background genes). **R4**:
`de_backend_requirements.lock.txt` (40 pinned packages) and `de_backend_phase_wrapper_v1.sh`
(`bash -n`-checked) make the execution packet's install sequence concrete; also found and fixed
that pydeseq2's own `DefaultInference` defaults to using every core on the node (verified
against its source) — `OMP_NUM_THREADS`/`OPENBLAS_NUM_THREADS` do not bound this separate
joblib pool — fixed via `_bounded_n_cpus()`, reading the real SLURM allocation. **Author
gene-set frozen**: `gene_sets/STEPHENSON_AUTHOR_MYELOID_SIGNATURES_v1.json` (IFN=99/TNF=200/
JAK=87 symbols, extracted by executing the pinned-commit notebook cell's own source, not
retyped) — freezing it is not adopting it. **The one flagged cross-group donor
(`sdrf::C-8821`) is now resolved**: Cindy chose directly ("进入目前分组少的那组" — assign it to
whichever group currently has fewer members), implemented as a one-time, hash-pinned decision
(`results/stephenson_icu_severity_grouping_resolve_mixed_visits.py`) — resolved donor-level
counts: comparison=25, reference=52, excluded-other=9 (of 86). **171 tests pass, 1 skipped**
(up from 164/1). Repeated-visit selection WITHIN an already-resolved group remains the one
still-open question. **Same standing caveat: no real DE run, no SLURM, no install, no push.**

**Follow-up review, `review:20261015_150231`, same day.** Independent re-test on HPCC
confirmed R1/R2 (171/1, unchanged) and accepted the C-8821 counts, then found 2 more real
defects (F1/F2), both fixed. **F2**: `stephenson_icu_severity_grouping_resolve_mixed_visits.py`
only checked that the input's own self-reported mixed-donor field equalled the expected donor —
a probe showed changing an UNRELATED donor's group, while leaving that stale field unchanged,
still produced an accepted (but wrong: 26/51 instead of 25/52) output and silently overwrote a
pre-existing output file (a sentinel placed at the output path was clobbered). Fixed: the script
now refuses unless both input files' sha256 match FROZEN, hardcoded values (never a hash
recomputed from "whatever was just read"), independently verifies the two named sample rows
belong to the expected donor and group, computes all conflicts before writing anything, and
refuses outright if its output path already exists. 6 new tests, including a genuine
before/after reproduction of the exact probe scenario against the unmodified prior script
(confirmed: it silently produced 26/51 and clobbered the sentinel; the fixed script rejects it
outright with 0 output written). **F1**: `de_backend_phase_wrapper_v1.sh` had 6 incomplete
points — an unused lock-file path placeholder, an unbound deployment location ("or wherever"),
no pytest dependency plan, treating `pip -c` as a package whitelist (it is not), an undecided
execution mode (`read -p` cannot work in a non-interactive `sbatch` job), and accepting a bare
pytest exit code as proof of a genuine pass (cannot distinguish passed from skipped). All fixed:
the wrapper now requires `REPO_CHECKOUT_DIR`/`EXPECTED_REPO_COMMIT`/a real SLURM allocation as
preflight checks before `conda create`; a new, separately-locked `de_backend_test_tooling.
lock.txt` (pytest==9.1.1 + its real transitive closure, verified via a REAL dry-run against the
main lock file as a constraint); a new `de_backend_verify_native_resolution.py` that exactly
compares the native resolver's package set against the lock file (replacing the human
confirmation prompt — functionally verified against real dry-run reports, both a match and a
deliberate mismatch); and a new `de_backend_verify_synthetic_test_passed.py` that parses pytest's
own JUnit XML and refuses unless the named test genuinely passed (functionally verified against
a real skipped run and a real passed run). **177 tests pass, 1 skipped** (up from 171/1). **Same
standing caveat: no real DE run, no SLURM, no install, no push.**

**Scope clarification, `review:20261015_224111`, same day.** The 25/52/9
independent-donor comparison view is confirmed to be ONE SPECIFIC analysis view, never a global
filter — repeated visits stay in the project's full `donor × visit × state` index for scFM
discovery and within-donor longitudinal work. New, metadata-only `results/
stephenson_longitudinal_feasibility_inventory.py` (3 tests, all passing, matching
independently-predicted counts): of 86 true donors, 15 have >1 visit (14 pairs + 1 triple);
**only 2 show ANY severity-tier transition and only 1 (`sdrf::C-8821`) crosses the
comparison/reference arm** — severity-tier change and arm change are different questions
(`sdrf::C-8896`'s Severe→Critical transition stays ICU→ICU, i.e. arm-stable). `Days_from_onset`
is now confirmed IDENTICAL across visits for ALL 15 multi-visit donors (not just some), with the
generation rule still unconfirmed. New `mono_LONGITUDINAL_FEASIBILITY_NOTE_v1.md` recommends (A)
full-index scFM discovery and (B) case-level description of the 2 transitioning donors as the
two well-supported near-term options — explicitly NOT a population-level progression model
(too few within-donor transitions in this cohort) and explicitly not a new decision on implementation assistant's
part. `de_design_providers.py`'s existing design-type guards are confirmed unchanged (no
weakening to force repeated visits through the unpaired backend); a future genuine
repeated-measures design type is named as a feasibility option only, not implemented.
`MONOCYTE_PILOT_PREREG_v1.md` §6 updated with this scope clarification. No new SLURM, install,
model run, real DE, push, or raw-expression-matrix read this round.

**Follow-up review, `review:20261015_230801`, 2026-10-16.** Found 3 real defects
in the longitudinal inventory and 2 real wording overstatements in its note, all fixed. **F1**
(overwrite): the inventory script's own test deleted the committed v1 results and `main()`
unconditionally overwrote its output path — fixed with explicit `roster_path`/`grouping_path`/
`out_tsv`/`out_summary` arguments, a refuse-to-overwrite guard, and a new `_v2` default output
path (v1 confirmed byte-identical before/after). **F2** (onset miscount): `sdrf::CV0231`'s
`days_from_onset` is the identical NON-NUMERIC string `"Healthy"` on both visits — v1 silently
folded this into "identical," reporting `unparseable=0` when it should have been 1. Fixed:
`_classify_onset` now separates raw-identity from numeric validity into 4 disjoint buckets
(`identical_numeric`=14, `identical_non_numeric`=1). **F4** (trajectory truncation): v1 only
ever compared a donor's first and last visit (discarding `sdrf::C-8896`'s middle visit) and
would silently sort an unparseable `collection_day` to the end. Fixed: `_order_donor_visits`
returns an explicit `order_resolved` flag, and `_any_adjacent_transition` checks every
consecutive step. Both F2 and F4 were mutation-killed directly against the real HPCC copy of the
script (reverting each fix broke the exact tests expected, then the clean file was re-synced and
re-hashed to confirm byte-identity before any further run). **New coverage evidence**: a bounded
`obs`-only read (never `X`/`layers`) of the frozen per-cohort raw object `stephenson_monocyte_
raw.h5ad` (`results/stephenson_sample_subtype_cell_counts.py`, sha256 of the source object
`11eb7450ff...4b`) found all 102 roster `sample_id`s match this object's 102 `Sample_ID`
categories exactly; joined via `results/stephenson_longitudinal_monocyte_coverage.py`, **all 12
symptomatic multi-visit donors share at least one monocyte subtype (`CD83_CD14_mono`) present
across every visit**, and the 2 clinically-transitioning donors share all 5 subtypes across
every visit — a real, named mapping gap is reported as `unknown`, never folded into a 0 count.
Also corrected the note's 12-pair severity breakdown (had omitted 2 `Critical→Critical` pairs)
and its §3 claim that a new THIRD DE-design type was required even for the simplest 2-visit
case (the EXISTING `paired_or_repeated_within_donor` design already covers it when `condition` is
defined as Initial/Follow-up rather than severity; a third type remains needed only for
`sdrf::C-8896`'s 3-visit case and genuine partial-overlap/covariate generalization). 15 new
tests across 3 files; **392 tests pass, 1 skipped**, independently re-confirmed this round (not
merely trusted from the prior round's reported count) by re-running the full relevant scfm/
stephenson suite in a fresh isolated HPCC snapshot. `MONOCYTE_PILOT_PREREG_v1.md` §6 and
`mono_LONGITUDINAL_FEASIBILITY_NOTE_v1.md` updated throughout (§1/§1a/§3/§4/§5). **Same standing
caveat: no new SLURM, install, model run, real DE, push, or raw expression-matrix (`X`) read.**

**Second follow-up review, `review:20261016_102834`, 2026-10-16.** Accepted the
raw-data coverage inventory but rejected the "subtype coverage is no longer a limiting factor"
framing as overstated (presence ≠ usable sample size). Found and fixed 2 more real defects, which
changed no already-committed number: **R1** — `stephenson_sample_subtype_cell_counts.py` indexed
categorical codes directly; anndata's `-1` ("missing") code would silently resolve to the LAST
category via Python negative indexing instead of raising (reproduced on a fixture: `sample_codes=
[0,-1]`/`subtype_codes=[-1,0]` silently counted; fixed with explicit code-range validation that
refuses any code outside `[0, len(categories))`; the real object has 0 missing/out-of-range
codes). **R2** — `stephenson_longitudinal_feasibility_inventory.py` gated
`severity_any_adjacent_transition` on `order_resolved` but not the first-to-last fields or the
summary's transition-count dicts, so an order-unresolved donor could still get a direction label
from arbitrary input row order (reproduced on a fixture: `order_resolved=False` yet
`severity_transition_first_to_last="worsening"`, counted under `{"Moderate->Critical": 1}`; fixed
— unresolved donors now read `"unknown"` and count under a separate `"order_unresolved"` summary
key; the real roster has 0 order-unresolved donors). Both mutation-killed directly against the
real HPCC script copy (injected, confirmed the exact expected failures, restored, re-hashed
byte-identical before any further run). New `results/stephenson_monocyte_pilot_design_readiness_
table.py` reuses the already-verified roster/grouping/counts/coverage tables (no new h5ad read)
to produce a donor×visit×subtype real-count table plus a per-subtype 1/20/50/100-cell sensitivity
listing across all 5 subtypes (not just `CD83_CD14_mono`) — explicitly descriptive, not a new
inclusion rule. Confirmed **no QC-admitted monocyte analysis input exists** anywhere in the
repository; only the raw, provisional `stephenson_monocyte_raw.h5ad` is available
(`qc_admitted_input_status` in the new summary JSON). Corrected
`mono_LONGITUDINAL_FEASIBILITY_NOTE_v1.md` §4(C)'s overstated "14 fully-paired donors (or any
subset)" (3 of the 14 are Asymptomatic, outside the current symptomatic comparison scope — the
real in-scope population is 11 symptomatic 2-visit donors + `sdrf::C-8896`'s 3 visits) and its
implicit pre-selection of `CD83_CD14_mono`. 5 new tests; **397 passed, 1 skipped** on the full
relevant scfm/stephenson suite, independently re-run this round (up from 392/1). **Same standing
caveat, still true: no new SLURM, install, model run, real DE, push, or raw expression-matrix
(`X`) read.**

**Third follow-up review, `review:20261016_120026`, 2026-10-16 — corrects two
items from the paragraph above.** (1) Repeated-visit selection is NOT an unanswered blocking
question: Cindy already asked to check sampling-time distribution/overlap first, then asked
directly why not use both samples — the correct status is **DEFERRED pending time evidence and
analysis design**, full visit index available and unblocked, not re-listed as an open blocker
alongside (A)/(B)/(C). (2) "No QC-admitted monocyte input exists" was itself overstated — that
search only covered the repository, not HPCC's actual QC data directories. Real, bounded,
read-only checks (cell-id joins only, no 12G/25G matrix read) found Gate 1 QC applied for
Stephenson (2026-09-21) and DecontX complete for Stephenson (2026-09-18, mean contamination
13.0%), both joining 100% by `cell_id` to all 128,503 monocyte-pilot cells (0 mismatches;
Gate 1 `flag_would_delete` rate and DecontX contamination both vary meaningfully by subtype —
see `results/stephenson_monocyte_pilot_design_readiness_table.py`'s `QC_SOURCE_EVIDENCE`).
Corrected status: **`qc_admission_not_verified_for_this_pilot`** — real QC/decontamination
artifacts exist and are cell-identity-joinable, but no decision exists that either one
constitutes admission for THIS monocyte pilot, and the decontaminated counts' actual values
(not just the per-cell score) were not read. 6 new/updated tests; v2 design-readiness outputs
(v1 left untouched, byte-identical before/after). **Owner decision, 2026-10-16:** Cindy replied
"允许" to the scoped proposal in `review:20261016_130032`. Existing atlas raw
counts are accepted for an input-to-state-discovery technical pilot only, retaining all visits
and attaching historical QC flags/scores without filtering or replacing counts. This resolves
provisional technical input use, not final scientific QC admission or disease/pathway claims.
Next: prepare the concrete execution package; model/checkpoint confirmation and job-specific
SLURM authorization remain outstanding. No job is authorized by this decision. Repeated-visit
selection remains DEFERRED. See `docs/decisions.md` and `review:20261016_130248`.

**Execution-packet round, `review:20261016_130032`, 2026-10-16.** Fixed the
review's 3 wording issues (historical-vs-current QC completion scope; DecontX ESTIMATED score
no longer characterized as proven contamination or used to imply a subtype preference;
evidence citation now a committed, re-runnable script + its real JSON output instead of
handoff prose). Delivered `scripts/34_scfm_discovery_framework/
STEPHENSON_MONOCYTE_SCPRINT2_PILOT_EXECUTION_PACKET_v1.md`: reuses the one proven-working
scPRINT-2 entry point (`run_scprint2_human_covid.sh`/`scprint2_env`, confirmed live and
responsive again this round) and `model_registry.py`'s small-v2.ckpt candidate (still a
candidate). Named and implemented the one real gap — the raw monocyte object mixes 192
Antibody-Capture features into its var axis, unlike the one cohort scPRINT-2 has actually
run on in this project — as `results/stephenson_monocyte_scprint2_pilot_prep.py` (GEX-only
subset + Gate 1/DecontX annotation ONLY, never filtering cells or substituting counts):
implemented and tested against a synthetic fixture (8 tests, 2 of them mutation-kills), NOT
run against the real 5.7 G object (no SLURM authorized; forbidden on the login node).
Wrote, did not submit, `run_scprint2_stephenson_monocyte_pilot.sh`. 22 new/updated tests this
round; full relevant suite independently re-run. **Remaining for the owner**: final
model/checkpoint confirmation, per-job SLURM authorization, and whether the reused 2h/128G
resource ceiling (evidenced on a smaller input: 47,977×13,066 vs this input's 128,503×24,737)
is sufficient — stated as UNVERIFIED, not assumed.

**Owner decision + packet v2, `review:20261016_141242`, 2026-10-16.** Cindy
replied "ok": **scPRINT-2 `small-v2.ckpt` is confirmed as the FIRST technical-pilot model**
(not the final manuscript backbone; UCE not excluded). Independent review found real defects
in the v1 packet, all fixed by reading the ACTUALLY INSTALLED library source directly (not
inferred from an old artifact's existence): **R1** `cli.py` writes `<--output_filename>
_embedded.h5ad`, not the literal path — v1's script checked the wrong path. **R2**
`scdataloader/preprocess.py` merges on Ensembl IDs when `is_symbol=False`; our var_names are
HGNC symbols — v1 wrongly left this unset. A real probe against the installed gene table (no
model weights touched) confirmed `--is_symbol true` maps 19,194/24,737 genes, clearing the
library's own 10,000-gene floor. **R3** the library unconditionally filters cells below 200
detected genes regardless of `skip_validate` — a real, bounded, chunked probe of this object's
own sparse structure (never densified, login-node-safe) found 0/128,503 cells fall below that
threshold. **R4** two of the packet's own "mutation-kill" tests were self-defeating (`try:
call(); assert False; except AssertionError: pass` catches its own assertion) and
`validate_prepared_adata()` trusted a caller-supplied report instead of re-deriving identity
from the real objects — both fixed and mutation-killed against the real production file (6
independent bypass scenarios from the review's own probe table, all now caught). **R5** the
job script lacked its own prep step, used `--requeue` (conflicts with this round's no-retry
scope), had no runtime checkpoint hash check, and ended with only `ls` — rewritten as one
bundled, fail-stop job (`--no-requeue`, runtime SHA256 check, a real validator
(`stephenson_monocyte_scprint2_pilot_validate.py`, new) as the job's own last step). **R6** the
new crosscheck script reintroduced the `-1`-categorical-code defect fixed elsewhere — fixed,
re-verified against real data (still 128,503/128,503 joined, 0 mismatches). Delivered
`STEPHENSON_MONOCYTE_SCPRINT2_PILOT_EXECUTION_PACKET_v2.md` (v1 left on disk, superseded, not
deleted). A known prerequisite gap, not fixed this round: the HPCC checkout the job script
points at is still stale (`35131cf`) with unrelated uncommitted changes. 72 tests
(new/updated) re-run this round, all passing. **Still not authorized: SLURM submission,
install, push, or final backbone selection.**

**Follow-up review + packet v3, `review:20261016_152556`, 2026-10-16.** Two real
gaps remained in v2, both closed: (1) `stephenson_monocyte_scprint2_pilot_validate.py` took
only `expected_n_obs` — a same-count output of unrelated ids, duplicate ids, or misrouted QC
annotations would have passed; rewritten to check exact cell-id coverage against the frozen
prepared-input manifest, the real primary `scprint_emb` key specifically (not any
`X_scprint*`-prefixed key — the installed `cell_emb.py` writes a UMAP to a separate
`X_scprint_umap` key), array dtype/shape, and QC-column correspondence by id; all 9 of the
review's probe scenarios now covered by tests, 2 independent mutation-kills performed and
reverted. (2) The job script's `PREP_RESULTS_DIR` pointed at the still-stale live checkout;
now points at `302006_updategit_deploy_scprint2_pilot`, an independent clean git export
(verified clean, verified at the authoring commit, hash-matched file-by-file against local) —
the stale checkout itself was never read, modified, or cleaned. Also fixed: the prep
validator's bare `assert` statements (stripped under `python -O`/`PYTHONOPTIMIZE`, verified by
reverting and re-running the new subprocess test) and an overclaiming mutation-kill test name.
Job script also now `cd`s into the per-run directory before `scprint2 embed` (the installed
library writes intermediate files to a relative `data/` path), verifies the three real input
files' sha256 before the run (not just after), and pins
`seed_everything`/`batch_size`/`num_workers`/`how`/`max_len` explicitly. 49 tests across the
touched files pass on HPCC. Delivered
`STEPHENSON_MONOCYTE_SCPRINT2_PILOT_EXECUTION_PACKET_v3.md` (v1/v2 left on disk, superseded).
**The only remaining item before this specific job is a per-job SLURM authorization request —
no further code gap is open.** Still not authorized: SLURM submission, install, push, or final
backbone selection.

**Follow-up review + packet v4, `review:20261016_202903`, 2026-10-16.** Four real
problems found in v3, all independently reproduced and closed, no new owner decision needed.
**R1**: `--time=3:00:00` exceeded `short_gpu`'s real `MaxTime=02:00:00` (confirmed directly via
`scontrol show partition short_gpu`) — lowered to `--time=2:00:00`; sufficiency for this larger
input remains genuinely unverified. **R2**: the validator's "recovered from the expected
manifest" claim for a missing required column was never backed by anything actually being
written — the review's own probe (a real 3-cell h5ad, `main()` called directly, output reopened)
showed the columns still absent while the validator returned a technical PASS; the fake-recovery
branch is now removed, missing is always a hard failure. **R3**: the validator checked only 2 QC
columns, never donor/sample identity — the review's `swapped_patient_id`/`swapped_Sample_ID`
probes (ids/embedding/QC all correct, only the identity value reversed) both returned
`hard_failures=[]`; `patient_id`/`Sample_ID` (confirmed real columns on
`stephenson_monocyte_raw.h5ad` via a bounded h5py read — there is no `visit` column on this
object) are now required and checked by id like the QC columns. **R4**: the job script recorded
whatever `HEAD` a clean deployment happened to have, so a fast-forward during queue wait would
still pass — the submitter must now supply `EXPECTED_DEPLOY_COMMIT` via `sbatch --export`,
compared against the deployment's live `HEAD`. Re-exported `302006_updategit_deploy_scprint2_pilot`
to the new commit `3de4bca` (same directory, same stale-checkout-untouched guarantee, verified
clean and hash-identical). 41 tests pass (the prior round's "49" claim for the 2 touched test
files was wrong; the real count for those 2 files alone is 38 before this round's 3 net-new
tests); 2 mutation-kills performed and reverted (hash-matched), plus the real review probes
re-run directly against the fix, all now correctly failing/passing as expected. The revision
guard was tested against both a disposable stand-in git directory and the real redeployed
directory. Delivered `STEPHENSON_MONOCYTE_SCPRINT2_PILOT_EXECUTION_PACKET_v4.md` (v1/v2/v3 left
on disk, superseded). **Update (2026-10-17): the single authorized run happened and failed on
identity, not on resources.** The user authorized exactly this one job; submitted as SLURM
`30044888` (pinned to `3de4bca`), it reached **FAILED** (ExitCode 1:0, Elapsed 00:52:05 — not a
timeout; MaxRSS ~35.4G of 128G — not OOM). `prep` succeeded; `embed`'s own output passed through
no crash but `validate` correctly hard-failed it: the installed `scdataloader==2.1.2`
`Preprocessor` unconditionally overwrites `obs_names` with fresh UUIDs before the model ever runs
(`preprocess.py:424`, no disabling flag in any of its 32 parameters), so the embedded output's
128,503 cell ids matched 0 of the prepared input's — VALIDATION_FAILED, classified and traced in
`implementation-note:20261017_164353`. All run artifacts preserved untouched. An identity-
recovery adapter (`stephenson_monocyte_scprint2_pilot_recover_cell_ids.py`, new, committed
`6a7601b`, 3 real defects found by independent review and fixed in `919cbcf`) relabels
`obs_names` from the preserved
`cell_id_stripped` column, which the UUID defect never touches; applied for real to job
30044888's artifacts (`recovered_v2/`, new versioned directory — `recovered_v1` kept as the
pre-fix record): 128,503/128,503 cells recovered, 0 QC/identity mismatches by id. This closes the
obs-identity question but NOT the separate question of whether `obsm["scprint_emb"]`'s rows stay
aligned with those now-correctly-labeled rows — that is addressed only by read-only source tracing
so far (no reordering operation found anywhere in the call chain), not by a live fixture test.
**Formal re-validation and resources/budget are not yet authorized; no new SLURM job, retry, or
resubmission has occurred.** Still not authorized: any new SLURM submission (including a
short CPU-only re-validation job), install, push, or final backbone selection.

The existing QC rebuild is unfinished, so historical downstream numbers remain provisional.
Its completion is a dependency for applicable real-data scientific claims, not a prerequisite for
framework interface development or small-fixture testing. No new HC2/Gate 2 job follows from this
scope change. The following table records inherited workstream state, not an instruction to resume
the former benchmark or manuscript plan.

| Phase | State |
|---|---|
| 8-cohort DecontX, initial pass | Complete |
| DecontX **refit** on the QC-retained set | Not started (waits on Gate 2) — this is a different run from the initial pass above |
| Gate 1 (basic QC, cohort-calibrated) | Rules frozen; applied; see caveat below |
| Gate 2 safety layer (P0.2–P0.5) | **Complete** — fail-closed config, status enum, fingerprinted checkpoints, runtime diagnostics; **190 tests** |
| Cohort runnability classification (Gate A) | **Conditionally approved** as a working classification — **not** a deletion permit |
| Sensitivity design (Gate B) | **Conditionally approved for scoring-only** — direct residual-rate stress grid; not prevalence or deletion authorization |
| Sensitivity-design review packet | Written AI methodology adjudication completed at the owner's request on 2026-09-24; current evidence PDF is 12 pages; no external human handover recorded |
| Gate 2 (Scrublet doublet scoring) | **Single-capture candidate audit passed** for GSM5102900_HC1; 11/11 parameters and deposit C1–C3 accepted. **Update (2026-10-11/12): permit + dry-run accepted, then the single authorized scoring pilot ran and completed.** Permit commit `b3b479f`; clean detached HPCC worktree `worktrees/b3b479fe...a9617/`; dry-run passed, automated review-reviewed (`review:20261011_234753`); user then authorized exactly one scoring job (short/2h/48G/4CPU), submitted as SLURM `29882644`, **COMPLETED 0:0, 00:03:16, MaxRSS 2.3G**. `SCORING_MANIFEST.json`: `status="complete"`, `n_cells=7666` (the full deposit, not a subsample), all 8 frozen priors examined (no `threshold_failed`), `outcome_columns_read=[]`, `deletion_authorised=false`, `production_parameter_ruling_written=false` (global slot still empty). Graph acceptance independently re-checked from the raw per-prior hashes (not just the manifest's own claim): all 8 constructions byte-identical (`k_adj=132`, shape `[22998,131]`, `acceptance_problems=[]`). Per-prior `n_called` (8/8/8/225/225/251/291/251) independently recomputed from the raw boolean vectors, exact match. This answers only whether the frozen scoring pipeline runs correctly on one capture — it is not a Harmony-vs-alternative decision, not a deletion authorization, and not a statement that Gate 2 is complete for the other seven cohorts. No second job submitted, no push, no outcome access, no cell deletion. Full evidence: the local implementation notes (this round's handoff).** |
| Stage 31 doublet re-audit | Full-capture-before-subsampling redesign **conditionally approved** on 2026-09-24; input/mapping validation and new dimension audit pending; old 30–60-cell/sample design not approved for scientific use |
| 5-FM full-scale embedding extraction | Not started |
| Batch-effect check | **In progress — 4 review rounds, 2026-10-08/09, 57 tests passing on HPCC.** Metadata (`scripts/33_batch_effect_audit/`, v4): all 12 current + 3 proposed Cambridge multi-donor dependency components enumerated with evidence status; collision-aware alias crosswalk; Stephenson capture/count units corrected to sample_id-level. **Stephenson/Cambridge pooling resolved for the Cambridge subset from the primary source** (Stephenson et al. 2021 *Nat Med*: pools of 4 combined pre-staining/loading, genotype-split) and cross-checked against this project's own `Site` field: 47 donors, currently 47 singletons in `canonical_donor_labels.csv`, collapse into 3 real components (21/22/4; two mix fatal/survived) — **proposed, not applied**. 83 non-Cambridge Stephenson donors remain unknown, not asserted independent. A known, documented (not silently resolved) limitation: the sample-to-pool pairing check is membership-only and cannot detect a swapped pool order for a donor in two pools. **Stage-33 diagnostic CLI implemented and validated end-to-end, twice** (`run_batch_diagnostics.py` + `make_synthetic_batch_fixture.py` + `prepare_provisional_diagnostic_input.py`, v5 2026-10-09, 4th review round): expression/PCA + optional FM-embedding mode, patient-balanced sampling, the shared scIB batch-ASW (`scripts/00_shared/code/batch_asw.py`; the second legacy caller now raises if invoked), a neighborhood-mixing diagnostic, and a real opt-in Harmony correction (transductive/exploratory, never out-of-sample). A 4th review found and this session fixed real bugs: both samplers used Python's per-process-randomized `hash()` (a fixed `--seed` did not actually reproduce a selection across processes); the real-data preparer capped cells per donor **across both lineages combined**, starving the rarer lineage (confirmed: re-running the real dry-run after the fix roughly doubled the usable cell count, 13,232→24,841, with Platelet recovering most dramatically, e.g. COMBAT 43→1,775); a requested-but-failed Harmony correction could leave the whole run `status="complete"`; `code_revision` could read `"UNKNOWN"` while still reporting success; evaluation-boundary leakage/group-integrity tests checked a parallel re-derivation instead of the actual harness's real fold/row identity. All fixed, re-validated: the latest synthetic run (`synthetic_fixture_v3`) completed with strict `status: "complete"` (not `"complete_with_partial_results"`) and a verified `code_revision`, showing the same qualitative result as before — the injected batch effect detected (neighbor-mixing ~0.09–0.18 pre-correction vs. ~0.75 expected under none) and reduced by Harmony (~0.48–0.69) in all 8 lineage×embedding×scope combinations; biological-conservation (KMeans+ARI) is near zero pre-correction when the batch effect dominates the naive clustering and rises sharply post-correction — a property of that specific metric, not a general "correction preserves biology" claim. **Expression/5-FM diagnostics on real cohort data are still not run** — blocked on resolved Gate-2 input + full-scale 5-FM embeddings; a bounded, two-allocation provisional job against the existing pre-Gate-2 input is specified (`scripts/33_batch_effect_audit/results/PROVISIONAL_JOB_SPEC_v3.md`, with real reviewable sbatch scripts) but **not submitted**. **5th review round (2026-10-09/10) found and fixed 5 execution-integrity defects independent of the scientific audit**: a stage-label free-text mismatch between the job spec's own two commands; `code_revision_verified` that was literally `code_revision is not None` (any string "verified"); a check-then-act output guard (the real root cause of this session's earlier login-node process race, not just the orphaned-process explanation given at the time); unknown-outcome cells withheld from the entire technical-integration population instead of just bio-labeled scoring; a failed corrected-UMAP not propagating to run status. All fixed with regression tests (`tests/test_execution_readiness.py`, 15/17 passing + 2 honest environment-dependent skips; `tests/test_preparer_to_diagnostic_integration.py`, a real end-to-end preparer→runner run against the exact job-spec v3 argument contract, passing). See `docs/execution-plan.md`, `scripts/33_batch_effect_audit/results/MANIFEST.md`, and `scripts/33_batch_effect_audit/results/FINDINGS_REPORT_v1.md` §10. **Addendum (2026-10-10): Job 1 executed.** After a 6th/7th review round found and fixed a stale/dirty deployment binding, a Job 1 acceptance-criteria mismatch, a premature `status="extracted"` publish, and a missing wrapper working-tree-clean check (`review:20261010_115720`, `20261010_120725.md`), the user authorized Job 1 only; it was submitted once (SLURM job `29804463`, `aea9a61` deployment) and **completed successfully**: 24,841 cells × 11,709 genes, matching the historical v2 estimate exactly cohort-by-cohort. A follow-up review (`review:20261010_191146`) found the completion handoff's claim that `anndata.read_h5ad(..., backed="r")` avoids loading sparse layers was **wrong** (it reads the installed anndata source and confirms `backed="r"` still materializes `layers` in full) — corrected in the local implementation notes (see the round's handoff) with a new bounded h5py-only validator, `scripts/33_batch_effect_audit/code/validate_job1_output.py`, whose metadata checks (all cells) and sampled, outcome-blind numeric checks (fixed small sample) both pass against the real output; full source-Matrix-Market numeric equality remains explicitly `NOT_VERIFIED`. **Addendum (2026-10-11): Job 2 authorized and executed.** The user authorized one run of `run_batch_diagnostics.py` (provisional cohort-level diagnostics, exploratory Harmony included) against the exact packet in `job2_preparation_packet_v2.md`. Submitted once via the staged wrapper (`stage33_wrappers/2ff4859/step2_run_diagnostics_v2.sh`, `aea9a61` application) as SLURM job `29868753`: **COMPLETED, ExitCode 0:0, Elapsed 00:03:47, MaxRSS ~3.25 GB** (well under the 16 GB allocation). `manifest.json`: `status="complete"`, `code_revision_verified=true` (method `git_head_clean_tree`), `input_sha256` matches the packet's expected value exactly, both lineages' `expression_pca` embedding completed with BOTH `uncorrected` and `harmony_corrected` scopes present (Harmony succeeded, not merely requested). Neighbor-batch-mixing rose from 0.061→0.440 (Monocyte_Myeloid) and 0.104→0.407 (Platelet) post-Harmony — consistent with the synthetic-fixture pattern already validated pre-submission — while `bio_ari` stayed near zero in both scopes for both lineages (same "near-chance clustering agreement" caveat already on record for this metric, not a biology-preserved claim). 11,959/12,998 Monocyte_Myeloid and 10,824/11,843 Platelet diagnostic cells have a known outcome and are bio-scoring-eligible; the rest (including all 40 of the disputed donor `CV0201`'s cells, confirmed outcome-`unknown`, split 20 Monocyte_Myeloid / 20 Platelet across the two lineages, NOT 40 in each — see the Stephenson row below) remain in the technical/neighbor-mixing diagnostic only. This is a provisional, pre-Gate-2 exploratory run — not a final correction decision, not a method comparison, not a manuscript result. Full detail: the local implementation notes (this round's handoff). **Correction (2026-10-11, same day): a review (`review:20261011_183017`) found the mixing-only framing above incomplete.** Outcome-stratified `batch_asw_neg_scib` improved post-Harmony in Monocyte_Myeloid (0.8444→0.8972) but got slightly WORSE in Platelet (0.8591→0.8492; both outcome strata individually worsened: fatal 0.1398→0.1583, survived 0.1420→0.1433 absolute silhouette, higher = worse by this metric's convention — `batch_asw_abs_silhouette` averages `abs(signed silhouette)` per cell, so this alone does not establish that cohorts became more geometrically separated rather than more strongly cross-cohort-inverted; no signed breakdown was computed). Neighbor-mixing's improvement also differs by lineage, not "roughly quadrupled in both": 7.17× (Monocyte_Myeloid) vs. 3.92× (Platelet). `bio_ari` near-zero in both scopes answers only whether a 2-cluster KMeans partition matches the fatal/survived label — not a biology-preservation claim in either direction. No re-tuning, no method decision, and no QC change follow from this. Full corrected per-lineage table: `scripts/33_batch_effect_audit/results/JOB2_PROVISIONAL_INTERPRETATION_v1.md`. |
| Stephenson source-identity crosswalk (4th pass correction) | **Open (2026-10-11, 4th pass — corrects the 3rd pass row below's "zero-consequence" overclaim and closes the audit tools' remaining concurrency gap; does not replace the confirmed facts above it).** A same-day follow-up review (`review:20261011_124253`) found: (1) the 3rd pass's claim that excluding `CV0201`/`CV0326` is "zero-consequence for every other donor" overstated the real result — removing 2 members from the 21-donor component does not fragment the remaining 19 (they stay one group), but this says nothing about how those 2 donors' OWN 73 cells should be treated in any future CV-dependency grouping; corrected to a conditional statement (see `stephenson_identity_decision_sheet_v3.md` §2–3) with an interim policy recommendation (keep both source mappings + conflict flag, do not auto-exclude, do not change canonical); (2) both new audit tools' `--out-dir` guard still had a real (if narrow) concurrency gap for two FRESH simultaneous starts — fixed by adopting the project's existing `claim_ownership()` primitive (same one `run_batch_diagnostics.py` already uses), with a direct atomicity test and a CLI-level companion test; (3) the BCR evidence numbers (1,816/645/1,230) were CONTIG ROW counts, not deduplicated cell counts — `stephenson_disputed_donor_impact.py` now reports both separately; (4) a prior claim of "never reading outcome values" was corrected to describe what actually happens (the column IS loaded to compute an aggregate known/unknown count using the SAME `bio_known_mask` definition `run_batch_diagnostics.py` uses — not merely non-null, since the string `"unknown"` is itself non-null — individual values are still never printed/returned/used to filter cells); (5) added a real full-CellID join between Job 1's 73 selected disputed-donor cells and GEX's own cell_donor map (matched/unmatched/mismatch counts, per pool×lineage breakdown) where the prior round only reported pool-prefix sets. The v2 wrapper fix is now actually deployed (not just committed) to `stage33_wrappers/2ff4859/step2_run_diagnostics_v2.sh`, sha256-verified against the git blob, `bash -n`-checked, and its `population_stage` h5py read re-verified against the real file — the frozen `aea9a61` application remains untouched. See `scripts/33_batch_effect_audit/results/stephenson_identity_decision_sheet_v3.md` and `job2_preparation_packet_v2.md` for full detail. |
| Stephenson source-identity crosswalk (3rd pass correction) | **Open (2026-10-11, 3rd pass — corrects the 2nd pass row below's wording, does not replace its confirmed facts).** A 2026-10-11 review (`review:20261011_113512`, F1) found the 2nd pass row overclaimed: GEX is NOT a fifth independent source (it is one more OBSERVED, processed per-cell label, possibly sharing BCR/TCR's own demultiplexing lineage — though its `patient_id` field IS now independently checked against its own `sample_id`-derived donor: 0 mismatches across all 130,644 real Cambridge cells); the "96,242 disagreements" figure wrongly summed two overlapping cell-set intersections (21,481 GEX∩BCR + 74,761 GEX∩TCR) without checking their overlap — the real distinct union is **94,457** cells (1,785 shared by both BCR and TCR), still 0 disagreements; "decision option 3 has zero consequence" was asserted, not tested — now tested (`component_diff_excluding_disputed_donors.csv`): excluding `CV0201`/`CV0326` entirely shrinks the 21-donor component to 19 (sizes `[22,19,4]`, no fragmentation), so it is ALSO zero-consequence for every other donor, now confirmed rather than assumed. **New evidence this pass**: GEX's own `Resample`/`Collection_Day` fields (not read in the 2nd pass) show `CV0201` has a genuine D0/D28 two-visit pattern — `Initial`/`BGCV08` then `Resample`/`BGCV06` — matching the same pattern already established for 11 other donors, and BCR directly OBSERVES `CV0201` in `BGCV06` (1,816 cells) and `BGCV08` (645 cells), never `BGCV13` — real per-cell evidence, not merely a subset/coverage-gap pattern. These 2 donors' cells ARE present in Job 1's output (40 + 33 cells; checked with the new `stephenson_disputed_donor_impact.py`), but Job 2's call path does not read pool/`cv_group` identity, so the conflict still does not block Job 2 (see `scripts/33_batch_effect_audit/results/job2_preparation_packet_v1.md`, new this pass, which also documents and fixes — in a new, undeployed wrapper file only — a `backed="r"` metadata-only-read defect inherited from the original R1 finding). Full detail: `scripts/33_batch_effect_audit/results/stephenson_identity_decision_sheet_v2.md`. |
| Stephenson source-identity crosswalk | **Open (2026-10-10, 2nd pass — corrects the first pass's wording).** `audit_stephenson_source_identity.py` now cross-references FIVE sources for all 15 Cambridge `BGCV*` pools: the official E-MTAB-10026 SDRF, this project's local `stephenson_pool_composition.csv`/`canonical_donor_labels.csv`, the REAL deposited/processed GEX h5ad's own per-cell `obs['sample_id']` (ground truth of what was actually sequenced, read via `h5py`, legacy AnnData 0.1.0 encoding), BCR, and TCR. **Exactly 2 pools are a real, unresolved disagreement**: `BGCV06`/`BGCV13` (`CV0201`/`CV0326`) — SDRF alone says `CV0201∈BGCV13`+`CV0326∈BGCV06`; the local table, BCR, AND the real GEX data all agree on the opposite (`CV0201∈BGCV06`+`BGCV08`, `CV0326∈BGCV13`). A full-CellID per-cell join (not pool-level sets) between GEX and BCR/TCR for the whole Cambridge subset found **zero** disagreements (21,481 GEX∩BCR cells, 74,761 GEX∩TCR cells, 100% agreement) — strong, not dispositive, evidence toward the local/BCR/GEX side; SDRF is the official deposit and this is not resolved by vote. Previously reported as a "3-way conflict" and "12 resample donors" — both corrected: TCR silence on `CV0326` is a coverage gap, not a third vote, and of the 12 Cambridge donors seen in >1 pool, 11 are real `Initial`/`Resample` (D0/D28) pairs and 1 (`CV0902`) is same-day/same-flag in two pools, not a resample. **The existing proposed 21/22/4-donor dependency components are UNAFFECTED by this conflict either way** (`component_diff_sdrf_vs_local.csv`: 0 donors change component — `BGCV06`/`BGCV08`/`BGCV13` are already bridged by other, unaffected donors). Job 2's wrapper/runner (`run_batch_diagnostics.py`, `step2_run_diagnostics.sh`) uses `--batch-col cohort --patient-col patient_id --lineage-col lineage --outcome-col outcome` — grepped directly, no `cv_group` reference anywhere in that call path, so this conflict does **not** block Job 2 on that basis (a prior handoff's claim to the contrary is retracted). See the local implementation notes for exact result paths. |
| LOCO / statistics / results integration | Not started |
| Manuscript recompile | Required after the above; the 2026-08-06 compile is stale |

**No Gate 2 job is running.** All eight cohorts' Gate 2 outputs are quarantined (see below);
the Liu run reached its walltime and was deliberately not rescued, since P0 requires a rerun
regardless. The archived `MEMORY.md` claim of "no pending HPCC runs" was nevertheless wrong
for the period it covered — the rebuild was active throughout.

---

## Cohort scope

Eight cohorts. **Schulte-Schrepping and Bernardes are permanently dropped** (Cindy, 2026-09-19).

`Qiu · Cillo · Chongqing · Stephenson · Liu · Wilk · Ren · COMBAT`

---

## Cell-count ladder (bind every number to its level)

Counts are only meaningful with their level named. Do not quote a number without one.

| Level | Meaning |
|---|---|
| source cohort population | all barcodes in the cohort as distributed |
| candidate barcodes | after demultiplexing/whitelist, before cell calling |
| called cells | after cell calling (Liu only; others ship called cells) |
| Gate 1 retained | after frozen basic QC |
| Gate 2 retained | after doublet exclusion — **does not exist yet** |
| benchmark cells | Platelet+Monocyte extracted from Gate 2 retained |
| donors / patients | biological patient, repeated draws collapsed |
| defined-outcome subset | donors with non-null mortality outcome |

### Gate 1 removal, per cohort

Artifact: `05_RESULTS/rebuild_upfront_qc/gate1_qc/<cohort>/audit_per_filter.csv` on HPCC.

| Cohort | Removed by generic rules | Note |
|---|---:|---|
| Ren | 9,264 (0.63%) | |
| Stephenson | 9,868 (1.52%) | |
| COMBAT | 17,656 (2.11%) | |
| Cillo | 2,890 (4.06%) | |
| Chongqing | 42,391 (4.69%) | |
| Qiu | 4,766 (7.42%) | |
| Liu | 161,527 (20.41%) | **superseded** — see below |
| Wilk | 141,344 (57.46%) | **superseded** — see below |

**Artifact hazard.** `gate1_qc/liu/` and `gate1_qc/wilk/` still contain the *generic-rule*
outputs (629,740 and 104,663 cells). Those are **not** the frozen results. The frozen
cohort-specific rules produced different sets, materialized separately:

- **Liu** — barcode-rank *inflection* cell calling: **644,679 / 791,267** candidate barcodes
  (81.5%). Eligible platelet donors fell 49 → 42 at the prespecified ≥20-platelet minimum.
- **Wilk** — platform-calibrated Seq-Well rules: **236,145 / 246,007** retained (4.0% removed;
  platelet removal 76.7% → 2.46%).

Both live in `05_RESULTS/rebuild_upfront_qc/gate2_input/<cohort>/`. Anything reading
`gate1_qc/{liu,wilk}/` is reading a superseded set.

### Assembled benchmark input (pre-rebuild)

`1,236,881` Platelet+Monocyte cells / `665` donors / `80` mortality events, from the
pre-Gate-1 DecontX path. **This will change** once Gate 2 completes and DecontX is refit.

---

## Gate 2 is blocked — two defects found 2026-09-21

An external code audit of the repository plus follow-up verification found:

1. **The expected-rate prior is unsound everywhere, by two different routes.**
   Because `expected_doublet_rate` feeds Scrublet's score calibration, **the scores themselves
   are affected** — these are not merely mis-reported numbers.

   The two Gate 2 entry points obtain the prior differently, so "`prior_counts.csv` is missing"
   does not decide the question. Evidence is the diagnostics themselves
   (`adhoc_gate2_prior_provenance.py`; a capture whose `n_cells_recovered_pre_qc` equals its
   `n_cells_scored_post_qc` had no real prior):

   | Cohort | Path | Captures with a real prior | Cause |
   |---|---|---:|---|
   | Qiu, Cillo, Stephenson | combined | 0 / 190 | `--prior-counts` file absent; `os.path.exists()` fails open |
   | COMBAT | split | 0 / 31 | `decontx_input/combat/split_per_channel` **does not exist**; `n_pre` silently `None` |
   | Wilk | combined | 23 / 23 | prior present |
   | Chongqing, Ren, Liu | — | — | cancelled, timed out, or incomplete at the time; no complete diagnostics available |

   **But Wilk's prior is not correct either.** Both Liu's and Wilk's priors come from the same
   line — `flags.groupby('channel_id').size()` in `materialize_liu_wilk_for_gate2.py` — which
   counts *cells entering Gate 1*, i.e. the already-processed candidate population, not droplets
   recovered from the instrument. For Liu that set is the demuxlet-SNG candidate population.
   So **no prior in use has been shown to be a valid loading count**, and the loading-rate
   formula `0.008 × n_recovered / 1000` is being fed a quantity it was not defined for.

   **Settled by the crosswalk (P0.1, complete).** `gate2_source_crosswalk.py` mapped 790 source
   artifacts onto all 746 Gate 2 captures — every capture has a located source, and
   `n_input_barcodes` is NA for all of them. No deposit holds an unfiltered droplet matrix or a
   Cell Ranger `metrics_summary`.

   **Graded by the T7 evidence layer (complete).** An earlier reading of the above as
   "UNAVAILABLE for all eight cohorts" was too coarse: it conflated *no instrument-level
   recovered count* with *no evidence at all*. The literature survey supplies nominal loading
   for four cohorts, so grading now uses three orthogonal fields
   (`loading_evidence_status` / `upstream_doublet_handling` / `rate_model_status`) in
   `reports/manuscript/revision/literature_audit/T7_doublet_prior_evidence/`:

   - **Loading known** — qiu (15,000/sample), wilk (15,000/array), combat (50,000/channel)
   - **Partial or unmappable** — ren (cohort-level range only), stephenson (98 pools, GEM wells
     not identifiable)
   - **Unavailable** — liu, cillo, chongqing
   - **No cohort is `DIRECT_RESIDUAL_ESTIMATE`**; the build script asserts it.

   Upstream doublet filtering has *already happened* in six of the eight cohorts — every one
   except Qiu and Chongqing, whose `upstream_doublet_handling` is `NONE_REPORTED`. That is
   what makes a residual prior unknown rather than merely missing. Ren is the sharpest case: our Gate 2 input
   is demonstrably post author-Scrublet and post manual removal (summing
   `Ncells_manual_removal` over 284 samples gives 1,462,702, matching our per-capture counts
   sample by sample), so their `expected_doublet_rate = 0.08` — an input prior, not a measured
   removed fraction — must not be reused.

   Capture coverage splits three ways and the kinds must not be merged: **580 mapped**, **23
   `cohort_protocol_inheritance`** (Wilk, `capture_specific=false`), **143
   `unmapped_to_gem_well`** (Stephenson). 580 captures have **capture-specific evidence
   linkage** — linkage, not measurement: a mapped row may hold a measured count, a published
   protocol value, a derived value, a proxy, or an upstream-processing record, and identity
   mapping does not establish a directly measured residual prior.

   T7 also holds the v3 adjudication and two scenario tables:
   `rate-sensitivity-scenarios-v3.csv` (21 rows, `runtime_eligible` false on **21/21**) and
   `pooled-same-donor-pairing-scenarios-v3.csv` (58 rows, false on **58/58**). They document
   possible assumptions and their stage mismatch; they **do not select a rate**. No T7 table is
   a runtime prior configuration, no Gate 2 code may consume one directly or through a wrapper,
   and a numeric value, a readable column or a `middle` label authorizes nothing. Promotion to a
   runtime parameter needs a separate, explicit, auditable decision in `docs/decisions.md`.

   **Gate 2 rate selection remains unresolved.** Nominal loading exists for some cohorts but is
   not the residual prevalence of the matrix Scrublet is given; those are different stages.
   Consequence, per `docs/data-contracts.md`: Gate 2 must run a **rate sensitivity grid**
   (scores recomputed at each rate, since the rate enters calibration), treat the loading rate
   as an uncertain parameter, and never force top-`r` deletion. None of that has been run.

   **All eight cohorts' Gate 2 outputs are quarantined**, under three distinct reasons
   (verified against `05_RESULTS/rebuild_upfront_qc/gate2_scrublet/`; no un-quarantined
   directory remains):

   | Suffix | Cohorts | Reason |
   |---|---|---|
   | `.INVALID_postqc_prior_20260921` | chongqing, cillo, combat, qiu, ren, stephenson | prior silently fell back to post-Gate-1 counts |
   | `.INVALID_prior_definition` | wilk | prior was present, but counts cells entering Gate 1 rather than droplets recovered |
   | `.INVALID_pre_status-enum_pre-checkpoint` | liu | run predates the status enum and checkpoint fingerprinting; reached walltime and was deliberately not rescued |

   Wilk survives the *fallback* defect but not the *definition* defect.

1b. **Split-path summaries are incomplete after a resume.** `gate2_scrublet_split.py` accumulates
   `summary_rows` only for channels processed in the current invocation. Observed: COMBAT's
   diagnostics hold **31 rows against 70 actual captures**.

2. **Silent threshold-detection failure recorded as a successful all-singlet call.**
   `scrublet_filter_combined.py` reads `predicted_doublet` unconditionally after its
   `try/except`. Scanpy warns rather than raising when auto-thresholding fails. Confirmed in
   **1 Wilk capture**: `status='ok'`, `threshold=NaN`, 0 flagged. This violates the frozen rule
   that failed/skipped/unreliable must stay distinct from singlet.

**No cells have been deleted.** Gate 2 remains scoring-only and awaits implementation/config
audit plus separate run authorization; the methodology review itself is complete.

### Withdrawn claim

An earlier analysis in this session concluded Scrublet's automatic threshold was "systematically
conservative," based on it sitting 2.3–3.4× above the expected-rate quantile. **That conclusion
is withdrawn.** `expected_doublet_rate` is a Bayesian calibration prior, not a deletion quota;
the detected fraction should approximate *detectable* doublet fraction × prevalence, not total
prevalence. The comparison does not establish that the threshold is wrong — and six of the
expected rates it rested on were themselves computed wrongly (defect 1).

---

## Cohort runnability classification (working, not a deletion permit)

Conditionally approved 2026-09-22. "Standard candidate" means only *suitable for designing a
scoring experiment*; it confers no rate and no deletion eligibility.

| Cohort | primary_class | secondary_action |
|---|---|---|
| Qiu, Chongqing | `standard_gate2_candidate` | — |
| Cillo | `upstream_doublet_filter_already_applied` | `standard_gate2_candidate`, limited to the residual types HTO did not cover |
| Liu, COMBAT | `upstream_doublet_filter_already_applied` | residual-input diagnostics only; no calibrated same-donor prevalence claim |
| **Ren** | `upstream_doublet_filter_already_applied` | **`residual_scoring_only`** |
| Stephenson, Wilk | `not_eligible_without_additional_evidence` | — |

Ren's input is downstream of the authors' own Scrublet (removed 131,962 cells, 7.6%) and a
further undefined manual removal (136,006, 8.5%). Per-capture counts match
`Ncells_manual_removal` for **284/284** samples — **counts-level, not barcode-level**; no
barcode identity comparison was performed. The authors' `expected_doublet_rate = 0.08` was
their *input prior*, not a measured removed fraction, and is not reusable.

## Sensitivity design — written scoring-only adjudication completed

The owner commissioned automated review to provide the methodological decision on 2026-09-24. The complete
ruling is `reports/manuscript/revision/methodology_review_2026-09-24/REVIEW.md`; it records prior
involvement and is an AI review, not external human sign-off. The original September 23 AI
opinion is unavailable; it was not reconstructed or treated as independently verified.

The approved design is the fixed algorithmic prior stress grid
`0.005, 0.01, 0.02, 0.04, 0.08, 0.12, 0.16, 0.20`. These settings are **not estimates or bounds
of cohort prevalence**. No vendor-to-post-Gate-1 conversion is supported by the current
evidence. `r`, loading-composition `s`, and class-specific survival remain unidentified;
`e in [0,1]` belongs in explanatory models only. Post-demux donor shares have no universal bias
direction as loading proxies.

Qiu/Chongqing are conditional standard-scoring candidates. Cillo/Liu/COMBAT are residual-input
diagnostics: upstream filtering changes the residual doublet mix, so unrestricted simulated
pairs do not establish calibrated same-donor probabilities. Ren remains
`residual_scoring_only`; Stephenson/Wilk remain `not_eligible_without_additional_evidence`.
All scoring units must satisfy the review's raw-count, physical-capture and parent-state
conditions; a cohort label alone does not prove those conditions.

Within a capture/seed, simulations, preprocessing, PCA and neighbours must be fixed across
priors. Rank and AUC invariance then follow from the score's monotone form and are not biological
validation. Threshold/call variation is a reportable sensitivity result, not grounds to choose
a favourable rate or discard a completed stress experiment.

The follow-up `reports/manuscript/revision/methodology_review_2026-09-24/GATE2_PARAMETER_RULING.md`
approves **all eleven methodological parameters** for a Gate 2 scoring-only pilot: simulation
ratio 2.0, UMI subsampling 1.0, log False, PC 30, k = round(0.5*sqrt(post-filter N)), prior SD
0.0, Euclidean distance, variance normalization True, centering True, approximate neighbors
False, neighbor-parent diagnostic False. Allowed values, scope and invalidation are fixed in
that ruling. This does not select Stage 31's PCA dimension.

The methodology mapping is complete. GSM5102900_HC1's complete deposit has passed C1–C3;
its scoring candidate passed the scoped final audit on 2026-09-25. See
`reports/manuscript/revision/gate2_capture_audit_2026-09-24/FINAL_REVIEW_2026-09-25.md`.
The global production slot stays empty; a validated permit supplies only a scoped ruling.
The graph is checked during pilot execution and acceptance, with failure stopping the run;
distance invariance is not an additional gate. **Update (2026-10-11/12): permit, clean deployment
dry-run, and the single authorized scoring pilot are all done** (permit commit `b3b479f`,
deployment `worktrees/b3b479fe...a9617/`, SLURM job `29882644` COMPLETED and accepted — see the
Gate 2 row above for the full acceptance evidence). The global production slot remains empty
(`production_parameter_ruling_written=false`, independently confirmed in the scoring manifest);
no production permit, additional cohort scoring job, or deletion authorization is implied by
anything above beyond that single scoped pilot on GSM5102900_HC1. **Update (2026-10-12, revised
same day per `review:20261012_124351`):** the pilot's own result is interpreted in
`scripts/32_rebuild_upfront_qc/GSM5102900_HC1_PILOT_INTERPRETATION_v1.md`. `n_called` is **not**
monotonic in the prior (291 at prior 0.16, falling to 251 at prior 0.20); **the threshold itself
is also not monotonic** (0.283172 at prior 0.02, falling to 0.188922 at prior 0.04) — Scrublet
fits each prior's threshold independently (`skimage.filters.threshold_minimum` on that prior's own
simulated-score histogram) and the prior is also an input to the per-cell score formula itself, so
nothing in the algorithm requires monotonicity. The 8 priors collapse into **four**, not three,
distinct candidate sets (machine-counted by value, not read off a table): `{0.005,0.01,0.02}`→8
cells, `{0.04,0.08}`→225, `{0.12,0.2}`→251, `{0.16}`→291 on its own. All 8 priors induce the
identical score *rank order* (`all_argsorts_identical=true`) and every candidate set nests in
every other (`all_pairs_nested=true`) — rank order being identical does not mean the score
*values* share one scale across priors, so threshold values from different priors are not
comparable cut points on one axis. This is a reportable sensitivity result about this one capture,
not a deletion rule, not a rate choice, and not evidence about any other capture.
Stage 31 separately requires raw-source, physical-capture and raw-barcode identity closure
before the reviewer can authorize a new dimension-only audit.

The current evidence PDF is
`reports/manuscript/revision/sensitivity_design_review/sensitivity_design_review_packet.pdf`
(12 pages; SHA-256
`1ea5507daf5775fbdf8da1c2482ec0a52e511309b3a663968680ab3602343c5a`). Its evidence remains frozen
at 2026-09-22; the 2026-09-24 review is a separate decision, including a correction to the
fourteen-fold residual-rate claim. The external handover log remains unfilled.

---

## Gate 2 mechanical safety layer (complete 2026-09-22)

The machinery that made the previous runs invalid has been replaced. **This does not unblock
Gate 2** — it makes an unauthorized run impossible rather than silently wrong.

| Module | Guarantees |
|---|---|
| `gate2_runtime_config.py` | the only route a rate can reach Scrublet; missing/ineligible/duplicate/mis-staged/out-of-range all raise; T7 paths are refused by name |
| `gate2_status.py` | six-value per-cell enum, separate from capture execution status; deletion is exactly `call_status == "predicted_doublet"` |
| `gate2_score_core.py` | validates Scrublet's return before accepting calls; a warning-only threshold failure becomes `threshold_failed`, never singlet |
| `gate2_checkpoint.py` | 18-field fingerprint, refuse-on-any-mismatch, atomic temp+fsync+replace, collision-resistant keys, no pickle |
| `gate2_run.py` | one runner for both layouts; diagnostics rebuilt for every capture including restored ones |
| `gate2_diagnostics.py` | score-distribution evidence computed while both arrays are in memory — AUC (simulated = positive class), 101 percentiles, shared-bin histograms, threshold tail fractions |

Both fallbacks that caused the invalid runs are gone, and the four scripts that contained them
(`gate2_scrublet_{score,split}.py`, `scrublet_filter_{combined,split}.py`) now exit with a
refusal rather than running old semantics.

Randomness is explicit: a base seed is required, per-capture seeds are derived deterministically
from base seed + cohort + raw channel, and both are recorded and fingerprinted.

**Stage 31 defect status (fixed in code; old artifacts remain stale).** The old script collapsed
unexamined cells to `False` and undefined patient rates to 0%. The replacement keeps a nullable
`doublet_call`, an explicit `examination_status`, examined-only denominators and NA for patients
with no examined cells; filtering, aggregation, export and statistics all validate the same
status↔call invariant. Threshold metadata is required, batched writes are staged before commit,
and outcome values are excluded from the integration-only analysis path (raw h5ad outcome columns
are acknowledged as incidentally loaded and then removed by a `Sample_ID` allowlist).

One explicitly authorized integration-only technical validation ran as SLURM job **29010326**
and completed in **2:27** with **12.2 GiB MaxRSS**. It did not run mortality analyses. It showed
the material size of the old semantic defect: among 240 patients, old-zero → new-NA for **84**,
old-zero → true measured zero for **1**, old-zero → new-positive for **1** (CV0104), and
old-positive → new-NA for **1**. The v1 call rate is **uncalibrated** for this material and does
not authorize deletion.

The v1 directory
`contamination_audit_v2_fixed_semantics_2026-09-22_v1` is retained as a technical validation
record but is **not in the evidence chain**: the run exposed an implicit/default PCA-dimension
problem. The code now requires `--n-prin-comps` with no default and provides a dimension-audit
preflight that attempts the Scanpy/Scrublet preprocessing measurement without making doublet
calls. A separately authorized dimension-audit run, SLURM **29016446**, subsequently completed
in **1:27** with about **12.24 GiB MaxRSS**: 301/301 samples were measured without error and all
161 samples eligible by the cell-count rule were dimension-feasible for
`n_prin_comps <= 29`. This establishes a computational feasible region only.

The 2026-09-24 written methodology review confirms a blocker upstream of PCA dimension:
Stage 31 first subsamples to at most 60 cells per patient and then runs Scrublet independently by
`Sample_ID`, leaving eligible groups of only 30–60 cells. No fixed `n_prin_comps` can restore the
expression manifold or neighbourhood information removed by that subsampling. Therefore no
dimension is selected and no v2 detection run is authorized. The required redesign is to score
the complete pre-subsample per-capture/sample raw-count matrices first, then map status/calls to
the 60-cell analytic subset by unique capture/barcode identity. Complete inputs mean the
cell-called/basic-QC population with adequate parent-state representation, not all empty
droplets. Sample_ID or sequencing lane may name a capture only where that equivalence is proven.
If those complete inputs or mappings do not exist, the Stage 31
doublet-rate outcome analysis must remain outside the evidence chain. The four legacy artifacts
remain registered stale and untouched.

**Existing September work was found and must be reused, not mistaken for the required input.**
HPCC contains the 2026-09-13 multi-cell-type DecontX build under
`05_RESULTS/bib_revision_tier2/contamination_audit/decontx_multicelltype/`: 11,939 cells,
144 patients and 161 Sample_IDs, comprising 6,729 canonical platelet cells, 2,545 Ren
multi-lineage cells and 2,665 Stephenson multi-lineage cells over 13,065 intersected genes. Its
upstream Ren raw matrix/annotation and Stephenson `with_raw.h5ad`/SDRF remain available under
`01_DATA/raw/monocyte_source/`, together with mapping code and logs.

That assembled object is **not** a full-capture Scrublet input: its platelet arm already came
from the 60-cells-per-patient canonical subset and its non-platelet arm was deliberately sampled
for the DecontX ambient-background question (including patient × lineage caps). It can supply
source loaders and crosswalk evidence, not doublet neighbourhoods. The remaining inventory must
use the actual Stage 31 pre-subsample population — **302 Sample_IDs across nine COVID `Data_NO`
groups** — not the eight Gate 2 cohorts. Of those, 301 survived the 60-cell patient-balanced
subsample and entered the dimension audit; Bernardes sample `004_TA` (16 cells) was the missing
unit. Ren and Stephenson are demonstrated overlapping source groups; Qiu sepsis raw data also
exists but is outside this COVID-only analysis.

Liu is now the strongest additional source candidate: 55 H5 partitions have the same 32,738
Gene Expression features and lane-ID sets matching 55 demuxlet and singlet outputs. All also
contain ADT features, which must be excluded from any RNA-only Scrublet input. The files use a
fixed 30,000-barcode boundary; its mechanism remains unknown, and integer UMI-count provenance is
strongly supported but awaits source metadata. `B*_10xlaneNN` remains a source partition key,
not a proven GEM-well key. Existing upstream-QC code already provides a reproducible
lane+barcode→demuxlet-donor bridge, but the Seurat-suffixed canonical IDs and donor×timepoint
Sample_IDs are not yet reconciled. Liu therefore remains
`pending_capture_identity + pending_cell_level_mapping`, not runnable.

**Two layers, and only the first exists.** `call_status == "predicted_doublet"` marks a
*mechanical candidate*; a final deletion authorization additionally requires cohort eligibility,
an approved runtime rate, a capture reliability criterion and a sensitivity analysis. **No
final-authorization artifact exists** — deliberately not even an all-false column. No cell in
any cohort has deletion authorization.

Vocabulary was corrected twice for the same reason — a name asserting more than the code can
support: `DELETION_ELIGIBLE`/`deletion_mask` → `MECHANICAL_DOUBLET_CANDIDATE_STATUSES`/
`mechanical_candidate_mask`, and `auto_singlet` → `scrublet_noncandidate` (because "singlet"
already meant demuxlet/genotype SNG — one *donor*, not one cell). Both are breaking renames with
no alias. Checkpoint schema is now `gate2-checkpoint/4`.

Interface compatibility with the installed Scanpy 1.11.3 was confirmed by one independent
synthetic smoke test (300×200 integer matrix, verdict COMPATIBLE). The repository suite itself
executes real `sc.pp.scrublet` **0** times; the two figures are about different things.

`method_unreliable` exists in the enum but **nothing assigns it**. No reliability threshold has
been selected; that remains a policy decision.

---

## Gate 2 capture GSM5102900_HC1 — C1–C3 and preflight passed (2026-09-24)

| condition | verdict |
|---|---|
| C1 un-normalised raw counts with source evidence | **PASS** (deposit) |
| C2 physical capture + raw-barcode namespace | **PASS** (uniform `-1` GEM suffix, 7,666 unique 16 nt) |
| C3 complete declared-stage population | **PASS on the deposit (7,666 cells) / FAIL on the container (65 = 0.85%)** |
| PC=30 strict feasibility | **PASS** — `N_eff=7,666`, `G_HVG=1,576`, 30 < 1,576 |
| kNN graph | **UNVERIFIED** — built inside `calculate_doublet_scores`; the score-free preflight cannot reach it |
| scoring permit | **none** — config non-runnable, `PRODUCTION_PARAMETER_RULING = {}` |

Preflight job **29073092** at commit `9aab7be`: COMPLETED 0:0, 2:18, MaxRSS 2.13 GB. Resolved
`n_neighbors=44`, `k_adj=132`, seed `10,137,777`. Evidence archived under
`reports/manuscript/revision/gate2_capture_audit_2026-09-24/run_evidence_29073092/`.

Scope: **one capture**. Nothing here extends to the other five conditional cohorts, whose C1–C3
remain unestablished. The earlier `gate2_cohort_eligibility_2026-09-24/ELIGIBILITY_DRAFT.md`
0-of-6 statement predates this one capture’s acceptance. Scoring the deposit — not the canonical container — is a condition of the C3 pass.

## Reproducibility state

`f0c1070` committed 24 scripts that had been running the pipeline from the HPCC mirror while
never existing in the repository — including every cohort-specific Gate 1 rule and the entire
Gate 2 path. Before that commit the repository could not reproduce the described runs.

The former Gate 2 checkpoint gaps are closed: checkpoints use the schema-4 fingerprint, bind the
runtime config and input identity, refuse mismatches, write atomically with fsync, and have
fresh-vs-resume coverage. Runtime diagnostics also record the relevant environment and
distribution evidence.

Remaining gaps:

- Gate B (Gate 2) has a conditional scoring-only design. **Update (2026-10-11/12):** one capture,
  `GSM5102900_HC1`, now has a reviewed runtime config, a scoped permit, and one completed,
  accepted authorized scoring run (SLURM `29882644`; see the Gate 2 row above). No other
  cohort/capture has a reviewed runtime config or run authorization yet, and no global production
  ruling has been activated (`production_parameter_ruling_written=false`). The stress-grid values
  must never be described as measured prevalence — confirmed non-monotonic across the prior grid
  on the one capture examined so far, see
  `scripts/32_rebuild_upfront_qc/GSM5102900_HC1_PILOT_INTERPRETATION_v1.md`.
- The external handover log is empty. The current complete AI methods ruling is
  `reports/manuscript/revision/methodology_review_2026-09-24/REVIEW.md`; the missing September 23
  original is disclosed and is not a reason to repeat delivery audits.
- Scrublet output is not yet wired into the DecontX/assembly chain;
  `assemble_platelet_monocyte_benchmark_input.py` still reads the original
  `*_decontx_input` directories.
- Stage 31 must first complete a read-only 302-Sample_ID/nine-`Data_NO` source inventory, reusing
  the existing Ren/Stephenson loaders and mappings; then redesign scoring on complete
  pre-subsample capture/sample inputs before any dimension choice or v2 run. All old
  outcome-facing artifacts remain stale.
- Local review/build/assembly edits and methodology documents are present. Git status is the
  source for their commit state; a working-tree document is not proof of a committed release.

---

## Manuscript

Last successful Overleaf compile: **2026-08-06** (`docs/submission-checklist.md`). Superseded by
2026-09-08→11 method/conclusion/attribution edits and by this QC rebuild. **Must be recompiled**
after the rebuilt benchmark numbers land; do not treat the 2026-08-06 state as submission-ready.

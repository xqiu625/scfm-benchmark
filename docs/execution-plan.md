# Execution Plan — 302006

**Current priorities updated 2026-10-12 from the owner's scope decision.** Earlier workstream
records below are retained with their original verification limits. Organized by dependency, not by date. Replaces
`docs/memory-bank/PLAN.md` (last updated 2026-05-27), archived under
`archive/superseded_2026-09-21/`.

`PLAN.md` was not merely stale — it described a roadmap that predates the QC rebuild entirely,
so following it would route work around the actual critical path. It must not be used.

Current state of each item lives in `docs/current-status.md`, not here.

---

## Active plan — scFM biological discovery framework

**2026-10-15 owner clarification:** implement the approved within-COVID contrast as
`{Severe, Critical}` vs `{Mild, Moderate}`, preserving original severity and independent
sampling-location evidence. **ICU-only eligibility for the severe/critical arm is now owner-approved:**
require confirmed ICU location at sampling, not later ICU admission or worst clinical status.
Report included, non-ICU, missing-location and conflicting records before freezing the real input;
resolve visit selection separately and do not reopen either approved decision. After the COVID
pilot, prepare a separate sepsis ICU vs COVID-19 ICU study design with explicit disease criteria,
sampling-time comparability and disease/site/batch identifiability checks. This future direction
does not authorize immediate pooled fitting.

The project no longer pursues a benchmark or mortality leaderboard. The primary scientific output
is cell-state, gene-program and pathway insight supported by expression evidence and validation.

**Use one primary scFM.** Before model-specific implementation, prepare a bounded shortlist of
at most 2–3 candidates using verified checkpoint/species, gene-interpretation interfaces and
deployment evidence. Recommend one to the owner; do not run a full benchmark to make the choice.
Shared framework contracts can proceed in parallel with this read-only selection work. Once chosen,
implement one production adapter. Other scFMs are optional targeted paper controls, not required
framework integrations; cross-model consensus is not an MVP acceptance condition.

1. **Framework foundation — implemented 2026-10-12, corrected twice same day**
   (`scripts/34_scfm_discovery_framework/`): dataset config schema
   (`code/dataset_schema.py`, type-checked, plus a required `counts_slot` field),
   model/checkpoint registry (`code/model_registry.py`, 2 entries from the
   shortlist, neither pinned as the backbone; scPRINT-2's species evidence is
   `checkpoint_metadata_confirmed`, independently re-verified against the
   checkpoint file itself), artifact/provenance schema with
   `planned`/`executed`/`validated` states, full structural AND type validation
   shared by build/write/read, and NaN/Infinity rejection on both write and read
   (`code/artifact_schema.py`), `validate`/`plan` CLI with real matrix-structure
   checks (counts-slot-aware gene axis, config-directory-relative path
   resolution, explicit failure on an unresolvable obs/var index) and reviewable
   plan provenance (`code/cli.py`), GSEA/GO-KEGG provider contracts with no live
   call implemented yet and a finite-score/non-empty-ranking validator
   (`code/enrichment_providers.py`), and a real runnable example
   (`examples/run_example.py`) exercising the whole
   validate→plan→read_manifest chain. Two further small-scope gaps fixed same
   day (`review:20261012_171916`): a JSON numeric overflow (`1e999`
   → `inf`) accepted on read but not write (fixed with a recursive
   non-finite-float check shared by build/write/read), and `validate`
   crashing with a raw `TypeError` on a type-invalid `expression_matrix_path`
   instead of a structured report (fixed by matching `plan`'s fail-fast
   order). **98 tests pass** on HPCC (up from 43 → 78 → 91 → 98 across the
   four correction rounds), including independent out-of-band mutation-kill
   checks (disabling the structural-check, gene-axis-mapping, path-resolution,
   and non-finite-check functions each make their targeted tests genuinely
   fail) plus inline mutation-kill tests for the manifest and gene-ranking
   validators. Old entry points untouched; `plan`/`validate` never submit a
   job, download a weight, or read a full real matrix (bounded `h5py`
   structural reads only).
   **Local gene-set cache + offline GSEA implemented** (`code/gene_set_cache.py`
   imports a local `.gmt` into a versioned, hash-verified cache, never a URL;
   `code/gsea_runner.py`'s `LocalPrerankGSEAProvider` runs the installed
   `gseapy.prerank` against it, with species/namespace consistency checks
   that fail closed and a `no_overlap` state kept distinct from a completed
   run with no significant terms). A second runnable example
   (`examples/run_gsea_example.py`) builds a synthetic gene set with a known
   enrichment direction and runs the real backend end to end, including a
   test that injects a sign flip and confirms the example's own consistency
   check catches it. **131 tests pass total.** This is model-agnostic
   scaffolding — no adapter yet, no live GO/KEGG annotation fetch
   (`AnnotationProvider` is still a stub). The specific production *adapter*
   is the only piece still gated on the backbone confirmation below.
   **GSEA identity/provenance gaps fixed same day** (`review:20261012_220407`
   F1-F5): the installed gseapy was found (by reading its real source) to
   silently case-fold gene IDs for matching, producing a self-contradictory
   result — fixed by pre-filtering every gene set to its exact,
   case-sensitive ranking intersection before the backend call, plus a
   defense-in-depth check on its reported leading edge. Every run now
   freezes the full ranking to `RANKING_SNAPSHOT.json` and binds source/
   dependency identity into the result; `code/gsea_replay.py` replays a run
   from those frozen artifacts alone, in a new process, and compares every
   scientific field. A shared `validate_gsea_result()` now rejects an
   incomplete "complete" record or an impossible probability field at both
   write and read; the gene-set cache's own file path is resolved relative
   to its manifest's directory, not the import-time cwd. All gene sets being
   filtered out by size is now a structured, audited failure, never an
   unhandled exception. **153 tests pass total** (up from 131), including
   two out-of-band mutation-kill checks (reverting the path-resolution and
   size-filter-precheck fixes in scratch mutants) that confirm the new tests
   genuinely fail without them. **This demonstrates the offline GSEA backend
   runs and reproduces correctly on synthetic data — it is not yet a frozen
   evidence package reviewed and accepted for a real dataset**, and does not
   change the status of item 2 below: the production adapter, live GO/KEGG
   fetch, and Docker packaging remain future increments, not immediate next
   steps, pending the owner's backbone and comparison-design decisions.
   **Replay-identity and failure-ordering gaps fixed next day**
   (`review:20261013_114111` R1/R2/C1): a numerically-identical
   replay was not proof it used the same ranking/cache/source identity as
   the original run — a swapped-but-internally-valid ranking snapshot, a
   swapped cache collection/version, and a tampered source-revision field
   were all previously accepted as "REPLAY MATCH". Fixed by comparing the
   replay's reconstructed ranking hash and the cache's full semantic
   identity (species/namespace/version, not just collection name/content
   hash) against what the ORIGINAL result recorded, by requiring every
   numeric run parameter present with the correct type, and by treating a
   source/dependency drift as a strict-replay failure by default (an
   explicit `--allow-identity-drift` flag is the separately-labeled
   cross-version comparison mode). An invalid (NaN/Inf) ranking previously
   made the snapshot writer raise mid-write before validation ever ran,
   leaving a bare `.inprogress` file and no result at all — fixed by
   validating first. A within-ranking case collision (two distinct IDs
   differing only by case) can still trigger the backend's internal
   case-folding even after the F1 pre-filter, which only prevents
   cross-set leakage; this stays correctly fail-closed (confirmed by a new
   direct regression), but the prior round's doc wording overclaimed the
   pre-filter's own scope and has been corrected. **165 tests pass total**
   (up from 153), including two genuine out-of-band scratch mutants
   (reverting the new replay ranking-identity check, and reverting the
   validate-before-snapshot ordering) that reproduce the exact original
   failure modes when deleted. Still the same standing caveat: this closes
   software evidence-chain gaps on synthetic data, not a claim about any
   real dataset or biological finding.
   **Replay-status scoping, replay-label fix, and first PCA/gene-ranking
   round** (`review:20261013_135634` R3/R4 + PCA/ranking
   interfaces): a re-review found the cache-identity check inside
   `gsea_replay.py` only ran for `status="complete"`, so a `"no_overlap"`
   original replayed against a swapped cache (different collection/version,
   also vacuously zero-overlap) still printed "REPLAY MATCH" (R3) — fixed
   by restricting strict replay to `"complete"` originals only, refusing
   any other status outright before the snapshot/cache/backend is ever
   touched. The printed label depended on whether a drift was *observed*
   rather than on whether `--allow-identity-drift` was actually requested
   (R4) — fixed by keying the label on the requested mode itself. **205
   tests pass total** (up from 165), including a real out-of-band HPCC
   scratch mutant confirming a genuinely new bug found this round (a
   bare uncaught `SystemExit` instead of a clean audited `failed` record
   for an oversized/unrecognized counts_slot read — see below) reproduces
   exactly when reverted. This round also implements
   `code/embedding_providers.PCAEmbeddingProvider` (a real
   `scikit-learn`-backed PCA CONTROL baseline, reusing `dataset_schema`/
   `cli.py`'s own input validation, label-free feature selection only,
   persisted hash-checked round-trip-readable artifacts) and
   `code/gene_ranking_providers.py`'s two evidence-stream contracts
   (`ExpressionDERankingProvider` with a donor-blocked-design
   estimability guard; `ModelMaskingRankingProvider`, whose masking
   operation is explicitly a gene-input operation, never a cell-selection
   mask — correcting the prior handoff's own overclaim). **Both ranking
   providers are type/validation/serialization contracts only this
   round** — no real statistical/model backend, no fabricated `"complete"`
   results. Same standing caveat: software-contract and offline-baseline
   work on synthetic data, not a real scFM/PCA comparison or biological
   finding.
   **F1-F4 design-contract and evidence-identity fixes**
   (`review:20261013_210131`): a bounded synthetic probe against
   commit `14a5037` found four gaps in the new PCA/ranking-interface round
   (none in the already-accepted GSEA R3/R4 work). **F1**: a donor-blocked
   estimability check applied unconditionally refused an ordinary unpaired
   3-vs-3-independent-donor design, while one cross-group donor out of
   several was wrongly enough to call a donor-blocked design "estimable" —
   fixed by an explicit, code-recognized `SampleStructure.design_type`
   that decides which check applies, checked against the real donor/group
   input; a `comparison.group_column`/`sample_structure.group_column`
   mismatch is now refused; declared covariates are flagged UNVERIFIED.
   **F2**: a duplicated source gene ID was silently resolved onto
   whichever column came last in the ID->position map — fixed by
   validating the gene axis (non-empty/string/unique) before any
   feature-selection branch, refusing a duplicate outright. **F3**: an
   all-zero counts block produced a NaN `explained_variance_ratio_` that
   reached the npz write before the JSON manifest's own validation caught
   it, leaving an orphaned npz with no result file — fixed by checking
   every array about to be persisted for finiteness up front, plus an
   explicit zero-total-variance refusal and record-then-reraise for
   genuinely unexpected backend exceptions. **F4**: two different
   `perturbable_genes` inputs serialized to an identical
   `ModelMaskingResult` — fixed by binding each result's actual input
   content identity (gene set + hash; group/donor value hashes); the
   shared `GeneRanking` type (used by the already-tested GSEA chain) is
   deliberately left unchanged, and `make_evidence_ranking`'s substring
   check is now documented as a lint-level tripwire only. **226 tests
   pass total** (up from 205: 21 new), with a dedicated mutation-kill test
   per fix. The standalone human+mouse PCA example was re-run end to end
   and produces numerically identical output to the prior round. Same
   standing caveat: software-contract fixes on synthetic data, no real
   dataset, no real DE/model backend.
   **R1-R3 entry-point/evidence-envelope fixes + first pseudobulk-aggregation increment**
   (`review:20261013_222450`): review of commit `0fec319` found F1/F4 only partially
   closed. **R1**: `run()`'s entry point only required ANY compared-group sample, never BOTH
   levels present, and an empty/`None` donor value reached `sorted()` uncaught — fixed by a
   `_validate_group_donor_values` check (length match, non-empty strings, both levels present)
   run before any donor grouping. **R2**: the partially-paired refusal message overclaimed
   mathematical inestimability (disproved with the review's own full-rank 6x5 counterexample) —
   reworded as a stated SUPPORT-SCOPE limit, not a math claim; straddling donors are never
   silently dropped. **R3**: `ExpressionDEResult`/`ModelMaskingResult` carried only a length+hash
   summary with no recoverable input and an identity hash using `json.dumps(..., default=str)` —
   fixed by adding a round-trippable `input_mapping`/`perturbable_genes` field with a
   recomputed-and-compared hash, dedicated `write_*_result`/`read_*_result` functions (validated
   identically at write and read), and a `_safe_json_identity` helper that never silently
   stringifies a non-JSON-native value. Also corrected the three outstanding doc-accuracy notes
   from the prior round (log1p/standardize phase of `input_counts_sha256`, the unsaved
   standardize mean/std, and the dense-read bound message no longer suggesting a
   `cell_ids_subset` workaround that does not work). **New `pseudobulk_providers.py`**: sample x
   (caller-selected) cell-type/state RAW count-sum aggregation — `sample_id_column` is always the
   aggregation unit (never donor/animal ID, so repeat visits never silently merge), conflicting
   per-sample metadata is refused, a counts-conservation check runs, and results share the same
   strict write/read/tamper-rejection contract as the other providers. Computes NO DE statistic,
   p-value, log-fold-change, or `GeneRanking`; `ExpressionDERankingProvider` itself is untouched.
   **266 tests pass total** (up from 226: 18 + 21 new), with a dedicated mutation-kill test per
   new guard. Both the PCA and GSEA standalone examples were re-run end to end on HPCC. Same
   standing caveat: software-contract fixes and a synthetic-data aggregation increment, no real
   dataset, no real DE statistic.
   **F1-F3 pseudobulk-identity/type-validation fixes** (`review:20261014_100436`): review
   of commit `32b22fa` found the npz identity binding and the evidence-envelope's status contract
   only partially closed. **F1**: `sample_rows`/`gene_ids` LABELS carried no binding to the npz
   array's actual CONTENT (reversing them, or editing `n_distinct_donors`/a row's `library_size`/
   `condition_value`, passed unchanged) — fixed with a new `identity_sha256` recompute-and-compare
   field plus `load_pseudobulk_counts` now recomputing each row's library size/each gene's total
   count directly from the loaded array; `dataset_identity` gained `selected_raw_counts_sha256`/
   `selected_cell_metadata_sha256`/`gene_ids_sha256` content hashes (a config hash is not a
   content hash of the expression matrix). **F2**: `cell_type_or_state_label=True` passed the
   entry truthiness check and wrote an orphan npz before the writer's type check rejected it; a
   sample with an all-missing or mixed-missing explicitly-requested `condition`/`visit` value
   either silently wrote `null` or crashed `sorted()` with an uncaught `TypeError` — fixed by
   TYPE-checking every entry argument and refusing (never raising) a missing/non-string value for
   an explicitly-requested column before any grouping, plus splitting `_write_pseudobulk_result`'s
   validation into a pre-npz check. **F3**: the envelope accepted an illegitimate `'complete'`
   record (`ranking={}`, `backend='not_implemented'`) while REJECTING a genuine `'failed'` record
   carrying a diagnostic `None` donor_value — fixed by re-validating `ranking` against
   `GeneRanking.validate()`, refusing complete-status placeholder backends, type-checking
   `MaskingScoreSpec`, and gating the strict per-row/per-field VALUE checks (not the key-presence/
   hash checks) to `status in ('complete', 'not_run')` so a `'failed'` record's raw diagnostic
   input can still be written/read. A new standalone `run_pseudobulk_example.py` demonstrates the
   real interface end to end for both species, including the failed-record read-back path.
   **288 tests pass total** (up from 266: 13 + 9 new). A read-only DE-backend dependency check
   confirmed the clean env is Python 3.10.18 with neither `pydeseq2` nor `rpy2` installed; this
   round did not install/upgrade anything. Same standing caveat: software-contract fixes on
   synthetic data, no real dataset, no real DE statistic.
   **R1-R3 equal-marginals/type-check/failed-record fixes** (`review:20261014_114357`):
   review of commit `0a40d4f` found F1/F2/F3 only partially closed. **R1**: a fixture with equal
   row/column totals (`[[1, 9], [9, 1]]`) could have `sample_rows`/`gene_ids` reversed, or a
   `condition_value` rewritten, with every recomputed total (and a self-consistently rehashed
   `identity_sha256`) still agreeing — a total is a many-to-one projection of content, never a
   content identity. Fixed by saving `sample_ids`/`gene_ids`/a normalized per-row identity
   snapshot INSIDE the npz itself (plain numpy string arrays, covered by `aggregation_sha256`);
   `load_pseudobulk_counts` now compares the manifest's claims DIRECTLY against this independent
   in-npz copy. **R2**: `n_samples`/`n_pseudobulk_samples`/`n_distinct_donors` had the same
   `isinstance(x, int) and ... != ...` shape with no `else` branch, silently skipping the
   consistency check for a non-int value instead of rejecting it — now type-checked
   unconditionally first in all three places. **R3**: the `n_samples == len(input_mapping)` and
   `perturbable_genes` value checks ran unconditionally, which made a genuinely `'failed'`
   length-mismatch/invalid-gene record (the real provider's own diagnostic) unwritable — fixed by
   new `input_lengths`/`input_mapping_buildable` fields (recorded always) gating the length check,
   and gating `perturbable_genes`' value check the same way `masking_score` already was.
   **300 tests pass total** (up from 288: 6 + 6 new); the 11 new tests independently confirmed to
   fail against the pre-fix code for the exact reason each fix addresses. Same standing caveat:
   software-contract fixes on synthetic data, no real dataset, no real DE statistic.

   **P2 fix + DE design-matrix-preparation increment ("B") + DE-backend environment plan**
   (`review:20261014_151747`): review of commit `3131146` confirmed R1-R3 closed and
   found one new gap — `input_mapping_buildable=False` was never tied to `status`, so a
   genuinely `'failed'` record could have only its `status` rewritten to `'not_run'`/`'complete'`
   and still be accepted. Fixed by rejecting that combination directly; independently confirmed
   via a before/after run against the unmodified `3131146` code. **New `de_design_providers.py`**
   (the "B" increment, gated on the P2 fix passing): pseudobulk-to-DE input/design-matrix
   PREPARATION only — consumes an already-validated pseudobulk aggregation plus the SAME
   `ComparisonSpec`/`SampleStructure` interfaces, supports exactly two synthetic designs
   (independent-unpaired, fully-paired), builds a real design matrix/contrast/rank/residual-df,
   aligns by `sample_id` lookup (never position), and refuses covariates/repeated
   samples/partial pairing/non-full-rank/insufficient residual df rather than silently pooling.
   Output status is `'prepared'`/`'failed'` only, deliberately distinct from
   `ExpressionDEResult`'s states. `ExpressionDERankingProvider` itself remains unchanged — still
   no real backend. **New `DE_BACKEND_ENVIRONMENT_PLAN_v1.md`** (plan only, nothing installed):
   recommends pinning `pydeseq2==0.5.2` (not PyPI-latest `0.5.4`, which needs Python >=3.11 vs.
   the env's 3.10.18); corrects a prior round's single-non-interactive-shell "no R" framing after
   confirming via a login shell that R IS available on this cluster; documents a draft
   `pip install --dry-run` command, an LFC-shrinkage plan, and a positive/null synthetic
   verification plan (paired-null preserving pairing structure; no exact-zero requirement on a
   null draw). **320 tests pass total** (up from 300: 3 gene-ranking + 17 de-design new —
   correcting an earlier "+6" gene-ranking mislabel to the actual "+3", `67 → 70`). Same
   standing caveat: synthetic-fixture/design-matrix work only — no real dataset, no real DE
   statistical backend, no new scientific claim.

   **F1/F2/F3 fixes + DE-design standalone example + isolated environment plan**
   (`review:20261014_203915`): review of commit `51d6b75` found two P1 gaps in
   `de_design_providers.py` and one P2 gap in the environment plan. **F1**: `prepare()` took a
   caller-supplied `counts` array independent of the manifest (same-shape reversed/negated
   counts, or a reordered manifest, were accepted) — fixed by having `prepare()` take a
   `pseudobulk_result_path` and call the real `read_pseudobulk_result`/`load_pseudobulk_counts`
   itself; `PreparedDEDesign` now also carries the manifest's path/hash and upstream
   `dataset_identity` forward. **F2**: declared column names were never checked against the
   pseudobulk aggregation's actual `params`, and read-time validation's "estimable" check is
   vacuous for any full-rank design, so a reversed/zero contrast or swapped reference/comparison
   level passed — fixed by a column-name binding check at `prepare()` and a read-time
   recompute-and-compare of every semantic field against the persisted identity. **F3**: the
   environment plan pointed its install draft at the shared `302006_env_clean` and misquoted
   `lfc_shrink`'s signature — fixed with a new isolated-env path (`302006_de_backend_v1`,
   Python 3.11 + `pydeseq2==0.5.4` primary, `0.5.2`/3.10 fallback), a real executed dependency
   dry-run, and a corrected backend-call draft. All seven of the review's reproduction mutations
   are now independently confirmed rejected (previously all accepted). New standalone example
   `examples/run_de_design_example.py` (human unpaired + mouse paired, real
   aggregate→prepare→write→read/load chain). **331 tests pass total** (up from 320: 11 new).
   Same standing caveat applies.
2. **Discovery and pathway chain:** connect state discovery to gene-level evidence and ranked genes,
   local enrichment, versioned GO/KEGG annotation providers and a report showing evidence and stage
   status. Include an expression baseline and sample-level validation hooks. Design a limited
   comparator/ablation study to test the backbone choice and the framework's added value separately;
   do not predeclare the selected scFM universally superior. Prepare and verify Docker
   packaging incrementally; report unbuilt images and untested adapters explicitly.
3. **Human monocyte application:** verify the monocyte definition rather than assuming the existing
   `Monocyte_Myeloid` label is equivalent; confirm eligible inputs, patient/sample identities,
   biological contrast and validation design. Prepare a concrete real-data pilot for owner review.

   **Owner decision recorded 2026-10-14** (`docs/decisions.md`): the human pilot's question is
   severe vs. mild/moderate COVID-19 within one cohort, monocyte cell states/programs — `confirmed`.
   Per `review:20261014_204536`, a feasibility/decision package was produced this round
   (metadata reads and literature checks only, no real DE/model run): `results/
   monocyte_pilot_feasibility_v1.tsv` (real, HPCC-read-only metadata — the monocyte-annotated
   universe is currently exactly {Ren, Stephenson, Qiu-sepsis}; Stephenson `proposed` primary, 90
   patients with mild/moderate/severe/fatal all populated; Ren `proposed` alternate, 53 patients,
   no moderate tier in this slice), `MONOCYTE_PILOT_PREREG_v1.md` (cohort/severity-mapping/subtype/
   repeated-visit/contrast-direction items marked `confirmed`/`proposed`/`pending_owner`),
   `POSITIVE_CONTROL_EVIDENCE_v1.md` (three literature candidates independently re-verified via
   `WebSearch`, 2026-10-14 — Schulte-Schrepping and Combes support an ISG-high-in-mild direction,
   Lee directly contradicts it with an ISG-high-in-severe finding in classical monocytes; no
   direction frozen), and `MODEL_DECISION_v1.md` (restates scPRINT-2's status, proposes a
   masking-score definition that avoids validating against its own target gene set, separates
   target/training/validation data roles). Cohort, severity-value handling (fatal/recovered),
   subtype scope, positive-control direction, checkpoint confirmation, and the validation split
   all remain `pending_owner`.

   **R1-R4 follow-up** (`review:20261014_221155`): **R1** — `write_prepared_de_design`/
   `load_prepared_de_counts` could be called directly (bypassing `prepare()`) with an internally-
   consistent but wrong artifact; fixed via a shared `_source_binding_problems` check that
   re-resolves the claimed upstream pseudobulk source LIVE on both paths — all 8 of the review's
   own reproduction cases now rejected (341 tests pass, up from 331). **R2** — Stephenson's real
   native severity field (`Status_on_day_collection_summary`, obs-only `h5py` extraction from the
   raw object) gives Asymptomatic=9/Mild=23/Moderate=30/Severe=13/Critical=15, superseding the
   earlier mild=32/moderate=28/severe=22 numbers (which came from a different, disagreeing
   cross-object join); confirmed a clean identity join (0 conflicts), 12/90 repeated-visit
   patients (severity-consistent), and Site as a real 3-way confound; Ren's equivalent check was
   attempted and found genuinely blocked locally (no severity field in its raw files). **R3** —
   gave the masking score a concrete `s`/`Δ(g)` contract, corrected "scPRINT-2 standing default"
   to "candidate," froze the real MSigDB Hallmark IFN-alpha gene set (97 genes, 96 covered in
   Stephenson's panel) as `gene_sets/`, and softened the three-paper literature disagreement from
   "direct contradiction" to an acknowledged cross-cohort difference. **R4** — re-ran the
   dependency dry-run with `--ignore-installed` (full 40-package closure), re-verified the
   backend call against `pydeseq2==0.5.4`'s own source, switched deployment to an explicit
   `--prefix`. New `results/monocyte_pilot_feasibility_v2.tsv` supersedes v1's Stephenson row.

   **Incremental follow-up** (`review:20261014_231259`): the counts-identity check's
   `np.allclose` was fixed to exact equality (343 tests pass, up from 341; genuine before/after
   against `546bee0`). The metadata extraction script's categorical decoder had a latent
   missing-code bug (fixed; confirmed it never affected this round's own real findings) and its
   diagnostics are now hard failures. A real further finding: Stephenson's 90 `patient_id`
   strings overcount true biological donors by at least 4 — a direct SDRF
   (`Characteristics[individual]`) crosswalk found 3 hidden same-donor groups invisible to
   string matching, 2 showing a real severity progression between visits. Corrected true-donor
   (n=86) severity: Asymptomatic=9/Mild=23/Moderate=30/Severe=13/Critical=11. New `results/
   monocyte_pilot_feasibility_v3.tsv` and `results/stephenson_donor_identity_roster_v1.tsv`
   supersede v2. `MODEL_DECISION_v1.md`'s masking contract sharpened (separate PCA/scFM fits,
   explicit cell→donor aggregation, one primary statistic chosen in advance, signed vs. absolute
   `Δ(g)` distinguished, pretrain-overlap `unknown`). New `DE_BACKEND_EXECUTION_PACKET_v1.md`
   consolidates the environment evidence for the next phase.

   **Round following `review:20261015_104226`/`_114223.md`:** donor-identity roster and
   summary now generated by one script (86 is the SDRF-confirmed true-donor count outright, not
   policy-pending); author-code investigation (`results/stephenson_author_code_investigation_v1.md`)
   pinned the real GO:0034340 99-gene IFN list used by the authors (distinct from this project's
   frozen 97-gene Hallmark set) and confirmed the Fig.2f/Table 7 severity-collapse rule is not in
   the public repo; `ExpressionDERankingProvider.run_pydeseq2_from_prepared_design(...)`
   implemented (honest `not_run` on the real, confirmed `pydeseq2` `ImportError` — the numerical
   call path itself not yet exercised); `DE_BACKEND_EXECUTION_PACKET_v1.md` made concrete
   (real paths, re-verified hashes, draft SLURM table, one scoped authorization request).
   Separately, Cindy confirmed directly (via automated review, `docs/decisions.md` 2026-10-15) the main
   comparison is now (Severe+Critical, confirmed ICU-at-sampling) vs (Mild+Moderate); implemented
   as `results/stephenson_icu_severity_grouping.py` (ICU location derived only from the real
   `Status_on_day_collection` string) plus a reviewable per-sample/per-donor report.

   **Round following `review:20261015_142903`, same day:** fixed the backend
   consumption entry point so it re-validates a caller-supplied manifest/counts via the SAME
   checks the write/load path already applies (closes a probe's content-tamper/fractional-count/
   flipped-contrast/zeroed-design-matrix findings), and passes `contrast_vector` as an `ndarray`
   (a plain `list` crashes pydeseq2 v0.5.4's own constructor for a 2-column design). Replaced the
   synthetic fixture's constant "null" gene (which distorted size-factor normalization) with a
   7-gene design including 4 stable background genes. Bound pydeseq2's own joblib worker count to
   the real SLURM allocation (`_bounded_n_cpus()` — `OMP_NUM_THREADS` alone does not do this).
   Made the execution packet's install sequence concrete (`de_backend_requirements.lock.txt`,
   `de_backend_phase_wrapper_v1.sh`). Froze the authors' own 99/200/87-symbol IFN/TNF/JAK gene
   lists as `gene_sets/STEPHENSON_AUTHOR_MYELOID_SIGNATURES_v1.json`. Cindy resolved the one
   flagged cross-group donor (`sdrf::C-8821`) directly, into the comparison arm — a one-time,
   hash-pinned decision (`results/stephenson_icu_severity_grouping_resolve_mixed_visits.py`);
   resolved counts: comparison=25, reference=52, excluded-other=9.

   **Round following `review:20261015_150231`, same day:** a probe found the C-8821
   resolution script only trusted the input's own self-reported state, not a frozen approved
   hash, and silently overwrote any pre-existing output — fixed (frozen-hash gate, verified
   sample/group checks, no-overwrite guard; genuine before/after against the unmodified prior
   script reproduces the exact defect, now rejected). Rewrote `de_backend_phase_wrapper_v1.sh`
   to require a verified repo checkout/commit and a real SLURM allocation as preflight (before
   `conda create`), added a separately-locked pytest closure, replaced the human confirmation
   prompt with an automated exact package-set comparison, and required structured (JUnit)
   proof the synthetic test genuinely passed rather than merely exiting 0. **177 tests pass, 1
   skipped** (up from 171/1).

   **Scope clarification, `review:20261015_224111`, same day:** the 25/52/9
   independent-donor comparison view is one specific analysis view, not a global filter — full
   `donor × visit × state` data stays available for scFM discovery and within-donor
   longitudinal work. New metadata-only `results/stephenson_longitudinal_feasibility_
   inventory.py` found: 15 of 86 donors have >1 visit, but only 2 show any severity-tier
   transition and only 1 crosses the comparison/reference arm — too few for a population-level
   progression claim, enough for case description. New `mono_LONGITUDINAL_FEASIBILITY_NOTE_v1.md`
   recommends (A) full-index scFM discovery and (B) case-level description of the 2
   transitioning donors as the two well-supported near-term options (not a new decision — for
   Cindy to choose). `de_design_providers.py`'s existing guards are unchanged. **Next step, not
   yet done:** ask Cindy the one still-open question from the independent-donor view —
   repeated-visit selection WITHIN an already-resolved group (already raised directly per
   `review:20261015_120317` item 5, awaiting reply) — plus, separately, whether to
   proceed with (A)/(B)/both from the new feasibility note.

   **Follow-up review, `review:20261015_230801`, 2026-10-16:** fixed 3 real
   defects (overwrite-on-rerun, an onset-day miscount that silently folded `sdrf::CV0231`'s
   non-numeric `"Healthy"` into "identical," and first/last-only trajectory comparison that
   discarded `sdrf::C-8896`'s middle visit) plus 2 wording overstatements in the note (the
   12-pair breakdown undercounted by 2; §3 wrongly claimed a new DE-design type was needed for
   the simplest 2-visit case, when `paired_or_repeated_within_donor` already covers it). Added
   real evidence: a bounded `obs`-only read of the frozen Stephenson monocyte raw object found
   all 12 symptomatic multi-visit donors share at least one monocyte subtype across every visit
   (`results/stephenson_sample_subtype_cell_counts.py`/`stephenson_longitudinal_monocyte_
   coverage.py`), so subtype coverage is not a limiting factor for options (A)/(B). The note now
   also names option (C) — a minimal, named, NOT-yet-run within-donor DE design the existing
   backend already supports — as a feasibility suggestion only. **Still the same open items:**
   repeated-visit selection WITHIN an already-resolved group, and which of (A)/(B)/(C)/combination
   to pursue — both remain for Cindy, not decided by this round.

   **Second follow-up review, `review:20261016_102834`, 2026-10-16:** accepted the
   real coverage inventory but rejected "subtype coverage is no longer a limiting factor, so a
   paired DE can proceed" as overstated — presence (>0 cells) and usable sample size are different
   claims. Found and fixed 2 more real defects, neither changing any already-committed real
   number: **R1**, `stephenson_sample_subtype_cell_counts.py` indexed categorical codes directly,
   so a `-1` ("missing," per anndata's on-disk format) code would silently wrap to the LAST
   category via Python negative indexing (reproduced on a fixture; fixed with explicit code-range
   validation; the real object has 0 missing/out-of-range codes). **R2**,
   `stephenson_longitudinal_feasibility_inventory.py`'s `order_resolved=False` gated the
   any-adjacent-transition field but not the first-to-last fields or summary counts, so an
   order-unresolved donor could still get a direction label from arbitrary input row order
   (reproduced on a fixture; fixed — unresolved donors now read `"unknown"` and count under a
   separate `"order_unresolved"` key; the real roster has 0 order-unresolved donors). Added a new
   `results/stephenson_monocyte_pilot_design_readiness_table.py` (donor×visit×subtype real counts
   joined with severity/ICU evidence, plus a per-subtype 1/20/50/100-cell sensitivity listing,
   explicitly descriptive, not a new threshold) and confirmed NO QC-admitted monocyte analysis
   input exists yet — only the raw, provisional `stephenson_monocyte_raw.h5ad`. Corrected
   `mono_LONGITUDINAL_FEASIBILITY_NOTE_v1.md` §4(C)'s overstated "14 fully-paired donors (or any
   subset)" framing (3 of those 14 are Asymptomatic, outside the current symptomatic comparison
   scope — 11 symptomatic 2-visit + `sdrf::C-8896`'s 3 visits is the real in-scope population) and
   its implicit pre-selection of `CD83_CD14_mono`. 5 new tests (27 total across the 2 fixed files'
   own suites); real HPCC mutation-kill directly reproduced both defects on isolated fixtures, then
   confirmed byte-identical against the clean source before any further run; **397 passed, 1
   skipped** on the full relevant scfm/stephenson suite (up from 392/1, independently re-run, not
   carried over from the prior round's count). **Still the same open items, unchanged by this
   round:** repeated-visit selection WITHIN an already-resolved group, which of (A)/(B)/(C) to
   pursue, and now also whether a provisional-QC monocyte input is acceptable for a feasibility
   pilot (`MONOCYTE_PILOT_PREREG_v1.md` §3) — all remain for Cindy.

   **Third follow-up review, `review:20261016_120026`, 2026-10-16 — supersedes the
   "awaiting reply"/blocking framing of repeated-visit selection above, and the prior round's
   "no QC-admitted input exists" claim.** Cindy has ALREADY addressed repeated-visit selection in
   this conversation: she asked to check the two groups' sampling-time distribution/overlap first,
   then asked directly why both samples of a repeated-visit donor can't just both be used. **The
   current, correct status is: the independent-donor cross-sectional view's visit-selection rule is
   DEFERRED pending time evidence and the specific analysis design — not an unanswered blocking
   engineering question, and not "Cindy declined early-visit selection" either.** The full
   qualifying-visit index stays available and is NOT blocked by this deferral; this does not
   re-open the already-approved ICU grouping or the `sdrf::C-8821` resolution, and no subtype/
   cutoff/visit-selection rule is decided by this entry. Separately, the prior round's "no
   QC-admitted monocyte input exists" was itself overstated — a repository-only search missed the
   actual HPCC data directories. Real, bounded, read-only checks (cell-id join only, no large
   matrix read) found: Gate 1 QC IS applied for Stephenson (2026-09-21) and DecontX IS complete for
   Stephenson (2026-09-18, mean contamination 13.0%), and BOTH join 100% by `cell_id` to every one
   of the monocyte pilot's 128,503 cells (0 mismatches) — but no decision exists that either one
   constitutes admission for THIS monocyte pilot specifically, and the actual decontaminated counts
   (not just the per-cell contamination score) were not read. Corrected status:
   `qc_admission_not_verified_for_this_pilot` (see `results/stephenson_monocyte_pilot_design_
   readiness_table.py`'s `QC_SOURCE_EVIDENCE`, v2 outputs). **Owner decision, 2026-10-16:**
   Cindy approved provisional raw-count use for the input-to-state-discovery technical pilot
   by replying "允许" to the scoped automated review question. Retain all visits and attach existing QC
   flags/scores; no filtering, corrected-count substitution, or disease/pathway conclusions
   follow from this approval. Prepare the exact execution package now, reusing existing
   adapters and verified QC evidence. Model/checkpoint selection and per-job SLURM authorization
   remain separate; no job is authorized yet. Repeated-visit selection remains DEFERRED and
   does not block preparation. See `docs/decisions.md` and
   `review:20261016_130248`; do not re-ask the resolved input-use question.

   **Execution-packet round, `review:20261016_130032` + `mono_20261016_130248.md`,
   2026-10-16.** Closed the review's §2 wording fixes: `STEPHENSON_QC_SOURCE_CROSSCHECK_v1.md`
   now says HISTORICAL Gate 1/DecontX artifacts exist and are joinable (not that a current
   admission pipeline is complete — Gate 2 for Stephenson is still only an `INVALID`-named
   file); `design_readiness_table.py`'s DecontX commentary no longer calls C1_CD16_mono/
   CD16_mono "coarser/less-specific" or claims they "carry" higher contamination — only that
   their historical DecontX ESTIMATED score is higher, with no subtype preference implied; the
   evidence citation now points to a new, committed, re-runnable script
   (`results/stephenson_qc_source_crosscheck.py` + its real output
   `stephenson_qc_source_crosscheck_v1.json`) instead of prose in a handoff. Delivered
   `STEPHENSON_MONOCYTE_SCPRINT2_PILOT_EXECUTION_PACKET_v1.md`: reuses the proven-working
   `run_scprint2_human_covid.sh`/`scprint2_env` entry point and `model_registry.py`'s
   scPRINT-2 small-v2.ckpt candidate (still a candidate, not confirmed); names the one real
   missing piece (the raw object mixes 192 Antibody-Capture features into its var axis, unlike
   the one proven-working input, which was already Gene-Expression-only) and implements it as
   `results/stephenson_monocyte_scprint2_pilot_prep.py` (GEX-only subset + Gate 1/DecontX
   annotation, never filtering cells or replacing counts) — implemented and tested against a
   synthetic fixture, NOT run against the real 5.7 G object (forbidden on the login node;
   no SLURM authorized). Wrote (not submitted) `run_scprint2_stephenson_monocyte_pilot.sh`.
   Independently re-confirmed via a live CLI call that `scprint2 embed --help` matches the
   flags the proven script already uses, and that the ONE proven real scPRINT-2 output in
   this project used gene-SYMBOL input with `--is_symbol` unset — replicated that exact,
   evidenced flag choice rather than guessing a "corrected" one. **Remaining for the owner**:
   final model/checkpoint confirmation, per-job SLURM authorization, and whether the reused
   2h/128G ceiling (evidenced on a smaller input) suffices for this larger one — named
   UNVERIFIED, not assumed.

   **Owner decision + packet v2, `review:20261016_141242`, 2026-10-16.** Cindy
   replied "ok" — **scPRINT-2 `small-v2.ckpt` confirmed as the FIRST technical-pilot model**
   (final backbone still unconfirmed; UCE not excluded; recorded in `docs/decisions.md`, do not
   ask again). That v1 packet's "replicate the proven call's flags" reasoning was itself wrong:
   independent review found the proven output's producing library VERSION cannot be confirmed
   to match what is installed now, so this round verified every preprocessing claim against
   the ACTUALLY INSTALLED source instead. **R1**: installed `cli.py` writes
   `<--output_filename>_embedded.h5ad`, not the literal path — fixed in the job script. **R2**:
   installed `scdataloader/preprocess.py` merges on Ensembl IDs when `is_symbol=False`; this
   object's var_names are HGNC symbols. A REAL probe against the installed gene table (no model
   weights touched) found `--is_symbol true` maps 19,194/24,737 genes — well over the library's
   own 10,000-gene floor; switched to `true`. **R3**: the library unconditionally filters cells
   under 200 detected genes; a real, bounded, chunked probe of this object's actual sparse
   structure (indices/indptr only, never densified, login-node-safe) found 0/128,503 cells
   below that threshold for THIS input. **R4**: 2 of the packet's own mutation-kill tests were
   self-defeating (their own `assert False` was caught by their own `except AssertionError`)
   and `validate_prepared_adata()` trusted a caller-supplied report instead of re-deriving
   identity from the real objects — both fixed, re-mutation-killed against the real production
   file (6 independent bypass scenarios from the review's probe table now all caught: reversed
   cell order, corrupted counts, empty GEX object, non-finite QC score, a stale "order
   preserved" report, a no-op validator). **R5**: bundled prep→embed→validate into one
   fail-stop job (new `run_scprint2_stephenson_monocyte_pilot.sh`), `--no-requeue` (was
   `--requeue`, conflicting with this round's no-retry scope), a RUNTIME checkpoint SHA256
   check (not just a recorded value), and a real validator
   (`results/stephenson_monocyte_scprint2_pilot_validate.py`, new) as the job's own last step —
   acceptance criteria are now wired into the run path, not only described in a document.
   Also fixed the prep script's overwrite guard (previously checked only the h5ad path, not
   the summary — an independent fixture had proven a pre-existing summary got silently
   overwritten). **R6**: the new `stephenson_qc_source_crosscheck.py` had reintroduced the
   `-1`-categorical-code defect fixed elsewhere in an earlier round — fixed, re-verified
   against real data (128,503/128,503 still joined, 0 mismatches). Delivered
   `STEPHENSON_MONOCYTE_SCPRINT2_PILOT_EXECUTION_PACKET_v2.md` (v1 superseded, left on disk).
   **Known prerequisite gap, not fixed this round**: the HPCC checkout the job script points
   at is still stale (`35131cf`, unrelated uncommitted changes) — must be resolved before any
   real submission. 72 tests re-run this round, all passing. **Still not authorized: SLURM
   submission, install, push, or final backbone selection.**

   **Follow-up review + packet v3, `review:20261016_152556`, 2026-10-16.** Two
   real gaps closed, no new owner decision needed (input scope and first-pilot model stay as
   already confirmed). **Blocker 1**: the output validator took only `expected_n_obs`, so a
   same-count output made of unrelated ids, duplicate ids, or misrouted QC annotations would
   have passed — rewritten to take the frozen prepared-input `.obs` and check exact cell-id
   coverage (no missing/extra/duplicate ids), the real primary `scprint_emb` key specifically
   (the installed `cell_emb.py:214,227-251` writes a UMAP to a SEPARATE `X_scprint_umap` key,
   which the old prefix match wrongly accepted alone), array dtype/shape (string or
   zero-column embeddings now fail), and QC-column correspondence BY ID. **Blocker 2**: the
   job script's `PREP_RESULTS_DIR` pointed at the still-stale live checkout — now points at
   `302006_updategit_deploy_scprint2_pilot`, an independent clean `git` export, verified clean
   and file-hash-identical to the authoring commit `d0aa764`; the stale checkout itself was
   never read, modified, or cleaned. Also fixed: `validate_prepared_adata()`'s bare `assert`
   statements (stripped under `python -O`/`PYTHONOPTIMIZE`; verified by reverting and
   re-running the new subprocess test, confirming the old form really was strippable) and an
   overclaiming mutation-kill test name. Job script now also `cd`s into the per-run directory
   before `scprint2 embed` (the installed library writes intermediate files to a relative
   `data/` path when no logger provides `save_dir`), verifies the three real input files'
   sha256 BEFORE the run (fail closed, not just recorded after), and pins
   `seed_everything`/`batch_size`/`num_workers`/`how`/`max_len` explicitly (confirmed as the
   installed CLI's own current flag names via `jsonargparse`'s `add_class_arguments`
   `nested_key=None` default). 49 tests pass; 2 further mutation-kills performed and reverted
   on the real validator. Delivered `STEPHENSON_MONOCYTE_SCPRINT2_PILOT_EXECUTION_PACKET_v3.md`
   (v1/v2 superseded, left on disk). **The only remaining item before this specific job is a
   per-job SLURM authorization request — no further code gap is open.** Still not authorized:
   SLURM submission, install, push, or final backbone selection.

   **Follow-up review + packet v4, `review:20261016_202903`, 2026-10-16.** Four
   real problems in v3, independently reproduced via the review's own probes, all closed; no
   new owner decision needed. **R1**: `--time=3:00:00` exceeded `short_gpu`'s real
   `MaxTime=02:00:00` (queried directly) — lowered to `--time=2:00:00`; whether 2h is actually
   enough for this larger input stays genuinely unverified. **R2**: the validator's
   "recovered_from_expected_manifest" claim for a missing required column was never backed by
   an actual write-back — the review's probe called `main()` on a real file, reopened the
   output, and found the columns still absent despite a technical PASS; the fake-recovery branch
   is removed, missing is now always a hard failure. **R3**: the validator never checked
   donor/sample identity (only 2 QC columns) — the review's `swapped_patient_id`/
   `swapped_Sample_ID` probes (ids/embedding/QC correct, only the identity value reversed) both
   returned `hard_failures=[]`; `patient_id`/`Sample_ID` (confirmed real columns on
   `stephenson_monocyte_raw.h5ad` via a bounded h5py read — no `visit` column exists on this
   object, so none was invented) are now required and checked by id. **R4**: the job script only
   recorded whatever `HEAD` a clean deployment had — a fast-forward during queue wait would
   still pass; the submitter must now supply `EXPECTED_DEPLOY_COMMIT` via `sbatch --export`,
   checked against the deployment's live `HEAD`. Re-exported the SAME deploy directory to the
   new commit `3de4bca`, verified clean and hash-identical; the stale checkout remained at
   `35131cf`, untouched. 41 tests pass (corrects the prior round's wrong "49" claim — the real
   count for those 2 files was 38 before this round's 3 net-new tests); 2 mutation-kills
   performed and reverted (hash-matched); the review's own probe scenarios re-run directly
   against the fix, now correctly failing/passing. The revision guard was tested against a
   stand-in git directory AND the real redeployed directory. Delivered
   `STEPHENSON_MONOCYTE_SCPRINT2_PILOT_EXECUTION_PACKET_v4.md` (v1/v2/v3 superseded, left on
   disk). **The only remaining item is a per-job SLURM authorization request, pinned to commit
   `3de4bca` — no further code gap is open.** Still not authorized: SLURM submission, install,
   push, or final backbone selection.

After each increment, hand off implemented files, tests and artifacts for review before expanding
scope. The already identified HC2 preflight defects may receive bounded code-only repairs
(`review:20261012_132816`); do not turn them into another deployment/scoring workstream.
No new SLURM submission, push, cell deletion or production Harmony choice is authorized by this plan.

## Earlier QC and benchmark plan — conditional reference, not current scheduling

The sections below preserve the previous dependency plan. They are not a mandate to resume the
benchmark. Apply relevant QC contracts before real biological interpretation; preserve existing
results. Framework development can proceed independently of the full historical QC rebuild.

### Legacy P0 — Unblock Gate 2 for applicable real-data analyses

Real-data analyses requiring a finalized Gate 2 retained set remain gated on it. Framework software
development is not gated on it. Nothing may be deleted from any cohort until Gate 2 is reviewed and frozen.

**P0 items 2–5 are implemented (the mechanical safety layer).** The complete 2026-09-24
methodology ruling is `reports/manuscript/revision/methodology_review_2026-09-24/REVIEW.md`.
It conditionally approves the scoring-only stress design and defines the input and inference
limits. Its follow-up `GATE2_PARAMETER_RULING.md` approves all eleven methodological parameters.
**Update (2026-10-11/12):** a reviewed runtime config and a scoped pilot permit now exist for
GSM5102900_HC1 (permit commit `b3b479f`, deployed and dry-run-passed; see `docs/current-status.md`
and `review:20261011_234753`). The single scoring job itself is tracked separately — see
item 6 below and `docs/current-status.md` for its latest status.

| Item | State |
|---|---|
| 1. expected-rate data contract | schema and validation built; direct residual stress grid conditionally approved for scoring-only; a reviewed runtime config exists for one capture only (`GSM5102900_HC1`, `gate2_scoring_config_GSM5102900_HC1.py`, permit commit `b3b479f`); no config exists yet for any other cohort/capture |
| 2. fail-closed input validation | **done** — `gate2_runtime_config.py` |
| 3. explicit call-status enum | **done** — `gate2_status.py` |
| 4. fingerprint-bound atomic checkpoints | **done** — `gate2_checkpoint.py` |
| 5. fresh-vs-resume equivalence tests | **done** — `tests/test_gate2_safety_layer.py`, 190 tests |
| 6. re-score the eligible cohorts | GSM5102900_HC1 pilot: permit/deployment/dry-run passed and automated review-reviewed; single authorized scoring job ran and completed (SLURM `29882644`), accepted — see `docs/current-status.md`. The other seven cohorts are not re-scored by this; Ren is `residual_scoring_only`, Stephenson and Wilk remain ineligible |
| 7. reliability policy | recorded in `docs/decisions.md`; `method_unreliable` is reserved and unassigned |

These had to land together. Fixing the prior alone and re-running would have produced a second
invalid pass, because threshold failures would still be recorded as singlet calls and stale
checkpoints would still be reused.

1. **Bind the scoring-stage data contract before generating any runtime CSV.**
   Per cohort, define and record: `channel_id`, `n_candidate_barcodes`, `n_called_pre_qc`,
   `n_post_demux_singlets`, `n_gate1_retained`, `n_donors`, the stage at which donor proportions
   are taken, the evidence grade, and the source artifact plus checksum. Distinguish a physical
   capture from a deposited file or analysis channel. The runtime prior comes from the fixed
   algorithmic stress grid; no cell-count column or donor proxy computes it. Loading/selection
   formulas are explanatory evidence only. Contract goes in `docs/data-contracts.md`.

2. **Make every input defect fail loudly.** Stop the capture or cohort; never fall back. Covers:
   `--prior-counts` absent; file missing; channel missing from the prior table; duplicate
   channel; prior below the scored cell count; prior non-integer or non-positive; raw split cell
   file missing; unknown platform. Both current fallbacks go —
   `if ... os.path.exists(...)` and `n_cells_recovered is None else n_cells`. This defect class,
   a missing input that raises nothing, has now produced three separate incidents here.

3. **Replace the boolean doublet call with an explicit status enum.**
   `predicted_doublet | scrublet_noncandidate | too_small | threshold_failed | method_unreliable |
   runtime_failed`. `scrublet_noncandidate` only where auto-thresholding actually succeeded; a keep mask
   is never on its own evidence of a confirmed singlet. `False` must not stand for the last four.
   Reclassify the known Wilk capture (`status='ok'`, `threshold=NaN`) as `threshold_failed`.

4. **Bind checkpoints to a full fingerprint**, and refuse — not warn — on any mismatch:
   schema version, cohort, raw channel string, ordered cell-ID checksum, matrix fingerprint,
   prior row and expected rate, donor proportions, platform, seed, all Scrublet parameters,
   scanpy/scrublet/scikit-image/sklearn versions, code commit, status and diagnostics.
   Also fix: `.npz` written to a same-directory temp then `os.replace`; the `safe` filename
   transform currently maps `A-B` and `A_B` to the same file; validate restored array lengths;
   store barcode IDs (or at least their ordered checksum); apply the same schema to the split
   path's CSV checkpoints; re-aggregate **all** channels' diagnostics on resume, not just the
   newly processed ones.

5. **Prove fresh-vs-resume equivalence** on a small real-structure fixture, comparing by cell ID
   rather than row position: scores, threshold, status, retained mask, filtered matrix,
   diagnostics. Cases to cover: normal capture; too-small capture; threshold-failed capture;
   multiplexed capture; single-donor capture; interruption leaving a half-written temp; a stale
   checkpoint; a changed prior; a changed cell order.

6. **Prepare, audit, then separately authorize the scoring experiment.** Use the fixed direct
   residual expected-rate stress grid `0.005, 0.01, 0.02, 0.04, 0.08, 0.12, 0.16, 0.20` with a
   common seed and frozen failure rules. Fix simulations, preprocessing, PCA and neighbours
   within capture/seed across the grid. The grid is a calibration sensitivity axis, not a
   prevalence interval; rank/AUC stability is not biological validation. Qiu/Chongqing are
   conditional standard-scoring candidates; Cillo/Liu/COMBAT are residual-input diagnostics
   with an unvalidated residual simulation distribution; Ren is residual-scoring-only;
   Stephenson/Wilk are not run. Do not add `e` as a runtime axis or derive priors from vendor
   tables. Apply the eleven values/rules in `GATE2_PARAMETER_RULING.md`; they are no longer open.
   Accept at least one capture's own C1–C3 evidence, then audit the exact config, including
   resolved k, strict PC feasibility, and direct cross-prior PCA/neighbour-graph identity.
   The reviewer then records a small-pilot permit under the owner's latest delegation, naming
   one complete capture, input/config hashes, resource limits and a new output location. Do not
   subsample cells to make the pilot small or activate the production slot before these gates.

7. **Record the reliability policy honestly — do not claim pre-registration.**
   Four cohorts' preliminary outputs have already been seen, so genuine prior blindness is not
   recoverable. The policy is specified in `docs/decisions.md` with the disclosure wording
   agreed there, and applies prospectively to the corrected-prior, version-bound rerun.

8. **Gate 2 review and freeze.** Review material must include: both expected rates and their
   source, threshold placement, observed vs simulated distributions, per-lineage candidate
   counts, status distribution, and the cumulative Gate 1 + Gate 2 donor × lineage retention
   table. Outcome-group retention is reported **once, for disclosure only**, and may not modify
   any frozen rule.

---

## P0b — Stage 31 boolean doublet fallback (fixed; integration-only validation complete)

Found during the Gate 2 read-only audit. It is **not** part of the Gate 2 safety layer. The code
semantics are fixed and one separately authorized integration-only validation has run, but the
legacy artifacts remain stale and the v1 validation outputs are not manuscript evidence.

`scripts/31_donor_level_revision/code/doublet_detection_and_reaudit.py`

| Original line | Historical behaviour |
|---|---|
| 115 | `a.obs["predicted_doublet"] = False` fills **every** cell in the fallback branch |
| 119 | samples with `<30` cells `continue` — stay `False` |
| 124 | per-sample exceptions are caught — those samples stay `False` |
| 127 | `sum(skipna=True)` counts those `False` values as measured zeros |
| 139 | `rate_by_patient = groupby("patient_id")["predicted_doublet"].mean()` — skipped and failed samples contribute **0% doublet rate as if measured** |
| **154** | `a_clean = a[~a.obs["predicted_doublet"].fillna(False)].copy()` — **the actual deletion**, using the prohibited `fillna(False)` |

Same defect class as the Gate 2 one: `False` asserts "examined and not a doublet", when the
truth is "never examined".

**Line 154 was missed in the first registration and is the most consequential**: it is the only
place in the repository that actually removes cells on a doublet call, and it does so with
`.fillna(False)`. The direction is conservative — an unexamined cell becomes `False`, so `~False`
keeps it — but `a_clean` then travels downstream as a "doublet-cleaned" set in which
never-examined cells are indistinguishable from examined singlets. Found during the Gate 2
terminology audit. This historical path has now been replaced; the table is retained to identify
why the legacy artifacts are stale.

**Affected artifacts**: `doublet_scores_per_cell.csv`, `doublet_rate_per_patient.csv`.
**Affected statistics**: doublet rate vs mortality (Spearman ρ), doublet rate vs library depth
(Spearman ρ).

> **Until this is fixed and re-audited, the patient-level 0% doublet statistics are not reliable
> evidence and must not be cited.** The mortality-correlation statistic is the more serious of
> the two, because unexamined samples entering as 0% is not outcome-neutral.

**Semantics FIXED 2026-09-22; integration-only validation performed; v2/outcome reruns not
authorized.** `doublet_call` is now a
nullable boolean plus an explicit `examination_status`; only confirmed doublets are removed and
retained unknown cells stay distinguishable; per-patient rates use examined cells as the
denominator and are `NA` (not 0%) when a patient has none; the Spearman analyses exclude
undefined rates explicitly and report the contributing patient count.

An enforced invariant binds the two fields — `examined` iff the call is True/False, unknown iff
the call is NA — because closing `fillna(False)` alone left a second collapse path: an
`examined` cell with an NA call entered the denominator without contributing to the numerator
(measured: 0.250 where the answer was 1.000). A sample is marked `examined` only after a finite
automatic threshold is read from its own Scrublet metadata **and** its calls and scores validate
— scanpy warns rather than raising when thresholding fails, so field-level checks alone pass
output resting on no threshold. `too_small` comes only from the cell-count rule, and the call
column must be a nullable boolean dtype rather than merely non-NA. The evidence-bound suite
records **188 passed** for the six declared Stage 31 inputs; see
`scripts/31_donor_level_revision/run_evidence/2026-09-22_test_run.txt`.

The four legacy artifacts are registered stale in
`scripts/31_donor_level_revision/STALE_ARTIFACTS.md`, with a copy placed beside them on HPCC as
`STALE_READ_BEFORE_USE.md`; none was overwritten, deleted or moved. One separately authorized
integration-only run (SLURM **29010326**) wrote to
`contamination_audit_v2_fixed_semantics_2026-09-22_v1`, executed no mortality analysis and is
retained only as a technical validation record. It exposed an unresolved PCA-dimension issue, so
v1 does not re-enter the evidence chain. The code now has no default `--n-prin-comps` and provides
a measurement-only dimension preflight. The authorized dimension audit (SLURM **29016446**)
measured 301/301 samples and showed full computational feasibility for the eligible 161 at
`n_prin_comps <= 29`; it did **not** select a dimension.

The 2026-09-24 ruling **conditionally approves workflow redesign and defers dimension selection**.
Stage 31 currently runs
Scrublet only after the patient-balanced cap of 60 cells and then splits by `Sample_ID`, leaving
30–60 cells per eligible group. Computational feasibility at a chosen PC count does not validate
the resulting neighbourhood model. Before any v2 detection run, refactor the analysis to run
Scrublet on complete pre-subsample physical-capture raw-count matrices after documented basic
cell calling/QC, and map its status/calls onto the analytic 60-cell subset afterward by unique
capture/barcode identity. Do not interpret "complete" as including empty droplets, or infer a
physical capture from Sample_ID/lane names alone. First complete the read-only inventory of full
inputs, capture mappings, estimated resources and a new versioned output path. If the complete
inputs/mappings are unavailable, remove the doublet-rate outcome analysis from the evidence
plan. Record every target unit and all mapping/scoring failures; permit no silent exclusion or
zero-filling. Obtain a new dimension audit on the admissible population before choosing a fixed
integer; the 301/161 coverage table remains historical. No Stage 31 detection or outcome job is
authorized by this document update. Under the owner's latest delegation, the reviewer may
authorize a new dimension-only audit only after raw-source, physical-capture and raw-barcode
identity closure on the complete input. Gate 2's PC=30 ruling does not apply to Stage 31.

Do not restart this inventory from zero. The 2026-09-13 multi-cell-type DecontX work already
provides Ren and Stephenson source loaders, raw-source locations, mappings and logs. Its assembled
`multicelltype_raw.h5ad` is deliberately subsampled (canonical 60-cell platelet arm plus sampled
non-platelet lineages) and is therefore **not** the full-capture scoring input. Reuse its mapping
evidence only. The inventory denominator is the Stage 31 pre-subsample COVID population: **302
Sample_IDs across nine `Data_NO` groups**, not the eight Gate 2 cohorts; 301 survived the old
patient-balanced subsample. Produce a per-Sample-ID
crosswalk recording source object, raw-integer/full-gene/multi-lineage status, capture key,
mapping cardinality and cell coverage. Ren and Stephenson start as located-but-not-yet-proven-
complete; the other seven groups remain unknown until joined. Do not write the refactor until
this inventory has passed static review.

For Liu, reuse `filter_liu_demuxlet_singlets.py`, `build_liu_decontx_input.py` and their existing
55-lane artifacts rather than recreating the lane+barcode→donor bridge. Audit the upstream
provenance that those scripts assumed: whether `10xlaneNN` is a GEM well, how the fixed 30,000
barcodes were selected, and how lane+barcode/donor assignments map to the canonical
Seurat-suffixed cell IDs and donor×timepoint Sample_IDs. A raw string intersection between the
two ID formats is not a mapping test; trace the suffix-generation/integration rule first.

---

## P1 — Rebuild the benchmark input

Depends on P0.

1. Exclude confirmed doublets; retain and mark `unknown` where detection was unreliable.
2. Refit DecontX on the QC-retained **full-lineage** set, per capture.
3. Extract Platelet and Monocyte lineages.
4. Re-assemble the benchmark matrix — and actually wire it to Gate 2 output;
   `assemble_platelet_monocyte_benchmark_input.py` still hardcodes the original
   `*_decontx_input` directories, so Gate 2 currently has no effect downstream.
5. Emit the cohort-flow table: every level of the count ladder in `docs/current-status.md`,
   each bound to an artifact path and producing commit.
6. Audit batch metadata for the rebuilt input: cohort/study, platform/chemistry, physical
   capture/library where documented, biological patient, draw/timepoint, lineage and outcome.
   Record missingness and mapping evidence; do not equate patient, Sample_ID or deposited file
   with a technical batch. Tabulate patient counts and outcome availability by batch, including
   batch × disease/outcome overlap and combinations that cannot be disentangled.

---

## P2 — Foundation model benchmark

Depends on P1. Full-scale, not subsampled (author decision, 2026-09-19).

1. Extract embeddings for all five FMs on the rebuilt input.
2. Complete the **batch-effect check requested on 2026-10-04** before interpreting model results:
   - **Done (2026-10-09), `scripts/33_batch_effect_audit/`:** the metadata/dependency-component
     half of this item. Capture-level inventory and donor-coverage crosswalk (with explicit
     `not_assessable` status and alias-normalization for Liu/Chongqing, Cillo's naming gap left
     open, not patched) in `build_batch_metadata_inventory.py`; patient-level `cv_group` x event
     tables in `build_patient_batch_confounding_table.py`, reporting **all 12** multi-donor
     dependency components (not just the 8 that mix fatal+survived donors) plus a found-but-
     unresolved contradiction between `stephenson_pool_composition.csv` (15 pools, 4 donors each)
     and `canonical_donor_labels.csv` (every one of those donors a singleton `cv_group`) — see
     `scripts/33_batch_effect_audit/results/MANIFEST.md`. These tables show CV-dependency
     structure and outcome overlap; they are not themselves a batch-biology association finding.
   - **Done (2026-10-09), implemented and validated on a synthetic fixture, not yet run on real
     cohort data:** `scripts/33_batch_effect_audit/code/run_batch_diagnostics.py` — the actual
     diagnostic runner. Expression/PCA mode (no FM dependency) and an optional FM-embedding mode
     (cell-ID-aligned, missing embeddings reported as `unsupported`, never silently dropped);
     per-lineage patient-balanced diagnostic sampling (seeded, reported eligible/excluded counts,
     never substituted for the full benchmark population); PCA/UMAP; the shared, tested
     `scib_batch_asw` (moved to `scripts/00_shared/code/batch_asw.py` so
     `embedding_quality_option_c.py` and this CLI use one implementation — the second legacy
     caller in `scripts/12_umap_visualization/` now raises if called, rather than silently
     returning its old non-scIB-correct metric); a neighborhood-batch-mixing diagnostic reported
     alongside its expected-under-random-mixing reference value; an opt-in real Harmony
     correction (harmonypy, explicitly labeled transductive/exploratory, never an out-of-sample
     transform — proven by `tests/test_evaluation_boundaries.py`); and a structured
     `manifest.json`/`metrics_long.csv` with complete/partial/failed status per lineage ×
     embedding, now enforced (a requested-but-failed piece of work fails the whole run and its
     exit code, not just an inner status field no caller checks). Validated end-to-end twice
     (write → compute → write → reload) on `make_synthetic_batch_fixture.py`'s output (known
     batch shift + known bio signal, 2 lineages, 4 cohorts, repeated patient draws, patient-level
     outcome now assigned once and reused across lineages, not re-derived per lineage): the
     injected batch effect was correctly detected and measurably reduced by correction in all 8
     lineage×embedding×scope combinations, with the latest run (`synthetic_fixture_v3`) reaching
     strict `status: "complete"`; see `docs/current-status.md` for the actual numbers and the
     honest bio-conservation caveat. Both samplers (this CLI and the real-data preparer below) now
     use a stable, cross-process-reproducible seed derivation
     (`scripts/00_shared/code/stable_seed.py`) instead of Python's per-process-randomized
     `hash()`, confirmed by a committed test. 57 tests passing across
     `scripts/33_batch_effect_audit/` + the new shared modules, run together on HPCC.
     **Running this against real cohort data is still blocked** on resolved Gate-2 input (P0) and
     item 1 above (5-FM embeddings) for the full-scale version; a bounded, two-allocation
     provisional run against the existing pre-Gate-2 input is specified with regenerated real
     numbers (`scripts/33_batch_effect_audit/results/PROVISIONAL_JOB_SPEC_v2.md` — v1 is
     superseded, kept for audit trail) but not submitted.
     **Correction (2026-10-11): the sentence above is stale for Job 1.** `prepare_provisional_
     diagnostic_input.py` (Job 1 of the two-step plan) WAS submitted and completed
     2026-10-10 (SLURM `29804463`, 24,841 cells × 11,709 genes) — see `docs/current-status.md`'s
     "Batch-effect check" row and `docs/experiment-log.md`. Only Job 2
     (`run_batch_diagnostics.py`, the actual diagnostics run on Job 1's output) remains
     not submitted, now with a read-only preparation packet
     (`scripts/33_batch_effect_audit/results/job2_preparation_packet_v1.md`). The full-scale
     5-FM version above is still blocked as described.
     **Correction (2026-10-11, later same day): the sentence above is now stale for Job 2 too.**
     User-authorized Job 2 (`run_batch_diagnostics.py`, via `job2_preparation_packet_v2.md`) WAS
     submitted once and completed (SLURM `29868753`, COMPLETED, ExitCode `0:0`). See
     `docs/current-status.md`'s "Batch-effect check" row addendum, `docs/experiment-log.md`, and
     `scripts/33_batch_effect_audit/results/JOB2_PROVISIONAL_INTERPRETATION_v1.md` for the
     per-lineage result (mixing improved both lineages; outcome-stratified batch ASW improved in
     Monocyte_Myeloid but slightly worsened in Platelet). This is a pre-Gate-2 provisional
     exploratory run on Job 1's existing (capped, pre-Gate-2) input — it does not satisfy or
     advance the full-scale 5-FM version above, which remains blocked as described.
   - Examine expression/PCA and each of the five FM embeddings on matched cells, separately
     for Platelet and Monocyte. Plot PCA/UMAP coloured by cohort, platform, capture, patient,
     disease/outcome and relevant QC covariates; preserve input version and plotting parameters.
   - Quantify batch association in the original PCA/embedding space with batch silhouette and
     a neighbourhood-mixing diagnostic where the design supports
     them. Use patient-balanced sampling for diagnostics, report patient/cell counts, and
     document metric parameters and unavailable comparisons. Diagnostic sampling does not change
     the full-scale benchmark.
   - Assess within comparable biological strata; report cohort association separately from
     technical-batch evidence. Flag complete batch–biology confounding as non-identifiable.
     UMAP appearance alone does not establish presence or absence of a technical batch effect.
   - Deliver stratified plots, quantitative diagnostics and a short limitations report, building
     on the metadata tables above. Legacy embedding-integrity checks do not satisfy this item.
3. LOCO runs and patient-grouped CV; assess their results alongside the batch diagnostics.
   If correction is warranted, specify a separate sensitivity comparison that preserves
   biological signal and fits learned preprocessing within training folds; document handling
   of unseen cohorts. Batch inspection does not automatically imply batch correction.
4. Statistics with patient-level grouping.
5. Baseline ladder.

---

## P3 — Manuscript

Depends on P2.

1. Integrate rebuilt numbers; update every figure and table.
2. Rewrite the QC section to describe the frozen two-gate protocol and its
   platform-calibrated thresholds.
3. Update limitations — including the doublet-detection limits found in P0 (homotypic doublets
   undetectable; low-RNA platelets systematically harder to call).
4. Re-verify claim traceability (`docs/claim-traceability.md`).
5. Recompile. The 2026-08-06 compile is stale.

---

## P4 — Submission logistics

Zenodo DOI · figure upload · author bios · ORCID · cover letter · AI/LLM disclosure ·
portal submission. Tracked day-to-day in `docs/submission-checklist.md`.

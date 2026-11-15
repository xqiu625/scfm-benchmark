# Decisions

Why important structural and analytical choices were made. Append-only.

## 2026-06-09 — LoRA "+0.006 null" is capacity-confounded; within-backbone effect is +0.113 (roadmap J)

**Context.** The manuscript reports Geneformer LoRA domain adaptation as a null (+0.006 AUC,
abstract / Conclusion 6 / `tab:domainadapt`), "refuting the claim that domain adaptation closes
the gap." Reviewers (R1 W8, R2 M4) flagged that this compares LoRA-fine-tuned Geneformer-**V1-10M**
(0.705) against zero-shot Geneformer-**V2-316M** (0.699) — a 30× larger *different* backbone, so
capacity is confounded with adaptation.

**Decision (author).** Build + run the missing within-backbone control: zero-shot Geneformer-V1-10M
(`scripts/30_review_response_stats/code/geneformer_v1_zeroshot.py`, job 25425050), reusing the LoRA
script's exact data loader / tokenizer / cohort / labels / 5-fold patient CV but freezing the base
encoder and mean-pooling (no LoRA, no head).

**Outcome.** Zero-shot V1-10M = **0.592 ± 0.094**. So within the same backbone, **LoRA gives +0.113**
(0.592 → 0.705), not +0.006. The "domain adaptation is a null" framing is wrong as stated. The honest
revised story is *more* interesting: LoRA-adapting the small V1-10M (+0.113) brings it level with the
30× larger zero-shot V2-316M (0.705 vs 0.699) — adaptation **is** effective, it just doesn't beat
using a bigger pretrained model. `tab:domainadapt` / abstract / Conclusion 6 must be revised; the
C003 "+0.006" sub-claim is marked `inconsistent` in `claim-traceability.md`. This is the one roadmap
analysis that changes a conclusion rather than softening a claim.

## 2026-06-09 — Donor-clustered re-analysis in response to simulated peer review (stage 30)

**Context.** A 5-reviewer simulated peer review (academic-paper-reviewer, decision: Major
Revision) converged on one CRITICAL theme: significance was overstated because DeLong ran on
cell-weighted OOF predictions (~42k correlated cells from 322 donors treated as independent),
inflating SE ~√(cells/donor)≈11× so any ΔAUC reads "p<10⁻³". Reviewers predicted most pairs
would not survive donor-level inference.

**Decision (author).** Run the roadmap P0 statistical re-analyses with the **donor** as the
unit of inference before any claim re-wording: donor cluster bootstrap (B=5000) for per-FM AUC
CIs and pairwise ΔAUC CIs; fold-mean t-CIs + Nadeau-Bengio corrected paired test; BH-FDR;
cross-disease CI overlap; severe-label sensitivity. New stage `scripts/30_review_response_stats/`
(local, on the scp'd OOF preds — no retraining). OOF prediction CSVs committed so a reviewer can
independently recompute (addresses the reproducibility ask).

**Outcome — the paper gets *stronger* where it counts and *honestly softer* where it overreached.**
- **KEEP "STATE leads":** STATE significantly beats all four other FMs on within-disease OOF AUC
  under donor clustering **after BH-FDR** (q≤0.017) — the headline survives proper inference.
- **SOFTEN middle-pack:** scPRINT-2≈Geneformer≈scGPT are statistically indistinguishable
  (donor p 0.37–0.58); 4/10 cell-DeLong pairs lose significance; stop ordering scGPT above UCE.
- **SOFTEN per-fold ranking:** fold-mean CIs overlap; Nadeau-Bengio finds only STATE>UCE (p=0.021).
- **SOFTEN "only bidirectional generalizer":** reverse-transfer CI for STATE overlaps UCE/Geneformer
  (separates only from scPRINT-2). The scGPT 0.352 anti-correlation IS robust (keep it).
- **Label-robust:** severe-excluded raises every FM +0.05 but preserves STATE's lead.

These results back the manuscript revision (claim softening), not a retraction — the experiments
and the central finding are sound; the inference and the over-precise rankings needed alignment.

## 2026-06-08 — Revive scGPT Option C (reverses the 2026-05-01 drop)

**Context.** On 2026-05-01 scGPT was dropped from the Option-C canonical-input set (tooling
pain — env/CUDA/GLIBCXX failures; `step2c` excludes it), and the manuscript shows scGPT `—` for
the within-disease Option-C tables with a "rerun pending" note. The cross-disease scGPT 0.352 used
the legacy embedding.

**Decision (author, 2026-06-08).** Generate the scGPT Option-C embedding properly and fill the `—`.
Path that worked (`scripts/foundation_models/scGPT/embed_scgpt_option_c.py` + `submit_scgpt_option_c.sh`):
run `scgpt.tasks.embed_data` inside `07_CONTAINERS/scgpt_0.1.7.sif` on the canonical raw-counts input,
writing `obsm['X_scgpt']`. Three documented failure modes had to be cleared:
1. `~/.local` pyarrow needing GLIBC_2.28 vs the container's older libc → `PYTHONNOUSERSITE=1`.
2. container `embed_data` asserts `cell_type_key` present → add constant `cell_type='platelet'`.
3. autocast dtype mismatch in the value-embedding Linear → cast input X to float32.

The result (47977×512) **passes the embedding sanity check** (cos 0.892, per-dim std 0.0143,
0 zero-var) — genuinely discriminative, unlike the degenerate mouse scGPT (which stays dropped).
scGPT registered in stage 24 `OBSM_KEYS` + stage 27 `EMBEDDINGS`; both benchmarks rerun to fill
the within-disease tables. The legacy `step2c` still excludes scGPT (the new path is a standalone
container job, not a step2c branch) — documented here so the exclusion isn't read as "still dropped".

## 2026-06-07 — Keep C003 = 0.705; treat the cold-retrain gap as a traceability note

**Context.** The HPCC rerun campaign established project-level reproducibility:
C004 (full-cohort) and C005 (adaptation bootstrap) reproduce **byte-identical** to the
committed evidence on a from-scratch rerun. The C003 Geneformer-LoRA number, however,
did **not** fully reproduce on a cold from-scratch retrain: mean AUC **0.675 ± 0.095** vs
the committed **0.705 ± 0.135**. The gap is isolated to **fold 0** (cold 0.763 @ best-epoch 2
vs cached 0.916 @ epoch 12); folds 1–4 are bit-identical. It is a localized deep-learning
training instability (best-epoch selection on a noisy per-fold validation curve, sensitive to
the short_gpu preemption/resume boundary), not a pipeline defect.

**Decision (author).** For submission, **keep the reported 0.705** and retain the cold-retrain
result as a **traceability note** in `claim-traceability.md`. Do **not** retract the number,
do **not** redo the benchmark, and do **not** rewrite the paper over a 0.705→0.675 fold-0
fluctuation. The qualitative conclusion (LoRA gives marginal-to-no benefit vs zero-shot ~0.699)
is unchanged — if anything 0.675 strengthens it.

**Why this is defensible.** (1) Project reproducibility is demonstrated by the two deterministic
pipelines (C004/C005, byte-identical). (2) The C003 variance is documented, bounded (one fold),
and explained. (3) A reviewer rerunning from scratch would see the same fold-0 behavior, now
disclosed. The cached JSON + fold checkpoints are preserved (`results_mortality_lora.cached_20260606.json`,
`05_RESULTS/checkpoints/geneformer_peft/_cold_repro_backup_20260606/`) so the 0.705 is fully traceable.

## 2026-06-06 — Make the embedding sanity-check variance-aware (Geneformer false positive)

**Context.** During the HPCC rerun of C004, `primary_benchmark_option_c.py` dropped
**all Geneformer rows** — `scripts/00_shared/code/embedding_sanity.py` failed Geneformer
on `cell_pair_cos_sim_mean = 0.9519 > 0.95`. But Geneformer is *not* collapsed: per-dim
std mean 0.1273 (vs TF-Metazoa collapse 3.4e-5), 0/1152 zero-var dims, and committed
mortality AUC 0.699–0.738 with real per-fold spread (0.58–0.85). High cosine similarity is
expected biology here — Geneformer's rank-value tokenization is dominated by a shared
ribosomal signature on platelets (the manuscript says so). The 0.95 threshold, calibrated on
UCE (0.758) and TF-Metazoa (0.985), over-triggered.

**Decision.** The standalone high-cosine criterion now **only hard-fails when it co-occurs
with collapsed per-dim variance** (`per_dim_std_mean < cos_sim_discriminative_std = 5e-3`).
High cosine similarity *with* discriminative per-dim spread is demoted to a **WARNING**
(log + "verify downstream AUC != 0.5"), not a failure. TF-Metazoa's true collapse still
fails — both via this combined gate and independently via the existing `per_dim_std_mean <
1e-4` and zero-var-fraction criteria (those are unchanged).

**Why this is safe.** True collapse (the bug the guard exists for) gives AUC 0.5±0.0 *and*
near-zero per-dim variance; it is still caught. The change only stops penalizing FMs that are
legitimately high-similarity yet discriminative. Verified on HPCC: real Geneformer → PASS+warn;
synthetic collapse (cos≈1, std 1e-5) → FAIL; synthetic high-sim+discriminative (cos 0.99, std
0.40) → PASS+warn. The committed manuscript Geneformer numbers were always valid; this only
lets reruns regenerate them. See `experiment-log.md` (2026-06-06).

## 2026-06-06 — Conform repo to Universal Data Science Pipeline Architecture

**Decision.** Restructured this git mirror to follow
`UNIVERSAL_DATA_SCIENCE_PIPELINE_ARCHITECTURE.md` ("full physical conform").

**What moved.**
- `analyses/NN_*` → `scripts/NN_*` (numbered stage folders are now under `scripts/`,
  matching the template's `scripts/NN_stage/` convention).
- `analyses/README.md` → `scripts/README.md`; `analyses/MANUSCRIPT_MAPPING.md` →
  `scripts/MANUSCRIPT_MAPPING.md`.
- Old `scripts/README.md` (tooling overview) → `docs/tooling-reference.md`.
- Loose root scripts (`download_results_from_hpcc.sh`, `inspect_cell_indices.py`,
  `generate_acc17_21_figures.py`) → `scripts/00_setup/code/`.
- `manuscript/` → `reports/manuscript/`.

**What was added.**
- `configs/{defaults,hpcc,local,cloud}.yaml`
- `docs/{pipeline_index,data_contracts,claim_traceability,figure_table_inventory,experiment_log,audit_report,decisions}.md`
- `reports/audit/` audit report set; `scripts/99_audit/` audit stage scaffold.
- `data/{raw,interim,processed,external}/`, top-level `results/`, `figures/`,
  `notebooks/`, `logs/` (gitkeep placeholders).
- `src/scfm_benchmark/{io,preprocessing,features,modeling,evaluation,validation,visualization,reporting,submission}/`
  empty subpackages for promoting reusable logic out of `scripts/*/code`.
- `environment.yml`, `requirements.txt`.

**Decisions within the conform.**
1. **Kept ACC numbering, did not renumber to template 00–10.** The accession
   numbers (06, 17, 24, 29, …) are referenced directly by the manuscript via
   `MANUSCRIPT_MAPPING.md` and `tab:*` labels. Renumbering would break those
   references right before submission. `docs/pipeline-index.md` maps each ACC
   stage onto the architecture's generic roles instead.
2. **`scripts/foundation_models|finetuning|tracking|agent|cross_species` stay under
   `scripts/`** as named reusable runners rather than being forced into a single
   numbered stage — they are runnable workflow tooling, which `scripts/` is for.
3. **All path references rewritten to `scripts/` — repo is the single source of truth.**
   Every in-repo reference (relative and absolute, including
   `/bigdata/godziklab/shared/Xinru/302006/analyses/...` defaults) was converted
   `analyses/ → scripts/` so the repo is internally consistent. The one exception is
   `legacy/analyses/`, which was deliberately preserved (those stages were never
   moved). **Consequence:** the HPCC working copy (`302006/`) must be renamed
   `analyses/ → scripts/` to match, or the absolute-path defaults will point at a
   non-existent dir on the next run. See "Migration follow-ups".
4. **`legacy/` kept at top level** (not folded into `archive/`) because `scripts/README.md`
   references it as a status tier; documented here instead of moving.

**Migration follow-ups (not done automatically).**
- **Required:** on the HPCC working copy (`302006/`), rename `analyses/` → `scripts/`
  so the rewritten absolute-path defaults resolve. Any external job scripts that
  `cd analyses/...` must be updated too.
- Backfill `docs/claim-traceability.md` and `docs/figure-table-inventory.md` with
  real evidence filenames + manifests, then run `scripts/99_audit`.
- Generate per-stage `results/manifest.json` going forward (see architecture §9).

## 2026-06-06 — Centralize the base path in one resolver (opportunistic migration)

**Decision.** Added `src/scfm_benchmark/config.py` as the single source of truth for
`project_root` / `data_root` / `logs_dir` and merged run parameters. It loads
`configs/defaults.yaml` overlaid by the env config (`hpcc|local|cloud`, auto-detected
via `/bigdata` presence or `SCFM_ENV`), with `PROJECT_ROOT` / `SCFM_DATA_ROOT` /
`SCFM_LOGS_DIR` env-var overrides. Code reads paths via
`from scfm_benchmark.config import paths, settings`.

**Why this and not a mass rewrite.** The automated audit flagged ~130 `portability`
majors — hardcoded `/bigdata/.../302006` literals across ~169 stage scripts. Those
scripts run on HPCC and cannot be validated locally (memory: *302006 no local
validation*); mechanically rewriting all of them risks breaking live SLURM jobs right
before the BIB submission, for low real-world benefit (the path is real on the cluster).

**Consequence / scope.**
- The numbered `scripts/NN_*/code/*` keep their hardcoded `/bigdata` absolutes as a
  working fallback. They are **migrated opportunistically**: the next time a stage
  script is edited for analysis reasons, switch its paths to `paths.resolve(...)`.
- The audit's `portability` count (130: 122 in `scripts/`, 8 in `legacy/`) will stay high
  until migration completes — this is **expected**. (The auditor's regex only matches
  *quoted* absolute paths in code, so the unquoted `configs/*.yaml` values are correctly
  **not** flagged.) Track real progress by "stages migrated to `config.paths`", not by the
  raw audit number. This is recorded in
  `../archive/superseded_2026-09-21/docs_legacy/audit_report_2026-06.md`.
- `configs/hpcc.yaml` remains the authoritative declaration of the cluster base path;
  `config.py` is the loader, not a second source of truth.

## 2026-09-21 — Durable decisions consolidated from the archived memory-bank

**Context.** An external documentation audit found `docs/memory-bank/` had become a mixed store
of expired plans, mutually-overturned analyses, transient HPCC job state and durable facts — and
that entry-point docs still named it authoritative. The directory was archived to
`archive/superseded_2026-09-21/`. The decisions below were extracted first, because they remain
binding; intermediate judgements that later results overturned were deliberately not carried
over.

**Scope and identity**
- Five-FM primary set: STATE, UCE, Geneformer, scPRINT-2, scGPT.
- Mouse work belongs to **302035**, not 302006.
- Eight cohorts; Schulte-Schrepping and Bernardes permanently dropped (2026-09-19).
- Full-scale extraction, not the 60-cells/donor subsample (author override of the 2026-09-19
  subsampling decision).

**Analysis**
- CV grouping key is **biological patient**, not draw or `Sample_ID`.
- No synthetic fallback data, ever — migrated here from `PLAN.md:248` so it survives archival.

**QC (frozen; see `scripts/32_rebuild_upfront_qc/qc_checklist.md`)**
- Thresholds are calibrated to platform and upstream input state, not harmonized to a single
  number. The aim is comparable *retained-cell quality*, not comparable deletion rates.
- QC is outcome-blind: thresholds are never tuned on mortality labels or predictive performance,
  and outcome-group retention is disclosure only.
- Retention is a trigger for investigation, never a target to reverse-engineer thresholds from.
- `skipped` / `failed` / `unreliable` are distinct statuses and must never be folded into
  `singlet`.
- The Scrublet expected doublet rate is derived from pre-QC recovered counts per capture — never
  from the merged matrix or post-QC cell counts. Residual rate after demultiplexing is
  `r·s/((1−r)+r·s)` with `s = Σ p_d²` over observed donor shares; `r/n_donors` is wrong.
- A fixed 500-counts / 200-genes floor must not be re-packaged as cell calling.

**Still open — deliberately not decided here**
The final doublet-handling scheme. The audit recommends status-enum handling with unreliable
captures retained and marked `unknown`, and warns that the reliability criterion carries
post-hoc risk because four cohorts' Gate 2 results have already been seen. Whoever sets that
criterion must record, here, whether it was chosen with knowledge of those results.

## 2026-09-21 — Withdrawn: "Scrublet's automatic threshold is systematically conservative"

Asserted during Gate 2 review on the basis that the automatic threshold sat at the 99.8th
percentile of observed scores and 2.3–3.4× above the expected-rate quantile in all four
completed cohorts. **Withdrawn.** `expected_doublet_rate` is a Bayesian calibration prior, not a
deletion quota; the detected fraction approximates *detectable* doublet fraction × prevalence,
not total prevalence, and Scrublet cannot detect homotypic doublets by construction. The
comparison therefore does not establish that the threshold is wrong. Compounding this, six of
the eight cohorts' expected rates were themselves computed from post-QC counts (see
`docs/current-status.md`), so the quantity being compared against was also wrong.

Recorded because the A/B/C doublet-handling options proposed on top of this claim were
circulated before it was withdrawn.

## 2026-09-21 — Gate 2 reliability policy (two layers), and its disclosure wording

**Disclosure, to be carried into the manuscript verbatim:**

> Reliability policy was specified after exploratory inspection of four preliminary cohort
> outputs but before the corrected-prior, version-bound eight-cohort rerun, and without
> reference to mortality outcomes.

This is **not** a pre-registration and must never be described as one. Genuine prior blindness
was lost once the Qiu / Cillo / Stephenson / Wilk Gate 2 summaries were inspected. An AI
assistant participated in that inspection and in drafting this policy, so it is likewise not an
independent adjudicator; it can help formalize criteria, nothing more.

### Layer 1 — mechanical validity (frozen now; unaffected by the four cohorts' results)

Objective states requiring no subjective threshold:

`too_small` · `runtime_failed` · `threshold_failed` · `valid_auto_threshold`

`valid_auto_threshold` requires all of: threshold present and finite; observed scores of correct
length and finite; simulated scores present and finite; predicted mask of correct length;
checkpoint fingerprint matching in every field; input passing the raw-count contract.

The first three retain their cells and are marked `unknown`. None may be written as
`False` / `singlet`.

### Layer 2 — method reliability (procedure specified now; cutoff deferred)

Scrublet publishes no validated universal cutoff, so inventing one (`AUC >= 0.7`,
`sim_median > obs_p95`, …) would be arbitrary. What is fixed now is the *evaluation procedure*,
not a number: multi-seed threshold success rate; threshold distribution; call-set Jaccard across
seeds; simulated recall; observed/simulated tail enrichment; histogram modality; lineage-pair
composition; spatial concentration of calls in expression space.

**Handling rule:**
- mechanical failure → retain, mark `unknown`;
- mechanically valid but reliability insufficient → retain, mark `unreliable`;
- mechanically valid and reliability adequate → automatic calls may enter the primary analysis;
- top-`r` expected-rate calling is a **sensitivity analysis only**, never the primary deletion
  rule;
- mortality/outcome plays no part in any reliability judgement.

If a binary reliable/unreliable cutoff is ultimately needed, it should be confirmed by a
methodological reviewer who has not seen these outputs and who works from the algorithm
definition, simulation design and anonymized diagnostic fields only.

### Consequence for a capture judged unreliable

Retain its cells. Do not delete the capture. Capture-level exclusion requires independent
capture-level QC failure evidence, not Scrublet's inability to resolve doublets there.

## 2026-09-23 — Gate 2 scoring-only sensitivity design and Stage 31 redesign

**Gate 2.** The sensitivity-design packet received a methods review recorded in
`codex_output/20260923_111811.md`. Use a direct residual expected-rate stress grid
(`0.005, 0.01, 0.02, 0.04, 0.08, 0.12, 0.16, 0.20`) rather than a base-rate × removal-efficiency
grid. These values test Scrublet calibration stability and are not cohort prevalence estimates.
Do not transform vendor recovered-stage rates into post-Gate-1 residual rates, and do not use
donor-proportion or removal-efficiency assumptions as runtime inputs. Qiu and Chongqing may be
standard scoring candidates; Cillo, Liu and COMBAT may be residual diagnostics; Ren remains
residual-scoring-only; Stephenson and Wilk remain ineligible without additional evidence.

This ruling permits preparation of a runtime-config/design diff for review only. It does not
authorize an HPCC run or deletion of any cell. Material call-set variation across the grid is a
method-instability result, not a reason to select a favourable grid point. Outcome information
must not enter configuration or reliability decisions.

**Stage 31.** The frozen dimension bundle establishes computational feasibility, not scientific
validity. The current pipeline caps the analytic data at 60 cells per patient before running
Scrublet independently by `Sample_ID`; eligible groups contain only 30–60 cells. No fixed
`n_prin_comps` is approved for evidence generation under that order of operations. Scrublet must
instead run on the complete pre-subsample per-capture/sample raw-count input, with calls/status
mapped to the 60-cell analytic subset afterward. If those complete inputs or mappings cannot be
recovered, the Stage 31 doublet-rate outcome analysis remains outside the evidence chain.

A 2026-09-23 read-only HPCC check found that the September multi-cell-type DecontX work can be
reused for provenance and mapping but not as the scoring matrix. The assembled object contains
11,939 deliberately sampled cells from the canonical platelet subset plus Ren and Stephenson
multi-lineage backgrounds. The raw Ren and Stephenson sources remain available; their loaders
and mappings should seed a new inventory. The authoritative inventory scope is the Stage 31
COVID pre-subsample set (302 Sample_IDs, nine `Data_NO` groups), not the eight Gate 2 cohorts;
301 survived the old patient-balanced subsample. A source is not counted as covered until a
Sample_ID/capture-level join proves it.

Liu's 55 source partitions are promising but not yet runnable. RNA features are consistent
across all 55 after excluding ADT; the fixed 30,000-barcode selection mechanism and the meaning
of `10xlaneNN` remain unresolved. Existing rebuild scripts already encode a lane+barcode to
demuxlet-donor bridge and must be reused, but the bridge to canonical Seurat-suffixed cell IDs
and donor×timepoint Sample_IDs still requires evidence. Exact intersection of the unnormalised ID
strings is not evidence that no mapping exists.

The methods review was produced by Codex at the project owner's request. It is a documented
methodological ruling for this workflow, not a human external-review handover and not a claim of
pre-registration; `BUILD_PROVENANCE.md` remains unfilled until the PDF is actually circulated.

## 2026-09-24 — Written Gate 2 and Stage 31 methodology adjudication

At the owner's explicit request, Codex completed the methodological review rather than another
delivery audit. The full, source-linked ruling is
`reports/manuscript/revision/methodology_review_2026-09-24/REVIEW.md`. It records Q1–Q5,
eight cohort decisions, six Stage 31 redesign decisions and the five deferred PCA answers.
This is an AI review with prior involvement in package verification and exposure to existing
opinions; it is not blind review, external human sign-off or pre-registration. The missing
September 23 original was not reconstructed. This new record governs the qualifications below;
the historical entries above remain unchanged.

**Gate 2: scoring/sensitivity only, conditionally approved.** Retain the fixed algorithmic
stress grid `0.005, 0.01, 0.02, 0.04, 0.08, 0.12, 0.16, 0.20`. Its endpoints are not bounds on
cohort prevalence. Qiu/Chongqing may be standard scoring candidates subject to proven count
sources and physical capture units. Cillo/Liu/COMBAT may support residual-input diagnostics;
their outputs are not calibrated same-donor doublet rates. Ren remains residual-scoring-only;
Stephenson and Wilk remain ineligible on the current evidence.

Post-demux/HTO filtering changes the mix of doublet classes. Unrestricted synthetic pairing of
the survivors need not represent that residual mix; changing the prior alone cannot repair
it. `e` remains an explanatory parameter in `[0,1]`, not an additional runtime axis. The general
two-cell selection model uses singlet/within-source/cross-source survival probabilities:
`D = r * (s*a_W + (1-s)*a_X)` and `r_res = D / ((1-r)*a_S + D)`. The old `e=1` expression is
only a special case. Post-demux donor shares have no identifiable universal bias direction.

Hold simulations, preprocessing, PCA and the neighbour graph fixed within each capture/seed
across priors. The Bayesian score is monotone in the simulated-neighbour fraction, so rank
stability and observed-versus-simulated AUC stability are implementation checks, not independent
biological validation. Threshold/call changes are legitimate sensitivity findings; they do not
invalidate a correctly conducted stress experiment. Unknown/failed/unmapped states remain NA.

**Stage 31: workflow direction conditionally approved; no fixed PCA value selected.** Score
the complete cell-called/basic-QC RNA population of a proven physical capture before the
60-cell patient subset or lineage selection. Map results back by unique barcode/capture identity.
"Complete" does not mean unfiltered empty droplets; Sample_ID, donor and sequencing lane are
not automatically capture identities. The source inventory covers the 302 pre-subsample
Sample_IDs; the old 301/161 audit cannot select dimensions for the new inputs. Require complete
status accounting, zero silent exclusions, and a new dimension audit on the admissible inputs.
Unknown source/mapping/parent-state support blocks the affected scientific analysis; it does
not authorize deleting the affected cells or patients.

The review corrects the prior runbook's vendor-to-residual instructions and its proposed reuse
of a 10x model for Wilk. Seq-Well is not intrinsically incompatible with expression-based
diagnostics; Wilk remains excluded because the current source/filter-stage/unit evidence is
insufficient. A platform-specific prevalence claim would need additional calibration.

**Execution boundary.** This review authorizes documenting the design and preparing a scoped
implementation difference after its input conditions are met. It does not create runtime
configuration, fill the stress runner's eleven unadjudicated production parameters, authorize
HPCC execution, authorize outcome analysis, or authorize cell deletion. No external handover
entry is filled on the basis of this in-session review.

### 2026-09-24 — All eleven Gate 2 parameters adjudicated; execution gates retained

At the owner's request, the reviewer adjudicated all eleven parameters together in
`reports/manuscript/revision/methodology_review_2026-09-24/GATE2_PARAMETER_RULING.md`
(with a CSV and tagged-source verification record). The parameter packet and cohort eligibility
draft were received during this work and reviewed; Stage 31/Gate 2 scope conflation and the
packet's prior-SD description were corrected against Scanpy 1.11.3 source.

Approved for the scoped Gate 2 scoring-only pilot: sim_doublet_ratio=2.0,
synthetic_doublet_umi_subsampling=1.0, log_transform=False, n_prin_comps=30,
n_neighbors=None resolved as round(0.5*sqrt(post-filter observed N)), stdev_doublet_rate=0.0,
knn_dist_metric='euclidean', normalize_variance=True, mean_center=True,
use_approx_neighbors=False, get_doublet_neighbor_parents=False. Each approved range is a
singleton or that frozen k rule. Rationale, scope and invalidation are in the ruling, including
strict dimension feasibility, no fallback changes and no prevalence/optimality claim.

This closes the eleven methodological decisions; it does not activate the production slot.
The runner now resolves k before passing and recording it, fixing a None-versus-integer echo
mismatch. Config audit remains NOT READY: capture eligibility and direct cross-prior PCA/graph
identity evidence are still outstanding. No cohort data or outcome was scored.

The owner's latest execution instruction is: **11/11 approved + at least one capture meeting
C1–C3 + config audit passed → reviewer records a scoped small-pilot authorization → run**.
That conditional delegation does not authorize a pilot today, whole-cohort expansion or deletion.
Stage 31 is separately held until raw-source, physical-capture and raw-barcode identity close;
the reviewer may then authorize a new **dimension-only audit**, not detection or outcome work.
Gate 2's PC=30 does not settle Stage 31's deferred parameter. Earlier statements that only the
grid is approved or all eleven parameters remain undecided are superseded by this entry.

## 2026-09-24 — GSM5102900_HC1 passes C1–C3 and the score-free preflight

**Capture `GSM5102900_HC1` (Qiu, GSE167363) is the first Gate 2 unit to close C1–C3, on the
complete 7,666-cell deposit input.** One capture only; nothing here generalises to another.

- **C1 PASS.** The deposit's MatrixMarket banner declares `integer`, 0 non-integer values across
  9,903,096 streamed entries, `feature_types = {Gene Expression}` only, and GEO records
  `cellranger count` producing a feature-barcode count matrix. Un-normalised RNA counts with the
  source stated, not assumed.
- **C2 PASS.** 7,666 barcodes, all unique, all 16 nt, GEM-well suffix `-1` throughout. 10x
  documents that suffix as the GEM-well identifier, so a uniform `-1` is positive evidence for one
  GEM-well namespace. The per-GSM file count alone is not the proof.
- **C3 PASS on the deposit, FAIL on the canonical container.** The deposit is the complete
  CellRanger-called population at the declared stage. The container keeps **65 of 7,666 cells
  (0.85%)**, lineage-selected; the 65 join back cleanly (65/65) but completeness fails.
  **A scoring experiment on this capture must take the deposit input.**

**Score-free preflight, job 29073092, COMPLETED 0:0, 2:18, MaxRSS 2.13 GB** (commit `9aab7be`):
`N_eff = 7,666`, `G_HVG = 1,576`, so **PC=30 is strictly feasible** (30 < 1,576). Observed
manifold (7666, 30) and simulated (15332, 30) are finite with full rank; `n_neighbors = 44`,
`k_adj = 132`, seed `10,137,777` from `capture_seed`. No score, call or rate was produced and no
outcome column was read.

**The kNN graph remains unverified.** It is built inside `calculate_doublet_scores`, so the
score-free preflight could not construct it. `graph_preflight_verified = false` is carried in the
scoring config with its reason, and this limitation must be restated wherever the config is cited.

**No scoring permit exists.** `PRODUCTION_PARAMETER_RULING` stays `{}`. The scoring config is
non-runnable by construction and its wrapper refuses without a scoped pilot-permit artifact whose
`config_sha256` and `input_sha256` match the audited values. Evidence for the above is archived at
`reports/manuscript/revision/gate2_capture_audit_2026-09-24/run_evidence_29073092/`.


## 2026-09-25 — Gate 2 single-capture candidate passes final review; commit/push approved

The owner limited this review to frozen bytes/evidence, permit satisfiability/enforcement, and
pilot fail-closed behavior. The decision and exact approved digests are in
`reports/manuscript/revision/gate2_capture_audit_2026-09-24/FINAL_REVIEW_2026-09-25.md`.
Focused 73 passed, real Scanpy integration 12 passed, complete suite 726 passed: all final runs
exit 0. M4–M12 were each injected and killed; all 529 source/test file hashes match before,
after, and following every byte-identical restoration. An initial full-suite import-path failure
is preserved separately; fixing the invocation PYTHONPATH required no frozen-file edits.

One substantive defect was fixed within scope: identity-mismatched scores can no longer be
published under the input cell order, and the entrypoint stops before complete on identity drift.
No outcome access, cell deletion, production-ruling activation, or cohort scoring was performed.
The candidate is approved for commit/push. Next is a tracked scoped permit naming this reviewed
implementation commit, then deployment to a clean worktree and dry-run. The reviewer then checks
only the concrete permit/dry-run to authorize the unique one-capture scoring pilot. Graph topology
acceptance is enforced during that pilot; distances are recorded without an invariance gate.
No further architecture expansion or general mutation campaign is required by this decision.

## 2026-10-04 — Record batch-effect checking as required follow-up

The owner requested: “关于这个课题记一下关于 batch effect 相关.然后要做 batch effect 检查.”
Batch-effect assessment is now an explicit pending task in `docs/current-status.md` and
`docs/execution-plan.md`, covering batch metadata/outcome confounding, expression and five-FM
embedding diagnostics, and interpretation alongside patient-grouped CV/LOCO. The checklist is
an implementation plan; no correction method or metric threshold has been selected.

Local inspection found that `scripts/01_embedding_integrity/` reports numerical integrity
(NaN/Inf, dimensions, duplicates and related checks), not completion of the rebuilt benchmark's
batch-effect assessment. No `.h5ad`, `.npy`, `.npz` or `.parquet` inputs were found by the local
tracked/unignored file search. The recorded QC rebuild and full-scale embedding dependencies
remain pending; HPCC state was not checked in this documentation update. No quantitative
batch-effect analysis was run, and no conclusion about its presence or absence is recorded.

### 2026-10-08/09 — Batch-effect metadata audit v1/v2; scIB batch-ASW fix; HPCC checkout found stale

`scripts/33_batch_effect_audit/` built and run on HPCC (a dedicated scratch copy of inputs, not
the `302006_updategit` checkout found there, which was stale at `35131cf` — dozens of commits
behind this repo's `644bcf4` — and carried uncommitted Gate-2-related files; that checkout was
left untouched, not reset or cleaned). v1 (2026-10-08) had real defects a same-day review found
and v2 (2026-10-09) fixed, confirmed independently before fixing: duplicated Cillo donor rows and
a literal `"nan"` donor ID from a NaN-handling bug; COMBAT's donor-coverage row read as "0/124
recovered" when the real situation is that cohort's file has no donor identity to join on at all
(now `coverage_status="not_assessable"`); v1's "8 unsplittable blocks" language overclaimed —
these are CV dependency components with outcome overlap, not a demonstrated batch-biology
association, and v1 silently omitted the 4 homogeneous-outcome multi-donor components (12 total,
not 8). Both versions are kept (`scripts/33_batch_effect_audit/results/MANIFEST.md`); use v2.

**New, unresolved finding:** all 15 recorded Stephenson pools (`stephenson_pool_composition.csv`)
contradict `canonical_donor_labels.csv`, which gives every pooled donor its own singleton
`cv_group`. Several donors recur across multiple pools. What a Stephenson "pool" physically
denotes (a real multiplexed capture vs. a sequencing/library pool with no GEM-well implication)
is not established here and needs a check against primary E-MTAB-10026 metadata before it
changes any CV grouping.

Separately, `scripts/24_option_c_canonical_input/code/embedding_quality_option_c.py`'s "Batch
ASW" computed one global, unstratified `silhouette_score(X, cohort)` — not scIB's definition,
which conditions on a biological label before averaging, and so could read real cohort-level
biology as a batch effect. Fixed (`scib_batch_asw()`), hardened against the sklearn
`n_labels == n_samples` boundary and non-finite input, and covered by a committed, HPCC-run
regression suite (`tests/test_batch_asw.py`) that a same-day review confirmed does detect the
specific `np.abs()`-omission mutation it is meant to catch. No FM embeddings or expression
diagnostics on real data were run; that still needs the Gate-2 input and 5-FM-embedding
dependencies this file already tracked as pending.

### 2026-10-09 — Stephenson/Cambridge pre-capture pooling confirmed; v3 metadata fixes

A second review (`codex_output/20261009_101831.md`) of the 2026-10-08/09 batch-effect
metadata audit found and this session confirmed independently: (1) Stephenson rows used a
semicolon-joined, multi-pool `capture_id` with a patient-aggregate cell count mislabeled as a
capture count (CV0068: 4,951 patient total vs 2,701/2,250 real per-sample_id counts in
`stephenson_sample_severity.csv`); (2) the alias-normalization dict comprehensions could
silently drop one of several IDs sharing a normalized key; (3) donor-uniqueness was checked
globally rather than per-cohort; (4) `classify_event` used a tolerance-based comparison instead
of exact equality. All fixed in v3 (`scripts/33_batch_effect_audit/`); see
`scripts/33_batch_effect_audit/results/MANIFEST.md`.

**The Stephenson pooling question is now resolved for the Cambridge subset.** Primary source
(Stephenson et al. 2021, *Nature Medicine* 27:904-916, Methods, fetched from PMC8121667):
Cambridge "pools of four samples were generated by combining 0.5 million live cells per
individual" BEFORE antibody staining and 10x Chromium loading, with downstream donor splitting
by genotype (not TotalSeq-C hashtag — that panel is a surface-protein assay at this site, not a
donor barcode). Cross-checked against this project's own data: every donor named in any of the
15 recorded `stephenson_pool_composition.csv` pools has `Site=="Cambridge"` in
`stephenson_patient_outcomes.csv`. Conclusion: these pools are real physical multi-donor
captures, not merely sequencing-library groupings, for these 47 donors specifically. Sanger and
Newcastle donors have no pool evidence in local metadata and remain unresolved, not asserted
independent.

`scripts/33_batch_effect_audit/code/resolve_stephenson_grouping.py` computes the proposed
connected components (same construction already used for Liu/Cillo/COMBAT's `cv_group`): the 47
Cambridge pool members collapse into 3 components (21/22/4 donors), the two large ones mixing
fatal and survived donors. **This is a proposal, not applied** — `canonical_donor_labels.csv`
is unchanged; promoting the proposal (and its consequence, that Stephenson would gain 2 more
mixed-known-outcome dependency components, 8→10 project-wide) is a separate decision.

Also noted, audited but not acted on: `build_stephenson_decontx_input.py` sets
`channel_id = sample_id` and its own docstring calls this the GEM-well/library unit; for the 47
Cambridge pool members, `sample_id` is actually one donor's share of a shared physical pool, so
the already-completed DecontX run for these donors was fit at a finer granularity than the real
ambient-RNA-sharing unit. This is a methodological finding about a completed run, not a new QC
rule; no historical output was touched and no rerun was performed or authorized by this entry.

### 2026-10-09 (continued) — Stage-33 diagnostic CLI implemented, validated on synthetic fixture

A third review (`codex_output/20261009_113232.md`) asked that the diagnostic CLI, a real Harmony
correction path, evaluation-boundary tests, and remaining v3 metadata gaps not be deferred again.
All implemented this session:

- **v4 metadata**: a recurrence of the NaN-to-literal-"nan" bug (fixed once for Cillo in v2)
  found in Stephenson's own `sample_ids`/`pools` parsing in v3, confirmed directly (AP1 had
  `pool_id="nan"`) and fixed with a shared `_split_semicolon_field()` helper used by every
  adapter. Added an explicit `pooling_status` enum (`confirmed_pooled`/`confirmed_single_donor`/
  `unknown`) so a missing pool record is never read as a confirmed-independent capture. Stephenson
  (sample_id, pool) pairs are now cross-validated against the independent
  `stephenson_pool_composition.csv` donor lists, not merely checked for equal list length. The
  alias crosswalk now covers every canonical/capture-file ID (not just the previously-unmatched
  subset) and gives an unmatched ID its own explicit row rather than silently dropping it.
- **Shared batch-ASW module**: `scib_batch_asw` moved to `scripts/00_shared/code/batch_asw.py`
  (hardened further: explicit 2-D/1-D shape checks, missing-label rejection) so
  `embedding_quality_option_c.py` and the new diagnostic CLI share one implementation. The second
  legacy caller (`scripts/12_umap_visualization/.../batch_effect_analysis.py::compute_batch_asw`)
  now raises `NotImplementedError` pointing to the shared module if called, rather than being
  merely documented as retired while still able to silently run.
- **Stage-33 diagnostic CLI** (`run_batch_diagnostics.py`) and a synthetic fixture generator
  (`make_synthetic_batch_fixture.py`) built and run end-to-end on HPCC. First fixture version had
  a real bug (Python's per-process-randomized `hash()` used for outcome assignment, which also
  concentrated every "fatal" patient into one cohort) — found via the diagnostic run's own
  manifest (`batch_asw` skipping the `fatal` stratum with `n_batches=1`), fixed with a properly
  seeded RNG, and re-verified deterministic (`build(seed)` called twice gives byte-identical
  output). On the corrected fixture: the injected batch shift was detected and measurably reduced
  by an opt-in real Harmony correction (harmonypy 0.0.10) in all 8 lineage×embedding×scope
  combinations (neighbor-batch-mixing ~0.09–0.18 pre-correction vs. ~0.75 expected under no batch
  effect, rising to ~0.51–0.69 post-correction). Harmony is explicitly labeled transductive/
  exploratory, never an out-of-sample transform — enforced by a harmonypy API bug fix (this
  version's `sigma` auto-broadcast silently breaks below `nclust==1`, i.e. N<45 — fixed by passing
  `nclust` explicitly) and proven by a dedicated no-persisted-state test.
- **Evaluation-boundary tests** (`tests/test_evaluation_boundaries.py`) added against the
  existing, already-used `repeated_cv.py` LOCO/grouped-CV harness (not a new implementation):
  group integrity under `StratifiedGroupKFold`/`LeaveOneGroupOut`, train-only fitting proven with
  a deliberately leaky negative control confirming the check can actually catch leakage, and the
  Harmony no-persisted-fit-state contract.
- **Provisional job specification** prepared and NOT submitted
  (`scripts/33_batch_effect_audit/results/PROVISIONAL_JOB_SPEC_v1.md`): a bounded,
  patient-balanced extraction from the real pre-Gate-2 assembled input
  (`05_RESULTS/rebuild_upfront_qc/platelet_monocyte_benchmark_input_v2/`, 22 GB/50 GB
  `counts.mtx` files, streamed via the already-production `read_mtx_columns()`), real `--dry-run`
  verified at 665 donors → 13,232 cells, followed by the diagnostic CLI in expression/PCA mode
  (no FM embeddings exist yet). Exact command, resource ceilings, and acceptance/failure criteria
  are in that file. No SLURM job was submitted by this entry.

41 tests pass together across `scripts/33_batch_effect_audit/` (`tests/test_batch_asw.py`,
`test_batch_metadata_inventory.py`, `test_batch_confounding.py`, `test_adapt_stephenson.py`,
`test_evaluation_boundaries.py`), run on HPCC in `302006_env_clean`.

### 2026-10-09 (continued, 4th review round) — reproducible sampling, enforced run status, real-harness tests

A 4th review (`codex_output/20261009_133450.md`) of `9c7767f` found and this session fixed, each
confirmed against the actual code/output before fixing:

- Both `run_batch_diagnostics.py` and `prepare_provisional_diagnostic_input.py` still used
  Python's built-in `hash((seed, id))` for per-patient sub-seeding -- `hash()` on a
  string-containing tuple is per-process-randomized (`PYTHONHASHSEED`), so a fixed `--seed` did
  not actually reproduce the same selection across separate process invocations. New
  `scripts/00_shared/code/stable_seed.py::derive_seed()` fixes this, proven stable across
  `PYTHONHASHSEED` values by a committed test (and the original bug proven real by a
  negative-control test against the same environment).
- The real-data preparer capped cells per bare `donor_id` **across both lineages combined** --
  a donor with cells in both Platelet and Monocyte_Myeloid could lose most of the rarer lineage
  to the more abundant one in the random draw. Fixed to cap within each
  `(cohort, donor_id, lineage)` combination separately. **Real, material impact, confirmed by
  rerunning the actual `--dry-run` against the real HPCC input**: estimated usable cells roughly
  doubled, 13,232 → 24,841, with Platelet recovering most dramatically (e.g. COMBAT's Platelet
  count: 43 → 1,775; Ren's: 481 → 3,444).
- A Harmony correction or an explicitly-requested embedding that failed or was unavailable could
  leave the enclosing embedding, and the whole run, `status="complete"` with exit 0. Required
  work failing now sets `status="failed"` and a non-zero exit; optional misses (e.g. an
  unsupported, not-requested embedding) do not fail the run. `manifest.json`'s `code_revision`
  could read `"UNKNOWN (...)"` while `status` still said `"complete"` -- `--code-revision` can now
  be passed explicitly (required for a loose-file scratch deployment, which is not a git
  checkout), and the run now fails rather than publish unresolved provenance silently.
- Saved representations (`.npy`) had no cell-ID binding; UMAP coordinate CSVs carried
  `patient_id` (which repeats across cells) but no `cell_id`. Both now save/verify cell-ID-bound
  `.npz`/CSV outputs, reloaded immediately to confirm shape/identity/values survive the round
  trip. Neither CLI had a no-clobber check; both now refuse to run into an existing output path.
- `tests/test_evaluation_boundaries.py`'s group-integrity and leakage tests were rewritten to
  instrument the REAL `repeated_cv.py` harness (via pandas `.iloc` row-index identity carried
  through `run_repeated_kfold`/`run_loco`), not a separately re-constructed splitter or an
  untouched-column comparison -- the prior leakage test's invariant held trivially whether the
  tested implementation was correct or leaky; the new one perturbs ONLY a real fold's actual
  held-out rows and proves a deliberately leaky implementation's captured fitted state changes
  while a correct one's does not.
- New `tests/test_prepare_provisional_diagnostic_input.py` exercises the real `read_mtx_columns`
  streaming path against a small, real, on-disk MatrixMarket fixture (not `--dry-run`, not a
  mock) -- write, stream, h5ad write/reload, and value-correctness all verified, as the review
  asked be done before presenting the real 72 GB scan for approval.
- `is_multi_donor_capture` is now always derived centrally from `pooling_status` in
  `build_batch_metadata_inventory.py`'s `_row()` (True/False/`pd.NA`), not passed ad hoc per
  adapter -- "unknown" pooling can no longer read as `False` ("confirmed independent") in any
  adapter. The synthetic fixture's outcome is now assigned once per `(cohort, patient)` and
  reused across lineages/draws (previously re-derived per lineage, so the same patient could get
  a different outcome in each lineage). `MANIFEST.md`'s Stephenson-proposal arithmetic (12→14)
  was corrected to 12→15. The Stephenson sample-to-pool pairing check's known limitation (cannot
  detect a swapped pool order for a donor who is a genuine member of two pools, since the check
  is membership-only) is now documented and characterized by a test, not silently presented as
  fully resolved.

Re-running the real `--dry-run` also surfaced an infrastructure finding: the command was
SIGKILL'd (exit 137) twice on the HPCC login node while reading `decontx/cells.csv`, even though
`free -h` showed 441 GB free on that node at the time -- consistent with this project's own
documented operational constraint that login-node process limits are not simply a function of
total node memory. Fixed by a memory-pattern change (reading both files' `cell_id` column alone
before the heavier full load), not a logic change; the dry-run then succeeded in ~11 seconds.

`scripts/33_batch_effect_audit/results/PROVISIONAL_JOB_SPEC_v2.md` replaces v1 (kept for the
audit trail) with the regenerated real numbers, strict `status: "complete"` acceptance (not
`complete_with_partial_results`), and an explicit two-allocation (prepare, then diagnose) design,
each needing its own per-job authorization. No SLURM job was submitted. A fresh synthetic
validation run referencing the real commit (`synthetic_fixture_v3`) completed with strict
`status: "complete"` and `code_revision_verified: true`, reproducing the same qualitative
before/after-correction pattern as the earlier runs. 52 tests pass together on HPCC
(`302006_env_clean`), confirmed from a run writing to shared `/bigdata` storage after two
identical runs' output was invisible from a follow-up check due to the documented login-node
`/tmp`-is-not-shared round-robin behavior.

## 2026-10-09/10 — execution-readiness repair, before/after plots, full Stephenson table

Closed the two items deferred at the end of the previous round (rendered before/after UMAP plots;
the full 130-donor Stephenson readiness table, not just the 47 Cambridge pool members) — see
`scripts/33_batch_effect_audit/results/FINDINGS_REPORT_v1.md` §9. Then a fresh review
(`codex_output/20261009_230728.md`) found the provisional-job packet itself was not yet safe to
execute as documented: `PROVISIONAL_JOB_SPEC_v2.md`'s own two commands' `--population-stage`
values did not match each other, which the runner only warned on; `code_revision_verified` was
`code_revision is not None`, so any string "verified"; both CLIs' output guards were
check-then-act, not atomic (the real underlying cause of the earlier process race, independent of
the orphaned-process explanation given at the time); unknown-outcome cells were withheld from the
entire technical-integration population, not just bio-labeled scoring; a failed corrected-UMAP
did not propagate to run status. All five fixed with regression tests, not just documented --
`scripts/33_batch_effect_audit/results/FINDINGS_REPORT_v1.md` §10 has the full detail. New shared
modules `scripts/00_shared/code/{population_stage,run_ownership}.py`.
`PROVISIONAL_JOB_SPEC_v3.md` replaces v2 (v1/v2 kept for the audit trail), with real sbatch
scripts (`scripts/33_batch_effect_audit/code/slurm/step{1,2}_*.sh`) carrying exact interpreter
paths -- neither submitted. `tests/test_execution_readiness.py` (new, 15/17 passing + 2 honest
environment-dependent skips) and `tests/test_preparer_to_diagnostic_integration.py` (new, a real
end-to-end preparer-to-runner run against the exact job-spec v3 argument contract, passing). The
six pre-existing test files this round did not touch were not re-run together this round -- two
attempts at a combined `pytest tests/` run were killed by genuinely heavy HPCC login-node load
that evening before completing; marked `UNVERIFIED (this round)` rather than assumed.

## 2026-10-12 — Owner changes the objective to scFM pathway discovery; framework first

**Authority:** explicit owner instructions in the current conversation, confirmed in successive
messages. This is a scientific scope decision, not an inference from an AI review.

**Decision:** build a reusable Docker-based framework that uses pretrained scFMs to identify cell
states and disease-associated gene programs, connects gene evidence to GSEA and GO/KEGG annotations,
and produces traceable biological insights. The project is no longer a benchmark/mortality-ranking
study. Baselines, reproducibility checks and batch diagnostics remain validation tools.

**Sequence and scope:** build the framework first. The human application is restricted to
monocytes. Reserve species, study/source and patient identity in its interfaces. Proceed one reviewable increment
at a time. The next increment is framework input/model contracts, validate/plan and provenance.

**Supersedes:** the primary benchmark objective, human platelet application in new work, and the
fixed five-model set as a compulsory framework scope. The old separation of mouse work into 302035
does not prevent mouse support in the framework; it remains relevant to data provenance and does
not authorize moving or merging original projects. Preserve historical code and results.

**Not decided:** specific biological contrasts, final checkpoint set, production batch correction,
or final data-retention policy. No new SLURM job, push, deletion or data upload is authorized here.
Existing full-capture QC contracts still apply where relevant before extracting the study lineage.

**Implementation status:** scope/status/plan documentation updated; framework implementation and
real-data execution are separate deliverables, not claimed complete by this decision.

**Handoffs:** parent-project `codex_output/20261012_133205.md`, `20261012_133528.md`,
and `20261012_135325.md`. The latest scope takes precedence where the earlier handoffs differ.

## 2026-10-12 — Owner narrows the framework to one primary scFM

**Authority:** the owner explicitly requests one useful scFM for the framework, with a reasoned
comparison in the paper explaining its suitability relative to alternatives.

**Decision:** a single primary backbone is sufficient. The framework need not support all scFMs;
multi-model consensus, a second production adapter, and a broad benchmark are not required
deliverables. Keep model identity and species contracts explicit without building an unnecessary
multi-model orchestration platform. The exact model/checkpoint remains unselected.

**Evidence boundary:** assess suitability for cell-state/gene-program/pathway discovery, human
monocytes, plus practical reproducibility and deployment. Small
targeted comparisons and an expression baseline may substantiate the paper; they do not reinstate
the old benchmark objective. Do not assert biological superiority from usability alone, or attribute
framework-level gains solely to the pretrained model. Prefer an independent validation stage after
model selection and retain negative findings.

**Supersedes:** earlier handoff requirements to add a second model adapter or produce cross-model
consensus as part of the MVP. Framework-first and human-monocytes-only remain unchanged. No new execution authorization is implied.

**Handoff:** parent-project `codex_output/20261012_140650.md`.

## 2026-10-12 — Limited scFM shortlist prepared (recommendation, not a decision)

scPRINT-2 was proposed as a candidate for minimal validation, with UCE as an
alternative. This was a recommendation, not a production model selection.
Species coverage does not establish downstream efficacy. The detailed source
record is preserved in the related repository
(https://github.com/xqiu625/302033_git/blob/master/archive/scfm_benchmark_context_20261113/docs/decisions.md).

## 2026-10-14 — Owner selects the within-cohort COVID severity question

**Authority:** explicit owner response in the active conversation:
“同一 COVID 队列内：重症 vs 轻/中度，关注细胞状态与程序变化（建议）”.

**Decision:** the human monocyte pilot will address severe versus mild/moderate COVID-19
within one cohort, focusing on cell-state and gene-program differences. This selects the
scientific question; it does not reinstate the former mortality benchmark or authorize pooling
COVID-19 and sepsis cohorts. Patient/donor remains the biological replicate.

**Still to resolve before real analysis:** which eligible cohort, the exact severity definitions
(including how moderate/critical categories map), sampling-time window, monocyte subtype scope,
repeated-visit handling, retained/QC input version, positive-control signature and its expected
direction, donor-level validation split, and the main model/checkpoint. These are not implicitly
approved by the question selection. A proposed CD14/classical-monocyte restriction is still a
recommendation, not part of the owner's literal choice.

**Execution boundary:** metadata feasibility checks, literature verification and already-scoped
synthetic framework repairs can continue. No new SLURM job, installation, model-weight download,
push, QC-retention change or real-data fitting is authorized by this decision.

## 2026-10-15 — Owner requires direct human-in-the-loop decision dialogue

**Authority:** owner supplied the “Agentic AI in Bioinformatics” diagram and stated:
“我希望我们遵守这个 然后 任何需要我裁决的都直接问我,我们沟通,或者我去读 paper 我去找数据库 裁决.”

**Decision:** use the goal -> plan -> authorized execution -> evaluation -> scientific-output
loop with Cindy as the scientific decision maker. Codex/Claude must surface decisions directly
in conversation, supported by evidence, uncertainty, recommendations, alternatives and impact.
Do not leave actionable scientific choices only in files marked `pending_owner`. Ask focused
questions as they become relevant, and give precise paper/database reading tasks when human
judgment or access is useful. Evidence gathering and independent authorized implementation continue.

**Scope:** this establishes how decisions are made, not the answers to outstanding scientific
questions. It does not choose a cohort, merge Critical with Severe, select a subtype, freeze a
positive control, approve a model/checkpoint, or authorize any new job. It does not add the
diagram's named software packages or wet-lab automation to the current MVP. Existing decisions
remain valid; receipt of a paper, an AI recommendation or silence is not approval.

**Next dialogue:** directly ask whether the primary severe group includes Critical, preserving
the already-confirmed Mild+Moderate reference group. At the time of this entry that question
has been presented and remains unanswered; do not implement a new grouping on this basis.

## 2026-10-15 — Owner merges Severe + Critical and adds a future ICU disease comparison

**Authority:** Cindy explicitly selected “主比较合并 Severe + Critical，对照 Mild + Moderate”
and explained “凡是 ICU 的取样 都是 severe 或者 critical.因为后续我还想放 sepsis 数据.
sepsis icu 和 covid-19 ICU 比较.”

**Confirmed current contrast:** within the selected COVID cohort, comparison group =
`{Severe, Critical}`, reference group = `{Mild, Moderate}`. Direction = pooled Severe/Critical
minus pooled Mild/Moderate. Preserve the original severity labels alongside the derived group.
This supersedes the earlier pending Critical decision and the recommendation to use Severe alone;
do not ask the owner to choose this grouping again. This is a decision, not an executed analysis.

**Confirmed future direction:** prepare a later comparison of monocyte states/programs between
sepsis ICU samples and COVID-19 ICU samples. This does not replace the current within-COVID
contrast, authorize immediate cross-cohort pooling, choose the sepsis dataset/definition, or
cancel separately maintained research. Do not reframe this as a benchmark.

**Clinical mapping boundary:** the owner's statement is about ICU samples belonging to the
severe/critical group, not proof of the inverse. Stephenson Methods defines Severe as ICU and/or
noninvasive ventilation; therefore severity grouping and actual ICU location at sampling must be
recorded separately. Keep source evidence and unknown values for sampling location. Never infer
ICU status solely from Severe/Critical or from a later/worst clinical status.

**One follow-up question is pending:** should the current COVID comparison additionally restrict
its Severe/Critical arm to confirmed ICU sampling, or keep the severity-based main comparison and
reserve strict ICU restriction for the later cross-disease analysis? This has been asked directly.
Until answered, do not invent an ICU eligibility filter. Evidence gathering and backend work can
continue independently of this question.

**Remaining design decisions:** final cohort, subtype scope, visit selection/time window, identity
admission/QC input, positive control, model/checkpoint and validation split remain unresolved unless
separately approved. Later sepsis work needs an explicit case definition (including COVID overlap),
sampling-time alignment, donor structure and disease-by-site/batch feasibility review. This entry
does not authorize new SLURM jobs, installation, downloads, real fitting, push or changes to QC.

## 2026-10-15 — Owner confirms ICU-at-sampling restriction for the COVID severe arm

**Authority:** Cindy explicitly chose “B：当前 COVID 重症组就只保留采样时在 ICU 的
Severe/Critical。是的这个比较好解释 还和临床合并.”

**Decision:** eligibility for the current COVID comparison arm is the conjunction of
Severe/Critical severity and confirmed ICU location at the same sample's collection time.
The reference arm remains Mild/Moderate. This resolves the preceding ICU-only question;
do not ask it again or substitute the previous recommendation to reserve ICU restriction for later.

**Implementation meaning:** keep original severity and care-setting fields plus source evidence.
Severe/Critical without ICU evidence cannot enter the primary comparison arm. Missing/ambiguous
location remains unresolved, not inferred as ICU or recoded as non-ICU. Future ICU admission,
worst clinical status and ventilation alone are not proof of location at the sampled time.
Do not move excluded severe/critical records into the Mild/Moderate reference arm. Unexpected
severity/location conflicts must be reported without silently overriding the source classification.

**Scope of interpretation:** the main contrast concerns ICU-sampled severe/critical COVID versus
mild/moderate COVID, not an isolated causal effect of disease severity or ICU care. This supports
the planned later ICU sepsis–COVID comparison but does not by itself make the cohorts comparable.

**Still open:** repeated-visit selection is a separate scientific decision. A donor initially
Moderate and later Critical in ICU cannot be assigned a visit by silently choosing Initial-only,
the most severe visit, or the earliest ICU visit. Prepare a donor/visit eligibility table and the
consequences of candidate rules before asking Cindy to decide. Other unapproved design choices
and per-job execution boundaries remain unchanged. No real analysis was executed by this entry.

## 2026-10-15 — Correction: the three entries above are real, not fabricated

**What happened:** while reconciling several parallel subagents' work this round, Claude found
the three entries above (plus a matching note in `CLAUDE.md`) already present in the working
tree, attributed to the user, with no record in its own conversation of anyone asking these
questions. Because none of its own subagents had been told to touch this file, and because one
subagent's completion happened to land at the same timestamp, Claude concluded a subagent had
hallucinated a fake dialogue with Cindy and stashed all four files (`git stash`, not discarded)
rather than commit them, and said so to the user.

**Correction:** this was a wrong inference about WHO wrote the content, not about whether it was
real. Per `codex_output/20261015_121928.md`: Codex (the separate reviewing agent, in its own
concurrent conversation with Cindy) wrote these four entries directly into the shared working
tree, and the decisions themselves are real — Cindy made them in that Codex conversation and
confirmed so directly, in her own words, in this session afterward. Independently verified
before restoring: the stash's commit hash (`dc69b75d...`) and diff stat (7/9/84/10 insertions)
matched exactly what `20261015_121928.md` cited, and both of Codex's own prior handoffs
(`20261015_120031.md`, `20261015_120317.md`) existed on disk with timestamps preceding the
stash. The four files were restored unchanged from the stash (`git checkout stash@{0} -- ...`);
nothing in them was rewritten. The stash entry itself is kept, unmodified, as the audit trail —
its label is now understood to describe Claude's mistaken initial read, not an actual
fabrication.

**Standing lesson:** a Codex session and a Claude Code session sharing this working tree is
itself a form of concurrent writer — same hazard class as two SSH sessions landing on different
HPCC login nodes (`CLAUDE.md`, Operational constraints). An unexplained change should be
preserved and the source checked (another session, or the user directly), never discarded or
asserted fabricated purely because "this session didn't do it." Per `20261015_121928.md` §4:
going forward, Codex writes only to `codex_output/`; Claude Code owns consolidating `CLAUDE.md`,
the status/plan docs, prereg files, and this decision log, so the two sessions stop double-
writing the same files without first saying so in a handoff.

## 2026-10-15 — Owner resolves the one cross-group donor (C-8821) found this round

**Authority:** Cindy's direct reply in this round's Codex conversation (via
`codex_output/20261015_142903.md` §1, recorded here by Claude per the division of labor in the
entry above): asked how to assign the one donor (`sdrf::C-8821`) whose visits land in BOTH the
ICU severe/critical arm and the mild/moderate reference arm under the current candidate roster,
Cindy replied (quoted verbatim): "如果只能进入一组, 就进入目前分组少的那组" (if it can only enter
one group, assign it to whichever group currently has fewer members).

**Decision:** `sdrf::C-8821` is assigned to the ICU severe/critical comparison arm (the
fewer-populated group under the roster as it stood when asked), using sample `newcastle004v2`
(D12, Critical, `ITU_intubated`). Its other sample, `MH9179821` (D0, Moderate, `Ward_O2`), is
excluded from this round's mild/moderate reference arm by this decision — the underlying sample
row, its identity, and its QC status are unchanged; it is simply not counted as an independent
reference-arm input for the current main comparison. Resulting donor-level counts, each donor now
counted in exactly one group: comparison=25 (unchanged), reference=52 (down from 53). Implemented,
reproducibly, in `scripts/34_scfm_discovery_framework/results/
stephenson_icu_severity_grouping_resolve_mixed_visits.py`, which pins the exact input file hashes
this decision was applied against and refuses to run again if the mixed-donor list ever differs
from exactly this one donor.

**Scope — explicitly NOT a general rule:** this is a one-time decision under the CURRENT candidate
roster and counts, not the Stephenson et al. authors' own analysis rule, and not derived from any
expression or pathway result. It must never be re-applied automatically if the roster, QC, or
candidate counts change. If a future round finds additional cross-group donors, a fresh, explicit
counting-basis/tie-break decision is required — never inferred from row order, and never decided
by re-running the same script blindly against a changed input (the script itself refuses that).

**Still open:** repeated-visit selection WITHIN an already-resolved group (e.g. a donor with two
Mild visits, or two ICU Severe visits) is a SEPARATE, still-unanswered question, already asked
directly of Cindy per a prior entry above — this decision does not answer it, and `sdrf::C-8821`
happens to have only one ICU-qualifying sample so is unaffected by it either way.

## 2026-10-16 — Stephenson monocytes: provisional raw counts accepted for a technical pilot (owner)

**Authority:** Cindy's direct reply in the Codex conversation, 2026-10-16: "允许".
The immediately preceding question asked whether the existing atlas raw counts may be used
for an input-to-state-discovery technical pilot, retaining all visits and attaching historical
QC flags, without admitting its outputs into disease/pathway scientific conclusions. The
question explicitly excluded SLURM authorization. Context and independent QC identity audit:
`codex_output/mono_20261016_130032.md`; decision handoff:
`codex_output/mono_20261016_130248.md`.

**Decision:** that limited provisional input-use scope is approved. Prepare the concrete pilot
against `02_EMBEDDINGS/option_c_monocyte/stephenson_monocyte_raw.h5ad` on HPCC. Retain all visits
and donor/sample/cell identities; attach the existing Gate 1 flags and historical DecontX
estimated scores as annotations. This approval does not apply retention filtering or replace
raw counts with historical corrected counts. Technical outputs remain provisional and are not
evidence for disease-associated programs or pathway conclusions.

**Effect:** the provisional-input-use question is resolved; do not ask it again without material
new evidence. The full-visit index remains available; earliest-visit selection remains deferred.
Existing ICU-at-sampling eligibility and C-8821 cross-sectional/longitudinal distinctions persist.

**Not decided or authorized:** final scientific QC admission, cell deletion, Gate 2/refit policy,
model/checkpoint selection, scientific DE/masking conclusions, any SLURM submission, dependency
installation, push, or automatic retry. Prepare an exact revision/input/model/output/resource/
command package before requesting the remaining model decision and job-specific execution
authorization. No real-data model job was run when recording this decision.

## 2026-10-16 — scPRINT-2 small-v2 selected for the first technical pilot only (owner)

**Authority:** Cindy's direct reply in the Codex conversation: "ok", answering the explicit
question whether scPRINT-2 small-v2 should be the first technical-pilot model, without locking
the final paper backbone or authorizing SLURM. Recorded during review of commit 603e213.

**Decision:** use scPRINT-2 `small-v2.ckpt` as the first technical-pilot model. Reuse the existing
environment/checkpoint evidence and repair the verified producer/consumer interface defects
before presenting the executable job package. Do not ask the first-pilot model question again.
Final backbone selection remains conditional on validation; UCE is not scientifically excluded.

**Scope:** the previous provisional-raw-count technical scope persists: retain all visits and
identities, QC flags/scores as annotations only, no disease/pathway conclusions. This reply does
not approve hidden model-preprocessing cell filtering, any SLURM submission, installation, push,
or automatic retry. No real-data inference was run when recording this decision. Execution
authorization still requires a corrected, concrete job-specific package.

### 2026-10-17 — Stephenson monocyte scPRINT-2 pilot v4: one SLURM submission authorized

**Source:** Cindy replied "授权" in the Codex conversation to the explicit request accompanying
`codex_output/mono_20261017_152624.md`: one v4 run, 1 GPU, 8 CPU, 128G, at most 2 hours,
prep → embed → validate only, without automatic retry. This records actual owner authorization,
not an authorization inferred from an AI review.

**Executor and scope:** Claude Code may perform the final drift/duplicate check and submit exactly
once using `sbatch --export=ALL,EXPECTED_DEPLOY_COMMIT=3de4bca47ad4ad07d501f557d7b9065540477655
/bigdata/godziklab/shared/Xinru/302006/302006_updategit_deploy_scprint2_pilot/scripts/foundation_models/scPRINT/run_scprint2_stephenson_monocyte_pilot.sh`.
Frozen deployment revision: `3de4bca47ad4ad07d501f557d7b9065540477655`; partition `short_gpu`;
1 GPU, 8 CPU, 128G, walltime ≤ 02:00:00, no requeue. Monitoring and output validation are included.
Use a fresh versioned run directory and preserve all 128,503 input cells and their identities;
Gate1 and historical DecontX scores remain annotations only.

**Conditions:** if this exact run has already been submitted, monitor its existing job ID instead
of submitting again. Stop on revision/hash/resource drift; do not silently substitute a deployment.
No automatic retry, replacement job, install, push, DE, masking, enrichment, or expansion to other
datasets. This does not select the final backbone or admit scientific findings. A walltime timeout
is resource-limited incompletion, not proof that the model is scientifically unsuitable.

**Status at record creation:** authorized; Codex has not submitted a job in this turn and has not
received a real submission receipt. Claude must record the actual job ID and terminal validation.
Execution handoff: `codex_output/mono_20261017_153849.md` (project-root relative).

### 2026-10-18 — Packet A single CPU validation authorized; overnight follow-through requested

**Source:** Cindy replied in the Codex conversation: "授权,可以多跑一些东西吗? overnight 跑",
directly answering the scoped Packet A authorization request in
`codex_output/mono_20261017_235549.md` (2 CPU, 32G, at most 30 minutes, original validator,
no GPU and no automatic retry).

**Authorized executor/scope:** Claude Code may submit Packet A once, after the final
drift/duplicate check, using
`/opt/linux/rocky/8.x/x86_64/pkgs/slurm/24.11.1/bin/sbatch /bigdata/godziklab/shared/Xinru/302006/codex_mono_packet_review_20261017_234357/scripts/foundation_models/scPRINT/run_scprint2_recovered_identity_validate_cpu.sh`.
Wrapper SHA-256: `11c9dce6bda8e17184ae6e4bdb98c1c5242ed6957d19624d7645d64f03753af3`.
Partition short; 2 CPU, 32G, walltime 00:30:00, no requeue. Includes monitoring and collection
of the original validator's report on the hash-pinned recovered_v2 artifact. The old validator
deployment remains pinned to clean `3de4bca47ad4ad07d501f557d7b9065540477655`.
If already submitted, monitor the existing job rather than submitting again. Stop on drift.

**Overnight intent:** continue useful implementation, bounded tests and preparation alongside
Packet A, without re-requesting approval for ordinary reversible steps. At this decision time,
Packet B still has unresolved execution/checker defects documented in the review above; this
message does not identify a corrected B revision/command or grant a retry, an additional GPU
submission, installation, push, changed scientific criteria or scientific-result admission.
Prepare concrete additional execution packets before seeking their job-specific authorization.

**Status at recording:** A authorized; Codex has not submitted it in this turn. A successful
identity_validation does not close producer_row_alignment. Overnight handoff:
`codex_output/mono_20261018_000122.md` (project-root relative).

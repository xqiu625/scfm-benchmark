# Model/checkpoint decision page v1 — short, per `review:20261014_204536` §5.A.5

This is a decision PAGE, not a new evaluation — it restates and organizes what
`MODEL_SELECTION_v1.md` already established, adds the three things that document does not cover
(a CONCRETE masking-score definition proposal, an explicit train/validation/target data
separation, and UCE's trigger condition), and marks every open item `pending_owner`. No new
checkpoint download, no new HPCC inference run, no model-weight mutation this round. **Correction
per `review:20261014_221155` R3**: §5 below previously called scPRINT-2 "the standing
single-model default" — `docs/decisions.md` itself (line 779 area) states scPRINT-2 is
**RECOMMENDED as primary, explicitly distinct from "single-model direction is approved"** — this
page now calls it a **candidate**, consistently, not a decided default; only the single-model
DIRECTION (one main model + PCA/expression baseline, no five-model benchmark) is itself decided.

## 1. Primary CANDIDATE: scPRINT-2 / `small-v2.ckpt` — status, restated from `MODEL_SELECTION_v1.md`

| Item | Status | Source |
|---|---|---|
| Checkpoint identity | Local file SHA-256 (`2b586c14...`) matches the official HuggingFace page exactly | `MODEL_SELECTION_v1.md` §"Resolved this round" |
| Human+mouse coverage | `checkpoint_metadata_confirmed` — read from the checkpoint's OWN saved hyperparameters (16-species list including human+mouse), not the README | `MODEL_SELECTION_v1.md`, `model_registry.py` |
| Environment | Dedicated `scprint2_env` conda env, already set up | `MODEL_SELECTION_v1.md` |
| Interface | Official repo (`cantinilab/scPRINT-2`) exposes embeddings and a gene-network output | `github.com/cantinilab/scPRINT-2` (official repo, not independently re-verified this round) |
| Architecture caveat | Two incompatible checkpoint architectures exist (4-layer default vs. 33-layer); this project's prior attempts used the 4-layer checkpoint | `MODEL_SELECTION_v1.md` |

**What this status does NOT mean** (explicit, per `review:20261014_204536` §2): a
confirmed checkpoint/species-metadata identity is NOT a claim that mouse data has actually been
run through it successfully, NOT a claim the gene-network interface's output has been
masking-attributed or validated for this project's purposes, and NOT a ranking of this model as
"better" than any alternative — no comparative run exists.

## 2. Alternate/control: UCE — explicit trigger condition, not a parallel track

Per the owner's standing single-model direction (`docs/decisions.md`), UCE is kept as a
documented alternative with an EXPLICIT trigger, never a parallel benchmark track:

- **Trigger**: UCE is activated only if scPRINT-2 fails a pre-defined runnability/applicability
  gate for this project's actual monocyte data (e.g. a genuine mouse-data run error, a species/
  gene-mapping rejection, or a resource limit this project cannot work around) — not activated
  because a comparative score looks better on one side.

## 3. Masking score: a CONCRETE proposed contract (`s`, `Δ(g)`) — per `review:20261014_221155` R3, sharpened per `review:20261014_231259` C

The prior draft of this page described a masking target only loosely ("whatever genes the
model's mechanism reports as most attribution-bearing"), with no actual scalar `s` and no
computable `Δ(g) = s_full − s_masked(g)`. This section gives the concrete method the review asked
for, with the review's own follow-up ambiguities resolved explicitly below. **This is a method
PROPOSAL for the owner/implementer to review, not an approved fitting run, not yet implemented,
and not an IFN-set-defined target** (avoiding the circularity `review:20261014_204536`
§3 warned about: the positive-control gene set is used ONLY for an after-the-fact overlap check,
per the last bullet below, never to define `s` itself).

1. **Development/validation split, donor-level** (per §4's table below): split the chosen
   cohort's donors into DEVELOPMENT and HELD-OUT sets, fixed BEFORE fitting anything, stratified
   by the severe/mild+moderate label so both sets contain both classes. A single, recorded random
   seed controls the split; the split itself is saved as an artifact (donor IDs in each set), not
   regenerated silently on a later run.
2. **Cell-level embedding → donor-level embedding, explicit, never flattening away the pilot's
   own cell-state goal**: the model (and, separately, PCA) produces a PER-CELL embedding first.
   Donor-level `e` is the MEAN of that donor's own cell embeddings, restricted to the SAME
   monocyte subtype scope `MONOCYTE_PILOT_PREREG_v1.md` §5 eventually freezes — a flat whole-donor
   average is a SUMMARY statistic used only for this linear severity-direction score, never a
   replacement for the pilot's primary cell-state/program discovery (which works at the cell or
   subpopulation level, not this donor-level summary). The masked and unmasked passes for a given
   donor MUST use the exact same cell set (no re-subsampling between the two), so any difference
   in `e` is attributable to the masking itself, not to a different population of cells.
3. **Fit a frozen linear severity-direction score on DEVELOPMENT donors only — ONE primary
   statistic, chosen in advance, not switched after seeing results**: the PRIMARY proposed
   statistic is the mean-difference vector between the severe-group and mild+moderate-group mean
   donor-level embeddings (simpler, fewer assumptions than a fitted classifier at this sample
   size). An L2-regularized logistic regression on the same two groups is a documented ALTERNATIVE,
   not an interchangeable option — if used instead of the mean-difference, that must be decided
   BEFORE fitting, not chosen after comparing which one "looks better." **Freeze `w`, `b`, and any
   feature-scaling constants (e.g. a per-dimension standardization mean/SD) immediately after this
   fit** — none of the three may be recomputed or adjusted using held-out-donor data at any later
   step.
4. **The scalar `s`**: for any donor-level embedding `e` (full or masked), `s(e) = w · e + b`
   (after applying the SAME frozen scaling from step 3). `s_full = s(e_full)`, where `e_full` is
   the aggregated embedding from the model's normal, unmasked forward pass on a held-out donor's
   own cell set.
5. **Masking and `Δ(g)`, via the REAL scPRINT-2 interface, never an assumed one**: before this
   step is implemented, scPRINT-2's OWN official masking/perturbation interface must be checked
   directly (its repo/docs, not assumed from other scFMs' conventions) and the exact, real
   mechanism it actually exposes (e.g. token dropout, zero-imputation, or whatever it genuinely
   supports at its native input resolution) written down here in place of this placeholder — a
   hand-rolled substitute is never acceptable. Once that real mechanism is confirmed: mask/ablate
   gene (or gene set) `g` in the input, recompute the embedding `e_masked(g)` through the SAME
   frozen model and the SAME cell set as step 2, and compute `s_masked(g) = s(e_masked(g))`.
   **`Δ(g) = s_full − s_masked(g)`**. `Δ(g)` is a SIGNED number describing `g`'s contribution TO
   this fitted score — **it is not itself evidence that `g` is transcriptionally up- or
   down-regulated in severe disease, and not a causal claim**; it only says the model's learned
   representation of `g` moved this particular linear projection by this amount.
6. **Signed `Δ(g)` vs. `|Δ(g)|` ranking are DIFFERENT statistics — pick ONE as primary**: ranking
   genes by signed `Δ(g)` answers "which genes push the score toward severe" (directional);
   ranking by `|Δ(g)|` answers "which genes the score depends on most, either direction"
   (magnitude-only). These answer different questions and must not be silently swapped depending
   on which produces a more convenient-looking top list. The PROPOSED primary statistic is signed
   `Δ(g)`, matching this pilot's own directional severe-vs-mild/moderate question; `|Δ(g)|`
   ranking may be reported as a secondary, explicitly-labeled view.
7. **Aggregation statistic, explicit**: for each held-out DONOR (not cell — donor is the
   aggregation unit, consistent with `MONOCYTE_PILOT_PREREG_v1.md`'s "validation unit is the
   donor" rule), compute `Δ(g)` once per donor (using that donor's own `e_full`/`e_masked(g)`),
   then report the MEAN `Δ(g)` across held-out donors AND its dispersion (SD or a bootstrap CI
   across donors) — never a single point estimate presented without spread. Donors are weighted
   equally (one donor, one vote), never weighted by cell count.
8. **Missing genes**: if `g` is outside the model's own vocabulary/gene-mapping, `Δ(g)` is
   `undefined`/excluded for that gene — NEVER silently set to 0 or dropped without being counted
   as a coverage gap in the report. **Gene-set coverage must be reported against TWO distinct
   universes, never conflated**: (a) coverage against the FINAL DE gene universe (whatever gene
   set the actual DE/expression analysis ends up using) and (b) coverage against the genes the
   MODEL can actually measure/attribute a `Δ(g)` for (its own vocabulary/masking-resolution) —
   `gene_sets/README.md`'s reported "96/97 present in Stephenson's raw panel" is neither of these;
   it is a preliminary check against one cohort's raw gene panel only, not the final DE universe
   and not the model's own measurable-gene set.
9. **Random seed**: one fixed, recorded seed controls BOTH the donor split (step 1) and any
   stochastic step in the model's own forward pass (if its masking/embedding is non-deterministic)
   — reported alongside results, not re-rolled between runs without being named as a new run.
10. **PCA comparator, same EVALUATION logic, its OWN fitted parameters**: steps 2–9 are repeated
    IDENTICALLY in STRUCTURE, substituting PCA components (fit on the SAME development donors, in
    PCA's OWN feature space) for the model's embedding — same split, same `Δ`/aggregation
    definition, same primary statistic choice — but PCA fits its OWN `w_pca`/`b_pca`/scaling in
    its own space; **it never reuses the model's `w`/`b`, which live in a different embedding
    space entirely.** The model-vs-PCA comparison differs only in which embedding feeds the SAME
    evaluation method, never in a shared coefficient vector.
11. **Positive-control overlap check, AFTER the fact only**: once `Δ(g)` is computed for every
    gene in-vocabulary, check whether the TOP-ranked genes by signed `Δ(g)` (direction per
    `POSITIVE_CONTROL_EVIDENCE_v1.md`'s eventual frozen direction) overlap with the frozen
    positive-control gene set — an enrichment/overlap statistic computed on top of an
    independently-derived ranking, never a definition that uses the positive-control set to decide
    what counts as high-`Δ(g)` in the first place. Per `POSITIVE_CONTROL_EVIDENCE_v1.md`'s own
    R3 correction: a null/non-overlapping result here is `indeterminate`, not automatically an
    engineering failure. **2026-10-15 update**: the frozen `HALLMARK_INTERFERON_ALPHA_RESPONSE`
    set used here is explicitly NOT the 99-symbol GO:0034340-derived list Stephenson et al.
    actually scored for Fig. 2f/Table 7 (only 24/99 symbols overlap — see
    `results/stephenson_author_code_investigation_v1.md` §1 and `gene_sets/README.md`); this step
    remains a related-program calibration check against a project-chosen asset, not a
    reproduction of the paper's own gene list, until Cindy decides otherwise.

None of steps 1–11 has been run or implemented this round — this is the method contract for a
future, separately-authorized fitting step.

## 4. Target / training / validation data — explicitly separated, each `pending_owner` until chosen

| Role | Proposed source | Status |
|---|---|---|
| Target data (what the real pilot analyzes) | Monocyte-annotated cells from the chosen cohort (`MONOCYTE_PILOT_PREREG_v1.md` §2) | `pending_owner` (cohort choice) |
| Model weights (training data behind the checkpoint) | scPRINT-2 `small-v2.ckpt`'s own pretraining corpus. Per `review:20261014_231259` C: frozen pretrained weights do NOT by themselves establish that this corpus has NO overlap with this project's own cohort data (e.g. public COVID atlases could in principle overlap a general pretraining corpus) — this is `unknown`, explicitly, never asserted as "no overlap" without having checked | unchanged from `MODEL_SELECTION_v1.md`; overlap status `unknown` |
| Validation data (for any model-selection/masking-score judgment call) | A donor-level held-out split WITHIN the target cohort, per `MONOCYTE_PILOT_PREREG_v1.md` §8's "validation unit is the donor" rule | `pending_owner` — no split has been constructed |

Keeping these three roles separate matters specifically because `review:20261014_204536`
§2.3 warns that tuning a model choice ON the positive control and then treating the SAME result
as independent validation overstates performance — if a held-out donor/dependency split cannot
be built (too few donors), the pilot must be labeled exploratory, not silently reuse the full
dataset for both roles.

## 5. What this page does NOT do

Does not install, download, or re-verify any checkpoint. Does not run inference. Does not pick
scPRINT-2 over UCE as a final decision — scPRINT-2 is a RECOMMENDED candidate, not a confirmed
default (`docs/decisions.md`'s own wording distinguishes "scPRINT-2 recommended as primary" from
"single-model direction is approved"; only the latter is actually decided). Does not implement
the concrete masking-score contract in §3. Does not construct the proposed validation split.
Checkpoint confirmation remains `pending_owner`, stated plainly, not implied as already settled.

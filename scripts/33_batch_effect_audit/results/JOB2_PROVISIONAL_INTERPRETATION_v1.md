# Job 2 (SLURM 29868753) — provisional result interpretation, v1 (2026-10-11)

Read-only interpretation of an already-executed, already-accepted run. Produced in response to a
review (`codex_output/20261011_183017.md`) of `claude-code_output/20261011_153652.md`, which found
that handoff's prose (a) reported only `neighbor_batch_mixing`'s improvement and omitted that
Platelet's stratified batch ASW got slightly *worse* post-Harmony, (b) rounded the mixing
improvement to "roughly quadrupled in both lineages" when the two lineages' ratios differ
materially, and (c) stated "all 40 of `CV0201`'s cells" in a way that could be misread as 40 per
lineage rather than 40 total across both. This document corrects all three by reporting the exact
numbers directly from `manifest.json`/`metrics_long.csv`, re-read from the live output this round
(not recomputed, not re-run).

No new SLURM job, no canonical/QC/DecontX change, no new allocation. Source identity (unchanged
from the accepted run): application `aea9a6197539ddd5086e46ddd34561617ad250d7`, input
`e2da0e4f44350a1ecb28c1894308812f9234c331420672a52177f9445305d6ba`, manifest
`code_revision_verified=true`.

**Same-round correction (2026-10-11, in place — no new version number):** a follow-up review
(`codex_output/20261011_190243.md` §2) found two phrasing issues here: §1's summary folded
`bio_ari`'s near-zero value into "improved on every metric," implying a biology-preservation
claim §3 explicitly disclaims; and "higher = more batch separation" asserted a geometric
interpretation the metric's own implementation (`mean(abs(signed silhouette))`, sign discarded)
cannot support on its own. Both fixed in place below; no number changed.

## 1. Per-lineage metrics, both directions reported

Source: `metrics_long.csv` and `manifest.json` under
`/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/batch_effect_audit/provisional_diagnostic_input_v3/diagnostic_run/`,
re-read in full this round.

### Monocyte_Myeloid (n_diagnostic_cells_bio_known = 11,959 / 12,998 technical)

| Metric | Uncorrected | Harmony | Change | Direction |
|---|---:|---:|---:|---|
| `neighbor_batch_mixing` (all 12,998 cells) | 0.061456 | 0.440458 | +0.379001 | improved (mixing up; 7.17× uncorrected; `expected_under_random_mixing`=0.799925) |
| `batch_asw_neg_scib` (1 − mean abs silhouette) | 0.844430 | 0.897151 | +0.052721 | improved |
| `batch_asw_abs_silhouette`, fatal stratum | 0.164911 | 0.096654 | −0.068257 | improved (lower = better by this metric's convention) |
| `batch_asw_abs_silhouette`, survived stratum | 0.146229 | 0.109044 | −0.037185 | improved |
| `bio_ari` (KMeans-2 vs fatal/survived) | −0.015894 | −0.004029 | +0.011865 | near zero both ways — see §3 |

**Monocyte_Myeloid's two batch-mixing metric families (neighbor mixing, stratified batch ASW)
both improved**, including both outcome strata of batch ASW individually. `bio_ari` moving closer
to zero is not itself evidence of improved biology preservation — see §3; it is reported here as
a direction, not folded into "improved."

### Platelet (n_diagnostic_cells_bio_known = 10,824 / 11,843 technical)

| Metric | Uncorrected | Harmony | Change | Direction |
|---|---:|---:|---:|---|
| `neighbor_batch_mixing` (all 11,843 cells) | 0.103746 | 0.407177 | +0.303431 | improved (3.92× uncorrected; `expected_under_random_mixing`=0.798414) |
| `batch_asw_neg_scib` | 0.859078 | 0.849220 | **−0.009858** | **worse** |
| `batch_asw_abs_silhouette`, fatal stratum | 0.139818 | 0.158269 | **+0.018451** | **worse** |
| `batch_asw_abs_silhouette`, survived stratum | 0.142027 | 0.143291 | **+0.001264** | **worse** (small) |
| `bio_ari` | −0.006810 | −0.008879 | −0.002069 | near zero both ways — see §3 |

**Platelet's neighborhood-mixing diagnostic improved, but its stratified batch-ASW diagnostic got
worse in both outcome strata** (fatal moved more than survived). These are two different
instruments measuring different structure (local kNN composition vs. global embedding silhouette
within an outcome stratum); they disagree here, and that disagreement is the finding, not an
error in either metric. `batch_asw_abs_silhouette` is `mean(abs(per-cell signed silhouette))`
(`scripts/00_shared/code/batch_asw.py:113`) — the sign is discarded before averaging, so this run
cannot distinguish "cells moved further from their own cohort" from "cells' nearest-cohort
assignment became more strongly inverted" as the cause of the increase; both would raise the
same abs-mean. Do not read "higher abs-silhouette" as proof that cohorts became more geometrically
separated — only that this specific metric's value moved in the worse direction by its own
convention. A signed-silhouette breakdown was not computed this round.

## 2. Corrected statement — no single verdict across lineages or metrics

**Harmony did not uniformly improve batch mixing across both lineages and both metric families in
this run.** The accurate statement is: neighborhood mixing improved substantially in both
lineages (differently — 7.17× for Monocyte_Myeloid vs. 3.92× for Platelet, not "roughly
quadrupled in both," which overstates uniformity); outcome-stratified batch ASW improved in
Monocyte_Myeloid but got slightly worse in Platelet. This is not evidence of over-correction by
itself — there is no uncertainty estimate (no bootstrap/repeat run) and the two metrics have
different structure (local-neighborhood vs. global per-stratum silhouette) and different cell
populations (`neighbor_batch_mixing` uses all technical cells; batch ASW uses only outcome-known
cells, 11,959/12,998 and 10,824/11,843 respectively). No parameter should be re-tuned, and QC
should not be touched, on the basis of this single descriptive comparison.

## 3. `bio_ari` answers a narrow question

`bio_conservation()`'s ARI compares a 2-cluster KMeans partition of the embedding against the
fatal/survived label. Both uncorrected and Harmony-corrected values are near zero in both
lineages (range −0.0159 to −0.0040), meaning this specific clustering has weak-to-no
correspondence with the outcome label either way. This does **not** mean "no biological signal,"
"Harmony preserved biology," or "Harmony destroyed biology" — it is one clustering algorithm
against one coarse binary label, not a marker-gene, cell-type, or survival-prediction validation.
No other correction method (scVI, RPCA) and no held-out-cohort transform were run this round;
Harmony here is fit jointly on this diagnostic subsample (transductive/exploratory, per the
manifest's own label), not an out-of-sample projection.

## 4. `expected_under_random_mixing` is a reference, not a target

Manifest-reported random-mixing reference: 0.799925 (Monocyte_Myeloid), 0.798414 (Platelet). Both
corrected values (0.4405, 0.4072) remain well below these references. This gap should not be read
as "correction is incomplete" or used to justify stronger correction — full random mixing across
cohorts is not necessarily the correct target if cohorts differ in real biological composition
(disease mix, timepoint, severity distribution). No action is implied by this gap alone.

## 5. Outcome-eligibility population, exact

| Lineage | Technical diagnostic cells | Outcome-known (bio-scoring eligible) | Outcome-unknown |
|---|---:|---:|---:|
| Monocyte_Myeloid | 12,998 | 11,959 | 1,039 |
| Platelet | 11,843 | 10,824 | 1,019 |

Disputed donor `CV0201`'s 40 selected Job-1 cells split **by lineage, not 40 in each**:
`Monocyte_Myeloid` 20 (2 from `BGCV06`, 18 from `BGCV08`), `Platelet` 20 (all from `BGCV06`) — per
`stephenson_disputed_donor_impact_v2_20261011_151807/disputed_donor_impact.json`'s
`job1_gex_cross_check`, re-read this round. All 40 are outcome-`unknown` (0/40 bio-scoring
eligible in either lineage), so they do not enter the ASW/ARI tables above, though they do
participate in the shared PCA/Harmony fit alongside the eligible cells. `CV0326`'s 33 cells
(`Monocyte_Myeloid` 13, `Platelet` 20) are all outcome-known.

## 6. What this does not do

Does not change the accepted/complete status of Job 2. Does not authorize or imply a production
correction-method decision, a method comparison (Harmony vs. scVI/RPCA), Gate 2 advancement, or
any canonical/`cv_group`/QC/DecontX change. Does not claim statistical significance — all
differences above are single-run descriptive deltas with no confidence interval.

## 7. Visual check (paired cohort/outcome UMAP, existing coordinates only)

Plotted directly from the 4 already-written `*_umap_coords.csv` files (no refit of
UMAP/Harmony/PCA, no count-matrix read) into a new, separate directory —
`diagnostic_run/` itself was not touched. Qualitatively consistent with §1: in
Monocyte_Myeloid, cohort clusters that were visually separated uncorrected (e.g. `Ren 2021`,
`Stephenson 2021`, `COMBAT 2022` as distinct islands) overlap substantially post-Harmony. In
Platelet, most cohorts also overlap post-Harmony, but `Ren 2021` (brown) remains a visually
distinct, only partially-merged cluster even after correction — a qualitative pattern consistent
with (not proof of) the quantitative finding that Platelet's stratified batch ASW did not uniformly
improve. The uncorrected and harmony panels are two independent UMAP fits; axis scale/orientation
is arbitrary and not a displacement comparison between them. `fatal`/`survived`/`unknown` outcome
panels show no obvious spatial clustering by outcome in either scope, consistent with the
near-zero `bio_ari` in §3 (again: this is a visual, non-quantitative observation, not a
statistical test).

Paths:
```
/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/batch_effect_audit/job2_review_artifacts_20261011_183817/
  Monocyte_Myeloid_expression_pca_umap_cohort_outcome.png
  Platelet_expression_pca_umap_cohort_outcome.png
```

## 8. Paths

```
/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/batch_effect_audit/provisional_diagnostic_input_v3/diagnostic_run/
  manifest.json, metrics_long.csv
/bigdata/godziklab/shared/Xinru/302006/06_LOGS/stage33_diag_run_29868753.log
/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/batch_effect_audit/stephenson_disputed_donor_impact_v2_20261011_151807/disputed_donor_impact.json
```

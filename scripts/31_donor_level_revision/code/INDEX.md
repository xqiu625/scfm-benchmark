# Index — `scripts/31_donor_level_revision/code/`

168 files, flat directory, grown iteratively over months. Written 2026-09-14 ahead of a possible
full rebuild (decontaminate upfront + add cohorts), so the reusable protocol pieces are easy to
find and the historical debugging trail doesn't get mistaken for something that needs re-running.
Not physically reorganized (would break every `.sh` wrapper's hardcoded paths and risk breaking
HPCC-side references) — this is a map, not a move.

## Core protocol / infrastructure (reusable, load-bearing)

- `repeated_cv.py` — the 10-seed `StratifiedGroupKFold` repeated-CV harness every patient-level
  baseline in this thread is built on. Start here for any new baseline.
- `masking_attribution_common.py` — shared masking-attribution logic (incl. `top_k_genes_per_cell`)
  used by all five per-model masking scripts.
- `transfer_utils.py` — cross-disease transfer helpers.
- `build_patient_map_v2.py` — canonical `sample_id -> patient_id` mapping.
- `donor_level_reanalysis.py` — patient/donor-level aggregation from cell-level OOF.

## Primary 5-FM masking attribution + GSEA (Section R6 / `main.tex` §3.6)

Per-model masking scripts, their SLURM wrappers, and test/dry-run variants:
`geneformer_masking_attribution*.py/.sh`, `scgpt_masking_attribution*.py/.sh`,
`state_masking_attribution*.py/.sh`, `uce_masking_attribution*.py/.sh`,
`scprint2_masking_attribution*.py/.sh`, `scprint2_multigene_occlusion*.py/.sh`,
`*_run_gsea*.sh`, `geneformer_run_gsea.py`, `compare_masking_attribution_5models.py`
(the cross-model comparison that produces `consensus_ranking.csv`),
`geneformer_gene_attribution.py/.sh`, `state_gene_attribution.py`.
`test_*_masking_only*.sh` are one-off dry-run checks, safe to ignore unless debugging a specific
model's masking step.

## Baseline ladder (this session, §Table 3 candidate additions)

- `pca_baseline_patient_oof_v2.py` (+ `.sh`) — the headline PCA-200 baseline (Table 3, AUC 0.7423).
- `pca_baseline_restricted_patient_oof.py`, `pca_restricted_vs_pca200.py` (+ `run_*.sh`) —
  non-platelet-lineage-marker-removed refit; answers the ambient-composition-confound question
  (0.7472 vs 0.7423, not significantly different).
- `pseudobulk_elasticnet_patient_oof.py` (+ `.sh`, `pseudobulk_en_chunk1/2.sh`),
  `pseudobulk_equivalence_test.py` — PseudobulkEN rung + TOST equivalence test.
- `hvg_marker_baselines_patient_oof_v2.py` (+ `.sh`) — HVG-mean and MarkerScore rungs.
- `platelet_tcell_ratio.py` (+ `.sh`), `pla_t_vs_pca200_paired.py` — composition-ratio rung.
- `baseline_ladder_auc.py`, `baseline_ladder_delong.py` — aggregation + exact DeLong/bootstrap
  across all rungs vs PCA-200.
- `resolution_tool.py` — the named `resolvable_delta_auc()` resolution-analysis function.
- `clinical_utility_and_consensus.py` (+ `.sh`) — PR-AUC, decision curve, consensus predictor.
- `consensus_attribution_signature.py` — builds the k-of-5 consensus gene set
  (531 at k=5, 1918 at k>=4) from `consensus_ranking.csv`.

## Contamination / ambient-RNA / doublet audit (this session, newest work)

- `consensus_gene_contamination_audit.py` (+ `.sh`) — the original 13-flagged-gene raw-count audit.
- `export_for_decontx.py`, `run_decontx.R` (+ `.sh`), `decontx_reaudit.py` (+ `.sh`) — DecontX on
  the platelet-only object.
- `decontx_ambient_profile.R` (+ `.sh`), `decontx_review_followup.py` (+ `.sh`) — genome-wide
  Δrho reference distribution + ambient-profile diagnostics, response to external review.
- `ren_multilineage_keepcols.py` (+ `run_ren_multilineage_*.sh`),
  `build_multicelltype_decontx_object.py` (+ `.sh`), `run_decontx_multicelltype.R` (+ `.sh`),
  `multicelltype_decontx_reaudit.py` (+ `.sh`) — the genuine multi-cell-type DecontX rerun
  (144-patient Ren+Stephenson subset); confirms platelet-only under-identification, shows CD3E's
  mortality signal is a pure artifact, CCL5's survives correction.
- `doublet_detection_and_reaudit.py` (+ `run_doublet_detection.sh`) — Scrublet per-sample doublet
  detection; ~30% rate is implausible (miscalibrated for platelets), excluding doublets barely
  changes the audit.
- `r6_purity_and_gene_audit.py` (+ `.sh`) — per-cell platelet purity, cross-lineage co-expression
  flagging, and the full 22-gene (not just 13) R6 depth/mortality/%mito audit.
- `install_decontx_biocmanager.sh` — the working celda/DecontX install recipe (BiocManager into a
  project-local R lib under `module load R/4.5.2` + `module load gcc/12.2.0`; conda/bioconda
  repeatedly failed, don't retry that route).
- `check_stephenson_categories.py` — one-off diagnostic (Stephenson `full_clustering` category
  list), not part of any pipeline.

## Cross-disease transfer / LOBO / adaptation curve

`cross_disease_96pool.py/.sh`, `cross_disease_fatal_vs_survived_patient*.py`, `jobN_lobo_v2.py/.sh`,
`jobS_cross_disease.py/.sh`, `jobS2_pooled_lobo.py/.sh`, `jobS2_round4_primary.py/.sh`,
`jobS_within_sepsis_cv.py/.sh`, `jobS_pca_baseline.py/.sh`, `jobS_embed` (dir),
`run_adaptation_curve_bootstrap_patient.py`, `jobH_step3/4` wrappers, `run_jobH*.sh`,
`jobQ_paired_inference_full.py`, `jobQ_qc_covid_only.py`, `jobR_winner_test_ensemble.py/.sh`,
`jobL_paired_inference.py/.sh`, `jobK_qc_baseline.py/.sh`.

## Monocyte positive control

`build_monocyte_canonical.py`, `monocyte_source_stratified_auc.py`,
`pca_baseline_monocyte_patient_oof.py`, `rerun_monocyte_mortality_patient.py`,
`jobO_monocyte_matched.py`, `jobO_investigate_flip.py`, `geneformer_monocyte_*`,
`step3_embed_5fm_monocyte.sh`.

## Chongqing sepsis cohort integration

`chongqing_platelet_monocyte_coexpr.py`, `integrate_chongqing_sepsis.py/.sh`,
`rerun_pca200_and_compile_chongqing.py/.sh`, `run_adaptation_curve_chongqing*.py/.sh`.

## Historical statistical-debugging trail (`round*`, most `job*` not listed above)

`round5_*` through `round17_*` (delong exact, incremental value, range statistic, cohort-only
null, max-LR permutation, nested AUC variants, harness unit test, primary-scaling check) and
`rerun_covid_sepsis_mortality_patient*.py/.sh`, `regen_primary_and_range_null_vectors.py`,
`null_model_cellcount.py`, `severity_consistency_audit.py`. These are the project's own
"never trust a commit message, always re-measure" audit trail (see docs/decisions.md) — mostly
superseded by later rounds, kept for provenance, not a starting point for new work.

## One-off verification scripts

`verify_gene_mapping.py`, `verify_patient_mapping.py`, `verify_state_attribution.py`,
`healthy_patient_count.py`, `inspect_scprint2*.py/.sh`. Point-in-time checks, not pipeline steps.

## Figures

`figures_v2/` — figure-generation scripts, separate from the analysis code above.

---

**If starting the rebuild**: the contamination/ambient-RNA section above is the direct
predecessor work — `install_decontx_biocmanager.sh`'s environment recipe, the multi-cell-type
object construction pattern in `build_multicelltype_decontx_object.py`, and the R6 gene audit in
`r6_purity_and_gene_audit.py` are the pieces most worth reusing rather than rebuilding from
scratch. `repeated_cv.py` and `masking_attribution_common.py` are the protocol backbone for
whatever the rebuilt pipeline's baseline ladder and attribution steps end up being.

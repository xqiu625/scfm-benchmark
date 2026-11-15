# Tables v3 (2026-09-09, generated from `figure_table_spec_v3.md` §3-4)

> All numbers read directly from the CSVs cited under each table (source path given); verified against `P0_onepager.md` (round-17 final) where the spec supplied an anchor. No hand-typed statistics except where explicitly marked.

## Table 3 — within-COVID main table

| model | repeat-mean +/- sd | ensemble-OOF AUC | vs PCA-200 (exact DeLong p) | vs null (exact DeLong p) |
|---|---|---|---|---|
| Geneformer | 0.786±0.026 | 0.8067 | 0.0360 | 0.0006 |
| STATE | 0.778±0.033 | 0.7946 | 0.0484 | 0.0014 |
| scPRINT-2 | 0.740±0.027 | 0.7561 | 0.6547 | 0.0113 |
| PCA-200 | 0.725±0.033 | 0.7423 | — | 0.0126 |
| scGPT | 0.720±0.047 | 0.7392 | 0.9089 | 0.0232 |
| UCE | 0.661±0.046 | 0.6765 | 0.0809 | 0.1571 |
| null (count+disease) | 0.554 | — | — | — |
| cohort-only | 0.486±0.038 | — | — | — |
| QC+cell count (as-fit) | 0.618 | — | — | — |

Pre-specified primary (best FM vs PCA-200, discrimination): winner-corrected permutation P=0.181; cross-check recentered bootstrap P=0.073. Confirmatory: leaders (STATE, Geneformer) vs UCE, range-statistic permutation P=0.011.

Source: `within_covid/donor_level_auc_repeated.csv`, `null_model_auc_repeated.csv`, `paired_inference_full/{winner_test_ensemble,delong_exact,each_model_vs_null,cohort_only_baseline}.csv`.

## Table 5 — cross-disease transfer table

| direction | model | LR AUC [95% CI] | RF AUC | severe-excluded LR AUC |
|---|---|---|---|---|
| forward (COVID→Sepsis, 91 test pts) | UCE | 0.684 [0.543, 0.805] | — | — |
| forward (COVID→Sepsis, 91 test pts) | STATE | 0.653 [0.490, 0.808] | — | — |
| forward (COVID→Sepsis, 91 test pts) | scGPT | 0.650 [0.495, 0.780] | — | — |
| forward (COVID→Sepsis, 91 test pts) | scPRINT-2 | 0.647 [0.496, 0.791] | — | — |
| forward (COVID→Sepsis, 91 test pts) | Geneformer | 0.644 [0.486, 0.824] | — | — |
| reverse (Sepsis→COVID, 96-pt train, 240 test pts) | STATE | 0.721 [0.621, 0.809] | 0.711 | 0.765 |
| reverse (Sepsis→COVID, 96-pt train, 240 test pts) | Geneformer | 0.712 [0.612, 0.814] | 0.636 | 0.775 |
| reverse (Sepsis→COVID, 96-pt train, 240 test pts) | UCE | 0.693 [0.596, 0.799] | 0.607 | 0.772 |
| reverse (Sepsis→COVID, 96-pt train, 240 test pts) | scGPT | 0.642 [0.525, 0.765] | 0.622 | 0.678 |
| reverse (Sepsis→COVID, 96-pt train, 240 test pts) | scPRINT-2 | 0.619 [0.493, 0.743] | 0.517 | 0.630 |

Source: `cross_disease_chongqing/cross_disease_chongqing_20260904_233237.csv` (forward), `cross_disease/cross_disease_patient_96pool_20260908_142954.csv` (reverse + severe-excluded sensitivity).

## Table 8 — incremental-value table (QC+count+cohort baseline)

| model | augmented AUC | delta AUC | exact DeLong SE | DeLong p | nested-LR P |
|---|---|---|---|---|---|
| UCE | 0.6430 | +0.0160 | 0.00889 | 0.07196 | P=0.0022 |
| scGPT | 0.6534 | +0.0264 | 0.00967 | 0.00632 | P<0.0001 |
| scPRINT-2 | 0.6733 | +0.0463 | 0.01665 | 0.00539 | P<0.0001 |
| PCA-200 | 0.6772 | +0.0502 | 0.01716 | 0.00342 | P<0.0001 (pre-specified primary) |
| Geneformer | 0.6802 | +0.0532 | 0.01395 | 0.00014 | P<0.0001 |
| STATE | 0.6919 | +0.0649 | 0.01711 | 0.00015 | P<0.0001 |

Baseline (QC+count+cohort) AUC = 0.627. The pre-specified primary is the PCA-200 row (no model selection); the five FMs are exploratory with both readings paired; nested-LR establishes existence, ΔAUC quantifies magnitude.

Source: `paired_inference_full/incremental_value_test.csv`; nested-LR P per P0_onepager.md round-10 (empirical permutation null, not the retired chi²(1) asymptotic reference).

## Table 9 — review-content table (goes into Introduction)

| study | scFM | single-cell input | task | cohort | events | from-scratch baseline | patient-level aggregation | correction | note |
|---|---|---|---|---|---|---|---|---|---|
| EGSP (Bioinformatics Advances 2026, PMC13032892) | scFoundation | embedding transfer (to bulk RNA-seq) | pan-cancer overall survival (TCGA, 25 cancer types) | 8,776 samples | 3,095 death / 5,681 censored (pooled) | yes | yes | no | peer-reviewed |
| Empirical Evaluation of scFMs for Cancer Outcomes (bioRxiv 2025.10.31.685892, PMC12637420) | 9 models | direct | 6 cancer clinical/biological outcome tasks | not disclosed per task (1,170 supervised experiments total) | not disclosed per task | yes | yes | no | preprint |
| TissueFormer (bioRxiv 2025.08.17.670735, PMC13466306) | Murine-adapted Geneformer | direct | COVID-19 severity (3-class) | 376 donors (3 cohorts) | severe: 46/83/28 per cohort | yes | yes | no | preprint |

Search strategy: PubMed (via Google, `site:pubmed.ncbi.nlm.nih.gov`) + bioRxiv/medRxiv, query `"single-cell foundation model" patient-level clinical outcome benchmark survival mortality 2025 2026`, searched 2026-09-09; 4 unique records recalled, 2 included (EGSP, 9-model), 2 excluded (scSurvival: not FM-based; Chen et al. hematology-CPD: not scRNA-seq). TissueFormer's inclusion path was different (identified in an earlier project literature pass, confirmed retrievable via a targeted search, not the general query). None of the three reports selection/multiplicity correction, and none discloses event counts at the granularity of its own comparisons.

## Supplementary tables

### tab:score_correlation

(see `paired_inference/score_correlation_pearson.csv`, 6 rows — reproduced verbatim below)

| Unnamed: 0 | STATE | UCE | Geneformer | scPRINT-2 | scGPT | PCA-200 |
|---|---|---|---|---|---|---|
| STATE | 1.0000 | 0.7944 | 0.9042 | 0.8972 | 0.8158 | 0.8591 |
| UCE | 0.7944 | 1.0000 | 0.7335 | 0.7547 | 0.8433 | 0.7070 |
| Geneformer | 0.9042 | 0.7335 | 1.0000 | 0.8770 | 0.7807 | 0.8303 |
| scPRINT-2 | 0.8972 | 0.7547 | 0.8770 | 1.0000 | 0.7593 | 0.7837 |
| scGPT | 0.8158 | 0.8433 | 0.7807 | 0.7593 | 1.0000 | 0.7756 |
| PCA-200 | 0.8591 | 0.7070 | 0.8303 | 0.7837 | 0.7756 | 1.0000 |

### tab:winner_ensemble

Source: `scripts/31_donor_level_revision/results_patient_v2/paired_inference_full/winner_test_ensemble.csv` (6 rows). See CSV directly for full precision; columns: row_type, model, auc_model, auc_pca, delta_auc, ci_lo, ci_hi, p_indiv_uncorrected, B, rng_seed.

### tab:se_check (DeLong vs bootstrap)

Source: `scripts/31_donor_level_revision/results_patient_v2/paired_inference_full/delong_exact.csv` (15 rows). See CSV directly for full precision; columns: pair, delta_auc, delong_exact_SE, delong_exact_p.

### tab:lobo_folds

Source: `scripts/31_donor_level_revision/results_patient_v2/lobo_v2/lobo_per_cohort_auc.csv` (25 rows). See CSV directly for full precision; columns: fm, cohort_data_no, held_out_auc, n_patients_held_out, n_fatal_held_out.

### tab:qc_strat

Source: `scripts/31_donor_level_revision/results_patient_v2/qc_baseline/patient_qc_distributions_covid_only.csv` (240 rows). See CSV directly for full precision; columns: patient_id, y, disease, n_cells, median_counts, median_genes, median_pct_mito.

### tab:fold_balance

Source: `scripts/31_donor_level_revision/results_patient_v2/paired_inference_full/fold_cohort_balance_seed42.csv` (5 rows). See CSV directly for full precision; columns: fold, n_test, n_fatal_test, cohorts_in_test.

### tab:qc_by_cohort

Source: `scripts/31_donor_level_revision/results_patient_v2/paired_inference_full/qc_by_cohort_association.csv` (9 rows). See CSV directly for full precision; columns: cohort, n, n_fatal, note, median_counts_fatal, median_counts_surv, mwu_p.

### tab:cohort_rates

Source: `scripts/31_donor_level_revision/results_patient_v2/paired_inference_full/item4_per_cohort_table.csv` (9 rows). See CSV directly for full precision; columns: cohort, n, n_fatal, raw_rate, laplace_rate.

### tab:incremental diagnostics (univariate)

Source: `scripts/31_donor_level_revision/results_patient_v2/paired_inference_full/item3_univariate_association.csv` (4 rows). See CSV directly for full precision; columns: feature, spearman_rho, spearman_p.

### tab:incremental diagnostics (coef signs)

Source: `scripts/31_donor_level_revision/results_patient_v2/paired_inference_full/item3_multivariate_coef_signs.csv` (4 rows). See CSV directly for full precision; columns: feature, mean_coef, pct_positive, n_fits.

### tab:severe_sens

Source: `scripts/31_donor_level_revision/results_patient_v2/cross_disease/cross_disease_patient_96pool_20260908_142954.csv`, `sens_severe_excluded_*` columns (see Table 5 above, which already includes the severe-excluded LR column).

### tab:aggregation_sens

Source: `scripts/31_donor_level_revision/results_patient_v2/within_covid/donor_level_auc_repeated.csv`, `mean`/`median` aggregation rows (both already loaded for Table 3; median rows available in the CSV but not reproduced here).

### tab:cellcount

Source: `scripts/31_donor_level_revision/results_patient_v2/within_covid/cellcount_association.csv` (6 rows). Columns: fm, spearman_rho, spearman_p, auc_low_count, n_low, n_fatal_low, auc_mid_count, n_mid, n_fatal_mid, auc_high_count, n_high, n_fatal_high.

### tab:nested_oof (3×5) and tab:nested_deviance

**No local CSV — transcribed from `P0_onepager.md` §3a round-16 section (already cross-checked against P0 in this session; not re-verified here).**

nested_oof, three reduced-model OOF AUCs (StandardScaler-corrected harness): QC+cohort+PCA-200 = 0.7098; PCA-200 alone = 0.7241; PCA-200+log1p(n_cells) = 0.7131. 0/5 models pass BH in every configuration (see P0 §3a for the full rank-by-rank table).

nested_deviance (in-sample, superseded to a supplementary cross-check per round-17): permutation null percentiles — median 0.188, 95th 0.790, 99th 1.230, max 3.819; observed statistics — STATE 4.3284, Geneformer 3.2227, scPRINT-2 3.1843, scGPT 1.5365, UCE 0.8749. Carries no headline weight (in-sample construction, not independently OOF-validated).

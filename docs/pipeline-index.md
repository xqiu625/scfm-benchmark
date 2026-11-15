# Pipeline Index

Central index of all pipeline stages. Stage folders live under `scripts/NN_*/`
and keep this project's **accession numbering** (ACC) rather than the template's
generic 00–10 numbering — the ACC numbers are referenced directly by the
manuscript (see `scripts/MANUSCRIPT_MAPPING.md`), so renumbering them would break
those references. The generic-role column below maps each ACC stage onto the
architecture's pipeline roles.

Status values: `planned` · `running` · `done` · `failed` · `needs-rerun` · `deprecated`

> **Status scope (2026-09-21):** `done` below records completion of the
> pre-rebuild stage. Stage 32 is rebuilding the benchmark input, so downstream
> stages will require refresh where P1/P2 changes their inputs. The active state
> and critical path live in `docs/current-status.md` and `docs/execution-plan.md`.

## Numbered analysis stages

| ACC | Folder | Generic role | Purpose | Manuscript role | Status |
|----:|---|---|---|---|---|
| 00 | `scripts/00_setup/` | 00_setup | HPCC result sync, cell-index inspection, ACC17–21 figure gen | infra | done |
| 00s | `scripts/00_shared/` | shared | Embedding sanity checks reused across stages | AUX | done |
| 01 | `scripts/01_embedding_integrity/` | 01_data_prep | Embedding integrity / batch-effect QC | SUPP | done |
| 02 | `scripts/02_embedding_classification/` | 05_modeling | Fine-tuned 5-FM classification benchmark | Fig 1 / Table 1 | done |
| 04 | `scripts/04_pca_baseline/` | 05_modeling | PCA baseline (elevated to main comparison) | `tab:baseline` | done |
| 05 | `scripts/05_sample_efficiency/` | 07_validation | Sample-efficiency robustness curves | Fig 4b | done |
| 06 | `scripts/06_batch_generalization/` | 06_integration | Cross-disease generalization (STATE bidirectional) | `tab:crossdisease` | done |
| 07 | `scripts/07_gene_dropout/` | 07_validation | Gene-dropout robustness | Fig 4c | done |
| 08 | `scripts/08_statistical_significance/` | 07_validation | Significance testing across FMs | SUPP | done |
| 09 | `scripts/09_embedding_ablation/` | 07_validation | Embedding ablation | SUPP | done |
| 10 | `scripts/10_clinical_metrics/` | 07_validation | Clinical metric computation | SUPP | done |
| 11 | `scripts/11_interpretability/` | 11_interpretability | Interpretability (SHAP/GSEA) | `tab:per_fm_*` | done |
| 12 | `scripts/12_umap_visualization/` | 02_eda | UMAP + cluster-quality metrics | figure | done |
| 15 | `scripts/15_differential_diagnosis/` | 05_modeling | 6-class differential diagnosis | Fig 2–3 / Table 2 | done |
| 17 | `scripts/17_domain_adapted_fm/` | 06_integration | Domain-adapted FM (Geneformer LoRA chain) | `tab:domainadapt` | done |
| 18 | `scripts/18_multimodal_fusion/` | 06_integration | Attention-based multimodal fusion | `tab:fusion` | done |
| 19 | `scripts/19_model_optimization/` | 05_modeling | Model/hyperparameter optimization | SUPP | done |
| 21 | `scripts/21_conditional_generation/` | 08_prediction | Conditional generation | EXPLORATORY | done |
| 22 | `scripts/22_esm_protein_multimodal/` | 04_feature_building | ESM protein multimodal features | SUPP | done |
| 24 | `scripts/24_option_c_canonical_input/` | 05_modeling | Canonical embedding pipeline + primary mortality benchmark | many `tab:*` | done |
| 26 | `scripts/26_figures/` | 10_reporting | Main + supplementary figure regeneration | all main figs | done |
| 27 | `scripts/27_covid_sepsis_only_rerun/` | 07_validation | COVID+Sepsis-only rerun | SUPP | done |
| 28 | `scripts/28_umap_covid_sepsis/` | 02_eda | UMAP on COVID+Sepsis subset | figure | done |
| 29 | `scripts/29_adaptation_curve/` | 07_validation | Adaptation-curve (Y-int/slope/asymptote) | `tab:adaptation_curve` | done |
| 30 | `scripts/30_review_response_stats/` | 07_validation | Clustered inference and reviewer-response reruns | revision evidence | done; refresh after rebuilt input |
| 31 | `scripts/31_donor_level_revision/` | 07_validation | Patient mapping, repeated CV, attribution and GSEA | revision evidence | done; refresh after rebuilt input |
| 32 | `scripts/32_rebuild_upfront_qc/` | 03_preprocessing | Eight-cohort Gate 1/Gate 2 rebuild, DecontX refit and assembly | upstream of current benchmark | **running / blocked at Gate 2** |
| 99 | `scripts/99_audit/` | 99_audit | Traceability / reproducibility audit (no new science) | audit | planned |

## Named tooling (reusable runners, not single stages)

| Folder | Generic role | Purpose |
|---|---|---|
| `scripts/foundation_models/` | 03_preprocessing → embeddings | Per-FM install + embedding generation (Geneformer, scGPT, UCE, scPRINT, TranscriptFormer) |
| `scripts/finetuning/` | 05_modeling | PEFT/LoRA/full fine-tuning submitters per FM |
| `scripts/tracking/` | reporting | W&B experiment tracking utilities |
| `scripts/agent/` | reporting | Benchmark agent helper |

## Maintenance

- When adding a stage, append a row here **before** writing its README.
- When a stage's numbers feed the manuscript, also add rows to
  `docs/claim-traceability.md` and `docs/figure-table-inventory.md`.
- Run `scripts/99_audit` to verify every path in this table exists.

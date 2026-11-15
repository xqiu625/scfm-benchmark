# scripts/ — Numbered workflow stages (Project 302006)

## Stress-Testing Single-Cell Foundation Models Across Disease and Species

Numbered stage folders (the architecture's `scripts/NN_stage/` layer). Stages keep
this project's **accession numbering**; see `../docs/pipeline-index.md` for the full
stage table and the mapping onto the architecture's generic pipeline roles.

**Project state:** Option C 6-FM pipeline; mortality binary endpoint on the
disease-matched COVID + Sepsis subset (42,095 cells, 322 donors).
**Living state:** `../docs/current-status.md`
**Submission punch list:** `../docs/submission-checklist.md`
**Manuscript ↔ stage mapping:** `MANUSCRIPT_MAPPING.md`
**Reusable runners** (not numbered stages): `foundation_models/`, `finetuning/`,
`tracking/`, `agent/` — described in `../docs/tooling-reference.md`.

---

## Folder convention

```
<ACC>_<stage_name>/
├── code/         # Python scripts, SLURM submitters (run_*.py, submit*.sh)
├── results/      # CSV tables, JSON outputs (committed to git per 2026-05-24 reorg)
├── figures/      # PDF + PNG figures (gitignored; uploaded to Overleaf separately)
├── logs/         # Per-stage SLURM logs (optional; HPCC is canonical)
├── manifest.json # Run provenance (architecture §9) — to be emitted going forward
└── README.md     # Per-stage documentation
```

---

## Status legend

| Marker | Meaning |
|---|---|
| **CORE** | Primary numbers feed the headline manuscript tables/figures. Should be the standard of evidence. |
| **SUPP** | Supports a supplementary table or referenced caveat in the manuscript. |
| **AUX** | Methodological infrastructure (sanity checks, shared code, figure generation). |
| **LEGACY** | Superseded by Option C pipeline; retained for reproducibility / "before correction" comparison. Numbers may be pre-Option-C and should not be used as new evidence. |
| **EXPLORATORY** | One-off side branches; not in the v1 submission. |

---

## Active analyses (Option C era)

| ACC | Folder | Manuscript role | Status |
|----:|---|---|---|
| 06 | `06_batch_generalization/` | `tab:crossdisease` (STATE bidirectional, STACK anti-correlated under severe=survived label, Job 24404316 2026-05-22) | **CORE** |
| 17 | `17_domain_adapted_fm/` | `tab:domainadapt` (Geneformer LoRA chain Jobs 23853347–23866507, May 6 2026, patient-level CV, no overlap; legacy continued-pretraining retained as documented data-overlap example) | **CORE** |
| 24 | `24_option_c_canonical_input/` | The standardized embedding pipeline + primary mortality benchmark (`tab:performance`, `tab:baseline`, `tab:delong`, `tab:calibration`, `tab:multiclassifier`, `tab:ablation`, `tab:fusion`, `tab:embedding_quality`, `tab:activation`, `tab:per_fm_gsea`, `tab:per_fm_shap`, `tab:per_fm_jaccard`) | **CORE** |
| 29 | `29_adaptation_curve/` | `tab:adaptation_curve` (Y-int / slope / asymptote on cross-disease sepsis; both single-seed legacy and 200-iter bootstrap variants in `results/`) | **CORE** |
| 00 | `00_shared/` | Shared embedding checks and workflow utilities | **AUX** |
| 01 | `01_embedding_integrity/` | Embedding NaN / dimensionality / cell-count checks before downstream evaluation | **AUX** |
| 26 | `26_figures/` | Figure generation scripts | **AUX** |

## Supporting analyses

| ACC | Folder | Notes |
|----:|---|---|
| 05 | `05_sample_efficiency/` | Subsampling at 5–100 % patient fractions | **SUPP** |
| 07 | `07_gene_dropout/` | Gene-masking robustness (10–70 %) | **SUPP** |
| 08 | `08_statistical_significance/` | DeLong / bootstrap CI infrastructure used by primary benchmark | **SUPP** |
| 09 | `09_embedding_ablation/` | PCA / normalization / top-PC removal | **SUPP** |
| 10 | `10_clinical_metrics/` | Sens@Spec, Cohen's κ, AUC-PR | **SUPP** |
| 11 | `11_interpretability/` | Top genes, biomarker analysis | **SUPP** |
| 12 | `12_umap_visualization/` | Per-FM UMAP | **SUPP** |
| 15 | `15_differential_diagnosis/` | 6-class disease classification. **2026-05-24**: silent cell-index-as-group fallback was removed; runs now require a patient-level group column or raise `RuntimeError`. | **SUPP** |
| 18 | `18_multimodal_fusion/` | Late-fusion ensemble | **SUPP** |
| 19 | `19_model_optimization/` | Quantization (FP16 / INT8) | **SUPP** |
| 27 | `27_covid_sepsis_only_rerun/` | Cohort-restriction rerun infrastructure (pre-Option C) | **SUPP** |
| 28 | `28_umap_covid_sepsis/` | COVID+Sepsis UMAP | **SUPP** |

## Legacy / superseded

| ACC | Folder | Reason |
|----:|---|---|
| 02 | `02_embedding_classification/` | Pre-Option-C severity benchmark; replaced by ACC 24 mortality binary on COVID+Sepsis subset. Numbers should not be cited as v1 evidence. | **LEGACY** |
| 03 | `03_zeroshot_benchmark/` | Empirically equivalent to ACC 02 (`diff` showed byte-for-byte identical scripts); not in v1. | **LEGACY** |
| 04 | `04_pca_baseline/` | Pre-Option-C PCA baseline; replaced by `tab:baseline` in ACC 24. | **LEGACY** |
| 13 | `13_shap_interpretability/` | Replaced by per-FM SHAP under ACC 24 (`tab:per_fm_shap`). | **LEGACY** |
| 14 | `14_scgpt_noise_analysis/` | scGPT was dropped from the active 6-FM set on 2026-05-01; retained as anti-correlated reference only. | **LEGACY** |
| 21 | `21_conditional_generation/` | Generative side branch; not in v1. | **EXPLORATORY** |
| 22 | `22_esm_protein_multimodal/` | ESM2 protein-language-model multimodal arm; not in v1 (figures retained). | **EXPLORATORY** |

---

## How to add a new analysis

1. Create `<NN>_<name>/code/`, `results/`, `figures/`, `logs/`, `README.md`.
2. Update this registry with the row + status marker.
3. If the analysis produces a manuscript table or figure, update `MANUSCRIPT_MAPPING.md`.

## Methodology floor

All analyses that produce manuscript-grade numbers must:

- Use **patient-level `StratifiedGroupKFold`** with `Sample_ID` (or equivalent donor column) as the group; cell-level CV is acceptable only when documented as a within-cohort separability check, not a cross-sample generalization test (see `tab:mouse_benchmark`, `tab:mouse_perturbation` infection contrasts for examples).
- Read embeddings via **explicit known `obsm` keys**, never via "first 2D obsm" fallback. The fallback was removed from `fusion_option_c.py`, `calibration_per_cohort.py`, `embedding_quality_option_c.py` on 2026-05-24 (external code-review finding); `primary_benchmark_option_c.py` already enforced this.
- Run an **embedding-sanity check** (cell-pair cosine + per-dimension variance) on every FM output before downstream evaluation. The TF-Metazoa near-miss case study (`ssec:mouse_benchmark` Methodology caveat) documents why this is required.
- For any "domain-adapted" claim, the **adapted weights must not have seen the evaluation cells during adaptation**. The legacy `Geneformer_adapted` (continued-pre-training) row in `tab:domainadapt` is retained explicitly as a documented data-overlap example, **not as a benchmark result**.

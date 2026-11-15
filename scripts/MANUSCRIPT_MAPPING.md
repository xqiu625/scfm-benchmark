# Analyses → Manuscript Mapping

## LMRL 2026 Workshop Submission

**Title:** Beyond Accuracy: A Robustness-Focused Benchmark of Single-Cell Foundation Models for Clinical Deployment

---

## Manuscript Section → Analysis Folder Mapping

```
MANUSCRIPT                          ANALYSES FOLDER
═══════════════════════════════════════════════════════════════

Section 3.1: Task 1 - COVID Severity
├── Figure 1: Model ranking         → 02_embedding_classification/
├── Table 1: Performance metrics    → 02_embedding_classification/
└── Supp: Zero-shot comparison      → 03_zeroshot_benchmark/

Section 3.2: Task 2 - Differential Diagnosis
├── Figure 2: 6×6 Confusion matrix  → 15_differential_diagnosis/
├── Figure 3: Per-disease AUC       → 15_differential_diagnosis/
└── Table 2: Model ranking (6-class)→ 15_differential_diagnosis/

Section 3.3: Robustness Analysis
├── Figure 4a: Distribution shift   → 06_batch_generalization/
├── Figure 4b: Sample efficiency    → 05_sample_efficiency/
├── Figure 4c: Gene dropout         → 07_gene_dropout/
└── Table 3: Robustness summary     → (aggregate)

Section 3.4: Embedding Quality
├── Figure 5: Ablation heatmap      → 09_embedding_ablation/
├── Table 4: PCA effects            → 09_embedding_ablation/
└── Supp: scGPT noise analysis      → 14_scgpt_noise_analysis/

Section 3.5: Interpretability (Discussion)
├── Figure 6: Top biomarkers        → 11_interpretability/
└── Supp: SHAP analysis             → 13_shap_interpretability/

Methods Section
├── Table M1: Dataset statistics    → 01_embedding_integrity/
├── Table M2: Model specifications  → 01_embedding_integrity/
└── Statistical significance        → 08_statistical_significance/
```

---

## Reorganized Analysis Structure

### Core Analyses (Required for Submission)

| ACC | Folder | Manuscript Section | Status | Priority |
|-----|--------|-------------------|--------|----------|
| 01 | `01_embedding_integrity/` | Methods | ✅ Done | - |
| 02 | `02_embedding_classification/` | Results 3.1 | ✅ Done | - |
| 05 | `05_sample_efficiency/` | Results 3.3 | ✅ Done | - |
| 06 | `06_batch_generalization/` | Results 3.3 | ✅ Done | - |
| 07 | `07_gene_dropout/` | Results 3.3 | ✅ Done | - |
| 08 | `08_statistical_significance/` | Methods | ✅ Done | - |
| 09 | `09_embedding_ablation/` | Results 3.4 | ✅ Done | - |
| **15** | `15_differential_diagnosis/` | **Results 3.2** | ⏳ **Pending** | **Critical** |

### Supplementary Analyses

| ACC | Folder | Use | Status |
|-----|--------|-----|--------|
| 03 | `03_zeroshot_benchmark/` | Supp: Zero-shot vs Fine-tuned | ✅ Done |
| 04 | `04_pca_baseline/` | Supp: Baseline comparison | ✅ Done |
| 10 | `10_clinical_metrics/` | Supp: Clinical utility | ✅ Done |
| 11 | `11_interpretability/` | Discussion: Biomarkers | ✅ Done |
| 12 | `12_umap_visualization/` | Supp: Embedding viz | ✅ Done |
| 13 | `13_shap_interpretability/` | Supp: SHAP | ⏳ Pending |
| 14 | `14_scgpt_noise_analysis/` | Supp: Noise characterization | ⏳ Pending |

---

## Figure Specifications

### Main Figures (4-5 for workshop paper)

#### Figure 1: Task 1 Model Ranking
**Source:** `02_embedding_classification/figures/`
**Content:** Bar chart showing 5 models × 2 tasks (Binary, 3-Class)
**Key Message:** STATE achieves best performance (0.944 AUC)
```
fig_01_task1_model_ranking.png
├── X-axis: Models (STATE, UCE, Geneformer, TranscriptFormer, scGPT)
├── Y-axis: AUC-ROC
├── Grouped bars: Binary (blue), 3-Class (orange)
└── Error bars: 95% CI from bootstrap
```

#### Figure 2: Task 2 Confusion Matrix (KEY FIGURE)
**Source:** `15_differential_diagnosis/figures/`
**Content:** 6×6 heatmap showing classification confusion
**Key Message:** Which diseases are confused with each other
```
fig_02_differential_diagnosis_confusion.png
├── Rows: True labels (Healthy, COVID, Sepsis, SLE, Flu, Pneumonia)
├── Cols: Predicted labels
├── Color: Prediction rate (0-1)
└── Annotations: Cell counts
```

#### Figure 3: Robustness Multi-Panel (KEY FIGURE)
**Source:** `05_sample_efficiency/`, `06_batch_generalization/`, `07_gene_dropout/`
**Content:** 3-panel figure showing robustness axes
**Key Message:** scGPT fails under distribution shift (0.524 AUC)
```
fig_03_robustness_multipanel.png
├── Panel A: Distribution Shift (COVID → Sepsis)
│   └── Bar chart: In-dist vs OOD AUC per model
├── Panel B: Sample Efficiency
│   └── Line plot: AUC vs % training data
└── Panel C: Gene Dropout
    └── Line plot: AUC vs % genes masked
```

#### Figure 4: Embedding Ablation Heatmap
**Source:** `09_embedding_ablation/figures/`
**Content:** Heatmap showing PCA effects across models
**Key Message:** scGPT improves +7% with PCA (noise in embeddings)
```
fig_04_embedding_ablation.png
├── Rows: Models
├── Cols: Ablation types (Original, PCA-64, PCA-128, PCA-256, PCA-512)
├── Color: AUC
└── Highlight: scGPT row showing improvement
```

#### Figure 5: Model Selection Flowchart (Discussion)
**Source:** Create new
**Content:** Decision tree for model selection
**Key Message:** Practical guidance for practitioners
```
fig_05_model_selection_guide.png
├── Decision nodes: Data size? OOD expected? Sequencing quality?
└── Leaf nodes: Recommended model
```

---

## Required Figure Updates

### Figures to Generate/Update

| Figure | Analysis Folder | Action | Priority |
|--------|----------------|--------|----------|
| fig_01 | 02_embedding_classification | Update with CI error bars | Medium |
| fig_02 | 15_differential_diagnosis | **Create new** | **Critical** |
| fig_03 | Multi-source | **Create combined panel** | **Critical** |
| fig_04 | 09_embedding_ablation | Update to heatmap format | Medium |
| fig_05 | New | **Create flowchart** | Low |

---

## Analysis Folder Updates Needed

### 02_embedding_classification/figures/
```
Required:
├── fig_01_task1_model_ranking.png     # Main figure
├── fig_01_task1_model_ranking.pdf     # Vector version
└── table_01_performance_metrics.csv   # For Table 1
```

### 15_differential_diagnosis/figures/ (NEW - CRITICAL)
```
Required:
├── fig_02_confusion_matrix_STATE.png  # Best model confusion
├── fig_02_confusion_matrix_all.png    # All models comparison
├── fig_03_per_disease_auc.png         # Bar chart per disease
└── table_02_differential_diagnosis.csv # For Table 2
```

### Combined Robustness Figure
```
Create new script to generate:
├── fig_03_robustness_multipanel.png
└── fig_03_robustness_multipanel.pdf

Combines data from:
├── 05_sample_efficiency/results/
├── 06_batch_generalization/results/
└── 07_gene_dropout/results/
```

---

## Scripts to Create

### 1. `create_manuscript_figures.py`
Generates all main figures in consistent style for LMRL submission.

### 2. `run_15_differential_diagnosis.py`
Runs the 6-class differential diagnosis analysis.

### 3. `create_robustness_multipanel.py`
Combines robustness results into single multi-panel figure.

---

## Submission Checklist by Analysis

| Analysis | Data Ready | Code Ready | Figure Ready | Table Ready |
|----------|------------|------------|--------------|-------------|
| 02_embedding_classification | ✅ | ✅ | ⚠️ Update | ✅ |
| 05_sample_efficiency | ✅ | ✅ | ⚠️ Combine | ✅ |
| 06_batch_generalization | ✅ | ✅ | ⚠️ Combine | ✅ |
| 07_gene_dropout | ✅ | ✅ | ⚠️ Combine | ✅ |
| 09_embedding_ablation | ✅ | ✅ | ⚠️ Update | ✅ |
| **15_differential_diagnosis** | ❌ **Run** | ✅ | ❌ **Create** | ❌ **Create** |

---

## Priority Actions

1. **[CRITICAL]** Run `15_differential_diagnosis` analysis
2. **[CRITICAL]** Create Figure 2 (confusion matrix)
3. **[HIGH]** Create Figure 3 (robustness multi-panel)
4. **[MEDIUM]** Update Figure 1 with error bars
5. **[MEDIUM]** Create Figure 4 (ablation heatmap)
6. **[LOW]** Create Figure 5 (flowchart)

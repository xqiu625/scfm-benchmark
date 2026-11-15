# 02 - Fine-tuned Benchmark

**Accession:** 02
**Status:** ✅ COMPLETE
**Last Run:** January 5, 2026

---

## Purpose

Evaluate fine-tuned classification performance of all 5 foundation models on:
- Binary task: Healthy vs COVID-19
- 3-Class task: Healthy vs Mild vs Severe

---

## Input

| File | Description |
|------|-------------|
| `embeddings/*_embeddings.h5ad` | All 5 model embeddings |
| Metadata | Disease state labels |

---

## Output

### Results
| File | Description |
|------|-------------|
| `results/table_02_finetuned_performance.csv` | AUC scores for all models |
| `results/table_02_confusion_matrices.csv` | Per-model confusion matrices |

### Figures
| File | Description |
|------|-------------|
| `figures/fig_02_model_rankings.png` | Bar chart of model AUCs |
| `figures/fig_02_roc_curves.png` | ROC curves for all models |
| `figures/fig_02_confusion_matrix.png` | Confusion matrix heatmaps |

---

## Code

| File | Description |
|------|-------------|
| `code/run_02_embedding_classification.py` | Main benchmark script |
| `code/submit_02_benchmark.sh` | SLURM job script |

---

## Key Results

### Binary Classification
| Rank | Model | AUC |
|:----:|-------|:---:|
| 1 | **STATE** | **0.944** |
| 2 | UCE | 0.910 |
| 3 | Geneformer | 0.845 |
| 4 | TranscriptFormer | 0.838 |
| 5 | scGPT | 0.735 |

### 3-Class Classification
| Rank | Model | AUC |
|:----:|-------|:---:|
| 1 | **STATE** | **0.948** |
| 2 | UCE | 0.912 |
| 3 | Geneformer | 0.851 |
| 4 | TranscriptFormer | 0.851 |
| 5 | scGPT | 0.738 |

**STATE achieves best performance on both tasks.**

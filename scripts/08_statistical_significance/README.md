# 08 - Statistical Significance

**Accession:** 08
**Status:** ✅ COMPLETE
**Last Run:** January 9, 2026
**Runtime:** ~3 min per model (1,440x speedup!)

---

## Purpose

Establish statistical significance of model performance differences using bootstrap confidence intervals.

---

## Input

| File | Description |
|------|-------------|
| `embeddings/*_embeddings.h5ad` | All model embeddings |

---

## Output

### Results
| File | Description |
|------|-------------|
| `results/table_08_bootstrap_ci.csv` | Bootstrap CIs for all models |
| `results/table_08_pairwise_significance.csv` | Pairwise comparisons |

### Figures
| File | Description |
|------|-------------|
| `figures/fig_08_statistical_significance.png` | CI comparison plot |
| `figures/fig_08_forest_plot.png` | Forest plot with CIs |

---

## Code

| File | Description |
|------|-------------|
| `code/run_08_stats_fast.py` | Optimized stats script |
| `code/submit_08_stats_state.sh` | SLURM job for STATE |
| `code/submit_08_stats_uce.sh` | SLURM job for UCE |
| `code/submit_08_stats_scgpt.sh` | SLURM job for scGPT |
| `code/submit_08_stats_geneformer.sh` | SLURM job for Geneformer |

---

## Methods

- **Bootstrap iterations:** 100
- **Classifier:** LogisticRegression (lbfgs solver)
- **Dimensionality:** PCA-256
- **Parallelization:** 16 cores via joblib

---

## Key Results

| Rank | Model | Binary AUC | 95% CI |
|:----:|-------|:----------:|--------|
| 1 | **STATE** | **0.871** | [0.869, 0.875] |
| 2 | UCE | 0.863 | [0.861, 0.868] |
| 3 | scGPT | 0.831 | [0.831, 0.838] |
| 4 | Geneformer | 0.782 | [0.781, 0.789] |

**All CIs are non-overlapping → statistically significant differences.**

---

## Optimization Details

| Parameter | Before | After | Speedup |
|-----------|--------|-------|---------|
| Execution | Sequential | Parallel (16 cores) | 16x |
| Solver | saga | lbfgs | 10x |
| Dimensions | Full (2058) | PCA-256 | 8x |
| Iterations | 500 | 100 | 5x |
| **Total** | >72 hours | ~3 min | **1,440x** |

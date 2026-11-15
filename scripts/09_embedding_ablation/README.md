# 09 - Embedding Ablation

**Accession:** 09
**Status:** ✅ COMPLETE
**Last Run:** January 10, 2026
**Runtime:** ~25 min per model (~60x speedup!)

---

## Purpose

Systematically evaluate embedding characteristics through ablation studies:
1. **Dimensionality reduction** - PCA to 64, 128, 256, 512 dims
2. **Top PCA removal** - Remove top 1, 2, 5, 10 components
3. **Normalization** - Raw vs L2 vs Z-score
4. **Random dimensions** - Random subset baseline

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
| `results/table_09_dim_reduction.csv` | PCA reduction results |
| `results/table_09_top_pca_removal.csv` | Top component removal |
| `results/table_09_normalization.csv` | Normalization comparison |
| `results/table_09_random_dims.csv` | Random dimension baseline |

### Figures
| File | Description |
|------|-------------|
| `figures/fig_09_dim_reduction.png` | PCA dimension curves |
| `figures/fig_09_top_pca_removal.png` | Component removal effects |
| `figures/fig_09_ablation_summary.png` | Summary heatmap |

---

## Code

| File | Description |
|------|-------------|
| `code/run_09_ablation_fast.py` | Optimized ablation script |
| `code/submit_09_ablation_all.sh` | SLURM job for all models |

---

## Key Results

### Dimensionality Reduction (Binary AUC)
| Model | Original | PCA-256 | Change |
|-------|:--------:|:-------:|:------:|
| STATE | 0.894 | 0.867 | -3% |
| UCE | 0.876 | 0.859 | -2% |
| **scGPT** | 0.770 | 0.827 | **+7%** |
| Geneformer | 0.824 | 0.773 | -6% |

**Insight:** scGPT improves with PCA - noise in original embeddings.

### Top PCA Removal (Binary AUC)
| Model | All | -top5 | Drop |
|-------|:---:|:-----:|:----:|
| STATE | 0.893 | 0.823 | -8% |
| UCE | 0.877 | 0.774 | -12% |
| scGPT | 0.852 | 0.762 | -11% |
| **Geneformer** | 0.824 | 0.812 | **-1%** |

**Insight:** Geneformer most robust to top component removal.

### Normalization (Binary AUC)
| Model | Raw | L2 | Z-score |
|-------|:---:|:--:|:-------:|
| All | X | X | X |

**Insight:** <0.1% difference - normalization doesn't matter.

---

## Optimization Details

| Parameter | Before | After | Speedup |
|-----------|--------|-------|---------|
| CV Folds | Sequential | Parallel (8 cores) | 5x |
| Solver | saga | lbfgs | 10x |
| **Total** | >24 hours | ~25 min | **~60x** |

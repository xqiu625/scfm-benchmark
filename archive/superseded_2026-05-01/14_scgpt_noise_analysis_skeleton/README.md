# 14 - scGPT Noise Analysis

**Accession:** 14
**Status:** ⏳ PENDING
**Priority:** Medium

---

## Purpose

Investigate why scGPT embeddings improve +7% after PCA dimensionality reduction.
- Characterize the "noise" that PCA removes
- Understand what information is encoded in removed components
- Provide insights into scGPT's embedding structure

---

## Background

From ablation analysis (09):
```
scGPT: Original 0.770 → PCA-256: 0.827 (+7.4%)
```

This is **anomalous** - typically PCA loses information. The improvement suggests:
1. scGPT embeddings contain structured noise
2. The noise may be task-irrelevant variance from pretraining
3. PCA acts as a denoiser, not just compressor

---

## Planned Analysis

### 1. Principal Component Characterization
- Compute PCA on scGPT embeddings (512 dims)
- Analyze variance explained by each component
- Compare to STATE/UCE/Geneformer PC structure

### 2. Noise Component Analysis
- Identify components that hurt classification (removed by PCA-256)
- Analyze what these components encode:
  - Batch effects?
  - Cell cycle?
  - Technical artifacts?
  - Pretraining task leakage?

### 3. Component-Performance Correlation
- Plot AUC vs number of PCs retained
- Find optimal dimensionality for scGPT
- Compare to other models' curves

### 4. Biological Signal Analysis
- Correlate PCs with known biological variables
- Test if removed PCs correlate with:
  - Total UMI counts
  - Mitochondrial %
  - Ribosomal gene %
  - Cell size proxies

---

## Expected Output

### Results
| File | Description |
|------|-------------|
| `results/pca_variance_explained.csv` | Variance per PC for all models |
| `results/pc_performance_curve.csv` | AUC vs # PCs retained |
| `results/noise_pc_correlations.csv` | Correlation of removed PCs with biological variables |

### Figures
| File | Description |
|------|-------------|
| `figures/fig_14_variance_explained.png` | Scree plot comparison |
| `figures/fig_14_optimal_dims.png` | AUC vs PCs curve |
| `figures/fig_14_noise_characterization.png` | What removed PCs encode |

---

## Key Questions to Answer

1. What is the optimal dimensionality for scGPT? (Answer: ~256)
2. What do the removed 256 dimensions encode?
3. Is this pattern unique to scGPT or shared with other generative models?
4. Can we predict which embeddings will benefit from PCA?

---

## Hypothesis

The "noise" in scGPT likely comes from its **generative pretraining objective**:
- scGPT is trained to reconstruct masked genes
- This encourages encoding information about individual gene values
- For classification, we need cell-level summary statistics
- PCA removes gene-level details, retaining cell-level patterns

This would explain why discriminative models (STATE, UCE) don't show this pattern.

---

## Code Template

```python
from sklearn.decomposition import PCA
import numpy as np

# Load embeddings
X_scgpt = adata_scgpt.obsm['X_scgpt']

# Full PCA
pca = PCA(n_components=512)
X_pca = pca.fit_transform(X_scgpt)

# Variance explained
var_explained = pca.explained_variance_ratio_

# Components removed by PCA-256 (257-512)
removed_components = pca.components_[256:]

# Correlate with biological variables
correlations = {}
for var in ['total_counts', 'pct_mito', 'pct_ribo']:
    correlations[var] = [np.corrcoef(X_pca[:, i], adata.obs[var])[0, 1]
                         for i in range(256, 512)]
```

---

## Expected Insight

This analysis should reveal whether scGPT's "noise" is:
- **Technical**: Batch effects, sequencing artifacts
- **Biological but irrelevant**: Cell cycle, metabolism
- **Pretraining leakage**: Gene-level reconstruction details

Understanding this will inform **best practices for using scGPT embeddings**.

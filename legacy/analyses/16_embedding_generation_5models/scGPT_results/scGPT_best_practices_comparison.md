# scGPT Best Practices vs. Project Implementation Comparison

**Date:** January 19, 2026
**Analysis:** Comparing official scGPT best practices with implementation in Project 302006

---

## 1. Preprocessing Pipeline Comparison

### Best Practice (Official Documentation)

```python
from scgpt.preprocess import Preprocessor

preprocessor = Preprocessor(
    use_key="X",                    # Key for raw data in adata.layers
    filter_gene_by_counts=3,        # Step 1: Filter genes
    filter_cell_by_counts=False,    # Step 2: Filter cells
    normalize_total=1e4,            # Step 3: Depth normalization
    result_normed_key="X_normed",
    log1p=data_is_raw,              # Step 4: Log1p transformation
    result_log1p_key="X_log1p",
    subset_hvg=n_hvg,               # Step 5: Select highly variable genes (1200-1800)
    hvg_flavor="seurat_v3",
    binning=n_bins,                 # Step 6: Value binning (51 bins)
    result_binned_key="X_binned"
)
preprocessor(adata, batch_key="batch")
```

### Original Project Implementation

```python
# From: scripts/foundation_models/scGPT/process_platelet_chunks_scgpt.sh
# Lines 521-538

# Binning only - no HVG selection, no normalization
X_binned = np.zeros_like(X, dtype=int)
for i in range(X.shape[1]):
    gene_expr = X[:, i]
    if np.max(gene_expr) > 0:
        non_zero = gene_expr[gene_expr > 0]
        if len(non_zero) > 50:
            bins = np.quantile(non_zero, np.linspace(0, 1, 51))
            X_binned[:, i] = np.digitize(gene_expr, bins)
        else:
            max_val = np.max(gene_expr)
            X_binned[:, i] = np.clip((gene_expr / max_val * 50).astype(int), 0, 50)
```

### Gap Analysis

| Preprocessing Step | Best Practice | Original Implementation | Status |
|-------------------|---------------|------------------------|--------|
| Filter genes by counts | `filter_gene_by_counts=3` | Not applied | **Missing** |
| Depth normalization | `normalize_total=1e4` | Not applied | **Missing** |
| Log1p transformation | `log1p=True` | Not applied | **Missing** |
| HVG selection | 1200-1800 genes | All genes in vocab (~15K) | **Missing** |
| Value binning | 51 bins via Preprocessor | Manual 51-bin quantile | Implemented |

---

## 2. Model Choice

### Best Practice
> "The whole-human model is recommended for most applications by default."

### Implementation
```bash
MODEL_DIR="$BASE_DIR/models/scGPT_human/scGPT_human"
```

**Status:** Correct - Using whole-human model

---

## 3. Zero-Shot vs Fine-tuned Performance

### Best Practice
> "scGPT and Geneformer in zero-shot configurations perform inconsistently compared to cell embeddings derived from HVG or generated using scVI or Harmony."
>
> "Fine-tuning is essential for optimal performance - don't rely on zero-shot alone"

### Implementation
- Using zero-shot embeddings (no fine-tuning)
- Benchmark runs logistic regression on frozen embeddings

**Status:** Intentional for benchmark design (comparing zero-shot performance across models)

---

## 4. Batch Integration Objectives

### Best Practice
> "For batch integration, turn on ESC, DAR and DSBN objectives specifically to facilitate batch integration."

```python
hyperparameter_defaults = dict(
    GEPC=True,              # Gene expression modelling for cell objective
    ecs_thres=0.8,          # Elastic cell similarity threshold
    dab_weight=1.0,         # Domain adversarial batch correction weight
    DSBN=True               # Domain-specific batch normalization
)
```

### Implementation
- None of these batch correction objectives enabled
- Using raw embeddings from pre-trained model

**Status:** **Missing** - Explains why scGPT embeddings encode batch effects

---

## 5. Processing Method

### Best Practice
- Use `Preprocessor` class for standardized pipeline
- Batch processing with appropriate batch_size (64-512)

### Implementation
```python
# Process cells one by one (memory efficient but slow)
for i in range(n_cells):
    gene_tensor = torch.tensor(gene_ids).unsqueeze(0)
    output = model(gene_tensor, value_tensor, ...)
```

**Status:** Suboptimal - Single-cell processing is very slow

---

## 6. Impact on Benchmark Results

### Observed Issues

1. **Batch effects encoded in embeddings**
   - UMAP visualization shows discrete crescents
   - Each crescent corresponds to a different batch (Data_NO)
   - Batch-shift AUC: 0.524 (near random)

2. **Performance metrics**
   - scGPT raw: 0.771 AUC (6-class severity)
   - scGPT + PCA-256: 0.831 AUC (+7.8% improvement)
   - PCA removes batch-dominant principal components

3. **Differential diagnosis**
   - scGPT raw: 0.893 Macro AUC
   - scGPT + PCA-256: 0.954 Macro AUC (+6.8% improvement)

### Root Causes

| Issue | Root Cause | Evidence |
|-------|-----------|----------|
| Batch effects | No HVG selection, no normalization, no DSBN | UMAP crescents match Data_NO batches |
| Poor generalization | Generative pre-training captures all patterns | 0.524 AUC on batch-shift task |
| PCA improvement | Top PCs capture batch variation | +7% with PCA-256 denoising |

---

## 7. Comparison with Other Models

| Model | Pre-training | Batch Integration | Batch-Shift AUC |
|-------|-------------|-------------------|-----------------|
| STATE | Perturbation prediction | Excellent | 0.789 |
| UCE | Contrastive (cross-species) | Good | 0.750 |
| Geneformer | Self-supervised (rank-based) | Good | 0.730 |
| TranscriptFormer | Cross-species | Good | N/A |
| **scGPT** | **Generative (masked expression)** | **Poor** | **0.524** |

**Key insight:** scGPT is the only model using generative pre-training, which learns to reconstruct ALL patterns including technical artifacts. Other models use objectives that inherently filter out batch effects.

---

## 8. Recommendations

### For Fair Comparison (Best Practices Implementation)

```python
from scgpt.preprocess import Preprocessor

preprocessor = Preprocessor(
    use_key="X",
    filter_gene_by_counts=3,
    normalize_total=1e4,
    log1p=True,
    subset_hvg=1200,
    hvg_flavor="seurat_v3",
    binning=51,
    result_binned_key="X_binned"
)
preprocessor(adata, batch_key="Data_NO")
```

### Expected Improvements
1. HVG selection reduces noise from non-informative genes
2. Proper normalization ensures consistent scaling across batches
3. Log transformation reduces skewness in expression values
4. Batch-aware preprocessing may reduce batch effects in embeddings

---

## 9. Files

- `scGPT_best_practices_comparison.md` - This document
- `../code/run_scgpt_best_practices.py` - New implementation following best practices
- `../code/submit_scgpt_best_practices.sh` - SLURM submission script

---

**Conclusion:** The original scGPT implementation missed several critical preprocessing steps recommended in official documentation. A new run with best practices will determine whether proper preprocessing can improve scGPT's batch integration and overall performance.

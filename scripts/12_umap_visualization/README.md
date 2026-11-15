# 12 - UMAP Visualization

**Accession:** 12
**Status:** ✅ COMPLETE
**Last Run:** January 5, 2026

---

## Purpose

Generate UMAP visualizations for all 5 foundation model embeddings to visualize:
- Disease state separation
- Severity gradients
- Model-specific embedding structure

---

## Input

| File | Description |
|------|-------------|
| `embeddings/*_embeddings.h5ad` | All 5 model embeddings |

---

## Output

### Figures
| File | Description |
|------|-------------|
| `figures/fig_12_umap_5models_6class.png` | All 5 models, 6-class coloring |
| `figures/fig_12_umap_5models_binary.png` | All 5 models, binary coloring |
| `figures/fig_12_state_severity_gradient.png` | STATE with severity gradient |
| `figures/fig_12_uce_severity_gradient.png` | UCE with severity gradient |

---

## Code

| File | Description |
|------|-------------|
| `code/run_12_umap_visualization.py` | UMAP generation script |

---

## Visualization Gallery

### 5-Model Comparison (6-Class)
![UMAP 5 Models](figures/fig_12_umap_5models_6class.png)

Shows embedding structure for:
- STATE (2,058 dims)
- UCE (1,280 dims)
- Geneformer (1,152 dims)
- TranscriptFormer (2,048 dims)
- scGPT (512 dims)

### Severity Gradient (STATE)
![STATE Severity](figures/fig_12_state_severity_gradient.png)

Shows continuous severity gradient in STATE embeddings.

---

## Key Observations

1. **STATE & UCE**: Clear disease state separation
2. **Geneformer**: Distinct clusters by cell state
3. **scGPT**: More diffuse embedding structure
4. **TranscriptFormer**: Moderate separation

---

## UMAP Parameters

```python
umap_params = {
    'n_neighbors': 15,
    'min_dist': 0.1,
    'metric': 'cosine',
    'random_state': 42
}
```

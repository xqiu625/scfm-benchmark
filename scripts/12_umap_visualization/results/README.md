# 12_umap_visualization Results

## UMAP Embedding Analysis for 5 Foundation Models

**Models:** STATE, UCE, Geneformer, TranscriptFormer, scGPT
**Data:** Platelet transcriptomes (~47,000 cells)
**Date:** January 18, 2026

---

## Key Finding: scGPT Encodes Batch Effects

### Batch Integration Comparison

| Model | UMAP Pattern | Batch Integration | Assessment |
|-------|--------------|-------------------|------------|
| **STATE** | Continuous, mixed | ✅ Excellent | Batches well-integrated |
| **UCE** | Continuous, mixed | ✅ Good | Minor batch clustering |
| **Geneformer** | Dispersed, mixed | ✅ Good | Batches well-integrated |
| **TranscriptFormer** | Continuous, mixed | ✅ Good | Batches integrated |
| **scGPT** | **Discrete crescents** | ❌ Poor | **Each crescent = 1 batch** |

### Evidence

The `umap_5models_by_batch.png` figure shows:
- STATE, UCE, Geneformer, TranscriptFormer: Colors (batches) are mixed throughout the embedding space
- scGPT: Each crescent-shaped cluster corresponds to a single batch (Data_NO)

---

## Implications

### 1. scGPT's Crescent Pattern Explained
- **NOT caused by data contamination** (verified with cleaned platelet-only data)
- **Caused by batch effects** encoded in scGPT embeddings
- Each discrete cluster represents a different data source

### 2. Explains Poor Generalization
| Model | Batch-Shift AUC (COVID→Sepsis) |
|-------|-------------------------------|
| STATE | 0.789 |
| UCE | 0.750 |
| Geneformer | 0.730 |
| scGPT | **0.524** (near random) |

scGPT fails on batch shift because its embeddings are batch-specific, not biologically meaningful.

### 3. Explains PCA Improvement
- scGPT +7% with PCA-256
- PCA likely removes batch-dominant principal components
- Other models don't benefit from PCA (they don't encode batch effects)

---

## Figures Generated

| Figure | Description |
|--------|-------------|
| `umap_5models_6class.png` | 6-class severity comparison |
| `umap_5models_binary.png` | Binary (severe vs non-severe) |
| `umap_5models_by_disease.png` | Colored by disease |
| `umap_5models_by_batch.png` | Colored by Data_NO (12 batches) |
| `scgpt_detailed_analysis.png` | scGPT detailed analysis |

---

## Conclusion

**scGPT's generative pre-training does not effectively remove technical variation from learned representations.** This fundamental difference in embedding structure explains its poor performance on:
- Batch-shift generalization (0.524 AUC)
- Overall severity prediction (0.735 AUC, lowest among models)

Other models (STATE, UCE, Geneformer, TranscriptFormer) produce batch-invariant embeddings that better capture biological signal.

---

**Date:** January 18, 2026
**Job IDs:** 22282372, 22282658

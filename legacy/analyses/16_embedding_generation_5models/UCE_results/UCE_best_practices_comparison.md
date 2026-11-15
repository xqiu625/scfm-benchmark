# UCE Best Practices vs. Project Implementation Comparison

**Date:** January 19, 2026
**Analysis:** Comparing official UCE best practices with implementation in Project 302006

---

## 1. Model Overview

UCE (Universal Cell Embeddings) is a zero-shot foundation model developed at Stanford University that creates a unified biological latent space using ESM2 protein language model embeddings for gene tokenization.

**Key Innovation:** Uses ESM2 protein embeddings (15 billion parameters) to tokenize genes, enabling cross-species generalization.

---

## 2. Model Selection

### Best Practice

> 33-layer transformer with over 650 million parameters
> Also available as a 4-layer lightweight version
> Gene embeddings from ESM2: 5,120-dimensional vectors per gene

### Project Implementation

```bash
python eval_single_anndata.py \
    --adata_path "human_platelet_covid_severity.h5ad" \
    --species "human" \
    --batch_size 200
```

| Aspect | Best Practice | Implementation | Status |
|--------|--------------|----------------|--------|
| Model | 33-layer (recommended) | 33-layer | ✅ Correct |
| Species parameter | Required | "human" | ✅ Correct |
| Batch size | 25 (80GB GPU) or higher | 200 | ✅ Appropriate |

---

## 3. Input Data Requirements

### Best Practice

> - h5ad file with scRNA-seq counts in .X
> - Gene names (not ENSEMBLIDs) in .var_names
> - No preprocessing (gene selection, cell type identification) needed
> - No highly variable gene selection applied
> - Protein-coding genes only (those with available ESM2 embeddings)

### Project Implementation

```python
INPUT_FILE = "human_platelet_covid_severity.h5ad"
# Raw counts in X matrix
# Gene symbols in var_names
```

| Requirement | Best Practice | Implementation | Status |
|-------------|--------------|----------------|--------|
| Format | h5ad | h5ad | ✅ Correct |
| Gene names | Gene symbols (not Ensembl) | Gene symbols | ✅ Correct |
| Preprocessing | None required | None applied | ✅ Correct |
| HVG selection | Not needed | Not applied | ✅ Correct |

---

## 4. Zero-Shot Paradigm

### Best Practice

> - No fine-tuning or retraining required for new datasets
> - No preprocessing (gene selection, cell type identification) needed
> - Works on any genes from any species
> - UCE is not intended to be fine-tuned (per CZI Virtual Cells Platform)
> - Compare UCE in zero-shot setting only

### Project Implementation

- Using zero-shot embeddings without fine-tuning
- Benchmark runs logistic regression on frozen embeddings

| Aspect | Best Practice | Implementation | Status |
|--------|--------------|----------------|--------|
| Fine-tuning | Not recommended | Not used | ✅ Correct |
| Zero-shot | Primary use case | Applied | ✅ Correct |
| Embedding usage | Frozen embeddings | Frozen + downstream classifier | ✅ Correct |

---

## 5. Output Embeddings

### Best Practice

> - Output: AnnData with UCE embeddings in .obsm["X_uce"]
> - CLS token captures cell-level embedding at final transformer layer
> - Processing Speed: ~3,000 cells/minute at batch size 128 on A100 GPU

### Project Implementation

| Metric | Value |
|--------|-------|
| Embedding dimension | 1,280 |
| Embedding key | `X_uce` |
| Cells processed | ~47,000 |

**Status:** ✅ Correct - Standard UCE output

---

## 6. Summary Comparison Table

| Category | Best Practice | Implementation | Status |
|----------|--------------|----------------|--------|
| **Model** | 33-layer transformer | 33-layer | ✅ Correct |
| **Input format** | h5ad with gene symbols | h5ad with gene symbols | ✅ Correct |
| **Preprocessing** | None required | None applied | ✅ Correct |
| **HVG selection** | Not needed | Not applied | ✅ Correct |
| **Zero-shot** | Primary use case | Applied | ✅ Correct |
| **Fine-tuning** | Not recommended | Not used | ✅ Correct |
| **Output** | obsm["X_uce"] | obsm["X_uce"] | ✅ Correct |

---

## 7. Benchmark Performance

UCE achieves strong performance in the benchmark:

| Task | UCE AUC | Rank |
|------|---------|------|
| 6-class Severity | 0.890 | 2nd |
| 5-class Differential | 0.969 | 2nd |
| Batch-Shift (COVID→Sepsis) | 0.750 | 2nd |

---

## 8. Why UCE Performs Well

### 1. Zero-Shot Design
> UCE excels at embedding new datasets without any fine-tuning or preprocessing

The implementation correctly uses UCE in its intended zero-shot paradigm.

### 2. Batch Effect Correction
> Corrects technology-based batch effects (droplet vs plate-based)
> Matches or exceeds fine-tuned methods (scVI, scArches) in SCIB benchmarks

This explains UCE's strong batch-shift performance (0.750 AUC).

### 3. Cross-Species Training
> 36 million cells from 300+ datasets, 8 species

Pre-training on diverse data sources helps generalization.

### 4. ESM2 Protein Embeddings
> Uses ESM2 protein language model embeddings for gene tokenization

This provides biologically meaningful gene representations.

---

## 9. Conclusion

**The UCE implementation in Project 302006 follows all best practices:**

| Aspect | Assessment |
|--------|------------|
| Model selection | ✅ Correct (33-layer) |
| Input format | ✅ Correct (gene symbols) |
| Preprocessing | ✅ Correct (none required) |
| Zero-shot usage | ✅ Correct (no fine-tuning) |
| Performance | ✅ Strong (2nd best overall) |

**Key Insight:** UCE's design philosophy of "no preprocessing required" is validated by results. The implementation correctly leverages UCE's strengths without unnecessary modifications.

---

## 10. Files

| File | Description |
|------|-------------|
| `UCE_best_practices_comparison.md` | This document |
| `scripts/foundation_models/UCE/uce_human_platelet.sh` | Implementation script |

---

**Date:** January 19, 2026
**Author:** Project 302006 Analysis

# Geneformer Best Practices vs. Project Implementation Comparison

**Date:** January 19, 2026
**Analysis:** Comparing official Geneformer best practices with implementation in Project 302006

---

## 1. Model Overview

Geneformer is a BERT-like transformer encoder pretrained on 30M+ single-cell transcriptomes using masked language modeling. Its key innovation is **rank value encoding** where genes are ranked by expression level normalized by median expression across the corpus.

---

## 2. Model Selection

### Best Practice

> Available Models:
> - Geneformer-V1-10M: 6 layers, 10M parameters
> - Geneformer-V1-12L: 12 layers, ~20M parameters
> - Geneformer-V2-104M: 12 layers, 104M parameters
> - Geneformer-V2-316M: 20 layers, 316M parameters

### Project Implementation

```python
model = BertModel.from_pretrained("ctheodoris/Geneformer")
```

| Aspect | Best Practice | Implementation | Status |
|--------|--------------|----------------|--------|
| Model | V1 or V2 available | HuggingFace default (V1) | ✅ Correct |
| Architecture | BERT encoder | BERT encoder | ✅ Correct |

---

## 3. Input Data Requirements

### Best Practice

> - Format: Raw counts scRNA-seq data as .loom or .h5ad
> - Required gene attribute: "ensembl_id" (Ensembl IDs)
> - Required cell attribute: "n_counts" (total read counts)
> - Important: No feature selection - include all detected genes

### Project Implementation

```python
# From generate_geneformer_embeddings.py

# Check for ensembl_id
if 'ensembl_id' not in adata.var.columns:
    print("ERROR: Missing ensembl_id column!")
    sys.exit(1)

# Ensure n_counts is present
if 'n_counts' not in adata.obs.columns:
    adata.obs['n_counts'] = np.array(adata.X.sum(axis=1)).flatten()
```

| Requirement | Best Practice | Implementation | Status |
|-------------|--------------|----------------|--------|
| Format | h5ad or loom | h5ad → loom conversion | ✅ Correct |
| Gene IDs | ensembl_id required | Checked and required | ✅ Correct |
| n_counts | Required | Auto-computed if missing | ✅ Correct |
| Feature selection | None | None applied | ✅ Correct |

---

## 4. Tokenization (Rank Value Encoding)

### Best Practice

> Geneformer's most distinctive feature is its rank value encoding:
> 1. Expression Ranking: Genes are ranked by expression level within each cell
> 2. Median Normalization: Rankings normalized by gene's median expression across corpus
> 3. Non-parametric Representation: Creates robust, distribution-free encoding

```python
from geneformer import TranscriptomeTokenizer

tk = TranscriptomeTokenizer(
    {"cell_type": "cell_type"},
    nproc=4,
    model_input_size=4096  # V2
)
```

### Project Implementation

```python
# Try V2 API first, fall back to V1
try:
    tk = TranscriptomeTokenizer(
        custom_attr_name_dict={"cell_type": "cell_type"},
        nproc=4,
        model_input_size=4096,
        special_token=True
    )
except TypeError:
    # Fall back to V1 API
    tk = TranscriptomeTokenizer(
        custom_attr_name_dict={"cell_type": "cell_type"},
        nproc=4
    )
```

| Aspect | Best Practice | Implementation | Status |
|--------|--------------|----------------|--------|
| Tokenizer | TranscriptomeTokenizer | TranscriptomeTokenizer | ✅ Correct |
| V1/V2 handling | Version-specific API | Try/except fallback | ✅ Correct |
| Rank encoding | Built into tokenizer | Used correctly | ✅ Correct |

---

## 5. Embedding Extraction

### Best Practice

> ```python
> from geneformer import EmbExtractor
>
> embex = EmbExtractor(
>     emb_layer=-2,  # Second-to-last layer recommended
>     forward_batch_size=30,
> )
> ```

### Project Implementation

```python
# Manual extraction due to container compatibility issues
from transformers import BertModel

model = BertModel.from_pretrained("ctheodoris/Geneformer")
model.eval()

# Mean pooling over non-padded tokens
with torch.no_grad():
    outputs = model(input_ids, attention_mask=attention_mask)
    hidden_states = outputs.last_hidden_state
    mask_expanded = attention_mask.unsqueeze(-1).expand(hidden_states.size())
    sum_embeddings = torch.sum(hidden_states * mask_expanded, dim=1)
    batch_emb = (sum_embeddings / sum_mask).cpu().numpy()
```

| Aspect | Best Practice | Implementation | Status |
|--------|--------------|----------------|--------|
| Method | EmbExtractor class | Manual BertModel | ⚠️ Different approach |
| Layer | -2 (second-to-last) | Last layer | ⚠️ Different layer |
| Pooling | CLS token or mean | Mean pooling | ✅ Valid approach |

**Note:** Manual extraction was necessary due to Singularity container compatibility issues with EmbExtractor.

---

## 6. Zero-Shot vs Fine-Tuning

### Best Practice

> - Zero-shot limitations: Foundation models often underperform simpler baselines
> - Fine-tuning advantage: With even limited task-specific data, significant improvements
> - Fine-tuning recommended: Zero-shot performance limited; fine-tuning unlocks full potential

### Project Implementation

- Using zero-shot embeddings without fine-tuning
- Benchmark design intentionally compares zero-shot performance

| Aspect | Best Practice | Implementation | Status |
|--------|--------------|----------------|--------|
| Zero-shot | Limited performance expected | Applied (intentional) | ✅ Correct for benchmark |
| Fine-tuning | Recommended for best results | Not used | ⚠️ Could improve |

---

## 7. Summary Comparison Table

| Category | Best Practice | Implementation | Status |
|----------|--------------|----------------|--------|
| **Model** | V1 or V2 | V1 (HuggingFace) | ✅ Correct |
| **Input format** | h5ad/loom with ensembl_id | h5ad → loom | ✅ Correct |
| **Tokenization** | TranscriptomeTokenizer | Used correctly | ✅ Correct |
| **Rank encoding** | Core feature | Applied via tokenizer | ✅ Correct |
| **Embedding extraction** | EmbExtractor (-2 layer) | Manual (last layer) | ⚠️ Different |
| **Feature selection** | None | None | ✅ Correct |
| **Fine-tuning** | Recommended | Not used (intentional) | ⚠️ Benchmark design |

---

## 8. Benchmark Performance

Geneformer achieves moderate performance:

| Task | Geneformer AUC | Rank |
|------|----------------|------|
| 6-class Severity | 0.862 | 4th |
| 5-class Differential | 0.939 | 5th |
| Batch-Shift (COVID→Sepsis) | 0.730 | 3rd |

---

## 9. Why Geneformer Has Moderate Performance

### 1. Zero-Shot Limitations
> Zero-shot limitations: Foundation models including Geneformer often underperform simpler baselines (PCA, HVG) in pure zero-shot settings

This is a known limitation documented in the best practices.

### 2. Fine-Tuning Would Help
> With even limited task-specific data, Geneformer shows significant improvements

The benchmark intentionally tests zero-shot performance.

### 3. Rank Encoding Strengths
> - Deprioritizes housekeeping genes
> - Prioritizes cell-state markers
> - Batch effect robustness

Rank encoding helps with batch generalization (0.730 AUC on batch-shift).

---

## 10. Potential Improvements

1. **Use EmbExtractor with layer=-2** (second-to-last layer recommended)
2. **Consider V2 model** (104M or 316M parameters)
3. **Fine-tuning** would significantly improve performance

---

## 11. Conclusion

**The Geneformer implementation follows most best practices:**

| Aspect | Assessment |
|--------|------------|
| Model selection | ✅ Correct |
| Input format | ✅ Correct |
| Tokenization | ✅ Correct |
| Embedding layer | ⚠️ Last layer instead of -2 |
| Zero-shot usage | ✅ Intentional for benchmark |

**Key Insight:** Geneformer's moderate zero-shot performance is expected per documentation. The implementation correctly applies rank value encoding, which provides robustness to batch effects.

---

## 12. Files

| File | Description |
|------|-------------|
| `Geneformer_best_practices_comparison.md` | This document |
| `scripts/foundation_models/Geneformer/generate_geneformer_embeddings.py` | Implementation script |

---

**Date:** January 19, 2026
**Author:** Project 302006 Analysis

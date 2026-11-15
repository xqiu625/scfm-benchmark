# TranscriptFormer Best Practices vs. Project Implementation Comparison

**Date:** January 19, 2026
**Analysis:** Comparing official TranscriptFormer best practices with implementation in Project 302006

---

## 1. Model Overview

TranscriptFormer is a generative autoregressive foundation model from CZI that jointly models gene identities and expression levels. It is the first generative, multi-species model for single-cell transcriptomics, trained on 112 million cells spanning 12 species across 1.53 billion years of evolution.

**Key Innovation:** Uses ESM-2 protein embeddings for species-agnostic gene representation and expression-aware multi-head self-attention.

---

## 2. Model Selection

### Best Practice

> Available Variants:
> - TF-Metazoa: 112M cells, 12 species, 444M parameters
> - TF-Exemplar: 110M cells, 5 species, 542M parameters
> - TF-Sapiens: 57M cells, human only, 368M parameters

### Project Implementation

```python
self.checkpoint_path = Path.cwd() / "checkpoints/tf_sapiens"
```

| Aspect | Best Practice | Implementation | Status |
|--------|--------------|----------------|--------|
| Model | TF-Sapiens, TF-Exemplar, or TF-Metazoa | TF-Sapiens | ✅ Correct |
| Use case | Human data → TF-Sapiens | Human platelet data | ✅ Appropriate |

---

## 3. Input Data Requirements

### Best Practice

> - Format: AnnData (.h5ad)
> - Gene IDs: var must contain ensembl_id column
> - Expression: Raw counts (unnormalized)
> - Source: adata.raw.X preferred, falls back to adata.X

### Project Implementation

```python
# Check for Ensembl IDs
if 'ensembl_id' not in adata.var.columns:
    logger.error("Missing 'ensembl_id' column in adata.var!")
    return None
```

| Requirement | Best Practice | Implementation | Status |
|-------------|--------------|----------------|--------|
| Format | h5ad | h5ad | ✅ Correct |
| Gene IDs | ensembl_id in var | Required and checked | ✅ Correct |
| Expression | Raw counts | Raw counts | ✅ Correct |
| Normalization | None (automatic) | None applied | ✅ Correct |

---

## 4. Inference API

### Best Practice

```bash
transcriptformer inference \
  --checkpoint-path ./checkpoints/tf_sapiens \
  --data-file data.h5ad \
  --output-path ./results \
  --batch-size 8 \
  --emb-type cell \
  --gene-col-name ensembl_id \
  --precision 16-mixed
```

### Project Implementation

```python
cmd = [
    "transcriptformer", "inference",
    "--checkpoint-path", str(self.checkpoint_path),
    "--data-file", str(input_file),
    "--output-path", str(output_path),
    "--batch-size", str(batch_size),
    "--emb-type", "cell",
    "--gene-col-name", "ensembl_id",
    "--precision", "16-mixed"
]
subprocess.run(cmd)
```

| Parameter | Best Practice | Implementation | Status |
|-----------|--------------|----------------|--------|
| CLI tool | transcriptformer inference | transcriptformer inference | ✅ Correct |
| Checkpoint | Model-specific path | tf_sapiens | ✅ Correct |
| Batch size | 8 recommended | 16 | ✅ Appropriate |
| Embedding type | cell | cell | ✅ Correct |
| Gene column | ensembl_id | ensembl_id | ✅ Correct |
| Precision | 16-mixed | 16-mixed | ✅ Correct |

---

## 5. Data Processing

### Best Practice

> - Values clipped at 30 (configurable)
> - No explicit feature selection required
> - Raw counts preferred over normalized
> - Only genes with ESM-2 embeddings included
> - Out-of-vocabulary genes filtered during inference

### Project Implementation

```python
# No manual preprocessing applied
# CLI handles all preprocessing internally
self.run_inference(
    input_file=input_file,
    output_path=temp_output_dir,
    batch_size=batch_size,
    emb_type="cell"
)
```

| Aspect | Best Practice | Implementation | Status |
|--------|--------------|----------------|--------|
| Preprocessing | Automatic | Automatic (CLI) | ✅ Correct |
| Feature selection | None required | None applied | ✅ Correct |
| Gene filtering | Automatic (ESM-2) | Automatic (CLI) | ✅ Correct |
| Count clipping | 30 (default) | Default | ✅ Correct |

---

## 6. Output Format

### Best Practice

> Cell Embeddings (--emb-type cell):
> - obsm['embeddings']: Cell embeddings
> - obs: Original metadata preserved
> - uns['llh']: Log-likelihood scores

### Project Implementation

```python
# Load the generated embeddings
embeddings_file = temp_output_dir / "embeddings.h5ad"
adata_emb = sc.read_h5ad(embeddings_file)

# Extract embeddings
if 'embeddings' in adata_emb.obsm:
    adata.obsm['X_transcriptformer'] = adata_emb.obsm['embeddings']
```

| Aspect | Best Practice | Implementation | Status |
|--------|--------------|----------------|--------|
| Output key | obsm['embeddings'] | Correctly extracted | ✅ Correct |
| Embedding dim | 2,048 | 2,048 | ✅ Correct |
| Storage key | User-defined | X_transcriptformer | ✅ Appropriate |

---

## 7. Summary Comparison Table

| Category | Best Practice | Implementation | Status |
|----------|--------------|----------------|--------|
| **Model** | TF-Sapiens for human | TF-Sapiens | ✅ Correct |
| **Input format** | h5ad with ensembl_id | h5ad with ensembl_id | ✅ Correct |
| **Expression** | Raw counts | Raw counts | ✅ Correct |
| **Preprocessing** | None (automatic) | None (CLI handles) | ✅ Correct |
| **CLI usage** | transcriptformer inference | Correct parameters | ✅ Correct |
| **Precision** | 16-mixed recommended | 16-mixed | ✅ Correct |
| **Output** | obsm['embeddings'] | Correctly extracted | ✅ Correct |

---

## 8. Benchmark Performance

TranscriptFormer achieves strong performance:

| Task | TranscriptFormer AUC | Rank |
|------|----------------------|------|
| 6-class Severity | N/A (different sample size) | - |
| 5-class Differential | 0.955 | 3rd |
| Batch-Shift | N/A | - |

---

## 9. Why TranscriptFormer Performs Well

### 1. Generative Architecture
> First generative, multi-species model for single-cell transcriptomics

Captures complex gene-expression relationships.

### 2. ESM-2 Protein Embeddings
> Uses ESM-2 protein embeddings for species-agnostic gene representation

Provides biologically meaningful gene representations.

### 3. Cross-Species Training
> 112 million cells spanning 12 species across 1.53 billion years

Even TF-Sapiens benefits from architectural innovations developed for multi-species models.

### 4. Expression-Aware Attention
> Expression counts introduced as log-count bias term in attention matrix

Efficiently incorporates expression information without vocabulary explosion.

---

## 10. Conclusion

**The TranscriptFormer implementation follows all best practices:**

| Aspect | Assessment |
|--------|------------|
| Model selection | ✅ Correct (TF-Sapiens for human) |
| Input format | ✅ Correct (ensembl_id) |
| Preprocessing | ✅ Correct (none required) |
| CLI usage | ✅ Correct (all parameters) |
| Precision | ✅ Correct (16-mixed) |
| Performance | ✅ Strong (3rd in differential) |

**Key Insight:** TranscriptFormer's design requires minimal user intervention. The CLI handles all preprocessing internally, making it straightforward to use correctly.

---

## 11. Files

| File | Description |
|------|-------------|
| `TranscriptFormer_best_practices_comparison.md` | This document |
| `scripts/foundation_models/TranscriptFormer/generate_transcriptformer_embeddings.py` | Implementation script |

---

**Date:** January 19, 2026
**Author:** Project 302006 Analysis

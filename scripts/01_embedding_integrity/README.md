# 01 - Embedding Integrity Check

**Accession:** 01
**Status:** ✅ COMPLETE
**Last Run:** January 2, 2026

---

## Purpose

Validate all foundation model embeddings before analysis:
- Check for NaN/Inf values
- Verify cell counts match metadata
- Confirm embedding dimensions
- Detect synthetic or random data

---

## Input

| File | Path | Description |
|------|------|-------------|
| STATE | `embeddings/state_embeddings.h5ad` | 47,808 cells × 2,058 dims |
| UCE | `embeddings/uce_embeddings.h5ad` | 47,808 cells × 1,280 dims |
| TranscriptFormer | `embeddings/transcriptformer_embeddings.h5ad` | 46,949 cells × 2,048 dims |
| scGPT | `embeddings/platelet_scgpt_combined.h5ad` | 47,808 cells × 512 dims |
| Geneformer | `embeddings/geneformer_embeddings.h5ad` | 46,949 cells × 1,152 dims |

---

## Output

### Results
| File | Description |
|------|-------------|
| `results/embedding_integrity_report.csv` | Validation results for all models |

### Figures
None (validation only)

---

## Code

| File | Description |
|------|-------------|
| `code/run_01_embedding_integrity.py` | Main validation script |

---

## Key Results

| Model | Status | Cells | Dims | Issues |
|-------|--------|-------|------|--------|
| STATE | ✅ VALID | 47,808 | 2,058 | None |
| UCE | ✅ VALID | 47,808 | 1,280 | None |
| TranscriptFormer | ✅ VALID | 46,949 | 2,048 | Re-generated Jan 2 |
| scGPT | ✅ VALID | 47,808 | 512 | None |
| Geneformer | ✅ VALID | 46,949 | 1,152 | Minor duplicates |

**All embeddings verified clean.**

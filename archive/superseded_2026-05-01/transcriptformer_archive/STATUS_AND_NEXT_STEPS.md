# TranscriptFormer Status & Next Steps

**Date:** November 26, 2025
**Project:** 302006 - COVID-19 Disease Severity Prediction

---

## Current Status

### ✅ Completed Successfully

1. **COVID Dataset Embeddings**
   - Status: ✅ Complete
   - Cells: 47,109
   - Embedding dimensions: 2,048
   - Processing time: ~15 minutes
   - Location: `/bigdata/godziklab/shared/Xinru/302006/02_EMBEDDINGS/transcriptformer/covid_transcriptformer.h5ad`

### ⚠️ Requires Fix

2. **HEK293T Dataset**
   - Status: ❌ Failed - Duplicate Ensembl ID issue
   - Error: `ValueError: Mismatch between expression matrix columns (24143) and gene names length (24142)`
   - Cause: 1 duplicate Ensembl ID in prepared data
   - Fix ready: Yes ✅

3. **HCT116 Dataset**
   - Status: ❌ Failed - Duplicate Ensembl ID issue
   - Error: `ValueError: Mismatch between expression matrix columns (24143) and gene names length (24142)`
   - Cause: 1 duplicate Ensembl ID in prepared data
   - Fix ready: Yes ✅

---

## The Problem

### What Happened

TranscriptFormer's `--remove-duplicate-genes` flag removes duplicate gene names but NOT the corresponding matrix columns, causing a dimension mismatch:

```
Prepared data:
- Expression matrix: 24,143 columns
- Ensembl IDs: 24,143 (with 1 duplicate)

TranscriptFormer loads and processes:
- Removes duplicate from gene list → 24,142 genes
- Matrix still has 24,143 columns
- ERROR: Dimension mismatch!
```

### Why COVID Worked

COVID dataset had **no duplicate Ensembl IDs** (12,729 unique genes), so it processed perfectly.

HEK293T and HCT116 each had **1 duplicate Ensembl ID** (out of 24,143 genes), triggering the error.

---

## The Solution

### What We Fixed

Updated `prepare_data_for_transcriptformer.py` to remove duplicate Ensembl IDs BEFORE saving:

```python
# Check for and remove duplicate Ensembl IDs
ensembl_ids = adata.var['ensembl_id'].values
unique_ids, unique_indices = np.unique(ensembl_ids, return_index=True)

n_duplicates = len(ensembl_ids) - len(unique_ids)
if n_duplicates > 0:
    # Remove BOTH the duplicate gene AND its matrix column
    adata = adata[:, unique_indices].copy()
```

This ensures perfect alignment between expression matrix and gene names.

---

## Next Steps to Complete

### Step 1: Fix HEK293T and HCT116 Data

Run the fix script on the cluster:

```bash
cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit
conda activate transcriptformer
bash scripts/foundation_models/TranscriptFormer/fix_and_regenerate_cell_lines.sh
```

**What it does:**
1. Reloads HEK293T and HCT116 datasets
2. Removes duplicate Ensembl IDs (1 per dataset)
3. Resaves cleaned datasets
4. Generates TranscriptFormer embeddings

**Expected time:** ~30-45 minutes total
- Data fixing: ~2 minutes
- HEK293T embeddings: ~5-7 minutes (13K cells)
- HCT116 embeddings: ~20-25 minutes (40K cells)

**Expected output:**
```
HEK293T: (13,166 cells × 24,142 genes) → 2048-dim embeddings
HCT116: (40,044 cells × 24,142 genes) → 2048-dim embeddings
```

### Step 2: Verify All Embeddings

After the fix script completes, verify all three datasets:

```python
import scanpy as sc

# Load and check each dataset
datasets = ['covid', 'hek293t', 'hct116']
base_path = '/bigdata/godziklab/shared/Xinru/302006/02_EMBEDDINGS/transcriptformer'

for dataset in datasets:
    adata = sc.read_h5ad(f'{base_path}/{dataset}_transcriptformer.h5ad')
    print(f"\n{dataset.upper()}:")
    print(f"  Cells: {adata.shape[0]:,}")
    print(f"  Genes: {adata.shape[1]:,}")
    print(f"  Embeddings: {adata.obsm['X_transcriptformer'].shape}")
    print(f"  Keys: {list(adata.obsm.keys())}")
```

**Expected output:**
```
COVID:
  Cells: 47,109
  Genes: 12,729
  Embeddings: (47109, 2048)
  Keys: ['X_transcriptformer']

HEK293T:
  Cells: 13,166
  Genes: 24,142
  Embeddings: (13166, 2048)
  Keys: ['X_transcriptformer']

HCT116:
  Cells: 40,044
  Genes: 24,142
  Embeddings: (40044, 2048)
  Keys: ['X_transcriptformer']
```

---

## All Challenges Encountered

We encountered and solved **7 major challenges** during TranscriptFormer implementation:

1. ✅ **Python Version Compatibility** - Python 3.6 → 3.11 via conda
2. ✅ **Package Dependency Conflicts** - zarr 3.x → 2.18.3, PyTorch ≤2.5.1
3. ✅ **TranscriptFormer API Misconception** - CLI-based, not Python API
4. ✅ **Ensembl ID Mapping Issues** - Duplicate mappings, length mismatches
5. ✅ **Gene Filtering for None Values** - Matrix dimension problems
6. ✅ **Missing Metadata Requirements** - Added assay/organism/tissue columns
7. ✅ **Duplicate Ensembl IDs After Preparation** - Fixed in data prep (current issue)

All solutions documented in: `scripts/foundation_models/TranscriptFormer/TROUBLESHOOTING.md`

---

## Files Created

### Scripts
1. `install_transcriptformer_conda.sh` - Installation with Python 3.11
2. `prepare_data_for_transcriptformer.py` - Data preparation (updated with duplicate removal)
3. `generate_transcriptformer_embeddings.py` - Embedding generation via CLI
4. `run_transcriptformer_pipeline.sh` - Full pipeline orchestration
5. `fix_and_regenerate_cell_lines.sh` - **Quick fix for HEK293T/HCT116** ← Run this!

### Documentation
1. `TROUBLESHOOTING.md` - Complete troubleshooting guide (7 challenges)
2. `STATUS_AND_NEXT_STEPS.md` - This file
3. `DAILY_PROGRESS.md` - Session summary

---

## Summary

### What Works ✅
- COVID dataset: **47,109 cells × 2048 embeddings** - Complete!
- Installation and environment setup
- Data preparation pipeline
- Embedding generation pipeline
- All 7 challenges documented and solved

### What Needs Fixing ⚠️
- HEK293T: Remove 1 duplicate Ensembl ID, regenerate embeddings (~7 min)
- HCT116: Remove 1 duplicate Ensembl ID, regenerate embeddings (~25 min)

### Quick Command
```bash
# Just run this one script to fix everything:
bash scripts/foundation_models/TranscriptFormer/fix_and_regenerate_cell_lines.sh
```

---

## After All Embeddings Complete

Once all three datasets have embeddings, proceed to:

### 1. Three-Model Comparison
Compare TranscriptFormer vs UCE vs scGPT:
- UMAP visualizations
- Clustering performance (ARI, NMI)
- COVID severity prediction accuracy
- Batch effect analysis

### 2. Model Performance Benchmarking
- Embedding generation time
- GPU memory usage
- Embedding quality metrics
- Downstream task performance

### 3. Final Analysis
- Which model best captures COVID severity signal?
- Ensemble strategies
- Publication-quality figures

---

## Quick Reference

### File Locations

**Data:**
```
/bigdata/godziklab/shared/Xinru/302006/01_DATA/transcriptformer_ready/
├── covid_transcriptformer_ready.h5ad       (47K cells × 12.7K genes)
├── hek293t_transcriptformer_ready.h5ad     (13K cells × 24.1K genes) ← Needs fix
└── hct116_transcriptformer_ready.h5ad      (40K cells × 24.1K genes) ← Needs fix
```

**Embeddings:**
```
/bigdata/godziklab/shared/Xinru/302006/02_EMBEDDINGS/transcriptformer/
├── covid_transcriptformer.h5ad             ✅ Complete (47K × 2048)
├── hek293t_transcriptformer.h5ad           ⚠️  Pending (will be 13K × 2048)
└── hct116_transcriptformer.h5ad            ⚠️  Pending (will be 40K × 2048)
```

### Key Parameters
- Model: TF-Sapiens (368M parameters, 57M human cells)
- Embedding dimensions: 2,048
- Batch size: 16
- Precision: fp16 (mixed precision)
- Processing speed: ~3.3 iterations/second

---

**Status:** 1/3 datasets complete, 2/3 require simple fix
**Time to complete:** ~35 minutes total
**Next action:** Run `fix_and_regenerate_cell_lines.sh`

---

*Last updated: November 26, 2025*

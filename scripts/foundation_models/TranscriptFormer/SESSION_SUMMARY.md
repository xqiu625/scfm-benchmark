# TranscriptFormer Implementation - Session Summary

**Date:** November 26, 2025
**Duration:** 5-6 hours
**Goal:** Generate TranscriptFormer embeddings for COVID, HEK293T, and HCT116 datasets

---

## TL;DR

✅ **COVID dataset successful:** 47,109 cells × 2,048 dimensions
❌ **HEK293T/HCT116 blocked:** Gene ID validation issues
📚 **8 challenges documented** in TROUBLESHOOTING.md (7 solved, 1 ongoing)

---

## What Worked

### COVID Dataset ✅
- **Status:** Complete
- **Output:** `/bigdata/godziklab/shared/Xinru/302006/02_EMBEDDINGS/transcriptformer/covid_transcriptformer.h5ad`
- **Specs:** 47,109 cells × 2,048 embedding dimensions
- **Time:** ~15 minutes processing
- **Why it worked:** Clean Ensembl IDs, no mixed gene ID types

### Installation Pipeline ✅
- Conda environment with Python 3.11
- Exact package version pinning
- TF-Sapiens model download
- **Script:** `install_transcriptformer_conda.sh`

### Documentation ✅
- Comprehensive troubleshooting guide (8 challenges)
- Daily progress log with recommendations
- Diagnostic tools for future debugging

---

## What's Blocked

### HEK293T & HCT116 Datasets ❌

**Issue:** TranscriptFormer reports "Found duplicate genes" despite our diagnostics showing no duplicates

**Root cause identified:**
- mygene.info returned mixed gene ID types
- 2 genes per dataset mapped to Entrez IDs (plain numbers) instead of Ensembl IDs
- Example: `126231` (Entrez) vs `ENSG00000000419` (Ensembl)

**Attempted solutions:**
1. Duplicate removal → Still failed
2. Version number stripping → Still failed
3. Entrez ID filtering → Prepared but not tested

**Mystery:** Our diagnostic script shows "no duplicates" but TranscriptFormer's internal validation still fails, suggesting opaque/undocumented validation logic

---

## Challenges Solved (7/8)

1. ✅ **Python Version** - Conda with Python 3.11
2. ✅ **Package Conflicts** - zarr 2.18.3, PyTorch ≤2.5.1
3. ✅ **API Misconception** - CLI-based, not Python API
4. ✅ **Ensembl Mapping** - mygene.info batching
5. ✅ **None Filtering** - Explicit None value removal
6. ✅ **Metadata Requirements** - assay/organism/tissue columns
7. ✅ **Duplicate IDs** - Version stripping + removal
8. ⚠️ **Mixed ID Types** - Identified, fix prepared, not tested

Full details: `TROUBLESHOOTING.md`

---

## Files Created

### Working Scripts (Keep)
- `install_transcriptformer_conda.sh` - Installation
- `prepare_data_for_transcriptformer.py` - Data prep (latest)
- `generate_transcriptformer_embeddings.py` - Embedding generation
- `diagnose_duplicates.py` - Diagnostic tool
- `final_fix_entrez_ids.sh` - Ready to test

### Failed Attempts (Archive)
- `install_transcriptformer.sh` - Superseded
- `fix_and_regenerate_cell_lines.sh` - Failed approach
- `regenerate_cell_lines.sh` - Failed approach
- `final_fix_and_run.sh` - Superseded

### Documentation
- `TROUBLESHOOTING.md` - Comprehensive (8 challenges)
- `DAILY_PROGRESS_1126.md` - Detailed session log
- `STATUS_AND_NEXT_STEPS.md` - Current status
- `SESSION_SUMMARY.md` - This file

---

## Recommendations

### Option 1: Continue Debugging 🔍
**Time:** 2-3 hours
**Approach:**
1. Test `final_fix_entrez_ids.sh`
2. Review TranscriptFormer source code
3. Contact developers if needed

**Pros:** Complete all 3 datasets
**Cons:** Time-intensive, uncertain outcome

### Option 2: Proceed with COVID Only ⚡ (RECOMMENDED)
**Time:** 30 minutes
**Approach:**
1. Run three-model comparison with COVID dataset
2. Note HEK293T/HCT116 as limitation
3. Focus on analysis

**Pros:** Move forward, COVID is primary dataset
**Cons:** Missing cell line comparisons

### Option 3: Start Fresh 🔄
**Time:** 1-2 hours
**Approach:**
1. Re-prepare HEK293T/HCT116 from raw data
2. Ensure all genes map to Ensembl only
3. Follow COVID's successful approach

**Pros:** Clean slate, higher success chance
**Cons:** Moderate time investment

---

## Key Learnings

1. **TranscriptFormer is exceptionally strict** about input format
2. **CLI tools can have opaque validation** - may need source review
3. **Gene ID mapping is complex** - multiple ID types exist
4. **Diagnostic tools are essential** - custom inspectors help
5. **Working examples are templates** - COVID shows requirements
6. **Version pinning is critical** - PyTorch ≤2.5.1, zarr 2.18.3
7. **Documentation saves time** - future debugging reference

---

## Next Session Checklist

- [ ] Review recommendations and choose approach
- [ ] If Option 2: Run three-model comparison with COVID
- [ ] If Option 1: Test `final_fix_entrez_ids.sh`
- [ ] If Option 3: Re-prepare from raw data
- [ ] Archive failed attempt scripts
- [ ] Update main project README if needed

---

## Quick Reference

### Working COVID Embeddings
```
Location: /bigdata/godziklab/shared/Xinru/302006/02_EMBEDDINGS/transcriptformer/
File: covid_transcriptformer.h5ad
Shape: (47109, 12729) → 2048 dims
```

### Test Diagnostic
```bash
python scripts/foundation_models/TranscriptFormer/diagnose_duplicates.py
```

### Installation
```bash
bash scripts/foundation_models/TranscriptFormer/install_transcriptformer_conda.sh
```

---

**Status:** Partial success - 1/3 datasets complete
**Recommendation:** Proceed with COVID, return to cell lines if time permits
**Documentation:** Complete and comprehensive

---

_Created: November 26, 2025, 5:00 PM_
_Project: 302006 - COVID-19 Disease Severity Prediction_

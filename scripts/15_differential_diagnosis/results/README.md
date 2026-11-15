# 15_differential_diagnosis Results

## 5-Class Differential Diagnosis Benchmark

**Task:** Distinguish between 5 inflammatory conditions from platelet transcriptomes
**Classes:** COVID-19, Sepsis, SLE, Pneumonia, Healthy
**Excluded:** Influenza (n=33 too small), Non_COVID (catch-all)

---

## Results Summary (All 6 Models)

| Rank | Model | Macro AUC | Weighted AUC | Notes |
|:----:|-------|:---------:|:------------:|-------|
| 1 | **STATE** | **0.978** | **0.962** | Best overall |
| 2 | UCE | 0.969 | 0.947 | |
| 3 | TranscriptFormer | 0.955 | 0.929 | |
| 4 | scGPT_PCA256 | 0.954 | 0.917 | PCA denoising required |
| 5 | Geneformer | 0.939 | 0.901 | |
| 6 | scGPT | 0.893 | 0.836 | Raw embeddings encode batch effects |

**Note:** scGPT improves +6% with PCA-256 due to batch effects encoded in raw embeddings.

---

## Per-Class AUC (One-vs-Rest)

| Model | COVID-19 | Healthy | Pneumonia | SLE | Sepsis |
|-------|:--------:|:-------:|:---------:|:---:|:------:|
| STATE | 0.958 | 0.961 | 0.973 | **0.999** | **0.998** |
| UCE | 0.941 | 0.948 | 0.963 | **0.999** | 0.995 |
| TranscriptFormer | 0.923 | 0.919 | 0.947 | 0.998 | 0.989 |
| scGPT_PCA256 | 0.909 | 0.889 | 0.980 | 0.999 | 0.994 |
| Geneformer | 0.892 | 0.915 | 0.908 | 0.996 | 0.984 |
| scGPT | 0.827 | 0.715 | 0.954 | 0.993 | 0.976 |

---

## Key Findings

1. **STATE achieves best performance** (0.978 Macro AUC)
2. **SLE easiest to classify** (0.993-0.999 for all models) - distinct autoimmune signature
3. **Sepsis well-separated** (0.976-0.998) - strong inflammatory signature
4. **scGPT improves +6% with PCA-256** (0.893 → 0.954) - removes batch effects
5. **scGPT struggles with Healthy** (0.715 AUC) - batch effects confound classification
6. **COVID-19 and Healthy most challenging** - overlapping with other conditions

---

## Sample Sizes

| Class | N Samples |
|-------|-----------|
| COVID-19 | ~38,673 |
| Sepsis | ~3,422 |
| Healthy | ~3,205 |
| Pneumonia | ~1,749 |
| SLE | ~169 |

---

## Files

- `differential_diagnosis_5class_results.csv` - Full results table
- `figures/` - Confusion matrices (pending)

---

**Date:** January 18, 2026
**Job IDs:** 22281132 (main), 22282656 (scGPT rerun)

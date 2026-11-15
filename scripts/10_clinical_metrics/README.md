# 10 - Clinical Metrics

**Accession:** 10
**Status:** ✅ COMPLETE

---

## Purpose

Evaluate clinical utility metrics beyond AUC-ROC:
- AUC-PR (Precision-Recall)
- Sensitivity at 90% Specificity
- Cohen's Kappa

---

## Key Results (Binary Task)

| Model | AUC-ROC | AUC-PR | Sens@90%Spec | Cohen's κ |
|-------|:-------:|:------:|:------------:|:---------:|
| **STATE** | **0.894** | **0.934** | **71.7%** | **0.611** |
| UCE | 0.876 | 0.923 | 67.4% | 0.568 |
| TranscriptFormer | 0.838 | 0.895 | 58.5% | 0.501 |
| Geneformer | 0.824 | 0.888 | 56.4% | 0.475 |
| scGPT | 0.776 | 0.836 | 40.8% | 0.402 |

**STATE has best clinical utility.**

---

## Output

### Figures
| File | Description |
|------|-------------|
| `figures/fig_10_clinical_metrics.png` | Clinical metrics comparison |
| `figures/fig_10_pr_curves.png` | Precision-Recall curves |

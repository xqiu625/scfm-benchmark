# 06 - Batch Generalization

**Accession:** 06
**Status:** ✅ COMPLETE

---

## Purpose

Test generalization under distribution shift: train on COVID-19, evaluate on sepsis.

---

## Key Results

| Rank | Model | Sepsis AUC |
|:----:|-------|:----------:|
| 1 | **STATE** | **0.789** |
| 2 | UCE | 0.750 |
| 3 | Geneformer | 0.730 |
| 4 | scGPT | 0.524 |

**STATE generalizes best. scGPT fails (near random).**

---

## Output

### Figures
| File | Description |
|------|-------------|
| `figures/fig_06_batch_generalization.png` | Transfer performance |

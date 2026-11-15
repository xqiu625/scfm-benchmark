# 07 - Gene Dropout

**Accession:** 07
**Status:** ✅ COMPLETE

---

## Purpose

Evaluate robustness to missing genes by randomly masking 10-70% of genes.

---

## Key Results

| Model | 0% | 70% | Drop |
|-------|:--:|:---:|:----:|
| **UCE** | 0.876 | 0.861 | **-1.7%** |
| STATE | 0.893 | 0.862 | -3.4% |

**UCE is most dropout-robust.**

---

## Output

### Figures
| File | Description |
|------|-------------|
| `figures/fig_07_gene_dropout.png` | Dropout curves |

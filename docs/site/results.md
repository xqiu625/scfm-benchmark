# Results

> **Provisional historical results.** This page predates the 2026-09
> eight-cohort QC rebuild. Do not use these values as current evidence; see
> `../current-status.md`. Regenerate this page after P2/P3 completes.

Headline numbers from the BIB submission (manuscript ~97% complete, May 2026).
All AUCs use patient-level `StratifiedGroupKFold` with `Sample_ID` (413
unique donors in the COVID cohort).

## Primary task: in-hospital mortality (binary)

47,808 platelets from 12 cohorts. `fatal = 1` (n = 6,237 COVID + 2,613 sepsis)
vs `survived = 0` (mild + moderate + recovered + severe; severe = ICU-discharged
survivor).

| Rank | Model | Binary AUC (LR) | LOBO mean AUC |
|---:|---|---:|---:|
| 1 | **STATE** | **0.825** | **0.789** |
| — | PCA-200 baseline | 0.814 | — |
| 2 | UCE | 0.794 | 0.750 |
| 3 | Geneformer | 0.778 | 0.730 |
| 4 | TranscriptFormer | 0.776 | — |
| 5 | scGPT | 0.663 | 0.494 |

**PCA-200 beats 4 of 5 foundation models.** Only STATE beats raw expression
with PCA. All DeLong pairwise tests significant except TF vs Geneformer.

## Robustness profile

| Axis | Best model | Value |
|---|---|---|
| Distribution shift (LOBO) | STATE | 0.789 |
| Gene dropout @ 70% missing | UCE | 0.0% drop |
| Sample efficiency @ 5% labels | UCE | -8.5% |
| FP16 quantisation | (all) | 0% AUC loss, 2× speedup |

## Cross-disease bidirectional transfer

Train COVID-mortality → test sepsis-mortality, and reverse. Same LR
hyperplane applied across diseases — a pure test of embedding-space
geometry.

| Model | COVID → Sepsis | Sepsis → COVID | Pattern |
|---|---:|---:|---|
| **STATE** | **0.800** | **0.706** | Bidirectional generalizer |
| STACK | 0.917 | 0.698 | Bidirectional generalizer |
| UCE | 0.798 | 0.514 | Forward-only |
| scPRINT-2 | 0.855 | 0.543 | Forward-only |
| TF | 0.775 | 0.629 | Moderate bidirectional |
| Geneformer | 0.742 | 0.581 | Weak bidirectional |
| **scGPT** | **0.352** ⚠ | 0.555 | **Anti-correlated** (CI [0.22, 0.50]) |

scGPT systematically inverts mortality prediction sign across diseases,
likely because its 512-dim embedding compresses platelet activation —
hyper-activated in COVID-fatal but exhausted in sepsis-fatal.

## Adaptation curve (the null finding)

`PCA-200+LR` from-scratch vs FMs, varying labelled sepsis cells `n` added
to the COVID training pool. 10 seeds per `n`, patient-stratified holdout.

| Method | `n=0` (zero-shot) | `n=100` | Asymptote |
|---|---:|---:|---:|
| **PCA-200+LR** | 0.891 | **0.936** | **0.964** |
| scPRINT-2 | **0.920** | 0.931 | 0.950 |
| STATE | 0.882 | 0.912 | 0.948 |
| TranscriptFormer | 0.811 | 0.860 | 0.912 |
| Geneformer | 0.804 | 0.833 | 0.878 |
| UCE | 0.646 | 0.756 | 0.876 |
| STACK | 0.643 | 0.887 | 0.890 |

**No FM has both higher Y-intercept AND steeper early slope than the
from-scratch baseline.** Only scPRINT-2 beats PCA on the zero-shot
intercept, but it is already flat (slope +0.010 to `n=100`). PCA's
asymptote beats every FM at full data.

## Multimodal fusion (the silver lining)

Late-averaging predicted probabilities across all 6 FMs beats any single
FM (AUC 0.849 vs STATE's 0.831). Suggests the FMs encode partially
non-redundant information — supported by the SHAP-Jaccard mechanism
diagnostic showing two distinct biological axes.

## What this means for the field

1. **Patient-level CV must be the floor** for any clinical scFM benchmark.
   Cell-level CV reports AUCs that are 7–17 percentage points inflated.
2. **PCA + LR is a required baseline.** Any FM that does not beat PCA on
   the target task does not justify its computational cost.
3. **High AUC is necessary but not sufficient.** Two FMs at the same AUC
   can encode incompatible biology — verify with the
   [Diagnostic Suite](diagnostic-suite.md) before deploying.
4. **Multimodal fusion is currently the best deployment strategy.**
   No single FM dominates across robustness, calibration, and generalisation.

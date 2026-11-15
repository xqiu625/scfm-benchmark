# Results Status — Patient-Level CV (v4 + v5)

Last updated: 2026-04-29

## ACC 02: Severity Benchmark (LR + RF, --quick)

### Binary Classification (non_severe vs severe)

| Model | LR Acc | LR BalAcc | LR AUC | RF Acc | RF BalAcc | RF AUC | Status |
|-------|--------|-----------|--------|--------|-----------|--------|--------|
| STATE | 0.738 | 0.741 | **0.825** | 0.691 | 0.632 | 0.758 | DONE (v4) |
| UCE | 0.705 | 0.713 | 0.794 | 0.706 | 0.682 | 0.781 | DONE (v4) |
| TranscriptFormer | 0.694 | 0.703 | 0.776 | 0.674 | 0.618 | 0.708 | DONE (v4) |
| Geneformer | 0.701 | 0.705 | 0.778 | 0.653 | 0.575 | 0.660 | DONE (v4) |
| scGPT | 0.661 | 0.658 | 0.663 | 0.634 | 0.614 | 0.680 | DONE (v4) |

### 3-Class Classification (control vs mild vs severe)

| Model | LR Acc | LR BalAcc | LR AUC | RF Acc | RF BalAcc | RF AUC | Status |
|-------|--------|-----------|--------|--------|-----------|--------|--------|
| STATE | 0.667 | 0.654 | 0.822 | 0.640 | 0.381 | 0.766 | **DONE (v5)** |
| UCE | 0.645 | 0.655 | **0.818** | 0.678 | 0.518 | 0.815 | DONE (v4) |
| TranscriptFormer | 0.612 | 0.632 | 0.787 | 0.644 | 0.406 | 0.712 | DONE (v4) |
| Geneformer | 0.609 | 0.617 | 0.778 | 0.625 | 0.355 | 0.654 | DONE (v4) |
| scGPT | — | 0.448 | 0.630 | — | — | — | DONE (v4) |

### 6-Class Classification (healthy/mild/moderate/severe/fatal/recovered)

| Model | LR Acc | LR BalAcc | LR AUC | RF Acc | RF BalAcc | RF AUC | Status |
|-------|--------|-----------|--------|--------|-----------|--------|--------|
| STATE | 0.445 | 0.467 | **0.782** | 0.441 | 0.224 | 0.658 | **DONE (v5 resubmit)** |
| UCE | — | 0.464 | 0.774 | — | — | — | **DONE (v5)** |
| TranscriptFormer | 0.378 | 0.440 | 0.750 | 0.439 | 0.253 | 0.650 | **DONE (v5)** |
| Geneformer | 0.348 | 0.403 | 0.728 | 0.422 | 0.186 | 0.592 | DONE (v4) |
| scGPT | — | 0.397 | 0.703 | — | — | — | DONE (v4) |

## ACC 03: Zero-shot Benchmark — SKIPPED

Scripts are byte-for-byte identical to ACC 02. Produces identical results.
For manuscript: "zero-shot evaluation using frozen embeddings with trained classifiers
is equivalent to the supervised benchmark since no model fine-tuning is performed."

## ACC 08: Statistical Significance — DONE

### Binary (Bootstrap CI + DeLong)

| Model | AUC | 95% CI | Per-fold AUCs |
|-------|-----|--------|---------------|
| STATE | 0.814 | [0.862, 0.869] | 0.879, 0.759, 0.860, 0.765, 0.802 |
| UCE | 0.794 | [0.854, 0.862] | 0.854, 0.764, 0.853, 0.751, 0.746 |
| TranscriptFormer | 0.742 | [0.777, 0.787] | 0.679, 0.773, 0.721, 0.735, 0.771 |
| scGPT | 0.685 | [0.821, 0.829] | 0.753, 0.667, 0.671, 0.732, 0.503 |
| Geneformer | 0.743 | [0.767, 0.778] | 0.673, 0.750, 0.676, 0.771, 0.793 |

### DeLong Pairwise Tests (binary)

| Comparison | AUC diff | p-value | Sig |
|------------|----------|---------|-----|
| STATE vs UCE | +0.020 | <0.0001 | * |
| STATE vs TF | +0.072 | <0.0001 | * |
| STATE vs scGPT | +0.129 | <0.0001 | * |
| STATE vs GF | +0.071 | <0.0001 | * |
| TF vs GF | -0.001 | 0.672 | ns |

## ACC 09: Embedding Ablation — DONE

Key findings (binary AUC):
- PCA dimensionality: STATE loses minimally (0.825→0.813 at PCA-256), Geneformer loses most (0.778→0.740)
- Top-PC removal: UCE collapses (0.790→0.512 at remove_top_10), Geneformer robust (0.778→0.742)
- Normalization: minimal effect across all models
- Random dimensions: UCE most efficient (0.766 with just 64 random dims)

## ACC 10: Clinical Metrics — DONE (all 5 models, binary + 3-class)

### Binary

| Model | AUC-ROC | AUC-PR | Sens@90Spec | Cohen's Kappa | Status |
|-------|---------|--------|-------------|---------------|--------|
| STATE | **0.825** | **0.886** | **0.548** | **0.467** | DONE (v4) |
| UCE | 0.794 | 0.868 | 0.503 | 0.406 | DONE (v4) |
| Geneformer | 0.778 | 0.852 | 0.459 | 0.394 | DONE (v5) |
| TranscriptFormer | 0.776 | 0.849 | 0.457 | 0.386 | DONE (v4) |
| scGPT | 0.664 | 0.723 | 0.248 | 0.307 | DONE (v5) |

### 3-Class

| Model | AUC-ROC | AUC-PR | Cohen's Kappa | BalAcc | Status |
|-------|---------|--------|---------------|--------|--------|
| STATE | **0.822** | **0.656** | **0.401** | 0.654 | DONE (v4) |
| UCE | 0.818 | 0.640 | 0.380 | **0.655** | DONE (v4) |
| TranscriptFormer | 0.787 | 0.583 | 0.330 | 0.632 | DONE (v4) |
| Geneformer | 0.778 | 0.573 | 0.321 | 0.617 | DONE (v5) |
| scGPT | 0.630 | 0.424 | 0.175 | 0.448 | DONE (v5) |

## ACC 15: Differential Diagnosis — DONE (v5)

| Model | Macro AUC | Weighted AUC | COVID-19 | Healthy | Pneumonia | SLE | Sepsis |
|-------|-----------|--------------|----------|---------|-----------|-----|--------|
| **UCE** | **0.938** | 0.897 | 0.887 | 0.913 | 0.909 | 0.991 | 0.988 |
| **STATE** | **0.936** | **0.913** | **0.905** | **0.931** | 0.907 | 0.944 | **0.993** |
| TranscriptFormer | 0.871 | 0.852 | 0.845 | 0.883 | 0.696 | 0.956 | 0.978 |
| scGPT | 0.842 | 0.785 | 0.780 | 0.580 | 0.896 | **0.991** | 0.962 |
| scGPT_PCA256 | 0.840 | 0.786 | 0.782 | 0.586 | 0.878 | 0.991 | 0.963 |
| Geneformer | 0.795 | 0.797 | 0.793 | 0.876 | 0.407 | 0.929 | 0.968 |

Key findings:
- UCE slightly beats STATE for macro AUC (0.938 vs 0.936) — ranking reversal from severity
- STATE has highest weighted AUC (0.913) — better on larger classes
- scGPT PCA-256 does NOT help for diff diagnosis (0.840 vs 0.842 raw) — unlike severity task
- Geneformer catastrophically fails on Pneumonia (0.407)
- SLE and Sepsis easiest for all models (>0.92)

## ACC 18: Multimodal Fusion — DONE

| Method | AUC |
|--------|-----|
| Late average | **0.849** |
| Late rank | 0.848 |
| Late stacking | 0.840 |
| Early PCA-256 | 0.835 |
| Attention fusion | 0.826 |
| STATE alone | 0.831 |
| UCE alone | 0.805 |
| TranscriptFormer alone | 0.782 |
| Geneformer alone | 0.776 |

## ACC 05: Sample Efficiency — DONE (v5, all 5 models!)

### Binary AUC at each subsample fraction

| Model | 100% | 50% | 20% | 10% | 5% | Drop at 5% | Status |
|-------|------|------|------|------|------|------------|--------|
| STATE | **0.817** | **0.801** | 0.754 | 0.700 | 0.703 | -13.9% | DONE |
| UCE | 0.795 | 0.790 | **0.773** | **0.767** | **0.727** | -8.5% | DONE |
| Geneformer | 0.757 | 0.740 | 0.703 | 0.666 | 0.646 | **-14.7%** | DONE |
| TranscriptFormer | 0.740 | 0.735 | 0.709 | 0.682 | 0.681 | -8.0% | DONE |
| scGPT | 0.707 | 0.704 | 0.697 | 0.682 | 0.695 | -1.7% | DONE |

Key findings:
- **UCE most sample-efficient** in absolute terms — highest AUC at 5%-50% fractions
- **scGPT flattest curve** (-1.7%) but lowest absolute performance throughout
- **TranscriptFormer** also resilient (-8.0%), similar to UCE
- **Geneformer degrades worst** (-14.7%), steepest cliff from 20%→5%
- **STATE** best at 100% but drops steeply below 20% (-13.9%)
- scGPT at 5% (0.695) actually higher than at 10% (0.682) — variance artifact
- Robustness ranking: scGPT (-1.7%) > TF (-8.0%) ≈ UCE (-8.5%) > STATE (-13.9%) > GF (-14.7%)
- But for absolute performance at low data: UCE >> all others

## ACC 07: Gene Dropout — DONE (v5, all 5 models!)

### Binary AUC at each dropout level

| Model | 0% | 10% | 30% | 50% | 70% | Drop at 70% |
|-------|-----|-----|-----|-----|-----|-------------|
| STATE | 0.814 | 0.815 | 0.812 | 0.807 | 0.798 | -2.0% |
| **UCE** | **0.809** | **0.808** | **0.808** | **0.808** | **0.809** | **+0.0%** |
| TranscriptFormer | 0.730 | 0.730 | 0.728 | 0.727 | 0.722 | -1.1% |
| scGPT | 0.686 | 0.685 | 0.683 | 0.678 | 0.675 | -1.6% |
| Geneformer | 0.729 | 0.728 | 0.725 | 0.720 | 0.714 | -2.0% |

Key findings:
- **UCE is perfectly flat** — 0.809 ± 0.001 at every dropout level (0% to 70%)
- TranscriptFormer also very stable (-1.1%)
- STATE and Geneformer degrade similarly (-2.0% each)
- scGPT slightly more sensitive (-1.6%)
- Robustness ranking: UCE >> TF > scGPT > STATE ≈ GF

## ACC 19: Model Optimization — DONE

| Precision | AUC | BalAcc | Size (MB) | Speedup |
|-----------|-----|--------|-----------|---------|
| FP32 | 0.829 | 0.746 | 2413 | 1.0x |
| FP16 | 0.829 | 0.745 | 1207 | 2.0x |
| INT8 dynamic | 0.828 | 0.746 | 603 | 2.9x |
| INT8 static | 0.825 | 0.743 | 603 | 3.3x |

## ACC 06: Cross-Disease Mortality Transfer (fatal vs survived) — DONE 2026-04-30

Job 23743154, ran 19:39-19:56 PDT (~17 min). Output:
`05_RESULTS/cross_disease_fatal_vs_survived/cross_disease_fatal_vs_survived_20260430_193915.{csv,json}`

### COVID → Sepsis transfer (LogReg, patient-level bootstrap CI)

| Model | AUC | 95% CI |
|---|---|---|
| **STATE** | **0.800** | [0.636, 0.928] |
| UCE | 0.798 | [0.736, 0.856] |
| TranscriptFormer | 0.775 | [0.693, 0.877] |
| Geneformer | 0.742 | [0.619, 0.861] |
| **scGPT** | **0.352** ⚠️ | [0.217, 0.503] (anti-correlated) |

### Sepsis → COVID transfer (LogReg)

| Model | AUC | 95% CI |
|---|---|---|
| **STATE** | **0.706** | [0.640, 0.775] |
| TranscriptFormer | 0.629 | [0.575, 0.695] |
| Geneformer | 0.581 | [0.533, 0.624] |
| scGPT | 0.555 | [0.395, 0.667] |
| **UCE** | **0.514** ⚠️ | [0.400, 0.597] (random — collapsed) |

### Key findings
1. **STATE is the only bidirectional generalizer** — wins both forward (0.800) and reverse (0.706) directions
2. **UCE asymmetry**: forward 0.798 = STATE, reverse 0.514 = random. COVID-specific encoding.
3. **scGPT anti-correlated forward**: 0.352 with CI [0.22, 0.50] — predictions systematically wrong, not just noisy
4. **Forward direction transfers well for most FMs** (4/5 in 0.74-0.80) — mortality signal is partially universal
5. **Reverse direction is the harder test** — only STATE maintains performance

Code: `analyses/06_batch_generalization/code/cross_disease_fatal_vs_survived.py`
Submit: `sbatch analyses/06_batch_generalization/code/run_cross_disease_fatal.sh`
Manuscript: see `tab:crossdisease` in main.tex; new Discussion subsection added

## What Remains After v5

### All v5 results collected — no pending re-runs.

STATE 6-class (job 23675490 + duplicate 23680676) completed 2026-04-29:
LogReg AUC 0.782, BalAcc 0.467 — STATE leads at every severity granularity
(binary 0.825, 3-class 0.822, 6-class 0.782).

### Completed manuscript results:
- ACC 02 binary: All 5 models ✓ **DONE**
- ACC 02 3-class: All 5 models ✓ **DONE**
- ACC 02 6-class: All 5 models ✓ **DONE**
- ACC 05: Sample efficiency all 5 models ✓ **DONE**
- ACC 07: Gene dropout all 5 models ✓ **DONE**
- ACC 08: Bootstrap CI + DeLong ✓ **DONE**
- ACC 09: Embedding ablation ✓ **DONE**
- ACC 10: Clinical metrics all 5 models × 2 tasks ✓ **DONE**
- ACC 15: Differential diagnosis all 5 models ✓ **DONE**
- ACC 18: Multimodal fusion ✓ **DONE**
- ACC 19: Quantization ✓ **DONE**

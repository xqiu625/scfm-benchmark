# Embedding Classification Results Summary

**Date:** 2026-01-26
**Method:** Downstream classification on frozen embeddings
**Strategies:** Linear Probe, DeepMLP, ResidualMLP, Attention
**Cross-validation:** 5-fold stratified

## Method Description

This analysis trains **classifiers on pre-computed, frozen embeddings** from foundation models. This is fundamentally different from fine-tuning:

```
Foundation Model Fine-tuning:
[Raw Data] → [Foundation Model (trainable)] → [Classifier] → [Prediction]

Our Method (Embedding Classification):
[Raw Data] → [Foundation Model] → [Frozen Embeddings] → [Trainable Classifier] → [Prediction]
                    ↓
              (pre-computed,
               not modified)
```

**Key distinction:**
- Foundation model weights are **never modified**
- Only the classification head is trained
- Embeddings are pre-computed and fixed
- This complies with model licenses (e.g., UCE restrictions)

## Results by Task

### Binary Classification (Severe/Fatal vs Others)

| Model | Zero-Shot | Linear | DeepMLP | ResMLP | Attention | **Best** | **Gain** |
|-------|-----------|--------|---------|--------|-----------|----------|----------|
| **STATE** | 0.895 | 0.884 | **0.942** | 0.940 | 0.913 | 0.942 | **+4.7%** |
| UCE | 0.877 | 0.846 | 0.906 | **0.908** | 0.884 | 0.908 | +3.1% |
| Geneformer | 0.813 | 0.766 | 0.847 | **0.849** | 0.817 | 0.849 | +3.6% |
| TranscriptFormer | 0.838 | 0.788 | **0.875** | 0.866 | 0.848 | 0.875 | +3.7% |
| scGPT | 0.775 | 0.729 | **0.732** | 0.731 | 0.728 | 0.732 | -4.3% |
| scGPT_BP | 0.766 | 0.576 | **0.771** | 0.749 | 0.521 | 0.771 | +0.5% |

### 3-Class Classification (Healthy vs Mild/Moderate vs Severe/Fatal)

| Model | Zero-Shot | Linear | DeepMLP | ResMLP | Attention | **Best** | **Gain** |
|-------|-----------|--------|---------|--------|-----------|----------|----------|
| **STATE** | 0.909 | 0.874 | 0.938 | **0.940** | 0.882 | 0.940 | **+3.1%** |
| UCE | 0.896 | 0.833 | **0.891** | 0.884 | 0.859 | 0.891 | -0.5% |
| Geneformer | 0.826 | 0.760 | **0.833** | 0.822 | 0.644 | 0.833 | +0.7% |
| TranscriptFormer | 0.851 | 0.776 | 0.818 | 0.817 | **0.820** | 0.820 | -3.1% |
| scGPT | 0.747 | **0.671** | 0.669 | 0.671 | 0.663 | 0.671 | -7.6% |
| scGPT_BP | 0.784 | 0.563 | **0.716** | 0.703 | 0.514 | 0.716 | -6.8% |

### 6-Class Classification (All severity levels)

| Model | Zero-Shot | Linear | DeepMLP | ResMLP | Attention | **Best** | **Gain** |
|-------|-----------|--------|---------|--------|-----------|----------|----------|
| **STATE** | 0.894 | 0.882 | **0.951** | 0.949 | 0.912 | 0.951 | **+5.7%** |
| UCE | 0.888 | 0.835 | **0.907** | 0.900 | 0.885 | 0.907 | +1.9% |
| Geneformer | 0.821 | 0.756 | **0.835** | 0.827 | 0.738 | 0.835 | +1.4% |
| TranscriptFormer | 0.849 | 0.781 | **0.860** | 0.851 | 0.827 | 0.860 | +1.1% |
| scGPT | 0.810 | 0.751 | 0.752 | **0.752** | 0.749 | 0.752 | -5.8% |
| scGPT_BP | 0.768 | 0.590 | **0.753** | 0.740 | 0.515 | 0.753 | -1.5% |

## Key Findings

### 1. STATE Consistently Best
- Achieves highest AUC across all tasks
- Best improvement: +5.7% on 6-class (0.894 → 0.951)
- Well-separated clusters in embedding space enable effective classifier training

### 2. DeepMLP is the Best Strategy
- Outperforms other strategies for most model/task combinations
- ResidualMLP is competitive but slightly worse
- Linear probe underperforms (often worse than zero-shot)
- Attention strategy struggles, especially on scGPT models

### 3. scGPT Models Hurt by Classification Training
- scGPT: Consistently negative gains (-4.3% to -7.6%)
- scGPT_BP: Mixed results, mostly negative (-6.8% to +0.5%)
- Root cause: Overlapping clusters in embedding space (low silhouette scores)
- Classifier training amplifies confusion between intermediate severity classes

### 4. Model Rankings (by best classifier AUC)

**Binary:**
1. STATE (0.942)
2. UCE (0.908)
3. TranscriptFormer (0.875)
4. Geneformer (0.849)
5. scGPT_BP (0.771)
6. scGPT (0.732)

**6-Class:**
1. STATE (0.951)
2. UCE (0.907)
3. TranscriptFormer (0.860)
4. Geneformer (0.835)
5. scGPT_BP (0.753)
6. scGPT (0.752)

## Relationship to Embedding Quality

| Model | Mean Silhouette | Classification Gain (6-class) |
|-------|-----------------|-------------------------------|
| STATE | -0.009 (best) | +5.7% (best) |
| UCE | -0.013 | +1.9% |
| Geneformer | -0.013 | +1.4% |
| TranscriptFormer | -0.021 | +1.1% |
| scGPT | -0.053 | -5.8% |
| scGPT_BP | -0.071 (worst) | -1.5% |

**Correlation:** Better cluster separation (higher silhouette) → better classification gains

## Output Files

- Results CSV: `scripts/02_embedding_classification/results/benchmark_embedding_classifier_*.csv`
- Checkpoints: `/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/checkpoints/finetune_v2/`

## Code

- Script: `scripts/finetuning/unified_finetune_v2.py`
- Submit: `scripts/02_embedding_classification/code/submit_finetune_v2.sh`

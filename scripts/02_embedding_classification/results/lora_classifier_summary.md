# LoRA-Style Classifier Results

**Date:** 2026-01-26
**Method:** LoRA-style classifier on frozen embeddings
**GPU:** Ada 6000 (short_gpu partition)

## Method Description

This analysis uses **LoRA-style (Low-Rank Adaptation) layers in the classification head**, trained on frozen embeddings. This is NOT fine-tuning of the foundation models.

```
Our Method:
[Pre-computed Embeddings] → [LoRA-style Classifier] → [Predictions]
       (frozen)                  (trainable)
```

**What LoRA-style means here:**
- Linear layers with low-rank adapter matrices (A and B)
- Main linear weights are frozen, only adapters are trained
- `output = Wx + (BA)x * scaling` where W is frozen, A and B are trainable
- Parameter efficient: ~1% of parameters are trainable

**Key distinction from true LoRA fine-tuning:**
- True LoRA: Injects adapters into transformer attention layers of the foundation model
- Our method: Uses LoRA-style layers only in the classification head
- Foundation model embeddings are **pre-computed and frozen**
- Complies with model licenses (e.g., UCE restrictions)

## Parameters
- LoRA rank (r): 8
- LoRA alpha: 16
- Dropout: 0.1
- Cross-validation: 5-fold stratified
- Loss: Focal Loss (gamma=2.0)
- Early stopping: patience=10

## Results - All Models, All Tasks

### Binary Classification

| Model | Zero-shot AUC | LoRA-Classifier AUC | Change |
|-------|---------------|---------------------|--------|
| STATE | 0.895 | 0.921 | **+2.6%** |
| UCE | 0.877 | 0.897 | **+2.0%** |
| TranscriptFormer | 0.838 | 0.876 | **+3.8%** |
| Geneformer | 0.813 | 0.846 | **+3.3%** |
| scGPT | 0.775 | 0.734 | -4.1% |
| scGPT_BP | 0.766 | 0.523 | -24.3% |

### 3-Class Classification

| Model | Zero-shot AUC | LoRA-Classifier AUC | Change |
|-------|---------------|---------------------|--------|
| STATE | 0.909 | 0.922 | **+1.3%** |
| UCE | 0.896 | 0.911 | **+1.5%** |
| TranscriptFormer | 0.851 | 0.887 | **+3.6%** |
| Geneformer | 0.826 | 0.849 | **+2.3%** |
| scGPT | 0.747 | 0.698 | -4.9% |
| scGPT_BP | 0.784 | 0.555 | -23.0% |

### 6-Class Classification

| Model | Zero-shot AUC | LoRA-Classifier AUC | Change |
|-------|---------------|---------------------|--------|
| STATE | 0.894 | 0.930 | **+3.6%** |
| UCE | 0.888 | 0.917 | **+2.9%** |
| TranscriptFormer | 0.849 | 0.896 | **+4.7%** |
| Geneformer | 0.821 | 0.853 | **+3.2%** |
| scGPT | 0.810 | 0.772 | -3.8% |
| scGPT_BP | 0.768 | 0.619 | -14.9% |

## Key Findings

1. **4/6 models improved** with LoRA-style classification across all tasks
2. **TranscriptFormer** consistently had the biggest gains (+3.6% to +4.7%)
3. **STATE** achieved highest overall AUC (0.930 on 6-class)
4. **scGPT models hurt by this approach** - consistent with other classification methods
5. **scGPT_BP severely degraded** (-14.9% to -24.3%) - LoRA-style layers don't help overlapping embeddings
6. **Best performers** (ranked by 6-class AUC):
   - STATE: 0.930
   - UCE: 0.917
   - TranscriptFormer: 0.896
   - Geneformer: 0.853

## Comparison: LoRA-Classifier vs DeepMLP Classifier (6-class)

| Model | Zero-shot | LoRA-Classifier | DeepMLP | Best Method |
|-------|-----------|-----------------|---------|-------------|
| STATE | 0.894 | 0.930 | **0.951** | DeepMLP |
| UCE | 0.888 | 0.917 | **0.907** | LoRA |
| TranscriptFormer | 0.849 | **0.896** | 0.860 | LoRA |
| Geneformer | 0.821 | 0.853 | **0.835** | LoRA |
| scGPT | 0.810 | 0.772 | **0.752** | LoRA (both bad) |
| scGPT_BP | 0.768 | 0.619 | **0.753** | DeepMLP |

**Observation:** LoRA-classifier outperforms DeepMLP for UCE, TranscriptFormer, and Geneformer, but DeepMLP is better for STATE and scGPT_BP.

## Output Files

Results saved to: `/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/lora_finetuning/`
- `{model}_lora_{task}_*.json` for each model and task combination

## Code

- Single model: `scripts/finetuning/geneformer_lora_finetune.py`
- All models: `scripts/finetuning/unified_lora_finetune.py`
- Submit scripts: `scripts/02_embedding_classification/code/submit_*.sh`

## Checkpoint Support

Checkpoints saved to: `05_RESULTS/checkpoints/lora/`
- Auto-saves every 5 epochs and when best AUC improves
- Auto-resumes from checkpoint if job restarts
- Enables running on short_gpu (2-hour limit) with auto-resume

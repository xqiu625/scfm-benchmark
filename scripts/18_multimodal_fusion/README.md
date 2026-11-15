# ACC 18: Multimodal Fusion

## Goal
Demonstrate multi-view learning by fusing heterogeneous FM representations.

## Strategies

### 1. Early Fusion (concatenate + PCA)
- Concat STATE (2058) + UCE (1280) + Geneformer (1152) + TranscriptFormer (2048) = 6538 dims
- PCA to 256 dims before classification

### 2. Late Fusion (average + rank + stacking)
- Average prediction probabilities across models
- Rank-based fusion for robustness
- Stacking meta-learner

### 3. Attention-Based Fusion (key differentiator)
- Project each model to 256 dims via learned linear layers
- Cross-attention: learned query attends over model embeddings
- Output: fused embedding + interpretable per-cell attention weights
- ~500K params, trained with classification loss
- **Key output:** attention heatmap showing which FM matters for which cells

### 4. Expression + Embedding Fusion
- Concatenate STATE embedding with top-100 raw gene expression features
- Tests whether FMs capture all info or raw features add value

## Execution
```bash
sbatch scripts/18_multimodal_fusion/code/submit_18_fusion.sh
```

## Expected Results
- Attention fusion beats best single model (STATE 0.944 AUC)
- Attention weights reveal per-cell model specialization (severe cells -> STATE, mild -> UCE)

## Env
- `state_env` (PyTorch 2.5.1) | GPU partition, 1 GPU, 128GB, 12h

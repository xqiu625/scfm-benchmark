# ACC 22: ESM Protein Multimodal Fusion

Integrates protein-level representations (ESM-2) with cell-level foundation model
embeddings for multimodal severity prediction.

## Concept

Current multimodal fusion combines 4 cell-level FMs (STATE, UCE, Geneformer, TranscriptFormer).
This module adds a 5th modality: **protein structural/functional embeddings** from ESM-2,
aggregated per cell by expression-weighted attention.

For each cell:
1. Identify which genes are expressed
2. Look up ESM-2 embeddings for corresponding proteins (1280-dim)
3. Aggregate protein embeddings weighted by expression level
4. Use as 5th input to attention fusion

## Pipeline

```
Step 1: identify_degs.py        — Find top DEGs across severity classes
Step 2: get_protein_sequences.py — Retrieve UniProt sequences for DEGs
Step 3: generate_esm_embeddings.py — Generate ESM-2 protein embeddings
Step 4: build_cell_protein_embeddings.py — Per-cell protein features via expression-weighted aggregation
Step 5: multimodal_fusion_5model.py — Attention fusion with 5 modalities
```

## Expected Output
- Per-protein ESM embeddings (top 500 DEGs x 1280-dim)
- Per-cell protein embeddings (47K cells x 256-dim after projection)
- 5-modality attention fusion AUC (compare vs 4-modality baseline)

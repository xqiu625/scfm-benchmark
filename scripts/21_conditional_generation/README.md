# ACC 21: Conditional Generation in Embedding Space

## Goal
Demonstrate generative modeling connecting disease states to therapeutic interventions.

## Architecture
- **Conditional VAE** on STATE embeddings (2058 dims) + severity condition (6-class one-hot)
- Encoder: [2058+6] -> 512 -> 256 -> 128 (mu, logvar)
- Decoder: [128+6] -> 256 -> 512 -> 2058
- Training: beta-VAE loss (recon + 0.5 * KL), 200 epochs

## Key Experiment: Severity Transition Generation
1. Encode severe cells with "severe" condition
2. Decode from same latent with "recovered" condition
3. Evaluate:
   - L2 distance from generated to nearest real recovered cells
   - Cosine similarity of transition vector to real recovery direction
   - UMAP overlay: generated vs real recovered cells
   - Classifier test: LogReg classifies generated embeddings

## Execution
```bash
sbatch scripts/21_conditional_generation/code/submit_21_cvae.sh
```

## Expected Results
- Generated "recovered" embeddings land near real recovered cluster on UMAP
- Cosine similarity > 0.5 with real recovery direction
- >70% of generated cells classified as recovered by LogReg

## Env
- `state_env` | GPU partition, 1 GPU, 64GB, 4h

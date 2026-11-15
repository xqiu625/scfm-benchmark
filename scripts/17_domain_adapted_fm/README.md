# ACC 17: Domain-Adapted Foundation Model

## Goal
Transform from "evaluator" to "builder" of foundation models by domain-adapting existing FMs and training a new one from scratch.

## Approaches

### A: Continue Pre-Training scGPT (scgpt_continue_pretrain.py)
- Load scGPT whole-human checkpoint
- Continue pre-training with MGEM (masked gene expression modeling) on 47,808 platelet cells
- 50 epochs, lr=5e-5, 15% masking
- Output: `X_scGPT_adapted` embeddings (512 dims)
- **Env:** `scpeft_env` | **GPU:** A100 x1

### B: Train Platelet Transformer from Scratch (platelet_transformer.py)
- 6-layer transformer, 384 hidden, 6 heads (~20M params)
- Pre-train on platelet (48K) + X-Atlas (500K) with MGEM
- 100 epochs, lr=1e-3, gradient accumulation (effective batch 256)
- Output: `X_platelet_tf` embeddings (384 dims)
- **Env:** `state_env` | **GPU:** A100 x1

### C: Benchmark (benchmark_adapted_models.py)
- Run 7-model comparison: 5 original + scGPT_adapted + PlateletTF
- Binary, 3-class, 6-class + sample efficiency

## Execution Order
```bash
# Run A and B in parallel
sbatch scripts/17_domain_adapted_fm/code/submit_17a_scgpt_pretrain.sh
sbatch scripts/17_domain_adapted_fm/code/submit_17b_train_scratch.sh

# After both complete: extract embeddings then benchmark
sbatch scripts/17_domain_adapted_fm/code/submit_17c_benchmark.sh
```

## Expected Results
- scGPT_adapted improves over raw scGPT (0.735 AUC), especially on OOD (0.524)
- PlateletTF competitive despite 20x fewer params than STATE (603M)

## Key Narrative
> "I rescued scGPT from OOD failure by domain adaptation and trained a 20M-param transformer from scratch that rivaled 600M-param models on platelet classification."

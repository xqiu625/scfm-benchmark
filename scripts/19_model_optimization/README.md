# ACC 19: Model Optimization & Distributed Training

## Goal
Demonstrate production ML engineering skills: quantization, FSDP, torch.compile.

## Components

### A: Quantization Benchmark (quantization_benchmark.py)
- STATE SE-600M (603M params) at FP32, FP16, INT8 dynamic, INT8 static
- Measures: inference latency, peak GPU memory, AUC
- **Expected:** FP16 maintains AUC within 0.1%, 2x speedup; INT8 ~4x memory savings

### B: FSDP Distributed Fine-Tuning (fsdp_finetune.py)
- Wrap SEClassifier with FullyShardedDataParallel
- SHARD_GRAD_OP strategy + mixed precision
- Compare 1-GPU vs 2-GPU scaling

### C: torch.compile Inference (torch_compile_benchmark.py)
- Apply torch.compile(mode='reduce-overhead') to all 5 FM inference pipelines
- Measure warmup vs steady-state throughput
- **Expected:** 20-40% speedup on steady-state inference

### D: Report (optimization_report.py)
- Consolidate into publication figures: precision vs accuracy, FSDP scaling curve, compile speedup chart

## Execution
```bash
# A + C can run in parallel on 1 GPU
sbatch scripts/19_model_optimization/code/submit_19a_quantization.sh
sbatch scripts/19_model_optimization/code/submit_19c_compile.sh

# B needs 2 GPUs
sbatch scripts/19_model_optimization/code/submit_19b_fsdp.sh
```

## Env
- `state_env` (PyTorch 2.5.1 supports FSDP + torch.compile)
- GPU partition, A100 x1-2, 128-256GB, 4-12h

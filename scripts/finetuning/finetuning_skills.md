# Generalizable Computational Skills
## Extracted from Project 302006 Memory Bank

---

## 1. Data Integrity & Validation

### Never Use Synthetic Fallbacks
- If a computation fails, **raise an error** instead of generating placeholder/dummy data
- Never use `np.random.randn()` or similar as fallback outputs
- Never assume missing labels or metadata -- raise an error and fix upstream
- All conclusions must be traceable to real, verified data

### Validate Outputs at Every Stage
- After generating any intermediate result, verify: correct dimensions, no NaN values, correct row/column counts, metadata preserved
- Create dedicated validation scripts to check data integrity after each pipeline step
- Compare input vs output cell/row counts to catch silent data loss

---

## 2. Experiment Tracking & Reproducibility

### Memory Bank System
Maintain a set of living documents for any multi-week computational project:

| Document | Purpose |
|----------|---------|
| **projectbrief.md** | Executive summary, objectives, dataset description, key results |
| **progress.md** | Detailed status of all tasks, results tables, job logs |
| **activecontext.md** | Current focus, immediate action items, latest results |
| **techcontext.md** | Environment details, library versions, directory structure |
| **PLAN.md** | Phased execution plan with tasks, dependencies, deliverables |

### Version Everything
- Track library versions explicitly (e.g., `scikit-learn 1.7.1`, `torch 2.0+`)
- Document conda environment names and paths
- Record SLURM job IDs and parameters for every computation
- Use git for all code; tag commits for key milestones

---

## 3. Performance Optimization Patterns

### Parallelization
- Use `joblib` for embarrassingly parallel workloads (bootstrap iterations, cross-validation folds)
- Set parallelism to match available cores (e.g., 16 cores for bootstrap, 8 for CV folds)
- Use SLURM array jobs for batch processing independent chunks (`--array=1-N%M` to limit concurrency)

### Algorithmic Speedups
- Choose faster solvers when accuracy trade-off is acceptable (e.g., `lbfgs` 10x faster than `saga` for logistic regression)
- Apply dimensionality reduction (PCA) before expensive downstream operations to reduce input size
- Reduce iteration counts where diminishing returns apply (e.g., 100 bootstrap iterations vs 500 when CIs are already tight)

### Achieved Speedups Reference
| Technique | Before | After | Speedup |
|-----------|--------|-------|---------|
| Parallel bootstrap + fast solver + PCA | 72h | 3 min | 1,440x |
| Parallel CV + fast solver | 24h | 25 min | 60x |

---

## 4. HPC / Cluster Job Management

### SLURM Best Practices
```bash
# Standard SLURM header template
#SBATCH --partition=<partition>
#SBATCH --mem=<memory>
#SBATCH --cpus-per-task=<cores>
#SBATCH --time=<walltime>
#SBATCH --output=/path/to/logs/%x_%j.out
#SBATCH --error=/path/to/logs/%x_%j.err
#SBATCH --mail-user=<email>
```

- **Use centralized log directory** -- all stdout/stderr in one place for easy monitoring
- **Set `PYTHONNOUSERSITE=1`** in SLURM scripts to avoid conflicts with user `~/.local` packages
- **Choose partition by need**: CPU-heavy (batch/highmem), GPU (gpu), quick iteration (short/short_gpu)
- **Use array jobs** for independent tasks: `--array=1-188%10` processes 188 chunks, max 10 concurrent
- **Skip completed work**: check for existing output files before re-running

### Checkpoint & Resume Pattern
- Save checkpoints periodically (every N epochs + on best metric improvement)
- Store: model state, optimizer state, scheduler state, epoch number, best metric
- Auto-detect and resume from checkpoint on job restart
- Essential for preemptible/time-limited partitions (e.g., 2-hour short_gpu jobs)

### Partition Selection Guide
| Partition | Use Case | Typical Limits |
|-----------|----------|---------------|
| **batch** | Standard CPU jobs | 128GB, 7 days |
| **highmem** | Large data loading, merging | 256-512GB |
| **short** | Quick CPU debugging | 256GB, short wall time |
| **gpu** | Training, embedding generation | A100, 100GB+, 7 days |
| **short_gpu** | Quick GPU iteration | 2-hour limit, fast queue |

---

## 5. Phased Project Execution

### Structure Work in Phases with Clear Dependencies
```
Phase 1: Data Preparation (can run in parallel with Phase 2)
  Task 1.1 -> Task 1.2 (sequential dependency)

Phase 2: Model Processing (independent of Phase 1)
  Task 2.1 (test) -> Task 2.2 (batch) -> Task 2.3 (combine)

Phase 3: Analysis (depends on Phase 1 + 2)
  Task 3.1 -> Task 3.2

Phase 4: Writing (can overlap with Phase 3)
```

### Identify Parallel Opportunities
- Tasks with no data dependencies can run concurrently
- Writing (methods section) can proceed while computation runs
- Use a critical path analysis to find the bottleneck

### Risk Mitigation
- Have a **backup plan** for each risky step (e.g., publish 2-model paper if 3rd model fails)
- **Cell ID / row matching issues**: document normalization strategy, fall back to common subset
- **Negative results are still publishable**: frame as "models achieve similar performance despite different architectures"

---

## 6. Benchmarking & Evaluation Framework

### Multi-Metric Evaluation
Never rely on a single metric. Use a comprehensive evaluation suite:

| Category | Metrics |
|----------|---------|
| **Accuracy** | AUC-ROC, AUC-PR, F1-macro, F1-weighted, accuracy |
| **Robustness** | Sample efficiency curves, feature dropout, batch-shift generalization |
| **Statistical** | Bootstrap confidence intervals, non-overlapping CIs for significance |
| **Clinical/Domain** | Sensitivity at fixed specificity, Cohen's kappa, calibration |
| **Ablation** | Dimensionality reduction, feature importance, normalization sensitivity |

### Baseline Comparison
- Always include a **simple baseline** (PCA + LogReg, XGBoost on raw features)
- Strong baselines prevent overstating model contributions (e.g., XGBoost at 0.897 vs foundation model at 0.944)
- Include zero-shot / no-training baselines where applicable

### Statistical Significance
- Use bootstrap confidence intervals (100-1000 iterations) with parallel computation
- Non-overlapping 95% CIs indicate statistically significant differences
- Report CIs alongside point estimates in all results tables

---

## 7. Environment & Dependency Management

### One Environment Per Model/Tool
- Maintain separate conda environments for each major tool to avoid dependency conflicts
- Document environment paths, Python versions, and key package versions
- Use setup scripts to recreate environments (`setup_<tool>_env.sh`)

### SLURM Script Environment Setup Template
```bash
#!/bin/bash
# Activate specific environment
source activate <env_name>

# Prevent user site-packages from interfering
export PYTHONNOUSERSITE=1

# Precision optimization for GPU
# Use torch.cuda.amp.autocast() for mixed precision
```

### Container Strategy
- Use Singularity/Docker containers for tools with complex dependencies
- Bind host directories into containers: `singularity exec --bind /data:/data container.sif python script.py`
- Keep container images versioned and documented

---

## 8. Terminology & Scientific Rigor

### Use Precise Method Names
Define clear terminology upfront and enforce it throughout:

| Method | What Happens | Correct Term |
|--------|-------------|--------------|
| Update model weights via backprop | Full model training | **Fine-tuning** |
| Train classifier on frozen embeddings | Only classifier trained | **Embedding classification** |
| Use model output directly, no training | Nothing trained | **Zero-shot / Feature extraction** |

### Decision Tree for Method Naming
```
Did you update the pretrained model weights?
  YES -> "Fine-tuning" or "End-to-end training"
  NO  -> Did you train any classifier head?
         YES -> "Embedding classification" / "Linear probing"
         NO  -> "Zero-shot" / "Feature extraction"
```

---

## 9. Debugging & Noise Diagnosis

### Dimensionality Reduction as a Diagnostic Tool
- If a model **improves with PCA reduction**, it likely has noisy embeddings (signal buried in noise dimensions)
- If a model **degrades with PCA**, it uses all dimensions effectively
- Removing top PCA components tests whether performance depends on dominant variance directions or distributed features
- Model **robust to top PCA removal** = distributed representation (good)
- Model **sensitive to top PCA removal** = concentrated representation (fragile)

### Normalization Sensitivity Test
- Apply L2, z-score, and raw normalization to features
- If all give similar results (<0.1% difference), normalization is irrelevant for this data
- Quick test to eliminate an entire category of hyperparameter tuning

### Random Feature Baseline
- Replace features with random subsets of same dimensionality
- Performance gap between real vs random features quantifies how much signal the model captures
- If random features perform well, the task may be too easy or leaking information

---

## 10. Results Documentation Pattern

### Standard Results Table Format
Always record results in structured markdown tables:

```markdown
| Rank | Method | Primary Metric | Secondary Metric | Notes |
|:----:|--------|:--------------:|:----------------:|-------|
| 1    | Best   | **0.944**      | 0.948            | Key insight |
| 2    | Second | 0.910          | 0.912            |             |
```

### Track What Was Tried and What Failed
- Document failed approaches with root cause analysis (e.g., "LoRA fine-tuning failed: AUC ~0.53, near random. Root cause: gene count limited to 512 by GPU memory")
- Record the **decision** (e.g., "Paused until more GPU memory available")
- This prevents repeating mistakes and informs future work

### File & Path Discipline
- Use a consistent directory structure (DATA / EMBEDDINGS / RESULTS / FIGURES / LOGS)
- Name output files with descriptive, sortable names
- Record exact file paths for all key outputs in documentation
- Use centralized results directories, not scattered across home folders

---

## 11. GPU Memory & Compute Optimization

### Memory-Constrained Training Strategies
- **Mixed precision** (`torch.cuda.amp.autocast()`) cuts memory ~50%
- **Gradient accumulation** simulates larger batch sizes
- **LoRA/PEFT adapters** train <2% of parameters while freezing backbone
- **Adaptive batch sizing** based on (samples x features) to stay within memory

### When Fine-Tuning Underperforms Pre-Computed Embeddings
Likely causes:
1. Memory limits force truncating input (e.g., 256 genes vs 12K)
2. Attention is O(n^2) -- full-length inputs may be infeasible
3. Pre-computed embeddings already encode full input via model's own optimized inference
4. Pre-training on large corpus (e.g., 167M samples) already captured task-relevant patterns

### Profiling
- Use `torch.profiler` + `nsight-sys` to identify bottlenecks
- Export runtime/memory per step as CSV for analysis
- Profile before optimizing -- measure, don't guess

---

## 12. Workflow Automation Patterns

### Idempotent Scripts
- Check for existing outputs before re-running expensive steps
- Scripts should be safe to re-run without side effects
- Use hash-based validation to verify cached results haven't been corrupted

### SLURM GPU Fallback
- Design scripts with CPU fallback when GPU unavailable
- Template: try A100 GPU partition first, fall back to highmem CPU with adjusted `--mem` and `--cpus-per-task`

### Efficiency Benchmarking Grid
- Automate parameter sweeps: `(model x precision x feature_count)` grid
- Save all results to structured CSV files
- Use shell scripts to orchestrate the full grid: `scripts/efficiency_bench.sh`

---

---

## Appendix A: Project Directory Structure (HPCC)

```
/bigdata/godziklab/shared/Xinru/302006/
├── 01_DATA/                  # Raw and processed data
│   ├── gpu_chunks_hvg/       # 188 HVG-filtered chunks (INPUT for Geneformer)
│   └── raw/                  # Original data files
│
├── 02_EMBEDDINGS/            # Foundation model embeddings
│   ├── state/
│   ├── uce/
│   │   ├── human_platelet_covid_severity_uce_adata.h5ad  (47,977 cells)
│   │   ├── hek293t_platelets_uce_adata.h5ad              (88,434 cells)
│   │   └── hct116_platelets_uce_adata.h5ad               (89,738 cells)
│   ├── transcriptformer/
│   ├── scgpt/
│   └── geneformer/
│
├── 03_MODELS/                # Trained models
│
├── 04_ANALYSIS/              # Analysis outputs
│
├── 05_RESULTS/               # Final results
│   ├── 15_differential_diagnosis/
│   ├── checkpoints/
│   │   ├── lora/             # LoRA fine-tuning checkpoints
│   │   └── finetune_v2/      # Unified fine-tuning checkpoints
│   ├── platelet_scgpt_analysis/
│   │   └── platelet_scgpt_combined.h5ad  (226,149 cells, 512-dim)
│   ├── geneformer_embeddings/
│   │   ├── covid_*_hvg_geneformer.h5ad   (9 files)
│   │   ├── hek293t_*_hvg_geneformer.h5ad (89 files)
│   │   └── hct116_*_hvg_geneformer.h5ad  (90 files)
│   ├── platelet_uce_scgpt_combined.h5ad
│   ├── platelet_geneformer_combined.h5ad
│   ├── platelet_three_models_combined.h5ad
│   ├── two_model_comparison/
│   └── three_model_comparison/
│       ├── figures/
│       └── tables/
│
├── 06_LOGS/                  # All SLURM job logs (stdout/stderr)
│   └── geneformer/           # Per-chunk logs: chunk_*.out, chunk_*.err
│
├── scfm/                     # Foundation model resources
│   ├── geneformer.sif        # Singularity container
│   └── UCE/                  # UCE framework
│
├── geneformer_models/
│   └── Geneformer/           # Pretrained model (316M params)
│       ├── config.json
│       ├── model.safetensors
│       └── geneformer/
│           └── token_dictionary_gc104M.pkl
│
├── scpeft_env/               # scPEFT conda environment (Python 3.9, PyTorch 2.1+cu118)
│
└── 302006_updategit/         # Git repository
    ├── scripts/
    │   ├── data/                  # Data download scripts
    │   ├── drug_discovery/        # Drug discovery pipeline
    │   ├── foundation_models/     # Embedding generation (per model)
    │   ├── preprocessing/         # Data preprocessing
    │   ├── download_results_from_hpcc.sh
    │   └── finetuning/            # Fine-tuning scripts
    │       ├── setup_scpeft.sh                # scPEFT env setup (→ scpeft_env)
    │       ├── scpeft_finetune.py             # Wrapper: --backbone scgpt|geneformer
    │       ├── scgpt_peft_finetune.py         # scGPT PEFT (LoRA/PREFIX/ENCODER/TOKEN)
    │       ├── geneformer_peft_finetune.py    # Geneformer LoRA
    │       ├── submit_scgpt_peft.sh           # SLURM: scGPT PEFT
    │       ├── submit_geneformer_peft.sh      # SLURM: Geneformer PEFT
    │       ├── state_true_finetune_v3.py      # STATE full fine-tuning (603M params)
    │       ├── state_lora_finetune.py         # STATE LoRA fine-tuning
    │       ├── submit_true_finetune.sh        # SLURM: STATE full fine-tuning
    │       ├── submit_lora_finetune.sh        # SLURM: STATE LoRA
    │       ├── submit_state_6class.sh         # SLURM: STATE 6-class only
    │       ├── verify_true_finetune_setup.py  # Verification: STATE setup
    │       ├── verify_scpeft_setup.py         # Verification: scPEFT setup
    │       └── finetuning_skills.md           # This document
    ├── scripts/                  # All analysis code (source of truth)
    │   ├── 01_embedding_integrity/
    │   ├── 02_embedding_classification/
    │   ├── ...
    │   ├── 15_differential_diagnosis/
    │   └── 16_embedding_generation_5models/
    ├── reports/manuscript/
    └── docs/
        └── runbooks/         # hpcc.md, gate1-gate2.md
```

**Git repo:** `/bigdata/godziklab/shared/Xinru/302006/302006_updategit/`
**Data root:** `/bigdata/godziklab/shared/Xinru/302006/`
**All SLURM logs:** `/bigdata/godziklab/shared/Xinru/302006/06_LOGS/`
**scPEFT conda env:** `/bigdata/godziklab/shared/Xinru/302006/scpeft_env`

---

## Appendix B: Conda Environments (Full Reference)

### Environment Summary

| Environment | Path | Python | Purpose |
|-------------|------|--------|---------|
| **302006_env** | `/bigdata/godziklab/shared/Xinru/302006_env` | 3.9.23 | Main project environment |
| **302006_env_clean** | `/rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean` | 3.9+ | **Recommended for analysis scripts** |
| **state_env** | `/bigdata/godziklab/shared/Xinru/302006/state_env` | 3.9+ | STATE model embedding generation |
| **geneformer_py310** | `/rhome/xinruqiu/bigdata/.conda/envs/geneformer_py310` | 3.10 | Geneformer embedding generation |
| **scgpt** | `/rhome/xinruqiu/bigdata/.conda/envs/scgpt` | 3.9+ | scGPT embedding generation |
| **transcriptformer** | `/rhome/xinruqiu/bigdata/.conda/envs/transcriptformer` | 3.9+ | TranscriptFormer embedding generation |
| **scpeft_env** | `/bigdata/godziklab/shared/Xinru/302006/scpeft_env` | 3.9 | scGPT/Geneformer PEFT fine-tuning (scPEFT) |
| **citeseq_analysis** | `/bigdata/nairlab/shared/XQiu/conda_envs/citeseq_analysis` | - | General analysis (legacy) |

### 302006_env_clean (Recommended for Analysis)

**Activation:** `source activate 302006_env_clean`

| Package | Version | Purpose |
|---------|---------|---------|
| matplotlib | 3.10.3 | Visualization |
| numpy | 2.2.6 | Numerical computing |
| pandas | 2.3.1 | Data manipulation |
| scanpy | 1.11.3 | Single-cell analysis |
| scikit-learn | 1.7.1 | Machine learning (classifiers, metrics, CV) |
| seaborn | 0.13.2 | Statistical visualization |
| h5py | - | HDF5 file I/O |
| joblib | - | Parallel processing |

### state_env (STATE Model)

**Activation:** `source activate /bigdata/godziklab/shared/Xinru/302006/state_env`

| Package | Version | Purpose |
|---------|---------|---------|
| numpy | 2.3.5 | Numerical computing |
| pandas | 2.3.3 | Data manipulation |
| scanpy | 1.11.5 | Single-cell analysis |
| scikit-learn | 1.7.2 | Machine learning |

### scpeft_env (scGPT/Geneformer PEFT Fine-Tuning via scPEFT)

**Activation:** `source activate /bigdata/godziklab/shared/Xinru/302006/scpeft_env`
**Path:** `/bigdata/godziklab/shared/Xinru/302006/scpeft_env`
**Setup script:** `scripts/finetuning/setup_scpeft.sh`

| Package | Version | Purpose |
|---------|---------|---------|
| python | 3.9 | Runtime |
| torch | 2.1.0+cu118 | GPU training (matches scPEFT requirements) |
| transformers | 4.45.2 | Geneformer BertForSequenceClassification |
| peft | 0.4.0 | Parameter-efficient fine-tuning adapters |
| accelerate | 0.21.0 | Training acceleration |
| scanpy | - | Single-cell data loading |
| scikit-learn | - | Cross-validation, metrics |
| torchtext | 0.16.0 | Required by scPEFT's bundled scGPT |
| loralib | - | LoRA adapter layers |
| ipython | - | Required by scgpt internally |
| flash-attn | - | Flash attention (optional) |

**Note:** Uses scPEFT's bundled scGPT (cloned to `03_MODELS/scPEFT/`). Set `PYTHONPATH=/bigdata/godziklab/shared/Xinru/302006/03_MODELS/scPEFT:$PYTHONPATH` to import scGPT.

### geneformer_py310 (Geneformer)

**Activation:** `source activate geneformer_py310`

| Package | Version | Purpose |
|---------|---------|---------|
| python | 3.10 | Runtime |
| geneformer | - | Geneformer package |
| torch | - | PyTorch |

### scgpt (scGPT Embeddings)

**Activation:** `source activate scgpt`

| Package | Version | Purpose |
|---------|---------|---------|
| scgpt | - | scGPT package |
| torch | - | PyTorch |
| flash-attn | - | Flash attention for efficient inference |

### transcriptformer (TranscriptFormer)

**Activation:** `source activate transcriptformer`

| Package | Version | Purpose |
|---------|---------|---------|
| transcriptformer | - | TranscriptFormer package |
| torch | - | PyTorch |

### Environment Usage Rules

1. **Always add to SLURM scripts:**
   ```bash
   source activate <env_name>
   export PYTHONNOUSERSITE=1   # Prevent ~/.local conflicts
   ```
2. **All packages stored in shared folder** -- do not install to home directory
3. **Geneformer and scGPT run via Singularity** for complex dependency isolation
4. **Use `302006_env_clean`** as default for analysis/visualization scripts
5. **Use model-specific envs** only when generating embeddings or fine-tuning that model

---

## Appendix C: Cluster Access & Interactive Sessions

### SSH Access
```bash
ssh -X xinruqiu@cluster.hpcc.ucr.edu
```

### SLURM Email Notifications
```bash
#SBATCH --mail-user=xqiu014@ucr.edu
```

### Interactive Session Commands

| Session Type | Command |
|-------------|---------|
| **High memory CPU** | `srun -p highmem --mem=512gb -c 4 --pty bash -l` |
| **Standard CPU** | `srun -p batch --mem=128gb -c 4 --pty bash -l` |
| **Short CPU** | `srun -p short --mem=256gb -c 4 --pty bash -l` |
| **Basic GPU** | `srun -p gpu --qos=gpu --gres=gpu:1 --mem=200G -c 4 --pty bash -l` |
| **A100 GPU** | `srun -p gpu --qos=gpu --gres=gpu:a100:1 --mem=64G -c 4 --pty bash -l` |
| **High-mem A100** | `srun -p gpu --qos=gpu --gres=gpu:a100:1 --mem=100G -c 4 --pty bash -l` |
| **Extended A100** | `srun -p gpu --qos=gpu --gres=gpu:a100:1 --mem=64G -c 4 --time=4:00:00 --pty bash -l` |

### SLURM Partitions

| Partition | Time Limit | Notes |
|-----------|-----------|-------|
| **gpu** | 7 days | Long wait (~30 hrs), A100/H100 GPUs |
| **short_gpu** | 2 hours | Fast start, private Ada 6000 GPUs |
| **preempt_gpu** | - | Not yet available for godziklab; available after Feb 2026 shutdown |
| **batch** | 7 days | Standard CPU |
| **highmem** | 7 days | 256-512GB RAM |
| **short** | Short | Quick CPU debugging |

### Singularity Container Usage
```bash
module load singularity
singularity exec --bind /bigdata:/bigdata \
  /bigdata/godziklab/shared/Xinru/302006/scfm/geneformer.sif \
  python script.py
```

---

## Appendix D: Embedding Specifications

| Model | obsm Key | Dimensions | Cells | File/Location |
|-------|----------|:----------:|:-----:|---------------|
| **STATE** | `X_emb` | 2,058 | 47,808 | 02_EMBEDDINGS/state/ |
| **UCE** | `X_uce` | 1,280 | 47,808 | 02_EMBEDDINGS/uce/ |
| **TranscriptFormer** | `X_transcriptformer` | 2,048 | 46,949 | 02_EMBEDDINGS/transcriptformer/ |
| **scGPT** | `X_scGPT` | 512 | 47,808 | `platelet_scgpt_combined.h5ad` |
| **scGPT_PCA256** | `X_scGPT` + PCA | 256 | 47,808 | Derived from scGPT |
| **Geneformer** | `X_geneformer` | 1,152 | 46,949 | 02_EMBEDDINGS/geneformer/ |

### UCE Embedding Generation Command
```bash
python eval_single_anndata.py \
  --adata_path "/path/to/dataset.h5ad" \
  --dir "/output/directory/" \
  --species "human" \
  --batch_size 200
```

Output: AnnData with UCE embeddings in `.obsm['X_uce']`

---

## Appendix E: Dataset Reference

### COVID-19 Platelet Single-Cell RNA-seq
- **Total cells**: 47,808 platelets
- **Severity levels**: 6 classes (healthy, mild, moderate, severe, fatal, recovered)
- **Distribution**: severe (41.4%), fatal (19.7%), mild (15.4%), moderate (9.0%), recovered (7.7%), healthy (6.7%)

### Classification Tasks
- **Binary**: COVID vs Healthy
- **3-class**: Healthy vs Mild/Moderate vs Severe/Fatal
- **6-class**: Full severity spectrum
- **4-class Differential Diagnosis**: COVID-19 (38,673), Sepsis (3,422), SLE (169), Influenza (33)

### Data Sources (11 Batches)

| Batch ID | Study | PMID |
|----------|-------|------|
| 302004data01 | Wilk et al. | 32514174 |
| 302004data02 | Lee et al. | 32651212 |
| 302004data04 | Schulte-Schrepping et al. | 32810438 |
| 302004data05 | Arunachalam et al. | 32788292 |
| 302004data15 | Ren et al. | 33657410 |
| 302004data16 | Liu et al. | 33713619 |
| 302004data17 | Stephenson et al. | 33879890 |
| 302004data18 | Jiang et al. | 32554932 |
| 302005data01 | Bernardes et al. | 33296687 |
| 302005data02 | Combes et al. | 33494096 |
| 302005data03 | Qiu et al. | 34558746 |

*302005data04 (Mistry et al., PMID: 31754027) is SLE data, excluded from COVID analysis.*

### Perturbation Datasets
- **HEK293T**: 88,434 cells (1,362 perturbations with >=100 cells)
- **HCT116**: 89,738 cells (1,500 perturbations with >=100 cells)
- **Mouse**: 97,047 cells (cross-species transfer target)
- **Ortholog mapping**: 85.6% success rate, 13,142 orthologs

---

## Appendix F: scGPT PEFT Fine-Tuning — Methods, Results & Progress

### Overview

Parameter-efficient fine-tuning (PEFT) of the scGPT foundation model for COVID-19 severity classification from platelet single-cell RNA-seq data. Uses the scPEFT framework (Nature Machine Intelligence, 2025) with LoRA adapters.

### Method

#### Model Architecture

| Component | Details |
|-----------|---------|
| **Base model** | scGPT whole-human (pretrained on 33M cells from CELLxGENE) |
| **Architecture** | 12-layer Transformer, d_model=512, 8 attention heads, 51M total params |
| **Pretrained weights** | `best_model.pt` (147/193 layers loaded; remaining are task-specific heads) |
| **PEFT method** | LoRA (Low-Rank Adaptation) via scPEFT framework |
| **LoRA config** | r=16, alpha=32, dropout=0.1, applied to attention layers |
| **Trainable params** | 635,906 / 51,430,915 (**1.24%**) |
| **Classification head** | Linear layer on CLS token output → n_classes |

#### Input Representation

| Parameter | Value |
|-----------|-------|
| **Input** | Raw gene expression counts per cell |
| **Tokenization** | Non-zero genes ranked by expression (descending), mapped to scGPT vocabulary (60,697 genes) |
| **Max sequence length** | 1,200 genes per cell |
| **Padding** | Zero-padded to max length; padding mask applied |
| **Vocabulary** | scGPT whole-human `vocab.json` |

#### Training Configuration

| Parameter | Value |
|-----------|-------|
| **Optimizer** | AdamW (lr=1e-4, weight_decay=0.01) |
| **Scheduler** | Cosine annealing (T_max=20, eta_min=1e-6) |
| **Loss** | Cross-entropy |
| **Batch size** | 16 |
| **Epochs** | 20 per fold |
| **Gradient clipping** | Max norm 1.0 |
| **Cross-validation** | 5-fold stratified (seed=42) |
| **Metrics** | AUC-ROC (primary), accuracy (secondary) |
| **Checkpointing** | Save on best validation AUC per fold |
| **Resume** | Checkpoint stores model, optimizer, scheduler state + epoch + best AUC |

#### Data

| Attribute | Value |
|-----------|-------|
| **Dataset** | COVID-19 platelet scRNA-seq (47,808 cells, 12,395 genes) |
| **Binary task** | COVID (severe+fatal) vs Non-COVID (healthy+mild+moderate+recovered) — 29,219 vs 18,589 cells |
| **3-class task** | Healthy vs Mild/Moderate vs Severe/Fatal |
| **6-class task** | Full severity spectrum (healthy, mild, moderate, severe, fatal, recovered) |
| **Source** | 11 batches from 11 published studies (see Appendix E) |

#### Freezing Strategy

Base model weights are frozen. Only LoRA adapter parameters and the classification head are trained via backpropagation. This follows the scPEFT `freeze_parameters(model, DownstreamTasks.Identification)` pattern.

### Results: Binary Classification (COVID vs Non-COVID)

#### 5-Fold Cross-Validation Summary

| Fold | Best AUC | Best Acc | Best Epoch | Sessions | Status |
|------|----------|----------|:----------:|:--------:|--------|
| 1 | 0.8937 | 0.8110 | 20 | 3 | Done |
| 2 | 0.8884 | 0.8065 | 20 | 3 | Done |
| 3 | 0.8878 | 0.8032 | 17 | 3 | Done |
| 4 | **0.8990** | **0.8139** | 19 | 5 | Done |
| 5 | 0.8872 | 0.7849 | 19 | 9 | Done |
| **Mean** | **0.8912** | **0.8048** | — | — | — |
| **Std** | **±0.0045** | **±0.0080** | — | — | — |

*Mean/Std computed across all 5 folds using current best AUC values.*

#### Context: Comparison with Embedding Classification

| Method | Binary AUC | Notes |
|--------|:----------:|-------|
| STATE embedding + XGBoost | **0.944** | Frozen embeddings, no model weight updates |
| UCE embedding + XGBoost | 0.910 | Frozen embeddings |
| **scGPT LoRA PEFT** | **0.891** | **End-to-end fine-tuning (1.24% params), 5-fold CV** |
| scGPT embedding + XGBoost | 0.878 | Frozen 512-dim embeddings |
| Raw gene expression + XGBoost | 0.897 | No foundation model |
| Geneformer embedding + XGBoost | 0.862 | Frozen embeddings |
| TranscriptFormer embedding + XGBoost | 0.847 | Frozen embeddings |

**Key observations:**
1. scGPT LoRA PEFT (0.891) **improves over frozen scGPT embeddings** (0.878) by +0.013 AUC — fine-tuning helps scGPT learn task-specific features
2. PEFT still **underperforms STATE embedding classification** (0.944) by -0.053 AUC — STATE's pretrained embeddings already capture sufficient discriminative signal without fine-tuning
3. PEFT is **competitive with raw gene expression baseline** (0.897) — fine-tuning scGPT with 1.24% of parameters approaches the performance of a task-specific classifier on raw features
4. Platelets are OOD for all models (absent from pretraining data), making this a genuine generalization test

### Training Progress

#### Status (2026-02-12)

| Step | Status | Notes |
|------|--------|-------|
| Clone scPEFT repo | **Done** | `/bigdata/.../03_MODELS/scPEFT/` |
| Create scpeft_env | **Done** | `/bigdata/.../scpeft_env` (Python 3.9, PyTorch 2.1+cu118) |
| Install dependencies | **Done** | torchtext 0.16.0, numpy<2, scanpy, scikit-learn |
| Download scGPT checkpoint | **Done** | `best_model.pt`, `args.json`, `vocab.json` via AWS S3 |
| Fix env issues | **Done** | LD_LIBRARY_PATH, GeneVocab API (see Issues table) |
| Run scGPT PEFT binary | **Done** | AUC = 0.8912 ± 0.0045, Acc = 0.8048 ± 0.0080 (5-fold CV) |
| Run scGPT PEFT 3class | **Done** | AUC = 0.8984 ± 0.0030, Acc = 0.7470 ± 0.0107 (5-fold CV). Folds: 0.9013, 0.8994, 0.8926, 0.8988, 0.9000. |
| Run scGPT PEFT 6class | **Done** | AUC = 0.9020 ± 0.0049, Acc = 0.6324 ± 0.0158 (5-fold CV). Folds: 0.8995, 0.9034, 0.9109, 0.8985, 0.8976. Batch job 22489013 completed Feb 12. |
| Run Geneformer PEFT | Pending | `sbatch submit_geneformer_peft.sh` |

#### Detailed Training Log: Fold 1 (LoRA, Binary)

| Epoch | Train Loss | Val AUC | Val Acc | Notes |
|-------|-----------|---------|---------|-------|
| 1 | 0.6519 | 0.7858 | 0.7039 | |
| 2 | 0.5783 | 0.8222 | 0.7364 | |
| 3 | 0.5256 | 0.8395 | 0.7599 | |
| 4 | 0.4988 | 0.8601 | 0.7845 | |
| 5 | 0.4763 | 0.8704 | 0.7886 | |
| 6 | 0.4637 | 0.8740 | 0.7931 | |
| 7 | 0.4537 | 0.8833 | 0.7990 | Checkpoint saved. Session 1 killed by 2hr limit. |
| 8 | 0.4656 | 0.8772 | 0.7827 | Session 2 resumed from epoch 7 |
| 9 | 0.4338 | 0.8828 | 0.8006 | |
| 10 | 0.4272 | 0.8873 | 0.8029 | Checkpoint saved (new best) |
| 11 | 0.4239 | 0.8866 | 0.8046 | |
| 12 | 0.4187 | 0.8923 | 0.8085 | Checkpoint saved (new best) |
| 13 | 0.4140 | 0.8933 | 0.8036 | Checkpoint saved (new best) |
| 14 | 0.4145 | 0.8931 | 0.8108 | Session 2 killed by 2hr limit. |
| 14 | 0.4181 | 0.8900 | 0.8062 | Session 3 resumed from epoch 13 |
| 15 | 0.4101 | 0.8924 | 0.8106 | |
| 16 | 0.4112 | 0.8919 | 0.8126 | |
| 17 | 0.4100 | 0.8930 | 0.8105 | |
| 18 | 0.4083 | 0.8933 | 0.8098 | |
| 19 | 0.4072 | 0.8932 | 0.8100 | |
| 20 | 0.4116 | 0.8937 | 0.8110 | **Fold 1 complete. Best AUC: 0.8937** |

#### Training Characteristics

- **~15 min/epoch** on NVIDIA A100 GPU
- **Convergence**: AUC plateaus around epoch 13-15 for most folds; remaining epochs show marginal fluctuation
- **Total wall time**: ~25 hours for binary task (20 epochs × 5 folds), run across ~20 interactive `short_gpu` sessions (2-hour limit each)
- **Checkpoint overhead**: Folds with best AUC before epoch 20 re-run trailing epochs on each resume (~1 hr/session for folds 3+4)

### Issues Encountered & Fixes

| Issue | Root Cause | Fix |
|-------|-----------|-----|
| `No module named 'torchtext'` | scPEFT's scGPT imports torchtext internally | `pip install torchtext==0.16.0 --index-url .../cu118` |
| `CXXABI_1.3.15 not found` | System libstdc++ too old for conda ICU 78 | `export LD_LIBRARY_PATH=.../scpeft_env/lib:$LD_LIBRARY_PATH` |
| `No module named 'numpy'` | `PYTHONNOUSERSITE=1` blocks ~/.local packages | `pip install "numpy<2"` into scpeft_env |
| `No module named 'scanpy'` | Same as above | `pip install scanpy scikit-learn` into scpeft_env |
| `NumPy 2.0 incompatible` | PyTorch 2.1 compiled against NumPy 1.x | Pin `numpy<2` |
| `GeneVocab has no attribute 'stoi'` | scPEFT GeneVocab uses `__getitem__` directly | Use `vocab[gene]` and `gene in vocab` |

### Environment Activation (copy-paste ready)

```bash
source activate /bigdata/godziklab/shared/Xinru/302006/scpeft_env
export PYTHONPATH=/bigdata/godziklab/shared/Xinru/302006/03_MODELS/scPEFT:$PYTHONPATH
export LD_LIBRARY_PATH=/bigdata/godziklab/shared/Xinru/302006/scpeft_env/lib:$LD_LIBRARY_PATH
export PYTHONNOUSERSITE=1
```

### Next Steps

1. **Complete binary fold 5** (2 epochs remaining) — next session
2. **Run 3-class classification** — `python scgpt_peft_finetune.py --task 3class --peft LORA --epochs 20 --resume`
3. **Run 6-class classification** — `python scgpt_peft_finetune.py --task 6class --peft LORA --epochs 20 --resume`
4. **Run Geneformer PEFT** — `sbatch submit_geneformer_peft.sh`

**Results directory**: `/bigdata/.../05_RESULTS/scgpt_peft/`
**Checkpoint directory**: `/bigdata/.../05_RESULTS/checkpoints/scgpt_peft/`

---

## Summary: Top 10 Actionable Rules

1. **Fail loudly** -- raise errors, never generate fake fallback data
2. **Validate every intermediate output** -- dimensions, counts, NaNs, metadata
3. **Parallelize aggressively** -- joblib for CPU, SLURM arrays for cluster
4. **Choose faster algorithms first** -- solver choice alone can give 10x speedup
5. **Checkpoint everything** -- save state periodically for preemptible jobs
6. **Always include a simple baseline** -- prevents overstating results
7. **Use bootstrap CIs for significance** -- non-overlapping CIs = significant
8. **One environment per tool** -- avoid dependency hell
9. **Document failures** -- record what didn't work and why
10. **Use precise terminology** -- say exactly what your method does

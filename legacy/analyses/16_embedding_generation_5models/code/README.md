# Embedding Generation Scripts

This directory contains all embedding generation scripts for the 5 foundation models benchmarked in Project 302006.

## Scripts by Model

| Model | Script | Type | Description |
|-------|--------|------|-------------|
| **STATE** | `generate_state_embeddings.py` | Python | STATE SE-600M model embeddings via Inference API |
| **UCE** | `generate_uce_embeddings.py` | Python | UCE 33-layer model via `eval_single_anndata.py` |
| **scGPT** | `generate_scgpt_embeddings.sh` | Shell | SLURM batch job submission for chunked processing |
| **scGPT** | `run_scgpt_chunk.py` | Python | Core scGPT inference with batched processing (~10x faster) |
| **scGPT (Best Practices)** | `run_scgpt_best_practices.py` | Python | scGPT with HVG selection, normalization, log1p |
| **Geneformer** | `generate_geneformer_embeddings.py` | Python | Geneformer V1 with TranscriptomeTokenizer |
| **TranscriptFormer** | `generate_transcriptformer_embeddings.py` | Python | TF-Sapiens via CLI inference |
| **scPRINT-2** | `generate_scprint2_embeddings.py` | Python | scPRINT-2 small-v2 via Embedder API |

## scGPT Improved Design

The scGPT scripts have been refactored for better maintainability and performance:

```
generate_scgpt_embeddings.sh  # Shell wrapper: job submission, path validation
         |
         v
run_scgpt_chunk.py            # Python core: batched inference, embedding extraction
```

**Key improvements:**
1. **Separated concerns** - Shell script handles SLURM, Python handles inference
2. **Batched processing** - ~10x faster than cell-by-cell (batch_size=64)
3. **Parameterized paths** - Pass base directory as argument for portability
4. **Error handling** - `set -e` for fail-fast behavior

**Usage:**
```bash
# Default path
./generate_scgpt_embeddings.sh

# Custom path
./generate_scgpt_embeddings.sh /path/to/project
```

## Supporting Scripts

| Script | Description |
|--------|-------------|
| `submit_scgpt_best_practices.sh` | SLURM submission script for scGPT best practices |
| `submit_uce_embeddings.sh` | SLURM submission script for UCE |
| `run_scgpt_bp_benchmark.py` | Benchmark comparing original vs best practices scGPT |

## UCE Script Design

The UCE scripts follow the same separation of concerns as scGPT:

```
submit_uce_embeddings.sh      # SLURM submission wrapper
         |
         v
generate_uce_embeddings.py    # Python core: validation, subprocess, error handling
         |
         v
eval_single_anndata.py        # UCE repository script (external)
```

**Usage:**
```bash
# Single dataset
python generate_uce_embeddings.py --base-dir /path/to/project --dataset covid

# All datasets
python generate_uce_embeddings.py --base-dir /path/to/project --dataset all

# Custom file
python generate_uce_embeddings.py --input data.h5ad --output embeddings/ --uce-dir /path/to/UCE
```

## Model Requirements

### STATE
- Package: `pip install arc-state`
- Model: SE-600M checkpoint (epoch4 or epoch15)
- Input: Raw counts with gene symbols

### UCE
- Repository: [UCE GitHub](https://github.com/snap-stanford/UCE)
- Model: 33-layer transformer (650M parameters)
- Input: h5ad with gene symbols (not Ensembl IDs)

### scGPT
- Container: `scgpt_0.1.7.sif`
- Model: scGPT_human
- Input: Raw counts with gene symbols
- Best practices: HVG selection, normalize_total, log1p, binning

### Geneformer
- Package: `pip install geneformer` or HuggingFace
- Model: `ctheodoris/Geneformer`
- Input: loom file with `ensembl_id` and `n_counts`

### TranscriptFormer
- Package: `pip install transcriptformer`
- Model: TF-Sapiens (368M parameters)
- Input: h5ad with `ensembl_id` column

### scPRINT-2
- Package: `uv pip install scprint2` (requires lamin: `lamin init --storage <path> --name <name> --modules bionty`)
- Setup: `scprint2 easy_setup` (populates ontologies)
- Model: small-v2 checkpoint from HuggingFace (`hf download jkobject/scPRINT small-v2.ckpt`)
- Input: Raw counts h5ad with Ensembl IDs or HUGO symbols, `obs['organism_ontology_term_id']`
- Supports 16 species including human and mouse
- Pre-trained on 350M+ cells

## Output Format

All models output embeddings in AnnData format:

| Model | Embedding Key | Dimension |
|-------|---------------|-----------|
| STATE | `X_state` | 2,058 |
| UCE | `X_uce` | 1,280 |
| scGPT | `X_scgpt` | 512 |
| Geneformer | `X_geneformer` | 1,152 |
| TranscriptFormer | `X_transcriptformer` | 2,048 |
| scPRINT-2 | `X_scprint2` | 256 (d_model of small-v2) |

---

**Date:** January 19, 2026
**Project:** 302006 - COVID Disease State Prediction

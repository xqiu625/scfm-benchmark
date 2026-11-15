## TranscriptFormer Implementation for Project 302006

### Overview

TranscriptFormer is a breakthrough foundation model that **outperforms UCE, scGPT, and Geneformer by up to 7%** in disease identification tasks. It's trained on 112 million cells across 12 species and excels at:
- Zero-shot disease identification (F1: 0.859 vs UCE: 0.805)
- Cross-species analysis
- Generative modeling

### Key Advantages

1. **Superior Performance**: 7% better than existing models in benchmarks
2. **No Manual Preprocessing**: Automatic normalization and feature selection
3. **Multi-species Training**: Robust cross-species generalization
4. **Large Scale**: Handles 200K+ cells efficiently

### Requirements

- **Python**: 3.11+ (strict requirement)
- **PyTorch**: ≤2.5.1 (versions ≥2.6.0 cause pickle errors)
- **GPU**: 16GB VRAM minimum (A100 40GB recommended)
- **Storage**: ~5GB for TF-Sapiens model

### Installation Steps

```bash
# 1. Run installation script
bash scripts/foundation_models/TranscriptFormer/install_transcriptformer.sh

# This will:
#   - Check Python version (must be 3.11+)
#   - Create virtual environment
#   - Install PyTorch 2.5.1
#   - Install TranscriptFormer
#   - Download TF-Sapiens model weights
```

### Data Preparation

TranscriptFormer requires:
- **Raw counts** (no normalization)
- **Ensembl gene IDs** in `adata.var['ensembl_id']`
- **H5AD format**

```bash
# 2. Prepare data
source /bigdata/godziklab/shared/Xinru/transcriptformer_env/bin/activate
python scripts/foundation_models/TranscriptFormer/prepare_data_for_transcriptformer.py

# This will:
#   - Load raw count data
#   - Add Ensembl IDs via mygene.info
#   - Perform basic QC (200-5000 genes/cell, <20% MT)
#   - Save to 01_DATA/transcriptformer_ready/
```

### Generate Embeddings

```bash
# 3. Generate embeddings on GPU
srun -p short_gpu --gres=gpu:1 --mem=256gb -c 4 --pty bash -l
source /bigdata/godziklab/shared/Xinru/transcriptformer_env/bin/activate
cd /bigdata/godziklab/shared/Xinru/302006

# Run with auto-detected batch size
python scripts/foundation_models/TranscriptFormer/generate_transcriptformer_embeddings.py

# Or specify batch size manually
python scripts/foundation_models/TranscriptFormer/generate_transcriptformer_embeddings.py --batch-size 16

# Or process single dataset
python scripts/foundation_models/TranscriptFormer/generate_transcriptformer_embeddings.py --dataset covid
```

### Batch Size Guidelines

Based on GPU memory:
- **A100 80GB**: batch_size=32-64
- **A100 40GB**: batch_size=16-32
- **V100 32GB**: batch_size=8-16
- **T4 16GB**: batch_size=2-4

### Expected Performance

Processing times (approximate):
- **COVID (47,977 cells)**: ~30-60 minutes on A100 40GB
- **HEK293T (88,434 cells)**: ~60-90 minutes
- **HCT116 (89,738 cells)**: ~60-90 minutes

### Output

Embeddings saved to: `02_EMBEDDINGS/transcriptformer/`

Files:
- `covid_transcriptformer.h5ad`
- `hek293t_transcriptformer.h5ad`
- `hct116_transcriptformer.h5ad`

Each contains `adata.obsm['X_transcriptformer']` with embedding dimensions (typically 512 or 1280).

### Troubleshooting

**PyTorch version error**:
```bash
pip install torch==2.5.1 --index-url https://download.pytorch.org/whl/cu121
```

**Python version error**:
```bash
# Load Python 3.11+ module on cluster
module load python/3.11
```

**Out of memory**:
```bash
# Reduce batch size
python generate_transcriptformer_embeddings.py --batch-size 4
```

**Missing Ensembl IDs**:
```bash
# Install mygene for automatic mapping
pip install mygene
```

### Next Steps

After generating TranscriptFormer embeddings:

1. **Merge embeddings** with UCE and scGPT
2. **Three-model comparison** analysis
3. **Evaluate performance** on COVID severity prediction

### Comparison with Other Models

| Model | F1 Score | Training Data | Parameters |
|-------|----------|---------------|------------|
| **TranscriptFormer** | **0.859** | 112M cells, 12 species | 368M-542M |
| UCE | 0.805 | 36M cells, 3 species | 1280-dim |
| scGPT | 0.798 | 33M cells, human | 512-dim |
| Geneformer | 0.797 | 30M cells, human | 256-dim |

### References

- TranscriptFormer paper: Nature (2025)
- GitHub: https://github.com/czbiohub-sf/transcriptformer
- Documentation: https://transcriptformer.readthedocs.io/

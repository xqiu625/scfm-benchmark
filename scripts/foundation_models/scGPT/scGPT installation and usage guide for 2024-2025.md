# scGPT installation and usage guide for 2024-2025

**scGPT**, the foundation model for single-cell genomics trained on over 33 million cells, has become significantly more accessible in 2024-2025 with **Docker containers providing the easiest installation path**, PyTorch 2.1.2 support, and multiple cloud-based alternatives that eliminate local setup entirely. The fastest way to get started is through the **Virtual Cells Platform Google Colab notebook** at virtualcellmodels.cziscience.com, which provides a complete working environment with GPU support in under 5 minutes.

## Quick start: Three paths to scGPT

The landscape of scGPT installation has evolved to accommodate different user needs and technical expertise levels. For immediate access without any installation, the **Virtual Cells Platform** offers a comprehensive Google Colab environment with pre-configured dependencies and automated model downloads. Researchers seeking the most reliable local installation should use **Docker** (`docker pull xueerchen/scgpt:0.1.7`), which eliminates dependency conflicts entirely. For traditional Python environments, a **conda setup with PyPI installation** provides the most flexibility, though it requires careful attention to PyTorch 2.1.2 and CUDA 11.7 compatibility.

### Docker installation eliminates dependency headaches

Docker has emerged as the preferred installation method for scGPT, solving the complex dependency management that previously plagued users. The official and community-maintained Docker images include all necessary components pre-configured:

```bash
# Pull and run the most stable Docker image
docker pull xueerchen/scgpt:0.1.7
docker run --gpus all --rm -it xueerchen/scgpt:0.1.7 bash

# Alternative community images with different configurations
docker pull kzkedzierska/scgpt:0.2  # Updated dependencies
docker pull wangbaonan/scgpt:latest  # Latest features
```

This approach bypasses the common R installation errors, flash-attention compilation issues, and CUDA version mismatches that frequently occur with manual installations. The container includes Python 3.10, PyTorch 2.1.2, and all bioinformatics dependencies pre-configured.

### Manual installation requires specific version alignment

For users requiring local installation, success depends on maintaining strict version compatibility. The **critical requirement** is PyTorch 2.1.2 (not 2.0 or 1.13 as in older tutorials), which represents a major update from earlier documentation. The installation sequence matters significantly:

```bash
# Create clean environment with Python 3.10
conda create -n scgpt python=3.10
conda activate scgpt

# Install PyTorch 2.1.2 first (critical version)
pip install torch==2.1.2 torchvision==0.16.2

# Install scGPT with optional optimizations
pip install scgpt
pip install "flash-attn<1.0.5"  # Optional but recommended
pip install wandb  # For experiment tracking

# Verify successful installation
python -c "import scgpt; print('scGPT version:', scgpt.__version__)"
```

The most common installation failure occurs when R is not present on the system, as scGPT depends on rpy2 for certain statistical functions. Installing R before scGPT resolves this: `sudo apt install r-base r-base-dev` on Ubuntu or `conda install -c conda-forge r-base` in conda environments.

## Working code for immediate productivity

The power of scGPT lies in its ability to generate biologically meaningful cell embeddings with minimal code. The simplest working example for cell annotation demonstrates the high-level API:

```python
import scgpt as scg
import scanpy as sc
from pathlib import Path

# Load your single-cell data
adata = sc.read_h5ad("your_data.h5ad")

# Download pre-trained model (one-time setup)
import gdown
folder_id = '1oWh_-ZRdhtoGQ2Fw24HP41FgLoomVo-y'
gdown.download_folder(f"https://drive.google.com/drive/folders/{folder_id}", 
                      output="./scGPT_human", quiet=False)

# Generate 512-dimensional embeddings
model_dir = Path("./scGPT_human")
embed_adata = scg.tasks.embed_data(
    adata,
    model_dir,
    gene_col="feature_name",  # or "index" for standard data
    batch_size=64,
    return_new_adata=True,
)

# Use embeddings for downstream analysis
adata.obsm["X_scgpt"] = embed_adata.X
sc.pp.neighbors(adata, use_rep="X_scgpt")
sc.tl.umap(adata)
sc.pl.umap(adata, color="cell_type", frameon=False)
```

This code handles the complete workflow from data loading to visualization, automatically managing gene vocabulary matching and value binning internally. For **CZ CELLxGENE data**, use `gene_col="feature_name"`; for other sources, `gene_col="index"` typically works.

### Batch integration requires minimal configuration

Batch effect correction, one of scGPT's strongest capabilities, requires only slight modifications to enable domain adaptation:

```python
from scgpt.preprocess import Preprocessor

# Configure preprocessing for batch integration
preprocessor = Preprocessor(
    use_key="X",
    filter_gene_by_counts=3,
    normalize_total=1e4,
    log1p=True,
    subset_hvg=1200,  # Use 1200 highly variable genes
    hvg_flavor="seurat_v3",
    binning=51,  # Value binning for transformer input
    result_binned_key="X_binned",
)

# Apply preprocessing with batch information
adata.obs["batch"] = adata.obs["batch"].astype("category")
preprocessor(adata, batch_key="batch")

# Fine-tune for integration (simplified)
config = dict(
    load_model="./scGPT_human",
    GEPC=True,  # Gene expression prediction
    ecs_thres=0.8,  # Elastic cell similarity
    dab_weight=1.0,  # Domain adaptation weight
    mask_ratio=0.4,
    epochs=15,
    batch_size=64,
    lr=1e-4,
)
```

## Cloud platforms eliminate installation entirely

The emergence of cloud-based scGPT platforms has dramatically lowered the barrier to entry. **Superbio.ai** offers browser-based, no-code interfaces for cell annotation, reference mapping, and gene regulatory network inference at app.superbio.ai. Their freemium model provides immediate access with compute credits, making it ideal for researchers without programming expertise.

For code-based cloud solutions, the **Virtual Cells Platform** (virtualcellmodels.cziscience.com) provides the most comprehensive Google Colab environment. Their notebook includes pre-configured PyTorch 2.1.2, automated model downloads, and direct integration with CZ CELLxGENE Census data containing over 50 million cells. This platform has become the recommended starting point for new users, requiring only a Google account and browser access.

**BioTuring Platform** offers professional-grade notebook environments at colab.bioturing.com with curated scGPT workflows and GPU support. While subscription-based, it provides enterprise features including data privacy controls and dedicated compute resources. For large-scale institutional deployments, **AWS SageMaker** integrates with AWS HealthOmics, offering scalable genomics workloads with pay-per-use pricing on instances ranging from ml.t2.medium for development to ml.p4d.24xlarge for production workloads.

## Hardware requirements shape deployment strategies

scGPT's transformer architecture demands significant computational resources, with **GPU memory being the primary bottleneck**. The model requires a minimum of 8GB VRAM for basic inference, but 16-24GB is recommended for fine-tuning tasks. NVIDIA GPUs with Ampere architecture or newer (RTX 3090, A100, H100) provide optimal performance through tensor core acceleration.

The critical technical requirement is **CUDA 11.7 compatibility** with flash-attention support for efficient attention computation. System RAM requirements scale with dataset size, with 32GB recommended for preprocessing large datasets and handling compilation of flash-attention. For CPU-only deployment, while possible, expect 10-50x slower inference compared to GPU execution.

Memory optimization strategies can extend hardware capabilities. Using mixed precision (FP16) reduces memory usage by 50% with minimal accuracy loss:

```python
import torch

# Enable mixed precision for memory efficiency
model = model.half()  # Convert to FP16
torch.backends.cudnn.benchmark = True

# Monitor GPU memory usage
print(f"Allocated: {torch.cuda.memory_allocated()/1024**3:.2f} GB")
```

## Common issues have clear solutions

Installation failures follow predictable patterns with established solutions. The **most frequent error** involves missing R installation, manifesting as "rpy2 in API mode cannot be built without R". Installing R before scGPT consistently resolves this issue.

Flash-attention compilation failures, the second most common problem, typically stem from CUDA toolkit absence or version mismatches. The solution involves ensuring CUDA 11.7 is properly installed with `nvcc` available in PATH:

```bash
# Verify CUDA installation
nvcc --version  # Should show 11.7.x

# Install with relaxed requirements if needed
pip install flash-attn==1.0.4 --no-build-isolation

# Alternative: Skip flash-attention (slower but works)
pip install scgpt --no-deps
pip install torch scanpy anndata numpy pandas
```

For **orbax version conflicts**, common with certain dependency combinations, pinning versions resolves the issue: `pip install scgpt "flash-attn<1.0.5" "orbax<0.1.8"`.

## Data preprocessing determines model performance

scGPT's performance critically depends on proper data preprocessing, particularly gene vocabulary matching and value binning. The standard preprocessing pipeline handles most single-cell datasets effectively:

```python
def preprocess_for_scgpt(adata, vocab_file="./scGPT_human/vocab.json"):
    """Standard preprocessing pipeline for scGPT"""
    
    # Filter low-quality cells and genes
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    
    # Select highly variable genes
    sc.pp.highly_variable_genes(adata, n_top_genes=3000, flavor='seurat_v3')
    adata_hvg = adata[:, adata.var['highly_variable']].copy()
    
    # Match genes with model vocabulary
    from scgpt.tokenizer.gene_tokenizer import GeneVocab
    vocab = GeneVocab.from_file(vocab_file)
    
    adata_hvg.var["in_vocab"] = [gene in vocab for gene in adata_hvg.var_names]
    adata_final = adata_hvg[:, adata_hvg.var["in_vocab"]].copy()
    
    print(f"Using {adata_final.shape[1]} genes from vocabulary")
    return adata_final
```

The binning step, unique to scGPT, discretizes continuous expression values into 51 bins, enabling the transformer to process gene expression as tokens. This approach maintains biological signal while enabling attention mechanisms to capture cell-state relationships.

## Fine-tuning unlocks domain-specific performance

While pre-trained scGPT models achieve strong zero-shot performance, fine-tuning on domain-specific data can boost cell-type annotation accuracy from ~85% to **99.5% F1-score**, as demonstrated in recent retinal cell studies. The fine-tuning process requires careful configuration of training objectives:

```python
# Fine-tuning configuration for cell annotation
config = dict(
    seed=42,
    load_model="./scGPT_human",
    CLS=True,  # Classification objective
    MVC=True,  # Masked value prediction
    mask_ratio=0.0,  # No masking for classification
    epochs=10,
    lr=1e-4,
    batch_size=32,
    dropout=0.2,
)

# Key hyperparameters for different tasks:
# - Cell annotation: CLS=True, mask_ratio=0.0
# - Batch integration: dab_weight=1.0, ecs_thres=0.8
# - Perturbation: MVC=True, mask_ratio=0.4
```

Popular datasets for fine-tuning include **PBMC datasets** (3K/8K/10K cells) available through scvi-tools, the **Multiple Sclerosis dataset** from scGPT tutorials demonstrating transfer learning, and the **Human Retina dataset** achieving state-of-the-art classification performance. These datasets, available as preprocessed .h5ad files, provide benchmarks for comparing model performance.

## Jupyter notebooks accelerate learning

The scGPT ecosystem includes extensive Jupyter notebook resources that demonstrate real-world applications. The official GitHub repository (github.com/bowang-lab/scGPT) contains tutorials covering zero-shot cell embedding, reference mapping to 33+ million cells, and multi-batch integration workflows.

For immediate experimentation, Google Colab notebooks provide the fastest path to results. The Virtual Cells Platform notebook walks through the complete workflow from data loading to UMAP visualization in approximately 20 cells of code. These notebooks handle environment setup automatically, including GPU configuration and dependency installation.

Creating custom notebooks benefits from understanding the execution flow:

```python
# Standard notebook setup
import warnings
warnings.filterwarnings("ignore")

# GPU configuration
import torch
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print(f"Using device: {device}")

# Install dependencies (Colab)
!pip install -q scgpt scanpy gdown

# Download models
import gdown
model_url = "https://drive.google.com/drive/folders/[MODEL_ID]"
gdown.download_folder(model_url, quiet=False)
```

## Performance optimization extends capabilities

Large datasets exceeding available memory require chunked processing strategies. Processing 1 million cells becomes feasible through batch iteration:

```python
def process_large_dataset(adata, model_dir, chunk_size=1000):
    """Process large datasets in memory-efficient chunks"""
    
    embeddings = []
    n_chunks = len(adata) // chunk_size + 1
    
    for i in range(n_chunks):
        start_idx = i * chunk_size
        end_idx = min((i + 1) * chunk_size, len(adata))
        
        chunk = adata[start_idx:end_idx].copy()
        chunk_embed = scg.tasks.embed_data(
            chunk, model_dir, 
            batch_size=64,
            return_new_adata=True
        )
        embeddings.append(chunk_embed.X)
        
        # Clear GPU memory between chunks
        torch.cuda.empty_cache()
    
    return np.vstack(embeddings)
```

Additional optimizations include gradient checkpointing for fine-tuning (trading compute for memory), using flash-attention for 2-4x speedup on compatible GPUs, and employing mixed precision training with automatic mixed precision (AMP) for 2x memory savings.

## The ecosystem continues rapid evolution

The scGPT ecosystem in 2024-2025 represents a mature, production-ready platform for single-cell analysis. The combination of Docker containers for reproducible environments, cloud platforms for zero-installation access, and comprehensive documentation has transformed scGPT from a research prototype into an accessible tool for the broader genomics community.

Future developments include **HuggingFace integration** (currently in beta) promising simplified model distribution, expanded pre-trained models for specific organs and disease states, and integration with spatial transcriptomics workflows. The active development community, evidenced by regular GitHub updates and expanding cloud platform support, ensures continued improvement in both functionality and accessibility.

For researchers beginning with scGPT today, the recommended path starts with the Virtual Cells Platform for initial exploration, progresses to Docker containers for reproducible local analysis, and culminates in custom fine-tuning for domain-specific applications. This progression, supported by the code examples and troubleshooting solutions provided above, enables researchers to leverage state-of-the-art foundation models for single-cell genomics regardless of their computational expertise level.
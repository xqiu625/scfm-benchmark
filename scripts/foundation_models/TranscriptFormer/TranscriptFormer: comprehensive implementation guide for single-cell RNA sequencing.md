# TranscriptFormer: A comprehensive scRNA-seq implementation guide

TranscriptFormer is a **generative foundation model** from the Chan Zuckerberg Initiative that enables cross-species single-cell RNA sequencing analysis without task-specific fine-tuning. Trained on **112 million cells** spanning 12 species and **1.53 billion years of evolution**, it represents the first multi-species generative model for single-cell transcriptomics. The tool excels at zero-shot cell type classification, cross-species annotation transfer, and gene regulatory inference—outperforming alternatives like scGPT and Geneformer on key benchmarks. **No official Docker or Singularity containers currently exist**, requiring native Python installation or custom container creation.

## How TranscriptFormer revolutionizes single-cell analysis

TranscriptFormer's architecture centers on a **12-layer transformer encoder** with 16 attention heads and a 2048-dimensional hidden layer, totaling approximately 300–540 million trainable parameters depending on variant. What distinguishes it from competitors is its expression-aware multi-head self-attention mechanism, which incorporates expression counts as a log-count bias term in the attention matrix rather than using token duplication like other models.

The model leverages **ESM-2 protein embeddings** as gene representations, enabling generalization to species not seen during training. Three pretrained variants serve different use cases:

| Variant | Training scope | Parameters | Best for |
|---------|---------------|------------|----------|
| **TF-Sapiens** | 57M human cells | 368M | Human-only studies |
| **TF-Exemplar** | 110M cells, 5 model organisms | 542M | Common lab species |
| **TF-Metazoa** | 112M cells, 12 species | 444M | Cross-species & rare organisms |

The species coverage spans vertebrates (human, mouse, zebrafish, chicken), invertebrates (C. elegans, Drosophila, sea urchin), yeast, and even the malaria parasite—providing unprecedented evolutionary breadth for transfer learning.

## Installation requires Python 3.11 and GPU hardware

TranscriptFormer demands **Python 3.11 or higher** and substantial GPU resources. An NVIDIA A100 with 40GB VRAM is recommended, though 16GB GPUs work with reduced batch sizes (1–4).

```bash
# Create environment using uv package manager
uv venv --python=3.11
source .venv/bin/activate

# Install from PyPI
uv pip install transcriptformer

# Download model weights (choose one)
transcriptformer download tf-sapiens      # ~4GB, human-only
transcriptformer download tf-exemplar     # ~5GB, multi-species
transcriptformer download tf-metazoa      # Largest, broadest coverage
```

For development installation from source:

```bash
git clone https://github.com/czi-ai/transcriptformer.git
cd transcriptformer
uv pip install -e .
```

**Critical dependency note**: PyTorch must be version **2.5.1 or earlier**—version 2.6.0+ causes pickle serialization errors. Core dependencies include PyTorch Lightning, anndata, scanpy, hydra-core, and h5py.

## Running inference follows a straightforward workflow

TranscriptFormer accepts **H5AD files** (AnnData format) with raw, unnormalized counts. The `var` dataframe must contain an `ensembl_id` column with Ensembl gene identifiers. No manual preprocessing is required—the model handles all normalization and transformation internally.

```bash
# Basic inference on human data
transcriptformer inference \
    --checkpoint-path ./checkpoints/tf_sapiens \
    --data-file your_data.h5ad \
    --output-path ./results \
    --batch-size 8

# Cross-species analysis (e.g., mouse with human model)
transcriptformer inference \
    --checkpoint-path ./checkpoints/tf_sapiens \
    --data-file mouse_data.h5ad \
    --pretrained-embedding ./checkpoints/all_embeddings/mus_musculus_gene.h5 \
    --batch-size 4

# Extract contextual gene embeddings for regulatory analysis
transcriptformer inference \
    --checkpoint-path ./checkpoints/tf_exemplar \
    --data-file your_data.h5ad \
    --emb-type cge \
    --num-gpus 4
```

Output files are saved as AnnData objects containing cell embeddings in `obsm['embeddings']`, contextual gene embeddings in `uns['cge_embeddings']` (when using `--emb-type cge`), and log-likelihood scores in `uns['llh']` for quality assessment.

For Python API usage:

```python
import anndata as ad
from transcriptformer.model.inference import run_inference
from omegaconf import OmegaConf

cfg = OmegaConf.load("./conf/inference_config.yaml")
cfg.model.checkpoint_path = "./checkpoints/tf_exemplar"

adata = ad.read_h5ad("your_data.h5ad")
adata_output = run_inference(cfg, data_files=[adata])
embeddings = adata_output.obsm["embeddings"]
```

## No official containers exist—custom builds required

TranscriptFormer currently lacks Docker images, Singularity containers, or presence on BioContainers, DockerHub, or GitHub Container Registry. Users requiring containerized deployment must build custom images.

**Custom Dockerfile template**:

```dockerfile
FROM nvidia/cuda:11.8-runtime-ubuntu22.04

RUN apt-get update && apt-get install -y \
    python3.11 python3-pip && \
    pip3 install transcriptformer

# Optional: pre-download model weights
RUN transcriptformer download tf-sapiens

WORKDIR /workspace
ENTRYPOINT ["transcriptformer"]
```

**Convert to Singularity for HPC clusters**:

```bash
singularity build transcriptformer.sif docker://your-registry/transcriptformer:latest
```

For native HPC deployment without containers:

```bash
module load python/3.11 cuda/11.8
python -m venv tf_env && source tf_env/bin/activate
pip install transcriptformer
transcriptformer inference --num-gpus 4 --data-file data.h5ad
```

## GitHub repository and official resources

The official repository at **https://github.com/czi-ai/transcriptformer** (MIT License) contains source code, documentation, and test data. Key resources include:

- **Model cards**: https://virtualcellmodels.cziscience.com/model/transcriptformer
- **Interactive tutorial**: https://virtualcellmodels.cziscience.com/quickstart/transcriptformer-quickstart (demonstrates macaque-to-marmoset cross-species annotation in ~15 minutes)
- **Preprint**: Pearce et al. (2025), "A Cross-Species Generative Cell Atlas," bioRxiv DOI: 10.1101/2025.04.25.650731
- **Third-party integration**: The Helical framework provides standardized comparisons with Geneformer

## Best practices for effective implementation

**Input data preparation** requires minimal effort—TranscriptFormer expects raw counts with Ensembl gene IDs. Count values are internally clipped at 30 by default, matching the training configuration. For datasets exceeding 100,000 cells, enable multi-GPU inference with `--num-gpus -1` (uses all available GPUs).

**Model selection** depends on your species: use TF-Sapiens for human-only analyses (fastest inference), TF-Exemplar for common model organisms (human, mouse, zebrafish, C. elegans, Drosophila), and TF-Metazoa for unusual species or maximum cross-species generalization.

**Memory management**: If encountering out-of-memory errors, reduce batch size to 1–4 for 16GB GPUs, or add `--oom-dataloader --n-data-workers 4` flags. The `--precision 16-mixed` setting (default) optimizes memory usage.

**For novel species** not in the training vocabulary, generate custom ESM-2 embeddings using scripts in the repository, then pass them via `--pretrained-embedding`. TranscriptFormer maintains strong performance even for species **685 million years** divergent from training organisms.

## Conclusion

TranscriptFormer represents a significant advancement in foundation models for single-cell analysis, offering state-of-the-art zero-shot performance and unprecedented cross-species capabilities. The straightforward CLI interface and minimal preprocessing requirements lower adoption barriers, though the absence of official containers creates friction for HPC and cloud deployment. For teams prioritizing cross-species annotation transfer or rare cell type identification, TranscriptFormer's generative architecture and evolutionary training breadth provide clear advantages over scGPT and Geneformer. The tool remains in active development at CZI, with the bioRxiv preprint awaiting formal peer review—users should monitor the GitHub repository for updates and community validation studies.
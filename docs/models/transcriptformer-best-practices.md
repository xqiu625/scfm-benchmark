# TranscriptFormer: Comprehensive Research Document

## A Cross-Species Generative Cell Atlas Across 1.5 Billion Years of Evolution

---

## Executive Summary

TranscriptFormer is a family of generative foundation models for single-cell transcriptomics developed by the Chan Zuckerberg Initiative (CZI). It represents a cross-species generative cell atlas trained on up to 112 million cells spanning 1.53 billion years of evolution across 12 species. Unlike previous single-cell foundation models that focus primarily on embedding generation, TranscriptFormer is the **first generative, multi-species model** for single-cell transcriptomics, enabling it to function as a "virtual instrument" for probing cellular biology through prompting and hypothesis generation.

**Key Innovations:**
- **Generative autoregressive architecture** jointly modeling genes and expression counts
- **Cross-species capability** spanning 1.53 billion years of evolutionary distance
- **ESM-2 protein embeddings** for species-agnostic gene representation
- **Expression-aware multi-head self-attention** with novel count likelihood
- **Zero-shot performance** on out-of-distribution species and disease states

**Citation:** Pearce, J. D., Simmonds, S. E., Mahmoudabadi, G., Krishnan, L., Palla, G., Istrate, A.-M., Tarashansky, A., Nelson, B., Valenzuela, O., Li, D., Quake, S. R., & Karaletsos, T. (2025). A Cross-Species Generative Cell Atlas Across 1.5 Billion Years of Evolution: The TranscriptFormer Single-cell Model. bioRxiv. DOI: 10.1101/2025.04.25.650731

---

## Table of Contents

1. [Model Overview](#model-overview)
2. [Architecture](#architecture)
3. [Model Variants](#model-variants)
4. [Training Data](#training-data)
5. [Training Methodology](#training-methodology)
6. [Capabilities & Applications](#capabilities--applications)
7. [Benchmark Performance](#benchmark-performance)
8. [Cross-Species Generalization](#cross-species-generalization)
9. [Generative Prompting](#generative-prompting)
10. [Comparison with Other Models](#comparison-with-other-models)
11. [Installation & Usage](#installation--usage)
12. [Limitations](#limitations)
13. [Ecosystem & Resources](#ecosystem--resources)
14. [Key Takeaways](#key-takeaways)

---

## Model Overview

### What is TranscriptFormer?

TranscriptFormer is a generative autoregressive foundation model that jointly models gene identities and their expression levels as "cell sentences." It captures the probability of each gene-expression tuple conditioned on preceding tokens, enabling it to represent cellular states and facilitate "virtual instruments" for biological inquiry.

### Key Differentiators

| Feature | TranscriptFormer | Traditional Models |
|---------|------------------|-------------------|
| **Paradigm** | Generative autoregressive | Masked language modeling (BERT-like) |
| **Cross-species** | Native (12 species, 1.53B years) | Requires ortholog conversion |
| **Gene representation** | ESM-2 protein embeddings | Rank encoding or learned tokens |
| **Expression modeling** | Explicit count likelihood | Implicit or binned |
| **Prompting** | Yes (for gene interactions) | Limited/None |

### Developers

- **Organization:** Chan Zuckerberg Initiative (CZI)
- **Authors:** James D. Pearce, Sara E. Simmonds*, Gita Mahmoudabadi*, Lakshmi Krishnan*, Giovanni Palla, Ana-Maria Istrate, Alexander Tarashansky, Benjamin Nelson, Omar Valenzuela, Donghui Li, Stephen R. Quake, Theofanis Karaletsos
- **Publication Date:** April 29, 2025 (bioRxiv preprint, updated October 2025)

---

## Architecture

### Core Architecture

TranscriptFormer employs a **transformer encoder** with several novel architectural components:

**Specifications:**
- **Layers:** 12 transformer encoder layers
- **Attention Heads:** 16 per layer
- **Hidden Dimension:** 2048
- **Core Parameters:** ~300 million trainable parameters in transformer layers
- **Total Parameters (TF-Metazoa):** 444 million trainable + 633 million non-trainable (frozen ESM-2 embeddings)

### Input Representation: "Cell Sentences"

Each cell is represented as a sequence of gene-expression tuples:

```
Cell = [(gene₁, count₁), (gene₂, count₂), ..., (geneₙ, countₙ)]
```

**Key Input Components:**
1. **Gene Embeddings:** Derived from ESM-2 (Evolutionary Scale Model 2) protein embeddings for each gene
2. **Assay Token:** Captures sequencing platform metadata
3. **Expression Counts:** Raw transcript counts (clipped at 30 by default)

### Novel Architectural Features

#### 1. Expression-Aware Multi-Head Self-Attention

Expression counts are introduced as a **log-count bias term** in the attention matrix, avoiding explicit token duplication:

```
Attention(Q, K, V) = softmax(QKᵀ/√d + log_count_bias) × V
```

This design:
- Efficiently incorporates expression information
- Avoids vocabulary explosion from count tokenization
- Preserves count-level variability

#### 2. Coupled Gene and Transcript Heads

TranscriptFormer uses a novel coupling between:
- **Gene Head:** Predicts gene identity at each position
- **Transcript/Count Head:** Predicts expression count for each gene

#### 3. Causal Masking

Autoregressive generation is enabled through causal masking, where each position can only attend to previous positions, enabling left-to-right generation.

#### 4. Count Likelihood (Zero-Truncated Poisson)

Expression counts are modeled using a **zero-truncated Poisson negative log-likelihood (NLL)**:

```
L = -log P(c_j | λ_j) = λ_j + log(1 - e^{-λ_j}) - c_j log(λ_j) + log(c_j!)
```

Where:
- λ_j = predicted rate parameter at position j
- c_j = observed ground truth count

This formulation accounts for the discrete, non-negative nature of gene expression and the exclusion of zero counts from input sequences.

### ESM-2 Integration

TranscriptFormer maps genes into a shared, species-agnostic embedding space using **ESM-2 amino acid embeddings** for coding genes:

**Why ESM-2?**
- Captures evolutionary and functional similarity between genes
- Allows orthologous genes across species to be represented similarly
- Enables cross-species analysis without explicit ortholog mapping
- Pre-computed embeddings available for 24+ species

**ESM-2 Model Used:** esm2_t33_650M_UR50D (33 layers, 650M parameters, 1024-dim embeddings)

---

## Model Variants

TranscriptFormer comes in three variants with identical core architecture but different training data:

### TF-Metazoa (Largest & Most Diverse)

| Attribute | Value |
|-----------|-------|
| **Training Cells** | 112 million |
| **Species** | 12 |
| **Trainable Parameters** | 444 million |
| **Non-trainable Parameters** | 633 million (frozen ESM-2) |
| **Vocabulary Size** | 247,388 genes |
| **Evolutionary Span** | 1.53 billion years |

**Species Coverage:**
- **Vertebrates (6):** Human (*Homo sapiens*), Mouse (*Mus musculus*), Rabbit (*Oryctolagus cuniculus*), Chicken (*Gallus gallus*), African clawed frog (*Xenopus laevis*), Zebrafish (*Danio rerio*)
- **Invertebrates (4):** Sea urchin (*Lytechinus variegatus*), *C. elegans*, Fruit fly (*Drosophila melanogaster*), Freshwater sponge (*Spongilla lacustris*)
- **Fungus (1):** Yeast (*Saccharomyces cerevisiae*)
- **Protist (1):** Malaria parasite (*Plasmodium falciparum*)

### TF-Exemplar (Model Organisms Focus)

| Attribute | Value |
|-----------|-------|
| **Training Cells** | 110 million |
| **Species** | 5 (human + 4 model organisms) |
| **Trainable Parameters** | 542 million |
| **Non-trainable Parameters** | 282 million |
| **Vocabulary Size** | 110,290 genes |

**Species:** Human, Mouse, Zebrafish, Fruit fly, *C. elegans*

### TF-Sapiens (Human-Only)

| Attribute | Value |
|-----------|-------|
| **Training Cells** | 57 million |
| **Species** | 1 (human only) |
| **Trainable Parameters** | 368 million |
| **Non-trainable Parameters** | 61 million |
| **Vocabulary Size** | 23,829 genes |

---

## Training Data

### Data Sources

TranscriptFormer's training corpus was assembled from multiple sources:

**Primary Sources:**
1. **CZ CELLxGENE Database:**
   - 57 million human cells (644 datasets)
   - 27 million mouse cells (88 datasets)
   - Filtered to exclude spatial assays and duplicate cells
   - Both normal (67M cells) and diseased (17M cells across 123 diseases) retained

2. **Tabula Sapiens** (human reference atlas)
3. **ZebraHub** (zebrafish developmental atlas)
4. **Bgee** (cross-species expression database)
5. **Other Public Datasets:** GEO, Human Cell Atlas Portal, Cell Portal

### Data Composition

| Category | Details |
|----------|---------|
| **Total Cells** | 112 million (TF-Metazoa) |
| **Species** | 12 |
| **Evolutionary Span** | 1.53 billion years |
| **Normal Cells** | ~67 million |
| **Diseased Cells** | ~17 million (123 disease types) |
| **Datasets** | 700+ |

### Data Processing

1. **Quality Control:**
   - UMI count filtering
   - Log1p normalization of expression values

2. **Count Processing:**
   - Values clipped at 30 (configurable)
   - No explicit feature selection required
   - Raw counts preferred over normalized

3. **Gene Filtering:**
   - Only genes with ESM-2 embeddings included
   - Out-of-vocabulary genes filtered during inference

4. **Species Balancing:**
   - Up-weighting of low-resource species
   - Prevents human/mouse dominance

---

## Training Methodology

### Training Objective

TranscriptFormer is trained via **log-likelihood maximization** on the joint probability of gene-expression tuples:

```
L = -Σ log P(gene_i, count_i | gene_<i, count_<i)
```

The model learns:
1. **Gene prediction loss:** Cross-entropy for gene identity
2. **Count prediction loss:** Zero-truncated Poisson NLL for expression counts

### Self-Supervised Learning

**Key Property:** TranscriptFormer is trained without any metadata labels (cell type, species, tissue, donor) except for assay technology identifiers.

This approach:
- Enables learning from massive unlabeled datasets
- Forces the model to discover biologically meaningful features from data
- Avoids overfitting to human annotation artifacts

### Training Configuration

| Parameter | Value |
|-----------|-------|
| **Optimizer** | AdamW |
| **Learning Rate Schedule** | Linear warm-up → cosine decay |
| **Global Batch Size** | ~4-5 million tokens |
| **Precision** | Mixed fp16/bf16 |
| **Infrastructure** | H100 cluster with DDP |
| **Total Tokens Processed** | ~3.5 trillion (~20 epochs) |

### Data Augmentation

**Gene Shuffling:** Expressed genes are randomly permuted each batch to remove positional bias (since genes in a cell don't have inherent order).

---

## Capabilities & Applications

### Core Capabilities

1. **Cell Type Classification**
   - Zero-shot classification across species
   - Transfer learning to unseen species
   - Robust to rare/challenging cell types

2. **Disease State Identification**
   - Zero-shot detection of disease states
   - SARS-CoV-2 infection detection without COVID-specific training
   - Host-pathogen interaction analysis

3. **Cross-Species Annotation Transfer**
   - Transfer cell type labels across species
   - Maintain performance at >600M years evolutionary distance
   - Spermatogenesis cell type transfer (mammals ↔ birds)

4. **Gene-Gene Interaction Prediction**
   - Predict transcription factor targets via prompting
   - Generate cell type-specific regulatory relationships
   - Point-wise mutual information (PMI) analysis

5. **Contextualized Gene Embeddings (CGE)**
   - Gene-specific representations per cell
   - Capture context-dependent gene behavior
   - Encode cell type, tissue, donor variations

### Virtual Instrument Paradigm

TranscriptFormer functions as a "virtual instrument" for biological inquiry:

```
Prompt: [marker genes for cell type X]
→ Model predicts: likely transcription factors, gene interactions
```

This enables:
- In silico hypothesis generation
- Computational experimentation before wet lab work
- Rapid screening of biological questions

---

## Benchmark Performance

### Human Cell Type Classification (Tabula Sapiens v2)

| Model | Macro F1 Score |
|-------|----------------|
| **TF-Metazoa** | **0.910 ± 0.001** |
| **TF-Exemplar** | **0.907 ± 0.001** |
| TF-Sapiens | 0.906 ± 0.001 |
| UCE | 0.906 ± 0.001 |
| scGPT | Comparable |
| Geneformer | Lower |

**Note:** Benchmark approaching saturation; most models achieve near-perfect accuracy on abundant cell types.

### Challenging Cell Types

TranscriptFormer excels on difficult cell types where other models struggle:
- Myeloid leukocytes
- T-cells and subtypes
- Innate lymphoid cells
- Rare cell populations

Multi-species pretraining (TF-Metazoa, TF-Exemplar) provides advantages over human-only training.

### Disease State Prediction (COVID-19 Lung)

| Model | F1 Score (Infected vs. Healthy) |
|-------|--------------------------------|
| **TranscriptFormer** | **~0.85-0.86** |
| Baselines | ~0.80 |

TranscriptFormer distinguishes SARS-CoV-2-infected cells from healthy cells **without any COVID-specific training**, demonstrating emergent disease detection capability.

---

## Cross-Species Generalization

### Out-of-Distribution Species Performance

TranscriptFormer was evaluated on species **never seen during training**:

| Species | Evolutionary Distance | TF-Metazoa F1 | UCE F1 |
|---------|----------------------|---------------|--------|
| Mouse lemur (*M. murinus*) | 74 Mya | High | Moderate |
| Zebrafish (*D. rerio*) | 435 Mya | High | Moderate |
| Sea lamprey (*P. marinus*) | 615 Mya | Moderate | Low |
| Stony coral (*S. pistillata*) | 685 Mya | **>0.7** | <0.6 |

**Key Finding:** TF-Metazoa outperforms UCE by ~15% on out-of-distribution species, particularly at larger evolutionary distances.

### Species Transfer Learning (Spermatogenesis)

Cross-species transfer of cell type annotations across mammals and birds:

| Transfer Direction | TF-Metazoa | UCE |
|--------------------|------------|-----|
| Human → Chicken (310 Mya) | 0.448 ± 0.006 | 0.251 ± 0.004 |
| Average across species pairs | Strong | Moderate |

### Why Multi-Species Training Helps

1. **Evolutionary Diversity:** More diverse training data improves generalization
2. **Conserved Patterns:** Model learns evolutionarily conserved gene expression programs
3. **Disproportionate Benefit:** Even small additions of diverse species (2M cells) significantly boost performance

---

## Generative Prompting

### Transcription Factor Prediction

TranscriptFormer can be prompted to predict cell type-specific transcription factors:

**Approach:**
1. Condition model on marker genes for specific cell type
2. Generate probability distributions for transcription factors
3. Analyze point-wise mutual information (PMI) for gene associations

**Results:**
- Recovers 70-87% of known regulatory relationships
- Predictions align with independent experimental observations
- PMI spectra match STRING v12.0 annotations

### Cell Type-Conditioned Gene Interactions

Model-generated heatmaps of transcription factor expression across 112 human cell types (from Tabula Sapiens 2.0) recapitulate known features:
- Diagonal pattern of cell type-specific TFs
- Vertical bands of broadly expressed regulators

### Validation

Gene-gene interaction predictions validated against:
- STRING database (protein-protein interactions)
- ChIP-seq data (transcription factor binding)
- Independent experimental observations

---

## Comparison with Other Models

### TranscriptFormer vs. Geneformer

| Aspect | TranscriptFormer | Geneformer |
|--------|------------------|------------|
| **Paradigm** | Generative autoregressive | Masked language modeling |
| **Cross-species** | Native (12 species) | Human-only (requires ortholog conversion) |
| **Gene encoding** | ESM-2 protein embeddings | Rank value encoding |
| **Expression modeling** | Explicit Poisson likelihood | Implicit in ranking |
| **Zero-shot performance** | Strong | Limited (per independent studies) |
| **In silico perturbation** | Via prompting | Core feature |
| **Training data** | 112M cells | 30M (V1) / 104M (V2) cells |
| **Publication** | bioRxiv 2025 | Nature 2023 |

### TranscriptFormer vs. UCE

| Aspect | TranscriptFormer | UCE |
|--------|------------------|-----|
| **Cross-species** | Strong (trained on 12 species) | Strong (ESM-2 based) |
| **Out-of-distribution species** | +15% over UCE | Baseline |
| **Generative** | Yes | No (embedding only) |
| **Cell type classification** | Comparable/superior | Comparable |
| **Batch correction** | Moderate | Strong |

### Zero-Shot Limitations (Field-Wide)

A 2025 Genome Biology study revealed that both Geneformer and scGPT underperform simpler baselines (HVG selection, scVI, Harmony) in zero-shot settings:
- Cell type clustering: Foundation models often worse than baselines
- Batch integration: Geneformer particularly struggles
- **TranscriptFormer addresses this:** Shows strong zero-shot performance on unseen species

---

## Installation & Usage

### Requirements

- **Python:** ≥3.11
- **GPU:** A100 40GB recommended (can run on 16GB with batch_size=2)
- **Multi-GPU:** Supported via DDP for large datasets

### Installation

**From PyPI:**
```bash
# Create virtual environment
uv venv --python=3.11
source .venv/bin/activate

# Install
uv pip install transcriptformer
```

**From Source:**
```bash
git clone https://github.com/czi-ai/transcriptformer.git
cd transcriptformer
uv venv --python=3.11
source .venv/bin/activate
uv pip install -e .
```

### Download Models

```bash
# Download specific model
transcriptformer download tf-sapiens
transcriptformer download tf-exemplar
transcriptformer download tf-metazoa

# Download all models and embeddings
transcriptformer download all

# Download only embedding files
transcriptformer download all-embeddings
```

### Running Inference

**Basic Cell Embedding Extraction:**
```bash
transcriptformer inference \
  --checkpoint-path ./checkpoints/tf_sapiens \
  --data-file data.h5ad \
  --output-path ./results \
  --batch-size 8
```

**Out-of-Distribution Species:**
```bash
transcriptformer inference \
  --checkpoint-path ./checkpoints/tf_sapiens \
  --data-file mouse_data.h5ad \
  --pretrained-embedding ./checkpoints/all_embeddings/mus_musculus_gene.h5 \
  --batch-size 8
```

**Contextual Gene Embeddings:**
```bash
transcriptformer inference \
  --checkpoint-path ./checkpoints/tf_sapiens \
  --data-file data.h5ad \
  --emb-type cge \
  --batch-size 8
```

**Multi-GPU Inference:**
```bash
transcriptformer inference \
  --checkpoint-path ./checkpoints/tf_sapiens \
  --data-file large_data.h5ad \
  --num-gpus 4 \
  --batch-size 32
```

### Input Data Format

- **Format:** AnnData (.h5ad)
- **Gene IDs:** `var` must contain `ensembl_id` column
- **Expression:** Raw counts (unnormalized)
- **Source:** `adata.raw.X` preferred, falls back to `adata.X`

### Output Format

**Cell Embeddings (`--emb-type cell`):**
- `obsm['embeddings']`: Cell embeddings
- `obs`: Original metadata preserved
- `uns['llh']`: Log-likelihood scores

**Contextual Gene Embeddings (`--emb-type cge`):**
- `uns['cge_embeddings']`: Gene embeddings (n_gene_instances × embedding_dim)
- `uns['cge_cell_indices']`: Cell index for each embedding
- `uns['cge_gene_names']`: Gene names

### Python API (via Helical)

```python
from helical.models.transcriptformer import TranscriptFormer, TranscriptFormerConfig

# Initialize model
config = TranscriptFormerConfig(model_name="tf-exemplar")
model = TranscriptFormer(config)

# Process data
dataset = model.process_data(adata, gene_names="ensembl_id")

# Get embeddings
embeddings = model.get_embeddings(dataset)
```

---

## Limitations

### Technical Limitations

1. **No Spatial Transcriptomics:** Does not handle spatial data directly
2. **No Multi-Modal Data:** Cannot integrate ATAC-seq, proteomics natively
3. **GPU Memory:** Large batch inference requires significant VRAM
4. **Perturbation Prediction:** Not specialized for zero-shot perturbation prediction (unlike Geneformer's focus)

### Biological Limitations

1. **Sampling Bias:** Overrepresentation of well-studied tissues/species (human, mouse)
2. **Demographic Underrepresentation:** Human data may reflect limited genetic ancestries
3. **Extreme Distances:** Performance degrades for very distant species or unusual tissues
4. **Batch Effects:** Does not explicitly handle batch effects (shared limitation with most foundation models)

### Model Risks

1. **Hallucinations:** May produce biologically implausible gene combinations with insufficient context
2. **Not Validated for Clinical Use:** Predictions should not substitute for experimental validation
3. **Homology-Based Data Leakage:** Conserved genes may leak information between training/validation sets

---

## Ecosystem & Resources

### Official Resources

| Resource | URL |
|----------|-----|
| **GitHub Repository** | https://github.com/czi-ai/transcriptformer |
| **Virtual Cells Platform** | https://virtualcellmodels.cziscience.com/model/transcriptformer |
| **Quickstart Tutorial** | https://virtualcellmodels.cziscience.com/quickstart/transcriptformer-quickstart |
| **bioRxiv Paper** | https://www.biorxiv.org/content/10.1101/2025.04.25.650731 |

### Model Checkpoints

| Model | Download Size |
|-------|--------------|
| TF-Sapiens | ~3 GB |
| TF-Exemplar | ~4.2 GB |
| TF-Metazoa | ~5 GB |
| All Embeddings | ~5 GB (24 species) |

### Available ESM-2 Embeddings (24 Species)

**In Training Data:**
- Human, Mouse, Zebrafish, Fruit fly, *C. elegans*
- Rabbit, Chicken, African clawed frog
- Sea urchin, Freshwater sponge, Yeast, Malaria parasite

**Out-of-Distribution:**
- Rat, Pig, Chimpanzee, Gorilla
- Rhesus macaque, Marmoset
- Western clawed frog, Platypus, Opossum
- Naked mole-rat, Sea lamprey, Coral

### Related CZI Models

1. **GREmLN:** Graph-aware gene regulatory network model
2. **rBio:** Reasoning model trained on virtual cell simulations (built on TranscriptFormer)

### Third-Party Integrations

1. **Helical Package:** Unified interface for foundation model comparison
2. **CZ CELLxGENE:** Pre-computed embeddings available
3. **Prism ToolSpace:** Contextualized gene embedding database

---

## Key Takeaways

### Scientific Contributions

1. **First generative, multi-species single-cell foundation model** spanning 1.53 billion years of evolution
2. **Novel architecture** combining ESM-2 embeddings, expression-aware attention, and count likelihood
3. **Demonstrates that evolutionary diversity enhances generalization** beyond dataset size alone
4. **Establishes virtual instrument paradigm** for computational biology

### Practical Implications

1. **Strong zero-shot performance** on unseen species and disease states
2. **Enables cross-species analysis** without explicit ortholog mapping
3. **Accessible via CLI and Python API** with comprehensive documentation
4. **Open source** under MIT license

### Comparison Summary

| Strength | Best Model |
|----------|------------|
| Cross-species generalization | **TranscriptFormer** |
| Zero-shot cell type classification | **TranscriptFormer** |
| In silico perturbation | Geneformer |
| Batch correction | UCE, scVI |
| Computational efficiency | Smaller models |
| Generative prompting | **TranscriptFormer** |

### Future Directions

1. **Multimodal expansion:** Incorporating spatial, proteomic, genomic data
2. **Batch effect handling:** Unified approach for dataset integration
3. **Perturbation prediction:** Specialized enhancements
4. **Species expansion:** More diverse organisms
5. **Integration with rBio:** Natural language interface for biological queries

---

## References

1. Pearce, J. D., et al. (2025). A Cross-Species Generative Cell Atlas Across 1.5 Billion Years of Evolution: The TranscriptFormer Single-cell Model. bioRxiv. DOI: 10.1101/2025.04.25.650731

2. Lin, Z., et al. (2023). Evolutionary-scale prediction of atomic-level protein structure with a language model. Science, 379(6637), 1123-1130. [ESM-2]

3. Theodoris, C. V., et al. (2023). Transfer learning enables predictions in network biology. Nature, 618, 616-624. [Geneformer]

4. Rosen, Y., et al. (2024). Universal Cell Embeddings: A Foundation Model for Cell Biology. Nature Methods. [UCE]

5. Cui, H., et al. (2024). scGPT: toward building a foundation model for single-cell multi-omics using generative AI. Nature Methods. [scGPT]

6. Bunne, C., et al. (2024). How to build the virtual cell with artificial intelligence: Priorities and opportunities. Cell, 187(25), 7045-7063.

---

*Document generated: January 2026*
*Research conducted using public information from bioRxiv, CZI, and related academic sources*

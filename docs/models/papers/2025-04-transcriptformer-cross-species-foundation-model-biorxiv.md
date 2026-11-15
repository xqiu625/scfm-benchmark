# TranscriptFormer: Cross-Species Foundation Model — bioRxiv (2025-04)

## 📊 Paper Metadata
- **Title:** A Cross-Species Generative Cell Atlas Across 1.5 Billion Years of Evolution: The TranscriptFormer Single-cell Model
- **Authors:** James D Pearce, Sara E Simmonds, Gita Mahmoudabadi, Lakshmi Krishnan, Giovanni Palla, et al.
- **Publication:** bioRxiv preprint (2025)
- **Institution:** Chan Zuckerberg Initiative, Stanford University
- **Paper Link:** https://doi.org/10.1101/2025.04.25.650731
- **Code/Data:** https://virtualcellmodels.cziscience.com/paper/transcriptformer-2025-04

## 🎨 Key Figures

### Figure 1: transcriptformer-overview.png
*The former local overview image is not tracked in this repository; consult the
cited source for the figure.*

**Why this figure is exceptional:**
- **Comprehensive workflow visualization:** Shows complete pipeline from data collection across 12 species to downstream applications
- **Cross-species scope clarity:** Effectively communicates the unprecedented evolutionary span (1.53 billion years) through phylogenetic tree
- **Technical architecture integration:** Seamlessly connects biological concepts with computational framework
- **Multi-modal output demonstration:** Illustrates both embedding generation and generative capabilities

**Design principles to mimic:**
- Hierarchical information flow from data → model → outputs → applications
- Species diversity visualization through phylogenetic representation
- Clear input/output specifications with dimensionality details
- Integration of biological and computational perspectives

### Figure 2: cross-species-performance.png
*The former local performance image is not tracked in this repository; consult
the cited source for the figure.*

**Evolutionary generalization excellence:**
- **Evolutionary distance correlation:** Demonstrates robust performance across 685 million years of evolutionary divergence
- **Model comparison framework:** Systematic benchmarking against UCE, ESM2-CE baselines
- **Species transfer validation:** Cross-species reference mapping with quantitative F1 metrics
- **Phylogenetic context integration:** Performance visualization aligned with evolutionary relationships

**Template value for your studies:**
- Multi-species validation protocol design
- Evolutionary distance as experimental variable
- Cross-validation across phylogenetic groups
- Quantitative generalization assessment metrics

## 🔄 Key Scientific Insights

```python
### 1. Conceptual Innovation
- **Cross-Species Foundation Model:** First generative model spanning 1.53 billion years of evolution across 12 species
- **Generative Cell Atlas:** Virtual instrument paradigm for biological inquiry through model prompting
- **Evolutionary Embedding Space:** Shared representation capturing conserved gene expression signatures

### 2. Methodological Framework
- **TranscriptFormer Architecture:** Autoregressive transformer with novel gene-count coupling
- **Expression-Aware Attention:** Count-biased attention mechanism for bag-of-transcripts representation
- **Training Variants:**
  1. **TF-Metazoa:** All 12 species (444M parameters)
  2. **TF-Exemplar:** 5 model organisms (542M parameters)
  3. **TF-Sapiens:** Human-only (368M parameters)

### 3. Validation Strategy
**Comprehensive Evaluation Across:**
- **112 million cells** spanning diverse phylogeny
- **Zero-shot generalization** to unseen species
- **Cell type classification** across evolutionary distances
- **Disease state prediction** in human COVID-19 data
```

## 🔬 Critical Technical Details
```python
### 1. Model Architecture

# Key architecture components:
- Transformer layers: 12 layers, 16 attention heads, 2048 model dimension
- Expression-aware attention: Count-based bias in attention computation
- Dual decoder heads: Categorical gene selection + Zero-truncated Poisson counts
- ESM-2 embeddings: Evolutionary protein embeddings for cross-species alignment

### 2. Training Methodology
- **Scale:** 3.5 trillion tokens, 15 epochs for TF-Metazoa
- **Hardware:** 1000 H100 GPUs with mixed precision training
- **Sequence processing:** 2048 token limit with random gene shuffling
- **Loss function:** Combined categorical + zero-truncated Poisson NLL

### 3. Performance Achievements
- **Cross-species F1:** 0.778 ± 0.002 average across out-of-distribution species
- **Human cell type:** 0.910 ± 0.001 macro F1 on Tabula Sapiens 2.0
- **Disease detection:** 0.859 F1 on COVID-19 infection classification
- **Evolutionary robustness:** >0.65 F1 even at 685M years divergence
```

## Baseline Models, Evaluation Metrics, and Datasets
```python
### Baseline Models (4 methods)
- **Foundation Models:** UCE, scGPT, Geneformer
- **Embedding Baseline:** ESM2-CE (protein-based averaging)
- **Traditional Methods:** scVI, Harmony (integration methods)

### Evaluation Metrics
- **Primary:** Macro F1 score, cell type classification accuracy
- **Cross-species:** Reference mapping via k-NN classification
- **Generative:** Pointwise mutual information (PMI) for gene interactions

### Datasets
- **Training Data:** 112M cells from CZ CELLxGENE + curated atlases
- **Evaluation Species:** Coral, zebrafish, sea lamprey, mouse lemur, tropical frog
- **Notable Benchmarks:** Tabula Sapiens 2.0, COVID-19 lung dataset, spermatogenesis atlas
```

## 💭 Critical Research Implications
```python
### 1. Methodological Impact
- **Cross-Species Foundation Models:** Establishes framework for evolutionary-scale modeling
- **Generative Cell Biology:** Virtual instrument paradigm for hypothesis generation
- **Zero-Shot Generalization:** Demonstrates transferability across vast evolutionary distances

### 2. Single-Cell Genomics Relevance
- **Atlas Integration:** Unified framework for multi-species cell atlas construction
- **Biological Discovery:** Model prompting for transcription factor-gene interaction inference
- **Disease Modeling:** Cross-species disease state recognition capabilities

### 3. Computational Biology Insights
- **Evolutionary Conservation:** Quantifies preserved gene expression patterns across phylogeny
- **Contextual Gene Embeddings:** Captures cell type, tissue, and donor-specific variation
- **Generative Modeling:** Probabilistic framework for cellular state prediction
```

## 💻 Computational Requirements
```python
### Hardware Specifications
- **GPU:** 1000 H100 GPUs for training (16-bit mixed precision)
- **RAM:** ~4-5M token global batch size across distributed setup
- **Storage:** Multi-TB for 112M cell dataset and model checkpoints

### Software Environment
- **OS:** Linux-based cluster environment
- **Frameworks:** PyTorch Lightning, ESM-2, PyTorch Flex Attention
- **Dependencies:** Scanpy, AnnData, CZ CELLxGENE API

### Training/Inference Times
- **Model Training:** 3.5 trillion tokens, ~weeks on 1000 H100s
- **Inference Speed:** Real-time cell embedding generation
- **Scalability Limits:** 247K gene vocabulary, 2048 sequence length
``` 

## 🚀 Future Directions & Limitations
```python
### Potential Extensions
- **Multi-Modal Integration:** Addition of spatial, proteomic, and epigenetic data
- **Perturbation Prediction:** Zero-shot modeling of genetic and chemical perturbations
- **Temporal Dynamics:** Integration of developmental and aging trajectories

### Current Limitations
- **Batch Effects:** Limited handling of technical variation across datasets
- **Perturbation Modeling:** Not specialized for intervention prediction
- **Storage Constraints:** Browser storage APIs not supported in artifacts

### Open Questions
- How do learned representations scale to entire phylogenetic trees?
- Can generative prompting replace experimental cell atlas construction?
- What is the optimal evolutionary sampling strategy for foundation models?
```

## 📋 Implementation Checklist
```python
### For Reproducing Results
- [ ] Access to CZ CELLxGENE API and curated datasets
- [ ] H100 GPU cluster or equivalent computational resources
- [ ] ESM-2 protein embedding preprocessing pipeline
- [ ] PyTorch Lightning distributed training setup

### For Adapting to New Domains
- [ ] Species-specific gene annotation and orthology mapping
- [ ] Custom genome reference integration for non-model organisms
- [ ] Evaluation dataset curation for target phylogenetic groups
- [ ] Cross-species validation protocol implementation
```

## 🔗 Related Work & Context
- **UCE (Universal Cell Embeddings):** Predecessor multi-species embedding model
- **scGPT/Geneformer:** Single-cell foundation models with different architectures
- **ESM-2:** Protein language model enabling cross-species gene representation

## 🧩 Relevance to My Work

```python
### Key Overlaps with Project
- **Shared Task Type:** Foundation model development for biological sequence data
- **Shared Methodology:** Transformer architectures, cross-species generalization, zero-shot evaluation
- **Shared Data Type:** Large-scale genomic datasets, multi-species validation

### What You Can Reuse or Build Upon
- **Cross-species training strategy:** Evolutionary sampling and orthology-based gene mapping
- **Generative evaluation protocol:** PMI-based interaction inference, cell type-specific prompting
- **Architecture innovations:** Expression-aware attention, dual decoder head design
- **Benchmark design:** Out-of-distribution species evaluation, evolutionary distance metrics

### What This Paper Lacks (That You Can Contribute)
- **Perturbation Modeling:** No specialized architecture for intervention prediction
- **Temporal Dynamics:** Limited modeling of developmental or disease progression
- **Multi-Modal Integration:** Single transcriptomic modality focus

### Potential for Citation or Collaboration
- **Cite for:** Cross-species foundation model precedent, evolutionary generalization validation
- **Contact for:** Model access, dataset sharing, collaborative benchmarking
```
---
*Note based on analysis of: Pearce et al. "A Cross-Species Generative Cell Atlas Across 1.5 Billion Years of Evolution: The TranscriptFormer Single-cell Model" bioRxiv 2025*

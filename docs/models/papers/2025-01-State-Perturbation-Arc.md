# STATE Perturbation Model — Arc Institute (2025-01)

## 📊 Paper Metadata
- **Title:** Predicting cellular responses to perturbation across diverse contexts with State
- **Authors:** Abhinav K. Adduri, Dhruv Gautam, Beatrice Bevilacqua, Alishba Imran, et al.
- **Publication:** Nature (2025)
- **Institution:** Arc Institute (primary), UC Berkeley, Stanford, UCSF
- **Paper Link:** [Nature Publication]
- **Code/Data:** https://github.com/ArcInstitute/State

## 🎨 Key Figures

### Figure 1: state-architecture-overview.png
> Architecture figure is not tracked locally; consult the source paper.

**Why this figure is exceptional:**
- **Multi-scale design clarity:** Shows how State operates across genes, cells, and populations with clear visual hierarchy
- **Transformer innovation:** Demonstrates set-based attention mechanism for handling cellular heterogeneity
- **Performance validation:** Includes empirical results showing optimal cell set size (256 cells) with clear performance curves
- **Architectural completeness:** Shows both State Transition (ST) and State Embedding (SE) components working together

**Design principles to mimic:**
- Hierarchical panel layout showing different scales of operation
- Integration of biological concepts with technical architecture
- Clear methodology flow from input to output
- Quantitative performance specifications with validation curves

### Figure 2: state-comprehensive-evaluation.png
> Evaluation figure is not tracked locally; consult the source paper.

**Evaluation excellence:**
- **Comprehensive benchmarking:** Compares against 6+ baseline methods across multiple datasets
- **Biological relevance:** Uses Cell-Eval framework with biologically meaningful metrics
- **Statistical rigor:** Shows consistent improvements across all evaluation dimensions
- **Practical significance:** Demonstrates 50%+ improvements in perturbation discrimination

**Template value for your studies:**
- Multi-metric evaluation framework design
- Biological validation alongside technical metrics
- Statistical reporting with error bars and significance testing

## 🔄 Key Scientific Insights

```python
### 1. Conceptual Innovation
- **Set-based perturbation modeling:** First transformer to use self-attention across cell sets for perturbation prediction
- **Multi-scale architecture:** Combines population-level learning (ST) with cell-level embeddings (SE)
- **Heterogeneity modeling:** Explicitly accounts for biological variation within cell populations

### 2. Methodological Framework
- **State Transition (ST):** Transformer with bidirectional attention across cell sets (100M+ perturbed cells)
- **State Embedding (SE):** Self-supervised cell representation learning (167M cells)
- **Training Strategy:** Maximum Mean Discrepancy (MMD) loss for distributional matching
- **Variant Strategies:**
  1. **ST+HVG:** Direct gene expression input (2000 HVGs)
  2. **ST+SE:** Cell embedding input (600M parameter SE model)
  3. **Zero-shot adaptation:** Cross-dataset transfer learning
  4. **Multi-modal support:** Chemical, genetic, and signaling perturbations

### 3. Validation Strategy
**Comprehensive Evaluation Across:**
- **100M+ cells** spanning chemical, genetic, and signaling perturbations
- **70+ contexts** including cell lines, donors, and tissue types
- **Multiple perturbation types** (drugs, CRISPR, cytokines)
- **Cross-dataset generalization** with held-out contexts and zero-shot evaluation
```

## 🔬 Critical Technical Details
```python
### 1. State Transition Architecture

# Key architecture components:
- Transformer backbone: LLaMA or GPT2 with bidirectional attention
- Cell set size: 32-512 cells (optimal: 256 for large datasets)
- Hidden dimensions: 128-1488 depending on dataset complexity
- Loss function: Maximum Mean Discrepancy with energy distance kernel

### 2. State Embedding Model
- **Architecture:** 600M parameter encoder-decoder transformer
- **Training data:** 167M cells from Arc scBaseCount, CZ CELLxGENE, Tahoe-100M
- **Gene representation:** ESM2-based protein embeddings + expression-aware soft binning
- **Dual-axis loss:** Gene-wise and cell-wise reconstruction objectives

### 3. Performance Metrics
- **Perturbation discrimination:** 54% improvement over baselines (Tahoe-100M)
- **DE gene prediction:** 184% AUPRC improvement (genetic perturbations)
- **Effect size correlation:** 70% improvement approaching r=0.8 (Tahoe-100M)
- **Zero-shot transfer:** 17%+ improvement with State embeddings
```

## Baseline Models, Evaluation Metrics, and Datasets
```python
### Baseline Models (6 methods)
- **Simple baselines:** Perturbation mean, Context mean, Linear model
- **Deep learning:** scVI, CPA (autoencoder-based)
- **Foundation models:** scGPT (transformer-based)

### Evaluation Metrics
- **Primary:** Perturbation discrimination, Pearson correlation, AUPRC for DE genes
- **Differential Expression:** Spearman correlation of log-fold changes, DE overlap accuracy
- **Effect Size:** Correlation of perturbation effect magnitudes across conditions

### Datasets
- **Tahoe-100M:** 100M cells, 1138 drug conditions, 50 cancer cell lines
- **Parse-PBMC:** 10M cells, 90 cytokine perturbations, 12 donors, 18 cell types
- **Replogle-Nadig:** 624K cells, 1677 genetic perturbations, 4 cell lines
```

## 💭 Critical Research Implications
```python
### 1. Methodological Impact
- **First reliable improvement over linear baselines:** Addresses fundamental challenge in perturbation modeling
- **Set-based attention for biology:** New paradigm for handling cellular heterogeneity
- **Cross-modal transfer:** Demonstrates learning across chemical, genetic, and signaling perturbations

### 2. Single-Cell Genomics Relevance
- **Foundation model for perturbations:** Scales to 100M+ cells with consistent improvements
- **Zero-shot prediction:** Enables perturbation prediction in novel cellular contexts
- **Multi-scale modeling:** Bridges molecular, cellular, and population-level effects

### 3. Drug Discovery Insights
- **Virtual screening potential:** Can predict drug effects without experimental testing
- **Context-specific responses:** Models cell-type and condition-dependent drug effects
- **Mechanism discovery:** Attention maps reveal cellular heterogeneity patterns
```

## 💻 Computational Requirements
```python
### Hardware Specifications
- **GPU:** Multiple H100 GPUs for SE pretraining, single A100 sufficient for ST training
- **RAM:** 500GB+ for large dataset loading and processing
- **CPU:** High-core count for data preprocessing and cell set construction
- **Storage:** 10TB+ for storing large-scale single-cell datasets

### Software Environment
- **OS:** Linux (Ubuntu/CentOS recommended)
- **Frameworks:** PyTorch Lightning, HuggingFace Transformers, Scanpy
- **Dependencies:** Flash Attention 2, geomloss (MMD), distributed training support

### Training/Inference Times
- **SE pretraining:** Days to weeks on multi-node GPU clusters (167M cells)
- **ST training:** Hours to days depending on dataset size and complexity
- **Inference:** Real-time prediction for new perturbations once trained
```

## 🚀 Future Directions & Limitations
```python
### Potential Extensions
- **Multi-modal integration:** Incorporate proteomics, metabolomics, morphology data
- **Temporal dynamics:** Model perturbation effects over time courses
- **Combinatorial perturbations:** Predict effects of multiple simultaneous interventions
- **Mechanistic interpretation:** Enhanced attention analysis for biological insights

### Current Limitations
- **Computational requirements:** Large-scale training requires significant resources
- **Data dependency:** Performance scales with available perturbation data
- **Context boundaries:** Limited evaluation on entirely novel tissue types or species

### Open Questions
- How does State performance scale to even larger datasets (1B+ cells)?
- Can the approach generalize to non-human organisms and primary tissues?
- What are the limits of zero-shot transfer across perturbation modalities?
```

## 📋 Implementation Checklist
```python
### For Reproducing Results
- [ ] Access to large-scale single-cell datasets (Tahoe-100M, scBaseCount)
- [ ] Multi-GPU compute infrastructure (preferably H100s for SE training)
- [ ] PyTorch Lightning environment with distributed training support
- [ ] Cell-Eval evaluation framework implementation

### For Adapting to New Domains
- [ ] Domain-specific cell type considerations and validation
- [ ] Perturbation-specific data requirements and preprocessing
- [ ] Model architecture adjustments for different data scales
- [ ] Evaluation metric adaptation for domain-relevant outcomes
```

## 🔗 Related Work & Context
- **CellOT/Neural OT:** State provides theoretical connection to optimal transport while being more flexible
- **scGPT/Foundation models:** State specifically designed for perturbation prediction vs. general cell modeling
- **CPA/scVI:** State's transformer architecture vs. VAE-based approaches, with superior generalization

## 🧩 Relevance to My Work

```python
### Key Overlaps with Project
- **Shared Task Type:** Perturbation prediction across cellular contexts
- **Shared Methodology:** Transformer architecture, foundation model pretraining, multi-scale learning
- **Shared Data Type:** Single-cell RNA-seq, perturbation screens, multi-modal omics

### What You Can Reuse or Build Upon
- **Set-based attention mechanism:** Novel approach to handling cellular heterogeneity
- **Multi-scale architecture:** Combining population and single-cell level modeling
- **Cell-Eval framework:** Comprehensive biological evaluation metrics
- **Zero-shot transfer strategy:** Cross-dataset and cross-modal generalization approach

### What This Paper Lacks (That You Can Contribute)
- **Limited temporal modeling:** No trajectory-aware perturbation dynamics
- **Single timepoint focus:** Lacks modeling of perturbation kinetics and time courses
- **Mechanistic interpretation:** Attention analysis could be more biologically interpretable
- **Primary tissue validation:** Most validation on cell lines rather than primary cells

### Potential for Citation or Collaboration
- **Cite for:** State-of-the-art perturbation prediction baseline, Cell-Eval metrics, set-based attention
- **Contact for:** Model weights, training pipelines, evaluation frameworks, potential collaboration on extensions
```

---
*Note based on analysis of: Adduri et al. "Predicting cellular responses to perturbation across diverse contexts with State" Nature (2025)*

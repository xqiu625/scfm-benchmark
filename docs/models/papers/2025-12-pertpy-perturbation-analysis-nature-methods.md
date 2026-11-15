# Pertpy: Perturbation Analysis — Nature Methods (2025-12)

## Paper Metadata
- **Title:** Pertpy: an end-to-end framework for perturbation analysis
- **Authors:** Lukas Heumos, Yuge Ji, Lilly May, Tessa D. Green, Stefan Peidli, et al.
- **Publication:** Nature Methods (2025)
- **Institution:** Helmholtz Center Munich, Technical University of Munich, Harvard Medical School, Wellcome Sanger Institute
- **Paper Link:** https://doi.org/10.1038/s41592-025-02909-7
- **Code/Data:** https://github.com/scverse/pertpy (Apache 2.0 license)
- **Tags:** `single-cell`, `perturbation-analysis`, `CRISPR-screens`, `drug-response`, `computational-framework`, `Python`

---

## One-Sentence Summary
Pertpy is a Python-based modular framework within the scverse ecosystem that provides end-to-end analysis of large-scale single-cell perturbation experiments, integrating metadata annotation, perturbation space learning, differential analysis, and compositional changes with fast, scalable implementations.

---

## Problem Statement

### Biological/Scientific Challenge
- **Primary Problem:** Understanding cellular responses to perturbations (genetic, chemical, environmental) is crucial for describing biological mechanisms, but current analysis tools are fragmented and lack integration
- **Current Limitations:** Existing perturbation analysis frameworks (MUSIC, ScMAGeCK, SCEPTRE, GSFA, FR-Perturb) primarily focus on CRISPR perturbation analysis, neglecting other perturbation data types
- **Gap in Knowledge:** No framework exists that scales to genome-scale datasets, contextualizes data with public annotations, and uses common data structures across tools

### Technical Challenges
- Challenge 1: Growing size and complexity of high-throughput perturbation screens pose interpretation challenges
- Challenge 2: Lack of meaningful lower-dimensional representations for perturbation data
- Challenge 3: Many existing tools suffer from maintenance issues or are confined to the R ecosystem
- Challenge 4: No unified framework for integrating external metadata with measured data

---

## Key Innovation

### Core Methodological Contribution
**Pertpy:** A Python-based modular framework for end-to-end analysis of large-scale single-cell perturbation experiments, providing access to harmonized datasets, metadata databases, and over 100 composable analysis functions organized into interoperable modules.

**Why It Matters:**
- **Advancement 1:** Integrates external metadata from databases (DepMap, GDSC, CMap, PubChem, ChEMBL) for unprecedented contextualization
- **Advancement 2:** Provides sparse and memory-efficient implementations leveraging JAX for GPU acceleration, substantially faster than original implementations
- **Advancement 3:** Part of scverse ecosystem enabling seamless interoperability with existing single-cell analysis workflows (scanpy, scvi-tools, squidpy)

### Technical Architecture
```python
# Core components and workflow
1. Data Transformation:
   - gRNA assignment (threshold-based or Poisson-Gaussian mixture model)
   - Batch effect correction via scverse ecosystem integration
   - Innovation: Supports both low and high MOI scenarios

2. Metadata Annotation:
   - Cell line annotation from DepMap, GDSC, CCLE, OncoTree
   - Drug annotation from CMap, PubChem, ChEMBL
   - Innovation: API-based automatic metadata enrichment

3. Perturbation Space Generation:
   - Pseudobulk, Centroid, MLP classifier, Logistic regression spaces
   - Innovation: Per-perturbation embeddings summarizing cellular responses

4. Downstream Analysis:
   - Differential expression (PyDESeq2, edgeR, Wilcoxon, t-test)
   - Compositional analysis (scCODA 2.0, tascCODA 2.0, Milo)
   - MCPs via DIALOGUE with fast linear programming solver
   - 18+ distance metrics for perturbation comparison
```

---

## Main Results

### 1. Method Validation
**Performance Metrics:**
- **Runtime:** Substantially faster than original R implementations (Extended Data Fig. 1)
- **Memory:** More efficient memory usage through JAX-based implementations
- **Mixscape Agreement:** 0.97 correlation with original Seurat implementation

**Key Achievements:**
- Over 100 composable analysis functions in organized modules
- 30+ public perturbational single-cell datasets available through data loaders
- Successful validation of all reimplemented methods against original implementations

### 2. Biological Discoveries

**Finding 1:** CRISPRa Screen Analysis (Norman et al. dataset)
- **Evidence:** 111,255 cells, 287 perturbations in K562 cells
- **Discovery:** TP73 (labeled as pioneer factor) clusters with G1 cell cycle perturbations; pro-growth cluster separates into KLF-targeting and MAPK-encoding subclusters
- **Significance:** Enables refinement of gene program annotations and discovery of unannotated clusters (neutrophil degranulation pathway)

**Finding 2:** Drug Response Deconvolution (MIX-Seq dataset)
- **Evidence:** 154,710 cells, 172 cell lines, 13 drug treatments
- **Discovery:** Identified viability-dependent (interferon signaling) and viability-independent (protein translation) responses to dabrafenib
- **Significance:** Protein translation upregulation mechanism distinct from BRAF inhibition was previously unknown

**Finding 3:** TNBC Treatment Response
- **Evidence:** 22 patients, pre/post treatment with NACT +/- anti-PD-L1
- **Discovery:** Compositional shifts in CD4/CD8 T cell populations; MCP2 associated with treatment response involving IL-7/IL-7R interaction and AP-1 complex genes (JUN, FOS, FOSB)
- **Significance:** Identified potential biomarkers for treatment response prediction

### 3. Application Demonstrations

**Use Case 1: CRISPRa Screen Analysis**
- **Setup:** Norman et al. dataset with MLP classifier perturbation space
- **Results:** Gene programs cluster by biological function; discovered previously unannotated gene programs
- **Insight:** Perturbation space learning can refine and extend existing annotations

**Use Case 2: Drug Response Analysis**
- **Setup:** MIX-Seq dataset with metadata annotation from multiple databases
- **Results:** 0.88 average Pearson correlation with CCLE bulk expression; successful viability deconvolution
- **Insight:** Integration of external databases enables sophisticated drug response analysis

**Use Case 3: Cancer Treatment Response**
- **Setup:** TNBC dataset with scCODA and DIALOGUE analysis
- **Results:** Identified cell type compositional changes and multicellular programs associated with treatment response
- **Insight:** Combined compositional and MCP analysis reveals treatment mechanisms

---

## Key Figures for Learning

### Figure 1: Framework Overview
**What Makes This Figure Excellent:**
- **Visual Design:** Clear modular architecture showing data flow from input through preprocessing to downstream analysis
- **Technical Content:** Shows all major components: data transformation, metadata annotation, preprocessing, perturbation space, and knowledge inference modules
- **Pedagogical Value:** Comprehensive overview of end-to-end perturbation analysis pipeline

**Design Principles to Adopt:**
1. Modular workflow visualization with numbered steps
2. Clear separation between preprocessing and analysis modules
3. Integration of multiple downstream analysis types in unified framework

**Teaching Points:**
- Perturbation analysis involves multiple interconnected steps
- Metadata enrichment is crucial for biological interpretation
- Multiple downstream analyses can be derived from same preprocessed data

### Figure 2: Perturbation Space Learning
**Educational Value:**
- **For Beginners:** Shows how perturbations with similar effects cluster together
- **For Advanced:** Demonstrates how classifier-based embeddings can reveal gene program structure
- **For Practitioners:** Provides workflow for discovering and annotating gene programs

**Reusable Templates:**
- UMAP visualization with gene program coloring
- Cluster disentanglement analysis workflow
- Pathway enrichment bar plots

### Figure 3: Drug Response Deconvolution
**Educational Value:**
- **For Beginners:** Linear regression approach to separate viability effects
- **For Advanced:** Integration of multiple drug sensitivity databases
- **For Practitioners:** Volcano plot visualization for viability-dependent/independent genes

### Figure 4: Treatment Response Analysis
**Educational Value:**
- **For Beginners:** Compositional analysis basics with scCODA
- **For Advanced:** MSE distance matrices for treatment comparison
- **For Practitioners:** Box plots for cell type proportion changes

---

## Technical Implementation

### Computational Pipeline
```python
### 1. Data Processing
Input: AnnData/MuData objects with perturbation annotations
Processing steps:
  - gRNA assignment (threshold or Poisson-Gaussian mixture)
  - Quality control and batch effect removal
  - Metadata annotation from public databases
Output: Annotated AnnData with perturbation labels

### 2. Core Algorithms
Mixscape: Classifies cells as successfully/unsuccessfully perturbed
  - Parameters: k=20 nearest neighbors, probability threshold=0.5
scCODA 2.0: Bayesian compositional analysis
  - Parameters: NumPyro no-U-turn sampler, FDR=0.1
DIALOGUE: Multicellular program discovery
  - Parameters: Linear programming solver, n_MCPs user-defined

### 3. Analysis & Validation
Distance metrics (18+ implemented):
  - MSE, Energy distance, Wasserstein, MMD, etc.
Statistical tests:
  - Monte Carlo permutation testing for distances
  - t-tests with Benjamini-Hochberg correction for MCPs
```

### Software & Hardware Requirements
**Minimal Setup:**
- **Compute:** CPU sufficient for small datasets (<50,000 cells)
- **Memory:** 16GB RAM recommended
- **Storage:** Varies by dataset size

**Optimal Setup:**
- **Compute:** GPU for JAX acceleration (CUDA-compatible)
- **Memory:** 64GB+ RAM for large datasets
- **Time:** Minutes to hours depending on analysis and dataset size

**Software Dependencies:**
```bash
# Key packages
- pertpy: 0.10.0+ (main framework)
- scanpy: single-cell analysis foundation
- anndata/mudata: data structures
- jax: GPU acceleration
- numpyro: Bayesian inference
- ott-jax: optimal transport
- pytorch: neural network training
```

---

## Benchmarking & Comparisons

### Baseline Methods Compared
**CRISPR Analysis Frameworks:**
- MUSIC: R-based, limited to CRISPR - lacks metadata integration
- ScMAGeCK: R/Bioconductor, CRISPR-specific - no compositional analysis
- SCEPTRE: R, guide assignment focus - limited downstream analysis
- GSFA: R, latent factors - no perturbation spaces
- FR-Perturb: R/Python, Euclidean distance only - limited distance metrics

### Performance Comparison Table
| Feature | pertpy | MUSIC | ScMAGeCK | SCEPTRE | GSFA | FR-Perturb |
|---------|--------|-------|----------|---------|------|------------|
| Language | Python | R | R | R | R | R, Python |
| Ecosystem | scverse | - | Bioconductor | - | - | - |
| Dataloaders | Yes | No | No | No | No | No |
| Metadata annotation | Yes | No | No | No | No | No |
| Guide assignment | Yes | No | No | Yes | No | No |
| Compositional analysis | Yes | No | No | No | No | No |
| MCPs | Yes | No | No | No | Yes | Yes |
| Distance metrics | 18+ | Topic-based | Selection scores | No | Euclidean | Euclidean |
| Perturbation embeddings | Multiple | Topics | Linear regression | No | Latent factors | Latent factors |

### When to Use This Method vs Alternatives
**Use pertpy when:**
- Working with diverse perturbation types (CRISPR, drugs, disease states)
- Need integration with scverse ecosystem (scanpy, scvi-tools)
- Require metadata contextualization from public databases
- Analyzing large-scale perturbation screens (>100,000 cells)
- Need multiple distance metrics or perturbation space methods

**Use alternatives when:**
- Working exclusively in R ecosystem
- Need specialized CRISPR guide assignment only (SCEPTRE)
- Simple differential expression sufficient

---

## Research Implications

### For Single-Cell Genomics
- **Methodological:** Unified framework for perturbation experiments across modalities
- **Analytical:** Novel perturbation space representations enable new downstream analyses
- **Biological:** Integration of metadata enables biology-informed discoveries

### For Foundation Models & AI
- **Training Data:** Curated perturbation datasets with harmonized ontologies for model training
- **Benchmarking:** Distance metrics provide evaluation framework for perturbation prediction models
- **Architecture Design:** Perturbation spaces inform multimodal model design

### For Drug Discovery
- **Target Identification:** Cell type prioritization via Augur identifies responsive populations
- **Perturbation Screens:** Scalable analysis of genome-scale CRISPR and drug screens
- **Patient Stratification:** Treatment response analysis enables biomarker discovery

### For Cancer Research
- **TME Analysis:** Compositional changes reveal treatment-induced immune remodeling
- **Resistance Mechanisms:** MCP analysis identifies potential resistance signatures
- **Combination Therapy:** CINEMA-OT enables synergistic effect analysis

---

## Limitations & Considerations

### Technical Limitations
- **Scalability:** Out-of-memory Dask implementations needed for >100M cells (Tahoe-100M)
- **Mixscape Assumptions:** Assumes perturbation effects are additive; may miss complex interactions
- **MCP Validation:** Multiple testing correction across MCPs may be needed for strict inference

### Biological Limitations
- **Cell Type Dependence:** Some analyses require pre-defined cell type labels
- **Perturbation Types:** Different perturbation types may require different analytical approaches
- **Temporal Dynamics:** Static snapshots may miss dynamic perturbation responses

### Practical Constraints
- **Cost:** Computational resources for large datasets; database API access
- **Expertise:** Requires familiarity with single-cell analysis and Python
- **Time:** Complex pipelines may require iterative parameter tuning
- **Sample Requirements:** Some methods require multiple replicates per condition

---

## Future Directions

### Immediate Extensions (< 1 year)
- Integration with spatial transcriptomics (Squidpy combination)
- Out-of-memory implementations via Dask for massive datasets

### Medium-term Goals (1-3 years)
- Perturbation atlas creation with harmonized data collection
- Foundation model integration for perturbation prediction evaluation
- Multimodal perturbation data analysis (transcriptomics + proteomics)

### Long-term Vision (3+ years)
- Real-time perturbation response prediction in clinical settings
- Comprehensive perturbation atlases characterizing inducible cell states
- Integration with generative models for perturbation design

### Open Questions
- How to optimally combine different perturbation space methods?
- What distance metrics best capture biological perturbation effects?
- How to extend compositional analysis to spatial contexts?

---

## Connections to AI4Bio Learning Hub

### Relevant Tutorial Categories
- **Single-Cell Analysis Fundamentals:**
  - Concepts to teach: AnnData structure, preprocessing, batch correction
  - Code examples: Standard scanpy workflow integration with pertpy

- **Perturbation Biology:**
  - Concepts to teach: CRISPR mechanisms, drug response, compositional changes
  - Code examples: Mixscape analysis, scCODA usage, distance calculations

- **Machine Learning for Biology:**
  - Concepts to teach: Classifier-based embeddings, optimal transport, Bayesian inference
  - Code examples: MLP perturbation space, CINEMA-OT counterfactual analysis

### Cross-References
- **Prerequisite Concepts:** scanpy basics, AnnData/MuData structures, single-cell QC
- **Related Methods:** scvi-tools (VAE-based), Squidpy (spatial), decoupler-py (functional)
- **Application Examples:** Drug discovery pipelines, cancer immunotherapy analysis

### Teaching Modules to Create
```python
Module 1: Introduction to Perturbation Analysis
- Learning objectives: Understand perturbation types, data structures, preprocessing
- Hands-on exercise: Load and preprocess CRISPR screen data
- Expected time: 2 hours

Module 2: Perturbation Space Learning
- Learning objectives: Compare space methods, interpret embeddings
- Hands-on exercise: Generate and compare pseudobulk vs MLP spaces
- Expected time: 3 hours

Module 3: Compositional Analysis
- Learning objectives: Bayesian compositional modeling, Milo neighborhoods
- Hands-on exercise: Apply scCODA to treatment dataset
- Expected time: 2 hours
```

---

## Related Work & Context

### Builds Upon
- **scverse ecosystem (2023):** Provides computational ecosystem foundation
- **Perturb-seq (2016):** Established pooled CRISPR screening with single-cell readouts
- **scPerturb (2024):** Harmonized perturbation dataset collection

### Complements
- **scvi-tools:** Probabilistic modeling for integration and perturbation prediction
- **Squidpy:** Spatial omics analysis for spatial perturbation studies
- **decoupler-py:** Functional inference from perturbation signatures

### Compared To
- **R-based frameworks:** Python alternative with scverse integration
- **Method-specific tools:** Unified framework vs specialized single-task tools

### Inspired Follow-ups
- Tahoe-100M: Giga-scale perturbation atlas leveraging pertpy infrastructure
- PerturBase: Database extending scperturb with pertpy compatibility

---

## Educational Resources

### For Students Learning This Topic
**Prerequisites:**
- Basic single-cell RNA-seq analysis (scanpy)
- Python programming fundamentals
- Understanding of experimental perturbation types

**Recommended Reading Order:**
1. scverse ecosystem paper (Virshup et al. 2023)
2. This paper (Heumos et al. 2025)
3. Specific method papers (scCODA, DIALOGUE, CINEMA-OT)

### Hands-On Learning
**Available Datasets:**
- Norman et al. CRISPRa screen (GSE133344)
- McFarland MIX-Seq drug screen (figshare)
- Zhang et al. TNBC treatment (GSE169246)

**Coding Tutorials:**
- pertpy documentation: https://pertpy.readthedocs.io
- 15+ use case tutorials covering diverse applications

**Example Notebooks:**
- CRISPR screen analysis
- Drug response deconvolution
- Treatment response analysis

---

## Practical Implementation Guide

### For Experimentalists
**Planning Your Experiment:**
- [ ] Define perturbation type (genetic, chemical, disease)
- [ ] Design appropriate controls for each condition
- [ ] Plan for sufficient replicates (n >= 3 per condition recommended)
- [ ] Consider metadata needs for downstream annotation

**Quality Control Checkpoints:**
- [ ] gRNA assignment quality (for CRISPR experiments)
- [ ] Mixscape classification of perturbation success
- [ ] Batch effect assessment and correction

### For Computational Biologists
**Setting Up Analysis:**
```bash
# Installation
pip install pertpy

# Basic usage
import pertpy as pt
import scanpy as sc

# Load data
adata = pt.dt.norman_2019()

# Metadata annotation
pt.tl.annotate_cell_lines(adata)

# Perturbation space
ps = pt.tl.PseudobulkSpace()
ps_adata = ps.compute(adata, target_col='perturbation')
```

**Common Pitfalls:**
- Not filtering cells before compositional analysis: Use consistent cell type filtering
- Ignoring batch effects: Always check and correct before perturbation analysis
- Wrong reference cell type in scCODA: Use automatic reference detection

**Troubleshooting:**
- **Issue:** Memory errors with large datasets
  - **Solution:** Use sparse representations, subset data for exploration

- **Issue:** Mixscape classifies all cells as not perturbed
  - **Solution:** Adjust perturbation signature parameters, check gRNA expression

---

## Key Takeaways for AI4Bio Hub Users

### For Beginners
1. **Main Concept:** Pertpy provides a unified toolkit for analyzing how cells respond to different treatments or genetic changes
2. **Why It Matters:** Understanding perturbation responses is key to drug discovery and disease understanding
3. **Start Here:** Begin with the Norman CRISPRa dataset tutorial

### For Intermediate Learners
1. **Technical Depth:** Master perturbation space generation and distance metrics
2. **Implementation:** Build custom pipelines combining multiple pertpy modules
3. **Next Steps:** Explore CINEMA-OT for causal inference, DIALOGUE for multicellular programs

### For Advanced Practitioners
1. **Novel Contributions:** Fast JAX implementations, unified metadata annotation, perturbation spaces
2. **Integration Opportunities:** Combine with scvi-tools for generative modeling, Squidpy for spatial
3. **Research Directions:** Contribute new distance metrics, perturbation space methods, or database integrations

---

## Citation & Attribution

```bibtex
@article{heumos2025pertpy,
  title={Pertpy: an end-to-end framework for perturbation analysis},
  author={Heumos, Lukas and Ji, Yuge and May, Lilly and Green, Tessa D and Peidli, Stefan and Zhang, Xinyue and Wu, Xichen and Ostner, Johannes and Schumacher, Antonia and Hrovatin, Karin and others},
  journal={Nature Methods},
  year={2025},
  publisher={Nature Publishing Group},
  doi={10.1038/s41592-025-02909-7}
}
```

**Related Resources:**
- **Code Repository:** https://github.com/scverse/pertpy
- **Documentation:** https://pertpy.readthedocs.io
- **Reproducibility:** https://github.com/theislab/pertpy-reproducibility
- **Author Contact:** fabian.theis@helmholtz-muenchen.de

---

*Summary prepared for AI4Bio Learning Hub | Last updated: 2025-12*
*Contributing to tutorials: Single-cell perturbation analysis, CRISPR screen analysis, Drug response prediction, Compositional analysis*

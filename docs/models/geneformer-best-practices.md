# Geneformer: Comprehensive Deep Research Report

## Executive Summary

Geneformer is a foundational transformer model developed by Christina Theodoris and colleagues at the Broad Institute of MIT and Harvard (later at Gladstone Institutes). It was pretrained on large-scale single-cell transcriptomic data to enable context-aware predictions in network biology with limited data. The model represents one of the pioneering single-cell foundation models, demonstrating that transfer learning principles from NLP can be successfully applied to gene expression data for drug target discovery and disease modeling.

---

## Model Overview

**Name:** Geneformer  
**Institution:** Broad Institute of MIT and Harvard, Gladstone Institutes, UCSF  
**Release:** Originally pretrained June 2021 (V1), expanded December 2024 (V2)  
**Publication:** Nature 2023 (V1), bioRxiv 2024 (V2)  
**License:** Apache 2.0  
**Primary Author:** Christina V. Theodoris, MD, PhD  

### Key Citations
- **V1:** Theodoris, C.V., et al. (2023). Transfer learning enables predictions in network biology. Nature 618: 616-624. DOI: 10.1038/s41586-023-06139-9
- **V2:** Chen, H., et al. (2024). Quantized multi-task learning for context-specific representations of gene network dynamics. bioRxiv 2024.08.16.608180

---

## Architecture

### Core Design: BERT-like Transformer Encoder

Geneformer uses a **BERT-like bidirectional transformer encoder** architecture, adapted for single-cell transcriptomics:

| Component | V1 Specification | V2 Specification |
|-----------|------------------|------------------|
| Architecture | Transformer Encoder | Transformer Encoder |
| Layers | 6 (original), 12 (scaled) | 12, 20 layers |
| Parameters | ~10M (6L) | 104M, 316M |
| Embedding dimensions | 256 | 256-512 |
| Attention heads | 4 per layer | Variable |
| Feed forward size | 512 | Variable |
| Input size | 2,048 genes | 4,096 genes |
| Vocabulary | ~25,424 genes | ~20,275 genes |

### Key Architectural Features

1. **Full Dense Self-Attention**: Captures complex gene-gene dependencies across the entire input context
2. **Multi-Head Attention**: Multiple attention heads per layer, each capturing distinct aspects of gene regulatory networks
3. **Transformer Encoder Blocks**: Each contains layer normalization, multi-head attention, and position-wise feed-forward networks
4. **[CLS] Token**: Special classification token for cell-level embeddings

---

## Rank Value Encoding: The Tokenization Strategy

### Core Innovation

Geneformer's most distinctive feature is its **rank value encoding** approach for representing single-cell transcriptomes:

1. **Expression Ranking**: Genes are ranked by their expression level within each cell
2. **Median Normalization**: Rankings are normalized by each gene's median expression across the entire pretraining corpus
3. **Non-parametric Representation**: Creates a robust, distribution-free encoding

### Mathematical Formulation

For each cell, genes are sorted in descending order by:
```
Rank_gene = Expression_in_cell / Median_expression_across_corpus
```

### Biological Rationale

This approach provides several advantages:
- **Deprioritizes housekeeping genes**: Ubiquitously highly-expressed genes are scaled to lower ranks
- **Prioritizes cell-state markers**: Transcription factors and other genes that distinguish cell states move to higher ranks
- **Batch effect robustness**: Relative rankings remain more stable than absolute counts across technical artifacts
- **Non-parametric**: No assumptions about expression distributions

### Token Vocabulary
- **V1**: ~25,424 protein-coding and miRNA genes detected in median of 173,152+ cells
- **V2**: ~20,275 protein-coding genes
- **Special tokens**: Padding, masking, CLS (V2), EOS (V2)

---

## Training Data

### Genecorpus-30M (V1)
- **Size**: 29.9 million human single-cell transcriptomes
- **Sources**: 561 publicly available datasets
- **Tissues**: Broad range of human tissues
- **Assembly date**: 2021
- **Quality control**: 
  - Cells with total reads within 3 SD of dataset mean
  - Mitochondrial reads within 3 SD of dataset mean
  - Minimum 7 detected Ensembl-annotated genes per cell
- **Exclusions**: Malignant cells, immortalized cell lines (high mutational burden)

### Genecorpus-104M (V2)
- **Size**: ~104 million human single-cell transcriptomes (non-cancer)
- **Assembly date**: 2024
- **Data balancing**: No tissue >25% of total data
- **Deduplication**: By DOI to prevent duplicate cells across databases
- **Additional**: Cancer corpus (~14M cells) for continual learning

### Data Filtering Criteria
- Organism: Homo sapiens
- Primary data only (is_primary_data=True)
- Non-diseased tissues (disease="normal" for general pretraining)
- No high mutational burden cells

---

## Pretraining Methodology

### Self-Supervised Masked Learning Objective

Geneformer uses a **masked language modeling (MLM)** approach adapted for transcriptomics:

1. **Masking**: 15% of genes within each transcriptome are randomly masked
2. **Prediction**: Model predicts which gene should occupy each masked position
3. **Context**: Uses remaining unmasked genes as context

### Masking Strategy (BERT-style)
- 80% of masked tokens → [MASK] token
- 10% of masked tokens → Random gene token
- 10% of masked tokens → Correct token (unchanged)

### Training Characteristics
- **Loss function**: Cross-entropy loss on masked token predictions
- **Optimization**: AdamW optimizer
- **Precision**: Full precision (fp32) for pretraining
- **Completely unsupervised**: No cell type labels required

### Emergent Knowledge

During pretraining, Geneformer learned:
- **Network hierarchy**: Encoded in attention weights
- **Gene regulatory relationships**: Context-dependent gene interactions
- **Cell state dynamics**: Understanding of transcriptional programs

---

## Model Versions and Variants

### Available Pretrained Models

| Model | Layers | Parameters | Input Size | Corpus | Date |
|-------|--------|------------|------------|--------|------|
| Geneformer-V1-10M | 6 | 10M | 2,048 | 30M cells | June 2021 |
| Geneformer-V1-12L | 12 | ~20M | 2,048 | 30M cells | June 2021 |
| Geneformer-V2-104M | 12 | 104M | 4,096 | 104M cells | Dec 2024 |
| Geneformer-V2-316M | 20 | 316M | 4,096 | 104M cells | Dec 2024 |
| Geneformer-V2-CLcancer | 12 | 104M | 4,096 | +14M cancer | Dec 2024 |

### Cancer-Tuned Model

The **Geneformer-V2-104M_CLcancer** variant was created through:
- **Continual learning**: Extended pretraining on ~14M cancer cell transcriptomes
- **Catastrophic forgetting prevention**: Included 1% of non-cancer cells during continual learning
- **Learning rate strategy**: Rewarm to max LR with cosine decay schedule
- **Purpose**: Domain-specific predictions for cancer biology

---

## Key Capabilities

### 1. Zero-Shot Learning

Geneformer can perform predictions without task-specific fine-tuning:
- **In silico perturbation analysis**: Predict effects of gene deletions/activations
- **Network connectivity inference**: Identify gene regulatory relationships
- **Cell state prediction**: Understand transcriptional programs

### 2. Transfer Learning / Fine-Tuning

With limited task-specific data, Geneformer excels at:
- **Cell state classification**: Distinguish disease vs. healthy cells
- **Gene classification**: Identify dosage-sensitive genes, TF targets
- **Chromatin dynamics prediction**: Bivalent vs. non-methylated genes
- **Network centrality prediction**: Central vs. peripheral network genes

### 3. In Silico Perturbation

Unique capability to simulate gene perturbations:
- **In silico deletion**: Remove genes from rank encoding, measure embedding shift
- **In silico activation**: Add genes to rank encoding
- **Treatment analysis**: Predict therapeutic targets that shift disease → healthy state
- **Reprogramming simulation**: Model cell fate transitions

---

## Demonstrated Applications

### 1. Dosage Sensitivity Prediction
- **Task**: Distinguish dosage-sensitive vs. insensitive transcription factors
- **Performance**: AUC 0.91 with only 10,000 training cells
- **Context-awareness**: Predictions vary by cell type (fetal vs. adult neurons)

### 2. Chromatin Dynamics
- **Bivalent vs. non-methylated genes**: Fine-tuned on only 56 loci
- **Long-range vs. short-range TFs**: Predict regulatory influence distance
- **Performance**: Consistently outperformed alternative methods

### 3. Gene Network Hierarchy
- **Central vs. peripheral genes**: NOTCH1-dependent network in endothelial cells
- **Performance**: AUC 0.81 with ~30K training cells
- **Self-supervised learning**: Network hierarchy encoded in attention weights

### 4. Cardiomyopathy Therapeutic Targets (Key Validation)

**Experimental validation of model predictions:**
- Fine-tuned to distinguish cardiomyocytes from non-failing, hypertrophic, and dilated cardiomyopathy hearts
- **Accuracy**: 90%, Precision: 82%, Recall: 87%
- Identified candidate genes whose deletion would shift disease → healthy state
- **CRISPR validation**: Knockout of predicted targets restored contractile force in iPSC-derived cardiac microtissues
- **Novel discovery**: TEAD4 identified as novel transcription factor critical for cardiomyocyte function

### 5. Synthetic Reprogramming
- Added OCT4, SOX2, KLF4, MYC to fibroblast signatures
- Model predicted shift to iPSC-like transcriptional state
- Consistent with known Yamanaka factor reprogramming

---

## Benchmark Performance

### Cell Type Clustering (Zero-Shot)

Multiple independent benchmarks show mixed results for zero-shot settings:

| Study | Finding |
|-------|---------|
| Kedzierska et al. | scGPT and Geneformer underperform HVG, Harmony, scVI in zero-shot cell type clustering |
| scEval | Geneformer benefits from fine-tuning; lower learning rate + more epochs improve performance |
| GenePT comparison | Geneformer (71% accuracy) vs. scGPT (88%) vs. GenePT (88%) on disease classification |

### With Fine-Tuning

Fine-tuning dramatically improves Geneformer's performance:
- **Cell type classification**: >90% accuracy on specific tasks (Crohn's disease dataset)
- **Dosage sensitivity**: Outperforms SVM, random forest, logistic regression
- **Chromatin dynamics**: Superior to non-pretrained attention models

### Comparison with Other Foundation Models

| Metric | Geneformer | scGPT | UCE |
|--------|------------|-------|-----|
| Zero-shot clustering | Weak | Moderate | Strong |
| Fine-tuned classification | Strong | Strong | Strong |
| In silico perturbation | Strong | Moderate | Limited |
| Cross-species | Requires orthologs | Requires orthologs | Native via ESM2 |
| Batch correction | Moderate | Moderate | Strong |

### Key Benchmark Insights

1. **Zero-shot limitations**: Foundation models including Geneformer often underperform simpler baselines (PCA, HVG) in pure zero-shot settings
2. **Fine-tuning advantage**: With even limited task-specific data, Geneformer shows significant improvements
3. **Complementary strengths**: Different models excel at different tasks; ensemble approaches can improve performance

---

## Installation and Usage

### Installation

```bash
# Ensure git-lfs is installed
git lfs install

# Clone repository
git clone https://huggingface.co/ctheodoris/Geneformer
cd Geneformer
pip install .
```

### Basic Tokenization

```python
from geneformer import TranscriptomeTokenizer

# Initialize tokenizer with custom attributes
tk = TranscriptomeTokenizer(
    {"cell_type": "cell_type", "organ_major": "organ"}, 
    nproc=4,
    model_version='V2'  # or 'V1' for original model
)

# Tokenize data
tk.tokenize_data(
    "data_directory",      # Directory with .loom or .h5ad files
    "output_directory",    # Output directory
    "output_prefix"        # Prefix for output files
)
```

### Input Data Requirements

- **Format**: Raw counts scRNA-seq data as .loom or .h5ad
- **Required gene attribute**: "ensembl_id" (Ensembl IDs)
- **Required cell attribute**: "n_counts" (total read counts)
- **Optional**: "filter_pass" for cell filtering
- **Important**: No feature selection - include all detected genes

### Embedding Extraction

```python
from geneformer import EmbExtractor

embex = EmbExtractor(
    model_type="CellClassifier",
    num_classes=n_classes,
    max_ncells=None,
    emb_label=["cell_id"],
    emb_layer=0,  # -1 for last layer, -2 for second-to-last (recommended)
    forward_batch_size=30,
    nproc=8,
)

embs = embex.extract_embs(
    model_directory=model_dir,
    input_data_file="tokenized.dataset",
    output_directory=output_dir,
    output_prefix="embeddings",
)
```

### Classification (Fine-Tuning)

```python
from geneformer import Classifier

cc = Classifier(
    classifier="cell",  # or "gene"
    cell_state_dict={"state_key": "cell_type"},
    filter_data=None,
    training_args=training_args,
    freeze_layers=2,  # Freeze first N layers
    num_crossval_splits=5,
    model_version='V2',
    nproc=4,
    ngpu=1,
)

# Prepare and train
cc.prepare_data(
    input_data_file="tokenized.dataset",
    output_directory=output_dir,
    output_prefix="classifier",
)

cc.train_classifier(
    model_directory=pretrained_model_dir,
    num_classes=num_classes,
    train_data=train_data,
    eval_data=eval_data,
    output_directory=output_dir,
)
```

### Via Helical Package

```python
from helical import Geneformer, GeneformerConfig
import anndata as ad

# Configure
config = GeneformerConfig(batch_size=10)
geneformer = Geneformer(configurer=config)

# Process data
ann_data = ad.read_h5ad("dataset.h5ad")
dataset = geneformer.process_data(ann_data, gene_names="gene_name")

# Get embeddings
embeddings = geneformer.get_embeddings(dataset)
```

---

## V2 Advancements

### Expanded Training
- **3.5x more data**: 104M vs 30M cells
- **Doubled context**: 4,096 vs 2,048 genes per cell
- **Larger models**: Up to 316M parameters

### Continual Learning
- Domain-specific adaptation without forgetting
- Cancer-tuned model for oncology applications
- Learning rate rewarm strategy

### Multi-Task Learning
- Joint learning of multiple cell state features
- Context-specific representations
- Improved generalization

### Quantization
- 4-bit and 8-bit quantized models available
- Resource-efficient fine-tuning via LoRA adapters
- Minimal accuracy loss with significant memory savings
- Enables deployment on consumer hardware

---

## Comparison: Geneformer vs UCE

| Feature | Geneformer | UCE |
|---------|------------|-----|
| **Tokenization** | Rank value encoding | ESM2 protein embeddings |
| **Species** | Human-only (requires orthologs for other species) | Cross-species native |
| **Architecture** | BERT encoder | 33-layer transformer |
| **Parameters** | 10M - 316M | 650M+ |
| **Input representation** | Gene expression ranks | Protein sequences |
| **Zero-shot capability** | Limited | Strong |
| **Fine-tuning strength** | Strong | Available but not primary use |
| **In silico perturbation** | Strong (core feature) | Limited |
| **Batch correction** | Moderate | Strong |
| **Network biology focus** | Primary design goal | Secondary |

### When to Use Each

**Use Geneformer when:**
- Focus is on gene network dynamics
- Need in silico perturbation analysis
- Have human data or can use ortholog conversion
- Want to fine-tune for specific downstream tasks
- Studying drug targets and disease mechanisms

**Use UCE when:**
- Need cross-species analysis
- Want zero-shot universal embeddings
- Focus on cell type annotation at scale
- Need batch correction without fine-tuning
- Working with non-human species

---

## Species Considerations

### Human-Centric Design
- Geneformer was trained exclusively on human data
- Vocabulary based on human genes

### Cross-Species via Ortholog Conversion
- Mouse-Geneformer: 21M mouse cells, dedicated pretraining
- Ortholog mapping enables human Geneformer on mouse data
- ~15% of mouse genes lack human orthologs → information loss
- Performance degrades for complex tasks and distant species

### Mouse-Geneformer
- Developed by separate research group (2024)
- 21 million mouse scRNA-seq profiles
- SiLU activation (vs ReLU in human version)
- Successfully identifies disease-causing genes in mouse models
- Bidirectional: Mouse → human analysis also possible

---

## Attention Weight Analysis

### Network Hierarchy Encoding
Geneformer's attention weights reveal learned biological knowledge:

1. **Early layers**: Survey diverse gene ranks broadly
2. **Later layers**: Focus on highly-ranked, context-critical genes
3. **Transcription factors**: Receive higher attention weights
4. **Central network nodes**: Prioritized by attention mechanism

### Interpretability Applications
- Identify key regulatory genes
- Understand context-dependent gene relationships
- Discover novel network connections

---

## Limitations

### Technical Limitations
1. **Human-centric**: Requires ortholog conversion for other species
2. **GPU requirement**: Efficient usage requires GPU resources
3. **Hyperparameter sensitivity**: No universal defaults; tuning required per task
4. **Zero-shot weakness**: Limited zero-shot performance vs fine-tuned

### Biological Limitations
1. **Cancer exclusion**: Original pretraining excluded cancer cells (addressed by V2 continual learning)
2. **High mutation burden**: Cannot interpret cells with extensive somatic mutations
3. **Static representation**: Doesn't model temporal dynamics
4. **No spatial information**: scRNA-seq only, no spatial context

### Comparative Limitations
1. **Cross-species**: Less elegant than UCE's protein-based approach
2. **Perturbation prediction**: Recent benchmarks show foundation models struggle to outperform PCA on perturbation tasks
3. **Batch effects**: Not as robust as methods specifically designed for batch correction

---

## Ecosystem and Integrations

### HuggingFace Integration
- Models hosted on HuggingFace Model Hub
- Compatible with HuggingFace Transformers library
- BERT fine-tuning objectives readily applicable

### NVIDIA BioNeMo
- Geneformer included in BioNeMo Framework
- Optimized for NVIDIA GPU infrastructure
- Part of Clara suite for drug discovery

### CELLxGENE Census
- Pre-computed Geneformer embeddings available
- Fine-tuned models for cell subclass prediction
- Integration with Census API

### Helical Package
- Unified interface for multiple foundation models
- Standardized embedding extraction
- Comparison tools (Geneformer vs UCE, TranscriptFormer)

### Virtual Cells Platform
- Listed on CZI Virtual Cells Platform
- Standardized model cards and documentation

---

## Future Directions

### From Authors
1. **Larger corpuses**: Continued expansion of pretraining data
2. **Domain-specific variants**: More continual learning for specific tissues/diseases
3. **Network dynamics modeling**: Better capture of temporal gene network changes

### Community Development
1. **Species-specific models**: Mouse-Geneformer as template for other organisms
2. **Multi-omics integration**: Combining with ATAC-seq, spatial data
3. **Improved perturbation prediction**: Addressing benchmark limitations

---

## Resources

### Primary Resources
- **HuggingFace Model**: https://huggingface.co/ctheodoris/Geneformer
- **HuggingFace Dataset**: https://huggingface.co/datasets/ctheodoris/Genecorpus-30M
- **Documentation**: https://geneformer.readthedocs.io/
- **Virtual Cells**: https://virtualcellmodels.cziscience.com/model/geneformer

### Code and Tools
- **GitHub (jkobject fork)**: https://github.com/jkobject/geneformer
- **NVIDIA BioNeMo**: https://nvidia.github.io/bionemo-framework/models/geneformer/
- **Helical Package**: https://pypi.org/project/helical/

### Tutorials
- CELLxGENE Census tutorial for Geneformer prediction
- Helical comparison notebooks (Geneformer vs UCE, TranscriptFormer)
- BioNeMo inference tutorials

---

## Key Takeaways

1. **Pioneer in single-cell foundation models**: One of the first to apply transfer learning to scRNA-seq
2. **Rank value encoding**: Novel tokenization preserving biological relevance
3. **In silico perturbation**: Unique strength for drug target discovery
4. **Experimental validation**: Cardiomyopathy therapeutic targets validated with CRISPR
5. **Network biology focus**: Designed specifically for understanding gene regulatory networks
6. **Transfer learning power**: Excellent performance with limited task-specific data
7. **Human-centric**: Requires adaptation for other species
8. **Fine-tuning recommended**: Zero-shot performance limited; fine-tuning unlocks full potential
9. **Active development**: V2 with 3x more data, larger models, quantization support
10. **Ecosystem integration**: Available across multiple platforms (HuggingFace, BioNeMo, Helical, CELLxGENE)

---

## Citation

```bibtex
@article{theodoris2023transfer,
  title={Transfer learning enables predictions in network biology},
  author={Theodoris, Christina V and Xiao, Ling and Chopra, Anant and Chaffin, Mark D and Al Sayed, Zeina R and Hill, Matthew C and Mantineo, Helene and Brydon, Elizabeth M and Zeng, Zexian and Liu, X Shirley and Ellinor, Patrick T},
  journal={Nature},
  volume={618},
  number={7965},
  pages={616--624},
  year={2023},
  publisher={Nature Publishing Group}
}

@article{chen2024quantized,
  title={Quantized multi-task learning for context-specific representations of gene network dynamics},
  author={Chen, Han and Venkatesh, Maanasa S and Gomez Ortega, Janet and Mahesh, Sai Vishnu and Nandi, Tarun and Madduri, Ravi and Pelka, Karin and Theodoris, Christina V},
  journal={bioRxiv},
  pages={2024--08},
  year={2024},
  publisher={Cold Spring Harbor Laboratory}
}
```

---

*Document compiled: January 2026*
*Research sources: Nature 2023, bioRxiv 2024, HuggingFace, NVIDIA BioNeMo, CELLxGENE Census, Helical documentation, multiple benchmark studies*

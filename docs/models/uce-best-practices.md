# UCE: Universal Cell Embeddings — Comprehensive Research Summary
Overview
UCE is a zero-shot foundation model for single-cell gene expression data developed at Stanford University by Yanay Rosen, Yusuf Roohani, and collaborators (Jure Leskovec lab, Stephen Quake, Tabula Sapiens Consortium). Released November 2023 via bioRxiv preprint, UCE creates a unified biological latent space that can represent any cell from any tissue or species without additional fine-tuning, retraining, or preprocessing.
Key Innovation: UCE uses ESM2 protein language model embeddings (15 billion parameters) to tokenize genes, enabling cross-species generalization to organisms never seen during training.

Architecture Details
Model Specifications:

33-layer transformer with over 650 million parameters
Also available as a 4-layer lightweight version
Gene embeddings from ESM2: 5,120-dimensional vectors per gene
CLS token captures cell-level embedding at final transformer layer

Cell Representation Pipeline:

Gene Sampling: Given gene expression for a single cell, UCE samples genes with replacement, weighted by expression level (only non-zero expressed genes)
Protein Tokenization: Each gene is converted to its protein product's ESM2 embedding
Chromosomal Organization: Gene tokens are sorted by genomic location and grouped by chromosome with special start/end tokens
CLS Token: A special cell-level token is prepended to capture the overall cell embedding
Transformer Processing: The combined representation passes through the 33-layer transformer
Output: Cell embedding is taken from the CLS token at the final layer

Training Objective: Masked binary prediction

During training, 20% of expressed genes are masked
Masked genes are combined with randomly sampled non-expressed genes
Model predicts whether each query gene was expressed (binary classification)
Binary cross-entropy loss function


Training Data
Integrated Mega-scale Atlas (IMA):

36 million cells from 300+ datasets
1,000+ uniquely named cell types
8 species: Human, mouse, zebrafish, mouse lemur, crab-eating macaque, rhesus macaque, tropical clawed frog, pig
33.9 million cells from CZ CELLxGENE (human + mouse)
2.3 million cells from 28 additional datasets

Data Preprocessing:

CELLxGENE datasets: Filter cells by minimum gene counts (200), genes by minimum cells count (10)
No highly variable gene selection applied for CELLxGENE
Other datasets: Top 8,000 highly variable genes (Seurat v3 formula via Scanpy)
Protein-coding genes only (those with available ESM2 embeddings)

Training Infrastructure:

40 days across 24 A100 80GB GPUs
Completely self-supervised (no cell type annotations used)
PyTorch with HuggingFace Accelerator


Key Capabilities
1. Zero-Shot Embedding

No fine-tuning or retraining required for new datasets
No preprocessing (gene selection, cell type identification) needed
Works on any genes from any species, even non-homologs

2. Cross-Species Generalization

Successfully embeds novel species not in training (green monkey, chicken)
Uses protein embeddings to bypass sequence homology requirements
Less effective for evolutionarily distant species (e.g., Drosophila)

3. Batch Effect Correction

Corrects technology-based batch effects (droplet vs plate-based)
Matches or exceeds fine-tuned methods (scVI, scArches) in SCIB benchmarks

4. Cell Type Organization

Emergent clustering by biological function, not experimental artifacts
Embedding distances correlate with Cell Ontology tree distances
Macrophages from different tissues align in UCE space despite distinct transcriptional identities


Benchmark Performance
Single-Cell Integration Benchmark (SCIB):

Outperforms Geneformer by 13.9% overall score
16.2% higher biological conservation score
10.1% higher batch correction score

Cell Type Silhouette Scores (Tabula Sapiens v2):

UCE has highest silhouette score for 67% of cell types
Outperforms Geneformer on 80% of cell types
Outperforms scGPT on 83% of cell types
Outperforms tGPT on 73% of cell types
B cells: 93% higher silhouette score vs scGPT, 25% higher vs Geneformer

vs Fine-Tuned Methods:

UCE zero-shot performs slightly better than scVI and scArches (which require dataset-specific training)

Perturbation Prediction (Important Limitation):

Recent benchmarks (Nature Methods, August 2025) show UCE and other foundation models do not outperform simple linear baselines for gene perturbation effect prediction
UCE was not designed for perturbation prediction; lacks built-in in silico perturbation functionality


Inference Workflow
Basic Embedding Generation:
bashpython eval_single_anndata.py \
  --adata_path {path_to_anndata} \
  --dir {output_dir} \
  --species {species} \
  --model_loc {model_loc} \
  --batch_size {batch_size}
Parameters:

adata_path: h5ad file with scRNA-seq counts in .X, gene names (not ENSEMBLIDs) in .var_names
species: Species of the dataset
model_loc: Path to model weights .torch file
batch_size: 25 for 33-layer model on 80GB GPU; 100 for 4-layer model

Output: AnnData with UCE embeddings in .obsm["X_uce"]
Processing Speed: ~3,000 cells/minute at batch size 128 on A100 GPU (bf16 mode)

Helical Library Integration
UCE is available through the Helical framework for easier use:
Basic Embedding:
pythonfrom helical.models.uce import UCE, UCEConfig
import anndata as ad

configurer = UCEConfig(batch_size=10)
uce = UCE(configurer=configurer)
ann_data = ad.read_h5ad("dataset.h5ad")
data_loader = uce.process_data(ann_data)
embeddings = uce.get_embeddings(data_loader)
Fine-Tuning for Classification (despite official guidance against fine-tuning):
pythonfrom helical.models.uce import UCEConfig, UCEFineTuningModel
import anndata as ad

ann_data = ad.read_h5ad("dataset.h5ad")
label_set = set(cell_types)

configurer = UCEConfig(batch_size=10)
uce_fine_tune = UCEFineTuningModel(
    uce_config=configurer,
    fine_tuning_head="classification",
    output_size=len(label_set)
)

dataset = uce_fine_tune.process_data(ann_data)
# Map labels to integers
class_id_dict = dict(zip(label_set, range(len(label_set))))
cell_types = [class_id_dict[ct] for ct in ann_data.obs.cell_type]

uce_fine_tune.train(train_input_data=dataset, train_labels=cell_types)

Case Study: Norn Cell Discovery
UCE demonstrated its utility through the Norn cell case study:

Background: Norn cells are recently discovered erythropoietin (Epo)-producing fibroblast-like cells in the kidney
Method: Zero-shot embedding of mouse renal cells identified Norn cell clusters
Cross-Tissue Search: Logistic classifier trained on UCE embeddings applied across all 36M cells in IMA
Findings:

Predicted Norn-like cells in gonad, heart, and lung tissues
Identified differences in lung Norn-like cells between COPD, IPF (idiopathic pulmonary fibrosis), and control patients
COPD Norn-like cells express hypoxia pathway genes at higher levels than IPF


Significance: Demonstrates UCE's ability to perform unbiased cross-tissue, cross-species biological discovery


Comparison with Other Models
FeatureUCEscGPTGeneformerSTATEParameters650M~50M~30M600M (SE)Training cells36M33M30M167M obs + 100M pertArchitecture33-layer transformerGPT-style decoderBERT-style encoderBidirectional + set-basedGene tokenizationESM2 protein embeddingsExpression binsRank-basedLog-normalizedZero-shotYes (designed for)Requires fine-tuningRequires fine-tuningRequires fine-tuningCross-speciesYes (via protein embeddings)LimitedLimitedNoPerturbation predictionNot designed forClaims capabilityClaims capabilityPrimary focusPreprocessingNone requiredRequiredRequiredRequired

Limitations & Considerations
Acknowledged Limitations:

Distant Species: Poor performance on evolutionarily distant species (e.g., Drosophila fly cells incorrectly annotated)
Coarse Cell Types: Limited by resolution of expert-annotated cell type labels
No Raw Transcript Information: Does not utilize detailed information in raw RNA transcripts
Perturbation Prediction: Not designed for and does not outperform linear baselines
Black Box: Internal decision-making not easily interpretable

Official Guidance:

UCE is not intended to be fine-tuned (per CZI Virtual Cells Platform)
Compare UCE in zero-shot setting only
Embeddings from 33-layer and 4-layer models are not compatible

Biases/Risks:

May reflect biases in training data
Potential for inaccurate outputs or hallucinations
Risk of incorrect biological interpretations


Resources
GitHub: https://github.com/snap-stanford/UCE
Preprint: bioRxiv 2023.11.28.568918 (v2: October 2024)
Model Weights:

33-layer: Available via Figshare
4-layer: Downloaded automatically
CZI Virtual Cells Platform: https://virtualcellmodels.cziscience.com/model/uce
Processed Datasets: Google Drive link in GitHub repo
Helical Integration: https://helical.readthedocs.io/


License
MIT License - Open for research and commercial use

Key Takeaways

Zero-shot champion: UCE excels at embedding new datasets without any fine-tuning or preprocessing
Cross-species via proteins: ESM2 embeddings enable species-agnostic gene representation
Strong integration: Outperforms competitors on batch correction and biological conservation
Not for perturbation: Unlike STATE, UCE is not designed for predicting perturbation effects
No preprocessing: Uniquely requires no gene selection or normalization
Emergent biology: Discovers biological organization not explicitly trained for
Practical: Fast inference (~3,000 cells/minute), accessible via multiple frameworks
Limitations: Struggles with distant species, limited to transcriptomics, black-box nature

UCE represents a significant advance in creating universal cell representations, but should be used for its strengths (zero-shot embedding, cross-species integration, cell type annotation) rather than perturbation prediction tasks where models like STATE are specifically designed.

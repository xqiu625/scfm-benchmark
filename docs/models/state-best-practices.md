# STATE: Single-Cell Foundation Model — Deep Research Report
Overview
STATE (State Transition and Embedding) is Arc Institute's first-generation virtual cell model, released in June 2025. It is specifically designed to predict how cells respond to perturbations (drugs, cytokines, genetic modifications). STATE is notable for being the first model to consistently beat simple linear baselines Arc Institute in perturbation prediction tasks—a significant achievement given that none of the prior deep learning models was able to consistently outperform the mean prediction or the linear model. Nature
Key distinguishing features:

Two-module architecture (State Embedding + State Transition)
Set-based transformer operating on cell populations rather than individual cells
Trained on 270+ million cells (167M observational + 100M+ perturbational)
First model to consistently outperform linear baselines for perturbation prediction


Model Architecture
Two-Module Design
STATE consists of two interlocking modules:
1. State Embedding (SE) Model
The SE model converts transcriptome data into a smooth multidimensional vector space:

Architecturally, the SE encoder is a dense, bidirectional transformer trained to predict log-normalized gene expression. bioRxiv
The SE decoder is a smaller, specialized MLP that predicts gene expression from a combination of the learned cell embedding and the target gene embedding. bioRxiv
This architectural asymmetry encourages the learning of generalizable representations of cell state in a single vector embedding. bioRxiv

Dual-axis loss function:

SE is trained with a loss computed along two axes: it predicts expression across genes within each cell, and for each gene across cells in each minibatch. bioRxiv
This dual-axis formulation encourages the model to capture relative variation in gene expression both within individual cells and across the population. bioRxiv

2. State Transition (ST) Model
The ST model predicts how cells transition between states after perturbation:

This model is built on a bidirectional transformer architecture that uniquely leverages self-attention over sets of cells. Arc Institute
ST can flexibly capture biological and technical heterogeneity (such as cell cycle state or biases in RNA-seq data) without relying on explicit distributional assumptions. Arc Institute


Training Data
STATE is trained on the largest perturbation dataset to date:
Training Datasets

Tahoe-100M: A giga-scale single-cell atlas of 100 million transcriptomic profiles measuring how each of 1,100 small-molecule perturbations impact cells across 50 cancer cell lines. bioRxiv
Parse-PBMC: Additional perturbational data from blood cells
Replogle-Nadig: CRISPR-based genetic perturbation data
Observational data: ~167 million cells from scBaseCount and CELLxGENE

Total scale:

State is trained on observational data from nearly 170 million cells and perturbational data from over 100 million cells across 70 cell lines. Arc Institute


Installation & Setup
Requirements

Python 3.10+
NVIDIA GPU recommended
Package manager: uv (from Astral)

Installation from PyPI
bashuv tool install arc-state
Installation from Source
bashgit clone git@github.com:ArcInstitute/state.git
cd state
uv run state
For Development
bashgit clone git@github.com:ArcInstitute/state.git
cd state
uv tool install -e .
CLI Access
bashstate --help
# Output: usage: state [-h] {emb,tx}

Pretrained Models
Available on HuggingFace (arcinstitute/):
ModelParametersUse CaseSE-600M600MCell embedding generationST (Tahoe-100M trained)-Perturbation prediction

Data Preprocessing
Standard Pipeline
Training Data Preprocessing:
bashstate tx preprocess_train \
  --adata /path/to/raw_data.h5ad \
  --output /path/to/preprocessed_training_data.h5ad \
  --num_hvgs 2000
This performs:

Total counts normalization (sc.pp.normalize_total)
Log1p transformation (sc.pp.log1p)
Highly variable gene selection (sc.pp.highly_variable_genes)
Storage in .obsm['X_hvg']

Inference Data Preprocessing:
bashstate tx preprocess_infer \
  --adata /path/to/real_data.h5ad \
  --output /path/to/control_template.h5ad \
  --control_condition "DMSO" \
  --pert_col "treatment" \
  --seed 42
Input Requirements

Format: AnnData (.h5ad)
Matrix format: CSR sparse matrix
Required columns: gene_name in var dataframe


Task-Specific Workflows
1. Perturbation Prediction (State Transition)
Configuration via TOML files:
toml# Example: fewshot.toml
[datasets]
replogle = "/path/to/replogle/dataset/"

[training]
replogle = "train"

[zeroshot]
"replogle.jurkat" = "test"
"replogle.k562" = "val"

[fewshot]
[fewshot."replogle.rpe1"]
val = ["AARS", "TUFM"]
test = ["NUP107", "RPUSD4"]
Training Command:
bashstate tx train \
  data.kwargs.toml_config_path="examples/fewshot.toml" \
  data.kwargs.embed_key=X_hvg \
  data.kwargs.num_workers=12 \
  data.kwargs.batch_col=batch_var \
  data.kwargs.pert_col=target_gene \
  data.kwargs.cell_type_key=cell_type \
  data.kwargs.control_pert=TARGET1 \
  training.max_steps=40000 \
  training.val_freq=100 \
  training.ckpt_every_n_steps=100 \
  training.batch_size=8 \
  training.lr=1e-4 \
  model.kwargs.cell_set_len=64 \
  model.kwargs.hidden_dim=328 \
  model=pertsets \
  output_dir="$HOME/state" \
  name="test"
Key Hyperparameters:
ParameterDefault/RecommendedDescriptiontraining.max_steps40000Total training stepstraining.batch_size8Batch sizetraining.lr1e-4Learning ratemodel.kwargs.cell_set_len64-256Cells per set (256 optimal for Tahoe-100M)model.kwargs.hidden_dim328Hidden dimension
Evaluation:
bashstate tx predict --output-dir $HOME/state/test/ --checkpoint final.ckpt
Inference on new data:
bashstate tx infer \
  --output $HOME/state/test/ \
  --output_dir /path/to/model/ \
  --checkpoint /path/to/model/final.ckpt \
  --adata /path/to/anndata/processed.h5 \
  --pert_col gene \
  --embed_key X_hvg
2. Cell Embedding (State Embedding)
Generate embeddings:
bashstate emb transform \
  --model-folder /path/to/SE-600M \
  --checkpoint /path/to/SE-600M/se600m_epoch15.ckpt \
  --input /path/to/data.h5ad \
  --output /path/to/output.h5ad
Train new SE model:
bashstate emb fit --conf ${CONFIG}
3. Vector Database Integration
Build database:
bashuv tool install ".[vectordb]"

state emb transform \
  --model-folder /path/to/SE-600M \
  --input /path/to/data.h5ad \
  --lancedb tmp/state_embeddings.lancedb \
  --gene-column gene_symbols
Query similar cells:
bashstate emb query \
  --lancedb tmp/state_embeddings.lancedb \
  --input tmp/query.h5ad \
  --output tmp/similar_cells.csv \
  --k 3
```

---

## Evaluation Paradigms

STATE supports two evaluation modes:

### Zero-shot Evaluation
- Tests generalization to completely unseen cell types
- Entire cell lines held out from training

### Few-shot Evaluation
- Tests generalization with limited examples
- Specific perturbations held out within cell types

---

## Benchmarking & Performance

### Key Results

During benchmarking on Tahoe-100M, State demonstrated a 50 percent improvement in distinguishing perturbation effects and achieved twice the accuracy in identifying true differentially expressed genes compared to existing models. 

### Comparison with Prior Models

Prior foundation models struggled against simple baselines:
- Surprisingly, we found that foundation models generally underperformed compared to a simple baseline model that uses the mean of the training samples. 
- In the differential expression space (Pearson Delta), even the simplest baseline model, Train Mean reached better correlation than foundation models (scGPT, scFoundation). 

STATE addresses this limitation by being the first to consistently outperform linear baselines.

### Evaluation Metrics

STATE uses the **cell-eval** framework (https://github.com/ArcInstitute/cell-eval):
- Pearson ΔC: Correlation between predicted/observed mean expression deltas
- DEG identification accuracy
- Perturbation separability in embedding space

---

## Best Practices

### Data Preparation
1. Use CSR sparse matrix format for h5ad files
2. Ensure `gene_name` column exists in var dataframe
3. For control matching, specify control condition (e.g., "DMSO")
4. Use 2000 HVGs as default for preprocessing

### Training Recommendations
1. **Cell set size**: Use 256 cells per set for optimal performance on large datasets
2. **Self-attention**: Critical for performance—removing it substantially degrades results
3. **Validation frequency**: Check every 100 steps for early detection of issues
4. **Checkpoint saving**: Save regularly (every 100-1000 steps)

### Evaluation Guidelines
1. Use held-out PERTURBATIONS, not held-out cells from training perturbations
2. Always compare against linear baselines
3. Check for mode collapse using rank metrics
4. Test zero-shot cross-cell-type generalization

### Configuration Tips
- Dataset paths should point to directories containing h5ad files
- Use format `"dataset_name.cell_type"` for cell type specification
- Ensure control condition is available across all splits
- No conflicts between zeroshot and fewshot assignments

---

## Limitations

### Current Scope
- While Tahoe-100M is broad and deep within the domain of in vitro cancer cell line models, it does not yet encompass the full range of primary cell types, tissue microenvironments, or in vivo perturbations. 

### Biological Validation
- Model released June 2025—independent validation studies still emerging
- Performance on rare cell types and novel perturbations needs evaluation

### Technical Constraints
- Requires GPU for efficient training/inference
- Large memory footprint for SE-600M model
- Non-commercial license restricts certain applications

---

## Common Issues & Solutions

| Issue | Solution |
|-------|----------|
| Missing gene names | Ensure `gene_name` column in adata.var |
| Memory errors | Reduce cell_set_len, use smaller batches |
| Slow training | Increase num_workers, verify GPU utilization |
| Poor zero-shot | Train on more diverse perturbations |
| Mode collapse | Check rank metrics, adjust architecture |

---

## Associated Tools & Resources

### Official Repositories
- **STATE**: https://github.com/ArcInstitute/state
- **cell-eval**: https://github.com/ArcInstitute/cell-eval
- **cell-load**: https://github.com/ArcInstitute/cell-load

### Pretrained Models
- HuggingFace: https://huggingface.co/arcinstitute

### Tutorials (Google Colab)
- Train ST for genetic perturbation: [Colab](https://colab.research.google.com/drive/1Ih-KtTEsPqDQnjTh6etVv_f-gRAA86ZN)
- Inference with Tahoe-100M trained model: [Colab](https://colab.research.google.com/drive/1bq5v7hixnM-tZHwNdgPiuuDo6kuiwLKJ)
- Embed and annotate with SE: [Colab](https://colab.research.google.com/drive/1uJinTJLSesJeot0mP254fQpSxGuDEsZt)
- Virtual Cell Challenge: [Colab](https://colab.research.google.com/drive/1QKOtYP7bMpdgDJEipDxaJqOchv7oQ-_l)

### Data Resources
- Arc Virtual Cell Atlas: https://arcinstitute.org/tools/virtualcellatlas
- Tahoe-100M preprint: https://www.biorxiv.org/content/10.1101/2025.02.20.639398v2

---

## Licensing

- **Code**: CC BY-NC-SA 4.0 (Creative Commons Attribution-NonCommercial-ShareAlike)
- **Model weights**: Arc Research Institute State Model Non-Commercial License
- **Usage**: Subject to Arc Research Institute State Model Acceptable Use Policy

---

## Key Takeaways

1. **STATE is purpose-built for perturbation prediction**—unlike general-purpose scFMs like scGPT or Geneformer
2. **Set-based architecture** captures population-level heterogeneity
3. **First to beat linear baselines** consistently—a critical benchmark
4. **Large training data** (270M+ cells) with emphasis on perturbational data
5. **Two-module design** allows flexible use of SE embeddings or full ST predictions
6. **TOML-based configuration** makes experimental setup straightforward
7. **Zero-shot and few-shot** evaluation paradigms built into the framework
8. **Non-commercial license** limits deployment in commercial settings

---

## Citation

For publications using STATE:
```
Preprint: https://arcinstitute.org/manuscripts/State
```

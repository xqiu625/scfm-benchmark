# STATE Fine-tuning Guide

---

## Best Practices for Fine-Tuning the State Model (Arc Institute)

### Overview
State is designed to predict how various stem cells, cancer cells, and immune cells respond to drugs, cytokines, or genetic perturbations. It consists of two interlocking modules:
- **State Embedding (SE) model**: Converts transcriptome data into a smooth multidimensional vector space
- **State Transition (ST) model**: Predicts how cells transition between states in response to perturbations

### Installation & Setup
From the [GitHub repository](https://github.com/ArcInstitute/state):

```bash
# Install via uv
git clone git@github.com:ArcInstitute/state.git
cd state
uv tool install -e .
```

### Data Preprocessing (Critical First Step)

Prepares training h5ad files by normalizing counts, applying log1p, and selecting highly variable genes (HVGs). The HVG matrix is stored in .obsm["X_hvg"].

```bash
state tx preprocess_train \
    --adata /path/to/raw_data.h5ad \
    --output /path/to/preprocessed_training_data.h5ad \
    --num_hvgs 2000
```

### Training Configuration

**Key parameters for fine-tuning:**
```bash
state tx train \
    data.kwargs.toml_config_path=examples/fewshot.toml \
    data.kwargs.embed_key=X_hvg \
    data.kwargs.pert_col=target_gene \
    data.kwargs.cell_type_key=cell_type \
    training.max_steps=40000 \
    training.batch_size=8 \
    model=state \
    output_dir="$HOME/state" \
    name="test"
```

### Best Practices

1. **Data Splits Configuration**
   - Use TOML config files for defining zero-shot, few-shot, and training splits
   - Use the "dataset.cell_type" format in [zeroshot] and [fewshot]. Anything not listed defaults to training.

2. **Embedding Preparation**
   - Create embedding profiles before training using `state emb`
   - The system auto-detects gene-name fields (checks `_index`, `gene_name`, `gene_symbols`, `feature_name`, `gene_id`, `symbol`)

3. **Evaluation**
   - Use the Cell-Eval framework for biologically relevant metrics
   - Three key metrics: differential expression score, perturbation discrimination score, and mean absolute error (MAE).

4. **Inference on New Data**
```bash
state tx infer \
    --model-dir /path/to/run \
    --checkpoint /path/to/run/checkpoints/final.ckpt \
    --adata /path/to/preprocessed_data.h5ad \
    --pert-col gene \
    --embed-key X_hvg \
    --output predictions.h5ad
```

### Training Resources
- **Google Colab tutorial**: Available for Virtual Cell Challenge training
- **Arc Virtual Cell Atlas**: Over 330 million cells for training data (gs://arc-institute-virtual-cell-atlas)
- **Tahoe-100M dataset**: Large-scale perturbation data across 70 cell lines

### Key Considerations

1. **Scale matters**: State is trained on observational data from nearly 170 million cells and perturbational data from over 100 million cells across 70 cell lines.

2. **Context generalization**: The model is designed to predict perturbation effects in cellular contexts not seen during training

3. **Licensing**: Code is CC BY-NC-SA 4.0; model weights have a non-commercial license

I'll search for more detailed information on hyperparameters, dataset preparation, and configuration specifics.Let me also fetch the bioRxiv preprint for more detailed hyperparameter information:Let me also look at the cell-eval repository for evaluation best practices:Now I have comprehensive information. Here's a detailed guide on **Best Practices for Fine-Tuning the STATE Model**:

---

## Comprehensive Guide: Fine-Tuning Arc Institute's STATE Model

### 1. Architecture Overview

State consists of two interlocking modules:
- **State Embedding (SE) Model**: Converts transcriptomes into smooth vector representations, trained on 167 million human cells
- **State Transition (ST) Model**: Predicts perturbation effects using bidirectional transformer with self-attention over sets of cells

---

### 2. Installation & Setup

```bash
# Install via PyPI
uv tool install arc-state

# Or from source
git clone git@github.com:ArcInstitute/state.git
cd state
uv tool install -e .

# Related packages
uv pip install -U cell-eval   # Evaluation framework
uv pip install -U cell-load   # Dataloaders
```

---

### 3. Data Preprocessing (Critical First Step)

#### Training Data Preprocessing
```bash
state tx preprocess_train \
  --adata /path/to/raw_data.h5ad \
  --output /path/to/preprocessed_training_data.h5ad \
  --num_hvgs 2000
```

This performs:
- Normalizes total counts per cell (`sc.pp.normalize_total`)
- Applies log1p transformation (`sc.pp.log1p`)
- Identifies highly variable genes (`sc.pp.highly_variable_genes`)
- Stores the HVG expression matrix in `.obsm['X_hvg']`

#### Inference Data Preprocessing
```bash
state tx preprocess_infer \
  --adata /path/to/real_data.h5ad \
  --output /path/to/control_template.h5ad \
  --control_condition "DMSO" \
  --pert_col "treatment" \
  --seed 42
```

---

### 4. Configuration File Setup (TOML)

Create a TOML configuration file to define datasets and splits:

```toml
# Dataset paths
[datasets]
replogle = "/path/to/replogle_dataset/"
jurkat = "/path/to/jurkat_dataset/"

# Training specifications
[training]
replogle = "train"
jurkat = "train"

# Zeroshot: Hold out entire cell types
[zeroshot]
"replogle.jurkat" = "test"
"jurkat.rpe1" = "val"

# Fewshot: Hold out specific perturbations
[fewshot]
[fewshot."replogle.rpe1"]
val = ["AARS"]
test = ["AARS", "NUP107", "RPUSD4"]  # can overlap with val

[fewshot."jurkat.k562"]
val = ["GENE1", "GENE2"]
test = ["GENE3", "GENE4"]
```

---

### 5. Training Hyperparameters

#### Full Training Command with All Parameters

```bash
state tx train \
    # Data configuration
    data.kwargs.toml_config_path="examples/fewshot.toml" \
    data.kwargs.embed_key=X_hvg \
    data.kwargs.num_workers=12 \
    data.kwargs.batch_col=batch_var \
    data.kwargs.pert_col=target_gene \
    data.kwargs.cell_type_key=cell_type \
    data.kwargs.control_pert=non-targeting \
    
    # Training hyperparameters
    training.max_steps=40000 \
    training.val_freq=100 \
    training.ckpt_every_n_steps=100 \
    training.batch_size=8 \
    training.lr=1e-4 \
    
    # Model architecture
    model.kwargs.cell_set_len=64 \
    model.kwargs.hidden_dim=328 \
    model=pertsets \
    
    # Logging
    wandb.tags="[test]" \
    output_dir="$HOME/state" \
    name="my_experiment"
```

#### Recommended Hyperparameter Ranges

| Parameter | Recommended Value | Notes |
|-----------|-------------------|-------|
| `learning_rate` | `1e-4` | Standard starting point |
| `batch_size` | `8` | Adjust based on GPU memory |
| `max_steps` | `40000` | For full training |
| `val_freq` | `100` | Validation frequency |
| `cell_set_len` | `64` | Cells per set in transformer |
| `hidden_dim` | `328` | Transformer hidden dimension |
| `num_hvgs` | `2000` | Highly variable genes to select |
| `num_workers` | `12-24` | Data loading parallelism |

---

### 6. Data Loading Parameters

From the [cell-load documentation](https://github.com/ArcInstitute/cell-load):

```python
from cell_load.data_modules import PerturbationDataModule

dm = PerturbationDataModule(
    # Required parameters
    toml_config_path="/path/to/config.toml",
    embed_key="X_hvg",
    num_workers=24,
    batch_col="gem_group",
    pert_col="gene",
    cell_type_key="cell_type",
    control_pert="non-targeting",
    
    # Optional parameters
    barcode=True,
    perturbation_features_file="/path/to/gene_embeddings.pt",
    output_space="gene",
    basal_mapping_strategy="random",
    n_basal_samples=1,
    should_yield_control_cells=True,
    batch_size=128,
)
```

---

### 7. Model Evaluation with Cell-Eval

#### Evaluation Command
```bash
cell-eval run \
    -ap /path/to/predictions.h5ad \
    -ar /path/to/ground_truth.h5ad \
    --num-threads 64 \
    --profile full
```

#### Python API
```python
from cell_eval import MetricsEvaluator

evaluator = MetricsEvaluator(
    adata_pred=adata_pred,
    adata_real=adata_real,
    control_pert="control",
    pert_col="perturbation",
    num_threads=64,
)
results, agg_results = evaluator.compute()
```

#### Key Metrics
- **Differential Expression Score**: How well the model identifies true DE genes
- **Perturbation Discrimination Score**: How well the model distinguishes perturbation effects
- **Mean Absolute Error (MAE)**: Expression-level accuracy across all genes

---

### 8. Training Strategies

#### Zeroshot Evaluation (Unseen Cell Types)
```toml
[zeroshot]
"replogle.jurkat" = "test"     # All jurkat cells to test
"replogle.rpe1" = "val"        # All rpe1 cells to validation
```

#### Fewshot Evaluation (Limited Perturbation Examples)
```toml
[fewshot]
[fewshot."replogle.k562"]
val = ["AARS"]
test = ["NUP107", "RPUSD4"]
```

#### Mixed Strategy
Combine both for comprehensive evaluation of generalization ability.

---

### 9. Inference on New Data

```bash
state tx infer \
    --model-dir /path/to/trained_model/ \
    --checkpoint /path/to/checkpoints/final.ckpt \
    --adata /path/to/preprocessed_data.h5ad \
    --pert-col gene \
    --embed-key X_hvg \
    --output predictions.h5ad
```

---

### 10. Best Practices Summary

1. **Data Quality**:
   - Ensure CSR matrix format for h5ad files
   - Include `gene_name` in the `var` dataframe
   - Normalize and log-transform before training

2. **Split Design**:
   - Control cells must be available across all splits
   - Control cell mapping is done within the same file only

3. **Training Tips**:
   - Start with recommended hyperparameters
   - Use validation frequency of 100 steps
   - Checkpoint regularly (`ckpt_every_n_steps=100`)

4. **Evaluation**:
   - Always evaluate using Cell-Eval for standardized metrics
   - Compare against baselines (State outperforms simple linear models)

---

### 11. Resources

- **Google Colab Tutorials**:
  - [Train ST on Replogle-Nadig](https://colab.research.google.com/drive/1Ih-KtTEsPqDQnjTh6etVv_f-gRAA86ZN)
  - [Inference with Tahoe-100M checkpoint](https://colab.research.google.com/drive/1bq5v7hixnM-tZHwNdgPiuuDo6kuiwLKJ)
  - [Embed with SE model](https://colab.research.google.com/drive/1uJinTJLSesJeot0mP254fQpSxGuDEsZt)
  - [Virtual Cell Challenge training](https://colab.research.google.com/drive/1QKOtYP7bMpdgDJEipDxaJqOchv7oQ-_l)

- **Data**:
  - Arc Virtual Cell Atlas: `gs://arc-institute-virtual-cell-atlas`
  - Tahoe-100M dataset (100M+ perturbed cells)

- **Pre-trained Models**: Available on [HuggingFace](https://huggingface.co/arcinstitute)

Would you like me to elaborate on any specific aspect, such as troubleshooting common issues or optimizing for specific hardware configurations?

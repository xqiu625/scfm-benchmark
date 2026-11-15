# scGPT Best Practices

---

## **scGPT Single-Cell Foundation Model: Best Practices Guide**

### **Overview**
scGPT is a foundation model for single-cell biology, based on a generative pretrained transformer across a repository of over 33 million cells. It supports multiple downstream tasks including cell type annotation, batch integration, multi-omic integration, perturbation prediction, and gene network inference.

---

### **1. Installation & Setup**

**Hardware Requirements:**
- Using scGPT requires an NVIDIA GPU with an up-to-date GPU driver
- Supported GPUs: Ampere (A100, RTX 3090), Ada (RTX 4090), or Hopper (H100) architectures for FlashAttention
- scGPT works with Python >= 3.7.13 and R >=3.6.1

**Installation Options:**

```bash
# Basic installation
pip install scgpt "flash-attn<1.0.5"  # flash-attn is optional but recommended

# For orbax package issues (as of 2023.09):
pip install scgpt "flash-attn<1.0.5" "orbax<0.1.8"
```

**Note:** Flash-attention is now optional. The pretrained weights can be loaded on pytorch CPU, GPU, and flash-attn backends using the same load_pretrained function.

**Docker alternative:**
```bash
docker pull xueerchen/scgpt:0.1.7
docker run --gpus all --rm -it xueerchen/scgpt:0.1.7 bash
```

---

### **2. Choosing the Right Pretrained Model**

The whole-human model is recommended for most applications by default. If your fine-tuning dataset shares similar cell type context with the training data of the organ-specific models, these models can usually demonstrate competitive performance as well.

Available checkpoints:
- **Whole-human model**: General purpose, recommended default
- **Organ-specific models**: Blood, kidney, etc. - use when your data matches the organ
- **Continual pretrained model (scGPT_CP)**: Better for zero-shot cell embedding tasks with batch effects

---

### **3. Data Preprocessing Pipeline**

The standard preprocessing follows these steps:

```python
from scgpt.preprocess import Preprocessor

preprocessor = Preprocessor(
    use_key="X",                    # Key for raw data in adata.layers
    filter_gene_by_counts=3,        # Step 1: Filter genes
    filter_cell_by_counts=False,    # Step 2: Filter cells
    normalize_total=1e4,            # Step 3: Depth normalization
    result_normed_key="X_normed",
    log1p=data_is_raw,              # Step 4: Log1p transformation
    result_log1p_key="X_log1p",
    subset_hvg=n_hvg,               # Step 5: Select highly variable genes
    hvg_flavor="seurat_v3" if data_is_raw else "cell_ranger",
    binning=n_bins,                 # Step 6: Value binning (typically 51 bins)
    result_binned_key="X_binned"
)

preprocessor(adata, batch_key="batch")
```

**Key preprocessing notes:**
- The standard pipeline includes depth normalization, log normalization, and highly variable gene (HVG) selection. Value binning is introduced to convert expression counts into relative values.
- Use `seurat_v3` flavor for raw count data
- Typical HVG count: 1200-1800 genes depending on task

---

### **4. Task-Specific Best Practices**

#### **Cell Type Annotation**

The fine-tuning pipeline summarized in the following steps can be used as a general recipe: 1. Specify hyper-parameter setup 2. Load and pre-process data 3. Load the pre-trained scGPT model 4. Finetune scGPT with task-specific objectives 5. Evaluate fine-tuned scGPT

**Recommended hyperparameters:**
- Learning rate: 0.0001, decaying by 10% after each epoch
- Batch size: 64
- Epochs: Typically 10-15
- Max sequence length: n_hvg + 1

#### **Batch Integration**

For batch integration, in addition to the default gene modelling objectives, turn on ESC, DAR and DSBN objectives specifically to facilitate batch integration.

Key settings:
```python
hyperparameter_defaults = dict(
    GEPC=True,              # Gene expression modelling for cell objective
    ecs_thres=0.8,          # Elastic cell similarity threshold
    dab_weight=1.0,         # Domain adversarial batch correction weight
    mask_ratio=0.4,
    n_bins=51,
    lr=0.0001,
    batch_size=64,
    dropout=0.2,
    DSBN=True               # Domain-specific batch normalization
)
```

#### **Perturbation Prediction**

For genetic perturbation response prediction:
- Use the GEARS framework integration
- The model can predict effects of single-gene and combination perturbations on gene expression profiles
- Fine-tune for 15+ epochs for perturbation tasks

#### **Gene Network Inference (GRN)**

Use the pre-trained blood model or task-specific models for GRN inference, particularly regarding gene program extraction and network visualization.

The GRN workflow extracts gene-gene interaction networks from learned attention weights and gene embeddings.

#### **Multi-omic Integration**

scGPT is the only foundation model that can integrate multiple single-cell omics including scATAC-seq data.

For multi-omic integration, scGPT uses condition tokens to encode modality information (scRNA-seq vs scATAC-seq).

---

### **5. Zero-Shot vs Fine-tuned Applications**

**Critical consideration:** Recent benchmarks reveal important limitations:

scGPT and Geneformer in zero-shot configurations perform inconsistently compared to cell embeddings derived from HVG or generated using scVI or Harmony.

**Recommendations:**
- **For exploratory analysis with unknown cell types:** Consider using scVI or Harmony as baselines; don't rely solely on zero-shot scGPT
- **For defined tasks with training data available:** Fine-tuning provides substantial improvements
- **Use the continual pretrained checkpoint (scGPT_CP)** for better zero-shot cell embedding performance

Zero-shot evaluation of these models exposes vulnerabilities that are not evident if evaluated with fine-tuning alone.

---

### **6. Validation & Quality Control**

The model's performance may degrade when analyzing cell types, tissues, or species not well represented in the training data. Users should validate model outputs against independent datasets to mitigate biases and inaccuracies. It is advised to use the model in conjunction with expert biological knowledge.

**Best practices:**
- Always validate against known marker genes
- Compare with established methods (scVI, Harmony, Seurat) on your specific dataset
- Be cautious with rare cell types or underrepresented tissues
- Use multiple evaluation metrics (ASW, NMI, ARI for clustering; avgBIO, avgBATCH for integration)

---

### **7. Training Tips & Hyperparameters**

The model is optimized using the Adam optimizer with a mini-batch size of 512. The initial learning rate is 0.0001 with a weight decay of 0.9 applied after each epoch.

**Fine-tuning recommendations:**
- Start with pretrained weights (faster convergence, better results)
- Use wandb for logging and visualization
- Enable mixed precision training (AMP) for speed
- Set appropriate dropout (0.2 is typical)

**For large datasets:**
- Use `per_seq_batch_sample = True` for balanced batch sampling
- Consider gradient accumulation if memory-limited

---

### **8. Common Issues & Solutions**

| Issue | Solution |
|-------|----------|
| Flash-attn installation fails | Use CUDA 11.7, try `flash-attn<1.0.5`, or skip flash-attn entirely |
| Out of memory | Reduce batch_size, reduce n_hvg, enable AMP |
| Poor zero-shot performance | Fine-tune on your specific task; compare with baselines |
| Gene vocabulary mismatch | Use the provided vocab.json; check gene_info.csv for ENSEMBL conversion |

---

### **9. Online Resources & Tools**

- **Documentation**: https://scgpt.readthedocs.io
- **GitHub**: https://github.com/bowang-lab/scGPT
- **Web apps**: Reference mapping, cell annotation, and GRN inference available at Superbio.ai
- **Virtual Cells Platform**: https://virtualcellmodels.cziscience.com/model/scgpt

A comprehensive end-to-end protocol for fine-tuning scGPT for cell-type classification is available, demonstrating the model's efficiency in handling complex data and improving annotation accuracy achieving 99.5% F1-score on retinal data.

---

### **Key Takeaways**

1. **Fine-tuning is essential** for optimal performance - don't rely on zero-shot alone
2. **Use whole-human model** as default; organ-specific when data matches
3. **Proper preprocessing** (normalization, HVG selection, binning) is critical
4. **Enable task-specific objectives** (ECS, DAB, DSBN) for batch integration
5. **Always validate** against established methods and biological knowledge
6. **Be aware of limitations** with rare cell types and underrepresented tissues

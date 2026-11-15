# Data Processing Pipeline: COVID Platelet Foundation Model Benchmark

**Dataset:** COVID-19 Platelet scRNA-seq Data
**Cells:** ~47,808 cells
**Batches:** 11 (Data_NO)
**Models:** STATE, UCE, Geneformer, TranscriptFormer, scGPT, scGPT (Best Practices)

---

## 1. Data Ingestion

### Input Data
```
File: human_platelet_covid_severity_raw_counts.h5ad
Location: /bigdata/godziklab/shared/Xinru/302006/01_DATA/raw/
```

### Data Sources (11 Batches)

| Batch ID | Study | PMID |
|----------|-------|------|
| 302004data01 | Wilk et al. | 32514174 |
| 302004data02 | Lee et al. | 32651212 |
| 302004data04 | Schulte-Schrepping et al. | 32810438 |
| 302004data05 | Arunachalam et al. | 32788292 |
| 302004data15 | Ren et al. | 33657410 |
| 302004data16 | Liu et al. | 33713619 |
| 302004data17 | Stephenson et al. | 33879890 |
| 302004data18 | Jiang et al. | 32554932 |
| 302005data01 | Bernardes et al. | 33296687 |
| 302005data02 | Combes et al. | 33494096 |
| 302005data03 | Qiu et al. | 34558746 |

*Note: 302005data04 (Mistry et al., PMID: 31754027) contains SLE data and is excluded from COVID analysis.*

### Metadata Columns
| Column | Description |
|--------|-------------|
| `covid_severity` | Severity labels: healthy, mild, moderate, severe, fatal, recovered |
| `Data_NO` | Batch identifier (11 unique batches) |
| `Sample_ID` | Individual sample identifier |
| `Disease` | Disease condition |

### Loading Code
```python
import scanpy as sc
adata = sc.read_h5ad("human_platelet_covid_severity_raw_counts.h5ad")
```

---

## 2. Quality Control (QC)

### Step 2.1: Calculate QC Metrics
```python
sc.pp.calculate_qc_metrics(adata, inplace=True)
```
**Output columns:**
- `n_genes_by_counts`: Number of genes detected per cell
- `total_counts`: Total UMI counts per cell
- `pct_counts_mt`: Percentage of mitochondrial reads

### Step 2.2: Cell Filtering
```python
# Filter by gene count
sc.pp.filter_cells(adata, min_genes=200)   # Remove cells with <200 genes
sc.pp.filter_cells(adata, max_genes=5000)  # Remove cells with >5000 genes (doublets)

# Filter by mitochondrial content
adata = adata[adata.obs['pct_counts_mt'] < 20]  # Remove cells with >20% MT reads
```

### Step 2.3: Gene Filtering
```python
sc.pp.filter_genes(adata, min_cells=10)    # Remove genes in <10 cells
sc.pp.filter_genes(adata, min_counts=3)    # Remove genes with <3 total counts
```

### QC Summary
| Filter | Threshold | Purpose |
|--------|-----------|---------|
| min_genes | 200 | Remove low-quality/empty cells |
| max_genes | 5,000 | Remove doublets |
| pct_counts_mt | <20% | Remove dying/stressed cells |
| min_cells | 10 | Remove rarely detected genes |
| min_counts | 3 | Remove noise genes |

---

## 3. Initial Processing

### Step 3.1: Store Raw Counts
```python
adata.layers["counts"] = adata.X.copy()  # Preserve raw counts
```

### Step 3.2: Depth Normalization
```python
sc.pp.normalize_total(adata, target_sum=1e4)
```
**Purpose:** Normalize each cell to 10,000 total counts to account for sequencing depth differences.

**Why 1e4?** Standard practice in scRNA-seq; makes expression values interpretable as "counts per 10,000".

### Step 3.3: Log Transformation
```python
sc.pp.log1p(adata)
adata.layers["X_normed"] = adata.X.copy()  # Store normalized data
```
**Purpose:**
- Variance stabilization (high-expression genes don't dominate)
- Makes data more normally distributed
- `log1p` = log(x + 1) to handle zeros

---

## 4. Feature Engineering (Crucial Steps)

### Step 4.1: Highly Variable Gene (HVG) Selection

**This is the KEY step for batch effect removal!**

```python
sc.pp.highly_variable_genes(
    adata,
    n_top_genes=1200,           # Select top 1,200 genes
    flavor='seurat_v3',         # Method: Seurat v3 algorithm
    layer='counts',             # Use raw counts for selection
    batch_key='Data_NO'         # ← CRITICAL: Batch-aware selection!
)

# Subset to HVGs only
adata = adata[:, adata.var.highly_variable].copy()
```

#### Why Batch-Aware HVG Selection Matters

| Without `batch_key` | With `batch_key='Data_NO'` |
|---------------------|----------------------------|
| Selects genes variable across ALL data | Selects genes variable WITHIN each batch |
| May select batch-effect genes | Removes batch-specific variation |
| Technical variation dominates | Biological variation preserved |
| → scGPT original (iLISI=1.43) | → scGPT Best Practices (iLISI=3.77) |

#### HVG Selection Parameters
| Parameter | Value | Reason |
|-----------|-------|--------|
| `n_top_genes` | 1,200 | Balance between signal and noise |
| `flavor` | 'seurat_v3' | Best for batch-aware selection |
| `layer` | 'counts' | Use raw counts, not normalized |
| `batch_key` | 'Data_NO' | Enable batch-aware selection |

### Step 4.2: Value Binning (scGPT-specific)

```python
n_bins = 51  # Number of expression bins

X_binned = np.zeros(X.shape, dtype=np.int32)
for i in range(X.shape[1]):  # For each gene
    gene_expr = X[:, i]
    non_zero_mask = gene_expr > 0

    if non_zero_mask.sum() > n_bins:
        # Quantile binning for genes with enough non-zero values
        bins = np.quantile(gene_expr[non_zero_mask], np.linspace(0, 1, n_bins))
        bins = np.unique(bins)
        X_binned[:, i] = np.digitize(gene_expr, bins)
    else:
        # Simple scaling for sparse genes
        X_binned[:, i] = (gene_expr / gene_expr.max() * (n_bins - 1)).astype(int)

adata.layers["X_binned"] = X_binned
```

**Purpose:** Convert continuous expression values to discrete tokens for scGPT's transformer architecture.

### Step 4.3: Gene ID Mapping (Model-specific)

#### For TranscriptFormer & Geneformer:
```python
import mygene
mg = mygene.MyGeneInfo()

# Map gene symbols to Ensembl IDs
results = mg.querymany(
    gene_symbols,
    scopes='symbol',
    fields='ensembl.gene',
    species='human'
)

# Store mapping
adata.var['ensembl_id'] = ensembl_ids

# Remove version numbers: ENSG00000123456.1 → ENSG00000123456
adata.var['ensembl_id'] = adata.var['ensembl_id'].str.split('.').str[0]
```

---

## 5. Complete Pipeline Summary

```
┌─────────────────────────────────────────────────────────────┐
│                    RAW DATA (h5ad)                          │
│                    ~50,000+ cells                           │
└─────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────┐
│                    QUALITY CONTROL                          │
│  • filter_cells(min_genes=200, max_genes=5000)             │
│  • filter pct_counts_mt < 20%                              │
│  • filter_genes(min_cells=10, min_counts=3)                │
└─────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────┐
│                  INITIAL PROCESSING                         │
│  • Store raw counts in layers["counts"]                    │
│  • normalize_total(target_sum=1e4)                         │
│  • log1p()                                                 │
└─────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────┐
│              FEATURE ENGINEERING (CRUCIAL)                  │
│  • highly_variable_genes(n_top_genes=1200,                 │
│                          batch_key='Data_NO')  ← KEY!      │
│  • Value binning (51 bins) for scGPT                       │
│  • Gene ID mapping (Symbol → Ensembl)                      │
└─────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────┐
│                 EMBEDDING GENERATION                        │
│  • STATE: 2,058 dims                                       │
│  • UCE: 1,280 dims                                         │
│  • Geneformer: 1,152 dims                                  │
│  • TranscriptFormer: 2,048 dims                            │
│  • scGPT: 512 dims                                         │
│  • scGPT (BP): 512 dims                                    │
└─────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────┐
│                    47,808 cells                             │
│                    Ready for analysis                       │
└─────────────────────────────────────────────────────────────┘
```

---

## 6. Key Parameters Table

| Step | Parameter | Value | Impact |
|------|-----------|-------|--------|
| QC | min_genes | 200 | Removes empty droplets |
| QC | max_genes | 5,000 | Removes doublets |
| QC | pct_mt | <20% | Removes dying cells |
| Norm | target_sum | 1e4 | Standardizes depth |
| HVG | n_top_genes | 1,200 | Feature selection |
| HVG | **batch_key** | **Data_NO** | **Batch effect removal** |
| HVG | flavor | seurat_v3 | Selection algorithm |
| Binning | n_bins | 51 | scGPT tokenization |

---

## 7. Why Batch-Aware HVG Selection Works

### The Problem (scGPT Original)
```
Without batch_key:
  Gene X highly variable because:
    - Batch 1: high expression
    - Batch 2: low expression
    - This is TECHNICAL variation, not biological!

  Result: Embedding captures batch differences
         → UMAP shows "banana shapes" by batch
         → kBET = 0, iLISI = 1.43
```

### The Solution (scGPT Best Practices)
```
With batch_key='Data_NO':
  Gene X selected only if:
    - Variable WITHIN Batch 1 AND
    - Variable WITHIN Batch 2 AND
    - Variable WITHIN Batch 3 ...

  This selects BIOLOGICALLY variable genes
  that vary across disease severity, not batch.

  Result: Embedding captures biology
         → UMAP shows disease patterns
         → iLISI = 3.77 (164% improvement)
```

---

## 8. Code Reference

| Step | Script Location |
|------|-----------------|
| Data Loading | `scripts/16_embedding_generation_5models/code/run_scgpt_best_practices.py:399-401` |
| QC | `scripts/foundation_models/TranscriptFormer/prepare_data_for_transcriptformer.py:138-176` |
| Normalization | `scripts/16_embedding_generation_5models/code/run_scgpt_best_practices.py:113-118` |
| HVG Selection | `scripts/16_embedding_generation_5models/code/run_scgpt_best_practices.py:122-178` |
| Binning | `scripts/16_embedding_generation_5models/code/run_scgpt_best_practices.py:181-206` |
| Embedding Gen | `scripts/16_embedding_generation_5models/code/run_scgpt_best_practices.py:265-375` |

---

**Generated:** January 2026

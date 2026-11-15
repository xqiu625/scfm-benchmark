# Batch Effect Analysis: 6 Foundation Models Comparison

## Method

### Dataset
The dataset consists of 47,808 platelet cells, integrating COVID-19 and Sepsis PBMC scRNA-seq data aggregated from 11 independent studies (Wilk et al., Lee et al., Schulte-Schrepping et al., Arunachalam et al., Ren et al., Liu et al., Stephenson et al., Jiang et al., Bernardes et al., Combes et al., Qiu et al.).

### Models Evaluated
- STATE
- UCE
- Geneformer
- TranscriptFormer
- scGPT
- scGPT with Best Practices (scGPT-BP)

### Sampling
5,000 cells were randomly subsampled per model (with fixed random seed=42) from cells with valid severity labels, maintaining representation across all 11 batches.

### Evaluation Metrics

We evaluated model performance using six complementary metrics spanning two categories. The metrics listed for each category are distinct and complementary measures used to evaluate different aspects of that specific objective.

#### Bio-conservation Metrics (Severity-Based)

These are independent metrics that quantify different facets of how well biological information (disease severity) is preserved in the embeddings. They provide a multi-faceted view of biological signal preservation.

| Metric | What it measures | Interpretation |
|--------|------------------|----------------|
| **Bio_ASW** (Biological Average Silhouette Width) | Cluster separation - how well severity groups are separated in embedding space | Higher = better separation of severity classes |
| **ARI** (Adjusted Rand Index) | Clustering agreement - overlap between k-means clusters and true severity labels | Higher = better agreement with true labels |
| **NMI** (Normalized Mutual Information) | Information-theoretic measure of clustering agreement with severity labels | Higher = more information shared between clusters and labels |
| **cLISI** (cell-type Local Inverse Simpson's Index) | Local biological purity - whether cells are surrounded by same-severity neighbors | Lower = better (cells near same-severity neighbors) |

*Note: For radar visualization, cLISI is inverted so higher = better for all metrics.*

#### Batch Mixing Metrics (11 Data Sources)

These are independent metrics that assess how effectively batch effects have been removed and how well different batches (data sources) are integrated in the embedding space.

| Metric | What it measures | Interpretation |
|--------|------------------|----------------|
| **Batch_ASW_neg** (Negative Batch ASW) | Overall batch separation (negated) - lack of batch clustering | Higher = batches more intermixed (less separated) |
| **kBET** (k-nearest neighbor Batch Effect Test) | Local batch composition - whether neighborhoods match global batch proportions (k=50) | Higher = local neighborhoods match global batch distribution |
| **iLISI** (integration Local Inverse Simpson's Index) | Local batch diversity - effective number of data sources in local neighborhoods (perplexity=30) | Higher = more batches represented locally; max=11 for perfect mixing |

All metrics were normalized to [0,1] scale for radar visualization, where higher values indicate better performance.

---

## Results

### Quantitative Metrics

| Model | Bio_ASW | ARI | NMI | cLISI | Batch_ASW_neg | kBET | iLISI | n_batches |
|-------|---------|-----|-----|-------|---------------|------|-------|-----------|
| STATE | -0.010 | 0.018 | 0.048 | 3.56 | 0.026 | 0.177 | 4.02 | 11 |
| UCE | -0.013 | 0.008 | 0.065 | 3.15 | 0.014 | 0.018 | 2.35 | 11 |
| Geneformer | -0.013 | 0.006 | 0.009 | 3.68 | 0.027 | **0.371** | **4.31** | 11 |
| TranscriptFormer | -0.021 | 0.006 | 0.028 | 3.74 | **0.222** | 0.108 | 3.84 | 11 |
| scGPT | **-0.056** | **0.050** | **0.119** | 2.87 | -0.198 | 0.000 | 1.43 | 11 |
| scGPT_BP | -0.063 | -0.000 | 0.023 | 3.48 | 0.214 | 0.090 | 3.77 | 11 |

---

## Result Interpretation

### 1. scGPT shows severe batch effects but best bio-conservation
- scGPT achieved the highest ARI (0.050) and NMI (0.119), indicating its embeddings best preserve disease severity distinctions
- However, scGPT has the worst batch mixing:
  - kBET = 0.000 (complete rejection - local neighborhoods don't match global batch composition)
  - iLISI = 1.43 (cells predominantly surrounded by same-batch neighbors)
  - Batch_ASW_neg = -0.198 (negative value indicates batch separation)
- This suggests scGPT embeddings are dominated by technical batch variation rather than biological signal

### 2. Geneformer achieves best batch integration
- Highest kBET acceptance rate (0.371) and iLISI (4.31 out of 11 possible)
- This indicates Geneformer embeddings successfully mix cells across the 11 data sources
- However, bio-conservation metrics are modest (ARI=0.006, NMI=0.009)

### 3. scGPT Best Practices dramatically improves batch mixing
- Compared to standard scGPT:
  - iLISI improved from 1.43 → 3.77 (+164%)
  - kBET improved from 0.000 → 0.090
  - Batch_ASW_neg improved from -0.198 → 0.214
- The key preprocessing change: **batch-aware HVG selection** using `batch_key='Data_NO'`
- Trade-off: bio-conservation decreased (ARI: 0.050 → -0.000), suggesting the best practices preprocessing removes some biological signal along with batch effects

### 4. STATE and UCE show moderate performance
- **STATE**: Reasonable batch mixing (iLISI=4.02) with modest bio-conservation
- **UCE**: Poor batch mixing (iLISI=2.35, kBET=0.018) suggesting batch effects remain prominent

### 5. TranscriptFormer balances both objectives
- Good batch mixing (Batch_ASW_neg=0.222, iLISI=3.84)
- Moderate bio-conservation metrics

---

## Key Insight

There exists a fundamental **trade-off between batch correction and biological signal preservation**:
- Models with aggressive batch mixing (Geneformer) may lose disease-relevant variation
- Models preserving biological signal (scGPT) retain batch confounding

The scGPT Best Practices approach demonstrates that **preprocessing strategies can substantially improve batch integration**, highlighting the critical importance of data preparation in foundation model applications.

---

## Figures

### Generated Figures
1. `radar_bio_conservation_6models.png` - Bio-conservation metrics (severity-based)
2. `radar_batch_mixing_6models.png` - Batch mixing metrics (11 data sources)
3. `batch_effect_radar_6models.png` - Combined radar (all metrics)
4. `batch_effect_heatmap_6models.png` - Heatmap of all metrics
5. `batch_effect_metrics_6models.png` - Bar charts by metric category

### Code Location
```
scripts/12_umap_visualization/code/selected_code/
├── batch_effect_analysis.py
└── submit_batch_effect_analysis.sh
```

---

**Generated:** January 2026

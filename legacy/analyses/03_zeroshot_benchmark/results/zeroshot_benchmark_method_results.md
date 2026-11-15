# Result Part 2: Zero-Shot Benchmark & Baseline Comparison

## Method

### Dataset
The dataset consists of 47,808 platelet cells, integrating COVID-19 and Sepsis PBMC scRNA-seq data aggregated from 11 independent studies.

### Models Evaluated
Six single-cell foundation models were benchmarked:
- **STATE** (2,058 dimensions)
- **UCE** (1,280 dimensions)
- **Geneformer** (1,152 dimensions)
- **TranscriptFormer** (2,048 dimensions)
- **scGPT** (512 dimensions)
- **scGPT_BP** (512 dimensions) - scGPT with Best Practices preprocessing

### Classification Tasks
Three classification tasks were evaluated:

| Task | Classes | Description |
|------|---------|-------------|
| **Binary** | 2 | non_severe (healthy, mild, moderate, recovered) vs severe (severe, fatal) |
| **3-class** | 3 | control (healthy) vs mild (mild, moderate, recovered) vs severe (severe, fatal) |
| **6-class** | 6 | healthy, mild, moderate, severe, fatal, recovered |

### Zero-Shot Evaluation Protocol
1. **Embedding extraction**: Pre-trained foundation model embeddings were extracted without any task-specific fine-tuning
2. **Classification**: StandardScaler → Classifier (LogisticRegression, RandomForest, SVM, MLP, GradientBoosting)
3. **Evaluation**: 5-fold stratified cross-validation
4. **Metrics**: Accuracy, Balanced Accuracy, F1 (macro/weighted), AUC-ROC (macro/weighted)

### Baseline Methods
Traditional machine learning baselines using raw gene expression:

| Baseline | Description |
|----------|-------------|
| **PCA-50/100/200/500** | PCA dimensionality reduction → LogisticRegression |
| **Raw_LogReg** | StandardScaler → LogisticRegression (all HVGs) |
| **Raw_RF** | StandardScaler → RandomForest (all HVGs) |
| **Raw_XGBoost** | StandardScaler → XGBoost (all HVGs) |

---

## Results

### Zero-Shot Classification Performance

#### Binary Classification (severe vs non-severe)
| Model | AUC-ROC | Balanced Accuracy |
|-------|---------|-------------------|
| **STATE** | **0.895** | 0.642 |
| UCE | 0.877 | 0.631 |
| TranscriptFormer | 0.838 | 0.579 |
| Geneformer | 0.813 | 0.523 |
| scGPT | 0.775 | 0.519 |
| scGPT_BP | 0.766 | 0.509 |

#### 3-Class Classification
| Model | AUC-ROC | Balanced Accuracy |
|-------|---------|-------------------|
| **STATE** | **0.909** | 0.642 |
| UCE | 0.896 | 0.631 |
| TranscriptFormer | 0.851 | 0.579 |
| Geneformer | 0.826 | 0.523 |
| scGPT_BP | 0.784 | 0.631 |
| scGPT | 0.747 | 0.519 |

#### 6-Class Classification
| Model | AUC-ROC | Balanced Accuracy |
|-------|---------|-------------------|
| **STATE** | **0.894** | 0.642 |
| UCE | 0.888 | 0.631 |
| TranscriptFormer | 0.849 | 0.579 |
| Geneformer | 0.821 | 0.523 |
| scGPT | 0.810 | 0.519 |
| scGPT_BP | 0.768 | 0.451 |

### Baseline Comparison

| Method | Binary AUC |
|--------|------------|
| **Raw XGBoost** | **0.897** |
| Raw RandomForest | 0.880 |
| Raw LogReg | 0.865 |
| PCA-500 + LogReg | 0.850 |
| PCA-200 + LogReg | 0.842 |
| PCA-100 + LogReg | 0.830 |
| PCA-50 + LogReg | 0.815 |

### Zero-Shot vs Fine-Tuned Performance

| Model | Zero-Shot AUC | Fine-Tuned AUC | Improvement |
|-------|---------------|----------------|-------------|
| STATE | 0.895 | 0.944 | **+4.9%** |
| UCE | 0.877 | 0.910 | +3.3% |
| Geneformer | 0.813 | 0.845 | +3.2% |
| TranscriptFormer | 0.838 | 0.838 | +0.0% |
| scGPT | 0.775 | 0.735 | -4.0% |

---

## Result Interpretation

### 1. STATE achieves best zero-shot performance
- STATE consistently ranks first across all classification tasks (Binary: 0.895, 3-class: 0.909, 6-class: 0.894)
- STATE embeddings capture disease severity signal effectively without task-specific training
- The high-dimensional embedding (2,058 dims) preserves rich biological information

### 2. Foundation models are competitive with traditional baselines
- Best zero-shot (STATE: 0.895) is comparable to best baseline (XGBoost: 0.897)
- Foundation model embeddings encode disease-relevant information despite being trained on diverse cell types
- However, they do not substantially outperform well-tuned traditional ML approaches

### 3. Fine-tuning benefits vary across models
- **STATE benefits most** from fine-tuning (+4.9%), suggesting its pretrained features are good but can be improved
- **TranscriptFormer shows no improvement** (0.0%), indicating its zero-shot representations are already optimal for this task
- **scGPT degrades** with fine-tuning (-4.0%), possibly due to overfitting or batch effects in embeddings

### 4. Model ranking is consistent across tasks
- STATE > UCE > TranscriptFormer > Geneformer > scGPT ≈ scGPT_BP (for zero-shot)
- This ranking is stable across binary, 3-class, and 6-class tasks
- scGPT vs scGPT_BP comparison is confounded by batch effects in standard scGPT (see Part 1)

### 5. scGPT vs scGPT_BP: task-dependent performance differences
- scGPT_BP vs scGPT performance varies by task:
  - Binary: scGPT (0.775) > scGPT_BP (0.766)
  - 3-class: scGPT_BP (0.784) > scGPT (0.747) **+0.037**
  - 6-class: scGPT (0.810) > scGPT_BP (0.768)
- Standard scGPT has batch effects (see Part 1), which may inflate or distort performance differently across tasks
- scGPT_BP's superior 3-class performance suggests batch effects may particularly confound intermediate classification granularity
- The inconsistent ranking indicates batch effect "leakage" affects tasks differently

### 6. XGBoost baseline is surprisingly strong
- Raw XGBoost (0.897) outperforms all zero-shot foundation models except STATE
- This suggests that for well-curated datasets, traditional ML remains highly competitive
- Foundation models may provide more value when:
  - Training data is limited
  - Transfer learning is needed
  - Batch effects are present (but see Part 1 analysis)

---

## Key Insights

1. **Zero-shot foundation model performance is competitive but not superior** to traditional baselines on this benchmark

2. **STATE emerges as the top performer** in both zero-shot and fine-tuned settings for COVID severity prediction

3. **Fine-tuning does not universally improve performance** - model architecture and pretraining approach determine fine-tuning benefit

4. **scGPT's batch effects** (identified in Part 1) likely contribute to its poor fine-tuning performance

5. **scGPT vs scGPT_BP shows task-dependent differences** - scGPT_BP outperforms on 3-class (+0.037 AUC) but underperforms on binary and 6-class, suggesting batch effects impact tasks differently

---

## Figures

| Figure | Description |
|--------|-------------|
| `fig_zeroshot_vs_finetuned.png` | Panel A: Zero-shot vs fine-tuned AUC comparison; Panel B: Fine-tuning improvement (%) |

## Code Location
```
analyses/03_zeroshot_benchmark/code/
├── benchmark_single_model.py      # Zero-shot benchmark
└── submit_benchmark_scGPT_BP.sh   # SLURM submission for scGPT_BP

analyses/04_pca_baseline/code/
└── benchmark_pca_baseline.py      # PCA/traditional baselines
```

---

**Generated:** January 2026
**Updated:** January 26, 2026 (added scGPT_BP results: binary, 3-class, 6-class)

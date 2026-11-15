# 13 - SHAP Interpretability Analysis

**Accession:** 13
**Status:** ⏳ PENDING — scoped under TODO §8.1 (per-FM gene attribution)
**Priority:** P1 (planned for BIB submission)
**Owner:** to be implemented after Option C embeddings (§9.3) complete

---

## Purpose

Use SHAP (SHapley Additive exPlanations) to explain why STATE achieves the best performance:
- Identify which embedding dimensions are most influential for severity classification
- Understand feature interactions in the embedding space
- Provide clinically interpretable explanations for model predictions

---

## Planned Analysis

### 1. SHAP Value Computation
- Compute SHAP values for STATE embeddings on binary classification
- Use TreeExplainer (for XGBoost) or KernelExplainer (for LogisticRegression)

### 2. Feature Importance
- Rank embedding dimensions by mean |SHAP| value
- Identify top 20 most influential dimensions

### 3. Visualization
- SHAP summary plot (beeswarm)
- SHAP dependence plots for top features
- Force plots for individual predictions (healthy vs severe)

### 4. Cross-Model Comparison
- Compare SHAP patterns across STATE, UCE, Geneformer
- Identify model-specific vs shared important dimensions

---

## Expected Output

### Results
| File | Description |
|------|-------------|
| `results/shap_values_STATE.npy` | SHAP values matrix |
| `results/feature_importance_STATE.csv` | Ranked feature importance |
| `results/shap_summary_stats.csv` | Summary statistics |

### Figures
| File | Description |
|------|-------------|
| `figures/fig_13_shap_summary.png` | SHAP beeswarm plot |
| `figures/fig_13_shap_dependence.png` | Top feature dependence |
| `figures/fig_13_shap_force_plot.png` | Individual prediction explanations |

---

## Code Template

```python
import shap
from sklearn.linear_model import LogisticRegression

# Train model
clf = LogisticRegression(solver='lbfgs', max_iter=1000)
clf.fit(X_train, y_train)

# Compute SHAP values
explainer = shap.LinearExplainer(clf, X_train)
shap_values = explainer.shap_values(X_test)

# Summary plot
shap.summary_plot(shap_values, X_test, show=False)
plt.savefig('figures/fig_13_shap_summary.png', dpi=300, bbox_inches='tight')
```

---

## Key Questions to Answer

1. Which embedding dimensions drive healthy vs COVID classification?
2. Are the same dimensions important for mild vs severe distinction?
3. Do different models rely on different embedding features?
4. Can we map influential dimensions back to biological processes?

---

## Dependencies

- `shap>=0.42`
- `matplotlib`
- `numpy`
- Pre-computed embeddings from `01_embedding_integrity`

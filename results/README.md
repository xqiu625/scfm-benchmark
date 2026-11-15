# results/

Cross-stage, machine-readable outputs that feed downstream stages or the
manuscript. Per-stage outputs stay in `scripts/NN_*/results/`; this top-level
`results/` is for **aggregated / promoted** artifacts shared across stages.

| Folder | Contents |
|---|---|
| `features/` | shared feature/signature tables |
| `models/` | trained model artifacts (versioned — see architecture §10) |
| `predictions/` | prediction outputs (versioned) |
| `validation/` | cross-stage validation metrics |
| `submissions/` | final submission-formatted outputs |

**Versioning (architecture §10).** For important model/prediction outputs, write
into dated run folders rather than overwriting:

```
results/predictions/2026-06-06_ridge_v1/{predictions.csv,metrics.csv,model.pkl,manifest.json,README.md}
results/predictions/latest.txt   # points at the current run folder
```

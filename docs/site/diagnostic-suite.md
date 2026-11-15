# scFM Diagnostic Suite

> Historical research design, not the portable package's feature list. The old
> `diagnose` dispatcher only supports parts of this workflow and requires the
> HPCC research checkout. Use `demo` / `benchmark` for the standalone package.

A three-view methodology for characterising a single-cell foundation model
on a clinical task. The three views are deliberately complementary: each
isolates one dimension along which a "well-performing" FM can still fail in
deployment.

| View | Question it answers | Failure mode it surfaces |
|---|---|---|
| Adaptation curve | How much labelled data does the FM actually need? | FM gives no advantage over from-scratch baseline |
| GSEA per-dim | What biology do the learned axes encode? | FM dimensions encode confounders, not the task |
| SHAP + Jaccard | Which genes drive predictions? Do FMs agree? | FMs reach the same AUC via incompatible mechanisms |

Research implementations live under `scripts/`. Their inputs and execution
requirements differ; the portable package does not unify these three views.

---

## 1. Adaptation curve

**Code:** `scripts/29_adaptation_curve/code/adaptation_curve.py`

For a target task, train a classifier on (`n` labelled cells from the target
domain) + (all available cells from the source domain), for
`n ∈ {0, 10, 50, 100, 500, 1000, all}`. Repeat with 10 seeds at each `n`.
Compare per-FM curves to a from-scratch PCA-200 + LR baseline.

A useful FM should satisfy **both**:

1. **Higher Y-intercept** (`n=0`) than from-scratch — i.e., the pre-trained
   representation transfers without any target-domain labels.
2. **Steeper early slope** — i.e., target-domain labels are absorbed faster
   than they are by a from-scratch baseline.

In our cross-disease test (COVID → sepsis mortality), **no FM satisfies
both**. PCA-200+LR's asymptote (0.964) beats every FM at full data, and only
scPRINT-2 beats PCA's Y-intercept (0.920 vs 0.891) — but scPRINT-2 is
already flat, so the labelled-data budget is wasted. This is the
publishable null finding the field needs to wrestle with.

## 2. GSEA per-FM dimension

**Code:** `scripts/24_option_c_canonical_input/code/gsea_per_fm_dim.py`

Treat each FM embedding dimension as a continuous score per cell. For each
dim, rank genes by Pearson correlation with the score and run GSEA against
MSigDB Hallmark + KEGG + Reactome. Report the strongest enrichments per dim
across all dimensions for each FM.

The result is a **two-axis biology** for the platelet dataset:

| Axis | Models | Pathways |
|---|---|---|
| Cytoskeleton / adhesion | STATE, UCE, TF | Actin, integrin, focal adhesion |
| Translation / coronavirus | scPRINT-2, STACK, Geneformer | Ribosome, COVID-19 disease |

This is how you discover that two FMs achieving the same AUC are using
fundamentally different biology to get there.

## 3. SHAP + Jaccard overlap

**Code:** `scripts/24_option_c_canonical_input/code/per_fm_shap_attribution.py`
and `per_fm_shap_gene_map.py`

Train the patient-CV classifier as in `evaluate`, then run SHAP on the
trained model to get per-feature (per-dim) importance. Map dim-level
importance back to genes via the GSEA correlation table (view 2). Take the
top-`k` genes per FM and compute pairwise Jaccard overlap.

The Jaccard matrix is diagnostic: two FMs in the same biological axis
(e.g., STATE and UCE) overlap at 0.4; FMs across the two axes overlap at
≤0.02. **High AUC similarity + low Jaccard = the FMs disagree on
mechanism** — clinically dangerous because failure modes are not
correlated, so ensembling cannot save you.

---

## Composite verdict

A model passes the suite if it has:

- adaptation curve dominating from-scratch baseline,
- a coherent biology in view 2 (not noise / batch confounders),
- high Jaccard with at least one other model with a different architecture.

On the platelet COVID-mortality task in our paper, **no FM passes all three**.
STATE comes closest. The framework is task-agnostic — point it at your own
clinical h5ad and you get the same diagnostic.

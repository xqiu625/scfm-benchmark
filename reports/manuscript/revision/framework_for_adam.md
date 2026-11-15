# Framework: Stress-testing single-cell foundation models on underrepresented platelets — a patient-level mortality benchmark

For: Adam | From: Xinru | 2026-09-09 | Venue: *Briefings in Bioinformatics* (positioned as a benchmarking review of five published models)

**One-line story.** At realistic clinical event counts, no foundation model buys detectable advantage over a strong from-scratch baseline (PCA + logistic regression), in any disease, direction, or design cell we tested. The benchmark's contribution is a measured resolution limit (~0.10 AUC at 25 fatal events; detecting 0.03–0.05 needs ~10³ patients) and a protocol any new model can be dropped into. This reframes small apparent gains in clinical single-cell benchmarks as a question of statistical resolution rather than model ranking, and provides a quantitative basis for designing adequately powered evaluations of future foundation models.

**Design.**
- Cohorts: 240 COVID-19 patients (25 fatal, in-hospital), 91 sepsis patients (17 fatal, 28-day), 70 healthy (reference only); positive control on monocytes (148 patients, 19 fatal, in-distribution).
- Everything is patient-level: repeated 10×5-fold cross-validation grouped by patient, every prediction aggregated to the patient before scoring; patient-level OOF ROC-AUC as the single primary metric.
- Inference: one pre-specified primary test (best FM vs from-scratch baseline, winner-corrected permutation test); all else exploratory with exact DeLong SEs; Benjamini–Hochberg within each family.
- Controls: PCA-200+LR baseline; null (cell count + disease); cohort-only; QC-metrics-only; nested incremental-value tests over QC+count+cohort; monocyte control.

**Results.**
1. Primary test: not met — the best FM's advantage (+0.064, Geneformer) is indistinguishable from winner selection itself (permutation P = 0.181; cross-check P = 0.073; stable across all 25 leave-one-fatal-out refits).
2. Signal exists and is not an artifact: mortality is learnable by 5 of 6 methods, and representations add information beyond QC + cell count + cohort (pre-specified nested test, P < 0.0001) — the direct answer to "could you predict from cell counts".
3. Within sepsis: complete null at late sampling windows (17 events; resolves only AUC ≳ 0.71), consistent with published signals being documented only in early windows.
4. Cross-disease: uninformative forward (the sepsis task is not learnable even in-domain), rank-indeterminate reverse; adaptation unresolvable at 32-patient test set.
5. Monocytes (in-distribution): no FM beats the baseline either — the null is not an OOD artifact.
6. One confirmatory separation: the two leaders (STATE, Geneformer) separate from UCE (range-statistic P = 0.011) — a statement about UCE being unsuited (the only method that fails the count+disease null), not about a pretraining advantage. The corpus-composition explanation was falsified by the corpus measurement itself.

**Per-model profile (the multi-axis texture the nulls sit on top of).** The nulls are all FM-versus-baseline; the models are not indistinguishable from each other.
- **STATE (SE-600M)** — tied top within-COVID (0.778); strongest reverse transfer (0.721); best calibrated (STATE-RF: Brier 0.095 / ECE 0.079, the deployment-relevant option); top pooled LOBO (0.753, tied); largest increment over QC+count+cohort (+0.065).
- **Geneformer (V2-316M)** — top point estimate within-COVID (0.786); its +0.064 over the baseline is the observed max the primary test evaluates; second reverse transfer (0.712, 0.009 behind STATE; top under the severe-excluded variant, 0.775); pooled LOBO 0.748 (tied top).
- **scPRINT-2 (small-v2)** — third within-COVID (0.740), statistically inseparable from the baseline; second pooled LOBO (0.751, tied top); mid-table on every other axis.
- **scGPT (whole-human)** — indistinguishable from the baseline within-COVID (0.720); most cohort-consistent (LOBO σ = 0.123) but lowest pooled LOBO (0.638) — low variance around a poor mean; second on monocytes.
- **UCE (650M)** — last within-COVID (0.661) and the weakest method throughout: the only method that fails the count+disease null (p = 0.157), the incremental-value test (p = 0.072), and every transfer configuration; the only model that prefers monocytes over platelets; the most dropout-resistant (70% dimension zeroing costs it nothing); least cohort-consistent (σ = 0.152).

No single model leads on every axis; the deployment guidance therefore falls to the from-scratch baseline, with STATE-RF noted only as the best-calibrated option if an FM is used anyway.

**What this resolves from your review.**
- Cell-level vs patient-level: fully patient-level; one metric, defined once, used everywhere; no bare "AUC".
- Imbalance: patient-paired inference with explicit power statements; sepsis reported as a resolution bound, not "no signal".
- Cell counts: null 0.554, cohort-only 0.486, QC-only 0.376 (a suppression artifact); representations add value over all three (pre-specified test).
- OOD dismissal: the in-distribution control returns the same null; title now says "underrepresented" (measured platelet content ~0.12% of bounding corpora).
- Two statistical QC gates are now part of Methods (a selection-correction floor invariant and a harness unit test); every number traces to a script and a CSV.

**Ask.** Approve this framework (all-null + resolution ruler, benchmarking-review positioning) and we proceed to full text. The Introduction will carry a small review table of the three comparable patient-level scFM evaluations in the literature — none of which reports selection correction or event counts at the granularity of their comparisons.

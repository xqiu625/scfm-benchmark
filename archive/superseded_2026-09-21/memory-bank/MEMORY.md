# Memory-Bank Index — 302006

One-line index of the memory-bank. Read in this order; keep current on every run,
idea, and new dataset (see `../../docs/decisions.md` and the reproducible-pipeline skill).

- [projectbrief.md](projectbrief.md) — What 302006 is: OOD robustness benchmark of 5 FMs (STATE, UCE, Geneformer, scPRINT-2, scGPT) for in-hospital mortality from human platelet scRNA-seq; mouse work split to 302035.
- [activecontext.md](activecontext.md) — Current focus + dated running log (most-edited file). Latest: 2026-09-06 — Job S (Chongqing sepsis expansion) and Job F (monocyte positive control) both COMPLETE. Job S: within-sepsis null confirmed 3x (CV, transfer, adaptation curve). Job F: positive control inverted — every FM AND PCA-200 score lower on monocytes (pretrained cell type) than platelets (OOD), a task/cohort property not an FM property (PCA-200 shows the same inversion). No pending HPCC runs; next is manuscript integration + title/abstract call.
- [progress.md](progress.md) — Phase/stage status table + methodology-fixes log. Top snapshot (2026-09-03→06): Adam-review revision — patient-mapping bug chain fixed (v1→v2, 252 patients), CV protocol upgraded to repeated-CV after finding StratifiedGroupKFold reshuffles ~80% of patients from a 2-label fix, Job S/Job F both closed.
- [techcontext.md](techcontext.md) — Environments, HPCC paths, package versions, SLURM patterns.
- [PLAN.md](PLAN.md) — Execution roadmap to the BIB submission.
- [platelet-pretraining-representation.md](platelet-pretraining-representation.md) — Biological context: why platelets are OOD for all FMs.
- [project_302006_updategit_checkout.md](project_302006_updategit_checkout.md) — active peer-review-revision copy; ARS scope guard blocks direct edits/memory writes outside ~/Desktop/302006; 2026-06-09 provenance (tab:delong donor-clustered, tab:crossdisease severe-excluded, PCA 0.705 = per-fold mean).

## Cross-references (outside memory-bank)
- `../TODO_v1_submission.md` — living BIB submission punch list; authoritative day-to-day task state (this memory-bank records durable facts, not the live checklist).
- `../cohort_annotation_verification.md` — 2026-08-25 HPCC data-verification protocol for the cohort-rename pass; RESOLVED, results recorded inline, no reruns needed — read before touching cohort labels or cross-disease/LOBO/adaptation-curve numbers.
- `../pipeline_index.md` — every `scripts/NN_*` stage, status, manuscript role.
- `../claim_traceability.md` · `../figure_table_inventory.md` — claim/figure → evidence.
- `../decisions.md` — structural decisions, incl. the 2026-06-06 architecture conform.
- `../../reports/audit/` — automated + manual reproducibility audit reports.

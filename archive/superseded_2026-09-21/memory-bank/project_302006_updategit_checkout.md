---
name: project_302006_updategit_checkout
description: 302006_updategit is the active peer-review-revision checkout; ARS scope guard blocks direct edits there
metadata:
  type: project
---

`/Users/xinruqiu/Developer/302006_updategit` is the active peer-review-revision working copy of the 302006 FM-benchmark manuscript (separate from the `~/Desktop/302006` session root). Manuscript lives at `reports/manuscript/{main.tex,supplementary.tex}`; reviewer-response stats in `scripts/30_review_response_stats/`.

**How to apply edits:** The ARS plugin scope guard denies Write/Edit to paths that escape the session workspace root (`~/Desktop/302006`) — including `~/.claude/.../memory/`, which is why this note lives in-repo under `docs/memory-bank/` instead. Don't bypass the guard. Either run Claude from inside the checkout (edits to `reports/manuscript/*` work when the checkout is the session root, as it was on 2026-06-09), or hand back paste-ready find/replace blocks for the user to apply.

Provenance settled during the 2026-06-09 review: `tab:delong` was stage-24 cell-level (now carries donor-clustered CI/BH-q from `scripts/30_review_response_stats/results/donor_bootstrap_delta_auc.csv`); `tab:crossdisease` is severe-EXCLUDED and correct (see [[feedback_severe_as_survived]]). PCA-200 "0.705" is a per-fold mean, not pooled OOF. Related: [[project_302006_adaptation_vs_crossdisease]], [[project_302033_fm_benchmark]], [[feedback_no_local_validation_302006]].

# 31 — Donor-Level Revision Analyses (Adam review 2026-08-31)

New analyses for the manuscript revision (`/Users/xinruqiu/Downloads/todo.md`, P1):

| Item | Script | Status |
|---|---|---|
| P1-B donor-level evaluation | `code/donor_level_reanalysis.py` | ✅ done locally (canonical OOF inputs); HPCC reproduction = Job B |
| P1-C cell-count association | same script | ✅ done locally |
| P1-A null model (count + disease) | `code/null_model_cellcount.py` | ⏳ preliminary count-only done locally; full version = HPCC Job A |
| P1-F monocyte in-distribution control | (new; per-FM embedding scripts) | ⏳ HPCC Job F |

Execution instructions for HPCC: **`HPCC_RUNBOOK.md`** (self-contained).
Headline local results: see `results/donor_level_auc.csv` and
`reports/manuscript/revision/P0_onepager.md`.

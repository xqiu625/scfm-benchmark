# archive/superseded_2026-05-22/

Material moved out of `docs/` on 2026-05-22 once the canonical living
documentation shifted to `docs/memory-bank/` and the public mkdocs site at
`docs/site/`. Files here are retained for historical context only — do not
consult them for current project state.

| Item | Why archived | Where to look instead |
|---|---|---|
| `RESULTS_STATUS_v5.md` | Pre-Option-C severity-task snapshot (2026-04-30); superseded by the Option C / mortality-task framing | `docs/memory-bank/progress.md` |
| `RESULTS_STATUS_OPTION_C.md` | 2026-05-05 snapshot; still references the Geneformer LoRA chain as pending (the chain in fact completed 2026-05-06), and predates the cleaned 6-FM COVID+Sepsis subset numbers | `docs/memory-bank/activecontext.md` and `docs/memory-bank/progress.md` |
| `TODO_before_publication.md` | 2026-05-01 checklist framed as "Clinical Severity Prediction" — wrong task framing now that the manuscript is mortality binary; per-item status drifted from manuscript reality | `docs/TODO_v1_submission.md` |
| `rerun_v5/` | SLURM `.err` / `.log` dump from the Apr 26–27 v5 re-runs; pre-Option-C era; the numerical results were integrated into the manuscript directly | `06_LOGS/` on HPCC for live logs; `manuscript/main.tex` for integrated numbers |
| `statistical_significance_optimization_report.md` | 2026-01-27 report on the DeLong / bootstrap-CI fix (Critical fixes #2 + #3); the fix is now in code and reflected in `tab:delong` | `manuscript/main.tex` (Methods + `tab:delong`); `analyses/08_statistical_significance/code/` |

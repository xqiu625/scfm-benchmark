# Audit Report

> **Archived 2026-09-21. Historical snapshot.** This report covers the 2026-06 architecture and
> traceability audit. It does not include the 2026-09 Gate 1/Gate 2 audit. Use
> `docs/current-status.md` for current blockers and `docs/execution-plan.md` for
> the active remediation sequence.

Top-level summary of independent-review findings and their resolution. Detailed
per-dimension findings live in `reports/audit/`. The audit answers whether every
result, figure, and claim is traceable, reproducible, and internally consistent —
it does **not** generate new scientific conclusions.

Severity labels: `critical` · `major` · `minor` · `note`

## Audit-ready invariants (architecture §12)

A reviewer with no prior context should be able to answer:
1. What command generated each important result?
2. What data files did it use?
3. What code version + config produced it?
4. What figure/table uses that output?
5. What manuscript claim depends on it?
6. Can the value be recomputed from recorded files?
7. Any orphan figures/claims/outputs?

## Current status (2026-06-06, post-restructure)

| Dimension | Report | Status |
|---|---|---|
| **Automated audit** (reproducible-pipeline `audit.py`) | `reports/audit/automated_audit_2026-06-06_refresh.md` (latest) · `…_2026-06-06.md` (prior) | **refresh (post-fix): 0 critical · 132 major · 115 note** (was 1/132/115 → 0/213/30 originally) |
| Code review | `reports/audit/code_review.md` | not started |
| Data integrity | `reports/audit/data_integrity_check.md` | not started |
| Result reproduction | `reports/audit/result_reproduction.md` | not started |
| Figure/table check | `reports/audit/figure_table_check.md` | not started |
| Claim check | `reports/audit/claim_check.md` | not started |
| Unresolved issues | `reports/audit/unresolved_issues.md` | open |

### Automated audit interpretation (2026-06-06 refresh)

| Area | Count | Sev | Reading |
|---|---|---|---|
| portability | 130 | major | Hardcoded `/bigdata/.../302006` in code (**122 `scripts/`, 8 `legacy/`, 0 `configs/`** — the auditor only matches quoted paths, so configs are correctly clean). Real fix = migrate stage scripts to `from scfm_benchmark.config import paths` (added 2026-06-06). **Opportunistic** — not a mass rewrite (would risk live HPCC jobs; unvalidatable locally). Track progress by stages migrated, not the raw count. |
| stale-output | ~82 | note | **False positives** — the restructure touched every code file's mtime, so committed figures/results look "older than newest code." (The refreshed auditor now classes these as `note`, not `major`.) Dismiss; it's an mtime artifact of a git-based move. |
| manifest | ~~30~~ **0 actionable** | note | **A001 closed (2026-06-07): all 23 numbered stages (01–29) have `results/manifest.json`** via `scripts/99_audit/code/emit_manifests.py --all-stages` — 4 rich (06/17/24/29, with jobs/params/claims) + 19 minimal/stub (content hashes + git pin + README-derived description; `metadata_complete:false`). 7 stub manifests for HPCC-only/gitignored stages (04/07/19/22/26/27/28) carry `no_results_reason`. The remaining auditor manifest notes (`00_setup`, `00_shared`, `agent`, `cross_species`, `finetuning`, `foundation_models`, `tracking`) are **non-stage shared-code dirs** — not analysis stages, no manifest expected (auditor over-flags every `scripts/` subdir). |
| traceability | ~~1 critical + 2 major~~ **0 open** | resolved | **All 5 claims resolved (2026-06-06).** C001 removed, C002–C005 verified. The 2 `needs-evidence` majors are cleared: C003/C004/C005 reconciled to 3 dp against in-repo evidence. Verification corrected **3 wrong evidence-file pins** (C003 → `results_mortality_lora.json`; C004 → `covid_sepsis_only/`; C005 → `point_auc_mean` column). See `claim-traceability.md`. |
| evidence | 3 | note | `data/processed/covid_sepsis_mortality.h5ad` referenced but gitignored (lives on HPCC — expected); one `benchmark_*.csv` glob placeholder to make exact. |

**Traceability is now fully resolved** (C001 removed; C002–C005 verified). The remaining
automated findings are non-blocking: portability is opportunistic (live-HPCC risk),
stale-output is an mtime artifact of the git move, and the 30 manifest + 3 evidence items are
`note`. Re-run `audit.py` after manifests (`results/manifest.json`) are emitted to clear A001.

**Verification caveat worth keeping:** the C003/C004/C005 pass found that a *plausibly-named*
evidence file is not the same as a *reproducing* one. Three pins were wrong — each pointed at a
sibling cohort/run (overlap vs LoRA; full-pool vs COVID+Sepsis subset; bootstrap-mean vs
point-mean) that did not reproduce the cited number. The corrected pins do. Stages 17/24/29 keep
multiple cohort variants side-by-side; the manuscript headline tables use the COVID+Sepsis subset
(`covid_sepsis_only/`, main.tex L317).

## Open items carried from the restructure

| ID | Severity | Item |
|---|---|---|
| A001 | ~~major~~ **resolved** | **Closed (2026-06-07): all 23 numbered stages (01–29) manifested** via `emit_manifests.py --all-stages` (4 rich + 19 minimal/stub; content-hashed, git-pinned, provenance post-hoc). 7 stubs document HPCC-only/gitignored results via `no_results_reason`. Remaining auditor manifest notes are non-stage code dirs (`agent`, `finetuning`, `00_*`, …) — out of scope. Re-run `--all-stages` after any stage re-run to refresh. |
| A002 | ~~major~~ **resolved** | **C003/C004/C005 verified (2026-06-06):** cited numbers reconciled to 3 dp against in-repo evidence; 3 wrong evidence-file pins corrected (C003→`results_mortality_lora.json`; C004→`covid_sepsis_only/`; C005→`point_auc_mean`). `claim-traceability.md` rows flipped to `verified`. (`figure-table-inventory.md` figure/table rows may still need the same re-pin sweep — follow-up.) |
| A003 | major | In-repo paths now say `scripts/`; the HPCC working copy (`302006/`) must be renamed `analyses/`→`scripts/` to match, or absolute-path defaults break. Central resolver `scfm_benchmark.config.paths` added to ease the eventual migration. |
| A004 | note | `data/`, `results/` top-level layers are empty placeholders (real data on HPCC). |
| A005 | ~~critical~~ **resolved** | **C001 removed** (2026-06-06): verified against `main.tex` that the claim does not exist — STATE is never fine-tuned (L226); `0.948` is STATE's adaptation-curve asymptote (C005, Job 24018607); no `0.944` in manuscript. Was a stale scaffold seed from superseded 2026-01-26 DeepMLP CSVs. |
| A006 | ~~major~~ **resolved (false alarm)** | **C002 verified** (2026-06-06): manuscript `tab:crossdisease` traces to `_20260504_152836.csv` (STATE/UCE/scPRINT-2/Geneformer, severe-excluded `n_train=21286`) + `_20260430.csv` (scGPT); all 5 point AUCs match to 3 dp. `_20260522` is the severe-*included* rerun (`n_train=38673`), not the source — no missing 05-26 file. *Nice-to-have:* one consolidated 5-FM severe-excluded CSV. |
| A007 | **note (resolved)** | **iCloud working-tree corruption (2026-06-06):** the repo lives in `~/Desktop/` (iCloud Drive). Sync created **315 `" 2.*"` conflict copies** and **deleted 27 tracked files** mid-session (incl. all `cross_disease_fatal_vs_survived/`). Recovered: restored the 27 from git, deleted the 315 copies. **Recurs until the repo is relocated off iCloud.** |

Run `scripts/99_audit` (or the skill's `audit.py`) to refresh these checks once the C001/C002
reconciliation lands and manifests exist.

# docs/

Project documentation. Start with the current-state documents below. Several
traceability and public-site files are retained historical snapshots and are
labelled accordingly.

## Document lifecycle

| Label | Meaning |
|---|---|
| **CURRENT** | Binding project state or plan; update in place. |
| **APPEND-ONLY** | Durable decisions or run history; never rewrite old entries. |
| **REFERENCE** | Background material; verify its date before operational use. |
| **REFRESH REQUIRED** | Useful structure, but numbers/status predate the 2026-09 QC rebuild. |
| **ARCHIVED** | Historical only; never use as current state. |

## Current control documents

| Path | Purpose |
|---|---|
| **`current-status.md`** | **CURRENT.** Current phase, cohort scope, count ladder, and blockers. |
| **`execution-plan.md`** | **CURRENT.** Remaining work in dependency order (P0–P4). |
| **`data-contracts.md`** | **CURRENT.** Artifact schemas, invariants, and Gate 2 prior contract. |
| **`runbooks/gate1-gate2.md`** | **CURRENT.** Upfront-QC operational protocol. |
| **`runbooks/hpcc.md`** | **CURRENT.** Verified HPCC paths, environments, and job conventions. |
| **`decisions.md`** | **APPEND-ONLY.** Binding choices and withdrawals. |
| **`experiment-log.md`** | **APPEND-ONLY.** Major runs and conclusions. |

## Traceability documents

| Path | Lifecycle | Purpose |
|---|---|---|
| `pipeline-index.md` | Current structure; qualified statuses | Stage registry. Pre-rebuild `done` does not mean the rebuilt benchmark is complete. |
| `claim-traceability.md` | **REFRESH REQUIRED after P2/P3** | Current rows describe the pre-rebuild manuscript. |
| `figure-table-inventory.md` | **REFRESH REQUIRED after P2/P3** | Current rows describe pre-rebuild figures and tables. |
| `submission-checklist.md` | Logistics checklist | August history plus P4 logistics; scientific work follows `execution-plan.md`. |

## Reference material (rarely changes)

| Path | Purpose |
|---|---|
| `cohort-sources.md` | Source-dataset checklist. Verify inclusion against the current eight-cohort scope. |
| `platelet-representation-evidence.md` | Active evidence note with per-model evidence grades. |
| `models/` | Per-FM best-practice notes + paper summaries (STATE, UCE, scGPT, Geneformer, TranscriptFormer, plus Pertpy and scPEFT reviews). Reference material — not regularly updated. |
| `tooling-reference.md` | Reference map for reusable tooling under `scripts/`. |

## Public site

`site/` contains documentation-site source. The home page and CLI reference now
describe the portable synthetic demo and embedding-evaluation package. Historical
results/design pages remain explicitly labelled; their numbers predate the
2026-09 QC rebuild and are not current scientific evidence. Packaging details
and the public/software boundary are in `package-demo.md`.

## Reading order for new collaborators

1. `../README.md` — what the project is
2. `current-status.md` — where the work actually stands right now, and what is blocked
3. `execution-plan.md` — what remains, in dependency order
4. `decisions.md` — why the binding choices were made
5. `runbooks/gate1-gate2.md` — the QC protocol currently on the critical path
6. `runbooks/hpcc.md` — environments, paths, SLURM patterns
7. `platelet-representation-evidence.md` — biological context, with its evidence graded
8. `site/diagnostic-suite.md` — the methodology framing of the public-facing tool

> `../archive/superseded_2026-09-21/` holds the former `memory-bank/`. It is history, not status;
> several of its headline results were later overturned. Do not start there.

## What was archived

- `../archive/superseded_2026-05-01/` — 6 superseded `submit_rerun_*.sh` iterations, scPRINT v1 wrappers, the lumped `cross_disease_generalization.py`, the broken step2 mouse pipeline, the empty scGPT noise stub, and old SLURM log dumps. See that folder's `README.md` for per-file rationale.
- `../archive/superseded_2026-05-22/` — pre-Option-C `RESULTS_STATUS_*.md` snapshots, the severity-task `TODO_before_publication.md`, the v5 SLURM log dump, and the Jan 27 statistical-significance report. See that folder's `README.md` for per-file rationale.
- `../archive/superseded_2026-09-21/` — former memory-bank, the 2026-06 audit summary,
  the resolved cohort-annotation checklist, and the superseded command snippets.

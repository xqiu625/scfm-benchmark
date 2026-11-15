# Superseded — `docs/memory-bank/`, archived 2026-09-21

> **Do not use anything in this directory as current project status, plan, or instruction.**
> It is retained for history and audit only. Several headline results here were overturned by
> later work and are still written as authoritative.

| | |
|---|---|
| Original path | `docs/memory-bank/` |
| Archived | 2026-09-21 |
| Last date any file can be trusted as current | 2026-09-11 (`activecontext.md`) |
| Archived at commit | see `git log -- archive/superseded_2026-09-21/memory-bank` |

The files are kept together and unmodified. They were **not** split across archive locations,
because each one's reliability depends on reading it against the others from the same period.

---

## Why it was archived

The directory had stopped describing the project and started misdirecting readers. Concretely:

- `PLAN.md` was last updated **2026-05-27** and describes a roadmap predating the QC rebuild
  entirely. Following it routes work around the actual critical path.
- `progress.md` stopped at 2026-09-10, `activecontext.md` is headed 2026-09-09 while containing
  a 2026-09-11 addendum, and none of them contain the 2026-09-21 Gate 1 / Gate 2 rebuild.
- `MEMORY.md` stated "No pending HPCC runs" — false since the rebuild began, and exactly the
  kind of claim that leads a new session to think computation is finished.
- Multiple files still present overturned results as headline findings.
- Cell and patient counts appear at many different levels (source population, benchmark cells,
  platelet-only, patients, defined-outcome subset) with no stable labelling of which is which,
  so numbers cannot be compared across files without re-deriving them.
- `activecontext.md` opens with ~47 lines of debugging narrative about rounds 4–17 overturning
  each other — valuable as audit history, actively misleading as "current context."

---

## Where each file's still-valid content went

| Archived file | Live replacement |
|---|---|
| `MEMORY.md`, `activecontext.md`, `progress.md` | `docs/current-status.md` |
| `PLAN.md` (roadmap) | `docs/execution-plan.md` |
| `PLAN.md:248` (no-synthetic-fallback policy) | `docs/decisions.md` |
| `projectbrief.md` (project identity) | `README.md` + `docs/current-status.md` |
| `techcontext.md` (environments, paths, SLURM) | `docs/runbooks/hpcc.md` |
| `platelet-pretraining-representation.md` | `docs/platelet-representation-evidence.md` — **reworked, not copied** |
| `project_302006_updategit_checkout.md` | tooling note, not project documentation; not carried forward |

---

## Specific claims in here that are known wrong

Flagged so a reader who opens these files anyway is not misled:

- **"No pending HPCC runs"** (`MEMORY.md:7`) — false from the start of the QC rebuild.
- **Geneformer LoRA "+0.006"** (`PLAN.md:21`, `projectbrief.md:164`) — capacity-confounded.
  The within-backbone effect is **+0.113**; see the 2026-06-09 entry in `docs/decisions.md`.
  Several files still carry +0.006 as a corrected headline.
- **"Platelets are absent from all five models' pretraining data"** — the inference chain runs
  from one CELLxGENE CellGuide statement to five separate corpora and is weakest exactly where
  it is stated most strongly (Geneformer's Genecorpus is not CELLxGENE-derived at all).
- **"Standard droplet scRNA-seq does not capture platelets"** — contradicted by this project's
  own 10x cohorts, which contain tens of thousands of platelets. Platelets are captured and then
  removed by conventional count/gene QC floors. See `docs/platelet-representation-evidence.md`.
- **All benchmark AUCs and cell/patient counts** — predate the 2026-09 QC rebuild, which changes
  the benchmark input. Current counts, bound to artifacts, are in `docs/current-status.md`.
- **Title / framing described as "locked"** — `activecontext.md` itself records a subsequent
  retitling and several re-framings.

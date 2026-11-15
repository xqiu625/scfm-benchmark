# Gate 2 resumption readiness — v1 (2026-10-11), read-only status check

Produced in response to a review (`review:20261011_183017` §5, item 4) asking whether
Gate 2's scoring-only pilot — paused 2026-09-25 (`docs/PAUSE_HANDOFF_2026-09-25.md`) — can be
resumed, and exactly what is missing to do so. **This is a status check, not a resumption.** No
permit was created, no worktree deployed, no dry-run or scoring job run, nothing in
`scripts/32_rebuild_upfront_qc/` was touched. All findings below are independently verified this
round (commands and real-function calls, not re-asserted from the handoff) — see the per-item
evidence.

**Same-round correction (2026-10-11, in place — no new version number):** a follow-up review
(`review:20261011_190243` §3) found the original text of §2/§3/§7 overclaimed —
"single missing precondition / nothing else missing" implied the permit alone completes
resumption, `sacct` returning empty was described as "no job of any kind submitted since the
pause" (false; Job 1/Job 2 ran in this window under a different filter) and as ruling out a
dry-run ever having happened (dry-run never calls `sbatch`, so `sacct` cannot see it either way).
§2/§3/§4/§7/§9 below are corrected; the underlying facts (no permit, no output directory, hashes
match) are unchanged.

## 1. Repository/implementation identity vs. the frozen pause commit (`07b5662`)

| Check | Result |
|---|---|
| `git log --oneline 07b5662..HEAD -- scripts/32_rebuild_upfront_qc/` | **empty** — of the ~26 commits since `07b5662` (through `f4e3352`), none touch this directory; all are stage-33 batch-effect-audit work |
| `git diff 07b5662..HEAD --stat -- scripts/32_rebuild_upfront_qc/code/` | **empty** — byte-identical to the pause commit |
| `git status --short -- scripts/32_rebuild_upfront_qc/` | clean |
| `git rev-parse HEAD` | `f4e3352` (descendant of `07b5662`, not equal to it — expected, since 26 unrelated commits landed) |
| `git rev-parse origin/main` | `644bcf4` — an ancestor of local `HEAD`, itself a descendant of `07b5662`; confirms `07b5662` was pushed, but none of the subsequent batch-effect-audit commits have been |

**VERIFIED: no drift.** The reviewed implementation files are still exactly the bytes reviewed at
`07b5662`. Resuming does not require re-auditing the implementation itself — only creating the
permit the pause handoff describes as the next step.

## 2. Does a permit already exist?

`grep -ril permit scripts/ docs/ reports/` returns only pre-existing design/review documents
dated ≤2026-09-25 (`gate2_parameter_ruling.py` and `gate2_scoring_config_GSM5102900_HC1.py`,
which *define* the permit schema and binding logic but are not a filled instance; the pause
handoff itself; `FINAL_REVIEW_2026-09-25.md`; a `PERMIT_SATISFIABILITY_CHECK.txt`, also
found duplicated inside two frozen deployment snapshots `_deploy_9ea4c6a`/`_deploy_aea9a61`, same
pre-pause date, not a new instance). `git log --all --oneline --grep=permit -i` and
`--grep=gate2 -i` since `07b5662`: both empty on every branch. **VERIFIED-ABSENT: no tracked
permit artifact exists.** This is the next authorization artifact to create — see §7, which also
lists what remains incomplete *after* the permit (deployment, permit validation, dry-run evidence,
per-job authorization).

**Scoping note (added per `review:20261011_201340` item 5):** the grep/git-log search
above was run against `HEAD=f4e3352` and is a statement about that revision, not a standing
guarantee. Once a permit commit is created (as authorized and recorded in
`GATE2_PILOT_PREPARATION_PACKET_v1.md` §6), its own commit message will contain "Gate2"/"permit"
by design, so a bare "`--grep=permit` is empty on every branch" claim made *after* that commit
exists would be false; any future re-check must be read against the revision it is run at, not
treated as perpetually re-confirming this one.

## 3. Does a dry-run, scoring job, or output directory already exist?

`sacct -u xinruqiu -S 2026-09-25 -E now` filtered for Gate2/`GSM5102900`/scoring: **empty** —
**no matching Gate 2 scoring job found in the scheduler record.** This does not mean "no job of
any kind since the pause" (Job 1 `29804463` and Job 2 `29868753`, both unrelated stage-33 work,
ran in this window and would not match this filter) and it does not establish that a dry-run never
happened: the wrapper's `--dry-run` mode runs entirely as a plain shell/Python process and never
calls `sbatch` (`reports/manuscript/revision/gate2_capture_audit_2026-09-24/run_capture_scoring.sh`
lines 112–130), so `sacct` cannot see it either way. A separate `find` across the project root for
`*dry*run*`/`*gate2*dry*`/`*capture_scoring*` found no saved dry-run log or output from any round —
**no saved dry-run evidence exists**, which is the honest claim. The configured scoring output
root (`/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/gate2_capture_scoring/`) **does not
exist**. The only post-pause-dated Gate 2 artifact under `05_RESULTS/` is the already-documented
score-free preflight (`gate2_capture_preflight/9aab7bef_20260924_225552/`, job `29073092`,
predates the pause and is already recorded in the handoff) — not a scoring run, not a dry-run.
**VERIFIED-ABSENT** (scoring job, output directory, saved dry-run evidence); **UNVERIFIED** (an
un-saved/unlogged dry-run having been run and discarded — absence of evidence, not proof of
absence, per `DEVELOPMENT.md` rule 3).

## 4. Config/implementation hash reproducibility

The handoff's "Reviewed pilot identity" table lists `config SHA-256
0f02df7c7e2a81da1efa2f371d6890878ad8f747146986ae5bc5db7a94f1dd41` and `implementation SHA-256
30831273d4d3030ef5f9b123ecefbb44af2ae37aa9b7a95aab18ffeb8de80d84`. Deployed current `HEAD`
(`f4e3352`, same `scripts/32_rebuild_upfront_qc/` bytes as `07b5662` per §1) to a fresh HPCC
scratch checkout and called the real production functions from
`gate2_scoring_config_GSM5102900_HC1.py` directly — `implementation_sha256()` and
`config_sha256(_config_body())`, not a reimplementation — under `302006_env_clean`. Both
reproduced the handoff's values exactly: **`impl` match, `cfg` match.** **VERIFIED-MATCH.** The
reviewed hashes are still live and correct; a permit written now can cite them unchanged.
Independently re-confirmed the six pinned `DEPENDENCY_VERSIONS` this round too (Python 3.10.18,
scanpy 1.11.3, anndata 0.11.4, numpy 2.2.6, scipy 1.15.3, sklearn 1.7.1) by importing each package
directly under `302006_env_clean` on HPCC — all six match exactly. **Scope note**: this version
check lives inside the scoring entry point (`gate2_score_capture_GSM5102900_HC1.py:224-237`), not
the wrapper's `--dry-run` path, which returns before that entry point ever runs (see
`run_capture_scoring.sh` lines 127–130). A dry-run alone will not exercise this check; it has been
confirmed separately here and should be re-confirmed in whatever environment an eventual real
deployment uses, not inferred from a passing dry-run.

## 5. Quarantine status of earlier Gate 2 outputs

`docs/current-status.md`: "All eight cohorts' Gate 2 outputs are quarantined", under three
distinct reasons, verified against `05_RESULTS/rebuild_upfront_qc/gate2_scrublet/`'s actual
`.INVALID_*`-suffixed directories — unchanged since the pause (`docs/decisions.md`'s last Gate 2
entry is dated 2026-09-25; every entry after that is stage-33 batch-effect-audit work).
**VERIFIED-UNCHANGED.** No un-quarantined Gate 2 output exists to be accidentally reused.

## 6. Intentionally-excluded untracked paths

Both `reports/manuscript/revision/capture_provenance_v1_2026-09-23/` and
`reports/manuscript/revision/methodology_review_2026-09-24/parameter_ruling_tests.txt` still show
as untracked (`??`) in `git status --short`; `git log --oneline 07b5662..HEAD` for both paths is
empty — never swept into a commit. **VERIFIED-PRESENT, still untracked**, as intended.

## 7. What is actually missing — not a re-audit, but more than one remaining step

Per `docs/PAUSE_HANDOFF_2026-09-25.md` §"Exact next workflow" and
`docs/execution-plan.md` P0 item 6, **the next authorization artifact to create is a tracked,
scoped pilot permit** — the implementation itself needs no further audit (§1). Creating the permit
does not by itself complete resumption: deployment to a clean worktree, the permit passing the
wrapper's real validation, saved dry-run evidence, a automated review dry-run review, and separate owner
authorization for the actual scoring job all remain open after the permit exists (§9 below lists
these explicitly). The permit must be a new commit (after `07b5662`) that:

- references `07b5662` as the reviewed implementation commit;
- binds exactly one job, one capture (`GSM5102900_HC1`), all eight approved prior settings
  (`0.005, 0.01, 0.02, 0.04, 0.08, 0.12, 0.16, 0.20`), scoring-only, no outcome access, no cell
  deletion, no retry/no second job;
- cites the reviewed hashes confirmed live in §4 (`cfg 0f02df7c...1dd41`,
  `impl 30831273...80d84`) and the reviewed output root;
- states the resource ceiling: `short` partition, `02:00:00`, `48G`, 4 CPUs, maximum one job.

After the permit commit exists: deploy to a **new clean detached HPCC worktree** (not the current
long-lived checkout), run the wrapper's **dry-run only**, get a automated review dry-run review, and only
then ask the owner for the single scoring-pilot authorization — a dry-run review is not itself
that authorization, and this round's Job 2 authorization does not carry over to it (different
script, different workstream, separately gated per `DEVELOPMENT.md`'s "no SLURM submission without
explicit, per-job authorization").

## 8. What this round does not do

Does not create a permit. Does not deploy a worktree. Does not run a dry-run or any SLURM job.
Does not re-litigate the eleven ruled Scrublet parameters (`GATE2_PARAMETER_RULING.md`) or the
scoring-only design — both are frozen and unaffected by anything checked here. Does not touch
`canonical_donor_labels.csv`, QC, or DecontX. Does not push.

## 9. Owner decision this enables

Whether to authorize creating the pilot permit, deploying it to a clean HPCC worktree, and running
the wrapper's dry-run (a reversible sequence — permit is a tracked commit, dry-run makes no
scoring call and writes no output per §3) is the next concrete item. Even after that sequence
succeeds, **the single scoring job itself is a separate, further authorization**, to be requested
only after a automated review dry-run review. This is procedurally independent of: (a) the `CV0201`/`CV0326`
identity conflict, (b) Cambridge canonical-grouping promotion, and (c) how to use this round's
Job 2 diagnostic results going forward — none of these three block or are blocked by Gate 2 permit
creation. It is not scientifically independent of the project's overall critical path: Gate 2's
retained set and a subsequent DecontX refit on it are the upstream input the final
Platelet/Monocyte 5-FM benchmark depends on (`docs/execution-plan.md` P0), so this single-capture
pilot's eventual success would not itself mean Gate 2 is complete for all eight cohorts — only
that the scoring implementation behaves correctly on one.

## 10. Evidence paths

```
docs/PAUSE_HANDOFF_2026-09-25.md
docs/execution-plan.md (P0 section)
reports/manuscript/revision/methodology_review_2026-09-24/GATE2_PARAMETER_RULING.md
scripts/32_rebuild_upfront_qc/code/gate2_scoring_config_GSM5102900_HC1.py
/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/gate2_capture_preflight/9aab7bef_20260924_225552/
```

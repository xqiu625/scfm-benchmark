# Gate 2 single-capture pilot preparation packet — v1 (2026-10-11)

**UNSIGNED DRAFT. No permit exists. Nothing in this document is authorization.** Produced per
`review:20261011_190243` §4–5, which reviewed `GATE2_RESUMPTION_READINESS_v1.md` and
asked for this preparation packet to be written now, separately from creating the actual permit —
so that if/when the owner authorizes the permit-creation-and-dry-run sequence, no further lookup
or re-derivation of hashes/versions/parameters is needed; everything cited below was already
independently verified this session (see `GATE2_RESUMPTION_READINESS_v1.md`).

This packet does **not** authorize creating the permit, deploying a worktree, or running the
wrapper. See §6 for the exact scope that would still require the owner's explicit message.

**Same-round correction (2026-10-11, in place — no new version number):** a follow-up review
(`review:20261011_201340`, items 1–5) found five operational gaps in §1/§4/§5 below —
no concrete permit/log file paths, an inaccurate rationale for why a shallow clone would be
unsafe, deployment verification scoped to one directory instead of all seven
`IMPLEMENTATION_FILES`, no explicit pin on the dry-run's working directory/`SLURM_SUBMIT_DIR`,
and "re-cited"/"confirmed live this round" evidence labels that don't say which prior report or
when. All five are fixed in place below; no scope or hash changed.

## 1. Draft permit fields (not a live permit)

| Field | Draft value | Status |
|---|---|---|
| `permit_id` | `gate2-GSM5102900_HC1-pilot-<actual authorization UTC timestamp>` | placeholder — generated only after real authorization, not invented now |
| `authorised_by` | — | **must be filled from the real authorization message, verbatim, after it is received** |
| `authorised_utc` | — | **must be the real UTC time of that authorization, not this document's write time** |
| `capture` | `GSM5102900_HC1` | confirmed |
| `implementation_commit` | `07b5662362ed71b7f7856dcdf4b86d27bf6a6e7e` | confirmed (`git rev-parse 07b5662`, this round) |
| `config_sha256` | `0f02df7c7e2a81da1efa2f371d6890878ad8f747146986ae5bc5db7a94f1dd41` | confirmed live this round via the real `config_sha256()` function |
| `implementation_sha256` | `30831273d4d3030ef5f9b123ecefbb44af2ae37aa9b7a95aab18ffeb8de80d84` | confirmed live this round via the real `implementation_sha256()` function |
| `scope` | one job, one capture, eight approved priors, scoring only; no outcome access, no cell deletion, no retry, no second job | matches `docs/PAUSE_HANDOFF_2026-09-25.md` verbatim |
| `resource_ceiling` | `{"partition":"short","time":"2:00:00","mem":"48G","cpus_per_task":4,"max_jobs":1}` | matches the pause handoff and the wrapper's own `#SBATCH` header |
| `output_root` | `/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/gate2_capture_scoring` | matches `OUTPUT_ROOT` in `gate2_scoring_config_GSM5102900_HC1.py` |
| `input_sha256` | see §2 | confirmed live this round |

Eight frozen priors (unchanged, re-cited not re-derived): `0.005, 0.01, 0.02, 0.04, 0.08, 0.12,
0.16, 0.20`. Binds the full 7,666-cell deposit, never the 65-cell canonical subset.

**Concrete paths (added per `review:20261011_201340` item 1 — the draft above only had
placeholders):**

- Permit file (tracked, committed, lives in the repo per
  `_assert_permit_is_repository_artifact`): `scripts/32_rebuild_upfront_qc/permits/GSM5102900_HC1_pilot_<authorisation-UTC-timestamp>.json`.
- Deployment worktree: `/bigdata/godziklab/shared/Xinru/302006/worktrees/<permit-commit-full-sha>/`
  (full 40-char hash, matching the convention of the two existing September worktrees).
- Dry-run evidence: saved to
  `/bigdata/godziklab/shared/Xinru/302006/06_LOGS/gate2_pilot_readiness/<permit-commit-full-sha>_<run-UTC-timestamp>/`
  — **outside** the deployment worktree, so the saved log/output never becomes an untracked file
  inside that clone.
- Every one of these three paths must be confirmed **not to already exist** immediately before
  it is created or written to; none may overwrite a prior artifact.

## 2. Input identity (re-cited, not re-queried)

Directory: `/bigdata/godziklab/shared/Xinru/302006/01_DATA/raw/GSE167363_qiu/extracted`

| File | SHA-256 (confirmed this round, matches audited config) |
|---|---|
| `GSM5102900_HC1_barcodes.tsv.gz` | `e6dbf9f587d9a6cf3ecff31cbf4d53b09c6faf9837c82d5c1c0a26227f377a02` |
| `GSM5102900_HC1_features.tsv.gz` | `c7408e82af270ca6308a8aa7ea659795a15ab76f2e3a85bd5db4cd3d53244199` |
| `GSM5102900_HC1_matrix.mtx.gz` | `4ea259bfc4de02cc9799df563779fce7c9bd4c889532d3ff1178e2caf545b9f8` |

## 3. Dependency environment (re-cited, not re-queried)

Pinned (`DEPENDENCY_VERSIONS` in `gate2_scoring_config_GSM5102900_HC1.py`): Python `3.10.18`,
scanpy `1.11.3`, anndata `0.11.4`, numpy `2.2.6`, scipy `1.15.3`, sklearn `1.7.1`. These six values
were independently confirmed live in `302006_env_clean` on HPCC by direct package import, as
recorded in `GATE2_RESUMPTION_READINESS_v1.md` §4, verified by implementation assistant on 2026-10-11 — a prior
round's check, re-cited here, not re-run as part of writing this packet. This check runs inside
the scoring entry point, not the dry-run path (`GATE2_RESUMPTION_READINESS_v1.md` §4 scope note,
itself confirmed by reading `run_capture_scoring.sh` lines 127–130) — it must be re-confirmed live
again at actual deployment time, in that deployment's own shell, rather than assumed from this
citation.

## 4. Deployment plan — no `git push` required

The pause handoff's original wording ("push the permit commit, create a new clean detached
worktree") assumes push access. Push remains the user's separate, unexercised decision this
session (`origin/main` is still `644bcf4`, an ancestor of local `HEAD`). The plan below achieves
the same guarantee — a new, clean, tracked checkout whose `HEAD` is the permit commit, with the 7
guarded implementation files byte-identical to `07b5662` — without a push:

1. Create the permit commit locally, on top of current `HEAD`, touching only the new permit file
   (no implementation file edited — the wrapper's own guard at
   `run_capture_scoring.sh` lines 37–41 would refuse a dirty/untracked/unstaged guarded file).
2. `git bundle create /tmp/gate2_pilot_<permit-commit-sha>.bundle HEAD` (full history, not a
   shallow clone).

   **Correction (`review:20261011_201340` item 2):** the original rationale here —
   "a shallow clone would break `git rev-parse HEAD`" — was imprecise. `git rev-parse HEAD`
   resolves fine from a shallow clone. What the wrapper's real check
   (`gate2_scoring_config_GSM5102900_HC1.py::_assert_implementation_commit`) actually requires is:
   the reviewed commit `07b5662` must exist as an object in the deployed repo so
   `git merge-base --is-ancestor 07b5662 HEAD` can run, and it must in fact be an ancestor — a
   shallow clone truncated after `07b5662` would fail exactly that check. A full-history
   `git bundle create ... HEAD` avoids the risk without needing a shallow-clone experiment to
   prove it.
3. `scp` the bundle to HPCC; `git clone /path/to/bundle.bundle
   /bigdata/godziklab/shared/Xinru/302006/worktrees/<permit-commit-full-sha>/` — a genuine new
   clone with real git history, matching the convention already used for the two existing
   September worktrees (`worktrees/9aab7befba0bd2a59ea1224aff95d261237611c6/`,
   `worktrees/codex_gate2_final_20260925T204838Z/`).
4. **Explicit detached checkout and ancestry verification** (`review:20261011_201340`
   item 2) — do not assume the clone's default branch HEAD already equals the permit commit:
   - `git checkout --detach <permit-commit-full-sha>` inside the clone; confirm with
     `git rev-parse HEAD`.
   - `git merge-base --is-ancestor 07b5662362ed71b7f7856dcdf4b86d27bf6a6e7e HEAD` — must exit 0.
   - `git status --short` must be empty.
   - **All seven** `IMPLEMENTATION_FILES` (not just the `code/` directory — the list includes
     `reports/manuscript/revision/gate2_capture_audit_2026-09-24/run_capture_scoring.sh`, per
     `review:20261011_201340` item 3) must have an identical git blob id between
     `07b5662` and this clone's `HEAD`, checked per-file with
     `git rev-parse 07b5662:<path>` vs. `git rev-parse HEAD:<path>` — a directory-level
     `--stat` diff is not sufficient evidence that every one of the seven files was checked.
5. Run the wrapper's `--dry-run` from that clone only, never from the long-lived, possibly-dirty
   main checkout — confirmed this round to in fact be dirty and on an unrelated commit
   (`35131cf3`, not a descendant of `07b5662`/`64d9d71`) when checked directly on HPCC, which is
   exactly the risk this bundle-based plan avoids.

## 5. Dry-run command and acceptance evidence to save

```sh
DEPLOY=/bigdata/godziklab/shared/Xinru/302006/worktrees/<permit-commit-full-sha>
LOGDIR=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/gate2_pilot_readiness/<permit-commit-full-sha>_<run-UTC-timestamp>
mkdir -p "$LOGDIR"   # must not already exist
cd "$DEPLOY"
SLURM_SUBMIT_DIR="$DEPLOY" \
GATE2_PILOT_PERMIT="$DEPLOY/scripts/32_rebuild_upfront_qc/permits/<permit-filename>.json" \
  bash reports/manuscript/revision/gate2_capture_audit_2026-09-24/run_capture_scoring.sh --dry-run \
  >"$LOGDIR/stdout.log" 2>"$LOGDIR/stderr.log"
echo $? >"$LOGDIR/exit_code.txt"
```

**Fixed per `review:20261011_201340` item 4:**

- `SLURM_SUBMIT_DIR` is explicitly pinned to the deployment path. The wrapper reads this variable
  preferentially; leaving it unset risks inheriting a value from the login shell that points at
  the old, dirty main checkout instead of this clean deployment.
- `GATE2_PILOT_PERMIT` is the absolute path to the permit file **inside this deployment**, not a
  path in the local repo or any other checkout.
- stdout/stderr are redirected directly to files in `$LOGDIR` — **outside** the deployment
  clone — with the wrapper's own `$?` captured immediately after and saved as the evidence of
  exit code. No pipe through `tee` or any other process is used for this invocation, because a
  piped exit code reflects the last command in the pipe, not the wrapper itself.
- The frozen wrapper (`run_capture_scoring.sh`) is not edited to add this redirection — redirection
  happens at the shell invoking it, not inside the script.

Per the review's correction (not just exit code 0, which a *missing* permit also returns — see
`run_capture_scoring.sh` lines 112–121): the saved log must show the literal lines `permit  :
present, tracked and matching` **and** `DRY RUN COMPLETE: validations passed, including the permit
gate.` (not the "DRY RUN COMPLETE: ... scoring is blocked for want of a valid permit" line, which
is the no-permit path). Also save/record alongside the log, all in `$LOGDIR`: the six dependency
versions actually present in that invocation's environment (§3, re-confirmed in that shell, not
assumed); the `config`/`impl` sha256 lines the wrapper itself prints (lines 55, 88–90); the
deployment's own `HEAD` commit; the permit file's own sha256; confirmation `git status --short`
was empty in the clone both before and after; confirmation that `PRODUCTION_PARAMETER_RULING`
remains empty; confirmation the output directory
(`output_root/GSM5102900_HC1__cfg0f02df7c7e2a__impl30831273d4d3`) does not exist or is empty both
before and after (a dry-run must not create scoring output).

## 6. Authorization status

**Authorization for permit-creation + clean deployment + dry-run was received** in this session on
2026-10-11/12 (recorded UTC time: see the permit file's own `authorised_utc`, which is this
authorization's record time — the exact dispatch time of the user's message is not independently
known to implementation assistant and is not invented here, per `DEVELOPMENT.md` rule 3). The message, verbatim:

> 我批准按这份审查创建并提交限定 permit、完成全新 HPCC 部署及 dry-run。完成后交 automated review 审查；不批准
> SLURM 评分作业、push 或细胞删除。

This is materially equivalent in scope to the draft approval text in
`review:20261011_201340` ("给你的单次批准文案"): create and locally commit the scoped
`GSM5102900_HC1` permit, deploy to a clean HPCC detached worktree, run the real wrapper's
`--dry-run`, and hand the evidence to automated review. It explicitly does **not** authorize the SLURM
scoring job, `git push`, or cell deletion — matching this packet's own non-scope in §7.

Per `DEVELOPMENT.md` ("No SLURM submission without explicit authorization, per job" — the same per-job
discipline extends to creating a permit that would make the production entry point accept a real
scoring job) and `docs/PAUSE_HANDOFF_2026-09-25.md`'s scoped-permit resumption sequence:

- **The actual scoring job remains a further, separate authorization**, not granted by this
  message, to be requested only after a automated review dry-run review, per the pause handoff.
- Outcome access, expanding the capture beyond `GSM5102900_HC1`, and automatic retry remain
  out of scope, consistent with the permit's own `scope` field
  (`REQUIRED_SCOPE` in `gate2_scoring_config_GSM5102900_HC1.py`).

## 7. What this document does not do

This plan document itself does not create a permit file, deploy a worktree, or run any command.
Following the authorization in §6, those actions were carried out separately and are recorded
with full evidence in a implementation handoff, not in this plan. Neither this document
nor the actions it authorizes include: `git push`; any SLURM scoring job submission; outcome
access; cell deletion; expanding the capture beyond `GSM5102900_HC1`; automatic retry; or any
change to `canonical_donor_labels.csv`, QC, or DecontX.

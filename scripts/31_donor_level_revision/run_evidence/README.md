# Run evidence — stage 31

`2026-09-22_scheduler_check.txt` holds the raw stdout+stderr of the queries below. It exists
because the claim *no unauthorized stage 31 job was submitted* was first a statement about runtime
history, which git cannot corroborate, and then relayed as filtered prose.

**The first version of this file did not deserve the word it used.** It called itself "verbatim,
untruncated, commands exactly as run" while hand-abbreviating three paths to `...`, piping the log
listing through `tail -20`, and letting `sacct`'s column widths truncate `WorkDir` to `…c+` and a
step id to `29010326.exte+`. The recapture fixes the mechanism, not just the text:

- **within one execution**, each command is printed by the same string that is then evaluated; the
  `run()` that does it is echoed into the transcript and lives in `capture_scheduler_evidence.sh`,
  committed beside the output. That is a property of the display mechanism — it does **not**
  certify the file was never edited afterwards, for which git history is the record. (An earlier
  version wrote "the transcript cannot drift from what ran", which claims more than the mechanism
  delivers, and asserted it while the mechanism lived only in an ad-hoc heredoc.);
- the stream is **unfiltered**. An earlier version piped the ssh output through
  `grep -v` to drop SSH banner lines and still called the result "raw stdout+stderr"; the banners
  are now present, and their presence is itself part of the provenance;
- the script runs under `set -euo pipefail`. Without `pipefail` a failed `ssh` whose pipeline
  ended in a successful command would exit 0 and leave a truncated transcript looking complete, so
  **check the exit code** — `bash capture_scheduler_evidence.sh > out.txt; echo $?` must be 0;
- no `head`/`tail`; no abbreviated paths;
- `sacct --parsable2`, so no field is cut to a column width;
- the log cross-check is a real `find`, not a directory listing.

## Replaying it

```bash
bash scripts/31_donor_level_revision/run_evidence/capture_scheduler_evidence.sh > out.txt
```

The script hard-codes the cluster, account, window and paths, so a rerun differs from the
committed transcript in the **clocks**, the **repo state** (both checkouts), the **scheduler
state**, and the **filesystem state** — it lists the results, legacy, v1 and logs directories,
whose contents and mtimes change.

Check the exit code. A zero-length or truncated transcript with exit 0 was possible before
`pipefail` was set.

## Recorded with the capture

| | |
|---|---|
| cluster | `cluster.hpcc.ucr.edu`, **slurm 24.11.1** (printed in the file) |
| account | `xinruqiu` |
| capture host + clock | in the header, with timezone and offset |
| cluster clock at query | printed before the first command |
| capture-script local repo HEAD | in the header — the checkout the **script** ran from |
| HPCC checkout HEAD + status | the **first two remote commands**, because the machine being inspected is not the machine the script ran on |
| window | `sacct --starttime=2026-09-22T00:00:00` |

### What `--starttime` actually selected here

Do **not** infer the window from the flag name. An earlier version of this README asserted "a
`sacct` window that starts later than a job's submission will not show it", and the captured
output contradicts it: `28996424` was **submitted 2026-09-21T22:57:20** and still appears, because
it was *running into* the window and ended 2026-09-22T00:45:53. State the parameter, read the
output; the flag's exact selection semantics on this Slurm build are not established here.

## What the four checks found

| Check | Result |
|---|---|
| `sacct` allocations, canonical name `doublet_detect_v2` | exactly one: **29010326**, COMPLETED, ExitCode 0:0 |
| `squeue` now | nothing of this stage queued or running |
| `bib_revision_tier2/` directories | one new one, matching the authorized `--out-dir`; legacy `contamination_audit/` still 16 files |
| `find logs -name 'doublet_detect_*'` | **two** job ids — 29010326 **and 28355823** |

### The second job id, and why it matters

`28355823` is `JobName=doublet_detect` (not `_v2`), submitted **2026-09-13T22:45:26** — the
original run that produced the stale artifacts registered in `../STALE_ARTIFACTS.md`. It is
**outside the queried window and carries a different job name**, so the canonical-name test would
not have found it. The previous capture hid it behind `tail -20`.

That is a live demonstration that **job name is not a complete identity test for a stage 31 run**,
not a hypothetical caveat.

## What this establishes — and what it does not

> Within the recorded account and scheduler window, exactly one job used the canonical stage-31
> job name `doublet_detect_v2` (29010326), and the captured filesystem checks found only the
> authorized v1 output directory alongside the pre-existing legacy directory. These checks do not
> exclude a job submitted with `--job-name` overridden, a python invocation bypassing the wrapper,
> a job that failed before producing output or logs, a submission under a different account, or a
> submission outside the queried window.

The `find` result is the concrete reason that last clause is not boilerplate.


## The two checkouts are not the same, and were not

Splitting the HEAD field (above) immediately showed it. At capture time:

| | commit |
|---|---|
| local checkout the script ran from | `974a5c1` |
| HPCC checkout being inspected | `35131cf` — **225 commits behind** |

and `git -C $REPO_ROOT status --short` reports `scripts/31_donor_level_revision/` as **untracked
(`??`)** on HPCC. That directory has only ever been `rsync`-ed there; HPCC has never pulled this
branch.

So a statement of the form "N tests pass at commit X" was never literally true of the HPCC runs.
What is true is that the files were byte-identical to the local ones at the time of the run — and
that had to be checked, not assumed. Checking it found a real gap:

```
tests/test_stage31_doublet_semantics.py   local ed76188230732d23   HPCC fa9c2a129f68408f
scripts/31_donor_level_revision/STALE_ARTIFACTS.md                 HPCC MISSING
```

The test file differed by exactly one trailing newline — the byte removed in `3fe30f3` to satisfy
`git diff --check`, after which the file was not re-synced. Semantically identical, so the 188
figure held; but *"188 tests at HEAD"* did not, and the difference was found by measuring rather
than by the change being large enough to notice. `STALE_ARTIFACTS.md` had never been copied there
at all.

All six files were then synced and the suite re-run. **The pass count binds to the six files
listed in the transcript at the commit the transcript names, not to a later HEAD.** A transcript
is committed after it is produced, so its `local HEAD` is always the commit before its own; what
must be checked is that the commits in between do not touch those six files —
`git diff --name-only <transcript HEAD>..HEAD` answers that. **That claim now has its own evidence file**
rather than living in a report: `2026-09-22_test_run.txt`, produced by
`capture_test_run_evidence.sh`, both committed here.

**Standing rule from this:** before reporting a test count, verify byte identity of the files the
run used against the committed ones. `rsync` is not a version-control relationship, and a remote
`git rev-parse HEAD` identifies the remote repo's baseline, **not** the code a stage 31 test used.

### `capture_test_run_evidence.sh`

Order is the point, and every link exits rather than reporting:

```
HEAD blob == working tree          (exit 74)
 -> sync -> remote == local        (exit 71)
 -> run  -> remote == local AGAIN  (exit 71)
```

with the remote query's own transport status checked (exit 72) and its row count and path set
validated (exit 73) — a failed or truncated query that still emitted plausible text would
otherwise read as agreement.

Step 0 exists because the header's `local dirty` counts the **whole repo** and says nothing about
these six files. An earlier version reported `dirty: 4` and left the reader to assume the six were
not among them; they happened not to be, which is luck, not a mechanism. Each file is now checked
against `git show HEAD:<path>` and both hashes are printed.

It records the six paths, full SHA-256 from both sides at each stage with timestamps, the rsync
commands, the interpreter and library versions, and the pytest command — echoed before it runs,
not merely defined in the script.

**Two bugs in this script were found by injection, not by reading it:**

- `remote_hashes` merged stderr, so the SSH banner lines entered the data stream and the diff
  failed on them **every time** — the check never once confirmed agreement. The comparison now
  runs on the structured hash rows while the raw merged output is still printed.
- the check was `diff … && echo "…"`. Under `set -e` a failing command on the left of `&&` is a
  condition, not an error, so a mismatch printed a diff and the script **still exited 0**. It is
  now an explicit `if`/`exit 71`.

- the comparison ran inside `diff <(local) <(remote | …)`. Bash does not propagate a process
  substitution's exit status, so an `ssh` that failed **after** printing six correct hashes still
  produced `AGREE`. The remote query now writes to a file whose status is checked.

A decorative check inside the tool built to stop decorative checks. Three injections, each
confirming a different link bites:

The injections are **not** a claim in this file. `verify_evidence_guards.sh` runs all five and
writes `2026-09-22_guard_verification.txt`; it exits non-zero if any guard lets its fault through,
so "the guards bite" is a replayable result rather than something reported here.

| Injection | Exit |
|---|---|
| remote shell exits 55 after printing six correct hashes | **72** transport failure |
| hash query truncated to five rows | **73** row-count |
| remote `cd` fails — pytest must not run in the wrong directory | **1** remote fail-fast |
| a tracked file differs from its HEAD blob | **74** commit binding |
| remote content drifts after the run | **71** post-run mismatch |

The harness had two reporting bugs of its own, both turning *the injection did not apply* into
*the guard is decorative*: the old/new pair was one argument split on `||`, which occurs inside the
injected shell code, and the injector's status was unchecked; and the replacement took the first
occurrence, which for the pytest command is the **echoed display block**, not the executed one.
Old and new are separate arguments now, the occurrence index is explicit, and a failed injection
is fatal.

## Failure propagation in `capture_scheduler_evidence.sh`

The remote heredoc ran under `set -uo pipefail` with no `-e`, so a middle `sacct`/`git`/`ls`
failure was masked as long as the last command succeeded — SSH returned 0 and the transcript
looked complete. Both remote blocks in the test capture also carry `set -euo pipefail`. Without it a failing
`cd $REMOTE_ROOT` would leave pytest running in the **wrong directory**, and ssh returns only the
last command's status, so the earlier failure would vanish — fail-closed at the hash query was not
enough. Injection C is the proof.

`run()` now checks each status and exits 70 on the first failure, naming the
command in the transcript, and the transcript echoes the function with `declare -f run` so the
shell itself reports the definition. An earlier version hand-wrote a one-line paraphrase under the
word "verbatim", dropping the failed command's text and the abort line. That is deliberately not remote `set -e`, which would abort before the
handler could say which command died.

Verified by injecting `run "false"` mid-sequence: exit 70, the transcript ends with
`*** capture ABORTED -- this transcript is incomplete`, and none of the later blocks print.

(The old comment credited `pipefail` with catching a failed ssh whose pipeline ended in a
successful command. That pipeline was the `grep -v` filter, already removed — the comment
described a version of the script that no longer existed.)

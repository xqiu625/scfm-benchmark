# Development and research verification

Every rule below exists because it was violated in this repository, and the incident is named so
the rule is not mistaken for generic advice.

---

## 1. A factual claim about the code must carry the command that produced it

Any statement of the form "it does X", "X is never passed", "the tests cover X", "all N samples
have property P" is a **measurement**, not a recollection. State the command and its output, or
mark the statement `UNVERIFIED`.

> **Why.** Three claims in this repo were written with full confidence and were false.
> `SCRUBLET_PARAMS` was documented in the manifest as a binding the code did not have — it was
> never passed to `sc.pp.scrublet`. A comment said the manifest was "written and fsynced";
> `write_manifest` has no `fsync` call, and the wording had been carried over from
> `gate2_checkpoint.py`, which does fsync, without the implementation. `min(n_samples,
> n_features)=17` was read as "17 usable features", which a *minimum* cannot tell you.
> Each was one `grep` away.

Corollary: **never infer model or code correctness from output that looks right.** Scanpy *warns*
rather than raising when automatic thresholding fails, then returns finite scores and an all-False
prediction — output that passes every field-level check while resting on no threshold at all.

## 2. A new test is not evidence until a mutation kills it

After writing a test, break the behaviour it claims to cover and confirm the test fails. Report
the mutation, not the passing run. A mutation that "survives" must first be shown to have actually
applied.

> **Why.** Repeatedly, tests here passed while asserting nothing.
> A test asserted `drop(columns=["y"])` was present — making *parse-everything-then-delete*, the
> exact defect, the passing condition.
> A test sliced `src[src.index('_by_status'):src.index('per_sample_diagnostics.csv')]`, but that
> filename also appears in a comment 500 lines earlier, so the slice ran backwards to `''` and the
> test asserted nothing while passing.
> A test searched `_ast.unparse(fn)` for `highly_variable_genes`, which the function's own
> **docstring** contains, so deleting every call still passed.
> A guard test used `'os.fsync(' in src` at file level, so one unrelated `fsync` anywhere would
> license every durability claim in the module.
> Two mutations were reported as "surviving" when the string replacement had silently matched
> nothing because a comment sat between the statements it was joining.

Prefer **AST assertions bound to the specific function** over substring searches. Where a text
scan is all that is available, say in the test docstring that it is a lint-level tripwire, not a
semantic proof.

## 3. Absence of a guarantee is not a guarantee of the opposite

Do not convert "unverified" or "not guaranteed" into a definite direction — in either direction.

> **Why.** While correcting the false claim that the manifest was fsynced, the replacement said
> "nothing here survives power loss", which asserts certain loss. Without `fsync` the data may
> still survive. What is absent is the guarantee, not the possibility.

Scope claims to what the code does: `After os.replace() returns, a subsequent process exit does
not remove the replaced manifest` — not "no crash can lose it". And state the first-write case:
"leaves the previous manifest" is false when the output guard requires a new empty directory.

## 4. Report filenames carry measured timestamps

Local reports are named from `date +%Y%m%d_%H%M%S`, run at write time. Never
computed, never estimated.

> **Why.** Four of the first five were invented, all running 11–22 minutes *ahead* of the real
> write time. An external auditor reads ordering from these names and was sent to a stale report
> that sorted after a newer one.

An in-place correction keeps the existing report's filename, which is part of the audit trail.

## 5. A proxy proves only its own boundary; test the real producer-consumer interface

Before promoting a check to an end-to-end guarantee, name the exact boundary it exercises, and
validate the schema emitted by the real producer and consumed by the real code. A double can
establish control flow; it cannot establish a third-party library's field names, dtypes,
numerical behaviour, filesystem durability, scheduler state, or ordering relative to work done
by a caller — unless that real boundary is exercised.

Never generalise **one file to all files**, **existence to coverage**, **matching counts to
matching identities**, **a clean single-commit diff to a clean cumulative diff**, **a checksum
to all random state**, or **a guard at function entry to a guarantee about the caller**.

> **Why.** The Gate 2 test double put `_calls` and `_scores` into `uns`, and the runner was
> written to match the double — every unit test passed while real scanpy stores calls and scores
> in `obs`, so a genuine successful run would have been recorded as
> `calls_not_strict_boolean`. The same category error recurred all through this repository:
> one Liu file's first 50 features were generalised to all 55 ("no antibody capture" — every
> file is CITE-seq); a total feature-count difference was read as an RNA feature-space mismatch
> when it was entirely the antibody panel; row counts were treated as identity agreement; a
> `git diff --check` on one commit missed a defect the cumulative range reported; a
> `simulated_checksum` that matched across a grid was taken as proof of common random numbers
> while a per-rate `random_state` still varied; and an in-function `assert_may_score` was
> described as refusing "before any I/O" when the caller had already loaded the capture.

State exactly what the evidence proves. Label the stronger claim `UNVERIFIED` until the real
interface is exercised, and say which test would upgrade it.

**Corollary for mutation harnesses.** Restore from a baseline you have verified, not from
whatever was on disk when the harness started. A harness that timed out mid-run left a mutated
file on HPCC; the next harness copied that mutant as its "clean" baseline, and both its
one failing test and its one no-op injection were artifacts of the drift rather than findings.
Check the starting hash against the committed file before injecting anything.

## 6. A guard the caller can switch off, or a check whose result changes nothing, is not a guard

Two ways a safety property is present in the source and absent in the behaviour:

**The caller supplies the answer.** If an argument lets a caller declare the dangerous condition
absent, the guard now protects only callers who did not want to bypass it. Derive the condition
from a frozen definition the caller cannot pass in. The test must exercise the real path, not a
path where the condition has been defined away.

**The result is recorded and then ignored.** A field named `*_verified` or `*_agree` that no
branch reads is provenance, not validation — and it reads as validation to everyone downstream.
Either it can fail the record it describes, or it must not be described as a check.

> **Why.** `assert_parameters_adjudicated()` was first defined and never called; the fix made it
> callable with `unadjudicated=`, and the tests passed `{}` — the caller asserting that nothing
> was open, which is the original defect with a new spelling. In the same module,
> `identity_verified=False` and `parameters_echo_agree.agree=False` were both computed, written
> to the row, and then allowed to proceed to `examined` and into the Jaccard whose positional
> alignment the first of them exists to establish. Separately, three Scrublet parameters sat as
> module constants under a comment calling them "fixed", with no ruling anywhere — a constant
> with a confident comment is the same defect as a recorded default.

Corollary: state a guard's scope as what it refuses, not what it hopes. `REQUIRED - ruled` is
checkable; "the parameters are adjudicated" is a claim about a document that may not exist.

## 7. Close the execution loop before reporting a change as complete

This rule applies to every repository and every language. A change is not complete merely because
the intended function exists, a field is present, or a test helper can exercise it. Before handoff,
trace one real value through the complete path:

`authority/input -> validation -> execution -> decision branch -> persisted output -> reload`

At each arrow, verify that the same value or identity is carried forward. In particular:

- **Configuration must drive execution.** Do not record one parameter object while the executed
  call uses constants, defaults, or a separately constructed dictionary.
- **Validation must control behaviour.** A computed `verified`, `valid`, `agree`, or `status` field
  is only an observation unless failure changes control flow or output eligibility.
- **All exits must satisfy the contract.** Success, expected failure, early return, retry, and
  exception paths must produce the same schema, serialization guarantees, provenance, and atomicity.
- **Tests must enter through the production boundary.** Calling a validator directly proves the
  validator, not that production invokes it. Supplying a test-only bypass proves the bypassed path,
  not the guard.
- **Provenance must describe the execution that occurred.** Record resolved versions, parameters,
  identities, and source revisions from the running path, not from comments or intended defaults.
- **Round trips must be strict.** Test persistence with the real serializer and without permissive
  fallbacks such as stringifying unknown objects, accepting NaN as JSON, or silently dropping fields.

Before saying "fixed", report the exact production entry point exercised, the failure mutation that
was detected, and the artifact or state read back after persistence. If any link was not exercised,
name that link `UNVERIFIED` instead of upgrading a partial test into an end-to-end claim.

---

## Operational constraints (standing)

- **`cluster.hpcc.ucr.edu` round-robins between two login nodes** (`skylark`, `bluejay`).
  Measured 2026-09-23: four successive `ssh … hostname` calls returned skylark, bluejay, bluejay,
  bluejay. So `/tmp` is NOT shared, and neither is the process table — a file written by one call
  can be missing from the next, and `pgrep` can report nothing while the job runs on the other
  node. Put working files under `/bigdata/...`, and never conclude a remote job has finished from
  an absent `/tmp` file or an empty `pgrep`.

  This produced a false diagnosis the same day: the mutation harness's baseline copy was in
  `/tmp`, a status check landed on the other node and found no baseline and no process, the
  shared `/bigdata` runner was mid-mutation, and the conclusion drawn was "the SSH was cut and
  left a mutant behind". The harness was running normally; the `rsync` issued to "restore" it
  overwrote a file the live harness was using. Read the log the job actually writes.

- **All code runs on HPCC.** No local syntax checks, no local script runs. Sync with `rsync`, run
  under `/rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean`.
  Do not run the real 5 GB input on the login node — it gets OOM-killed, and the resulting exit
  code is a fact about the probe, not about the code.
- **No SLURM submission without explicit authorization**, per job. A code-review pass is not an
  execution authorization. Never modify or cancel an existing job.
- **Never delete, overwrite or move an existing checkpoint or result.** A rerun writes to a new
  versioned directory; the old one stays as the record of what the old code produced.
- **Scoped commits.** Stage only the files the task touched — never sweep in unrelated edits,
  build intermediates or `.DS_Store`.
- **Push is a separate decision from commit**, and belongs to the user.

## Outcome-data boundary (`--integration-only`)

Two inputs, two different guarantees; claiming one for both is false.

| Input | Guarantee |
|---|---|
| patient list CSV | **strict** — `y` excluded by `usecols` at parse time, never in memory |
| raw h5ad | **incidental load, acknowledged** — `read_h5ad` cannot select obs columns; only `REQUIRED_OBS_COLUMNS` survives the loader, the rest are dropped before use and recorded |

The manifest separates **authorization** (`outcome_analysis_requested`, known up front, the only
field that may be predicted) from **facts** (`known_outcome_values_loaded`,
`outcome_values_inspected`, `outcome_analyses_executed`, `container_obs_columns`), each starting
at "has not happened" and flipped by the code that does the thing, with an immediate atomic write.

Never use mortality/outcome information to design, select or tune a QC, doublet or DecontX rule.

## Open decisions — do not close them silently

**Owner's collaboration instruction (2026-10-15):** ask Cindy directly
when scientific adjudication is needed; a `pending_owner` note in a handoff is insufficient.
Bring evidence, a recommendation and consequences; investigate factual gaps before transferring
them to her. Continue independent authorized work. The supplied agent-loop diagram establishes
a collaboration pattern, not an authorization to expand the MVP technology stack or submit jobs.

- `docs/scrublet_dimension_parameter.md` is **OPEN**. There is no default `--n-prin-comps`; the
  script and wrapper both refuse to run without one. Do not pick a value by back-deriving it from
  the samples that failed or by maximising coverage.
- `scripts/31_donor_level_revision/STALE_ARTIFACTS.md` registers artifacts that must not enter the
  manuscript evidence chain. Check it before citing any stage 31 number.
- A mechanical doublet candidate is **not** a deletion authorization. No final-authorization
  artifact exists, deliberately — not even an all-false column.

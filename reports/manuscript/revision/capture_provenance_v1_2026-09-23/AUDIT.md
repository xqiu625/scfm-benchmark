# Bundle audit v2 — 2026-09-24

Answers to the four blocking findings, the provenance gap and the explanatory point, each with the
command result rather than a recollection. v1's five self-audit answers are superseded where the
numbers changed.

**Read-only. Not staged, not committed. No builder, DecontX output, Gate 2 use class, runtime
config or living document touched. No job submitted, no rate or parameter selected.**

---

## P0 — "no outcome column was read" was false

**Confirmed, and the claim is now a mechanism rather than a sentence.**

v1 read every input with `csv.DictReader`, which materialises the whole row. Four inputs carry
outcome columns:

| input | outcome columns present |
|---|---|
| `stephenson_patient_outcomes.csv` | `Outcome`, `outcome_unique_within_patient` |
| `canonical_donor_labels.csv` | `event`, `event_definition`, `disease` |
| `ren_gse158055_channel_inventory.csv` | `outcome` |
| `qiu_gse167363_channel_inventory.csv` | `disease_state` |

so the values were parsed and held in memory. Not indexing `r['Outcome']` is not not-reading it —
the same defect this repository already recorded for `pd.read_csv` followed by a column drop, which
makes this a repeat rather than a new mistake.

No outcome-free authoritative source exists for the Stephenson patient↔sample↔pool mapping: the two
candidates (`stephenson_patient_outcomes.csv`, `patient_sample_map.csv`) both carry outcome. So the
boundary is enforced at parse time instead:

- `read_columns()` reads via `pandas.read_csv(usecols=...)`, so other columns are excluded **at
  parse time**, which is the same mechanism this repository already audited for the patient-list
  CSV;
- it **refuses** any requested column matching `OUTCOME_PATTERN`, which is deliberately broader
  than the four known names so a new input with a differently spelled label column is refused
  rather than silently read;
- it records, per input, the columns read and the columns excluded, into
  **`outcome_column_boundary.csv`**. The boundary is now an artifact a reader can check.

Measured: `inputs read with usecols: 9; inputs where an outcome column was excluded at parse time: 4`.

Claim now made, and nothing stronger: **no outcome column is parsed.** Values are not read,
inspected, propagated, exported or analysed. The inputs themselves still contain outcome columns,
and that is stated rather than glossed.

## P0 — 15 coverage-only placeholder rows in the membership table

**Confirmed and removed.** They are coverage facts, not memberships, and they had inflated
`n_membership_rows`, `n_biological_samples` and the capture-status distribution. Cillo read as 83
rows / 82 samples when 12 of those rows recorded no membership at all.

They now live in `donor_coverage_reconciliation.csv` only, which is the cleaner of the two options
offered. An invariant refuses `membership_kind='unresolved'` in the long table, and a test refuses
any donor listed as `canonical_without_membership` appearing among the table's donors.

```
before: 2,537 rows   Chongqing 134  Cillo 83   Wilk 25
after : 2,522 rows   Chongqing 133  Cillo 71   Wilk 23
coverage-only donors, reconciliation table only: 15
```

## P1 — empty string counted as one analysis channel

**Confirmed.** Two independent fixes, because either alone would leave the trap open:

1. the placeholder rows that carried the empty key are gone;
2. every summary count now filters empties through `nonempty()`, so a future empty cannot form a
   category.

```
                 v1 reported    non-empty    v2 reports
Chongqing            134           133          133
Cillo                 45            44           44
Wilk                  24            23           23
```

v1's "the summary is mechanically rebuildable" was true and useless: the summary was rebuilt with
the same wrong logic, so agreement showed consistency, not correctness. The v2 test asserts the
count equals the **non-empty** cardinality, and additionally that it differs from the
with-empties cardinality wherever empties exist — so the assertion cannot be satisfied by
reproducing the defect.

## P1 — Liu `resolved` contradicted the bundle's own definition

**Confirmed. All 877 Liu rows are now `provisional`.**

`SCHEMA.md` defines `resolved` as *the deposit identifies the capture unit directly*, while
`FINDINGS.md` simultaneously stated that reading `B1_10xlane1` as one 10x chip well is this
repository's annotation. Both cannot hold. No independent, locatable authority was found that the
lane token is the physical capture, so the status moved rather than the definition.

The capture key is kept — it is the deposit's own GSM-title token and the unit demuxlet ran per —
with the reason recorded in `ambiguity_reason`.

Wilk is `provisional` for the adjacent reason: its `channel_id` is builder-derived. After this
change, **every `resolved` capture key comes from deposit metadata**, and the two keys that rest on
repository annotation or builder derivation are both `provisional`:

| cohort | capture key | source | status |
|---|---|---|---|
| Qiu / Chongqing / Ren / Cillo | GSM / sampleID | GEO deposit metadata | resolved |
| COMBAT | `obs/Channel_ID` | the deposited h5ad | resolved |
| Liu | lane token in the GSM title | deposit token, **repository gloss** as a well | provisional |
| Wilk | build-summary `channel_id` | **builder-derived** | provisional |
| Stephenson | — | none available | unresolved |

## Provenance gap — raw → extract was not replayable

**Closed.** `extract_hpcc_membership.py` is the read-only step from the raw HPCC sources to the two
extracts, and it writes `RAW_INPUT_MANIFEST.txt` listing every raw input actually read: sha256 for
each of the 55 Liu `.best` files, and for the 6.4 GB COMBAT h5ad the `<size>:<partial-hash>`
identity an earlier audit already recorded for that exact path
(`gate2_source_crosswalk.csv :: source_checksum`, cohort `combat`). The size is re-measured at run
time and must match, which makes the reuse checkable rather than assumed.

Re-running it reproduces both extracts: `COMBAT: 836148 cells, 70 channels, 124 donors, 979
membership pairs` and `Liu: 55 lanes, 50 SNG donors, 877 membership pairs`.

Outcome boundary of the extraction itself: `h5py` opens **named datasets**, so only
`obs/Channel_ID` and `obs/COMBAT_ID` are materialised from a deposit whose obs also holds clinical
fields. This is strict, not incidental — unlike `anndata.read_h5ad`, which would materialise all of
obs. The Liu `.best` files contain genotype calls and no outcome field.

## Explanatory — 40 `exact member-set match` vs 39 resolved Cillo captures

Not a contradiction. `GSM5465163` is an exact match of an **empty** member set
(`n_donors_in_lib=0`, `members` blank): it has an `h5ad_suffix` but names no donor, so it
contributes no membership row and therefore no capture. Stated in `FINDINGS.md` F3.

```
join_status='exact member-set match'  suffix=True  members=True   39
join_status='exact member-set match'  suffix=True  members=False    1   <- GSM5465163
join_status='PARTIAL ...'             suffix=False members=True     3
join_status='PARTIAL ...'             suffix=False members=False    3
join_status='no h5ad cells'           suffix=False members=True     2
```

---

# v3 round — 2026-09-24

Four further blocking items. Each answer is a command result; the mutation section below carries
the **real unified diffs**, not labels.

## P0 — the `usecols` test could not detect a dropped `usecols`

**Correct, and my v2 report overstated the kill.** `read_columns()` ends with
`return df[list(usecols)].to_dict('records')`, so removing `usecols=` from the parser still yields
only the requested keys. The v2 mutation I ran replaced **both** lines — the parser call *and* the
projection — so the kill never isolated the parser argument. The test checked returned keys, which
cannot see the difference.

Fixed by capturing what the parser actually received. `read_columns(path, usecols, reader=None)`
takes an injectable reader; `test_the_parser_itself_receives_usecols` records every call, separates
the header read (`nrows=0`) from the data read, and asserts the **data read** got
`usecols=wanted`.

Proof it now bites: the **minimal** mutation — only `usecols=usecols` deleted from the data read,
projection left in place — is KILLED. Its diff is one line (see below).

## P0 — the COMBAT identity only verified size

**Correct.** v2 compared `os.path.getsize(...)` against the size half of
`6409089483:d4435c68dd07c170` and never recomputed the digest, so a same-size file with replaced
content would have passed while the bundle claimed a verified identity.

`cheap_checksum()` is now a byte-for-byte port of `gate2_source_crosswalk.py:62`
(`sha256(str(size) + first 1 MiB + last 1 MiB)`, first 16 hex, prefixed by the size), and
`assert_raw_identity()` requires the **whole** string to match.

Three tests: the port is checked by `exec`-ing the original function out of
`gate2_source_crosswalk.py` and comparing outputs on a real file — an identity recorded by one
algorithm and verified by a different one is not a verification — plus same-size/head-changed and
same-size/tail-changed probes, each confirmed to produce a different identity and to be refused.

## P1 — publication was not transactional

**Correct.** v2 overwrote the COMBAT csv, then the Liu csv, then the manifest.

`publish()` now builds every output in full, writes each to a `.inprogress` sibling, fsyncs, re-reads
it byte-identical, and only then `os.replace`s all three. A failure leaves all three previous
versions and removes every temp. An empty payload is refused before anything is replaced.

`dest` is injectable and the fault-injection tests run in a temporary directory. That is not
convenience: driving `publish()` into a failure while its guard is mutated away would otherwise
leave a genuinely empty extract in the bundle — which is exactly what happened during the v3
mutation round (below).

## P1 — COMBAT was `intended` on post-QC observed co-occurrence

**Correct.** The pairs are donors observed among a channel's **retained** cells after the
depositors' QC. No pool-design artifact in this repository lists the intended donors per channel
(`combat_channel_inventory.csv` carries only an `n_donors` count), so `intended` claimed loading
design from co-occurrence.

Changed to **`observed_assigned`** — a new vocabulary entry, deliberately distinct from
`observed_singlet` (genotype demultiplexing) because this is the deposit's own per-cell donor
assignment. An invariant now requires **every** `observed_*` row to carry a stated limitation, so
neither observed kind can appear without "absence of a row is not absence from the capture".

## Claims narrowed

**Hashing bytes is not reading a value.** Checksumming opens each input with
`open(path, 'rb').read()`, so every byte of the file — an outcome field's bytes included — passes
through a hash. The claim is now exactly: *no outcome column is parsed into a field, inspected,
propagated, exported or analysed*, with the byte-level reality stated beside it in both the script
docstring and `SCHEMA.md`.

**Cited evidence is frozen.** The Stephenson rows rest on `T7_doublet_prior_evidence/memo.md`,
which the builder never parses. Its sha256 is now in `INPUT_CHECKSUMS.txt`, marked
`(cited, not parsed)`, and a test requires it.

---

## Two environments

| | env A — project | env B — fresh CI-like venv |
|---|---|---|
| python | 3.10.18 | 3.10.18 |
| pandas | 2.3.1 | 2.3.3 |
| h5py | 3.14.0 | 3.16.0 |
| numpy | 2.2.6 | 2.2.6 |
| pytest | 9.1.1 | 9.1.1 |
| command | `python -m pytest test_capture_provenance.py -q -p no:cacheprovider` | same, plus `-c /dev/null` |
| **exit code** | **0** | **0** |
| result | 26 passed | 26 passed |

Env B is `python -m venv` + `pip install pandas h5py numpy pytest` — nothing from the project
environment. Artifacts are byte-identical after both runs and no `.inprogress` residue remains.

Extraction re-run from raw reproduces both extracts byte-identically
(`combat f949f8ea9997`, `liu 5822848b0e75`).

---

## Negative tests and mutations

`test_capture_provenance.py`: **26 tests**, each written to fail if a specific earlier defect
returns. Baselines `build 748f4e74` / `extract d80ef296` verified before injection and after
restore. Kills key on the test process **exit code**.

| mutation | file | rc | verdict |
|---|---|---:|---|
| **MINIMAL**: only `usecols=usecols` removed from the data read; projection left in place | build | 1 | KILLED |
| COMBAT membership back to `intended` | build | 1 | KILLED |
| the observed-limitation text cleared on every observed row | build | 1 | KILLED |
| cited-but-unparsed memo dropped from the frozen checksums | build | 1 | KILLED |
| coverage-only rows appended to the membership table again | build | 3 | KILLED (invariant `SystemExit`) |
| summary counts stop excluding empty keys | build | 1 | KILLED |
| Liu back to deposit-`resolved` on a repository gloss | build | 1 | KILLED |
| outcome-column request guard removed | build | 1 | KILLED |
| raw identity check reverts to SIZE ONLY | extract | 1 | KILLED |
| `cheap_checksum` diverges from the repository algorithm (tail dropped) | extract | 1 | KILLED |
| publication becomes direct sequential writes | extract | 1 | KILLED |
| empty-payload refusal removed | extract | 1 | KILLED |
| temp files no longer cleaned up on failure | extract | 1 | KILLED |
| the observed-limitation **invariant** removed, alone | build | 0 | equivalent — see below |

The two diffs that carry the most weight:

```diff
--- build_capture_provenance.py (baseline)
+++ build_capture_provenance.py (mutated)
@@ -98,3 +98,3 @@
 })
-    df = reader(path, usecols=usecols, dtype=str, keep_default_na=False)
+    df = reader(path, dtype=str, keep_default_na=False)
 return df[list(usecols)].to_dict('records')
```

```diff
--- extract_hpcc_membership.py (baseline)
+++ extract_hpcc_membership.py (mutated)
@@ -193,4 +193,5 @@
                      raise SystemExit(f'{name}: temp file does not match its payload')
-        for name in outputs:                       # only now is anything visible
-            os.replace(temps[name], os.path.join(dest, name))
+        for name, payload in outputs.items():
+            with open(os.path.join(dest, name), 'wb') as fh:
+                fh.write(payload)
     finally:
```

The first is the isolation that was missing: one line, the projection untouched, and the kwargs
test fires. The second reverts the transactional publish and the fault-injection test fires.

### The one non-kill, and why it is not a gap

Removing the observed-limitation **invariant alone** changes nothing: no row lacks a reason, so
there is nothing for it to catch. Reported as a survivor it would imply a test gap that is not
there. Its defect-producing form — clearing the reason on every observed row — is KILLED by two
tests. A surviving mutation has to be shown to reproduce the defect first.

### A process error this round, mine

The v3 run reported `rc=2 :: 1 error` for the temp-cleanup mutation. That was **not** that
mutation's own kill. An earlier mutation in the same run removed the empty-payload refusal, so the
fault-injection test published an empty payload **into the real bundle**, leaving
`combat_donor_by_channel.csv` empty; the harness restored only the two `.py` files, so the next
iteration failed at collection with `EmptyDataError`. A contaminated baseline again — the same class
CLAUDE.md rule 5's corollary names.

Fixed at the source rather than in the harness: `publish(dest=...)` is injectable and every
fault-injection test runs in a `TemporaryDirectory`, with an AST test asserting that no `publish()`
call in the test file omits `dest=`. Re-run after the fix: every verdict above, and the real
extracts byte-identical afterwards.

## Still not asserted

- whether a `provisional` capture is in fact one physical capture;
- whether any existing DecontX or doublet result should be recomputed;
- whether Stephenson can be made runnable at any granularity;
- how the `docs/data-contracts.md` contradiction should be resolved;
- whether Liu's missing 56th lane was never deposited, never run, or merged.

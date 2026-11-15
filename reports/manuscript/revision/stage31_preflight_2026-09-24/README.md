# Stage 31 input preflight — 2026-09-24

**Read-only. No doublet score, call or rate was produced. No `n_prin_comps` was selected. No job
was submitted. No outcome value was read.**

Answers the question the 2026-09-24 methodology ruling put ahead of every parameter: *can the
redesigned workflow's inputs be identified at all?* It does not score anything, and it does not
approve the redesign.

Produced by `scripts/31_donor_level_revision/code/stage31_input_preflight.py`, run from a detached
worktree at commit `bffadf8827c88266d27164288ed42adb1c4b5179`.

## What was measured

| | |
|---|---|
| COVID patients | 240 |
| **Pre-subsample `Sample_ID`s** | **302** |
| **`Data_NO` groups** | **9** |
| Pre-subsample cells | 36,076 |
| Cells in the 60-cell analytic subset | 10,737 |
| Container, all diseases | 413 `Sample_ID`s, 47,977 cells × 13,066 genes |

302 / 9 match the project record exactly. The container is `X` only — dense, `int64`, integral on
probe, no layers.

## The four findings

### 1. The 302-vs-301 gap is the subsample dropping a whole `Sample_ID`

`004_TA` has **16 cells pre-subsample and 0 in the 60-cell subset**. The cap is applied **per
patient**, so a patient with several `Sample_ID`s can lose one entirely — patient `004` did.

That is the documented 302-vs-301 discrepancy, now measured rather than described. It is also a
concrete instance of why the ruling rejected the subsample-first design: a unit can vanish before
it is ever scored, and nothing downstream records that it existed.

### 2. There is no single barcode/suffix rule — 29 shapes across 9 groups

A first pass assumed `<barcode>-<gem>_<suffix>` and left **76% of cell ids unparsed**. Measured
instead, the canonical namespace has **29 distinct id shapes, up to 14 inside a single `Data_NO`
group** (`cell_id_namespace_by_data_no.csv`). The shapes encode different things per group:

```
<BC>-<n>_<n>                       Seurat merge suffix
Platelet_<BC>-<n>                  lineage token prefix
Megakaryocytes_Ki<n>_Ki<n>TA<n>_<BC>-<n>   nested sample/visit prefixes
BGCV<n>_<BC>-<n>                   sequencing pool prefix
<BC>-MH<n>  /  <BC>-newcastle<n>   site or donor suffix
GSM<n>_<n>.<n>_<BC>                accession prefix, barcode last
```

**C2 cannot be satisfied by one rule.** Each `Data_NO` group needs its own documented
prefix/suffix rule before a raw-barcode namespace can be fixed.

### 3. Container-internal identity is clean; raw-input identity is NOT established

The one element common to every id is a nucleotide run. Extracting it:

- **0** cells with no extractable nucleotide token;
- **0** `(Sample_ID, token)` duplicates **within this container**;
- **0** subset cells ambiguous, **0** unmatched.

**What that does and does not show.** The token is pulled out of the canonical cell-id *string* in
this same h5ad. No raw matrix or barcode file was opened, and the token was never joined against a
10x barcode universe. So it is a `nucleotide_token_candidate`, **not a verified raw barcode**, and
`raw_barcode_identity_verified` is `false` on every row.

What is shown: each of the 10,737 subset cells links deterministically, by its full canonical
`cell_id`, back to one pre-subsample row **of the same container**. That is a container-internal
identity candidate. It is **not** a link to a raw source, and the raw-input identity mapping the
redesign needs has not been performed.

Nor is it a mapping to a physical capture. `physical_capture_key` is `unknown` for **302 of 302**
samples; nothing here establishes that a `Sample_ID` is one capture, and C2 forbids assuming it.

### 4. Two C1/C3 conditions remain unresolved

- **`raw_counts_un_normalised` = unknown.** `int64` and integral-on-probe are *necessary*, not
  sufficient: an integral matrix can still be the output of a transformation. Provenance to an
  un-normalised RNA source is not established here.
- **`parent_states_retained` = unknown.** The container is named for platelets and most ids carry
  no lineage token, but 6,030 carry a `Megakaryocytes` prefix and a few carry `T_cell` / `B_cell` /
  `Other_cell`. Those tokens are inherited from upstream object names and are **not** evidence of
  current lineage composition. The question is open in both directions.

## Files

| file | contents |
|---|---|
| `source_inventory.csv` | 302 rows, one per pre-subsample `Sample_ID` |
| `capture_barcode_crosswalk.csv` | 36,076 rows, one per pre-subsample cell; the token column is a candidate, `raw_barcode_identity_verified=false` |
| `cell_id_namespace_by_data_no.csv` | the 29 observed id shapes, per `Data_NO` |
| `subset_mapping_summary.csv` | 301 rows — the 302nd contributes no subset cell |
| `PRECHECK_MANIFEST.json` | inputs, identities, outcome boundary, unknowns |

## Outcome boundary

`h5py` opens **named datasets**, so only `Sample_ID`, `Data_NO`, `Cell_ID`, `Cell_ID2`,
`orig.ident`, `nCount_RNA`, `nFeature_RNA` and the index were materialised. `Death`,
`covid_severity`, `mild_vs_severe`, `severity_score`, `three_group`, `progression_group`, `label`
were **never opened** — the script refuses to run if any of them appears in the opened set. The two
CSVs are read with `usecols`.

`disease` **is** inspected, from the patient list, because the COVID population cannot be defined
without it. It is recorded in the manifest as a **population selector**, not an outcome, and no
value of it reaches any output — only the resulting id set does.

No output file carries mortality, Death, survival, severity or any other clinical value.

## What this does not do

It does not select `n_prin_comps`, does not score, does not approve the redesign, does not
establish a capture mapping, and does not establish raw-barcode identity. A new dimension audit belongs on the complete inputs **after** the
capture mapping and raw-count provenance above are resolved — not on this container.

# Gate 2 doublet rate — sensitivity-design review packet

**For independent methodological review. Nothing in this package selects a rate, authorises a run,
or authorises deleting a cell.**

| | |
|---|---|
| PDF | `sensitivity_design_review_packet.pdf` — 12 pages |
| PDF SHA-256 | `1ea5507daf5775fbdf8da1c2482ec0a52e511309b3a663968680ab3602343c5a` |
| Built | 2026-09-24 (typography rebuild; see *What changed* below) |
| Evidence frozen at | **2026-09-22**, commit `45aeb9f7c4ce0c3d64f7c5939020c8b954ed2c4e` |
| Attachments | 8 CSVs in `attachments/`, listed with SHA-256 in `SHA256SUMS.txt` |
| Build record | `BUILD_PROVENANCE.md` (both builds) |

## Read this before reading the packet

### External human review has NOT happened

No external methodological reviewer has seen this material. **You are the first.**

The project record contains a methodological summary dated **2026-09-23** in
`docs/decisions.md`. Two things about it must be stated plainly:

1. it is an **internal AI (Codex) opinion**, not an external human review;
2. its **original file is not available**. It was referenced as
   `codex_output/20260923_111811.md`; a full search of the project directory on 2026-09-24
   did not find it. It has not been reconstructed, and nothing in this package treats it as
   verified.

That summary proposes a **direct residual expected-rate grid**. Treat it as an unreviewed internal
proposal. You are asked to verify or reject it independently, not to build on it.

### The evidence is frozen at 2026-09-22; the project has moved since

Every number here reflects the evidence as of 2026-09-22. Decisions recorded on 2026-09-23 and
later are **not** incorporated. In particular the Stage 31 review task has changed — see the
Stage 31 addendum shipped separately; do not treat the older Stage 31 bundle's request to choose a
fixed `n_prin_comps` as the current question.

### What the packet is, and what it is not

It is an **evidence-gap matrix** across eight cohorts plus the reasoning behind it: what is known
per cohort about chemistry, loading, recovered counts, upstream doublet handling and donor
proportions, each carried with an explicit evidence grade.

It is **not** a rate proposal. No cohort has an obtainable recovered count; a quantity graded
`assumed` or `unknown` is never presented as an estimate of a real prevalence.

## What changed in this build (2026-09-24)

**Typography only. No value, rate, cohort classification or matrix cell was changed** —
`cohort_evidence_gap_matrix.csv` is byte-identical to the version committed in `45aeb9f`, and the
builder's 12 enumerated invariant checks pass unchanged.

The 2026-09-22 build ran text off the **physical page edge** on three pages. Measured with
`pdftotext -bbox` against the 612pt page width:

| page | worst xMax | content |
|---|---|---|
| 2 | 613.6pt | `scripts/32_rebuild_upfront_qc/code/build_…` |
| 3 | 615.1pt | `not_eligible_without_a…` |
| 4 | 615.1pt | `not_eligible_without_a…` |

Two causes, both unbreakable words that no amount of line-stretching can fix:

- pandoc escapes `_` as `\_`, at which TeX will not break, so a 55-character snake_case value is
  one word;
- the matrix tables used a short separator row, so pandoc emitted non-wrapping columns.

Fixed at the source: a zero-width break after each underscore, wrapping table columns, and a
smaller face for the matrix. `hyphenat[htt]` was tried and **rejected** — it fixed the width but
inserted a real hyphen inside file paths (`scripts/32_rebuild_upfront_-qc`), which a reader could
copy as part of the path.

After the rebuild the worst overhang anywhere is **552.7pt**, i.e. inside the 612pt page and inside
the printable area on every page. All 12 pages were re-rendered and checked for page-boundary
violations; pages 2, 3 and 4 were additionally inspected visually. The packet is one page shorter
than the 2026-09-22 build because the matrix is set more compactly.

The superseded build (13 pages, 90,556 bytes, SHA-256 `5c08b685…`) remains recorded in
`BUILD_PROVENANCE.md`. **Use this build; do not circulate the 2026-09-22 PDF.**

## Attachments

The PDF is self-contained — the full matrix is rendered inline and you do not need to open a CSV to
read it. The CSVs are attached for machine use and for checking the narrative against the rows.

| file | rows | what it is |
|---|---:|---|
| `cohort_evidence_gap_matrix.csv` | 8 | the matrix, one row per cohort |
| `doublet_prior_capture_evidence.csv` | 605 | per-capture evidence rows |
| `capture_evidence_join_key.csv` | 746 | the materialised join key per capture |
| `doublet_prior_source_audit.csv` | 23 | the sources behind the matrix |
| `ren_gse158055_capture_qc.csv` | 284 | Ren per-sample QC |
| `pooled-same-donor-pairing-scenarios-v3.csv` | 58 | same-donor pairing scenarios |
| `liu_demuxlet_per_lane.csv` | 55 | Liu per-lane demuxlet SNG/DBL/AMB |
| `tenx_multiplet_rate_tables.csv` | 54 | vendor multiplet tables, reference data |

Verify with `shasum -a 256 -c SHA256SUMS.txt` from this directory.

`ATTACHMENT_MANIFEST.md` maps every file in this package back to the repository path it was copied
from, with its SHA-256, so the sent copy is traceable to version control. Both that manifest and
this README are generated from versioned sources by
`scripts/32_rebuild_upfront_qc/code/assemble_review_delivery.py`; nothing in this package is
hand-assembled.

## Handover status

**Not yet sent.** No external reviewer has received this package. The handover log records it as
unsent.

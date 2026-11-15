# Build provenance — sensitivity-design review packet

The PDF in this directory is a **frozen evidence package prepared for circulation to
independent methodological review**. It becomes *circulated* only once it is sent and the
handover log below is filled in. It is tracked in git, against the repository-wide `*.pdf` ignore rule,
so that "which version did the reviewer see" has an exact answer. Reproducibility is not a
substitute for a recorded correspondence between the circulated file and a commit.

**The PDF must never be edited by hand.** Any content change goes through
`scripts/32_rebuild_upfront_qc/code/build_sensitivity_review_packet.py`, followed by a new entry
in the build log below.

---

## Current build

| | |
|---|---|
| File | `sensitivity_design_review_packet.pdf` |
| SHA-256 | `1ea5507daf5775fbdf8da1c2482ec0a52e511309b3a663968680ab3602343c5a` |
| Size | 89,315 bytes |
| Pages | 12 |
| Built | 2026-09-24 (typography rebuild; evidence still frozen at 2026-09-22) |
| Supersedes | `5c08b6857f92bf4665163c4fac2a0b8967246fe3f9b55302fcdd2c631d45371b` (13 pages, 90,556 bytes, 2026-09-22) — that build ran text off the physical page edge on pages 2-4 and must not be circulated |
| Source | `README.md` and `cohort_evidence_gap_matrix.csv` **as committed alongside this PDF** — source, build script, ignore exception and artifact are one commit, so there is no window in which the PDF and its inputs disagree |
| Commit | cannot be written inside the file it names. `git log --follow -- reports/manuscript/revision/sensitivity_design_review/sensitivity_design_review_packet.pdf` gives it; `git log -1 --format=%H` at the time of sending goes in the handover log below |
| pandoc | 3.9 |
| XeTeX | 3.14159265-2.6-0.99999 (TeX Live 2018) |

## How it is built

```bash
python3 scripts/32_rebuild_upfront_qc/code/build_sensitivity_review_packet.py --pdf
```

That is the whole command. The pandoc invocation lives in `PANDOC_ARGV` inside the script and
**nowhere else** — an earlier version of this file transcribed it, and the two had already
drifted: the script's printed form was missing `title`, `subtitle` and `date`, so anyone
rebuilding from it would have produced a different title page. The script now runs pandoc
itself and prints the resulting SHA-256, size, page count and tool versions to be pasted below.

Without `--pdf` the script regenerates the matrix and validates it, but renders nothing.

### What the build actually depends on

Four layers, because "reproducible from two tracked inputs" understated it and could not explain
why a change in T7 makes the build fail:

| Layer | Contents | Tracked? |
|---|---|---|
| **Direct render inputs** | `README.md` (source of record), the generated `_packet_full.md` (README + matrix inline), `_header.tex` (page furniture) | README yes; the other two are intermediates, deleted after each build |
| **Upstream evidence inputs** | the T7 artifacts — `capture_evidence_join_key.csv`, `rate-sensitivity-scenarios-v3.csv`, `pooled-same-donor-pairing-scenarios-v3.csv`, `gate2_source_crosswalk.csv`, `doublet_prior_capture_evidence.csv`, `doublet_prior_source_audit.csv`, `ren_gse158055_capture_qc.csv`, `liu_demuxlet_per_lane.csv`, `tenx_multiplet_rate_tables.csv` — plus `docs/current-status.md`, which the classification check reads | yes |
| **Build implementation** | `build_sensitivity_review_packet.py`: the manual adjudication in `FACTS`, the expected-count constants, the twelve revalidation assertions, `MATRIX_FIELDS`, `HEADER_TEX`, `PANDOC_ARGV` | yes |
| **Build environment** | pandoc, XeTeX, and the fonts named in `PANDOC_ARGV` (Helvetica Neue, Menlo) | no — recorded by version in the table above, not pinned |

Scope of that protection, stated precisely: **a change affecting any of the twelve asserted
invariants causes the build to fail** before anything is rendered. Other upstream changes — a
reworded `count_definition`, a new column, an edited free-text field — may still alter the
rendered packet without tripping an assertion, and therefore require human review. The
assertions are a tripwire on the countable claims, not a proof that the packet is unchanged.

### Determinism: it is frozen and checksummed, not byte-reproducible

Two builds of identical inputs produce **different bytes**. Measured, not assumed:

```
without SOURCE_DATE_EPOCH   e4e765f4a20ca8e9…  vs  5bf6dbfa6b163276…
with    SOURCE_DATE_EPOCH   95587df3ae75c42f…  vs  bce2a1a0d8b9085b…
```

`SOURCE_DATE_EPOCH` does pin the PDF `CreationDate`, but something else in the XeTeX/pandoc path
still varies, so it is not used — pinning one field while the output still differs would only
make the file look more reproducible than it is.

What this artifact therefore is: **a frozen file with a recorded checksum**. The guarantee is
that the bytes circulated to a reviewer are identified exactly, not that a third party can
regenerate those same bytes. Anyone rebuilding gets an equivalent document with a different
hash; comparison should be on content, and the checksum here identifies the copy that was sent.

**The PDF must never be edited by hand.** Content changes go through the script and get a new
build-log row.

**The write is atomic.** pandoc renders to `_sensitivity_design_review_packet.inprogress.pdf`;
the tracked path is replaced only after the render succeeds and the output starts with `%PDF-`
and is a plausible size. A failed build therefore leaves the existing frozen artifact intact --
verified two ways: a pandoc that fails before writing anything, and a stub pandoc that
truncates its target, writes junk and *then* returns non-zero. In both the previous PDF's
checksum was unchanged and no temp file survived. The second case is the one a font error never
reaches, and it is what exposed a cleanup gap: the in-progress PDF was not being removed on
failure, and a first attempt to remove it deleted the file on the SUCCESS path too, because the
`finally` ran before the rename. The render, the validation and the rename now all sit inside
the same `try`, so the cleanup only ever sees a temp that did not become the artifact.

## Build log

| Date | SHA-256 (first 16) | Pages | Source | Note |
|---|---|---|---|---|
| 2026-09-24 | `1ea5507daf5775fb` | 12 | **current — circulate this one** | TYPOGRAPHY ONLY; no value, rate, cohort classification or matrix cell changed, and `cohort_evidence_gap_matrix.csv` is byte-identical to `45aeb9f`. The 2026-09-22 build ran text off the PHYSICAL page edge on pages 2-4, measured with `pdftotext -bbox` against the 612pt page width: 613.6pt (a build-script path), 615.1pt and 615.1pt (`not_eligible_without_a...`). Both causes were unbreakable words: pandoc escapes `_` as `\_`, at which TeX will not break, and the matrix tables used a short separator row so pandoc emitted non-wrapping columns. Fixed at source with a zero-width break after each underscore, wrapping table columns and a smaller face for the matrix. `hyphenat[htt]` was tried and REJECTED: it fixed the width but inserted a real hyphen inside file paths (`scripts/32_rebuild_upfront_-qc`), a transcription hazard in a document sent outside the project. After the rebuild the worst overhang anywhere is 552.7pt, inside the page on all 12 pages. 89,315 bytes. One page shorter than 2026-09-22 because the matrix is set more compactly. |
| 2026-09-22 | `5c08b6857f92bf46` | 13 | SUPERSEDED by the 2026-09-24 build | Ren's granularity corrected. Granularity is now an explicit per-cohort field; the previous substring test matched Ren's "... NOT per-capture" and filed it as capture-addressable, contradicting Q2. Q2's own chain diagram and the "cleanest per-capture case" phrase also updated. Atomic-write cleanup fixed. |
| 2026-09-22 | `2509e35960161999` | 13 | superseded | Granularity vocabulary unified in the matrix and summary table only; Ren still misclassified. Never circulated. |
| 2026-09-22 | `9492f3d0f70e0b5f` | 12 | superseded | Built by `--pdf`, the single definition of the command. Never circulated. |
| 2026-09-22 | `adece6c35c513faf` | 12 | superseded | First tracked build. Added the per-page footer. Never circulated. |

## Handover log

Recorded when the packet is sent, so the reviewer's copy is identifiable afterwards.

| Date sent | Recipient | Commit released | PDF SHA-256 | Notes |
|---|---|---|---|---|
| *(not yet sent)* | | | | |

## What the reviewer is asked to return

Per cohort: whether a scoring experiment may be designed at all; the permitted parameterisation;
the evidence grade and label it must carry; the post-hoc invalidation conditions; and an explicit
"insufficient evidence, do not proceed" option. Plus a decision on each of Q1–Q5.

Until that comes back: **no runtime configuration is created and no pilot is run.**

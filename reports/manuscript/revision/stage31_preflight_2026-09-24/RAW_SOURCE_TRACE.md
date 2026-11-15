# Stage 31 raw-count source, physical capture and raw barcode — read-only trace

**2026-09-24. Read-only. Nothing was scored, selected, submitted or deleted.**

Traces the three things the redesign needs, for the 9 `Data_NO` groups of the 302 pre-subsample
COVID `Sample_ID`s. Every row is either measured or `unknown`; nothing is inferred from a matching
cell count or a plausible-looking string.

## Headline: the container is an inherited artifact with no build script here

`01_DATA/raw/human_platelet_covid_severity_raw_counts.h5ad` is dated **2025-08-19**, which predates
every per-cohort raw download in `01_DATA/raw/` (all 2025-09). **No script in this repository
builds it.** `scripts/24_option_c_canonical_input/...` consumes it as a given input.

So the raw-count source of the Stage 31 input is **not reconstructible from this repository**. That
is the finding, not an obstacle to work around: C1 asks for un-normalised RNA counts *with source
evidence*, and the source evidence does not exist here.

`obs['orig.ident']` is `SeuratProject` for **all 47,977 cells** — the one field that would normally
name the originating object was flattened to the Seurat default, so it carries no provenance.

## What the container does tell us, measured

| `Data_NO` | cells | samples | token len | dominant cell-id shape | source attribution |
|---|---:|---:|---|---|---|
| 302004data01 | 180 | 7 | **12nt** | `GSM<n>_<n>.<n>_<BC>` | GSM accession embedded → a GEO series; **12nt = not 10x** (Seq-Well-like) |
| 302004data02 | 1,050 | 19 | 16nt | `<BC>-<n>` | unknown |
| 302004data04 | 822 | 26 | 16nt | `<n>_<BC>-<n>` | unknown |
| 302004data05 | 4,490 | 7 | 16nt | `<BC>-<n>_<n>` | unknown (bare Seurat merge suffix) |
| 302004data15 | 9,052 | 79 | 16nt | `Platelet_<BC>-<n>` | unknown; prefix is a lineage token, not a source |
| 302004data16 | 2,240 | 65 | 16nt | `<BC>-<n>_<n>` | unknown (bare Seurat merge suffix) |
| 302004data17 | 13,562 | 131 | 16nt | `<BC>-MH<n>`, `<BC>-newcastle<n>`, `BGCV<n>_<BC>-<n>` | **site and pool tokens match Stephenson E-MTAB-10026's own naming** |
| 302005data01 | 6,512 | 33 | 16nt | `Megakaryocytes_Ki<n>_Ki<n>TA<n>_<BC>-<n>`, `Coc<n>` | `Ki`/`Coc` read as Kiel/Cologne → Schulte-Schrepping; **token semantics, not a stated mapping** |
| 302005data02 | 6,331 | 24 | 16nt | `<BC>-<n>_<n>` | unknown (bare Seurat merge suffix) |

**5 of 9 groups carry no source token at all.** Their cell ids are a barcode plus a Seurat merge
index, which identifies position in a merge, not an origin.

The two attributions above are **not** equal in strength: `302004data17` embeds tokens that appear
in the deposit's own sample naming, while `302005data01`'s rests on reading `Ki`/`Coc` as city
abbreviations. Both are recorded as attributions to be confirmed against the deposits, not as
established mappings.

## Physical capture

**Unresolved for all 302 `Sample_ID`s.** Nothing in the container names a GEM well, lane or array.
The only capture-like tokens are `BGCV<n>` in `302004data17` — and those are *sequence pools*, which
the Gate 2 evidence already records as **not** GEM wells.

## Raw barcode

The nucleotide token is a **candidate**, not a verified raw barcode: it is extracted from the
canonical cell-id string, and no raw matrix or barcode file was opened.

Two measurements bound what it can be used for:

- **0** duplicate `(Sample_ID, token)` pairs — within a sample the token separates cells;
- **723 tokens appear in more than one `Data_NO` group** — so a barcode alone is **not** a global
  identity key. Any identity work must carry the sample/capture alongside it.

Token length also splits the groups: `302004data01` is 12nt and every other group is 16nt, so the
container mixes at least two barcode namespaces.

## What this means for C1–C3

| condition | status | why |
|---|---|---|
| **C1** un-normalised RNA counts with source evidence | **not satisfied** | `int64` and integral on probe, but the producing script does not exist here and `orig.ident` is flattened; provenance is absent |
| **C2** fixed physical capture + raw-barcode namespace | **not satisfied** | no capture named anywhere; 29 id shapes; two token lengths; 723 cross-group token collisions |
| **C3** complete cell-called/QC population, not pre-selected | **unresolved** | the container is named for platelets and most ids carry no lineage token, but 6,030 carry `Megakaryocytes` and a few `T_cell`/`B_cell`/`Other_cell`; those are inherited object names, not evidence of current composition |

## The next read, if this is to be resolved

Attribution has to come from the deposits, not from the container: for each `Data_NO` group, join
its tokens against the candidate raw source's own barcode list and check the overlap. That is a
real join against files in `01_DATA/raw/`, and it is **not** performed here — this trace stops at
what the container itself can support.

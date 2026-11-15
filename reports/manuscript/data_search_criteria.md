# What kind of data to look for (non-FASTQ, for consistent decontamination)

Scope decision (2026-09-14): decontaminate upfront and apply the correction the same way to
every cell in the pooled dataset. That means ambient/doublet correction from count matrices
across the board, for existing AND new cohorts alike — not CellRanger+CellBender-from-FASTQ for
some cohorts and matrix-only correction for others. Mixing two different processing pipelines
across cohorts is itself a reviewer-attackable inconsistency, worse than picking one method
uniformly. So: **only look for cohorts that ship already-aligned count matrices** — skip
anything that requires raw FASTQ + alignment to get there.

**Correction (2026-09-14):** `raw_feature_bc_matrix` is fine and should NOT be avoided —
"raw" here means unfiltered (empty droplets included), not unaligned. It is still a *matrix*
file, produced after alignment, not a FASTQ file. Having the raw/unfiltered matrix is actually
better than filtered-only where available: it enables either a CellBender-style empty-droplet
ambient estimate or DecontX on the filtered subset, whichever is used for the uniform
correction step. The only thing to actually avoid is FASTQ/raw sequencing reads that require
running alignment (CellRanger or equivalent) to produce a matrix at all.

This also means the four cohorts in the current study that happen to have open FASTQ in SRA
(GSE161918, GSE163668, GSE150728, GSE167363) should be treated the same as everything else —
DecontX on their existing filtered matrices, not CellRanger+CellBender just because the raw
reads happen to be available for them and not others.

---

## What "the right kind of data" looks like

Look for a GEO/ArrayExpress/CELLxGENE/Zenodo deposit whose **supplementary files** are one of:

- A 10x-style triplet: `matrix.mtx.gz` + `barcodes.tsv.gz` + `features.tsv.gz` (or `genes.tsv.gz`)
  — either `filtered_feature_bc_matrix` or `raw_feature_bc_matrix` are both fine; the raw/
  unfiltered version is actually preferred when available, since it includes empty droplets and
  supports a stronger ambient-RNA estimate
- A single `.h5` (10x HDF5 format) — `filtered_feature_bc_matrix.h5` or `raw_feature_bc_matrix.h5`,
  same as above
- A processed `.h5ad` (AnnData) or `.rds`/`.Rdata` (Seurat/SingleCellExperiment object) with
  counts already in `adata.X` / `assay(sce, "counts")` — this is the easiest case, often already
  merged across samples with metadata attached
- A `.loom` file with a counts layer

**Avoid / do not pursue for this track:**
- SRA/ENA accessions (`SRPxxxxxx`, `PRJNAxxxxxx`) or `.fastq.gz` files as the *only* access
  route — these are raw sequencing reads that require running alignment (CellRanger or
  equivalent) to produce any matrix at all, which is the actual inconsistency to avoid
- DGE/count files that are ambiguous about filtering status — check the methods section of the
  source paper so you know whether what you have is raw or filtered before using it, but either
  is usable

**Quick way to check a GEO series**: open the series' "Supplementary file" table — if you see
`*_feature_bc_matrix.h5` (raw or filtered) or a `barcodes/features/matrix` triplet, it's a
match, whichever filtering state it's in. If the only files are `.sra`/`.fastq.gz` or the GEO
record links out to SRA with no matrix file listed at all, skip it for this track (it's fine for
a future FASTQ-based rebuild, just not this one).

---

## Candidate cohorts, filtered to processed-matrix access

From the earlier cohort inventory, tiered by what's actually needed here:

| Cohort | Access | Status |
|---|---|---|
| GSE158055 (Ren) | Filtered matrices, BioProject only (no SRA relation) | **Use** — already in the study at 54/196 donors; recovering more donors here needs only outcome-label research, not new access |
| GSE279452 | Per-sample mtx triplets (666 files) | **Use** — 91 sepsis donors, largest sepsis arm, filtered matrices confirmed |
| GSE208337 | Per-sample mtx triplets | **Use** — 105 samples, Myl9/microthrombi, directly platelet-relevant |
| GSE174072 | rds count matrices | **Use** — processed R object, 70 samples |
| GSE303333 | Per-sample mtx triplets | **Use** — pediatric sepsis, 24 samples, small but processed-matrix |
| GSE151263 | Processed UMI txt | Already excluded from the current study (no per-patient outcome) — skip |
| GSE161918 (Liu) | Per-sample matrices + Seurat rds, *also* has open FASTQ in SRA | **Use the rds/matrix route only** — ignore the FASTQ option for this track |
| GSE163668 (Combes) | Per-sample matrices, *also* open FASTQ in SRA | **Use the matrix route only** |
| GSE150728 (Wilk) | rds count matrices, *also* open FASTQ in SRA | **Use the matrix route only** |
| GSE167363 (Qiu sepsis) | Per-sample mtx triplets, *also* open FASTQ in SRA | **Use the matrix route only** — already in the study this way |
| GSE180578 | RAW.tar + annotated h5ad | **Use either** — the RAW.tar matrix bundle and the annotated h5ad are both matrix data, not FASTQ |
| E-MTAB-10026 (Stephenson) | Processed matrices (already in use) | **Use** — recovering more donors (90→120) needs only outcome-label research |
| Schulte-Schrepping | EGA, processed | **Use if EGA access is pursued** — processed, not raw |
| GSE165080, GSE196198 | Series mtx triplet / per-celltype matrices | **Use** — both processed-matrix |
| GSE290679 | h5ad (raw + filtered) | **Use either** — both are matrix data; raw+filtered both being offered is a good sign for ambient-RNA work |
| COMBAT (EGAS00001005493) | Processed on CELLxGENE + Zenodo, controlled access | **Use** — CELLxGENE deposit is a processed h5ad; the EGA route is for the controlled-access version, still processed, not FASTQ |
| Teichlab paediatric+adult (EGAD00001007718) | Processed on CELLxGENE, controlled access | **Use** — same pattern as COMBAT |

Every candidate above already clears the "matrix, not FASTQ" bar. The four SRA-listed
ones (Liu/Combes/Wilk/Qiu-sepsis) are flagged specifically so the FASTQ option is not
accidentally taken for those just because it happens to be available — pull the matrix/rds
files instead, same as everything else, for pipeline consistency.

## Still your research task for each

1. **Per-patient outcome (mortality) label** — in the paper, its supplement, GEO sample
   characteristics, or by contacting the authors. This is the binding constraint, not access to
   cells.
2. **Note whether the matrix is raw (unfiltered) or filtered (cell-called)** for each cohort —
   both are usable, but the correction script needs to know which it's getting, and raw/
   unfiltered is preferred where both exist since it supports a stronger ambient estimate.
3. **Platelet annotation or markers present** — either the deposit already labels platelets, or
   the gene panel includes PPBP/PF4/GP9/ITGA2B/TUBB1/GP1BA so a marker-based call is possible
   post-hoc, same approach already used for this project's own cohorts.

Build the event ledger as you go (cohort → donors → confirmed fatal events, running total),
same target as before: **≥90 pooled fatal events** to meaningfully improve on the current 25
COVID + 17 sepsis = 42.

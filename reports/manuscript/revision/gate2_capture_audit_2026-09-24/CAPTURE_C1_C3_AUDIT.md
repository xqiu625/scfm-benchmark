# Gate 2 capture audit — one capture, C1–C3 with row-level source evidence

**2026-09-24. Read-only. No job submitted, no cell deleted, no outcome read,
`PRODUCTION_PARAMETER_RULING` not filled.**

## The capture

**`GSM5102900_HC1`** — Qiu, GSE167363. Exactly one capture; no other cohort is touched.

Chosen because it is the case most likely to *close*, and the reasons are checkable:

- the deposit ships a per-GSM CellRanger triplet, so C1/C2/C3 can be evidenced at row level;
- its name matches the canonical container exactly (6 of Qiu's 12 GEO names do; `P25/P50/NSES`
  vs `NS1/NS2/S1` do not, so a mismatched sample would have made the identity join ambiguous);
- one donor, one draw, no timepoint ambiguity;
- a healthy control, so there is no outcome coupling of any kind.

## C1 — un-normalised raw-count provenance: **PASS**

Evidence is the deposit's own statements, not inference:

| evidence | value |
|---|---|
| MatrixMarket header | `%%MatrixMarket matrix coordinate **integer** general` |
| dimensions / nnz | 33,538 genes × 7,666 cells, 9,903,096 non-zeros (streamed and counted; declared = observed) |
| values | integral, min 1.0 |
| `features.tsv.gz` feature types | `{'Gene Expression'}` — **only**; no Antibody Capture |
| GEO `!Sample_library_strategy` | `RNA-Seq` |
| GEO `!Sample_library_source` | `transcriptomic` |
| GEO `!Sample_data_processing` | "Cell Ranger software pipeline (v3.1.0) … We obtained the **feature-barcode matrix that contains gene expression counts**" |
| genome | GRCh38, `refdata-cellranger-GRCh38-3.0.0` |

File identities — full SHA-256, immutable paths and the GEO source for each:

```
GSM5102900_HC1_barcodes.tsv.gz      39,604 bytes
  sha256 e6dbf9f587d9a6cf3ecff31cbf4d53b09c6faf9837c82d5c1c0a26227f377a02
GSM5102900_HC1_features.tsv.gz     304,728 bytes
  sha256 c7408e82af270ca6308a8aa7ea659795a15ab76f2e3a85bd5db4cd3d53244199
GSM5102900_HC1_matrix.mtx.gz    38,785,404 bytes
  sha256 4ea259bfc4de02cc9799df563779fce7c9bd4c889532d3ff1178e2caf545b9f8
```

| | |
|---|---|
| local path | `/bigdata/godziklab/shared/Xinru/302006/01_DATA/raw/GSE167363_qiu/extracted/` |
| GEO series | <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE167363> |
| GEO sample | <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5102900> |
| GEO supplementary | `ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5102nnn/GSM5102900/suppl/` |
| canonical container | `/bigdata/godziklab/shared/Xinru/302006/01_DATA/raw/human_platelet_covid_severity_raw_counts.h5ad`, 5,027,039,977 bytes (not hashed; identity is path + size) |

**Replayable.** `verify_capture_evidence.py` recomputes every number in this document —
dimensions, nnz, integrality, feature types, barcode count/uniqueness/length/suffix distribution,
the 65/65 join and the retained fraction — and **exits non-zero if any of them disagrees**. Its
output is `EVIDENCE_MANIFEST.json`. The numbers below are therefore checkable, not asserted.

The matrix declares itself integer, the depositor names it a CellRanger feature-barcode count
matrix, and the feature file contains no non-RNA modality. This is the un-normalised RNA count
source C1 asks for, with the source stated rather than assumed.

## C2 — physical capture and raw-barcode namespace: **PASS**

| evidence | value |
|---|---|
| barcodes | 7,666, **all unique** |
| token length | **16 nt for every barcode** — one namespace, no mixing |
| GEM-well suffix | **`-1` for all 7,666** |
| processing | GEO `!Sample_data_processing`: **`cellranger count`**, Cell Ranger v3.1.0 |
| deposited files | one per-GSM triplet (GEO `supplementary_file_1..3`) |
| chemistry | `Chromium Next GEM Single Cell 3ʹ Reagent Kits v3.1` (GEO) |
| instrument | Illumina NovaSeq 6000 |

The uniform `-1` is the load-bearing observation. 10x documents the barcode suffix as the **GEM
well identifier**, and it is `-1` throughout a single-GEM-well sample; an aggregated matrix carries
`-1 … -N`. So a single `-1` across all 7,666 barcodes is positive evidence for **one GEM-well
barcode namespace**.

Stated precisely, because the earlier wording overreached: GEO reports `cellranger count` for this
sample and provides one per-GSM triplet; together with the uniform `-1` suffix, that is evidence
for a single GEM-well namespace. **The file count alone is not the proof** — a triplet is simply
how CellRanger writes one matrix, and on its own it would say nothing about how many GEM wells
went into it.

**Recorded, not resolved:** the deposit states 3′ v3.1 while the project record notes the paper
describes 5′. The ruling holds this blocks any implementation depending on a vendor rate; it does
not block a pure algorithmic-stress design, which is what the approved grid is.

## C3 — complete declared-stage population and parent-state support

**PASS for the deposit. FAIL for the canonical container.** The distinction is the whole result.

| input | cells | verdict |
|---|---:|---|
| deposit `GSM5102900_HC1` | **7,666** | complete CellRanger-called population at the declared stage, all lineages |
| canonical container, `Sample_ID=HC1` | **65** | **0.85%** of the deposit — platelet-selected |

Row-level join, container → deposit:

```
container HC1 cells                      65
deposit barcodes                      7,666
container tokens found in deposit     65/65   (100.0%)
container is a subset of the deposit   True
fraction of the deposit retained       0.85%
```

The identity join is clean — every container cell is a real deposit barcode, so the mapping
itself is sound. What fails is completeness: 0.85% is a lineage-selected slice, and C3 forbids an
input "pre-selected by target lineage".

Note C3's own wording: *complete* means the available cells and parent states of the known stage,
**not** feeding empty droplets. The CellRanger-called set is exactly that stage, so the deposit
satisfies C3 while the container cannot.

**Consequence:** a scoring experiment on this capture must take the **deposit's 7,666-cell
matrix** as input. The canonical container is not an admissible input for it.

## Verdict

| condition | verdict |
|---|---|
| C1 un-normalised raw counts with source evidence | **PASS** |
| C2 physical capture + raw-barcode namespace | **PASS** |
| C3 complete declared-stage population, not lineage-preselected | **PASS (deposit) / FAIL (container)** |

C1–C3 are closable for this capture **on the deposit input**.

## Where this audit stops, and why

Instruction 4 asks the config audit material to include **strict PC feasibility**. It is not
complete:

```
strict rule        n_prin_comps < min(n_cells_effective, n_usable_features)
n_cells_effective  7,666   MEASURED (filter_cells(min_genes=3): 7,666/7,666 pass)
genes passing      17,226  MEASURED (filter_genes(min_cells=3): 17,226/33,538)
n_usable_features  NOT MEASURED
approved value     n_prin_comps = 30
verdict            CANNOT BE CONFIRMED THIS ROUND
```

`n_usable_features` is the highly-variable-gene count, which requires materialising the matrix and
running the HVG step. That attempt was **OOM-killed on the login node** — the documented failure
mode for real inputs — and a SLURM job is not authorised this round.

30 < 7,666 holds on the cell side, and 30 < 17,226 holds on the pre-HVG gene side, so the
constraint can only bind through HVG. It is very unlikely to bind. **It is still not measured, so
it is recorded as unconfirmed rather than assumed** — that is exactly the substitution instruction
3 forbids.

The remaining config-audit items (config hash, the eleven approved values bound to this capture,
dependency versions, fixed seed and threshold rule, and evidence that the simulated parent pairs,
PCA and neighbour graph are identical across the eight priors) are **not produced here**, because
they would rest on a feasibility check that has not been completed. They are ready to produce in
one authorised job that measures HVG on this single 7,666-cell capture.

## Boundaries honoured

No SLURM job submitted. `PRODUCTION_PARAMETER_RULING` untouched (`{}`). No outcome field read —
the capture is a healthy control and only `!Sample_characteristics` text was read from GEO
metadata, never a project outcome column. No cell deleted. No doublet score, call or rate
produced. Exactly one capture examined.

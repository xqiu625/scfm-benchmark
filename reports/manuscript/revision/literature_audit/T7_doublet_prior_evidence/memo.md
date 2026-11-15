# Capture/library-level doublet-rate loading priors: evidence survey

**Scope.** Whether a defensible capture/library-level loading prior can be obtained from public
sources for eight cohorts. This is an evidence survey only: it does not select a prior, derive a
removal quota, or use outcome data. Every factual line below is tied to a source in
`doublet_prior_source_audit.csv`; the per-capture rows are in
`doublet_prior_capture_evidence.csv`.

**Stage vocabulary** is kept strictly separate throughout: `target_cells_loaded`,
`raw_droplet_barcodes`, `estimated_or_recovered_cells`, `called_cells_pre_study_qc`,
`called_cells_post_study_qc`, `post_demultiplex_singlets`, `post_doublet_filter_cells`,
`lineage_subset_cells`. Where a source gives a number at one stage, it is not promoted to another.

---

## 1. Executive summary

### Qiu — GSE167363 — VALIDATED (loading), PROXY_ONLY (per-capture cells)

Best evidence: Qiu et al., *J Leukoc Biol* 2021 Methods — "A total of 15,000 cells per sample were
loaded onto the 10x genomics platform… (Chromium Next GEM Single Cell V3.1)", with Cell Ranger
v3.1.0 `cellranger count`. The deposit ships per-sample barcodes/features/matrix triplets.

Supports a loading-based prior: **yes**, at a uniform protocol level (15,000 loaded per sample);
the corresponding recovered-cell figure is not reported, so the recovered-cell entry of the 10x
table has to be reached either through the deposited called-cell count or through the guide's
loaded→recovered column, not both at once. Deposited called cells per sample, counted from the
barcode files: median 5,136, range 2,853–8,345 over 12 samples, 12 distinct values — natural
variation consistent with default Cell Ranger calling rather than a forced count.

Largest uncertainty: a chemistry conflict. GEO's extract protocol states "Chromium Next GEM Single
Cell 3ʹ Reagent Kits v3.1"; the paper's processing section says Cell Ranger performed "single-cell
5′ unique molecular identifier (UMI) counting". Both cannot be right, and the 3′ v3.1 and 5′ v1.1
tables differ slightly. Unresolved from public sources.

### Chongqing — GSE279452 — UNAVAILABLE (loading), PROXY_ONLY (per-capture cells)

Best evidence: GEO only. NCBI GDS lists **no linked PubMed ID** for this series, so there is no
paper or supplement to search. The extract protocols name 10x Single Cell 3′ v3.1 and 5′ VDJ kits;
data processing states "Cellranger was utilized to calculate the raw counts and cell whitelist by
default parameter" — i.e. default (data-driven) cell calling, not a forced count. Per-sample
barcode files give a called-cell count per capture.

Supports a loading-based prior: **no** for loading; the deposit supports a called-cell number only
— counted from the barcode files: median 7,417, range 2,192–15,654 over 133 samples, with 131
distinct values, again consistent with default (not forced) cell calling.

Largest uncertainty: no publication means no loading figure can be obtained at all from the
specified source hierarchy, and the Cell Ranger version is not stated.

### Cillo — GSE180578 — UNAVAILABLE (loading), PROXY_ONLY (per-capture cells)

Best evidence: *Cell Rep Med* 2021 STAR Methods — "Two samples (each stained with unique
TotalSeq-C antibodies) were then pooled and loaded per lane of the 10X chip"; GEO adds "loaded on a
10X Genomics Chip A for use with the 5prime V1 chemistry kit". **No cell number is given for the
load.** Cell Ranger v3.1.0; deposited GEX matrices are filtered feature-barcode matrices
(verified by inspection: GSM5465125 contains 4,038 barcodes over 33,538 features).

Supports a loading-based prior: **no**. Only called cells per library are obtainable.

Largest uncertainty: 27 of 48 GEX libraries carry two samples and 10 pool two different donors, so
any prior must account for cross-donor doublets; the authors removed only HTO double-positive
cells, and no doublet rate is reported.

### Liu — GSE161918 — INELIGIBLE (the 30,000 figure), UNAVAILABLE (loading)

This is the cohort where the distinction the brief asks for matters most, and the answer is
unambiguous. The per-lane deposits are named `*_filtered_feature_bc_matrix.h5`, but they are **not
data-driven cell calls**: **all 56 GEX lanes contain exactly 30,000 barcodes** (verified by opening
every deposited h5 and reading `/matrix/shape`), and for lane `B1_10xlane1` the h5 barcode set and
the demuxlet `.best` barcode set are **identical (30,000 of 30,000 overlapping)**. The `.best` files carry 30,000 rows per lane (51 lanes exactly 30,000, three
29,999, one 29,998). That is the signature of a fixed candidate set (a forced cell count), used
both for the deposited matrix and as the demuxlet input.

So, against the definitions in the brief: the 30,000 is a **fixed candidate whitelist**, not
`estimated_or_recovered_cells` and not `called_cells_pre_study_qc`. The SNG/DBL/AMB counts are
`post_demultiplex_singlets` and cross-donor multiplet calls computed on that whitelist.

- Per-lane target loading: **not reported** in the paper, the GEO protocols, or the series design.
- Per-lane recovered-cell estimate: **not obtainable** — no raw droplet matrix is deposited (only
  the forced 30,000-barcode h5 files), so cell calling cannot be re-run from public data.
- Chemistry is recorded in the h5 attributes as "Single Cell 5′ R2-only" (Next GEM 5′ v1.1 per the
  GEO protocol), Cell Ranger 3.1.0 per GEO.

Supports a loading-based prior: **no**. Any per-lane number currently available is either a fixed
whitelist or post-demultiplex.

### Ren — GSE158055 — PROXY_ONLY

Best evidence: Table S1 of the *Cell* 2021 paper has a dedicated sheet, "Quality control metrics
for individual sample", with per-sample `inflection_UMI`, `Ncells_Original`, then stepwise
`Ncells_gene_lower/higher`, `Ncells_UMI_lower/higher`, `Ncells_Mt`, `Ncells_Scrublet`,
`Ncells_manual_removal` for all 284 samples (median `Ncells_Original` 6,566; range 364–35,616;
sum 2,001,714 → 1,462,702 after manual removal, matching the paper's headline count).

Two stage caveats that block promotion of these numbers:

1. Cell calling was **not** Cell Ranger. Gene expression was quantified with kallisto|bustools and
   cells were called at a per-sample UMI inflection point, so `Ncells_Original` is a
   `called_cells_pre_study_qc` under a different cell-calling definition than the 10x tables assume.
2. Loading is reported only as a cohort-level range: "The loaded cell numbers range from
   300-500,000 aiming for 300-14,000 single cells per reaction." That is not a per-capture number
   and must not be spread across captures.

The deposited series matrix is post-QC **and** post-author-doublet-removal (Scrublet with
`expected_doublet_rate = 0.08` cohort-wide, plus UMI/gene upper cutoffs).

Largest uncertainty: chemistry is mixed within the cohort — 223 samples 10x 5′, 61 samples 10x 3′
(Table S1 / GEO characteristic), spanning 3′ v2, 3′ v3 and 5′ v2 reagent kits, so a single
chemistry table does not apply cohort-wide.

### COMBAT — VALIDATED (loading), UNRESOLVED (model applicability)

Best evidence, verbatim from the *Cell* 2022 STAR Methods: cells were "counted on the BioRad TC20
Automated Cell Counter and loaded onto the 10X Genomics Chip G at **50,000 cells per channel**.
Each pool was loaded across seven channels." The Zenodo CITE-seq object gives `Channel_ID` and
`Pool_ID` per cell (10 pools × 7 channels = 70 GEM wells, 14 donors per channel).

Supports a loading-based prior: **the loading number is the strongest of any cohort here**, but the
10x expected-rate tables do not cover it — 50,000 cells loaded is roughly three times the top row
of the Next GEM tables (16,500 loaded / 10,000 recovered) and above the GEM-X v4 table's top row
(29,000 / 20,000). Applying either table at this load requires extrapolation the guides do not
sanction.

Other stages: recovered/called cells per channel are **not public**. The paper reports only the
post-QC total (`ngenes > 300`, `pct_mito < 10%`, n = 836,148 cells), and the Zenodo object is that
post-QC set. Cross-donor demultiplexing used demuxlet v2 and Vireo; doublets were screened using
Scrublet scores plus cross-modality congruence (a ~26k-cell doublet group was identified during
clustering), but no doublet rate per channel is reported.

### Stephenson — E-MTAB-10026 — PROXY_ONLY, capture mapping UNRESOLVED

Best evidence: *Nat Med* 2021 Methods gives site-specific loading — Newcastle/London
"Approximately 20,000–30,000 cells per sample were loaded onto the 10x Chromium controller"
(NextGEM 5′ v1.1 with Feature Barcoding), and Cambridge pooled four samples (0.5 M live cells each)
with "Two lanes of 25,000 cells were loaded per pool". Cell Ranger v4.0; 1,141,860 cells sequenced,
781,123 after QC (cohort level). Doublet handling differs by site: Scrublet per sample
(Newcastle), emptyDrops + genotype doublet tagging (Cambridge), souporcell v2 (London).

The ArrayExpress SDRF resolves samples to pools — 143 samples across 98 `sequence pool` values (83
pools of one sample, 15 pools of four) — but **not to GEM wells**: the Cambridge pools each span
two lanes, and which cells came from which lane is not recorded publicly.

Supports a loading-based prior: partially — loading is known per sample or per pool by site, but
the capture (GEM well) is not identifiable for the multiplexed pools, so a per-capture prior cannot
be assigned without assuming an even split, which the brief forbids.

Largest uncertainty: a chemistry conflict. The paper describes 5′ NextGEM v1.1 V(D)J kits with
Feature Barcoding; the SDRF annotates all 143 samples as `Comment[library construction] = 10xV3`.

### Wilk — GSE174072 — VALIDATED (loading), platform model available but not from 10x

Wilk's scRNA-seq arm is **Seq-Well, not 10x**, so the 10x tables are inapplicable by construction.
The loading density is documented identically in the paper and in GEO: "Immediately after Ficoll
separation, 50,000 PBMCs were resuspended in RPMI + 10% FCS at a concentration of 75,000 cells/ml.
200 µl of this cell suspension (**15,000 cells**) was then loaded onto Seq-Well arrays preloaded
with mRNA capture beads (ChemGenes)."

The authors report **no** doublet/multiplet rate. Their multiplet handling is a heuristic filter:
"To remove putative multiplets (where more than one cell may have loaded into a given well on an
array), cells that expressed >75 genes per 100 UMIs were also filtered out."

Platform-level expected-rate model: **it exists** — see section 5. Deposits are
`*_cell.counts.matrices.rds.gz` per array (45 scRNA deposits; the 25 fragment files belong to the
10x scATAC arm, which is a different platform and different chemistry, Next GEM ATAC v1.1).

---

## 4. 10x expected-multiplet-rate model — conclusions

Extracted directly from four official 10x user guides (tidy version in
`tenx_multiplet_rate_tables.csv`):

| document | chemistry | table range (recovered) | rate per 1,000 recovered | recovered ÷ loaded |
|---|---|---|---|---|
| CG000204 Rev D | Next GEM Single Cell 3′ v3.1 | 500–10,000 | 0.76–0.80 % | 62.5 % |
| CG000315 Rev B | Next GEM Single Cell 3′ v3.1 dual index | 500–10,000 | 0.80 % exactly | 60.6 % |
| CG000208 Rev G | Next GEM Single Cell V(D)J v1.1 (5′, Feature Barcode) | 500–10,000 | 0.76–0.80 % | 56.6–57.3 % |
| CG000731 Rev A | GEM-X Single Cell 3′ v4 | 500–20,000 | 0.40 % exactly | 69.0 % |

**1. What is the denominator?** The guides tabulate three columns together — Multiplet Rate (%),
\# of Cells Loaded, \# of Cells Recovered — so each rate is defined against a
(loaded, recovered) pair. The rate scales linearly with **cells recovered**; "cells loaded" is the
input that produces that recovery at the guide's assumed capture efficiency. It is not defined
against called cells after user QC, nor against droplets. Zheng et al. 2017 (the platform paper)
states the same relationship empirically: a titration "showed a linear relationship between the
multiplet rate and the number of recovered cells ranging from 1,200 to 9,500", consistent with
Poisson loading, with ~50% capture (cells detected ÷ cells loaded).

**2. Is 0.8% per 1,000 an exact formula?** No — it is an **empirical linear approximation** that
happens to reproduce the CG000315 table exactly (0.4% at 1,000 recovered … 8.0% at 10,000) and the
CG000204/CG000208 tables to within rounding (0.76–0.80%). Every value in the guides is written with
a "~". There is no closed-form formula in the user guides; there is a table.

**3. Do chemistries differ?** Yes, materially. GEM-X Single Cell 3′ v4 is **0.4% per 1,000
recovered — half the Next GEM rate** — and its assumed capture efficiency is higher (69% vs
57–63%). Using the 0.8%/1,000 figure on a GEM-X v4 library would double the expected rate.

**4. Is extrapolation beyond the table allowed?** The guides give no basis for it. Next GEM tables
stop at 10,000 recovered / ~16,500 loaded; GEM-X v4 stops at 20,000 recovered / 29,000 loaded.
Superloaded captures in this set — COMBAT at 50,000 cells per channel, Stephenson Cambridge at
25,000 per lane, Wilk/Qiu at 15,000 (Qiu on 10x) — sit at or beyond the top of the applicable
table. Zheng et al. validated linearity only over 1,200–9,500 recovered cells.

**5. When can a filtered matrix cell count proxy for recovered cells?** Only when the filtered
matrix is the unmodified output of the standard cell caller for that pipeline — no `--force-cells`
or `--expect-cells` override, and no user filtering applied before deposit. In this set that
condition **fails for Liu** (fixed 30,000 barcodes per lane, identical to the demuxlet candidate
set) and **fails for Ren** (kallisto|bustools inflection call, then study QC and Scrublet removal
before deposit). It is plausible but unverified for Cillo, Qiu and Chongqing, where the deposits
look like standard Cell Ranger filtered output; the diagnostic is whether per-capture counts vary
naturally rather than sitting at a round constant.

---

## 5. Seq-Well conclusion (Wilk)

A platform-appropriate prior **is** available, and it is not the 10x number.

- **Version.** Wilk et al. cite both Gierahn et al. 2017 and Hughes et al. 2020 for the method;
  Hughes 2020 is the Seq-Well S³ protocol paper.
- **Array capacity.** The Seq-Well protocol (Protocol Exchange companion, 10.1038/protex.2017.006a)
  describes "a PDMS array of 86,000 subnanoliter wells", with single-bead loading efficiency 95%
  and cell capture efficiency ~80%, and states explicitly that "dual occupancy rate … can be tuned
  by adjusting the number of cells loaded". Hughes 2020 instead writes of "a loaded Seq-Well device
  containing 85,000 wells". **This 86,000 vs 85,000 discrepancy is left unresolved here**; both are
  primary statements from the platform's own authors.
- **Expected-rate model.** Hughes 2020 (STAR Methods, "Doublet Removal") is the usable anchor:
  "Doublet removal was performed for each array individually using DoubletFinder. For each array,
  the expected doublet rate was estimated based on the cell loading density. A total of 20,000 cells
  were loaded to a loaded Seq-Well device containing 85,000 wells, which resulted in an expected
  doublet rate of 2.37%." That is a published loading-density → expected-rate mapping with an
  explicit worked example, i.e. the platform analogue of the 10x table.
- **Applicability to Wilk.** Wilk's loading density is documented and uniform (15,000 cells per
  array), so the Hughes mapping has the input it needs. What is *not* available for Wilk: any
  per-array occupancy measurement, any author-reported doublet rate, and confirmation that every
  array received the nominal 15,000 cells.

Stated plainly for the record: **do not borrow the 10x 0.8%-per-1,000 figure for Wilk.** The
Seq-Well geometry (fixed well count, gravity loading, tunable dual occupancy) is a different
process from Poisson GEM encapsulation, and the platform's own literature supplies its own model.

---

## 6. Missing evidence

Resolvable only with author contact, controlled-access data, or local HPCC files:

1. **Cell Ranger `metrics_summary.csv` is not deposited for any cohort in this set.** Its
   "Estimated Number of Cells" field is precisely the quantity the 10x table is indexed on; one CSV
   per capture would move Cillo, Qiu, Chongqing and Stephenson from PROXY_ONLY toward VALIDATED.
2. **Liu**: confirmation that `--force-cells 30000` (or equivalent) was used; per-lane target
   loading; any retained raw droplet matrices. None of this is public.
3. **COMBAT**: per-channel recovered/called cells before QC — likely present in the EGA-controlled
   Cell Ranger outputs (EGAD00001007963), not in the open Zenodo object.
4. **Stephenson**: the sample → GEM-well assignment for the 15 four-sample Cambridge pools (each
   split across two lanes), and per-lane loaded counts.
5. **Ren**: per-sample loaded/target cells (only a cohort-wide range is published), and whether any
   sample deviated from the stated range.
6. **Chongqing**: any manuscript or preprint; Cell Ranger version; loading protocol.
7. **Wilk**: per-array cell counting/occupancy records, and which array count (85,000 vs 86,000) is
   correct for the arrays used.
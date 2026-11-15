# Stephenson/COVIDPBMC author-code investigation — Fig.2f/Table 7, GO:0034340, clinical fields

Per `codex_output/20261015_114223.md` §5 item 3 and `codex_output/20261015_104226.md` §四 item 3:
a targeted, read-only investigation of the authors' own public repository,
[`scCOVID-19/COVIDPBMC`](https://github.com/scCOVID-19/COVIDPBMC), to resolve (or precisely scope
the remaining gap in) five questions about the Fig. 2f / Supplementary Table 7 myeloid
IFN/TNF/JAK-STAT scoring analysis. No author pipeline was run; no file in this project's own
expression results was consulted while reading the author code, so no grouping was picked by
looking at which one reproduces our own numbers.

**Repository state inspected**: default branch `master`, pinned to commit
`9f51fe9a9ffb45bf3fbd4c302fe6e3f8979f5deb` (the tip of `master`; `pushed_at: 2021-12-21T10:42:15Z`,
confirmed via `curl -s https://api.github.com/repos/scCOVID-19/COVIDPBMC/commits?per_page=1`). The
full recursive tree was listed via the GitHub Trees API
(`.../git/trees/9f51fe9a9ffb45bf3fbd4c302fe6e3f8979f5deb?recursive=1`); file contents fetched via
`raw.githubusercontent.com/.../9f51fe9a9ffb45bf3fbd4c302fe6e3f8979f5deb/...`; repo-wide text search
via `gh api search/code -f q='<term> repo:scCOVID-19/COVIDPBMC'` (authenticated `gh`, since the
unauthenticated GitHub code-search endpoint returns 401).

## 1. The IFN/TNF/JAK-STAT scoring notebook — found, read in full

**File**: `myeloid/myeloid_signature_scoring_publish.ipynb`, cells 0–9 (sha256 of the fetched file
at this commit: `6f5b470f044b22b1ff8cd89f2d18707400c0b3f4b0b33dd2931f7ee5c9089211`).

- **Cell 2**: subsets the full `covid_total` object to myeloid `full_clustering` categories
  `CD83_CD14_mono, CD14_mono, CD16_mono, C1_CD16_mono, DC1, DC2, DC3, ASDC, pDC, DC_prolif,
  Mono_prolif` — **all sites, no site filter applied in this cell.**
- **Cell 5**: concatenates this blood myeloid subset with the reannotated BAL dataset
  (`GSE145926`) into one `myeloid` object — confirms Codex's own observation that Table 6 mixes
  peripheral-blood and BAL labels.
- **Cell 7**: three gene lists are hardcoded verbatim in the notebook, each with its own MSigDB
  card URL in a comment:
  - `TNF_genes` (200 symbols) — comment cites
    `https://www.gsea-msigdb.org/gsea/msigdb/cards/HALLMARK_TNFA_SIGNALING_VIA_NFKB.html`.
  - `JAK_genes` (87 symbols) — comment cites
    `https://www.gsea-msigdb.org/gsea/msigdb/cards/HALLMARK_IL6_JAK_STAT3_SIGNALING`.
  - `IFN_genes` (99 symbols) — comment cites
    `https://www.gsea-msigdb.org/gsea/msigdb/cards/GO_RESPONSE_TO_TYPE_I_INTERFERON` — **this is
    the GO:0034340-derived set the paper's Methods p.15 names**, frozen as a literal Python list
    in this cell, not a live GO/MSigDB API call at analysis time.
- **Cell 8**: `sc.tl.score_genes(myeloid, <gene_list>, ctrl_size=50, gene_pool=None, n_bins=25,
  score_name='<X>_score', random_state=0, copy=False, use_raw=None)` for each of the three sets —
  confirms the paper's own scoring method exactly as Methods p.15 describes (`tl.score_genes`).
- **Cell 9**: exports `myeloid.obs[['full_clustering', 'severity', 'TNF_score', 'JAK_score',
  'IFN_score']]` to `TNF_IFN_JAK_scores_v2.csv` — **cell-level** scores, one row per cell, not yet
  aggregated or tested.

**Finding — GO:0034340's actual gene file, resolved**: the paper used a **frozen 99-symbol list**
checked into this notebook, not a live GO-annotation query at analysis time, and not today's
current MSigDB membership for the same GO term (not independently re-pulled this round — the
notebook's own hardcoded list is the authoritative source for what the paper actually scored).
**This list is a different gene set from the project's currently frozen
`HALLMARK_INTERFERON_ALPHA_RESPONSE_msigdb.json`**: direct symbol-set comparison (both loaded
locally) gives only **24 of 99 GO-set symbols overlapping with the 97-symbol Hallmark set** (e.g.
`HLA-A/B/E`, `IRF1-9`, `MAVS`, `TBK1`, `JAK1` are GO-set-only; `B2M`, `BATF2`, `CASP1/8`, `CD74`,
`CMPK2` are Hallmark-only). These are related but clearly non-identical interferon-response gene
sets, confirming the review's caution that the currently frozen Hallmark asset is not a precise
reproduction of what Fig. 2f/Table 7 actually scored.

## 2. The Healthy/Mild/Severe three-group recoding and the Table 7 statistical test — NOT found in this repository

Repo-wide authenticated code search (`gh api search/code`) for the terms that would appear in a
script computing Table 7's reported `Δ` and `p adj` values:

| Search term | Hits in `scCOVID-19/COVIDPBMC` |
|---|---|
| `IFN_score` | 1 — `myeloid/myeloid_signature_scoring_publish.ipynb` only (the export cell above, no test) |
| `TNF_IFN_JAK` | 1 — same file, same export cell |
| `GO_RESPONSE_TO_TYPE_I_INTERFERON` | 1 — same file, same gene-list comment |
| `wilcox.test` | 2 — both in `TcellAnalysis/notebooks/`, unrelated to myeloid IFN scoring |
| `p.adjust` | 4 — all in `TcellAnalysis/notebooks/`, none in `myeloid/` |
| `Fig2f`, `Fig_2f` | 0 |

The notebook that reads `TNF_IFN_JAK_scores_v2.csv` back in, collapses the six
`Status_on_day_collection_summary`/`severity` levels into the three shown in Fig. 2f
(Healthy/Mild/Severe), and performs the statistical test behind Table 7's `Δ`/`p adj` columns is
**not present in this public repository** — neither in `myeloid/R_plotting_publish.ipynb` (read in
full, 17 cells; it plots Fig. 1C/2B/2J proportions and a PAGA pseudotime figure, no IFN/TNF/JAK
cells) nor anywhere else the code search reached. **Single remaining question**: which script (not
checked into this repo, or run outside version control) performs the 6→3 severity collapse and
the Wilcoxon/t-test-with-FDR-correction that produced Table 7 — this cannot be resolved from the
author's committed code alone, and the next step (if needed) is a direct question to the
corresponding author or Nature Medicine's code-availability contact, not a guess from this
project's own data.

A related, partially-resolved sub-question: cell 9's export includes a column literally named
`severity` (lowercase), distinct in spelling from `Status_on_day_collection_summary` (the field
name used throughout `R_plotting_publish.ipynb` and `Notebooks/DA_analysis-symptom_onset*.Rmd`).
No cell in any notebook read this round assigns or derives `myeloid.obs['severity']` — it is
already present on the `covid_total` object the notebooks load (`covid.h5ad`), which is **not
itself published in this GitHub repository** (the processed object lives on the authors'
infrastructure / the ArrayExpress deposit, not in this code repo). **Single remaining question**:
whether `severity` is simply a duplicate/alias of `Status_on_day_collection_summary` or an
already-collapsed 3-level encoding used specifically for the myeloid scoring figure — not
resolvable from the files committed to this repository.

## 3. Sample list per comparison (Initial-only vs. all visits) — not resolved by this repository

Cell 2 of the scoring notebook subsets by `full_clustering` only; no cell restricts by
`Collection_Day`/visit/resample flag before scoring or export. Whether the authors excluded
resample visits before computing Table 7 cannot be determined from this notebook — the filtering,
if any, would happen in the still-missing downstream script identified in §2. **Single remaining
question**: same missing script as §2 — not independently resolvable here.

## 4. The Newcastle/London-only site restriction — found, and does NOT transfer to the IFN scoring notebook

**File**: `myeloid/R_plotting_publish.ipynb`, cell 10 (sha256:
`3a9c5e17a0b3130f28a657711a33f0be871895d14e0366ed820d398d1a301110`):

```r
#Exclude Cambridge samples
prolif<-prolif[!prolif$Site == "Cambridge",]
prolif<-prolif[!prolif$Status_on_day_collection_summary == "Asymptomatic",]
```

This is the Methods p.15 site restriction Codex flagged — **confirmed**, and confirmed to apply
specifically to the myeloid **proportions/differential-abundance** figure (Fig. 2B: `celltype`
frequency by `diagnosis`, feeding `myeloid_proportions.pdf`), which also separately excludes
`Asymptomatic`. **It is a different cell, a different input CSV
(`myeloid_clustering_for_stats.csv`), and a different analysis from the IFN/TNF/JAK scoring
notebook** (`myeloid/myeloid_signature_scoring_publish.ipynb`, §1 above), whose cell 2 subsets by
cell type only and does not filter by `Site` or exclude `Asymptomatic`. **Confirmed finding, not
an assumption**: the site restriction does not transfer to the Fig. 2f/Table 7 IFN-score
comparison as committed in this repository — consistent with the review's own caution not to
assume it does, and this project should not apply a blanket Cambridge exclusion to the monocyte
IFN positive-control analysis on the strength of the Methods' differential-abundance sentence
alone.

## 5. Clinical-field provenance (`Days_from_onset`, `Collection_Day`, `Status_on_day_collection_summary`, `Worst_Clinical_Status`) — partially resolved

**File**: `Notebooks/DA_analysis-symptom_onset.Rmd` and `...-symptom_onset - Updated.Rmd` (both
read; the "Updated" version's relevant cells, lines 44–76 and 126–134, are representative of both).

- `covid.meta$OrderedSeverity <- ordered(covid.meta$Status_on_day_collection_summary, levels=c(...))`
  — confirms `Status_on_day_collection_summary` is the author's own in-code name for the six-level
  sampling-day category (`Healthy, Asymptomatic, Mild, Moderate, Severe, Critical`), matching
  Supplementary Table 2's `Catergory on collection day` column exactly in spelling and level set —
  corroborating (not newly discovering) what `stephenson_sdrf_identity_crosswalk.py`'s sibling work
  already assumed.
- `Days_from_onset` is plotted directly, one boxplot point per row, against `OrderedSeverity`
  (`aes(x=OrderedSeverity, y=Days_from_onset)`), filtered only by excluding `NA` and a few
  diagnosis categories — **no cell in either Rmd addresses repeat-visit rows sharing an identical
  `Days_from_onset` value**, and no cell recomputes or redefines `Days_from_onset` per visit. This
  neither confirms nor refutes this project's own repeat-visit duplicate-onset observation; it
  shows the authors' own downstream analysis uses the field as a flat per-row value without
  flagging the issue, which is consistent with (but does not prove) "the saved value is a baseline
  captured once per individual." **Single remaining question, unchanged from before this
  investigation**: the exact per-row generation rule for `Days_from_onset` is not written down
  anywhere in this public repository.
- `Collection_Day` and `Worst_Clinical_Status` do not appear as assigned/derived variables in any
  notebook or script read this round (repo-wide `gh api search/code` for `Worst_Clinical_Status`
  and `Collection_Day` returns 0 hits in this repository — both are metadata-table columns
  consumed as-is, never computed by committed code). **Single remaining question**: their exact
  generation rule, if one exists beyond manual clinical annotation, is not in this repository
  either; Supplementary Table 2 itself is the only available definition (clinical case-note
  abstraction at the time of data entry), and this project has no basis to look further within
  the author's code for it.

## What this investigation did NOT do

Did not run any notebook. Did not fetch `covid.h5ad`, `full_bal.h5ad`, or any other data file (none
are in this public repository; they are large processed objects on the authors' own
infrastructure/deposit, out of scope for a code-only investigation). Did not contact the
corresponding author. Did not pick the Healthy/Mild/Severe grouping definition by checking which
one would make this project's own eventual DE/GSEA result look better — the grouping remains
unresolved for this project's purposes, and `POSITIVE_CONTROL_EVIDENCE_v1.md` marks the IFN-down
candidate accordingly.

# Files you were looking for, and the corpus-overlap table (2026-09-15)

## Both "missing" files exist as project artifacts

**`canonical_donor_labels.csv` with `well_block`/`cv_group`** — saved, 721 rows × 13 columns.
Columns: `cohort, accession, donor_id, disease, event, event_definition, platelet_capable,
label_source, join_basis, channel_unit, notes, cv_group, well_block`. Both grouping columns are
fully populated (721/721), so callers can switch their grouping key straight away. `cv_group` is
the connected component of the "shares a GEM well" graph — the minimum unit that can move between
folds without splitting a well; `well_block` is the raw pool/batch/library a donor was loaded in.

**The corpus-overlap plan** — `plan_verify-scfm-pretraining-corpus-platelet_dd0417db.json` is in
the project (artifact `6cb4632d-1fe9-46d4-897c-2a0b5e4fdc62`). It is a *plan*, not a result, but
two of its deliverables were completed and are also on file:
`census_platelet_megakaryocyte_by_version.csv` (CELLxGENE census platelet/megakaryocyte fractions
across 10 census releases) and `deliverable1_evidence_quotes.csv` (43 verbatim corpus-provenance
quotes covering Geneformer, scGPT, UCE, scPRINT-2, STATE).

## New: cohort × model corpus overlap

`scfm_corpus_overlap.csv` (long format: cohort, events, model, overlap_status, basis) and
`benchmark_cohort_cellxgene_presence.csv` (the DOI-level CELLxGENE check behind it).

Method: matched each benchmark cohort's publication DOI against all 390 CELLxGENE collections,
then combined with each model's own corpus statement from `deliverable1_evidence_quotes.csv`.
Surname/title matching was too noisy to use — DOI match only.

| cohort | events | scGPT | UCE | scPRINT-2 | Geneformer | STATE |
|---|---|---|---|---|---|---|
| Ren (GSE158055) | 12 | controls only | present | present | unverified | unlikely |
| Stephenson (E-MTAB-10026) | 8 | controls only | present | present | unverified | unlikely |
| COMBAT | 16 | controls only | present | present | unverified | unlikely |
| Liu (GSE161918) | 4 | controls only | present | present | unverified | unlikely |
| Chongqing (GSE279452) | 22 | absent | absent | absent | unverified | unlikely |
| Cillo (GSE180578) | 12 | absent | absent | absent | unverified | unlikely |
| Wilk 2021 JEM (GSE174072) | 5 | absent | absent | absent | unverified | unlikely |
| Schulte-Schrepping | 3 | absent | absent | absent | unverified | unlikely |
| Qiu (GSE167363) | 2 | absent | absent | absent | unverified | unlikely |
| Bernardes | 2 | absent | absent | absent | unverified | unlikely |

Reading of each column:

- **scGPT** used CXG census 2023-05-15 but **filtered to samples without disease conditions**, so
  the diseased donors — every event in this benchmark — are excluded by construction. Only the
  healthy controls from those four collections could be in its corpus. That is the most favourable
  overlap profile of the three census-derived models.
- **UCE** (IMA, CXG census 2023-07-10) and **scPRINT-2** (cellxgene mid-2023, primary cells) report
  no disease filter, so diseased cells from the four CXG collections are eligible and should be
  treated as present.
- **Geneformer V1 is now resolved by evidence** — see the section below. V2 is not yet resolvable.
- **STATE** is a perturbation model (Tahoe-100M / Replogle-Nadig); its 167M-cell observational
  embedding corpus is not enumerated in the quotes on file, so "unlikely but unverified".

## Geneformer, resolved for V1

Supplementary Table 1 of Theodoris et al. 2023 (`41586_2023_6139_MOESM4_ESM.xlsx`) is the full
Genecorpus-30M composition: 561 datasets with a `Source` column carrying complete literature
citations, so membership is auditable by reference rather than inferred from dates.

- **Ren 2021 is in Genecorpus-30M** — 2 datasets, 1,402,746 cells passing QC, i.e. essentially the
  whole cohort. Geneformer V1 was pretrained on Ren's cells.
- **Wilk is a partial, donor-level case.** The corpus contains Wilk et al.'s 2020 atlas
  (GSE150728, 121,311 cells). Our benchmark uses the 2021 JEM deposit (GSE174072), a different
  accession from the same Stanford biobank, and the JEM paper states it includes patients profiled
  in the 2020 cohort — so donor-level overlap is likely. Whether the 5 fatal GSE174072 donors are
  among them is **not** established and is worth checking against the 2020 donor list, because
  donor-level overlap is the strongest form of leakage available to a model.
- **Absent from Genecorpus-30M**: Stephenson, COMBAT, Liu, Cillo, Chongqing, Schulte-Schrepping,
  Bernardes, Qiu. Also worth noting: Lee 2020 (GSE149689) *is* in the corpus — no events, but it
  is on the severity-strata list.

So for Geneformer V1 the split is **17 overlapping events (Ren 12 + Wilk 5) vs 69 naive** — a
different partition from UCE/scPRINT-2 (40 vs 46), which is exactly why the arm has to be assigned
per model rather than once for the whole benchmark.

**V2 (Genecorpus-104M) cannot be resolved yet.** Its dataset list is cited to Chen & Venkatesh et
al., *Nature Computational Science* 2026, which is not retrievable via Crossref or PubMed at time
of writing, and the HuggingFace repo ships only tokenized Arrow shards with no source manifest.
Since V2 explicitly draws on CELLxGENE, the working assumption should be that it behaves like
UCE/scPRINT-2 (all four CXG cohorts present) — but if you benchmark a V2 checkpoint, that
assumption needs stating as such, or the corpus paper chasing down when it appears.

## Why this is a design asset rather than a caveat

**46 of the 86 events sit in cohorts absent from every census-derived corpus** (Chongqing 22,
Cillo 12, Wilk 5, Schulte 3, Qiu 2, Bernardes 2); the other 40 are in the four CXG collections.
That splits the leave-one-cohort-out evaluation into a corpus-naive arm and a corpus-overlapping
arm, which turns the standard objection to scFM benchmarks into a measurement: if a model's
advantage is larger on held-out cohorts that were in its pretraining data than on cohorts that
were not, that gap *is* the leakage effect, reported rather than assumed. Two things follow for
the paper: report LOCO results split by arm as a main result, and resolve the Geneformer rows
before doing so, since Geneformer is the one model whose arm assignment is currently unknown.
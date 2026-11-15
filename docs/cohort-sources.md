# Cohort Sources (for demographic audit / Supp Table 1)

**Last updated:** 2026-08-25 (corrected sepsis-source mis-citation and stale 4-cohort
structural caveat against a direct HPCC crosstab; see
`../archive/superseded_2026-09-21/docs_legacy/cohort-annotation-verification.md`)
**Used by:** Manuscript Methods (cohort description) + Supp Table 1 +
`inspect_cohort_demographics.py`
**Companion docs:** `submission-checklist.md` for submission-state context;
`../README.md` and `current-status.md` for benchmark-design rationale.

---

## Scope of this reference

The 302006 BIB benchmark covers the **human-platelet COVID + Sepsis mortality
scenario** (single source-cohort structure):

| Scenario | Species | Cells (full) | Primary task subset | Cohort framework |
|---|---|---:|---|---|
| Mortality | Human platelet | 47,977 | **42,095 cells / 322 donors / COVID + Sepsis only** | 11 source studies, integrated from Qiu et al.~2024 IJMS (PMID 38892136); see "Human cohorts" below |

Primary in-hospital mortality task: COVID + Sepsis subset only (42,095 cells, 322 donors, 8,850 fatal / 33,245 survived), explicitly disease-matched to avoid off-disease confounding by Flu, Lung_infect, Non_covid, SLE, and Healthy Controls. The 2026-05-22 label-definition correction includes severe COVID as survived (`y=0`) — semantically equivalent to severe sepsis-survived — which changes the survived class count from 15,049 (severe excluded) to 32,436 (severe included) and is the basis for the corrected `tab:crossdisease` and `tab:adaptation_curve`.

---

## Human cohorts

Integration source: `scripts/data/integrate_platelet_data_harmony.R` from Qiu et al. 2024 IJMS. The `Data_NO` column in the canonical h5ad maps to these accessions.

| Accession | First author | Journal/Year | Disease | Notes for demographic audit |
|---|---|---|---|---|
| GSE150728 | Wilk et al. | Nat Med 2020 | COVID-19 | Severe COVID atlas; demographics in Supplementary |
| GSE149689 | Lee et al. | Sci Immunol 2020 | COVID-19 + influenza | IFN-response focus; PBMC |
| EGAS00001004571 | Schulte-Schrepping et al. | Cell 2020 | severe COVID | German cohort, controlled access |
| GSE155673 | Arunachalam et al. | Science 2020 | mild vs.~severe COVID | systems immunology |
| GSE158055 | Ren et al. | Cell 2021 | COVID atlas (largest) | Chinese cohort, ~196 patients |
| covid19cellatlas (multiple) | Bernardes / Hashimoto et al. | various | COVID-19 fatal | late juncture; megakaryocyte focus |
| E-MTAB-10026 | Stephenson et al. | Nat Med 2021 | COVID-19 multi-omics | UK cohort |
| GSE151263 | Jiang et al. | JCI Insight 2020, PMID 32554932 | Sepsis (integrated `Disease` label) | source study framed as pneumonia/early-ARDS-vs-not in a sepsis population; re-labeled `Sepsis` in this integration. `Data_NO=302004data18`: 7 donors, 476 cells, `covid_severity=severe` for all 7 (zero fatal) |
| Bernardes_2020 (FastGenomics) | Bernardes et al. | Immunity 2020 | COVID-19 longitudinal | megakaryocyte response |
| GSE163668 | Combes et al. | Nature 2021 | COVID-19 (verified 2026-08-25: `Data_NO=302005data02` is 100% `Disease=COVID-19`/`Lung_infect` in the canonical h5ad, no Sepsis label) | protective immune states |
| **GSE167363** | **Qiu et al.** | **J Leukoc Biol 2021, PMID 34558746** | **Sepsis** | **The primary sepsis cohort (`Data_NO=302005data03`): 10 donors (4 fatal + 6 severe) + 2 HC. Anchors the cross-disease transfer test. Corrects a prior mis-citation as "Reyes et al." — Reyes et al. is a different paper (Nat Med 2020, `reyes2020immune`), cited elsewhere in the manuscript as background, not as the sepsis source.** |
| GSE142016 | Mistry et al. | PNAS 2019 | SLE | neutrophil diversity; off-disease (excluded from primary task) |

**Structural caveat for benchmark design (corrected 2026-08-25 against a direct HPCC crosstab of the canonical h5ad — see `../archive/superseded_2026-09-21/docs_legacy/cohort-annotation-verification.md`).** The disease-matched 322-donor / 42,095-cell primary-task cohort spans **eleven** source-study `Data_NO` groups, not four: Wilk (7 donors), Lee (11), Schulte-Schrepping (23), Arunachalam (7), Ren (59), Liu (52), Stephenson (102), Jiang (7), Bernardes (28), Combes (16), and Qiu (10) — sums to 322. Of these, only **seven** have at least one fatal donor and can form a standalone leave-one-cohort-out fold (AUC is undefined without both outcome classes): Schulte-Schrepping, Ren, Liu, Stephenson, Bernardes, and Combes (COVID-19), plus Qiu (sepsis). The other four — Wilk, Lee, Arunachalam, and Jiang (32 donors, 6,090 cells total) — consist entirely of severe/moderate ICU-survivor cases (zero fatal each); they remain in the pooled 5-fold `StratifiedGroupKFold` training and evaluation but are absent from the LOBO/per-cohort-AUC tables and figures specifically. Jiang contributes 7 of the study's 17 total sepsis donors, so the sepsis arm is genuinely two source studies (Qiu + Jiang), even though only Qiu's cohort appears in the per-cohort breakdown.

---

## Demographic info to extract (human cohorts)

For each human cohort, look for these in the source paper (Methods, Supplementary Table 1, or GEO/EGA SOFT/MAGE-TAB metadata):

1. **Sex** — count of M/F (and "not reported")
2. **Age** — median + IQR or range; per-patient if available
3. **Ancestry / ethnicity** — usually inferred from country of recruitment if not explicitly reported. Mark "not reported" if absent.
4. **Recruitment site / country** — used as an ancestry proxy when explicit ancestry data is missing
5. **N patients per severity** — sanity check against the `Data_NO × covid_severity` crosstab

---

## Audit performed 2026-05-01 — finding: most studies did NOT provide demographic info

**Finding from manual review of source papers and GEO records:** the majority of contributing human cohorts do not report sex, age, or ancestry at the patient or cohort level. This is a structural feature of how public scRNA-seq COVID-19 cohorts were curated in 2020–2021 (a global pandemic with rushed publication timelines), not a quirk of our integration.

**Decision (preserved through 2026-05-24):**

1. Document the audit explicitly in the manuscript Limitations.
2. Cite this as evidence that prospective scRNA-seq studies need standardized demographic reporting (field-level recommendation, aligned with CZI Virtual Cells Workshop recommendations).
3. Skip Supp Table 1 demographic columns and replace with a simple cohort-summary table (accession + N patients + N cells + severity distribution), which is what reviewers will actually verify.

This audit itself is a contribution — turns a vague "limitations" handwave into a concrete documented finding that the field needs to address.

---

## Method for filling in the cohort summary table

For each row in the canonical h5ad:

1. Identify which source dataset the `Data_NO` value corresponds to using the tables above + the integration script. The R script `scripts/data/integrate_platelet_data_harmony.R` (from the IJMS paper) holds the mapping.
2. Look up the source paper Table 1 / Supplementary Table 1 / cohort description.
3. Fill in `N cells`, `N patients`, severity distribution, country/site.
4. Report ancestry if the paper reports it; otherwise write `not reported (recruitment site: X)`.

The completed table becomes `Supp_Table_1_cohorts.csv` for the manuscript. Cite it in the Methods cohort description and in the relevant TODO_v1_submission section.

# Deliverable 6 — direction-framing evidence for the Discussion

Every quote below was checked programmatically against the retrieved source text: 25 of 27 rows verify as exact (whitespace-normalised) substrings, the remaining two are cell values read out of supplementary workbooks. Evidence table: [deliverable6_evidence_quotes.csv]({{artifact:art_6f1c1866-bf8a-4f95-9cbd-a10abbc938b9}}); the 6b landscape is [deliverable6b_platelet_endpoint_landscape.csv]({{artifact:art_0b1a6814-28b6-4960-9bf6-b4625566125a}}).

---

## 6a. Is platelet-transcriptome outcome prediction established for sepsis — and in which window?

Two primary human studies make an outcome claim from platelet transcriptomes in sepsis, and both sample early.

**Qiu 2021** — the study you already have — draws at hours: "PBMC collected from septic patients in the intensive care unit (ICU) at 0 and 6 h post sepsis diagnosis" (n = 5 septic patients, within-subject fatal-vs-surviving design).

**Middleton 2019 (Blood)** is the only paper that states a predictive claim outright, in its Key Points: "Sepsis-induced upregulation of ITGA2B predicts mortality." — restated in the abstract as "ITGA2B is upregulated and functional in sepsis due to trafficking from megakaryocytes and de novo synthesis in platelets and is associated with increased mortality." **Its sampling window is not reported in anything I could access**: the Europe PMC record is `isOpenAccess = N` and holds only front matter, so the human cohort's draw timing and size could not be verified. Flagged as unverified.

**The same group's 2021 platelet-sepsis cohort (n = 59 septic patients vs 46 donors) states its window explicitly** — "All septic patients were enrolled within 72 h of intensive care unit (ICU) admission." — and, importantly for us, the mortality signal there was *not* transcriptional: "While both the pro- and active forms of HPSE protein increased in platelets during sepsis, only the active form of HPSE protein significantly correlated with sepsis-associated mortality."

Everything else that looks like sepsis-outcome prediction from platelets is not platelet transcriptomics. The two 2024 "platelet-related gene" prognostic models are built on public bulk cohorts — "Sepsis patient data from the Gene Expression Omnibus (GEO) database and PRGs from the Molecular Signatures Database (MSigDB) were analyzed." — with no platelet isolation and no stated draw window. The only other prospective platelet-transcriptome sepsis study is veterinary and has no outcome endpoint ("6 client-owned dogs with sepsis and 6 healthy dogs were enrolled from March through August 2023").

![Published platelet-transcriptome outcome claims in sepsis all come from draws inside 72 h]({{artifact:art_6496401b-aed7-4512-92bb-5447fd4b0a98}})

> **Verdict 6a:** Outcome-predictive platelet-transcriptome signal in sepsis is documented only in early windows — 0–6 h post diagnosis, or enrolment within 72 h of ICU admission — and no human study reports platelet-transcriptome outcome prediction from draws beyond the first days, so our null at a median of 16 days post admission falls in untested territory rather than contradicting an established positive.

---

## 6b. Which clinical endpoints ARE established from platelet transcriptomics?

Full table with cohort sizes in `deliverable6b_platelet_endpoint_landscape.csv`. Ranked by strength of evidence:

1. **Cancer detection (tumour-educated platelets, thromboSeq)** — the field's signal-positive anchor and the only line with large validation cohorts. "We distinguished 228 patients with localized and metastasized tumors from 55 healthy individuals with 96% accuracy" (2015); "accurate TEP-based detection of early- and late-stage non-small-cell lung cancer (n = 518 late-stage validation cohort, accuracy, 88%…" (2017); the 2022 follow-up detected cancer in two-thirds of 1,096 stage I–IV samples at 99% specificity in asymptomatic controls, but specificity fell to 78% in symptomatic controls.
2. **SLE disease activity / vascular phenotype** — n = 89 (71 SLE, 18 controls) with repeat draws: "We outline a platelet transcriptomic signature for SLE activity that is significantly associated with measures of platelet activity and vascular health, implying it may serve as an indicator of cardiovascular risk." Association, explicitly not a validated predictor.
3. **Aplastic-anaemia treatment response / long-term prognosis** — a positive claim (J Thromb Haemost 2025) that I could not verify beyond the abstract; no full text available. **Flagged UNVERIFIED.**
4. **COVID-19 care level and platelet function** — n = 25 COVID patients (ICU and non-ICU) + 5 donors, drawn "within 72 hours of hospital or intensive care unit (ICU) admission"; the endpoint is function by care level ("Critically ill ICU COVID-19 patients had increased aggregation in response to low-dose agonist compared with non-ICU patients"), not mortality, and no classifier or AUC is reported. The 10-patient platelet RNA-seq study adds mechanism, not prediction.
5. **Sepsis mortality** — as in 6a: one single-transcript claim and one protein-level association.
6. **Descriptive only** — myocardial infarction in women (n = 40 vs 38; "do not significantly change over time"), immune thrombocytopenia (2026, unverified).
7. **Dengue** — nothing qualifying. A 69-hit search for platelet + dengue + transcriptomics returned only megakaryoblast/cell-line infection work and whole-blood severity studies. Do not cite dengue as a platelet-transcriptomics anchor.

> **Verdict 6b:** The only endpoint with large, replicated, prediction-grade evidence from platelet RNA is cancer detection; COVID severity, SLE activity and sepsis mortality exist as associations in cohorts of 5–89 patients without validated patient-level classifiers — so a benchmarking paper should treat cancer detection (and, at a stretch, COVID severity) as its signal-positive anchor and cannot be scooped on a negative sepsis-mortality result, because no positive one has been published.

---

## 6c. Does the Chongqing source paper report platelet/megakaryocyte findings?

**Yes for annotation, no for any traceable finding.** Their atlas does contain a platelet-lineage population: Supplementary Table 6 lists a `Megakaryocyte` cluster among 47 annotated populations, with 956 marker rows headed by PPBP (avg_log2FC 8.11; pct.1 0.991 vs pct.2 0.082), NRGN, PF4, MYL9, GP9 — and ITGA2B and GP1BA mark that cluster and no other. So megakaryocyte/platelet transcripts are present in the same PBMC libraries we embed.

But no platelet- or MK-centred result is traceable in any accessible material. The abundance-statistics workbooks behind Figs. 2b/3d/4d/5e contain no megakaryocyte or platelet row at all (Fig. 2b's tests cover the seven CD4 T-cell clusters); the 55-page peer review file contains zero occurrences of "platelet" or "megakaryocyt"; and in the clinical tables "Platelet Count" appears only as a laboratory variable. **Caveat:** the article is closed access with no PMCID, so I cannot verify what the body text says, and one supplementary PDF is image-only (0 extractable characters per page) and was therefore not searchable.

> **Verdict 6c:** The source study annotates a megakaryocyte population defined by canonical platelet genes (PPBP, PF4, GP9, ITGA2B) but reports no platelet- or megakaryocyte-focused analysis in any accessible supplementary material, so our null is specifically about 28-day mortality prediction from these cells — not a claim that platelet-lineage biology is absent or unstudied in this cohort.

---

## 6d. What do scFM evaluations report for donor-level clinical outcomes?

The closest evaluation is **TissueFormer (BMC Bioinformatics 2026)**, and it is close indeed: it predicts COVID-19 severity from blood scRNA-seq in three cohorts totalling 376 donors — including the Ren and Stephenson cohorts already in our benchmark — under the protocol we use: "Each dataset was evaluated independently using 5-fold stratified group k-fold cross-validation with all cells from a donor assigned to the same fold". Its margin over simple baselines is real but cohort-size dependent: "at N=64 TissueFormer achieved 70% balanced accuracy compared to 54% for the best classical baseline (logistic regression on cell type histograms) and 52% for random forest on pseudobulk", whereas "On the smaller COMBAT dataset, TissueFormer and logistic regression on pseudobulk were closer in performance (87% vs. 85% at N=256)". Two things to note before citing it as an FM result: its per-cell encoder is Geneformer-architecture, but the **foundation-model comparison is in its mouse brain-region task**, not the human COVID task, where the comparators are pseudobulk LR/RF, cell-type histograms, CellCnn, scAGG and ScRAT; and ScRAT matched it.

Where pretrained embeddings *have* been tested against plain expression for a patient-level clinical endpoint, they did not win. A 2026 survival-prediction paper derived patient-level scFoundation embeddings for 25 TCGA cancer types and reports of the prior attempts: "trained Ridge-regularized Cox proportional hazards models (Ridge–Cox) using embeddings from CancerFoundation and scGPT, yet these embeddings did not consistently outperform models trained directly on gene expression." Note this is bulk transcriptome input, not patient scRNA-seq cohorts.

At cell level the negative literature is now explicit — "our findings demonstrate that scGPT and Geneformer in zero-shot configurations perform inconsistently compared to cell embeddings derived from HVG or generated using scVI or Harmony" (Genome Biology 2025) — and a 2026 preprint reports that "by carefully choosing pre-processing steps and normalization procedures, it is possible to achieve SOTA or near-SOTA performance using simple representations where cells are viewed as vectors in gene expression space."

Meanwhile purpose-built supervised aggregators do reach high patient-level accuracy without any FM: scPhase reports "a mean AUC of 0.895 on receiver operating characteristic (ROC) curves" for COVID-19 severity, and scMILD reports sample-level AUROCs on COVID and lupus datasets — neither benchmarks against FM embeddings.

> **Verdict 6d:** Donor-level clinical prediction from scRNA-seq is demonstrably achievable — 70–87% balanced accuracy for COVID severity under donor-grouped cross-validation — but in every evaluation I could retrieve, the winning models are supervised multi-cell architectures rather than pretrained foundation-model embeddings, and the two studies that tested FM embeddings against expression or pseudobulk baselines found no consistent advantage; so "consistent with the emerging negative literature on foundation-model embeddings" is defensible, and I found no prior patient-level scFM evaluation on sepsis mortality at all.

---

## Claims I could not verify

- **Middleton 2019's human sampling window and cohort size** — body paywalled (Blood; Europe PMC front matter only). The mortality claim itself is quoted from the paper's own Key Points and abstract.
- **The Chongqing paper's main-text treatment of megakaryocytes** — closed access, no PMCID; one supplementary PDF is image-only and unsearchable.
- **Aplastic-anaemia platelet-transcriptomics prognosis (J Thromb Haemost 2025)** and **ITP platelet transcriptomics (2026)** — abstract only, no full text retrieved; cohort sizes and performance not checked.
- **Absence claims** are search-bounded, not proofs: 6a rests on 95 + 45 PubMed hits for platelet/sepsis/transcriptomics screened to 27 with platelet or megakaryocyte in the title; the dengue negative rests on 69 hits; 6d on 101 + 20 + 24 hits across three query formulations. A study that names neither platelets nor its endpoint in title/abstract could have been missed.

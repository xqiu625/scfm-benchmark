# Evidence Note — Platelet Representation in Foundation Model Pretraining Corpora

**Status:** active evidence note. Reworked 2026-09-21 from
`docs/memory-bank/platelet-pretraining-representation.md` (2026-05-27), which is now archived
under `archive/superseded_2026-09-21/`.

**What changed in this rework.** The archived version asserted, without per-model verification,
that platelets are absent from all five models' pretraining corpora, and asserted that droplet
scRNA-seq does not capture platelets. The second claim is contradicted by this project's own
data. Both are restated below with their actual evidential strength, and the benchmark numbers
quoted in the old "Implications" section have been removed rather than carried forward, because
they predate the 2026-09 QC rebuild.

---

## 1. The claim, graded by evidence strength

The headline claim is that platelets are out-of-distribution for all five foundation models.
That claim is *probably* right, but the evidence behind it is not uniform across models, and the
archived note presented one inference chain as if it were five independent verifications.

| Model | Primary corpus | Evidence for "no platelets" | Strength |
|---|---|---|---|
| scGPT | CELLxGENE (33M) | CELLxGENE CellGuide reports no collections for CL:0000233 | **Corpus-level, indirect** |
| UCE | CELLxGENE (33.9M) + 28 datasets | Same CellGuide inference; the 2.3M non-Census cells were not checked | **Partial** |
| scPRINT-2 | CELLxGENE Census + lamin DB | Same CellGuide inference; lamin DB portion not checked | **Partial** |
| STATE | CELLxGENE + scBaseCount + Tahoe-100M | CellGuide for the Census portion; scBaseCount (~230M cells, 72 tissues) **not** checked | **Weak for the scBaseCount portion** |
| Geneformer | Genecorpus-30M / V2 (561 GEO+HCA datasets) | Not CELLxGENE-derived at all; rested on the "droplet-based only" curation rule, which section 2 shows does not imply platelet-free | **None — record as `unknown / not established`, not `~0`** |

**Three distinct statements must be kept apart**, and the archived note conflated them:

1. *explicitly excluded* by a documented curation rule — not demonstrated for any of the five;
2. *unlabeled in the cell ontology* — this is what the CellGuide result actually shows;
3. *possibly present but unannotated* — not excluded by any evidence we hold.

The CellGuide statement ("Neither this cell type nor its ancestors and descendants are available
in any collections") is evidence for (2). Treating it as evidence for (1) is the inferential gap.

**Not yet done:** the verification snippet in section 5 has never been run, and no Census release
or version was recorded in the archived note. Until it is run against a pinned Census version,
the per-model "~0" figures should be cited as *expected*, not *measured*.

---

## 2. Correction: droplet scRNA-seq does capture platelets

The archived note stated that "standard droplet-based scRNA-seq protocols (10x Chromium,
Drop-seq) capture nucleated cells; platelets are systematically excluded." **This project's own
data contradicts that.** Every cohort here is standard droplet or Seq-Well PBMC sequencing, and
every one of them yields platelets:

| Cohort | Platform | Platelets present | Artifact |
|---|---|---|---|
| Qiu | 10x | 5,276 | `gate1_qc/qiu/` lineage audit |
| Stephenson | 10x | 16,911 | `gate1_qc/stephenson/` lineage audit |
| Wilk | Seq-Well | 10,502 | Seq-Well-rule retained set |
| Cillo | 10x | 692 | `gate1_qc/cillo/` lineage audit |

The accurate mechanism is different, and more useful to this project: **platelets are captured,
then removed by conventional QC.** Platelets carry very little RNA, so a generic `total_counts
>= 500` / `n_genes >= 200` filter deletes them preferentially. Measured on Wilk, the generic
thresholds removed **76.7% of platelets** versus 57.5% of cells overall; under the frozen
platform-calibrated Seq-Well rules that fell to **2.46%** (`qc_checklist.md` §6e).

So the correct framing is not "the assay cannot see platelets" but "the standard QC pipeline
discards them." That is a statement about analysis convention, not about droplet chemistry.

**This does not make the OOD claim stronger.** It supplies a *plausible mechanism* for
underrepresentation — if a corpus was built with conventional per-cell thresholds, platelets
would be removed at that step — but a mechanism is not a measurement. Actual representation
still has to be established separately for each model. Wording that can be defended:

> Platelets are captured by droplet scRNA-seq but their low RNA content makes them vulnerable to
> conventional complexity thresholds. This provides a plausible mechanism for underrepresentation
> in curated pretraining corpora; actual representation must be established separately for each
> model.

Note also that this mechanism *removes* the support that section 1 recorded for Geneformer.
Geneformer's "~0" rested entirely on its droplet-only curation rule, and droplet-only does not
imply platelet-free. **Geneformer should be listed as `unknown / not established`, not `~0`,**
until its Genecorpus is audited directly.

**Still true from the archived note** (low-level biology, independently sourced):
platelets are anucleate megakaryocyte fragments carrying ~2.2 fg RNA versus ~10 pg for a
nucleated cell (~5,000× less), and under 10% of circulating platelets contain detectable RNA.

---

## 3. Megakaryocyte signal (unchanged, still unverified)

Megakaryocytes are nucleated and are plausibly present in pretraining corpora. Shared
transcriptional programs could give indirect signal, moderated by the RNA remodeling platelets
undergo during thrombopoiesis. This was speculation in the archived note and remains speculation;
it has not been tested here.

---

## 4. What this supports in the manuscript

Framing only: the benchmark tests generalization to a clinically relevant cell type that is
underrepresented in pretraining corpora. It does **not**, on current evidence, support a claim
that any specific model's corpus provably contains zero platelets.

Benchmark AUC numbers that the archived note used to argue this point have been deliberately
omitted here. They predate the QC rebuild and must be re-derived from the rebuilt input; quoting
them in an evidence note is how superseded numbers survive. Current numbers belong in
`docs/current-status.md`, bound to artifacts.

---

## 5. Verification to run before citing per-model figures

Pin and record the Census release before quoting any number this produces.

```python
import cellxgene_census

CENSUS_VERSION = "..."  # record the exact release; results are version-dependent
census = cellxgene_census.open_soma(census_version=CENSUS_VERSION)
for ct in ("platelet", "megakaryocyte"):
    obs = cellxgene_census.get_obs(
        census, "homo_sapiens",
        value_filter=f"cell_type == '{ct}' and is_primary_data == True",
        column_names=["cell_type"])
    print(CENSUS_VERSION, ct, len(obs))
```

This covers only the CELLxGENE-derived portion. Geneformer's Genecorpus and STATE's scBaseCount
need separate treatment, and there is currently no plan on file for either.

---

## References

- CZ CELLxGENE CellGuide — Platelet: https://cellxgene.cziscience.com/cellguide/CL_0000233
- Viggiani et al. (2025). Single-cell RNA sequencing of platelets: challenges and potential. *J Thromb Thrombolysis*. https://doi.org/10.1007/s11239-025-03153-8
- Wolfsberger et al. (2024). The First Comprehensive Description of the Platelet Single Cell Transcriptome. https://pubmed.ncbi.nlm.nih.gov/39464035/
- Cui et al. (2024). scGPT. *Nature Methods*. https://doi.org/10.1038/s41592-024-02201-0
- Theodoris et al. (2023). Transfer learning enables predictions in network biology. *Nature* 618:616-624. https://doi.org/10.1038/s41586-023-06139-9
- Rosen et al. (2024). Universal Cell Embeddings. *Nature Methods*. https://github.com/snap-stanford/UCE
- Arc Institute STATE: https://arcinstitute.org/news/virtual-cell-model-state
- Genecorpus-30M: https://huggingface.co/datasets/ctheodoris/Genecorpus-30M

# Platelet Representation in Foundation Model Pretraining Data

**Last Updated:** May 27, 2026

---

## Key Finding

**Platelets are absent or negligibly represented in the pretraining data of all five primary foundation models evaluated in this project (STATE, UCE, Geneformer, scPRINT-2, scGPT).** This makes our 47,808-platelet COVID + sepsis mortality benchmark an out-of-distribution generalization test for every model. TranscriptFormer (legacy) and STACK (removed 2026-05-26) shared the same property.

---

## Why Platelets Are Missing from scRNA-seq Corpora

Platelets are **anucleate cell fragments** (no nucleus) derived from megakaryocytes. They have extremely low RNA content:

- **~2.2 fg RNA per platelet** vs **~10 pg per nucleated cell** (~5,000x less)
- Less than 10% of circulating platelets contain detectable RNA
- Standard droplet-based scRNA-seq protocols (10x Chromium, Drop-seq) capture nucleated cells; platelets are systematically excluded
- Specialized isolation protocols (e.g., acridine orange staining + FACS sorting) are required for platelet scRNA-seq

The [CZ CELLxGENE CellGuide page for platelets (CL:0000233)](https://cellxgene.cziscience.com/cellguide/CL_0000233) explicitly states:

> "Neither this cell type nor its ancestors and descendants are available in any collections."

---

## Per-Model Analysis

### STATE (Arc Institute)

| Attribute | Details |
|-----------|---------|
| **Pretraining cells** | ~270M (167M observational + 100M+ perturbational) |
| **Observational sources** | CELLxGENE (~107M cells) + scBaseCount (~230M cells) |
| **Perturbational sources** | Tahoe-100M (50 cancer cell lines), Parse-PBMC, Replogle-Nadig |
| **Platelet representation** | **~0**. CELLxGENE has no platelets. Tahoe-100M is cancer cell lines. scBaseCount spans 72 tissues but platelets are not captured by standard scRNA-seq. |

### scGPT (Wang Lab, University of Toronto)

| Attribute | Details |
|-----------|---------|
| **Pretraining cells** | 33M normal human cells |
| **Source** | CELLxGENE, 51 organs, 441 studies |
| **Variants** | Whole-human (33M), Blood (10.3M), Brain (13.2M), Kidney (814K) |
| **Platelet representation** | **~0**. Sourced from CELLxGENE which has no platelet entries. The blood variant (10.3M cells) likely contains nucleated blood cells only (PBMCs, not platelets). |

### Geneformer (Broad Institute / Gladstone)

| Attribute | Details |
|-----------|---------|
| **Pretraining cells** | 29.9M (V1, Genecorpus-30M) / 104M (V2) |
| **Sources** | 561 datasets from NCBI GEO, Human Cell Atlas, EMBL-EBI, etc. |
| **Curation** | Droplet-based platforms only; excluded malignant cells and immortalized lines; no tissue >25% of data |
| **Platelet representation** | **~0**. Restricted to droplet-based sequencing, which does not capture platelets. |

### UCE (Stanford, Leskovec Lab)

| Attribute | Details |
|-----------|---------|
| **Pretraining cells** | 36M cells (Integrated Mega-scale Atlas) |
| **Sources** | 33.9M from CELLxGENE + 2.3M from 28 additional datasets, 8 species |
| **Platelet representation** | **~0**. Primary source is CELLxGENE which has no platelets. |

### scPRINT-2 (Kalfon et al., Lamin / scPRINT consortium)

| Attribute | Details |
|-----------|---------|
| **Pretraining cells** | ~50M cells, ontology-aware via lamin (NCBITaxon, Cell Ontology) |
| **Sources** | CELLxGENE Census + curated lamin DB |
| **Platelet representation** | **~0**. Primary source is CELLxGENE which has no platelets. |

### TranscriptFormer (Chan Zuckerberg Initiative) — legacy, not in 5-FM primary set

| Attribute | Details |
|-----------|---------|
| **Pretraining cells** | 57M human (TF-Sapiens) / 112M across 12 species (TF-Metazoa) |
| **Sources** | CELLxGENE (644 human datasets, 88 mouse datasets), Tabula Sapiens, ZebraHub, Bgee |
| **Platelet representation** | **~0**. Primary source is CELLxGENE which has no platelets. |

---

## Summary Table (5-FM primary set + legacy)

| Model | Pretraining Scale | Primary Data Source | Platelet Cells | OOD for Platelets? |
|-------|:-----------------:|--------------------:|:--------------:|:-------------------:|
| **STATE** | 270M | CELLxGENE + scBaseCount + Tahoe-100M | ~0 | Yes |
| **UCE** | 36M | CELLxGENE + 28 datasets | ~0 | Yes |
| **Geneformer** | 30M (V1) / 104M (V2) | Genecorpus-30M (561 datasets) | ~0 | Yes |
| **scPRINT-2** | ~50M | CELLxGENE Census + lamin DB | ~0 | Yes |
| **scGPT** | 33M | CELLxGENE | ~0 | Yes |
| TranscriptFormer (legacy) | 57–112M | CELLxGENE + atlases | ~0 | Yes |

---

## Implications for the Benchmark

1. **Out-of-distribution test**: Our 47,808 platelet dataset is entirely OOD for all five primary foundation models. None were pretrained on platelets.

2. **Stronger generalization claim**: Performance on an unseen cell type is a more rigorous test of foundation model utility than evaluation on well-represented cell types (e.g., PBMCs, T cells).

3. **Explains baseline competitiveness**: PCA-200+LR on raw expression achieves 0.705 mortality AUC (patient-level CV), ties Geneformer (0.738) and STATE (0.724) within fold noise, and beats UCE (0.598). On the cross-disease adaptation curve PCA-200+LR asymptotes to 0.964 — beats every FM. This is consistent with the models having no platelet-specific representations to leverage; the OOD setting makes the FM advantage *conditional, not automatic*.

4. **Per-FM strengths differ on different axes**: STATE is the only bidirectional cross-disease generalizer (forward 0.798, reverse 0.707, severe-excluded); scPRINT-2 has the best within-disease zero-shot AUC (0.762); Geneformer + LoRA gives +0.006 marginal AUC over zero-shot under patient-level CV; scGPT exhibits an anti-correlated cross-disease forward transfer (0.352).

5. **Megakaryocyte signal**: Megakaryocytes (platelet precursors, nucleated) are likely present in pretraining data. Shared transcriptional programs between megakaryocytes and platelets could provide indirect signal, though platelets undergo significant RNA remodeling during thrombopoiesis.

6. **Manuscript framing**: This finding supports framing the benchmark as testing foundation model generalization to a clinically relevant but underrepresented cell type, rather than performance on a well-characterized cell population.

---

## Verification

To confirm platelet absence in CELLxGENE programmatically:

```python
import cellxgene_census

census = cellxgene_census.open_soma()
platelet_obs = cellxgene_census.get_obs(
    census, "homo_sapiens",
    value_filter="cell_type == 'platelet' and is_primary_data == True",
    column_names=["cell_type"]
)
print(f"Platelet cells in CELLxGENE: {len(platelet_obs)}")

# Also check megakaryocytes for comparison
mk_obs = cellxgene_census.get_obs(
    census, "homo_sapiens",
    value_filter="cell_type == 'megakaryocyte' and is_primary_data == True",
    column_names=["cell_type"]
)
print(f"Megakaryocyte cells in CELLxGENE: {len(mk_obs)}")
```

---

## References

- CZ CELLxGENE CellGuide — Platelet: https://cellxgene.cziscience.com/cellguide/CL_0000233
- Viggiani et al. (2025). Single-cell RNA sequencing of platelets: challenges and potential. *J Thromb Thrombolysis*. https://doi.org/10.1007/s11239-025-03153-8
- Wolfsberger et al. (2024). The First Comprehensive Description of the Platelet Single Cell Transcriptome. https://pubmed.ncbi.nlm.nih.gov/39464035/
- Cui et al. (2024). scGPT: toward building a foundation model for single-cell multi-omics. *Nature Methods*. https://doi.org/10.1038/s41592-024-02201-0
- Theodoris et al. (2023). Transfer learning enables predictions in network biology. *Nature* 618: 616-624. https://doi.org/10.1038/s41586-023-06139-9
- Rosen et al. (2024). Universal Cell Embeddings. *Nature Methods*. https://github.com/snap-stanford/UCE
- Pearce et al. (2025). TranscriptFormer. *bioRxiv*. https://doi.org/10.1101/2025.04.25.650731
- Arc Institute STATE announcement: https://arcinstitute.org/news/virtual-cell-model-state
- Genecorpus-30M: https://huggingface.co/datasets/ctheodoris/Genecorpus-30M

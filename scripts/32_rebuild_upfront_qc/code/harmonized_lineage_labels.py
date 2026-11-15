"""
Harmonized, marker-score-based lineage labeling used identically across every
cohort and every GEM-well channel as the z argument to decontX(). Never a
cohort's own published annotation -- see project_302006_rebuild_architecture
memory for why (per-cohort z makes ambient-removal magnitude a function of
each paper's annotation granularity on the same underlying soup).

Fixed at 8 categories: Platelet, T_cell, B_cell, NK_cell, Monocyte_Myeloid,
DC, HSPC, Erythroid. Deliberately excludes a "Pan_leukocyte" bucket (PTPRC
alone is not a lineage -- it is high in every nucleated category, so scoring
it as a competing cluster would be wrong, not just uninformative).

Platelet and non-platelet-lineage panels below are the exact panels already
used in scripts/31_donor_level_revision/code/r6_purity_and_gene_audit.py and
pca_baseline_restricted_patient_oof.py -- reused unchanged for continuity
between the retrospective audit and this upfront pipeline. HSPC is new here
(no prior panel existed; the audit never needed a progenitor bucket).
"""
import numpy as np
import scanpy as sc

LINEAGE_MARKERS = {
    "Platelet": ["PPBP", "PF4", "ITGA2B", "ITGB3", "GP9", "GP1BA", "TUBB1", "NRGN", "CLU", "SPARC"],
    "T_cell": ["CD3D", "CD3E", "CD3G", "CD2", "TRBC1", "TRBC2", "IL7R", "GZMA", "GZMB", "GZMK",
               "GZMH", "CD8A", "CD8B", "CD4", "CCL5", "IL32", "CD7"],
    "B_cell": ["CD79A", "CD79B", "MS4A1", "CD19", "IGHM", "IGHD", "CD74"],
    "NK_cell": ["GNLY", "NKG7", "KLRD1", "KLRB1", "PRF1", "NCAM1", "KLRF1"],
    "Monocyte_Myeloid": ["CD14", "LYZ", "FCN1", "S100A8", "S100A9", "S100A12", "FCER1G", "TYROBP",
                          "AIF1", "MNDA", "C1QA", "C1QB", "C1QC", "AXL", "VCAN"],
    "DC": ["FCER1A", "CLEC10A", "ITGAX", "CD1C"],
    "HSPC": ["CD34", "KIT", "GATA2", "SPINK2", "MPO", "ELANE"],
    "Erythroid": ["HBB", "HBA1", "HBA2", "SLC25A37", "ALAS2"],
}

MIN_SCORE_MARGIN = 0.0  # top score must exceed the runner-up by more than this to avoid "Unclear"
UNCLEAR_LABEL = "Unclear"


def label_channel(adata, layer=None, min_genes_present=3, random_state=0):
    """
    Assign one lineage label per cell using sc.tl.score_genes for each of the
    8 fixed panels, then argmax. Designed to run on ONE channel's filtered
    count matrix at a time (a few hundred to a few thousand cells) -- no
    fitting is shared across channels, so this is safe to call independently
    and in parallel for every channel in every cohort.

    Returns a pandas Series of labels indexed like adata.obs_names.
    """
    var_upper = {g.upper(): g for g in adata.var_names.astype(str)}
    scores = {}
    for lineage, genes in LINEAGE_MARKERS.items():
        present = [var_upper[g] for g in genes if g.upper() in var_upper]
        if len(present) < min(min_genes_present, len(genes)):
            scores[lineage] = np.full(adata.n_obs, -np.inf)
            continue
        key = f"_score_{lineage}"
        sc.tl.score_genes(adata, gene_list=present, score_name=key,
                           layer=layer, random_state=random_state, use_raw=False)
        scores[lineage] = adata.obs[key].to_numpy()
        del adata.obs[key]

    score_mat = np.column_stack([scores[k] for k in LINEAGE_MARKERS])
    order = np.argsort(-score_mat, axis=1)
    top_idx = order[:, 0]
    top_val = score_mat[np.arange(len(top_idx)), top_idx]
    second_val = score_mat[np.arange(len(top_idx)), order[:, 1]]

    names = list(LINEAGE_MARKERS.keys())
    labels = np.array([names[i] for i in top_idx], dtype=object)
    labels[(top_val - second_val) <= MIN_SCORE_MARGIN] = UNCLEAR_LABEL
    labels[~np.isfinite(top_val)] = UNCLEAR_LABEL

    import pandas as pd
    return pd.Series(labels, index=adata.obs_names, name="harmonized_lineage")


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description="Label one channel's filtered matrix with harmonized lineage z-labels")
    ap.add_argument("--h5ad", required=True, help="Filtered, cell-called count matrix for ONE channel")
    ap.add_argument("--out-csv", required=True)
    ap.add_argument("--channel-id", required=True, help="GEM-well/channel identifier, stamped into the output for the batch= argument downstream")
    args = ap.parse_args()

    a = sc.read_h5ad(args.h5ad)
    labels = label_channel(a)
    df = labels.to_frame()
    df["channel_id"] = args.channel_id
    df.to_csv(args.out_csv)
    print(f"{args.h5ad}: {a.n_obs} cells, label counts:")
    print(labels.value_counts())

"""
Prerank GSEA on the Geneformer gradient-attribution ranked gene list, against
the three pre-downloaded Enrichr-style genesets (Hallmark 2020, KEGG 2021,
Reactome 2022). Model-agnostic: reads MODEL_NAME's ranked CSV, writes GSEA
results under the same per-model output dir. Rerun with a different
MODEL_NAME once the other 4 models' attribution scripts exist.
"""
import os
import sys
import json
import pandas as pd
import gseapy as gp

MODEL_NAME = sys.argv[1] if len(sys.argv) > 1 else "geneformer"
BASE = "/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/gene_attribution_gsea"
RANKED_CSV = f"{BASE}/{MODEL_NAME}/{MODEL_NAME}_gene_attribution_ranked.csv"
GENESETS = {
    "Hallmark": f"{BASE}/genesets/MSigDB_Hallmark_2020.json",
    "KEGG": f"{BASE}/genesets/KEGG_2021_Human.json",
    "Reactome": f"{BASE}/genesets/Reactome_2022.json",
}
OUT_DIR = f"{BASE}/{MODEL_NAME}"

df = pd.read_csv(RANKED_CSV)
print(f"Loaded ranked list: {len(df)} genes", flush=True)
rnk = df[["gene", "mean_attribution"]].sort_values("mean_attribution", ascending=False)
rnk = rnk.drop_duplicates(subset="gene", keep="first")
print(f"After dedup: {len(rnk)} genes, range [{rnk.mean_attribution.min():.4f}, "
      f"{rnk.mean_attribution.max():.4f}]", flush=True)

for gs_name, gs_path in GENESETS.items():
    print(f"\n=== {gs_name} ===", flush=True)
    with open(gs_path) as f:
        gene_sets = json.load(f)
    res = gp.prerank(
        rnk=rnk,
        gene_sets=gene_sets,
        min_size=10,
        max_size=500,
        permutation_num=1000,
        outdir=None,
        seed=20260910,
        threads=4,
    )
    out = res.res2d.copy()
    out = out.sort_values("NOM p-val")
    out_path = f"{OUT_DIR}/{MODEL_NAME}_gsea_{gs_name.lower()}.csv"
    out.to_csv(out_path, index=False)
    print(f"Saved {out_path} ({len(out)} terms)", flush=True)
    sig = out[out["FDR q-val"] < 0.25]
    print(f"Terms with FDR q < 0.25: {len(sig)}", flush=True)
    cols = ["Term", "NES", "NOM p-val", "FDR q-val"]
    print(out[cols].head(15).to_string(index=False), flush=True)

print("\nDone.", flush=True)

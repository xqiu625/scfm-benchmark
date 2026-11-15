import sys
import pandas as pd
import pickle
import h5py

print("starting", flush=True)
sym2ens = pd.read_csv(
    "/bigdata/godziklab/shared/Xinru/302006/02_EMBEDDINGS/option_c_monocyte/gene_symbol_to_ensembl.csv"
).set_index("gene_symbol")["ensembl_id"].to_dict()
print(f"loaded sym2ens, n={len(sym2ens)}", flush=True)

with open("/bigdata/godziklab/shared/Xinru/302006/03_MODELS/scPEFT/geneformer_peft/geneformer/token_dictionary.pkl", "rb") as f:
    token_dict = pickle.load(f)
print(f"loaded token_dict, n={len(token_dict)}", flush=True)

with h5py.File("/bigdata/godziklab/shared/Xinru/302006/01_DATA/raw/human_platelet_covid_severity_raw_counts.h5ad", "r") as f:
    var_names = [x.decode() if isinstance(x, bytes) else x for x in f["var"]["_index"][:]]
print(f"loaded var_names, n={len(var_names)}", flush=True)

resolved = sum(1 for g in var_names if g in sym2ens)
mapped = sum(1 for g in var_names if sym2ens.get(g) in token_dict)
print(f"n_genes={len(var_names)}, resolved_to_ensembl={resolved}, mapped_to_token={mapped}")

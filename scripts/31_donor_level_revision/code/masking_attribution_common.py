"""
Shared cohort loading, deterministic top-K candidate-gene selection, and
output/aggregation logic for masking-based ("leave-one-gene-out") gene
attribution -- imported identically by geneformer_masking_attribution.py
(scpeft_env) and state_masking_attribution.py (state_env) so the
perturbation definition and downstream gene ranking are comparable by
CONSTRUCTION, not just by convention. Written per the user's explicit
2026-09-10 requirement that differences between models' pathway results be
interpretable as model-specific biological emphasis, not attribution-method
artifacts.

Design (supersedes the earlier gradient x input approach for both models):
  For a given cell, pick the SAME top-K expressed genes (by raw expression,
  descending) from a gene universe usable by BOTH architectures (i.e. the
  intersection of Geneformer's Ensembl-token-mappable genes and STATE's
  ESM2-protein-embedding-mappable genes). Feed this identical K-gene set
  into each model's own native input construction. For each gene g in the
  set, mask g out (architecture-specific mechanism, implemented in each
  model's own script) and measure:
      delta(g, cell) = s_full(cell) - s_masked_g(cell)
  where s = (pooled cell embedding) . w + b, w/b from ONE in-sample L2
  logistic regression fit per model on that model's own s_full embeddings
  across the whole cohort (same C=0.3 convention as the prior gradient
  pipeline). delta is aggregated (mean) per gene across all cell
  occurrences, filtered at MIN_GENE_OCCURRENCE, in the SAME output schema
  the existing GSEA runner (geneformer_run_gsea.py) already consumes.

Why K=256, not each model's native context length (512 for Geneformer,
2048 for STATE): a masking sweep costs O(K) full forward passes per cell
(one per candidate + one baseline), each O(seq_len^2) in attention -- at
native lengths this is intractable for STATE (2048^2 attention x 2049
passes x thousands of cells) and unnecessarily expensive for Geneformer.
Capping both models to the SAME K also removes "one model natively sees
512 genes/cell, the other sees 2048" as a confound on the comparison itself
-- this is a comparability improvement, not just a cost workaround.
"""
import json
import os
import pickle

import numpy as np
import pandas as pd
import scanpy as sc
import torch

RAW_H5AD = "/bigdata/godziklab/shared/Xinru/302006/01_DATA/raw/human_platelet_covid_severity_raw_counts.h5ad"
PATIENT_LIST_CSV = ("/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/"
                     "31_donor_level_revision/results_patient_v2/qc_baseline/patient_qc_distributions.csv")
SAMPLE_MAP_CSV = ("/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/"
                   "31_donor_level_revision/results/donor_meta_with_patient_id_v2.csv")
SYMBOL_TO_ENSEMBL_CSV = ("/bigdata/godziklab/shared/Xinru/302006/02_EMBEDDINGS/"
                          "option_c_monocyte/gene_symbol_to_ensembl.csv")
GENEFORMER_TOKEN_DICT = ("/bigdata/godziklab/shared/Xinru/302006/03_MODELS/geneformer/"
                          "geneformer-V2-316M-native/geneformer/token_dictionary_gc104M.pkl")
STATE_PROTEIN_EMBEDDINGS = "/bigdata/godziklab/shared/Xinru/302006/03_MODELS/SE-600M/protein_embeddings.pt"
SCGPT_VOCAB_JSON = ("/bigdata/godziklab/shared/Xinru/302006/03_MODELS/scgpt/"
                     "scGPT_human/scGPT_human/vocab.json")
UCE_HUMAN_PE = ("/bigdata/godziklab/shared/Xinru/302006/03_MODELS/uce/UCE/model_files/"
                 "protein_embeddings/Homo_sapiens.GRCh38.gene_symbol_to_embedding_ESM2.pt")
UCE_CHROM_CSV = "/bigdata/godziklab/shared/Xinru/302006/03_MODELS/uce/UCE/model_files/species_chrom.csv"

CELLS_PER_PATIENT = 60  # same convention as the gradient-based scripts
RNG_SEED = 20260910
TOP_K_GENES = 256
MIN_GENE_OCCURRENCE = 30


def load_canonical_cohort(max_cells=None):
    """Identical cohort restriction + stratified subsample used by both the
    prior gradient-based scripts and this masking-based pipeline, for
    direct comparability with earlier work. max_cells is a testing-only
    hook (None = full cohort, unchanged default behavior) to slice down to
    a small subset for cheap standalone validation before a full GPU run."""
    patient_list = pd.read_csv(PATIENT_LIST_CSV)
    patient_list = patient_list[patient_list["disease"] == "COVID-19"]
    assert len(patient_list) == 240 and int(patient_list["y"].sum()) == 25, \
        f"cohort mismatch: n={len(patient_list)}, fatal={int(patient_list['y'].sum())}"
    y_by_patient = patient_list.set_index("patient_id")["y"].to_dict()
    covid_patient_ids = set(patient_list["patient_id"].astype(str))
    print(f"Canonical cohort: {len(covid_patient_ids)} patients, "
          f"{sum(y_by_patient.values())} fatal", flush=True)

    sample_map = pd.read_csv(SAMPLE_MAP_CSV)
    sample_map["patient_id"] = sample_map["patient_id"].astype(str)
    sample_map = sample_map[sample_map["patient_id"].isin(covid_patient_ids)]
    sample_to_patient = sample_map.set_index("sample_id")["patient_id"].to_dict()

    print(f"Loading {RAW_H5AD} ...", flush=True)
    a = sc.read_h5ad(RAW_H5AD)
    a.obs["patient_id"] = a.obs["Sample_ID"].astype(str).map(sample_to_patient)
    a = a[a.obs["patient_id"].notna()].copy()
    a.obs["y"] = a.obs["patient_id"].map(y_by_patient).astype(int)
    n_patients_present = a.obs["patient_id"].nunique()
    print(f"Cells mapped to canonical cohort: {a.n_obs} cells, "
          f"{n_patients_present} / 240 patients present in this h5ad", flush=True)

    rng = np.random.default_rng(RNG_SEED)
    keep_idx = []
    for pid, grp in a.obs.groupby("patient_id"):
        n = len(grp)
        take = min(n, CELLS_PER_PATIENT)
        chosen = rng.choice(grp.index.to_numpy(), size=take, replace=False)
        keep_idx.extend(chosen.tolist())
    a = a[keep_idx].copy()
    print(f"After stratified per-patient subsample (cap={CELLS_PER_PATIENT}/patient): "
          f"{a.n_obs} cells, {a.obs['y'].sum()} from fatal patients' cells "
          f"(fatal patients contribute {a.obs.loc[a.obs.y == 1, 'patient_id'].nunique()} / "
          f"{a.obs['patient_id'].nunique()} patients)", flush=True)
    if max_cells is not None and a.n_obs > max_cells:
        # keep at least one fatal-patient cell if possible, for a meaningful smoke test
        fatal_idx = np.where(a.obs["y"].to_numpy() == 1)[0]
        other_idx = np.where(a.obs["y"].to_numpy() == 0)[0]
        n_fatal_keep = min(len(fatal_idx), max(1, max_cells // 4))
        keep = np.concatenate([fatal_idx[:n_fatal_keep], other_idx[:max_cells - n_fatal_keep]])
        a = a[np.sort(keep)].copy()
        print(f"[TESTING] sliced down to {a.n_obs} cells for validation "
              f"({int((a.obs['y']==1).sum())} from fatal patients)", flush=True)
    a.obs_names_make_unique()
    return a, y_by_patient


def shared_gene_universe(var_names):
    """
    Boolean mask over var_names: True where the gene symbol is usable by
    BOTH Geneformer (symbol -> Ensembl -> token_dictionary_gc104M.pkl) and
    STATE (symbol -> protein_embeddings.pt directly). Only genes in this
    intersection are eligible to be selected into any cell's top-K set, so
    the two models are always attributing the exact same candidate gene per
    occurrence, not just the same count of candidates.
    """
    sym2ens = pd.read_csv(SYMBOL_TO_ENSEMBL_CSV).set_index("gene_symbol")["ensembl_id"].to_dict()
    with open(GENEFORMER_TOKEN_DICT, "rb") as f:
        token_dict = pickle.load(f)
    protein_embeds = torch.load(STATE_PROTEIN_EMBEDDINGS, weights_only=False, map_location="cpu")

    geneformer_ok = np.array([sym2ens.get(g) in token_dict for g in var_names])
    state_ok = np.array([g in protein_embeds for g in var_names])
    both_ok = geneformer_ok & state_ok
    print(f"Shared gene universe: {both_ok.sum()} / {len(var_names)} genes usable by "
          f"both Geneformer ({geneformer_ok.sum()}) and STATE ({state_ok.sum()})", flush=True)
    assert both_ok.sum() > 5000, "shared universe too small -- mapping broke somewhere"
    return both_ok, sym2ens, token_dict, protein_embeds


def shared_gene_universe_3way(var_names):
    """
    Extends shared_gene_universe() to a 3-way (Geneformer & STATE & scGPT)
    intersection, for scGPT onward (and future models). NOT used by
    geneformer_masking_attribution.py / state_masking_attribution.py --
    those two already-committed, already-run (or in-progress) scripts use
    the original 2-way shared_gene_universe() above and are NOT being
    retroactively rerun for this: checked empirically (2026-09-10) that the
    3-way intersection drops only 7 genes out of 11,984 (0.06%: CYB561D2,
    MATR3, RGS5, SLFN12L, SPATA13, TBCE, TMSB15B), none of which appear in
    either model's top/bottom attribution lists so far -- judged negligible
    relative to discarding hours of completed/in-progress GPU compute.
    New models (scGPT here, UCE/scPRINT-2 next) should call THIS function,
    which will keep shrinking (properly) as each new model's own vocabulary
    is intersected in.
    """
    valid_mask, sym2ens, token_dict, protein_embeds = shared_gene_universe(var_names)
    with open(SCGPT_VOCAB_JSON) as f:
        scgpt_vocab = json.load(f)
    scgpt_ok = np.array([g in scgpt_vocab for g in var_names])
    three_way = valid_mask & scgpt_ok
    print(f"3-way shared gene universe (+scGPT): {three_way.sum()} / {len(var_names)} "
          f"(was {valid_mask.sum()} in the 2-way Geneformer/STATE universe, "
          f"scGPT alone covers {scgpt_ok.sum()})", flush=True)
    assert three_way.sum() > 5000, "3-way shared universe too small -- mapping broke somewhere"
    return three_way, sym2ens, token_dict, protein_embeds, scgpt_vocab


def shared_gene_universe_4way(var_names):
    """
    Extends shared_gene_universe_3way() to a 4-way (Geneformer & STATE &
    scGPT & UCE) intersection, for UCE onward. Not used by the three
    already-run scripts (Geneformer/STATE on the 2-way, scGPT on the
    3-way) -- check the actual delta each time before deciding whether an
    earlier model needs a rerun (see shared_gene_universe_3way's docstring
    for the precedent: a small delta was judged not worth discarding
    completed GPU compute for).

    UCE eligibility requires BOTH (a) the gene symbol (uppercased) exists
    in UCE's own human ESM2 protein-embedding dict, and (b) the gene
    symbol (uppercased) exists in UCE's species_chrom.csv reference table
    for species=="human" (needed for chromosome/genomic-start ordering --
    UCE's own tokenization ORDERS tokens by chromosome then genomic
    position, not by expression rank like the other 3 models). Both are
    genuine requirements verified by reading UCE's own
    adata_path_to_prot_chrom_starts() (data_proc/data_utils.py), which
    raises if either lookup fails -- a gene missing from either file is not
    usable by UCE at all, not just under-preferred.
    """
    four_way_pre, sym2ens, token_dict, protein_embeds, scgpt_vocab = shared_gene_universe_3way(var_names)
    uce_pe_raw = torch.load(UCE_HUMAN_PE, weights_only=False, map_location="cpu")
    uce_pe_upper = {k.upper() for k in uce_pe_raw.keys()}

    chrom_df = pd.read_csv(UCE_CHROM_CSV)
    human_chrom_genes_upper = set(
        chrom_df.loc[chrom_df["species"] == "human", "gene_symbol"].str.upper()
    )

    uce_ok = np.array([
        g.upper() in uce_pe_upper and g.upper() in human_chrom_genes_upper
        for g in var_names
    ])
    four_way = four_way_pre & uce_ok
    print(f"4-way shared gene universe (+UCE): {four_way.sum()} / {len(var_names)} "
          f"(was {four_way_pre.sum()} in the 3-way universe, UCE alone covers "
          f"{uce_ok.sum()})", flush=True)
    assert four_way.sum() > 5000, "4-way shared universe too small -- mapping broke somewhere"
    return four_way, sym2ens, token_dict, protein_embeds, scgpt_vocab, uce_pe_raw, chrom_df


def shared_gene_universe_5way(var_names, scprint2_human_ensembl_ids):
    """
    Extends shared_gene_universe_4way() to a 5-way (+ scPRINT-2) intersection,
    the final model. Unlike the other 4 models, scPRINT-2's gene vocabulary
    (`model.genes`) is only obtainable by loading the full checkpoint
    (~10 minutes on this cluster, confirmed empirically -- no lightweight
    vocab file ships separately), so this function takes the already-loaded
    model's human Ensembl gene ID set as an explicit argument (extracted by
    the calling script right after it loads the model for its own use)
    rather than loading anything itself -- avoids a redundant ~10-minute
    reload just to check membership.

    scprint2_human_ensembl_ids should be the SET of standard human Ensembl
    gene IDs (regex ^ENSG\\d+$ -- verified empirically that a naive
    startswith("ENSG") check wrongly also matches other species whose
    prefixes happen to start with the same 4 letters, e.g. chicken's
    ENSGALG...) present in model.genes.
    """
    five_way_pre, sym2ens, token_dict, protein_embeds, scgpt_vocab, uce_pe_raw, chrom_df = \
        shared_gene_universe_4way(var_names)
    scprint2_ok = np.array([
        sym2ens.get(g) in scprint2_human_ensembl_ids for g in var_names
    ])
    five_way = five_way_pre & scprint2_ok
    print(f"5-way shared gene universe (+scPRINT-2): {five_way.sum()} / {len(var_names)} "
          f"(was {five_way_pre.sum()} in the 4-way universe, scPRINT-2 alone covers "
          f"{scprint2_ok.sum()})", flush=True)
    assert five_way.sum() > 5000, "5-way shared universe too small -- mapping broke somewhere"
    return five_way, sym2ens, token_dict, protein_embeds, scgpt_vocab, uce_pe_raw, chrom_df


def top_k_genes_per_cell(X, valid_mask, k=TOP_K_GENES):
    """
    X: dense or sparse (n_cells, n_genes) raw-count matrix.
    valid_mask: boolean (n_genes,) from shared_gene_universe().
    Returns: list of length n_cells, each a 1D int array of up to k
    var_names column-indices, sorted by descending expression, restricted
    to valid_mask and to genes with nonzero expression in that cell
    (deterministic -- no randomness, so both model scripts independently
    recomputing this get identical candidate sets without needing to share
    a file).

    Uses kind='stable' (not the numpy default 'quicksort') explicitly:
    platelet single-cell data has heavy ties at low counts right at the
    top-K boundary, and quicksort's tie-break order is not guaranteed
    identical across numpy builds/versions -- confirmed a real risk here,
    since scpeft_env (numpy 1.26.4) and state_env (numpy 2.3.5) are
    genuinely different major versions. A stable sort breaks ties by
    original input order, which for `nz` (from np.where, ascending) means
    ties resolve by var_names column index -- reproducible across any
    numpy version/platform, since it depends only on the stable-sort
    guarantee, not on quicksort's internal partitioning.
    """
    n_cells = X.shape[0]
    out = []
    for i in range(n_cells):
        row = X[i].toarray().ravel() if hasattr(X[i], "toarray") else np.asarray(X[i]).ravel()
        row = row.copy()
        row[~valid_mask] = -1.0
        nz = np.where(row > 0)[0]
        if len(nz) == 0:
            out.append(np.array([], dtype=np.int64))
            continue
        order = nz[np.argsort(-row[nz], kind="stable")][:k]
        out.append(order)
    return out


def save_checkpoint_atomic(path, **arrays):
    """np.savez to a temp file then os.replace -- a raw np.savez(path, ...)
    can be killed mid-write by SLURM's walltime signal, leaving a truncated
    file that crashes the NEXT resume attempt (and the sbatch wrapper's
    set -e would then let the chain die there, needing manual intervention).
    The time-budget check that triggers this save leaves ~1800s of margin
    before the actual SLURM walltime limit in this project's own chained
    jobs, so the risk window is small in practice, but the fix is cheap.
    Passes an open file OBJECT (not a string path) to np.savez -- numpy
    only auto-appends '.npz' for string paths, and path.tmp doesn't end in
    .npz so a string-path call would silently write to
    '<path>.tmp.npz' instead, breaking the os.replace below."""
    tmp_path = f"{path}.tmp"
    with open(tmp_path, "wb") as f:
        np.savez(f, **arrays)
    os.replace(tmp_path, path)


def aggregate_and_save(gene_sum, gene_sumsq, gene_count, gene_symbols, out_path,
                        min_occurrence=MIN_GENE_OCCURRENCE):
    """Same ranked-list schema as the gradient-based scripts, so the
    existing model-agnostic geneformer_run_gsea.py works unchanged."""
    keep = gene_count >= min_occurrence
    mean_attr = np.divide(gene_sum, gene_count, out=np.zeros_like(gene_sum), where=gene_count > 0)
    var_attr = (np.divide(gene_sumsq, gene_count, out=np.zeros_like(gene_sum), where=gene_count > 0)
                - mean_attr ** 2)
    sd_attr = np.sqrt(np.clip(var_attr, 0, None))

    ranked = pd.DataFrame({
        "gene": [gene_symbols[i] for i in np.where(keep)[0]],
        "mean_attribution": mean_attr[keep],
        "sd_attribution": sd_attr[keep],
        "n_occurrences": gene_count[keep],
    }).sort_values("mean_attribution", ascending=False)

    print(f"Genes in ranked list (>= {min_occurrence} occurrences): {len(ranked)}", flush=True)
    print(f"Attribution (delta-score) range: [{ranked.mean_attribution.min():.4f}, "
          f"{ranked.mean_attribution.max():.4f}]", flush=True)
    ranked.to_csv(out_path, index=False)
    print(f"Saved {out_path}", flush=True)
    print("\nTop 20 genes by mean delta-score:", flush=True)
    print(ranked.head(20).to_string(index=False), flush=True)
    print("\nBottom 10 genes by mean delta-score:", flush=True)
    print(ranked.tail(10).to_string(index=False), flush=True)
    return ranked

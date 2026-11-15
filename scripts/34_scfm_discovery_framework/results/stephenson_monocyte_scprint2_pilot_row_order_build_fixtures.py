"""Bounded (h5py-only) fixture builder for
run_scprint2_row_order_smoke_test.sh (Packet B).

`anndata.read_h5ad(..., backed="r")` OOM-killed (exit 137, signal 9) on this
cluster's login node even for a 150-row slice of the real 5.7G frozen input --
confirmed empirically 2026-10-18, consistent with this project's standing rule
not to open a multi-GB h5ad via anndata on the login node (DEVELOPMENT.md). This
reads only the bounded slices actually needed: obs columns (small), var
(small), and the FIRST (n_cells+1) rows of the CSR X matrix (bounded by
indptr, never the full ~700MB data array), then builds a fresh, small AnnData
purely in memory. Verified against the real frozen input (2026-10-18): exit
0, ~1GB peak RSS, ~50s, correct 151-cell fixtures with a verified content
duplicate and a verified non-reversal/non-identity permutation.

Writes, into `run_dir`:
- `genelist.json`: real `ensembl_gene_id` values (NOT symbols -- see the
  gene-namespace note below) that `scprint2 embed --how some --genelist=...`
  will actually match against.
- `fixture_natural.h5ad`: the first `n_cells` real cells plus one injected
  content-duplicate of cell `DUP_SOURCE_POS`, appended at the end with a
  distinct `cell_id_stripped` (and hence a distinct recovered id) but
  IDENTICAL expression content -- an independent, within-run identity oracle
  (see `check_duplicate_row_self_consistency` in
  `stephenson_monocyte_scprint2_pilot_row_order_verify.py`).
- `fixture_order_b.h5ad`: the SAME (n_cells+1) cells, reordered by a fixed
  random derangement (seeded by `perm_seed`), explicitly asserted to differ
  from both the identity and the pure reversal permutations -- a pure
  reversal is the one specific permutation a "producer always reverses
  regardless of input order" bug would exactly cancel against when compared
  to a natural-order run (demonstrated in
  review:20261018_000122's independent review).
- `duplicate_pair.json`: `{"original_id": ..., "duplicate_id": ...}` for the
  injected duplicate, in final (post-`steph_` prefix) recovered-id form.

Gene-namespace note (B3): the installed `scdataloader.Preprocessor`, with
`--is_symbol true` (this project's pinned setting), REPLACES `adata.var.index`
with `ensembl_gene_id` via
`genesdf.drop_duplicates("symbol").set_index("symbol")`
(`preprocess.py:268-277`) BEFORE the Collator ever runs; `Collator._setup`
matches `genelist` against that SAME ensembl-indexed `load_genes()` table via
`df.index.isin(genelist)` (`collator.py:102-140`, both read directly,
2026-10-17/18). Passing raw HGNC symbols as `--genelist` would match ZERO
genes post-preprocessing. This reproduces the EXACT same resolution
(including pandas' `drop_duplicates(keep="first")` default, matched exactly
rather than approximated) so the genelist this script emits is the real,
verified-present set the installed Collator will actually accept.
"""
import json
import sys

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix

DUP_SOURCE_POS = 5


def read_obs_col(group, key):
    x = group[key]
    if isinstance(x, h5py.Dataset):
        return x.asstr()[:] if x.dtype.kind in "OSU" else x[:]
    cats = x["categories"].asstr()[:]
    codes = x["codes"][:]
    return np.array([cats[c] if c >= 0 else None for c in codes], dtype=object)


def build(frozen_input_path, n_cells, n_symbols, min_hits, perm_seed, run_dir):
    with h5py.File(frozen_input_path, "r") as f:
        obs_group = f["obs"]
        obs_index_name = obs_group.attrs["_index"]
        obs_cols = [c for c in obs_group.attrs["column-order"] if c != obs_index_name]
        obs_ids_full = read_obs_col(obs_group, obs_index_name)
        obs_data = {c: read_obs_col(obs_group, c) for c in obs_cols}

        var_group = f["var"]
        var_index_name = var_group.attrs["_index"]
        var_names = read_obs_col(var_group, var_index_name)
        var_cols = [c for c in var_group.attrs["column-order"] if c != var_index_name]
        var_data = {c: read_obs_col(var_group, c) for c in var_cols}

        n_vars = len(var_names)
        indptr_full = f["X"]["indptr"][: n_cells + 2]  # bounded: only the rows we need
        row_start, row_end = int(indptr_full[0]), int(indptr_full[n_cells])
        x_data = f["X"]["data"][row_start:row_end]
        x_indices = f["X"]["indices"][row_start:row_end]
        x_indptr = indptr_full[: n_cells + 1] - row_start

    obs_df = pd.DataFrame({c: v[:n_cells] for c, v in obs_data.items()}, index=obs_ids_full[:n_cells])
    var_df = pd.DataFrame({c: v for c, v in var_data.items()}, index=var_names)
    X = csr_matrix((x_data, x_indices, x_indptr), shape=(n_cells, n_vars))
    natural = ad.AnnData(X=X, obs=obs_df, var=var_df)

    from scdataloader.utils import load_genes

    table = load_genes(organisms="NCBITaxon:9606")
    table = table.copy()
    table["ensembl_gene_id"] = table.index
    resolved = table.drop_duplicates("symbol").set_index("symbol", drop=False)

    symbol_genelist = list(var_names[:n_symbols])
    ensembl_genelist = [resolved.loc[s, "ensembl_gene_id"] for s in symbol_genelist if s in resolved.index]
    unmapped = [s for s in symbol_genelist if s not in resolved.index]
    assert len(ensembl_genelist) >= min_hits, (
        f"only {len(ensembl_genelist)} of {len(symbol_genelist)} symbols resolved to a real "
        f"ensembl_gene_id (unmapped: {unmapped}) -- below the required minimum {min_hits}"
    )
    print(f"Gene namespace preflight: {len(ensembl_genelist)}/{len(symbol_genelist)} symbols resolved "
          f"to a real ensembl_gene_id ({len(unmapped)} unmapped: {unmapped[:5]}...)")

    with open(f"{run_dir}/genelist.json", "w") as fh:
        json.dump(ensembl_genelist, fh)

    dup_row = natural[DUP_SOURCE_POS : DUP_SOURCE_POS + 1].copy()
    original_stripped = str(dup_row.obs["cell_id_stripped"].iloc[0])
    dup_row.obs["cell_id_stripped"] = original_stripped + "_DUP_FOR_ROWORDER_TEST"
    dup_row.obs_names = [f"steph_{dup_row.obs['cell_id_stripped'].iloc[0]}"]

    natural_with_dup = ad.concat([natural, dup_row], join="outer")
    n_total = natural_with_dup.n_obs

    with open(f"{run_dir}/duplicate_pair.json", "w") as fh:
        json.dump(
            {"original_id": f"steph_{original_stripped}", "duplicate_id": str(dup_row.obs_names[0])},
            fh,
        )

    natural_with_dup.write_h5ad(f"{run_dir}/fixture_natural.h5ad")

    rng = np.random.default_rng(perm_seed)
    perm = rng.permutation(n_total)
    reversal = np.arange(n_total)[::-1]
    assert not np.array_equal(perm, reversal), (
        f"seed {perm_seed}'s permutation coincided with the reversal -- pick a different seed"
    )
    assert not np.array_equal(perm, np.arange(n_total)), "permutation coincided with identity"
    order_b = natural_with_dup[perm].copy()
    order_b.write_h5ad(f"{run_dir}/fixture_order_b.h5ad")

    print(f"Fixtures ready: {n_total} cells each ({n_cells} real + 1 injected content-duplicate), "
          f"order_b is a fixed derangement (seed={perm_seed}), {len(ensembl_genelist)} real "
          f"ensembl gene ids pinned for --genelist")


if __name__ == "__main__":
    if len(sys.argv) != 7:
        print(
            "usage: stephenson_monocyte_scprint2_pilot_row_order_build_fixtures.py "
            "<frozen_input.h5ad> <n_cells> <n_symbols_to_try> <min_gene_hits_required> "
            "<permutation_seed> <run_dir>",
            file=sys.stderr,
        )
        sys.exit(2)
    build(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]), sys.argv[6])

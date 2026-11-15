"""Tests for score_activation_markers.py — module score logic + dim correlation.

Pure-math unit tests; no HPCC embeddings required.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(
    0,
    str(Path(__file__).parent.parent / "scripts/24_option_c_canonical_input/code"),
)
from score_activation_markers import (
    PLATELET_ACTIVATION_GENES,
    compute_module_score,
    find_embedding,
)


class TestPlateletActivationGenes:
    """The hard-coded marker dictionary should contain canonical platelet markers."""

    def test_core_markers_present(self):
        core = PLATELET_ACTIVATION_GENES["core_activation"]
        # Canonical platelet activation markers
        for m in ("PF4", "PPBP", "SELP", "ITGA2B", "GP1BA"):
            assert m in core, f"{m} missing from core_activation set"

    def test_all_groups_non_empty(self):
        for name, genes in PLATELET_ACTIVATION_GENES.items():
            assert len(genes) > 0, f"{name} has no genes"
            assert all(isinstance(g, str) for g in genes)

    def test_no_duplicate_genes_within_group(self):
        for name, genes in PLATELET_ACTIVATION_GENES.items():
            assert len(genes) == len(set(genes)), f"{name} has duplicates"

    def test_groups_are_disjoint_or_intentional(self):
        """Sanity check — overlap between groups should be zero or documented."""
        seen = {}
        for name, genes in PLATELET_ACTIVATION_GENES.items():
            for g in genes:
                if g in seen:
                    pytest.fail(f"Gene {g!r} duplicated in {name!r} and {seen[g]!r}")
                seen[g] = name


class TestComputeModuleScore:
    """compute_module_score() runs sc.tl.score_genes; returns (n_obs,) array."""

    def _make_adata(self, n_cells=100, gene_names=None):
        import anndata

        if gene_names is None:
            gene_names = [f"GENE{i}" for i in range(50)]
        rng = np.random.default_rng(0)
        X = rng.poisson(2, size=(n_cells, len(gene_names))).astype(np.float32)
        obs = pd.DataFrame(index=[f"c{i}" for i in range(n_cells)])
        var = pd.DataFrame(index=gene_names)
        return anndata.AnnData(X=X, obs=obs, var=var)

    def test_score_shape_matches_n_obs(self):
        gene_names = ["PF4", "PPBP", "SELP", "GENE1", "GENE2", "GENE3"]
        adata = self._make_adata(n_cells=50, gene_names=gene_names)
        scores = compute_module_score(adata, "core", ["PF4", "PPBP", "SELP"])
        assert scores.shape == (50,)
        assert np.isfinite(scores).all()

    def test_score_returns_zeros_when_no_markers_in_var(self):
        adata = self._make_adata(n_cells=20, gene_names=["GENE1", "GENE2"])
        scores = compute_module_score(adata, "core", ["PF4", "SELP"])  # neither in var
        assert scores.shape == (20,)
        assert (scores == 0).all()

    def test_score_filters_to_valid_markers(self):
        """Markers absent from var_names are silently skipped, not errored on."""
        gene_names = ["PF4", "GENE1"]
        adata = self._make_adata(n_cells=10, gene_names=gene_names)
        # Only PF4 will be matched
        scores = compute_module_score(adata, "test", ["PF4", "MISSING_GENE"])
        assert scores.shape == (10,)
        assert np.isfinite(scores).all()


class TestFindEmbedding:
    """find_embedding() in score_activation_markers also handles scprint2 key."""

    def test_finds_scprint2_keys(self):
        import anndata

        n = 30
        d = 1600
        for key in ("X_scprint2", "X_scprint", "X_scPRINT", "X_scPRINT2"):
            a = anndata.AnnData(
                X=np.zeros((n, 50)),
                obsm={key: np.random.randn(n, d).astype(np.float32)},
            )
            X, found = find_embedding(a, "scprint2")
            assert found == key
            assert X.shape == (n, d)

"""Tests for embedding processing and cell protein embedding logic."""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).parent.parent / "scripts/22_esm_protein_multimodal/code"))


class TestExpressionWeightedEmbeddings:
    """Test the expression-weighted protein embedding aggregation (step 4 logic)."""

    def test_weighted_mean_basic(self):
        """Expression-weighted mean produces correct output."""
        # 3 cells, 2 genes, embedding dim 4
        expr_matrix = np.array(
            [
                [1.0, 0.0],  # cell 0: only gene A expressed
                [0.0, 1.0],  # cell 1: only gene B expressed
                [0.5, 0.5],  # cell 2: equal expression
            ]
        )
        esm_matrix = np.array(
            [
                [1.0, 0.0, 0.0, 0.0],  # gene A embedding
                [0.0, 1.0, 0.0, 0.0],  # gene B embedding
            ]
        )

        # Normalize expression
        row_sums = expr_matrix.sum(axis=1, keepdims=True)
        row_sums = np.maximum(row_sums, 1e-8)
        weights = expr_matrix / row_sums

        result = weights @ esm_matrix

        # Cell 0 should get gene A's embedding
        np.testing.assert_allclose(result[0], [1.0, 0.0, 0.0, 0.0])
        # Cell 1 should get gene B's embedding
        np.testing.assert_allclose(result[1], [0.0, 1.0, 0.0, 0.0])
        # Cell 2 should get average
        np.testing.assert_allclose(result[2], [0.5, 0.5, 0.0, 0.0])

    def test_zero_expression_handled(self):
        """Cells with zero expression don't cause NaN/Inf."""
        expr_matrix = np.array(
            [
                [0.0, 0.0],  # no expression
                [1.0, 2.0],  # normal
            ]
        )
        esm_matrix = np.array(
            [
                [1.0, 0.0],
                [0.0, 1.0],
            ]
        )

        row_sums = expr_matrix.sum(axis=1, keepdims=True)
        row_sums = np.maximum(row_sums, 1e-8)
        weights = expr_matrix / row_sums
        result = weights @ esm_matrix

        assert not np.isnan(result).any()
        assert not np.isinf(result).any()
        # Zero-expression cell should get near-zero embedding
        assert np.linalg.norm(result[0]) < 1e-4

    def test_negative_expression_clipped(self):
        """Negative expression values are clipped to 0."""
        expr_matrix = np.array([[-1.0, 2.0], [1.0, 1.0]])
        expr_matrix = np.maximum(expr_matrix, 0)
        esm_matrix = np.array([[1.0, 0.0], [0.0, 1.0]])

        row_sums = expr_matrix.sum(axis=1, keepdims=True)
        row_sums = np.maximum(row_sums, 1e-8)
        weights = expr_matrix / row_sums
        result = weights @ esm_matrix

        # First cell: only gene B contributes (gene A clipped to 0)
        np.testing.assert_allclose(result[0], [0.0, 1.0])

    def test_output_shape(self):
        """Output shape is (n_cells, esm_dim)."""
        n_cells, n_genes, esm_dim = 100, 50, 1280
        expr = np.random.poisson(3, size=(n_cells, n_genes)).astype(float)
        esm = np.random.randn(n_genes, esm_dim).astype(np.float32)

        row_sums = expr.sum(axis=1, keepdims=True)
        row_sums = np.maximum(row_sums, 1e-8)
        weights = expr / row_sums
        result = weights @ esm

        assert result.shape == (n_cells, esm_dim)

    def test_weights_sum_to_one(self):
        """Normalized expression weights sum to 1 per cell."""
        expr = np.random.poisson(5, size=(50, 20)).astype(float)
        expr = np.maximum(expr, 0)
        row_sums = expr.sum(axis=1, keepdims=True)
        row_sums = np.maximum(row_sums, 1e-8)
        weights = expr / row_sums

        np.testing.assert_allclose(weights.sum(axis=1), 1.0, atol=1e-6)


@pytest.mark.integration
class TestDEGIdentification:
    """Test DEG identification utilities (step 1 logic)."""

    def test_binary_label_creation(self):
        """Binary labels correctly map severity to 0/1."""
        import anndata
        from step1_identify_degs import create_binary_labels

        n = 100
        obs = pd.DataFrame(
            {"covid_severity": ["severe"] * 30 + ["mild"] * 40 + ["fatal"] * 30},
            index=[f"c{i}" for i in range(n)],
        )
        adata = anndata.AnnData(
            X=np.zeros((n, 10)),
            obs=obs,
        )
        labels = create_binary_labels(adata)
        assert labels.sum() == 60  # severe + fatal
        assert (1 - labels).sum() == 40  # mild

    def test_binary_label_case_insensitive(self):
        """Label mapping handles mixed case."""
        import anndata
        from step1_identify_degs import create_binary_labels

        obs = pd.DataFrame(
            {"covid_severity": ["Severe", "FATAL", "Mild", "moderate"]},
            index=[f"c{i}" for i in range(4)],
        )
        adata = anndata.AnnData(X=np.zeros((4, 5)), obs=obs)
        labels = create_binary_labels(adata)
        np.testing.assert_array_equal(labels, [1, 1, 0, 0])


@pytest.mark.integration
class TestProteinSequences:
    """Test protein sequence utilities (step 2 logic)."""

    def test_fasta_write_and_read(self):
        """FASTA write/read roundtrip preserves data."""
        import tempfile

        from step2_get_protein_sequences import write_fasta
        from step3_generate_esm_embeddings import read_fasta

        proteins = [
            {
                "uniprot_id": "P12345",
                "gene": "TP53",
                "protein_name": "Tumor protein p53",
                "sequence": "MEEPQSDPSVEPPLSQETFSDLWK" * 5,
            },
            {
                "uniprot_id": "Q67890",
                "gene": "BRCA1",
                "protein_name": "BRCA1",
                "sequence": "MDLSALREVE" * 10,
            },
        ]

        with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as f:
            write_fasta(proteins, f.name)
            sequences = read_fasta(f.name)

        assert len(sequences) == 2
        # Check sequences are preserved (minus newlines)
        assert sequences[0][1] == proteins[0]["sequence"]
        assert sequences[1][1] == proteins[1]["sequence"]

    def test_fasta_read_multiline(self):
        """FASTA reader handles multi-line sequences."""
        import tempfile

        from step3_generate_esm_embeddings import read_fasta

        content = ">P12345|TP53|test\nMEEPQSDP\nSVEPPLSQ\n>Q67890|BRCA1|test\nMDLSALREVE\n"
        with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as f:
            f.write(content)
            f.flush()
            sequences = read_fasta(f.name)

        assert len(sequences) == 2
        assert sequences[0][1] == "MEEPQSDPSVEPPLSQ"
        assert sequences[1][1] == "MDLSALREVE"

"""Shared fixtures for scfm-benchmark tests."""

import numpy as np
import pandas as pd
import pytest


@pytest.fixture
def synthetic_embeddings():
    """Create synthetic cell embeddings for testing (no HPCC data needed)."""
    np.random.seed(42)
    n_cells = 200
    models = {
        "STATE": np.random.randn(n_cells, 2058).astype(np.float32),
        "UCE": np.random.randn(n_cells, 1280).astype(np.float32),
        "Geneformer": np.random.randn(n_cells, 1152).astype(np.float32),
        "TranscriptFormer": np.random.randn(n_cells, 2048).astype(np.float32),
    }
    return models


@pytest.fixture
def binary_labels():
    """Create synthetic binary severity labels."""
    np.random.seed(42)
    n_cells = 200
    labels = np.zeros(n_cells, dtype=int)
    labels[:80] = 1  # 40% severe
    return labels


@pytest.fixture
def severity_labels():
    """Create synthetic multi-class severity labels."""
    np.random.seed(42)
    n_cells = 200
    classes = ["healthy", "mild", "moderate", "severe", "fatal"]
    return np.random.choice(classes, size=n_cells)


@pytest.fixture
def mortality_labels():
    """Binary in-hospital mortality (fatal=1 vs survived=0).

    Mirrors the primary task definition from main.tex Methods §"Primary Task:
    In-Hospital Mortality (Binary)": fatal=1; survived (mild+moderate+recovered+severe)=0;
    healthy and unknown are EXCLUDED from this fixture (since they're excluded
    from the actual analysis).
    """
    np.random.seed(42)
    n_cells = 200
    return np.random.choice([0, 1], size=n_cells, p=[0.7, 0.3])  # ~30% fatal


@pytest.fixture
def mortality_severities():
    """Synthetic 'covid_severity' string column for testing label-mask construction.

    Used to verify that build_binary() / create_binary_labels() correctly
    map fatal->1, survived-tokens->0, healthy/unknown->excluded.
    """
    np.random.seed(42)
    n = 200
    classes = ["fatal", "severe", "mild", "moderate", "recovered", "healthy", "unknown"]
    weights = [0.13, 0.40, 0.15, 0.10, 0.08, 0.07, 0.07]
    return np.random.choice(classes, size=n, p=weights)


@pytest.fixture
def patient_groups():
    """Synthetic Sample_ID groups (50 patients, ~4 cells each) for StratifiedGroupKFold."""
    np.random.seed(42)
    n_cells = 200
    return np.array([f"P{i // 4:03d}" for i in range(n_cells)])


@pytest.fixture
def mock_adata():
    """Create a minimal AnnData object for testing."""
    import anndata

    np.random.seed(42)
    n_cells, n_genes = 200, 500
    X = np.random.poisson(2, size=(n_cells, n_genes)).astype(np.float32)
    obs = pd.DataFrame(
        {"covid_severity": np.random.choice(["severe", "mild", "moderate"], n_cells)},
        index=[f"cell_{i}" for i in range(n_cells)],
    )
    var = pd.DataFrame(index=[f"GENE{i}" for i in range(n_genes)])
    adata = anndata.AnnData(X=X, obs=obs, var=var)
    return adata

"""Tests for benchmark evaluation utilities."""

import numpy as np
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler


class TestCrossValidation:
    """Test cross-validation evaluation pipeline."""

    def test_cv_returns_valid_auc(self, synthetic_embeddings, binary_labels):
        """CV evaluation returns AUC in [0, 1]."""
        X = synthetic_embeddings["STATE"]
        y = binary_labels
        cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
        # 5-fold CV on 200 cells -> train fold = 160 samples; PCA n_components
        # must be <= min(n_samples_per_fold, n_features). 64 is well below
        # both and keeps the test fast.
        clf = make_pipeline(
            StandardScaler(),
            PCA(n_components=64, random_state=42),
            LogisticRegression(max_iter=2000, random_state=42, class_weight="balanced"),
        )
        y_prob = cross_val_predict(clf, X, y, cv=cv, method="predict_proba")
        auc = roc_auc_score(y, y_prob[:, 1])
        assert 0.0 <= auc <= 1.0

    def test_cv_all_models_run(self, synthetic_embeddings, binary_labels):
        """All 4 models can be evaluated without error."""
        y = binary_labels
        for name, X in synthetic_embeddings.items():
            cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
            # PCA must fit min(train-fold n_samples, n_features); 5-fold on
            # 200 cells -> 160 train samples per fold.
            clf = make_pipeline(
                StandardScaler(),
                PCA(n_components=64, random_state=42),
                LogisticRegression(max_iter=2000, random_state=42, class_weight="balanced"),
            )
            y_prob = cross_val_predict(clf, X, y, cv=cv, method="predict_proba")
            auc = roc_auc_score(y, y_prob[:, 1])
            assert 0.0 <= auc <= 1.0, f"{name} produced invalid AUC: {auc}"


class TestEarlyFusion:
    """Test early fusion (concatenation) logic."""

    def test_concatenation_shape(self, synthetic_embeddings):
        """Concatenated embedding has sum of all dims."""
        X_concat = np.hstack([synthetic_embeddings[m] for m in sorted(synthetic_embeddings.keys())])
        expected_dim = sum(v.shape[1] for v in synthetic_embeddings.values())
        assert X_concat.shape == (200, expected_dim)

    def test_pca_reduces_dims(self, synthetic_embeddings):
        """PCA reduction produces target dimensionality."""
        X_concat = np.hstack([synthetic_embeddings[m] for m in sorted(synthetic_embeddings.keys())])
        # 200 cells x ~6538 concatenated features; PCA n_components capped at
        # min(n_samples, n_features) = 200. Use 128 to stay safely below.
        pca = PCA(n_components=128, random_state=42)
        X_reduced = pca.fit_transform(StandardScaler().fit_transform(X_concat))
        assert X_reduced.shape == (200, 128)


class TestLateFusion:
    """Test late fusion (probability averaging) logic."""

    def test_average_probs_valid(self):
        """Averaged probabilities stay in [0, 1]."""
        probs = {
            "A": np.random.uniform(0, 1, 100),
            "B": np.random.uniform(0, 1, 100),
            "C": np.random.uniform(0, 1, 100),
        }
        avg = np.mean(list(probs.values()), axis=0)
        assert avg.min() >= 0.0
        assert avg.max() <= 1.0

    def test_rank_fusion_valid(self):
        """Rank-based fusion produces valid scores."""
        from scipy.stats import rankdata

        probs = {
            "A": np.random.uniform(0, 1, 50),
            "B": np.random.uniform(0, 1, 50),
        }
        ranked = np.mean([rankdata(p) / len(p) for p in probs.values()], axis=0)
        assert ranked.min() > 0.0
        assert ranked.max() <= 1.0


class TestBootstrapCI:
    """Test bootstrap confidence interval logic."""

    def test_ci_contains_point_estimate(self):
        """95% CI should contain the point estimate for most seeds."""
        np.random.seed(42)
        y_true = np.array([1, 1, 1, 0, 0, 0, 1, 0, 1, 0] * 10)
        y_prob = y_true * 0.8 + (1 - y_true) * 0.2 + np.random.normal(0, 0.1, len(y_true))
        y_prob = np.clip(y_prob, 0, 1)

        point_auc = roc_auc_score(y_true, y_prob)

        # Bootstrap
        n_boot = 500
        aucs = []
        for _ in range(n_boot):
            idx = np.random.randint(0, len(y_true), len(y_true))
            if len(np.unique(y_true[idx])) < 2:
                continue
            aucs.append(roc_auc_score(y_true[idx], y_prob[idx]))

        ci_low = np.percentile(aucs, 2.5)
        ci_high = np.percentile(aucs, 97.5)

        assert ci_low <= point_auc <= ci_high
        assert ci_high - ci_low < 0.3  # CI shouldn't be absurdly wide

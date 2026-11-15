"""Tests for the attention-based multimodal fusion module."""

import sys
from pathlib import Path

import numpy as np
import pytest

torch = pytest.importorskip("torch")

# Add module path
sys.path.insert(0, str(Path(__file__).parent.parent / "scripts/18_multimodal_fusion/code"))
from attention_fusion_module import AttentionFusion, GatedFusion

MODEL_DIMS = {"STATE": 2058, "UCE": 1280, "Geneformer": 1152, "TranscriptFormer": 2048}


class TestAttentionFusion:
    """Test suite for AttentionFusion module."""

    def test_forward_shape(self):
        """Output logits have correct shape [batch, n_classes]."""
        model = AttentionFusion(MODEL_DIMS, common_dim=256, n_heads=4, n_classes=2)
        batch_size = 16
        inputs = {name: torch.randn(batch_size, dim) for name, dim in MODEL_DIMS.items()}
        logits = model(inputs)
        assert logits.shape == (batch_size, 2)

    def test_return_attention_weights(self):
        """Attention weights returned when requested, sum to ~1 per cell."""
        model = AttentionFusion(MODEL_DIMS, common_dim=256, n_heads=4)
        # eval() disables attention-weight dropout — softmax sum-to-1 only
        # holds without the dropout zeroing pattern (training-mode dropout
        # rescales surviving weights by 1/(1-p), breaking the invariant).
        model.eval()
        batch_size = 8
        inputs = {name: torch.randn(batch_size, dim) for name, dim in MODEL_DIMS.items()}
        logits, fused, attn = model(inputs, return_attention=True)

        assert attn.shape == (batch_size, len(MODEL_DIMS))
        # Attention weights should sum to approximately 1 (softmax output)
        sums = attn.sum(dim=1).detach().numpy()
        np.testing.assert_allclose(sums, 1.0, atol=1e-5)

    def test_fused_embedding_shape(self):
        """Fused embedding has correct common_dim."""
        common_dim = 128
        model = AttentionFusion(MODEL_DIMS, common_dim=common_dim)
        inputs = {name: torch.randn(4, dim) for name, dim in MODEL_DIMS.items()}
        _, fused, _ = model(inputs, return_attention=True)
        assert fused.shape == (4, common_dim)

    def test_missing_modality(self):
        """Model handles missing modality input (zeros fallback)."""
        model = AttentionFusion(MODEL_DIMS, common_dim=256)
        # Only provide 2 of 4 modalities
        partial_inputs = {
            "STATE": torch.randn(4, 2058),
            "UCE": torch.randn(4, 1280),
        }
        logits = model(partial_inputs)
        assert logits.shape == (4, 2)
        assert not torch.isnan(logits).any()

    def test_5_model_extension(self):
        """Fusion works with 5 modalities (ESM_protein added)."""
        dims_5 = {**MODEL_DIMS, "ESM_protein": 1280}
        model = AttentionFusion(dims_5, common_dim=256, n_heads=4)
        model.eval()  # disable attention dropout; required for sum-to-1 check
        inputs = {name: torch.randn(8, dim) for name, dim in dims_5.items()}
        logits, fused, attn = model(inputs, return_attention=True)

        assert logits.shape == (8, 2)
        assert attn.shape == (8, 5)
        np.testing.assert_allclose(attn.sum(dim=1).detach().numpy(), 1.0, atol=1e-5)

    def test_gradient_flow(self):
        """Gradients flow through all projection layers."""
        model = AttentionFusion(MODEL_DIMS, common_dim=256)
        inputs = {name: torch.randn(4, dim, requires_grad=True) for name, dim in MODEL_DIMS.items()}
        logits = model(inputs)
        loss = logits.sum()
        loss.backward()

        for name in MODEL_DIMS:
            assert inputs[name].grad is not None
            assert inputs[name].grad.abs().sum() > 0

    def test_deterministic_eval(self):
        """Model produces same output in eval mode with same input."""
        model = AttentionFusion(MODEL_DIMS, common_dim=256, dropout=0.1)
        model.eval()
        inputs = {name: torch.randn(4, dim) for name, dim in MODEL_DIMS.items()}

        with torch.no_grad():
            out1 = model(inputs)
            out2 = model(inputs)
        torch.testing.assert_close(out1, out2)


class TestGatedFusion:
    """Test suite for GatedFusion module."""

    def test_forward_shape(self):
        """Gate fusion output has correct shape."""
        model = GatedFusion(MODEL_DIMS, common_dim=256, n_classes=2)
        inputs = {name: torch.randn(8, dim) for name, dim in MODEL_DIMS.items()}
        logits = model(inputs)
        assert logits.shape == (8, 2)

    def test_gates_sum_to_one(self):
        """Gating weights sum to 1 per cell."""
        model = GatedFusion(MODEL_DIMS, common_dim=256)
        model.eval()  # consistency with AttentionFusion tests; gate has no dropout but harmless
        inputs = {name: torch.randn(8, dim) for name, dim in MODEL_DIMS.items()}
        _, _, gates = model(inputs, return_gates=True)
        sums = gates.sum(dim=1).detach().numpy()
        np.testing.assert_allclose(sums, 1.0, atol=1e-5)

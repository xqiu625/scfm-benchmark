#!/usr/bin/env python3
"""
ACC 18: Attention-Based Multimodal Fusion Module

Learn per-cell attention weights over foundation model embeddings.
Key differentiator: interpretable attention reveals which FM matters for which cells.

Architecture:
    Input: N model embeddings per cell [X_1, X_2, ..., X_N]
    -> Project each to common dim (256) via learned linear layers
    -> Cross-attention: learned query attends over model embeddings
    -> Output: fused 256-dim embedding + attention weights

Author: Xinru Qiu
Date: February 2026
"""

import torch
import torch.nn as nn
import torch.nn.functional as F
import numpy as np


class AttentionFusion(nn.Module):
    """
    Attention-based fusion of multiple foundation model embeddings.

    For each cell, learns to weight different model representations
    based on their relevance to the classification task.
    """

    def __init__(self, model_dims, common_dim=256, n_heads=4, n_classes=2, dropout=0.1):
        """
        Args:
            model_dims: dict mapping model_name -> embedding_dim
                e.g., {'STATE': 2058, 'UCE': 1280, 'Geneformer': 1152, 'TranscriptFormer': 2048}
            common_dim: projection dimension for all models
            n_heads: number of attention heads
            n_classes: classification output dimension
            dropout: dropout rate
        """
        super().__init__()

        self.model_names = sorted(model_dims.keys())
        self.n_models = len(self.model_names)
        self.common_dim = common_dim

        # Per-model projection layers
        self.projections = nn.ModuleDict({
            name: nn.Sequential(
                nn.Linear(dim, common_dim),
                nn.LayerNorm(common_dim),
                nn.GELU(),
                nn.Dropout(dropout),
            )
            for name, dim in model_dims.items()
        })

        # Learned query for cross-attention
        self.query = nn.Parameter(torch.randn(1, 1, common_dim) * 0.02)

        # Multi-head cross-attention
        self.cross_attention = nn.MultiheadAttention(
            embed_dim=common_dim,
            num_heads=n_heads,
            dropout=dropout,
            batch_first=True,
        )

        # Layer norm after attention
        self.attn_norm = nn.LayerNorm(common_dim)

        # Classification head
        self.classifier = nn.Sequential(
            nn.Linear(common_dim, common_dim // 2),
            nn.LayerNorm(common_dim // 2),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(common_dim // 2, n_classes),
        )

    def forward(self, embeddings_dict, return_attention=False):
        """
        Forward pass.

        Args:
            embeddings_dict: dict mapping model_name -> tensor [batch, model_dim]
            return_attention: if True, also return attention weights

        Returns:
            logits: [batch, n_classes]
            fused_embedding: [batch, common_dim] (optional)
            attention_weights: [batch, n_models] (optional)
        """
        batch_size = next(iter(embeddings_dict.values())).shape[0]

        # Project each model to common dimension
        projected = []
        for name in self.model_names:
            if name in embeddings_dict:
                proj = self.projections[name](embeddings_dict[name])
                projected.append(proj)
            else:
                # Missing model: use zeros
                projected.append(torch.zeros(batch_size, self.common_dim,
                                           device=next(self.parameters()).device))

        # Stack as key/value: [batch, n_models, common_dim]
        kv = torch.stack(projected, dim=1)

        # Query: [batch, 1, common_dim]
        query = self.query.expand(batch_size, -1, -1)

        # Cross-attention
        fused, attn_weights = self.cross_attention(
            query=query,
            key=kv,
            value=kv,
        )
        # fused: [batch, 1, common_dim], attn_weights: [batch, 1, n_models]

        fused = fused.squeeze(1)  # [batch, common_dim]
        fused = self.attn_norm(fused)

        attn_weights = attn_weights.squeeze(1)  # [batch, n_models]

        # Classification
        logits = self.classifier(fused)

        if return_attention:
            return logits, fused, attn_weights

        return logits

    def get_model_attention(self, embeddings_dict):
        """Get per-cell attention weights for interpretability."""
        _, _, attn_weights = self.forward(embeddings_dict, return_attention=True)
        return {
            name: attn_weights[:, i].detach().cpu().numpy()
            for i, name in enumerate(self.model_names)
        }


class GatedFusion(nn.Module):
    """
    Gated fusion: learn model-specific gates for each cell.
    Simpler alternative to attention fusion.
    """

    def __init__(self, model_dims, common_dim=256, n_classes=2, dropout=0.1):
        super().__init__()

        self.model_names = sorted(model_dims.keys())
        n_models = len(self.model_names)
        total_dim = sum(model_dims.values())

        # Per-model projection
        self.projections = nn.ModuleDict({
            name: nn.Linear(dim, common_dim)
            for name, dim in model_dims.items()
        })

        # Gate network: takes all embeddings, produces per-model weights
        self.gate = nn.Sequential(
            nn.Linear(total_dim, 128),
            nn.ReLU(),
            nn.Linear(128, n_models),
            nn.Softmax(dim=-1),
        )

        self.classifier = nn.Sequential(
            nn.LayerNorm(common_dim),
            nn.Linear(common_dim, common_dim // 2),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(common_dim // 2, n_classes),
        )

    def forward(self, embeddings_dict, return_gates=False):
        batch_size = next(iter(embeddings_dict.values())).shape[0]

        # Concatenate all embeddings for gating
        all_embs = torch.cat([embeddings_dict[name] for name in self.model_names], dim=1)
        gates = self.gate(all_embs)  # [batch, n_models]

        # Project and gate
        projected = []
        for i, name in enumerate(self.model_names):
            proj = self.projections[name](embeddings_dict[name])  # [batch, common_dim]
            gated = proj * gates[:, i:i+1]
            projected.append(gated)

        fused = sum(projected)  # [batch, common_dim]
        logits = self.classifier(fused)

        if return_gates:
            return logits, fused, gates

        return logits

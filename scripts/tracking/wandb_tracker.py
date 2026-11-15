"""
Weights & Biases experiment tracking for scfm-benchmark.

Provides a unified interface for logging benchmark results, model comparisons,
and attention fusion training to W&B.

Usage:
    from scripts.tracking.wandb_tracker import BenchmarkTracker

    tracker = BenchmarkTracker(project="scfm-benchmark", experiment="attention_fusion")
    tracker.log_model_result("STATE", auc=0.871, balanced_acc=0.85, ci_low=0.869, ci_high=0.875)
    tracker.finish()
"""

from typing import Optional

try:
    import wandb

    HAS_WANDB = True
except ImportError:
    HAS_WANDB = False

import numpy as np


class BenchmarkTracker:
    """W&B tracker for foundation model benchmark experiments."""

    def __init__(
        self,
        project: str = "scfm-benchmark",
        experiment: str = "benchmark",
        config: Optional[dict] = None,
        tags: Optional[list] = None,
        offline: bool = False,
    ):
        self.project = project
        self.experiment = experiment
        self.run = None

        if not HAS_WANDB:
            print("[WARNING] wandb not installed. Tracking disabled. pip install wandb")
            return

        mode = "offline" if offline else "online"
        self.run = wandb.init(
            project=project,
            name=experiment,
            config=config or {},
            tags=tags or [],
            mode=mode,
        )

    @property
    def active(self):
        return self.run is not None

    def log_model_result(self, model_name: str, **metrics):
        """Log single model evaluation metrics."""
        if not self.active:
            return
        prefixed = {f"{model_name}/{k}": v for k, v in metrics.items()}
        wandb.log(prefixed)

    def log_comparison_table(self, results: dict):
        """Log model comparison as W&B Table.

        Args:
            results: dict of {model_name: {metric: value, ...}}
        """
        if not self.active:
            return

        columns = ["model", "auc", "balanced_accuracy", "ci_low", "ci_high"]
        data = []
        for model, metrics in results.items():
            row = [
                model,
                metrics.get("auc", 0),
                metrics.get("balanced_accuracy", 0),
                metrics.get("ci_low", 0),
                metrics.get("ci_high", 0),
            ]
            data.append(row)

        table = wandb.Table(columns=columns, data=data)
        wandb.log({"model_comparison": table})

    def log_fusion_epoch(self, epoch: int, fold: int, train_loss: float, val_auc: float):
        """Log attention fusion training progress."""
        if not self.active:
            return
        wandb.log(
            {
                "epoch": epoch,
                "fold": fold,
                "train/loss": train_loss,
                "val/auc": val_auc,
            }
        )

    def log_attention_weights(self, attention_dict: dict, severity_labels: np.ndarray):
        """Log attention weight distribution as W&B bar chart."""
        if not self.active:
            return

        # Mean attention per model
        for name, weights in attention_dict.items():
            wandb.log({f"attention/{name}_mean": float(np.mean(weights))})

        # Per-severity breakdown
        for sev_val, sev_name in [(0, "non_severe"), (1, "severe")]:
            mask = severity_labels == sev_val
            for name, weights in attention_dict.items():
                wandb.log({f"attention/{sev_name}/{name}": float(np.mean(weights[mask]))})

    def log_robustness_curve(self, axis_name: str, x_values: list, model_results: dict):
        """Log robustness degradation curves (sample efficiency, gene dropout, etc.).

        Args:
            axis_name: e.g. "dropout_pct" or "training_fraction"
            x_values: list of x-axis values
            model_results: {model_name: [auc_at_x0, auc_at_x1, ...]}
        """
        if not self.active:
            return

        for model, aucs in model_results.items():
            for x, auc in zip(x_values, aucs):
                wandb.log(
                    {
                        f"robustness/{axis_name}": x,
                        f"robustness/{model}_auc": auc,
                    }
                )

    def log_artifact(self, name: str, path: str, artifact_type: str = "result"):
        """Log a file as W&B artifact."""
        if not self.active:
            return
        artifact = wandb.Artifact(name=name, type=artifact_type)
        artifact.add_file(path)
        self.run.log_artifact(artifact)

    def log_summary(self, summary: dict):
        """Log final summary metrics."""
        if not self.active:
            return
        for k, v in summary.items():
            wandb.run.summary[k] = v

    def finish(self):
        """Finish the W&B run."""
        if self.active:
            wandb.finish()
            self.run = None


class FusionTracker(BenchmarkTracker):
    """Specialized tracker for multimodal fusion experiments."""

    def __init__(self, n_models: int = 4, common_dim: int = 256, **kwargs):
        config = {
            "n_models": n_models,
            "common_dim": common_dim,
            "architecture": "cross_attention",
            "n_heads": 4,
            "epochs": 100,
            "batch_size": 256,
            "optimizer": "AdamW",
            "lr": 1e-3,
        }
        config.update(kwargs.get("config", {}))
        kwargs["config"] = config
        kwargs.setdefault("tags", []).extend(["fusion", f"{n_models}-model"])
        super().__init__(**kwargs)

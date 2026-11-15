#!/usr/bin/env python3
"""
Backfill existing benchmark results into W&B.

Reads completed result JSON/CSV files and logs them as historical W&B runs.
Run this once after creating your W&B account to populate the dashboard.

Usage:
    # First: wandb login
    python scripts/tracking/backfill_experiments.py --results-dir /path/to/results

    # Dry run (offline mode):
    python scripts/tracking/backfill_experiments.py --offline
"""

import argparse
import json
from pathlib import Path

import pandas as pd

try:
    import wandb
except ImportError:
    raise ImportError("pip install wandb")


WANDB_PROJECT = "scfm-benchmark"


def backfill_severity_benchmark(results_dir: Path, offline: bool = False):
    """Backfill COVID severity benchmark results."""
    # Look for benchmark result files
    patterns = [
        "benchmark_results_*.json",
        "covid_severity_benchmark_*.json",
        "fusion_results_*.json",
    ]

    for pattern in patterns:
        for f in sorted(results_dir.glob(f"**/{pattern}")):
            print(f"  Processing: {f.name}")
            with open(f) as fh:
                data = json.load(fh)

            wandb.init(
                project=WANDB_PROJECT,
                name=f"backfill/{f.stem}",
                config={"source_file": str(f)},
                tags=["backfill", "benchmark"],
                mode="offline" if offline else "online",
            )

            # Log each model's metrics
            for key, val in data.items():
                if isinstance(val, dict) and "auc" in val:
                    wandb.log({f"{key}/auc": val["auc"]})
                    if "balanced_accuracy" in val:
                        wandb.log({f"{key}/balanced_accuracy": val["balanced_accuracy"]})

            wandb.finish()


def backfill_bootstrap_ci(results_dir: Path, offline: bool = False):
    """Backfill statistical significance (bootstrap CI) results."""
    ci_files = list(results_dir.glob("**/bootstrap_ci*.csv")) + list(
        results_dir.glob("**/statistical_significance*.json")
    )

    for f in ci_files:
        print(f"  Processing: {f.name}")
        wandb.init(
            project=WANDB_PROJECT,
            name="backfill/bootstrap_ci",
            tags=["backfill", "statistics", "bootstrap"],
            mode="offline" if offline else "online",
        )

        if f.suffix == ".csv":
            df = pd.read_csv(f)
            for _, row in df.iterrows():
                model = row.get("model", row.get("Model", "unknown"))
                wandb.log(
                    {
                        f"{model}/auc": row.get("auc", row.get("AUC", 0)),
                        f"{model}/ci_low": row.get("ci_low", row.get("CI_low", 0)),
                        f"{model}/ci_high": row.get("ci_high", row.get("CI_high", 0)),
                    }
                )
        elif f.suffix == ".json":
            with open(f) as fh:
                data = json.load(fh)
            for model, metrics in data.items():
                if isinstance(metrics, dict):
                    for k, v in metrics.items():
                        if isinstance(v, (int, float)):
                            wandb.log({f"{model}/{k}": v})

        wandb.finish()


def backfill_robustness(results_dir: Path, offline: bool = False):
    """Backfill robustness results (batch shift, gene dropout, sample efficiency)."""
    robustness_files = (
        list(results_dir.glob("**/batch_shift*.json"))
        + list(results_dir.glob("**/gene_dropout*.json"))
        + list(results_dir.glob("**/gene_dropout*.csv"))
        + list(results_dir.glob("**/sample_efficiency*.json"))
    )

    for f in robustness_files:
        analysis = f.stem.split("_")[0:2]
        tag = "_".join(analysis)
        print(f"  Processing: {f.name} (tag={tag})")

        wandb.init(
            project=WANDB_PROJECT,
            name=f"backfill/{tag}",
            tags=["backfill", "robustness", tag],
            mode="offline" if offline else "online",
        )

        if f.suffix == ".json":
            with open(f) as fh:
                data = json.load(fh)

            # Flatten and log
            def _flatten(d, prefix=""):
                for k, v in d.items():
                    key = f"{prefix}/{k}" if prefix else k
                    if isinstance(v, (int, float)):
                        wandb.log({key: v})
                    elif isinstance(v, dict):
                        _flatten(v, key)

            _flatten(data)
        elif f.suffix == ".csv":
            df = pd.read_csv(f)
            # Log as table
            table = wandb.Table(dataframe=df)
            wandb.log({tag: table})

        wandb.finish()


def backfill_fusion(results_dir: Path, offline: bool = False):
    """Backfill multimodal fusion results."""
    fusion_files = list(results_dir.glob("**/fusion_*.json"))

    for f in fusion_files:
        print(f"  Processing: {f.name}")
        with open(f) as fh:
            data = json.load(fh)

        n_models = data.get("n_models", 4)
        wandb.init(
            project=WANDB_PROJECT,
            name=f"backfill/fusion_{n_models}model",
            config={"n_models": n_models},
            tags=["backfill", "fusion", f"{n_models}-model"],
            mode="offline" if offline else "online",
        )

        for key, val in data.items():
            if isinstance(val, dict) and "auc" in val:
                wandb.log({f"{key}/auc": val["auc"]})

        wandb.finish()


def main():
    parser = argparse.ArgumentParser(description="Backfill W&B experiments")
    parser.add_argument(
        "--results-dir", type=str, default="/bigdata/godziklab/shared/Xinru/302006/05_RESULTS"
    )
    parser.add_argument(
        "--offline", action="store_true", help="Run in offline mode (no W&B account needed)"
    )
    args = parser.parse_args()

    results_dir = Path(args.results_dir)
    mode_str = "OFFLINE" if args.offline else "ONLINE"
    print(f"Backfilling experiments from: {results_dir}")
    print(f"Mode: {mode_str}")
    print("=" * 60)

    print("\n[1/4] Severity benchmark results...")
    backfill_severity_benchmark(results_dir, args.offline)

    print("\n[2/4] Bootstrap CI results...")
    backfill_bootstrap_ci(results_dir, args.offline)

    print("\n[3/4] Robustness results...")
    backfill_robustness(results_dir, args.offline)

    print("\n[4/4] Fusion results...")
    backfill_fusion(results_dir, args.offline)

    print("\n" + "=" * 60)
    print("Done. Check your W&B dashboard at: https://wandb.ai")
    if args.offline:
        print("Offline runs stored locally. Sync with: wandb sync")


if __name__ == "__main__":
    main()

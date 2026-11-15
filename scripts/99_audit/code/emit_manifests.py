#!/usr/bin/env python3
"""Emit (or refresh) ``results/manifest.json`` for one or more pipeline stages.

Audit item **A001**: no stage emits a `results/manifest.json` (architecture §9).
This generator closes that gap. It is the canonical way to produce a stage manifest
going forward — point it at a stage and it walks ``results/``, content-hashes every
output, records the git commit, and merges in the stage's declared metadata
(SLURM jobs, parameters, the claims it backs, the canonical evidence files).

IMPORTANT — provenance honesty
------------------------------
The stages this repo ships were run on the HPCC weeks before the repo was
restructured; the SLURM logs, run-time package versions and wall-clock timings
live on the cluster, **not in this checkout**. So a manifest produced here is a
*post-hoc reconstruction*, and it says so: ``provenance.reconstructed_post_hoc =
true`` and ``runtime_seconds``/``package_versions`` are recorded as ``null`` with a
note rather than back-filled with values we cannot verify. Sizes + SHA-256 of the
committed outputs and the git commit ARE verifiable, and those anchor the manifest.

When a stage is *re-run* in-repo, write the manifest from inside the stage instead
(architecture §7 `write_manifest`) so these fields are captured for real.

Usage
-----
    python scripts/99_audit/code/emit_manifests.py --stage 17_domain_adapted_fm
    python scripts/99_audit/code/emit_manifests.py --all          # every stage in the registry
    python scripts/99_audit/code/emit_manifests.py --all --check  # don't write; report drift
"""
from __future__ import annotations

import argparse
import datetime as _dt
import hashlib
import json
import subprocess
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]

# --- per-stage declared metadata -------------------------------------------------
# Extend this registry one entry per stage to clear A001 across all 30 stages.
# `claims` / `primary_outputs` tie the manifest to docs/claim-traceability.md so the
# manifest doubles as the canonical evidence-file pin.
STAGE_META: dict[str, dict] = {
    "17_domain_adapted_fm": {
        "description": "Geneformer LoRA domain-adaptation vs zero-shot on the platelet "
                       "mortality task, under 5-fold patient-level StratifiedGroupKFold CV "
                       "(no train/test patient overlap).",
        "slurm_jobs": ["23853347-23866507"],
        "run_dates": ["2026-05-06"],
        "generating_code": [
            # The verified 0.705 patient-CV LoRA result is produced by the finetuning
            # PEFT script + its submit wrapper — NOT the stage-17 continued-pretraining
            # scripts (geneformer_continue_pretrain.py reproduces the superseded 0.879
            # overlap result that lives in the 7model_benchmark_*.csv files below).
            "scripts/finetuning/geneformer_peft_finetune.py",
            "scripts/finetuning/submit_geneformer_finetune_mortality.sh",
        ],
        "parameters": {
            "backbone": "Geneformer-V1-10M",
            "peft": "LoRA", "lora_r": 16, "lora_alpha": 32, "lora_dropout": 0.1,
            "epochs": 25, "learning_rate": 8e-4,
            "cv": "5-fold StratifiedGroupKFold (patient-level, 339 patients)",
        },
        "claims": [
            {"id": "C003", "table": "tab:domainadapt",
             "value": "Geneformer LoRA mean test AUC 0.705 +/- 0.135",
             "evidence": "results/geneformer_peft/results_mortality_lora.json"},
        ],
        "primary_outputs": ["geneformer_peft/results_mortality_lora.json"],
        "superseded_outputs": {
            "7model_benchmark_20260309_222122.csv":
                "legacy continued-pretraining/overlap run (Geneformer_adapted binary LR=0.879); "
                "NOT the patient-level-CV LoRA result — does not back C003.",
            "7model_benchmark_20260311_103812.csv":
                "legacy continued-pretraining/overlap run (Geneformer_adapted binary LR=0.879); "
                "NOT the patient-level-CV LoRA result — does not back C003.",
        },
        "inputs": [
            "data/processed/covid_sepsis_mortality.h5ad",
        ],
    },
    "24_option_c_canonical_input": {
        "description": "Option C standardized within-disease mortality benchmark on the "
                       "disease-matched COVID+Sepsis subset (42,095 cells / 322 donors): "
                       "primary LR/RF AUC, DeLong pairwise, calibration, plus auxiliary "
                       "ablation/GSEA/SHAP/activation analyses.",
        "slurm_jobs": ["23921397", "23823131", "23823133", "23851383"],
        "run_dates": ["2026-05-09", "2026-05-04", "2026-05-05"],
        "generating_code": [
            "scripts/24_option_c_canonical_input/code/primary_benchmark_option_c.py",
            "scripts/24_option_c_canonical_input/code/calibration_per_cohort.py",
        ],
        "parameters": {
            "cohort": "COVID+Sepsis subset (disease-matched)",
            "n_cells": 42095, "n_donors": 322, "n_fatal": 8850, "n_survived": 33245,
            "cv": "5-fold patient-level", "classifiers": ["LR (L2, balanced)", "RandomForest"],
            "ece_bins": 10,
        },
        "claims": [
            {"id": "C004", "table": "tab:performance",
             "value": "scPRINT-2 0.762 / Geneformer 0.738 / STATE 0.724 / UCE 0.598 (per-fold mean LR AUC)",
             "evidence": "results/option_c_primary_benchmark/covid_sepsis_only_mortality.csv"},
            {"id": "C004", "table": "tab:delong",
             "value": "STATE 0.735 > scPRINT-2 0.679 > Geneformer 0.670 > UCE 0.608 (OOF cell-weighted)",
             "evidence": "results/option_c_primary_benchmark/covid_sepsis_only/delong_pairwise.csv"},
            {"id": "C004", "table": "tab:calibration",
             "value": "STATE-LR AUC 0.732 / Brier 0.179 / ECE 0.156; STATE-RF Brier 0.149 / ECE 0.034",
             "evidence": "results/option_c_primary_benchmark/covid_sepsis_only/calibration.csv"},
        ],
        "primary_outputs": [
            "option_c_primary_benchmark/covid_sepsis_only_mortality.csv",
            "option_c_primary_benchmark/covid_sepsis_only/delong_pairwise.csv",
            "option_c_primary_benchmark/covid_sepsis_only/calibration.csv",
        ],
        "superseded_outputs": {
            "option_c_primary_benchmark/mortality_summary.csv":
                "non-disease-matched run (STATE-LR 0.712); main.tex L317 excludes from headline tables.",
            "option_c_primary_benchmark/delong_pairwise.csv":
                "non-disease-matched run (STATE 0.720); headline DeLong is covid_sepsis_only/delong_pairwise.csv.",
            "option_c_calibration/calibration_summary.csv":
                "non-disease-matched run (STATE 0.719); headline calibration is covid_sepsis_only/calibration.csv.",
        },
        "inputs": [
            "data/processed/covid_sepsis_mortality.h5ad",
        ],
    },
    "29_adaptation_curve": {
        "description": "Cross-disease adaptation curve: AUC on a fixed held-out sepsis test "
                       "split vs number of labelled sepsis cells added to the constant COVID "
                       "training pool, 10 seeds/point, 200-iteration bootstrap CIs.",
        "slurm_jobs": ["24018607"],
        "run_dates": ["2026-05-13"],
        "generating_code": [
            "scripts/29_adaptation_curve/code/run_adaptation_curve_bootstrap.py",
        ],
        "parameters": {
            "covid_pool_cells": 38673,
            "n_sepsis_grid": [0, 10, 50, 100, 500, 1000, "all (2882)"],
            "seeds_per_point": 10, "bootstrap_iters": 200,
            "classifier": "LR (L2, balanced) on standardized embeddings; PCA-200 for from-scratch baseline",
            "headline_column": "point_auc_mean (NOT boot_auc_mean)",
        },
        "claims": [
            {"id": "C005", "table": "tab:adaptation_curve",
             "value": "PCA-200+LR asymptote 0.964; scPRINT-2 Y-int 0.920; per-method Y-int/n100/n1000/asymptote",
             "evidence": "results/adaptation_curve_bootstrap/summary.csv"},
        ],
        "primary_outputs": ["adaptation_curve_bootstrap/summary.csv"],
        "superseded_outputs": {
            "adaptation_curve/summary.csv":
                "single-seed legacy variant; the manuscript uses the 200-iter bootstrap "
                "(adaptation_curve_bootstrap/summary.csv).",
        },
        "inputs": [
            "data/processed/covid_sepsis_mortality.h5ad",
        ],
    },
    "06_batch_generalization": {
        "description": "Bidirectional cross-disease generalization (fatal vs survived): train on one "
                       "disease, test on the other (COVID<->Sepsis), severe-excluded subset. Backs the "
                       "C002 cross-disease claim.",
        "slurm_jobs": ["23761950"],
        "run_dates": ["2026-05-04", "2026-04-30"],
        "generating_code": [
            "scripts/06_batch_generalization/code/cross_disease_fatal_vs_survived.py",
        ],
        "parameters": {
            "cohort": "severe-EXCLUDED cross-disease (n_train=21286)",
            "directions": "forward (COVID->Sepsis) + reverse (Sepsis->COVID)",
            "classifier": "LogisticRegression (auc)",
        },
        "claims": [
            {"id": "C002", "table": "tab:crossdisease",
             "value": "STATE fwd 0.798 / rev 0.707 (only primary FM clearing both directions)",
             "evidence": "results/cross_disease_fatal_vs_survived/cross_disease_fatal_vs_survived_20260504_152836.csv"},
        ],
        "primary_outputs": [
            "cross_disease_fatal_vs_survived/cross_disease_fatal_vs_survived_20260504_152836.csv",
            "cross_disease_fatal_vs_survived/cross_disease_fatal_vs_survived_20260430_193915.csv",
        ],
        "superseded_outputs": {
            "cross_disease_fatal_vs_survived/cross_disease_fatal_vs_survived_20260522_151856.csv":
                "severe-INCLUDED rerun (n_train=38673) — not the manuscript source (A006 false alarm).",
        },
        "inputs": [
            "data/processed/covid_sepsis_mortality.h5ad",
        ],
    },
}

GENERATED_BY = "scripts/99_audit/code/emit_manifests.py"

# Numbered stages whose results are NOT committed in-repo — a stub manifest documents
# where they actually live, so A001 is closed for the whole numbered pipeline.
NO_RESULTS_REASON: dict[str, str] = {
    "04_pca_baseline": "Results on HPCC (05_RESULTS) / folded into stage 24 option_c PCA baseline; not committed.",
    "07_gene_dropout": "Results on HPCC (05_RESULTS, gene-dropout robustness); not committed (backs tab:dropout).",
    "19_model_optimization": "Quantization/FSDP/torch.compile benchmarks — HPCC GPU runs; results not committed.",
    "22_esm_protein_multimodal": "ESM protein-multimodal — HPCC GPU run; results not committed.",
    "26_figures": "Generates manuscript figures (PDF/PNG, gitignored); no CSV results — see figure-table-inventory.md.",
    "27_covid_sepsis_only_rerun": "Writes to HPCC 05_RESULTS/.../covid_sepsis_only/ (committed copy lives under stage 24 results).",
    "28_umap_covid_sepsis": "UMAP coordinates on HPCC (05_RESULTS); not committed (large, regenerable).",
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def git_commit() -> str | None:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=REPO, text=True
        ).strip()
    except Exception:
        return None


def discover_stages() -> list[str]:
    """All numbered stage dirs under scripts/ that have a results/ directory."""
    out = []
    for d in sorted((REPO / "scripts").iterdir()):
        if d.is_dir() and d.name[:1].isdigit() and (d / "results").is_dir():
            out.append(d.name)
    return out


def derive_description(stage: str) -> str:
    """Best-effort description from the stage README's first '# ' heading."""
    rd = REPO / "scripts" / stage / "README.md"
    if rd.is_file():
        for line in rd.read_text(errors="replace").splitlines():
            line = line.strip()
            if line.startswith("# "):
                return f"{line[2:].strip()} (auto-derived from README; metadata not curated)"
    return f"{stage} (no README; metadata not curated)"


def build_manifest(stage: str, meta: dict) -> dict:
    stage_dir = REPO / "scripts" / stage
    results_dir = stage_dir / "results"
    outputs = []
    for f in sorted(results_dir.rglob("*")):
        if f.is_file() and f.name != "manifest.json":
            outputs.append({
                "path": str(f.relative_to(stage_dir)),
                "size_bytes": f.stat().st_size,
                "sha256": sha256(f),
            })

    inputs = []
    for rel in meta.get("inputs", []):
        p = REPO / rel
        inputs.append({
            "path": rel,
            "exists": p.exists(),
            "size_bytes": p.stat().st_size if p.exists() else None,
            "note": None if p.exists() else "gitignored / on HPCC (data/ not committed)",
        })

    metadata_complete = "description" in meta
    description = meta.get("description") or derive_description(stage)

    now = _dt.datetime.now(_dt.timezone.utc).replace(microsecond=0).isoformat()
    return {
        "stage": stage,
        "status": "success",
        "metadata_complete": metadata_complete,
        "description": description,
        "slurm_jobs": meta.get("slurm_jobs", []),
        "run_dates": meta.get("run_dates", []),
        "generating_code": meta.get("generating_code", []),
        "parameters": meta.get("parameters", {}),
        "claims": meta.get("claims", []),
        "primary_outputs": meta.get("primary_outputs", []),
        "superseded_outputs": meta.get("superseded_outputs", {}),
        "inputs": inputs,
        "outputs": outputs,
        "n_outputs": len(outputs),
        "in_repo_outputs": len(outputs) > 0,
        "no_results_reason": (None if outputs else
                              NO_RESULTS_REASON.get(stage,
                              "No committed results in this checkout (HPCC-only / gitignored / utility stage).")),
        "git_commit": git_commit(),
        # run-time environment NOT captured (results predate this checkout; see module docstring)
        "runtime_seconds": None,
        "package_versions": None,
        "config": "configs/hpcc.yaml",
        "environment_ref": "environment.yml",
        "provenance": {
            "reconstructed_post_hoc": True,
            "generated_at": now,
            "generated_by": GENERATED_BY,
            "verified_in": "docs/claim-traceability.md (2026-06-06 traceability audit)",
            "note": "Manifest reconstructed during the 2026-06-06/07 audit (item A001). Outputs were "
                    "produced on the HPCC before the repo restructure; sizes + sha256 + git_commit "
                    "are verifiable, but run-time wall-clock and package versions were not captured "
                    "in-repo and are recorded as null rather than back-filled. When "
                    "metadata_complete=false, slurm_jobs/parameters/claims were not curated for this "
                    "stage — the output inventory + hashes + git pin are still authoritative.",
        },
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--stage", help="stage dir name under scripts/ (e.g. 17_domain_adapted_fm)")
    ap.add_argument("--all", action="store_true", help="emit for every stage in the registry (rich metadata only)")
    ap.add_argument("--all-stages", action="store_true",
                    help="emit for EVERY stage with a results/ dir — rich if in registry, "
                         "minimal (hashes+git pin+README description) otherwise. Closes A001.")
    ap.add_argument("--check", action="store_true",
                    help="do not write; report whether the on-disk manifest would change")
    args = ap.parse_args()

    if args.all_stages:
        stages = discover_stages()
    elif args.all:
        stages = list(STAGE_META)
    elif args.stage:
        if args.stage not in STAGE_META and not (REPO / "scripts" / args.stage / "results").is_dir():
            ap.error(f"stage {args.stage!r} has no registry entry and no results/ dir")
        stages = [args.stage]
    else:
        ap.error("pass --stage <name>, --all, or --all-stages")

    for stage in stages:
        manifest = build_manifest(stage, STAGE_META.get(stage, {}))
        out = REPO / "scripts" / stage / "results" / "manifest.json"
        text = json.dumps(manifest, indent=2) + "\n"
        if args.check:
            prev = out.read_text() if out.exists() else ""
            # ignore the volatile generated_at line when diffing
            def strip_ts(s: str) -> str:
                return "\n".join(l for l in s.splitlines() if '"generated_at"' not in l)
            status = "unchanged" if strip_ts(prev) == strip_ts(text) else "WOULD CHANGE"
            print(f"[check] {stage}: {status} ({manifest['n_outputs']} outputs)")
        else:
            out.write_text(text)
            tag = "rich " if manifest["metadata_complete"] else "min  "
            print(f"[write {tag}] {out.relative_to(REPO)}  "
                  f"({manifest['n_outputs']} outputs, {len(manifest['primary_outputs'])} primary)")


if __name__ == "__main__":
    main()

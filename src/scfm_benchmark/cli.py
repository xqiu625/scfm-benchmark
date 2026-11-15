"""scfm-bench CLI: portable demo/benchmark plus historical research dispatchers.

Thin dispatcher that wraps the numbered scripts/ stage pipelines.
Heavy work (FM embedding generation, large evaluations) is intended to run
on HPCC where the per-FM conda envs and model weights live; the CLI here
shells out to the appropriate scripts so users get a uniform entry point.

Repair notes (2026-08-06):
  - ``evaluate`` stages the single input h5ad into a temp dir as
    ``<stem>_covid_canonical.h5ad`` so it matches the
    ``*_covid_canonical.h5ad`` glob of primary_benchmark_option_c.py, and
    remaps its flags to that script's real CLI (--embedding-dir/--output-dir/
    --models/--analyses). Only the mortality task is wired up; the other
    tasks exit with a clear pointer to the owning script.
  - ``embed`` maps the primary 5 FMs to the submitters that actually exist
    in this clone and fails loudly for anything else.
  - ``diagnose`` points at run_adaptation_curve_bootstrap.py (the only flag
    it accepts is --models); GSEA from a single embedding is not supported
    in this clone and exits with a clear error; SHAP runs only when
    --raw-input is supplied.
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

from scfm_benchmark import __version__

REPO_ROOT = Path(__file__).resolve().parents[2]
OPTION_C = REPO_ROOT / "scripts" / "24_option_c_canonical_input" / "code"
ADAPT_CURVE = REPO_ROOT / "scripts" / "29_adaptation_curve" / "code"
FM_SCRIPTS = REPO_ROOT / "scripts" / "foundation_models"

# Primary 5-FM manuscript set -> (canonical-COVID submitter, how to launch it).
# Verified against the dirs/scripts present in this clone (2026-08-06):
#   STATE      : no per-FM dir exists; step2c_human_covid_canonical_input.sh is
#                the batch driver that submits the STATE canonical job (run it
#                with bash — it calls sbatch internally).
#   UCE        : foundation_models/UCE/uce_human_platelet.sh
#   Geneformer : foundation_models/Geneformer/run_geneformer_covid_option_c.sh
#   scPRINT-2  : foundation_models/scPRINT/run_scprint2_human_covid.sh
#   scGPT      : foundation_models/scGPT/submit_scgpt_option_c.sh
EMBED_SUBMITTERS = {
    "state": (OPTION_C / "step2c_human_covid_canonical_input.sh", "bash"),
    "uce": (FM_SCRIPTS / "UCE" / "uce_human_platelet.sh", "sbatch"),
    "geneformer": (FM_SCRIPTS / "Geneformer" / "run_geneformer_covid_option_c.sh", "sbatch"),
    "scprint2": (FM_SCRIPTS / "scPRINT" / "run_scprint2_human_covid.sh", "sbatch"),
    "scgpt": (FM_SCRIPTS / "scGPT" / "submit_scgpt_option_c.sh", "sbatch"),
}
VALID_MODELS = tuple(EMBED_SUBMITTERS)
VALID_TASKS = ("mortality", "tissue", "cross-disease")

# stem prefix (lowercase) -> display name used by run_adaptation_curve_bootstrap.py
_ADAPT_DISPLAY = {
    "state": "STATE",
    "uce": "UCE",
    "geneformer": "Geneformer",
    "scprint2": "scPRINT-2",
    "scgpt": "scGPT",
}


def _run(cmd: list[str]) -> int:
    print(f"$ {' '.join(cmd)}", file=sys.stderr)
    return subprocess.run(cmd, check=False).returncode


def _stage_embedding(emb: Path, staging_dir: Path) -> Path:
    """Link ``emb`` into ``staging_dir`` under a name matching the
    ``*_covid_canonical.h5ad`` glob used by the option-C scripts.

    The model key is derived downstream from the first ``_``-token of the
    stem, so ``<model>_covid_canonical.h5ad`` (or an already-canonical name)
    resolves to a registered OBSM key; anything else fails loudly inside the
    benchmark script's find_embedding.
    """
    name = emb.name if emb.name.endswith("_covid_canonical.h5ad") else f"{emb.stem}_covid_canonical.h5ad"
    link = staging_dir / name
    try:
        link.symlink_to(emb.resolve())
    except OSError:
        shutil.copy2(emb, link)
    return link


def cmd_embed(args: argparse.Namespace) -> int:
    """Submit the FM-specific embedding job to SLURM.

    Each FM has its own env / weights / CLI; we don't unify them in Python.
    Only the primary 5 FMs have canonical-COVID submitters in this clone.
    """
    entry = EMBED_SUBMITTERS.get(args.model)
    if entry is None:
        print(
            f"ERROR: no canonical-COVID submitter wired up for model={args.model!r}.\n"
            f"Supported models: {', '.join(VALID_MODELS)}.\n"
            f"Other FM code lives under {FM_SCRIPTS} but has no verified\n"
            f"canonical-COVID submitter in this clone; submit those manually.",
            file=sys.stderr,
        )
        return 2

    submit, launcher = entry
    if not submit.exists():
        print(
            f"ERROR: submitter for model={args.model!r} not found at {submit}.\n"
            f"This clone may be incomplete; see scripts/foundation_models/.",
            file=sys.stderr,
        )
        return 1

    print(
        f"# {submit.name} has HPCC-pinned paths; verify inside before submitting:\n"
        f"#   input : {args.h5ad}\n"
        f"#   output: {args.out}\n",
        file=sys.stderr,
    )
    return _run([launcher, str(submit)])


def cmd_evaluate(args: argparse.Namespace) -> int:
    """Run the primary patient-level benchmark on a pre-computed embedding."""
    if args.task != "mortality":
        owners = {
            "cross-disease": "scripts/06_batch_generalization/code/cross_disease_fatal_vs_survived.py",
            "tissue": "(no tissue-task benchmark exists in this clone)",
        }
        print(
            f"ERROR: --task {args.task} is not supported by `scfm-bench evaluate` in this clone.\n"
            f"Only the mortality task is wired up (via primary_benchmark_option_c.py).\n"
            f"See: {owners[args.task]}",
            file=sys.stderr,
        )
        return 2
    if args.groups_col != "Sample_ID":
        print(
            f"ERROR: --groups-col {args.groups_col!r} is not supported: "
            f"primary_benchmark_option_c.py hardcodes Sample_ID patient groups in this clone.",
            file=sys.stderr,
        )
        return 2

    script = OPTION_C / "primary_benchmark_option_c.py"
    if not script.exists():
        print(f"ERROR: benchmark script not found: {script}", file=sys.stderr)
        return 1

    with tempfile.TemporaryDirectory(prefix="scfm_eval_") as td:
        staged = _stage_embedding(args.emb, Path(td))
        model_key = staged.stem.split("_")[0].lower()
        return _run(
            [
                sys.executable,
                str(script),
                "--embedding-dir",
                td,
                "--output-dir",
                str(args.out),
                "--models",
                model_key,
                "--analyses",
                "mortality",
            ]
        )


def cmd_diagnose(args: argparse.Namespace) -> int:
    """scFM Diagnostic Suite: adaptation curve + GSEA + SHAP attribution.

    Only parts that can be driven correctly from this clone are executed:
      - adaptation curve: run_adaptation_curve_bootstrap.py --models <FM>.
        NOTE: that script has HPCC-pinned ROOT/OUT_DIR; results land in its
        own OUT_DIR, not in --out.
      - SHAP: per_fm_shap_attribution.py, only when --raw-input is given.
      - GSEA: NOT supported from a single embedding in this clone (it
        consumes per-FM gene-ranking CSVs, not an embedding) — a clear
        error is printed and the subcommand exits non-zero.
    """
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    model_key = args.emb.stem.split("_")[0].lower()
    failed = False

    # 1. Adaptation curve (bootstrap-corrected)
    adapt = ADAPT_CURVE / "run_adaptation_curve_bootstrap.py"
    if not adapt.exists():
        print(f"ERROR: adaptation-curve script not found: {adapt}", file=sys.stderr)
        failed = True
    elif model_key not in _ADAPT_DISPLAY:
        print(
            f"ERROR: cannot map {args.emb.name!r} to an FM known to the adaptation-curve "
            f"script ({', '.join(sorted(_ADAPT_DISPLAY))}). Rename the embedding to "
            f"<model>_....h5ad or run {adapt} directly.",
            file=sys.stderr,
        )
        failed = True
    else:
        print(
            "# NOTE: run_adaptation_curve_bootstrap.py uses HPCC-pinned ROOT/OUT_DIR;\n"
            "# results are written to its own OUT_DIR, not to --out.",
            file=sys.stderr,
        )
        failed |= _run([sys.executable, str(adapt), "--models", _ADAPT_DISPLAY[model_key]]) != 0

    # 2. GSEA per-dim — consumes top_genes_<fm>_dim<d>.csv rankings produced by
    #    biological_interpretation_option_c.py; cannot run on a bare embedding.
    print(
        "ERROR: the GSEA view is not supported by `scfm-bench diagnose` in this clone:\n"
        "gsea_per_fm_dim.py needs per-FM gene-ranking CSVs, not an embedding.\n"
        "Run scripts/24_option_c_canonical_input/code/biological_interpretation_option_c.py\n"
        "first, then scripts/24_option_c_canonical_input/code/run_gsea_per_fm_dim.sh.",
        file=sys.stderr,
    )
    failed = True

    # 3. SHAP attribution — needs the raw-counts h5ad for dim->gene mapping.
    shap_script = OPTION_C / "per_fm_shap_attribution.py"
    if not shap_script.exists():
        print(f"ERROR: SHAP script not found: {shap_script}", file=sys.stderr)
        failed = True
    elif args.raw_input is None:
        print(
            "ERROR: the SHAP view requires the raw-counts h5ad for Cell_ID-aligned\n"
            "dim->gene mapping; re-run with --raw-input <raw.h5ad>, or see\n"
            "scripts/24_option_c_canonical_input/code/run_shap_attribution.sh.",
            file=sys.stderr,
        )
        failed = True
    else:
        with tempfile.TemporaryDirectory(prefix="scfm_shap_") as td:
            _stage_embedding(args.emb, Path(td))
            failed |= _run(
                [
                    sys.executable,
                    str(shap_script),
                    "--embedding-dir",
                    td,
                    "--raw-input",
                    str(args.raw_input),
                    "--output-dir",
                    str(out / "shap"),
                    "--fms",
                    model_key,
                ]
            ) != 0

    if failed:
        print("diagnose: one or more views did not run (see messages above).", file=sys.stderr)
        return 2
    return 0


def cmd_report(args: argparse.Namespace) -> int:
    """Bundle results_dir artifacts into a single HTML report."""
    from scfm_benchmark.api import build_report

    return build_report(results_dir=Path(args.results_dir), out_path=Path(args.out))


def cmd_demo(args: argparse.Namespace) -> int:
    from scfm_benchmark.demo import run_demo

    run_demo(args.out)
    print(f"Synthetic demo complete: {args.out / 'results' / 'report.html'}")
    return 0


def cmd_benchmark(args: argparse.Namespace) -> int:
    from scfm_benchmark.benchmark import benchmark_from_csv

    def assignments(values, convert):
        result = {}
        for value in values:
            name, sep, item = value.partition("=")
            if not sep or not name or not item or name in result:
                raise ValueError("use unique NAME=VALUE assignments")
            result[name] = convert(item)
        return result

    benchmark_from_csv(args.metadata, assignments(args.embedding, Path), args.out,
                       n_splits=args.folds, seeds=args.seeds, pca=assignments(args.pca, int))
    print(f"Benchmark complete: {args.out / 'report.html'}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="scfm-bench",
        description="Single-cell foundation model robustness benchmark.",
    )
    p.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    sub = p.add_subparsers(dest="command", required=True)

    demo = sub.add_parser("demo", help="Offline synthetic CPU demo; no model downloads")
    demo.add_argument("--out", required=True, type=Path, help="New output directory")
    demo.set_defaults(func=cmd_demo)

    bench = sub.add_parser("benchmark", help="Portable donor-aware evaluation of embedding CSVs")
    bench.add_argument("--metadata", required=True, type=Path)
    bench.add_argument("--embedding", required=True, action="append", metavar="NAME=CSV")
    bench.add_argument("--pca", action="append", default=[], metavar="NAME=DIMENSIONS")
    bench.add_argument("--folds", type=int, default=5)
    bench.add_argument("--seeds", type=int, nargs="+", default=[42, 43, 44])
    bench.add_argument("--out", required=True, type=Path, help="New output directory")
    bench.set_defaults(func=cmd_benchmark)

    e = sub.add_parser("embed", help="Submit FM embedding job to SLURM (HPCC)")
    e.add_argument("--h5ad", required=True, type=Path, help="Input raw-counts h5ad")
    e.add_argument("--model", required=True,
                   help=f"FM to embed with. Supported: {', '.join(VALID_MODELS)}")
    e.add_argument("--out", required=True, type=Path, help="Output embedding h5ad")
    e.set_defaults(func=cmd_embed)

    v = sub.add_parser("evaluate", help="Patient-level mortality benchmark on an embedding")
    v.add_argument("--emb", required=True, type=Path)
    v.add_argument("--task", required=True, choices=VALID_TASKS)
    v.add_argument("--out", required=True, type=Path)
    v.add_argument("--groups-col", default="Sample_ID")
    v.set_defaults(func=cmd_evaluate)

    d = sub.add_parser("diagnose", help="scFM Diagnostic Suite (adapt + GSEA + SHAP)")
    d.add_argument("--emb", required=True, type=Path)
    d.add_argument("--out", required=True, type=Path)
    d.add_argument("--raw-input", type=Path, default=None,
                   help="Raw-counts h5ad required for the SHAP view")
    d.set_defaults(func=cmd_diagnose)

    r = sub.add_parser("report", help="Bundle results into a single HTML report")
    r.add_argument("--results-dir", required=True, type=Path)
    r.add_argument("--out", required=True, type=Path)
    r.set_defaults(func=cmd_report)

    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.command in {"demo", "benchmark"}:
        try:
            return args.func(args)
        except (ValueError, OSError) as exc:
            print(f"ERROR: {exc}", file=sys.stderr)
            return 2
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())

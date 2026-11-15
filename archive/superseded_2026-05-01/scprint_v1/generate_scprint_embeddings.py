#!/usr/bin/env python3
"""Generate scPRINT embeddings for an AnnData input.

Loads a scPRINT checkpoint, runs the Embedder over the input AnnData, and
writes a new h5ad with the scPRINT cell embeddings stored in `.obsm`.

Usage:
    python generate_scprint_embeddings.py \
        --input  /bigdata/godziklab/shared/Xinru/302006/02_EMBEDDINGS/state/covid_state.h5ad \
        --output /bigdata/godziklab/shared/Xinru/302006/02_EMBEDDINGS/scprint/covid_scprint.h5ad \
        --ckpt   auto \
        --batch-size 64 \
        --max-len 2000

`--ckpt auto` picks the largest *.ckpt under
$HF_HOME/hub/models--jkobject--scPRINT/snapshots/*. Pass an absolute path to
override.

Per project policy: no synthetic / random fallback. Any failure raises.
"""
from __future__ import annotations

import argparse
import glob
import logging
import os
import sys
from pathlib import Path

import numpy as np
import scanpy as sc

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)],
)
log = logging.getLogger("scprint_embed")


def find_default_ckpt() -> str:
    hf_home = os.environ.get("HF_HOME", "/bigdata/godziklab/shared/Xinru/.huggingface")
    pattern = f"{hf_home}/hub/models--jkobject--scPRINT/snapshots/*/*.ckpt"
    ckpts = sorted(glob.glob(pattern), key=os.path.getsize, reverse=True)
    if not ckpts:
        raise FileNotFoundError(
            f"No scPRINT checkpoints found at {pattern}. "
            f"Run `hf download jkobject/scPRINT` first."
        )
    chosen = ckpts[0]
    log.info("Auto-selected checkpoint: %s (%.2f GB)", chosen, os.path.getsize(chosen) / 1e9)
    return chosen


def load_model(ckpt_path: str):
    from scprint import scPrint  # noqa: WPS433  (import inside fn for faster --help)

    log.info("Loading scPRINT from %s", ckpt_path)
    model = scPrint.load_from_checkpoint(
        ckpt_path,
        precpt_gene_emb=None,
        transformer="normal",
        strict_loading=False,
    )
    model.eval()
    return model


def run_embedder(model, adata, batch_size: int, max_len: int):
    from scprint.tasks import Embedder  # noqa: WPS433

    log.info(
        "Running Embedder (batch_size=%d, max_len=%d) on %d cells x %d genes",
        batch_size,
        max_len,
        adata.n_obs,
        adata.n_vars,
    )
    embedder = Embedder(
        how="random expr",
        max_len=max_len,
        add_zero_genes=0,
        pred_embedding=[],
        keep_all_cls_pred=False,
        doplot=False,
        batch_size=batch_size,
        output_expression="none",
    )
    out = embedder(model, adata)
    if isinstance(out, tuple):
        adata_out = out[0]
    else:
        adata_out = out
    return adata_out


def attach_embedding(adata_in, adata_out):
    """Find scPRINT embedding in adata_out.obsm and attach to adata_in as X_scprint."""
    candidate_keys = ("scprint", "X_scprint", "scPRINT", "embedding", "X_embedding")
    found_key = None
    for key in candidate_keys:
        if key in adata_out.obsm:
            found_key = key
            break
    if found_key is None:
        # Fall back: take any obsm with matching n_obs and >0 width
        for key, arr in adata_out.obsm.items():
            if arr.shape[0] == adata_out.n_obs and arr.ndim == 2 and arr.shape[1] > 1:
                found_key = key
                break
    if found_key is None:
        raise RuntimeError(
            f"No scPRINT embedding found in adata_out.obsm. Keys: {list(adata_out.obsm.keys())}"
        )

    emb = np.asarray(adata_out.obsm[found_key])
    if emb.shape[0] != adata_in.n_obs:
        raise RuntimeError(
            f"Embedding row count {emb.shape[0]} != input cells {adata_in.n_obs}. "
            f"scPRINT may have filtered cells; refusing to mismatch obs."
        )
    if np.isnan(emb).any():
        raise RuntimeError("scPRINT embedding contains NaN values — aborting.")

    log.info("scPRINT embedding shape: %s (from key '%s')", emb.shape, found_key)
    adata_in.obsm["X_scprint"] = emb
    return adata_in


def process_chunked(
    model,
    adata,
    chunk_size: int,
    batch_size: int,
    max_len: int,
    checkpoint_dir: Path,
):
    """Process AnnData in chunks of N cells, saving each chunk's embeddings to checkpoint_dir.

    On rerun, completed chunks are skipped. After all chunks done, embeddings
    are reassembled in original cell order.
    """
    checkpoint_dir.mkdir(parents=True, exist_ok=True)
    n_chunks = (adata.n_obs + chunk_size - 1) // chunk_size
    log.info("Chunked mode: %d chunks of %d cells (total %d)", n_chunks, chunk_size, adata.n_obs)

    embedding_dim = None
    for chunk_idx in range(n_chunks):
        chunk_path = checkpoint_dir / f"chunk_{chunk_idx:05d}.npy"
        if chunk_path.exists():
            arr = np.load(chunk_path)
            log.info("Chunk %d/%d already done (%s, shape=%s) — skipping",
                     chunk_idx + 1, n_chunks, chunk_path.name, arr.shape)
            embedding_dim = arr.shape[1]
            continue

        start = chunk_idx * chunk_size
        end = min(start + chunk_size, adata.n_obs)
        log.info("Chunk %d/%d: cells [%d, %d)", chunk_idx + 1, n_chunks, start, end)

        adata_chunk = adata[start:end].copy()
        adata_out = run_embedder(model, adata_chunk, batch_size, max_len)

        emb = None
        for key in ("scprint", "X_scprint", "scPRINT", "embedding", "X_embedding"):
            if key in adata_out.obsm:
                emb = np.asarray(adata_out.obsm[key])
                break
        if emb is None:
            for key, arr in adata_out.obsm.items():
                if arr.shape[0] == adata_out.n_obs and arr.ndim == 2 and arr.shape[1] > 1:
                    emb = np.asarray(arr)
                    break
        if emb is None:
            raise RuntimeError(f"No embedding found in chunk {chunk_idx} obsm: {list(adata_out.obsm.keys())}")
        if emb.shape[0] != (end - start):
            raise RuntimeError(f"Chunk {chunk_idx}: got {emb.shape[0]} embs for {end - start} cells")
        if np.isnan(emb).any():
            raise RuntimeError(f"Chunk {chunk_idx}: NaN in embeddings")

        np.save(chunk_path, emb.astype(np.float32))
        log.info("  Saved chunk %d (shape=%s) to %s", chunk_idx + 1, emb.shape, chunk_path)
        embedding_dim = emb.shape[1]

    # Reassemble
    log.info("All %d chunks complete. Reassembling...", n_chunks)
    if embedding_dim is None:
        raise RuntimeError("No embedding dim recorded — cannot reassemble")
    full = np.zeros((adata.n_obs, embedding_dim), dtype=np.float32)
    for chunk_idx in range(n_chunks):
        start = chunk_idx * chunk_size
        end = min(start + chunk_size, adata.n_obs)
        arr = np.load(checkpoint_dir / f"chunk_{chunk_idx:05d}.npy")
        full[start:end] = arr
    if np.isnan(full).any():
        raise RuntimeError("Reassembled embeddings contain NaN")
    return full


def main() -> int:
    p = argparse.ArgumentParser(description="Generate scPRINT embeddings")
    p.add_argument("--input", required=True, help="Input AnnData h5ad path")
    p.add_argument("--output", required=True, help="Output AnnData h5ad path")
    p.add_argument(
        "--ckpt",
        default="auto",
        help='scPRINT checkpoint path, or "auto" to pick the largest under HF_HOME',
    )
    p.add_argument("--batch-size", type=int, default=64)
    p.add_argument("--max-len", type=int, default=2000)
    p.add_argument(
        "--checkpoint-dir",
        default=None,
        help="If set, process in chunks of --chunk-size cells, saving each chunk's "
             "embeddings here. Resumable across SLURM time limits.",
    )
    p.add_argument("--chunk-size", type=int, default=5000,
                   help="Cells per chunk in --checkpoint-dir mode (default 5000)")
    args = p.parse_args()

    in_path = Path(args.input)
    out_path = Path(args.output)
    if not in_path.exists():
        raise FileNotFoundError(in_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    ckpt = find_default_ckpt() if args.ckpt == "auto" else args.ckpt
    if not Path(ckpt).exists():
        raise FileNotFoundError(f"Checkpoint not found: {ckpt}")

    log.info("Reading %s", in_path)
    adata = sc.read_h5ad(in_path)
    log.info("Loaded %d cells x %d genes", adata.n_obs, adata.n_vars)

    model = load_model(ckpt)

    if args.checkpoint_dir:
        ckpt_dir = Path(args.checkpoint_dir)
        emb = process_chunked(model, adata, args.chunk_size, args.batch_size, args.max_len, ckpt_dir)
        if emb.shape[0] != adata.n_obs:
            raise RuntimeError(f"Reassembled {emb.shape[0]} != n_obs {adata.n_obs}")
        adata.obsm["X_scprint"] = emb
        log.info("scPRINT embedding shape: %s", emb.shape)
    else:
        adata_out = run_embedder(model, adata, args.batch_size, args.max_len)
        adata = attach_embedding(adata, adata_out)

    log.info("Writing %s", out_path)
    adata.write_h5ad(out_path, compression="gzip")
    log.info("Done. obsm keys on output: %s", list(adata.obsm.keys()))
    return 0


if __name__ == "__main__":
    sys.exit(main())

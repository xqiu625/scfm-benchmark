import torch
from scprint2 import scPRINT2

CKPT = "/bigdata/godziklab/shared/Xinru/.huggingface/scprint2_ckpts/small-v2.ckpt"
print("Loading model...", flush=True)
model = scPRINT2.load_from_checkpoint(CKPT, precpt_gene_emb=None, gene_pos_file=None)
print("Loaded.", flush=True)

print("type(model.expr_encoder):", type(model.expr_encoder), flush=True)
print("normalization:", getattr(model, "normalization", None), flush=True)

genes = list(model.genes)
human = [g for g in genes if g.startswith("ENSG")]
print("n human (ENSG) genes:", len(human), flush=True)
print("human sample:", human[:10], flush=True)

# index of first human gene, to understand ordering / whether gene_pos is a
# global flat index or needs an organism-relative offset
first_human_idx = next(i for i, g in enumerate(genes) if g.startswith("ENSG"))
print("first human gene global index:", first_human_idx, flush=True)

print("num classes+1 (cell_embs_count-ish):", len(model.classes) + 1, flush=True)

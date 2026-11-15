import torch
from scprint2 import scPRINT2

CKPT = "/bigdata/godziklab/shared/Xinru/.huggingface/scprint2_ckpts/small-v2.ckpt"
print("Loading model...", flush=True)
model = scPRINT2.load_from_checkpoint(CKPT, precpt_gene_emb=None, gene_pos_file=None)
print("Loaded.", flush=True)

print("type(model.transformer):", type(model.transformer), flush=True)
print("expr_emb_style:", getattr(model, "expr_emb_style", None), flush=True)
print("n_input_bins:", getattr(model, "n_input_bins", None), flush=True)
print("organisms:", getattr(model, "organisms", None), flush=True)
print("classes:", getattr(model, "classes", None), flush=True)
print("use_metacell_token:", getattr(model, "use_metacell_token", None), flush=True)
print("cell_transformer:", getattr(model, "cell_transformer", None), flush=True)
print("attn_bias:", getattr(model, "attn_bias", None), flush=True)
print("d_model:", getattr(model, "d_model", None), flush=True)

genes = model.genes
print("n genes:", len(genes), flush=True)
print("gene sample:", list(genes)[:10], flush=True)
print("gene type:", type(genes), flush=True)

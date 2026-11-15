import re
import torch
from scprint2 import scPRINT2

CKPT = "/bigdata/godziklab/shared/Xinru/.huggingface/scprint2_ckpts/small-v2.ckpt"
print("Loading model...", flush=True)
model = scPRINT2.load_from_checkpoint(CKPT, precpt_gene_emb=None, gene_pos_file=None)
print("Loaded.", flush=True)

genes = list(model.genes)
# human Ensembl gene IDs: ENSG followed directly by digits (ENSG00000012345),
# NOT ENSGALG... (chicken), ENSGACG... etc -- a prior naive startswith("ENSG")
# check wrongly matched chicken genes too.
human_re = re.compile(r"^ENSG\d+$")
human = [g for g in genes if human_re.match(g)]
print("n TRUE human genes:", len(human), flush=True)
print("human sample:", human[:10], flush=True)
human_idx = [i for i, g in enumerate(genes) if human_re.match(g)]
print("first human global index:", human_idx[0] if human_idx else None, flush=True)
print("last human global index:", human_idx[-1] if human_idx else None, flush=True)
print("contiguous block?:", (human_idx[-1] - human_idx[0] + 1) == len(human_idx) if human_idx else None, flush=True)

print("pos_encoder:", model.pos_encoder, flush=True)
print("use_metacell_token:", model.use_metacell_token, flush=True)

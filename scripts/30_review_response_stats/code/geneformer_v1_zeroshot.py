#!/usr/bin/env python3
"""
Zero-shot Geneformer V1-10M control (roadmap J).

The manuscript's "LoRA adds only +0.006 AUC" null compares LoRA-fine-tuned
Geneformer-**V1-10M** against zero-shot Geneformer-**V2-316M** — different backbones,
so capacity is confounded with adaptation. This generates the missing WITHIN-backbone
control: zero-shot Geneformer-V1-10M on the SAME cohort / labels / patient-level CV the
LoRA run used, by reusing the LoRA script's data loader + tokenizer and swapping
fine-tuning for frozen mean-pooled embedding extraction.

Result isolates the true LoRA adaptation effect:
    LoRA-V1-10M (committed 0.705)  -  zero-shot-V1-10M (this)

Output: 05_RESULTS/geneformer_peft/geneformer_v1_zeroshot_mortality.csv
        05_RESULTS/geneformer_peft/geneformer_v1_zeroshot_emb.npy
"""
import sys
import os
import numpy as np
import pandas as pd
import torch
from torch.utils.data import DataLoader
from transformers import BertModel
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.metrics import roc_auc_score, balanced_accuracy_score

# reuse the exact LoRA data/tokenizer machinery so the control is apples-to-apples
sys.path.insert(0, '/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/finetuning')
from geneformer_peft_finetune import (  # noqa: E402
    Config, load_covid_data, load_geneformer_token_dict, GeneformerDataset,
)

OUT = '/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/geneformer_peft'


def main():
    cfg = Config()
    print(f'V1-10M base model: {cfg.PRETRAIN_CKPT}', flush=True)
    print(f'N_FOLDS={cfg.N_FOLDS}  SEED={cfg.SEED}  MAX_LEN={cfg.MAX_INPUT_LEN}', flush=True)

    adata, labels, groups = load_covid_data(cfg, task='mortality')
    token_dict = load_geneformer_token_dict(cfg)
    ds = GeneformerDataset(adata, labels, token_dict, max_len=cfg.MAX_INPUT_LEN)
    y = np.asarray(labels)
    print(f'cells={len(y)}  donors={len(np.unique(groups))}  '
          f'fatal={int(y.sum())}  survived={int((y == 0).sum())}', flush=True)

    dev = 'cuda' if torch.cuda.is_available() else 'cpu'
    # frozen base encoder (no LoRA, no classifier head) -> zero-shot representation
    model = BertModel.from_pretrained(str(cfg.PRETRAIN_CKPT)).to(dev).eval()
    print(f'loaded BertModel on {dev}; hidden={model.config.hidden_size}', flush=True)

    loader = DataLoader(ds, batch_size=cfg.BATCH_SIZE, shuffle=False, num_workers=4)
    embs = []
    with torch.no_grad():
        for bi, batch in enumerate(loader):
            ids = batch['input_ids'].to(dev)
            mask = batch['attention_mask'].to(dev)
            hs = model(input_ids=ids, attention_mask=mask).last_hidden_state
            m = mask.unsqueeze(-1).float()
            pooled = (hs * m).sum(1) / m.sum(1).clamp(min=1)   # masked mean pool
            embs.append(pooled.float().cpu().numpy())
            if bi % 200 == 0:
                print(f'  batch {bi}/{len(loader)}', flush=True)
    X = np.concatenate(embs).astype(np.float32)
    print(f'zero-shot embeddings: {X.shape}', flush=True)

    skf = StratifiedGroupKFold(n_splits=cfg.N_FOLDS, shuffle=True, random_state=cfg.SEED)
    aucs, baccs = [], []
    for tr, te in skf.split(X, y, groups):
        clf = LogisticRegression(C=1.0, solver='lbfgs', class_weight='balanced', max_iter=1000)
        clf.fit(X[tr], y[tr])
        p = clf.predict_proba(X[te])[:, 1]
        aucs.append(roc_auc_score(y[te], p))
        baccs.append(balanced_accuracy_score(y[te], clf.predict(X[te])))

    m_auc, s_auc = float(np.mean(aucs)), float(np.std(aucs))
    print(f'\n=== Zero-shot Geneformer-V1-10M mortality LR: AUC = {m_auc:.4f} +/- {s_auc:.4f} ===',
          flush=True)
    print('folds: ' + ';'.join(f'{a:.3f}' for a in aucs), flush=True)
    print(f'LoRA-V1-10M (committed) = 0.705 +/- 0.135', flush=True)
    print(f'=> within-backbone LoRA adaptation effect = {0.705 - m_auc:+.4f} '
          f'(positive = LoRA helped over zero-shot, same V1-10M backbone)', flush=True)

    os.makedirs(OUT, exist_ok=True)
    np.save(os.path.join(OUT, 'geneformer_v1_zeroshot_emb.npy'), X)
    pd.DataFrame([{
        'model': 'Geneformer-V1-10M-zeroshot', 'classifier': 'LR',
        'auc_mean': m_auc, 'auc_std': s_auc, 'bacc_mean': float(np.mean(baccs)),
        'fold_aucs': ';'.join(f'{a:.3f}' for a in aucs), 'n_cells': int(X.shape[0]),
        'n_donors': int(len(np.unique(groups))),
    }]).to_csv(os.path.join(OUT, 'geneformer_v1_zeroshot_mortality.csv'), index=False)
    print(f'SAVED -> {OUT}/geneformer_v1_zeroshot_mortality.csv', flush=True)


if __name__ == '__main__':
    main()

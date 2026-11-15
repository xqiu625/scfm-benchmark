# scDrugMap（Nat Commun 2025）阅读笔记与引用要点

文献：`41467_2025_Article_67481.pdf`（doi:10.1038/s41467-025-67481-2；
Wang Q, Pan Y, Zhou M 等，Song Q 通讯）。
归档位置：`revision/literature_audit/` 同级笔记（本文件）。

## 它是什么

最大规模的临床终点 scFM benchmark：8 个 scFM + 2 个 LLM（GPT4o-mini、
LLaMa3-8B），495,000 细胞 / 60 数据集，药物反应预测（敏感/耐药二分类），
三种组织/多癌种/多治疗方案。headline：scFoundation pooled 最强；UCE
fine-tune 后 cross-data 最强；scGPT zero-shot cross-data 最强。

## 协议弱点（= 我们的卖点，Discussion 用）

1. **pooled-data 主评估 = 细胞级随机 10-fold CV**（原文："randomly
   partitioning the data within each biological category"）——kapoor2023 定义的
   泄漏模式，F1~0.9+ 数字大概率被抬高。
2. **无 from-scratch 简单基线**：基线为 MLP / scDEAL / DREEP（本身即复杂模型）。
3. 无校准（Brier/ECE）、无 repeated-CV 不确定性（单次 10-fold + s.e.m.）。
4. LoRA fine-tune > layer-freezing 的结论同样建立在细胞级 CV 上。

## 可引用的事实

- 肝癌数据上 MLP/scDEAL/DREEP 与 scGPT/CellLM/LLaMa3/Geneformer/scBERT
  全部打平（F1=0.924，Supplementary Note 3）——**常规方法在其自身协议下
  仍不输 FM**，支持本文主线。
- cross-data（跨研究）设定存在且被称为 "zero-shot"——与我们 cross-disease
  的分布偏移**不同量级**（同癌种跨研究 vs 跨疾病×未见细胞类型），引用时
  必须带这个限定。

## 建议的 Discussion 句式（relation to prior evaluations 节）

> The largest clinical-endpoint scFM benchmark to date (scDrugMap, 495k cells,
> 60 datasets) still evaluates pooled data with cell-level random cross-validation
> — the leakage pattern documented by Kapoor & Narayanan — and omits
> from-scratch simple baselines; in its own liver-cancer arm, MLP, scDEAL and
> DREEP match every foundation model. Its cross-study "zero-shot" transfer is
> within-cancer-type, a far milder shift than our cross-disease, unseen-cell-type
> setting, so their scGPT-strong zero-shot result does not contradict ours.

—— 注意措辞：说"protocol 不同所以数字不可比"，不贬低对方工作。

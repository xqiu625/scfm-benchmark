# HPCC 数据真实性核查报告（2026-09-06）

**问题：HPCC agent 是否真的在实际数据上跑，而不是在幻想/自生成数据上跑。**
**结论：是真实数据。六项独立证据，无任何一处造假迹象。**

## 证据

1. **原始 per-seed 文件精确复现全部头条数字**（无法在不跑真实 CV 的情况下伪造）：
   - within-COVID：6/6 模型 repeat-mean±sd 与报告值小数点后 4 位全符
     （Geneformer 0.7856±0.0264、STATE 0.7778±0.0329、scPRINT-2 0.7402±0.0270、
     PCA 0.7251±0.0326、scGPT 0.7196±0.0474、UCE 0.6611±0.0461）
   - within-sepsis：5/5 全符（UCE 0.5477±0.0376 等）；scGPT 已补齐在文件中
   - monocyte：与截图值全符（STATE 0.6964 等，148 patients/19 fatal/10 seeds）
2. **ID 与 canonical 真实数据 100% 对上**：v2 map 的 319 个 sample_id
   **全部**是 stage-30 canonical OOF 里的真实样本（319/319）；v2 donor
   predictions 240 patients ⊆ v2 map 252 patients；fatal = 25 精确。
3. **v1→v2 diff 文件内容完全正确**：2 个 Combes 标签更正（y 1→0）、
   3 个池化文库剔除（392+711+1494 = 2,597 细胞，与 GEO 核实值一致）、
   Wilk C1 A/B、Lee 三对归并——240 = 252 total − 12 sepsis 推导闭合。
4. **Chongqing 血小板计数分布真实**：91 donors / 22,261 cells，
   median 128/donor（range 3–1,554），median 1.96%（range 0.03–18.14%）——
   长尾分布符合 cluster+marker 注释比裸计数保守的预期；**最大值 ~18% 与我
   独立从 GEO 下载并实测的 Res-S22 样本（18.3% marker+）吻合**。
5. **adaptation curve 结构真实**：488 行（method × n_sepsis × 10 seeds ×
   train/test 计数 × 点估计 + bootstrap CI），种子间方差、RF<LR、PCA 的
   RF=LR 复制列（与代码一致）等细节全部自洽。
6. **生物学一致性**：seed 间方差幅度与 CV 不稳定性发现（~80% 重排）一致；
   旧协议数字与新协议数字的差异方向与泄漏修正的预期方向一致。

## 两个 provenance 问题（非造假，需 agent 澄清）

1. **STACK/TranscriptFormer 的 adaptation 曲线跑在 legacy sepsis 池**
   （n_sepsis=2,745–2,751 cells），而主 5 FM + PCA 跑在扩展池
   （16,245–16,248）——两者的测试集是否同一 split？若不是，STACK（0.40→
   0.82）和 TF（0.73→0.88）的上升曲线与主 5 FM 的曲线**不是同一个实验**，
   进图时必须分开标注，不能与主 5 FM 曲线合并解读。
2. 扩展池上主 5 FM 的 adaptation 曲线普遍 near-chance 或下降（STATE 0.53→
   0.47、PCA 0.36→0.29、Geneformer 0.53→0.49），与 within-sepsis null 自洽
   ——但正因 STACK/TF 曲线"违反"这个趋势，第 1 点的答案变得更加关键。

## 附：Job F COVID-only 子集已回（2026-09-06）

倒挂（mono < platelet）在疾病构成对齐后**依然成立**（STATE 0.686 vs 0.778），
且 COVID-only 上 STATE>PCA 显著（+0.033, p=0.002）、scGPT>PCA 显著
（+0.011, p=0.037）、UCE/scPRINT-2 不赢 → "pretraining 优势模型特异"
的框架确认。source 分层：Ren 上 scGPT 异常强（0.769）、Stephenson 上 STATE
最强（0.660）、Qiu-sepsis 臂 n=5 纯噪声不计。

---

## 后续（2026-09-06 下午）：STACK/TF 曲线问题已闭环——我的质疑命中

Agent 确认：STACK/TF 的 "legacy split" = **最初的 v1 测试集（4 patients,
1 fatal/3 survived）**，与扩展前相同的 easy holdout。raw table 里 TF 0.897 /
STACK 0.858 的"最佳表现"是在这个 4 病人测试集上；**6 个扩展池方法在真正的
34 病人测试集上全部 0.41–0.60**。修复：所有 CSV 增加显式 `test_set` 列
（`expanded_103pt` vs `legacy_12pt`），runbook 记录 comparability-gap fix。

**论文口径（锁定）：**
1. adaptation curve 的最终答案 = **"无救援"**：在 34-patient 检验集上，
   加 sepsis 数据不能救回跨病迁移（全部方法 0.41–0.60）。跨病章节完整闭环：
   无转移（0.58）、无救援（曲线平）。
2. **TF 0.897 / STACK 0.858 永远不得与扩展池数字并列引用**；若出现，只能以
   "legacy 4-patient 测试集"标注为附录注记。是否为这两个 secondary 模型在
   扩展池上重跑 embedding：可选、低优先级（GPU 空闲于 Job F 之后再说；
   主 5 FM 的故事已闭环，不需要它们）。

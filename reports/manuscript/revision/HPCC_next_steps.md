# HPCC 接下来的顺序与修复（2026-09-04 定稿版）

> 依据：代码审查（`code_review_2026-09-04.md`）+ Job S 首批结果与 adaptation
> curve 决策。本文件是当前唯一行动清单；细节规格仍以
> `scripts/31_donor_level_revision/HPCC_RUNBOOK.md` 为准。

---

## 一、先修（🔴 3 项必须修）——不修则结果表不算数

### Fix-1：v2 映射文件补 Task-4 标签更正（10 分钟，最先做）

`code/build_patient_map_v2.py` 缺 `SEVERITY_FIX` → v2 map 里 GSM4995436/442
仍是 y=1，污染计数报告和 null model。

- 加同一个 dict：`GSM4995436→moderate、GSM4995442→severe`（y 1→0）；
  `GSM4995432→moderate`（y 不变）
- 重新生成 `donor_meta_with_patient_id_v2.csv`，diff 报告同步更新
- 预期：COVID 25 fatal patients / 27 fatal sid；全队列 fatal patients 27

### Fix-2：transfer 评分改 patient-level 聚合（小编辑）

`code/jobS_cross_disease.py` 和 `code/cross_disease_fatal_vs_survived_patient_v2.py`
的 AUC 点估计目前是 cell-weighted（bootstrap 按 patient 聚类但分数在细胞级）。
与 P0 estimand（一切评分 patient-level）不一致。

- 主口径：每 patient 的预测概率取均值 → patient-level AUC（bootstrap 同现
  200 次 patient 聚类重采样，在聚合分数上做）
- cell-weighted 保留为 secondary 并在输出列名里标清（如 `auc_cell_weighted`）

### Fix-3：transfer 主协议改 severe=survived（一个开关）

两个 transfer 脚本 COVID 侧当前是 severe-EXCLUDED。P0 已定 **severe=survived
为主协议**。

- severe-included（COVID survived = mild/moderate/recovered/severe）作为主口径
- severe-excluded 保留为敏感性，两版同报、列名标主次（`primary_` / `sens_severe_excluded_`）
- 不写任何解释性文字，只标版本

---

## 二、同批修（🟡 2 项应该修）

### Fix-4：Chongqing 聚合体比例报告

`code/integrate_chongqing_sepsis.py`：在 platelet 簇中统计单核 marker
（CD14/LST1/FCN1）共表达比例，写进 `chongqing_platelet_counts_per_donor.csv`
（设计文档的硬性要求）。

### Fix-5：数字同步（全文档）

- Chongqing = **91 donors / 17 fatal**（93 个 Adult Sepsis − 2 个 mort28 缺失）。
  改掉 runbook Job S、`sepsis_arm_expansion.md`、`P0_onepager.md` 里所有
  "108 patients / 20 fatal"
- COVID fatal = 25 patients（Task-4 更正后），全队列 fatal patients = 27

---

## 三、执行顺序（修复落位后）

| # | 任务 | 说明 | 依赖 |
|---|---|---|---|
| 1 | **Fix-1**（v2 map 重建） | 10 分钟 | 无 |
| 2 | **Fix-2/3**（transfer 脚本改完） | 小编辑 | 无，可与 1 并行 |
| 3 | **Job H2 收尾** | repeated-CV（10×5）within-COVID：5 FM + PCA-200 + null，主表 = repeat-mean±sd + 同 seed 配对检验 | Fix-1 |
| 4 | **scGPT embedding 落地后**（Chongqing，已在队列） | within-sepsis 表补 scGPT 行；transfer 三分析同补 | GPU 队列 |
| 5 | **transfer 双向（Chongqing）按 Fix-2/3 重出** | patient-level 主口径 + severe-included 主协议 | Fix-2/3 + scGPT（可选并入） |
| 6 | **adaptation curve（已批准，option 1）** | sepsis 池 103（12+91）patients，10 seeds 训练池子采样设计不变；与旧曲线同协议可直接比较 | Job S Step 2–4 |
| 7 | **Job F monocyte** | Ren + Stephenson 为主（Qiu 可选）；5 FM embedding（GPU）；评估从第一天就 patient-level + repeated-CV | GPU 排在 H2 后；标签继承修正后的 v2 map |
| 8 | **Fix-4/5** | 随上述同批完成 | — |

**GPU 调度**：H2（CPU）期间 S 的下载/scGPT embedding 并行；F 的 embedding 排在
S 之后。

---

## 四、已定协议与决策（勿再变更）

- **Repeated-CV（10 seeds × 5-fold）只适用于 fold 型步骤**（within-COVID /
  within-sepsis / PCA / null）。cross-disease = 固定 split + 200 次 donor
  bootstrap（LR 确定性，seeds 无效；RF 可加 10 个 random_state，可选）。
  adaptation curve 自带 10-seed 训练池子采样设计，不改不重跑。
- **adaptation curve 按原计划跑（option 1）**：它回答"加 sepsis 数据能否救回
  跨病迁移"，within-sepsis null 不预设其答案。
- **within-sepsis 结果按真实 null 报告**（FM 0.51–0.55、PCA 0.43），时点混杂
  （fatal 采血更早，median 10 vs 18.5 天）写进 Discussion；**FM > PCA 的反转
  （4/5 FM，q=0.0037；scGPT n.s. q=0.40；UCE 显著领先 STATE q=0.010，
  均近随机水平）如实报，不解释**。
- 旧的 12-patient transfer 数字（STATE 0.80+）vs 新 0.58：两个候选解释
  （small-n artifact / Qiu-T0 与 Chongqing-d16 时点错配）并录，不裁定。
- 写作规则：数字如实摆，不写更正过程，不给自洽解释。

## 五、回传清单

1. 修正后的 `donor_meta_with_patient_id_v2.csv` + diff
2. within-COVID 主表（5 FM + PCA + null，repeat-mean±sd + 配对检验）
3. within-sepsis 主表（含 scGPT）
4. transfer 双向表（patient-level 主口径 + severe 两协议，含 scGPT）
5. adaptation curve（103-patient 池，per-seed 行）
6. monocyte（Job F）：各 FM + PCA 的 in-distribution vs platelet 对比表

---

## 六、Job F 结果与跟进（2026-09-06）

**monocyte（in-dist, 148 patients/19 fatal）vs platelet（OOD, within-COVID 240/25），repeated-CV：**

| FM | Monocyte | Platelet | Δ |
|---|---|---|---|
| STATE | 0.696±0.037 | 0.778±0.033 | −0.081 |
| scPRINT-2 | 0.672±0.058 | 0.740±0.027 | −0.068 |
| scGPT | 0.677±0.037 | 0.720±0.047 | −0.043 |
| UCE | 0.648±0.056 | 0.661±0.046 | −0.014 |
| PCA-200 | 0.654±0.043 | 0.725±0.033 | −0.072 |
| Geneformer | N/A（container issue） | 0.786±0.026 | — |

**判读（三层，分开写）：**
1. **设计的问题**：FM 在 in-distribution 是否赢 PCA——**混合答案**（STATE 和 scGPT 显著赢，
   scPRINT-2 不赢，UCE 打平或更低）→ "pretraining 优势"即使在 in-distribution 也是
   模型特异的。与我们"conditional, model-specific"的总框架一致。
2. **倒挂（全体 mono < platelet）**：**PCA-200 同样掉 0.072**——这是关键内控：
  掉幅与 pretraining 无关，是任务/队列属性（monocyte mortality 在这个队列就是更难）。
   不是 FM 的失败，也不能当 FM 的功劳。biological note（"consistent with"句式）：
   血小板在 COVID 死亡机制上比单核细胞更下游更近（immunothrombosis 文献，intro 已引）。
3. **遗留混杂**：monocyte 队列 = 143 COVID（Ren+Stephenson curated）+ 5 sepsis
   （Qiu marker 阈值法，最低 provenance），platelet 对比是 COVID-only →
   **已批准跑 COVID-only 子集**（Ren+Stephenson 143 patients，免新 CV 仅重聚合）。

**跟进请求（已发给 agent）：**
1. COVID-only monocyte 子集结果（同 repeated-CV 报告口径）
2. Geneformer monocyte embedding 若 container issue 修复便宜则补跑，否则表注 N/A
3. 按 source（Ren / Stephenson / Qiu-marker）分层报 monocyte AUC——检验
   provenance/质量担忧（Qiu marker 阈值法那 5 个病人）
4. Geneformer 之外的 FM 数字定稿后回传，进 within-COVID vs monocyte 对比终表

---

## 七、Job J — GSE151263 (Jiang) 剔除重跑（2026-09-07 用户决定 B）

**原因**：GSE151263 无 per-patient 结局记录（原文仅报文献死亡率），其 7 个
"survivor" 标签是研究设计推断，不可验证。决策：**GSE151263 保留在 platelet 数据集与队列描述中（它确实是 sepsis 数据），
仅从一切带结局标签的 mortality 分析中剔除**（Data_NO `302004data18`，
7 donors / 476 cells 不进任何 labeled 训练/测试池），使所有进入分析的标签均可溯源。

**范围（CPU，<1 h）**：
1. 仅在 mortality 分析管线中过滤 `Data_NO == '302004data18'`（476 cells）
   （数据集本身保留）。mortality 分析基数：sepsis 96 patients
   （5 Qiu GSE167363 + 91 GSE279452）/ 19 fatal；GSE167363 测试集
   5 patients / 2 fatal。
2. **重跑 reverse transfer**（sepsis→COVID）：训练池 103→96 patients。
3. **重跑 adaptation curve**：sepsis 训练池同步剔除 Jiang 细胞。
4. within-COVID（R1）与 within-sepsis（R2，GSE279452-only）**不受影响，不重跑**；
   forward transfer（测试集=91 GSE279452）不受影响。
5. 报告：reverse transfer 新旧对照（96 vs 103 训练池）、adaptation 曲线新终点；
   预期结论方向不变，仅汇报实际值。

---

## 八、外部评审后续任务（Jobs K–P，2026-09-07；依据 platelet_fm_benchmark_review.md）

### Job K — QC-metrics-only 基线（评审 A1，最高优先，可推翻主结果）
1. 从 canonical h5ad 计算/读取 per-cell QC：n_counts（nUMI）、n_genes、
   pct_counts_mt（X 是 raw counts，自行计算即可）。
2. 聚合到 patient：mean/median nUMI、nGene、%mito、n_cells。
3. 用这些 QC 特征做 LR，跑与主分析**完全相同**的 repeated-CV
   （10 seeds × 5-fold patient 分组），报 patient-level AUC repeat-mean±sd。
4. 分层报告：per-patient QC（nUMI/nGene/%mito/platelet cells）按 cohort
   （Data_NO）和按结局（y）的中位数分布表。
5. **决策规则（同 null model）：QC-only AUC ≥ 0.7 within-COVID → 立即停手报告**
   （主结果是技术假象，框架重构）；<0.7 → 作为 supplementary null 写入。
6. 文档化 platelet calling 规则：legacy = 2024 IJMS 整合（cluster+marker）；
   GSE279452 = cluster+marker 12 基因 panel（`integrate_chongqing_sepsis.py`），
   并补 Fix-4 的聚合体比例报告（platelet 簇中 CD14/LST1/FCN1 共表达）。

### Job L — 患者级配对推断正式版（评审 A3，定头条措辞）
1. 对 10 个 seed 各自的 `oof_patient_<FM>_seed<NN>.csv`（05_RESULTS/
   revision_donor_level_v2/）：per-patient 聚合（mean）→ 患者级配对
   bootstrap（B=10,000，重采样 patient）计算头条对比 ΔAUC：
   每个 FM − PCA-200、STATE − Geneformer。
2. 同时报：(a) 每 seed 的 Δ 与 p；(b) seed 集成（每 patient 10 个 OOF 分数取
   均值）后的配对 bootstrap；(c) 各 AUC 的 Hanley–McNeil SE（25 fatal/215
   survivors）；(d) 模型间分数相关系数。
3. 预期管理：seed-42 预演结果 STATE−PCA p=0.059、Geneformer−PCA p=0.002——
   正式版以此为准写头条（可能"Geneformer 显著、STATE 边缘"）。

### Job M — GSE279452 年龄/位点分层（评审 A2）
1. 从 GEO characteristics 取 per-sample age：91 个分析病人的年龄分布
   （直方图/中位/范围）；确认全部成人（group=='Adult Sepsis'，≥18 岁）；
   若有 pediatric 混入，报告并剔除。
2. 按年龄层（如 <40 / 40–60 / >60）报 fatal 计数。
3. 把 site（单中心）与 chemistry（10x）作为队列级常量写入 cohort 表，
   供 transfer 混杂声明引用。

### Job N — LOBO on v2（评审 C12，降为 secondary）
在 v2 数据上做 leave-one-cohort-out（within-COVID；仅 cohorts 同时含
≥1 fatal 与 ≥1 survivor 的可计算 AUC；其余队列留训练池）。
报 per-cohort AUC + 每 FM 的 mean/σ；与 pooled 值并列。

### Job O — monocyte 病人匹配对照（评审 C13）
取同时有 platelet 和 monocyte 的病人交集，按每病人细胞数匹配（下采样），
在**同一病人集**上重算 monocyte vs platelet 的 AUC 对比（含 PCA-200）。
交集过小则如实报告并声明非匹配设计的局限。

### Job P — 清理 legacy 4-patient 数字（评审 C11）
从 `results_patient_v2/adaptation_curve/` 的全部 CSV 中**删除**
`test_set == legacy_12pt` 的行（先备份一份到 `05_RESULTS/attic/`）。
TF 0.897 / STACK 0.858 不再存在于任何下游可引用的记录中。

**回传**：Job K 的 QC 基线 AUC（含决策规则结果）、分层 QC 表；Job L 全套
配对推断表；Job M 年龄分布与分层 fatal 表；Job N per-cohort AUC 表；
Job O 匹配对照结果（或交集过小的声明）；Job P 清理确认。

---

## 九、Job Q — 配对 bootstrap 补全（供图标注用，CPU <30 min）

Job L 只跑了 6 对（5 FM−PCA + STATE−Geneformer，within-COVID）。为让 F3/S2/F4
的显著性标注全部使用诚实推断（seed-ensemble 患者级配对 bootstrap，
B=10,000，与 Job L 同法），补跑：
1. **within-COVID 全部 15 对**（6 模型两两，含 PCA-200；null 不进）。
2. **within-sepsis 全部 15 对**（当前 S2/F3 用的是 per-seed Wilcoxon，
   属反保守口径，需替换）。
3. **monocyte FM−PCA 全部 5 对**（当前 F4 森林图的 p 同为 per-seed
   Wilcoxon：STATE 0.002 / scGPT 0.037，需替换）。
4. **QC 基线 within-COVID 变体**（评审 A1 口径核对）：当前
   `qc_baseline_auc.csv` 为 pooled 252 patients/27 fatal；补跑 COVID-only
   （240 patients/25 fatal）的 QC-metrics-only AUC，输出
   `qc_baseline/qc_baseline_auc_covid_only.csv`。
输出格式与 `paired_inference/seed_ensemble_paired_bootstrap.csv` 相同
（pair, n_patients, delta_auc, ci_lo, ci_hi, p_boot），分任务各一份
（`paired_inference_full/within_covid_15pairs.csv`、`within_sepsis_15pairs.csv`、
`monocyte_vs_pca_5pairs.csv`）。

---

## 十、Job R — seed-ensemble 版的 winner 校正与脆弱性复核（CPU <30 min）

评审 round-3 的决定项：主对比必须用 winner-corrected max-statistic 检验。
本地 seed-42 预演（donor_predictions，STATE 最弱 seed）结果：
familywise P=0.0003（偏乐观）、Geneformer−PCA fragility 0/25、
STATE−PCA fragility 17/25。请在 **10-seed 集成 OOF 分数**上正式复算：
1. **max-statistic familywise test**：B=10,000 患者重采样，P = 无一 FM
   AUC > PCA 的重采样比例（"does any FM beat PCA"）。
2. **leave-one-fatal-out fragility**（25 refits）：Geneformer−PCA 与
   STATE−PCA 各自 p 范围与 p>0.05 的 refit 数。
3. **分数相关矩阵（Pearson）+ 共识预测器**（5 FM 均值）的 AUC
   ——本地 seed-42 值：相关 0.61–0.88、consensus AUC 0.6995
   （低于最佳单模型与 PCA，ensembling 稀释）。
输出 `paired_inference_full/winner_test_ensemble.csv`、
`fragility_ensemble.csv`、`score_correlation_ensemble.csv`、
`consensus_ensemble.csv`。

**STATUS（2026-09-08 13:48，sbatch 28203972，16m54s，全部通过 sanity gate
n=240/fatal=25）：**

- **familywise P = 0.0045**（预演 0.0003，同量级 <0.01 → 按规则**不改正文**，
  只填表）。obs max Δ = 0.0644 [0.0145, 0.1283]（Geneformer 单模型即为 max）。
- 逐模型 vs PCA（未做多重校正的点估计，供参考）：STATE +0.0523(p=0.0328)、
  Geneformer +0.0644(p=0.0264)、scPRINT-2 +0.0138(n.s.)、scGPT −0.0032(n.s.)、
  UCE −0.0659(n.s.，方向为负)。
- **fragility**：Geneformer−PCA 6/25 refits 失显著（预演 0/25，比预演脆弱
  但仍 19/25 稳健）；STATE−PCA 7/25 失显著（预演 17/25，**比预演更稳健**）。
  两条口径方向一致：ensemble 版下 STATE 的边际显著性其实比单 seed 预演更可靠。
- **consensus**：5-FM 均值 AUC=0.7678，**反超 PCA-200 (+0.0255)**（预演
  consensus 0.6995 曾低于 PCA，ensemble 版反转——预演的"consensus 被稀释"
  结论不成立于正式版，如实记录不解释）。仍低于最佳单模型 Geneformer
  (0.8067, Δ=−0.0389)。
- 相关矩阵：0.707–0.904（预演 0.61–0.88，范围上移但量级一致）。
- 本地已同步 4 个 CSV 至 `results_patient_v2/paired_inference_full/`。


---

## 十一、Job J-补 — reverse transfer 96 池补齐 + legacy 参考值确认（2026-09-08）

背景：§七 item 2 指定 reverse transfer 训练池 103→96（剔 Jiang），但回传的
`cross_disease_patient_jobJ_20260907_154957.csv` 其 reverse 行为
n_train=2946 cells = 5-patient GSE167363-only 池，未实现 96 池
（91 GSE279452 + 5 Qiu = 25,207 cells）。当前 P0/R3 的逆向数字以 91 池
（chongqing 2026-09-04，STATE-LR 0.689 [0.583, 0.785]）为准；jobJ 文件的
5-patient 池两个方向已转 SI 小样本参考注（sn:gse167363_testset）。

1. **先查** `05_RESULTS/revision_donor_level_v2/` 是否已有 96 池 reverse 输出
   （跑过未回传）；有则直接回传，不必重跑。
2. **若无则补跑** reverse transfer（sepsis→COVID）：训练池 = 96 patients
   （GSE279452 91 + GSE167363 5，剔除 GSE151263），协议与 chongqing 完全
   一致（severe=survived 主协议 + severe-excluded 敏感性，LR+RF，200 次
   patient bootstrap），CPU <1h，输出格式同 `cross_disease_chongqing/`。
3. **顺带确认**：剔除 Jiang 后 legacy 池（5 patients）adaptation curve 各模型
   值（供 SI 小样本参考注更新；SI 注的 transfer 值已按 jobJ 文件更新为
   1.00/0.833/0.667）。

**回传**：96 池 reverse 表（或"已存在"确认 + 路径）、legacy adaptation 值；
runbook §十一 追加 STATUS 块。
**决策规则**：96 池回传后即为逆向 canonical 值，R3/P0/Abstract 同步更新；
若与 91 池值（0.689）差 <0.02 则只更新表格不动文字。

**STATUS（2026-09-08 14:44，sbatch 28205555，15m26s，exit 0）：**

- **步骤1（先查）**：`05_RESULTS/revision_donor_level_v2/` 无既有 96 池
  reverse 输出，确认需补跑。
- **步骤2（补跑）**：新脚本 `code/cross_disease_96pool.py`（合并
  `jobS_cross_disease.py` 的 Chongqing 91-patient 加载逻辑与
  `cross_disease_fatal_vs_survived_patient_v2.py` 的 Jiang 剔除 legacy
  Qiu 5-patient 加载逻辑，输出 schema 与 `cross_disease_chongqing/` 一致）。
  验证 pool 组成：Qiu 2,946 cells/5 patients + Chongqing 22,261 cells/
  91 patients = 25,207 cells/96 patients（与预期完全一致）。输出：
  `05_RESULTS/revision_donor_level_v2/cross_disease_patient_96pool_20260908_142954.{csv,json}`。
  **Sepsis→COVID（96 池，240 patients/25 fatal 测试集）主协议 LR/RF：**

  | FM | LR | LR 95% CI | RF |
  |---|---|---|---|
  | STATE | 0.7209 | [0.621, 0.809] | 0.7113 |
  | Geneformer | 0.7118 | [0.612, 0.814] | 0.6358 |
  | UCE | 0.6928 | [0.591, 0.788] | 0.6074 |
  | scGPT | 0.6417 | [0.525, 0.765] | 0.6221 |
  | scPRINT-2 | 0.6186 | [0.493, 0.743] | 0.5167 |

  vs. 91 池旧值（STATE-LR 0.689）：新值 0.7209，差 **0.032 ≥ 0.02 阈值**
  → **触发文字更新**（非仅换表）。STATE 仍最强，但 Geneformer LR 点估计
  （0.7118）与 STATE 相差仅 0.009，不再支持"唯 STATE"表述——P0 §3d 已
  改写为 STATE/Geneformer 并列陈述，CI 全重叠如实注明。
- **步骤3（legacy 参考值确认）**：`05_RESULTS/attic/adaptation_curve_bootstrap_chongqing_pre_jobP_backup/summary.csv`
  证实只有 STACK/TranscriptFormer 有 `legacy_5pt` 行（5 个主 FM+PCA-200
  均在 `expanded_96pt` 上跑，无独立 legacy-only 版本）。终点（n=all）值：
  STACK 0.7289、TranscriptFormer 0.7646——与此前 Job P 报告的数值一致，
  无需变更。
- **已更新**：`P0_onepager.md` §3d（96 池 reverse 表 + STATE/Geneformer
  并列措辞）。`main.tex` 未含此数字（尚未进入 revision 数字），无需改。
- 4 个输出文件（csv/json）已同步至本地
  `scripts/31_donor_level_revision/results_patient_v2/cross_disease/`。


---

## 十二、Job S — 主对比重建为 permutation max-statistic null（round-4 决定项，CPU <1h）

round-4 评审 item 1（决定性）：Job R 的 familywise P=0.0045 是**未校正的
max-tail**——重采样患者再取 max 不构造零分布（winner's curse 被重估而非
移除）；评审模拟显示该程序在全局零下 type I ≈ 0.36，且会产生比未校正
pairwise p 更小的"p 值"（正是我们 0.0045 vs 0.029 的模式）。必须重建：

1. **Permutation null（主）**：用 Job L/R 同款 seed-ensemble 患者分数
   （240 patients/25 fatal）。患者级打乱 mortality 标签 nperm=10,000 次；
   每次计算 5 FM + PCA-200 的 AUC，记 max_j(AUC_j − AUC_PCA)；
   P = (1 + #{max_perm ≥ max_obs})/(1+nperm)。只置换标签、不重训
   （scoring-step null，正文注明较弱 null）。输出 obs 每模型 delta 与 obs max。
2. **Recentered bootstrap（交叉验证）**：B=10,000 患者重采样；每次取
   max_j(δ_j^b − δ_j^obs)；P = fraction(recentered max ≥ max_j δ_j^obs)。
3. **LOFO on primary**：25 个 fatal 逐一剔除，每次在 239 人上重算 #1 的
   permutation P（nperm=10,000）；输出每 refit 的 P、范围、P>0.05 个数。
4. **SE 一致性检查（item 6）**：15 对 within-COVID 对比，输出 bootstrap SE
   （B=10,000 分布的 std）与 DeLong 配对 SE（用 ensemble 分数的实际相关），
   标记 bootstrap/DeLong < 0.8 的对比。
5. **Concatenated-pooled LOBO（item 4）**：5 个 LOBO fold 的 held-out 患者
   预测拼接，每模型算一个 pooled AUC；输出 per-cohort n/fatal 表。

sanity gate：n=240/fatal=25。输出至 `paired_inference_full/`：
`primary_permutation_null.csv`、`primary_recentered_bootstrap.csv`、
`fragility_primary_lofo.csv`、`se_check_delong.csv`、`lobo_pooled.csv`。
回传 + runbook §十二 追加 STATUS 块。

**决策规则**：permutation P < 0.05 → 头条成立（正当防守），R1/SI/P0 替换
P 值与 fragility；P ≥ 0.05 → 转全-null 上界框架（标题/摘要/R1 重写，
round-3 bottom line 适用），**立即上报，不要自行润色数字**。

**STATUS（2026-09-08 20:48-20:53，sbatch 28210993 主分析 33m50s exit 0 +
28205525 pooled LOBO 13m03s exit 0，两次因路径 bug 失败重跑后成功，sanity
gate n=240/fatal=25 通过）：**

- **⚠ 结论：permutation P = 0.18098 ≥ 0.05 — 触发升级，不自行润色，如实上报。**
  Job R 的 familywise P=0.0045（缺陷方法：resample-再取 max 不构造零分布）
  与本次的诚实 permutation null 之间差 40 倍量级，正是评审预警的模式。
- **(1) permutation null（主）**：obs max Δ = 0.0644（Geneformer − PCA-200）；
  **P = 0.18098**（nperm=10,000）。逐模型 obs delta：STATE +0.0523、
  Geneformer +0.0644、scPRINT-2 +0.0138、scGPT −0.0032、UCE −0.0659。
- **(2) recentered bootstrap（交叉验证）**：**P = 0.07280**——方向与 (1) 一致，
  同样不显著，非孤立结果。recentered-max 95% range [−0.0316, 0.0828]（跨零）。
- **(3) LOFO on primary（25 refits）**：**25/25 refits 全部 P>0.05**
  （范围 [0.14109, 0.28777]）——无一次剔除单个 fatal 病人后转为显著，
  结论不依赖任何单一病人，是稳定的不显著。
- **(4) SE 一致性检查**：15 对中 **8 对** bootstrap/DeLong ratio < 0.8
  （STATE-UCE 0.736、STATE-scPRINT-2 0.754、STATE-scGPT 0.657、
  STATE-PCA 0.874 达标未标、UCE-Geneformer 0.623、UCE-scPRINT-2 0.703、
  UCE-scGPT 0.751、Geneformer-scGPT 0.694、scPRINT-2-scGPT 0.577 等）——
  多数配对的 bootstrap SE 系统性小于 DeLong 近似 SE，即本轮修订全程使用的
  patient-level bootstrap 可能低估了方差，不只是本项分析的问题，值得追溯
  之前所有 bootstrap 显著性结论。
- **(5) Concatenated-pooled LOBO**：STATE 0.7525、Geneformer 0.7481、
  scPRINT-2 0.7514、UCE 0.6813、scGPT 0.6376（n=206 pooled patients/25
  fatal，5 个可计算 cohort，0 例重复病人）。与 within-COVID pooled 值方向
  一致但普遍略低（合理，LOBO 是更严格的跨队列泛化检验）。
- 4 个 CSV 已同步至本地 `results_patient_v2/paired_inference_full/`。
- **本任务未修改任何正文/SI/P0 文字**（按规则，升级情形下不自行润色）——
  等待用户对"头条论断是否成立"的决定。


---

## 十三、Job T — QC 增量价值检验 + range-statistic + 方法 vs null + DeLong 全面替换（round-5，CPU <2h）

依据 round-5 评审 order of work 1–4。全部用与 Job L/R/S 相同的 seed-ensemble
患者分数、v2 映射、seeds 42–51。

1. **QC 增量价值检验（决定性；同时回应 Adam 原始质疑与 round-5 item 1）**：
   基线协变量 = 患者级 mean nUMI + nGene + %mito + cell count + cohort one-hot。
   repeated-CV（10×5，patient 分组，预处理只在训练折内 fit）下跑 7 组 LR：
   (a) 仅基线协变量；(b) 基线 + 表示（5 个 FM embedding + PCA-200，共 6 组）。
   ensemble OOF 后，DeLong 配对检验 base vs base+rep（每组一次）。
   输出：每组 AUC、Δ、DeLong p。
2. **QC 诊断三件套**：(a) QC 基线 sign-oriented AUC（取反）+ 患者级 bootstrap CI；
   (b) 每 cohort 内 QC 特征与 mortality 的点二列相关表；(c) repeated-CV 每折的
   cohort 组成表（n patients / n fatal per cohort per fold）。
3. **Range-statistic permutation（item 3，leaders−UCE 升级为 confirmatory）**：
   与 Job S 相同的 10,000 次标签置换（rng seed 20260908，未保存则同种子重跑）：
   每次记录 5 个 FM 的 max − min AUC；观测值 0.1302（Geneformer−UCE）；
   P = (1 + #{range_perm ≥ 0.1302})/(1+10000)。输出 null 分布 50/90/95/99 百分位。
4. **方法 vs count+disease null（item 4）**：6 个方法（5 FM + PCA-200）各自与
   null 模型的 ensemble OOF 分数做 DeLong 配对检验（null 的 per-seed 患者 OOF
   若无则按同协议重生成）。输出 6 行 Δ/p。
5. **DeLong 全面对比替换（item 6）**：within-COVID 15 对、within-sepsis 15 对、
   monocyte 5 对（FM−PCA），全部报 DeLong SE 与双尾 p（ensemble 分数），
   作为正文报告的 SE 来源（bootstrap 仅保留给 max/range 统计量）。

sanity gate：n=240/25（within-COVID）、91/17（sepsis）、148/19（monocyte）。
输出至 `paired_inference_full/` 与 `qc_baseline/`，回传 + runbook §十三 STATUS。

**决策规则**：#1 若 6 组增量全不显著 → "信号可被 QC+计数+队列解释"，框架再降
一级，**立即上报**；若增量显著（尤其 FM 组）→ QC 混淆正式关闭，写入正文。
#3 P<0.05 → leaders−UCE 升为 confirmatory secondary result 写入 R1/摘要。

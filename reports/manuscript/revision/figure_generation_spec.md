# Figure Generation Spec v2 — 修订版全部数据图（Python/matplotlib/seaborn）

> 交接给绘图 agent。视觉语法参照 FigureYa 模板库
> （`/Users/xinruqiu/Documents/Skill-Stack/小丫画图/INDEX.html`，R/ggplot2 模板）——
> **只借格式，用 Python 复刻，不装 R**。每图标注对应的 FigureYa 模板名。
> **铁律**：数字必须与 `P0_onepager.md` §3 最终值一致——只许从指定 CSV 重算，
> 不许手打、不许引用旧版图。输出：PDF（矢量）+ PNG（≥300 dpi）→
> `reports/manuscript/figures/v2/`。英文图面。

## 全局风格

- 字体 Arial/Helvetica；浅灰虚线网格或无网格；direct label 优先；面板字母 A/B/C/D。
- 疾病色：COVID-19 `#C0392B`、Sepsis `#E67E22`、Healthy `#7F8C8D`。
- 模型固定色：STATE `#2C3E50`、UCE `#8E44AD`、Geneformer `#2980B9`、
  scPRINT-2 `#27AE60`、scGPT `#D35400`、PCA-200 `#7F8C8D`、null `#BDC3C7`。
- 显著性只标 q<0.05（值从 pairwise CSV 取，不自算 p）。
- **不允许整篇都是 bar chart**——按下列指定格式出图。

## 数据清单（repo 内路径，已验证，只读）

- `scripts/31_donor_level_revision/results_patient_v2/within_covid/`：
  `per_seed_auc_repeated.csv`、`donor_level_auc_repeated.csv`、
  `donor_level_pairwise_repeated.csv`、`null_model_auc_repeated.csv`、
  `donor_predictions_<FM>.csv`（y、n_cells、mean/median/trim、mean_rf）
- `.../within_sepsis/`：同名三件套（Chongqing）
- `.../monocyte/`：`donor_level_auc_repeated.csv`、`indistribution_vs_platelet.csv`、
  `monocyte_source_stratified_auc.csv`、`pairwise_vs_pca200_covid_only.csv`、
  `per_seed_auc_repeated.csv`
- `.../cross_disease_chongqing/cross_disease_chongqing_20260904_233237.csv`
- `.../adaptation_curve/all_rows.csv`、`summary.csv`、`split_documentation.txt`

---

## Figure 2 — within-COVID（240 patients / 25 fatal）

- **A（主）：Lollipop 图**（仿 FigureYa19Lollipop 的"点+茎"语法）——
  x = 模型（含 PCA、null），y = patient-level repeat-mean AUC（LR），
  茎端圆点 = mean，每模型后方叠加 10 个 per-seed 小 jitter 点（alpha≈0.5），
  茎长或误差帽 = repeat-sd。锚点：Geneformer 0.786、STATE 0.778、
  scPRINT-2 0.740、PCA 0.725、scGPT 0.720、UCE 0.661、null 0.554。
  显著性括注：STATE vs PCA（q=0.003）；Geneformer vs STATE（n.s.）。
- **B：多面板 patient-level ROC 曲线**（仿 FigureYa102multipanelROC）——
  从 `donor_predictions_<FM>.csv` 的 y 与 mean（LR 分）绘制各模型 ROC
  （patient 级），单面板多曲线 + 各 AUC 图例（值 = repeat-mean）；
  附一条 PCA 与 null 的对比虚线。对角线灰色。
- **C：box + density 组合图**（仿 FigureYa227boxdensity）——x = 模型，
  y = 10 个 per-seed AUC：上层 density 曲线、下层 boxplot+jitter，
  直观呈现 fold 不稳定性（repeated-CV 动机）。
- 输出 `Figure2_within_covid.pdf/png`。

## Figure 3 — within-sepsis null（91 patients / 17 fatal）

- 单面板：**violin + 内嵌 jitter + mean 点**（不用 bar）——x = 模型
  （5 FM + PCA），y = per-seed AUC（LR）；y 轴 [0.25, 0.75]，
  0.5 随机参考线（灰虚线）。PCA 显著低于 4/5 FM（q=0.0037）用
  星号注在 PCA 组上方，并在 caption 注明 scGPT vs PCA 不显著（q=0.40）。图内小字："28-day mortality; draws median 16 d
  post admission"。锚点：UCE 0.548、scPRINT-2 0.523、Geneformer 0.515、
  STATE 0.511、PCA 0.428。
- 输出 `Figure3_within_sepsis.pdf/png`。

## Figure 4 — in-distribution vs OOD（monocyte vs platelet）

- **A（主）：dumbbell/slope 图**——每模型一行，左点 monocyte
  （143 patients/17 fatal，COVID-only cut）、右点 platelet（240/25），
  水平连线并按 platelet 值降序排。数据 `indistribution_vs_platelet.csv`
  的 `COVID_only_Ren_Stephenson` 行（含 Geneformer）。
- **B：森林图**（仿 FigureYa90森林图）——monocyte 上 FM vs PCA-200 的
  ΔAUC + 95% CI（数据 `pairwise_vs_pca200_covid_only.csv`）；
  参照线 Δ=0；显著者（STATE p=0.002、scGPT p=0.037）加粗/标色，其余灰。
- **C（可选小图）**：source 分层点图（`monocyte_source_stratified_auc.csv`，
  Ren / Stephenson / Qiu-sepsis 三组；Qiu n=5 标注 "n=5, noise"）。
- 输出 `Figure4_indist_vs_ood.pdf/png`。

## Figure 5 — cross-disease transfer（Chongqing 扩展池）

- **森林图**（仿 FigureYa90森林图，双面板）：面板 A = COVID→Sepsis
  （91 patients/17 fatal），面板 B = Sepsis→COVID（240/25）；
  每模型一行，点 = patient-level AUC，横线 = bootstrap 95% CI；
  LR 圆点、RF 三角（同面板两种 marker）。数据
  `cross_disease_chongqing_20260904_233237.csv` 的 primary_* 列。
  0.5 参照线（灰虚线）；0.65 参考线（点线，旧稿阈值，便于对照）。
  锚点：正向 0.644–0.684 全重叠；逆向 STATE-RF 0.745 [0.655, 0.834] 最高。
- **禁止出现** legacy 4-patient 测试集数字（TF 0.897 / STACK 0.858）。
- 输出 `Figure5_cross_disease.pdf/png`。

## Figure 6 — adaptation curve（32-patient 检验集）

- **线 + ribbon**（仿 FigureYa37correlationV2）：x = sepsis 训练细胞数
  （log 标度，0→all），y = AUC，每方法一条线 + 10-seed mean±sd 半透明带；
  方法 = 5 FM + PCA-200。**只取 `all_rows.csv` 中
  `test_set == expanded_103pt` 的行**；`summary.csv` 校验聚合。
  0.5 参照线。终点（n=all）全体落 0.41–0.60。
- **不画 STACK/TranscriptFormer**（legacy 测试集，不可比）。
- 输出 `Figure6_adaptation_curve.pdf/png`。

## Figure 7 — calibration（within-COVID，patient-level）

- 双面板 reliability diagram：A=LR（`donor_predictions_<FM>.csv` 的 mean）、
  B=RF（mean_rf）；每模型一条曲线（10 等宽分箱：箱均预测概率 vs 实际死亡率，
  点面积∝箱内人数），对角虚线；图例标各模型 ECE（10 箱人口加权 |p−y|，
  与 Methods 公式一致）。
- 输出 `Figure7_calibration.pdf/png`。

---

## Supplementary figures

- **S1 计数关联（P1-C）**：散点 + loess + 边缘密度（仿 FigureYa37correlationV2
  / FigureYa163twoVarCor 语法）——donor n_cells(x, log) vs 预测分(y)，
  每 FM 一小图，标注 Spearman ρ（recompute from CSV; verified range 0.183–0.271，
  UCE 最低 0.183，均 p<0.005——不是"每个 FM 0.22–0.29"，这个旧锚点值有误，
  已在生成时用真实值覆盖）。`FigureS1_cellcount.pdf/png`。
- **S2 配对 ΔAUC 热图**（仿 FigureYa200pairwiseAUC 的矩阵语法）——
  within-COVID 全配对的 per-seed paired ΔAUC（数据
  `donor_level_pairwise_repeated.csv`），下三角热图 + 星号标显著。
  `FigureS2_pairwise_heatmap.pdf/png`。
- **S3 cohort 采样时间线**：已有成品
  `reports/manuscript/revision/literature_audit/T1T2_fm_corpora_and_cohort_windows/cohort_sampling_timeline.png`，
  校验即可；风格不统一再按同数据重绘。
- ~~S4 映射审计瀑布图~~ **已取消（2026-09-07 用户决定）**：病人数只在正文
  （Methods 一句归并规则）和 Figure 1（240/25）中说明，不做审计过程图。

## 验收清单

1. 全部点值与对应 CSV 重算一致（<0.0005）；F2/F3 与 P0 §3a/3b 完全一致；
2. F5/F6 无 legacy 数字；F6 只含 expanded_103pt；
3. 显著性只来自 pairwise CSV；无 bar chart 滥用（除允许的分组场景外，
   主图均为 lollipop/forest/dumbbell/violin/line 语法）；
4. PDF 矢量 + PNG ≥300dpi，Arial/Helvetica，疾病色与 Figure 1 一致；
5. `figures/v2/INDEX.md`：每图一行（文件、数据源、关键数字、脚本路径、
   仿用的 FigureYa 模板名）。

## 明确不做

不重画旧版图；不画 GSEA/SHAP/activation-axis 生物学图（需 v2 上重跑，
另立任务）；不把 legacy test_set 数字混进任何图；不改任何 CSV；不装 R。

---

# 修订附录 v2.1 — Job J–P 数据落地后的图更新清单（2026-09-08）

> 以下修改基于 Job J（Jiang 剔除）、Job L（配对推断终版）、Job N（LOBO v2）、
> Job K（QC 基线）、Job M（年龄分层）、Job O（matched 对照）、Job P（legacy 清除）。
> 规则同 v2：数字只从指定 CSV 重算；显著性只用 seed-ensemble 配对 bootstrap
> （Job L/Q 口径），**per-seed Wilcoxon 的 q 值全部作废**。

## F2（within-COVID）——显著性标注换推断口径

- 括注改为（数据源 `paired_inference/seed_ensemble_paired_bootstrap.csv`）：
  Geneformer − PCA = +0.064 [0.005, 0.126], **p=0.029**；
  STATE − PCA = +0.052 [0.003, 0.106], **p=0.038**；
  STATE − Geneformer：n.s.（p=0.60，写作 "could not be distinguished"）；
  scPRINT-2/scGPT vs PCA：n.s.；UCE vs PCA：marginal（p=0.077，不标星）。
- caption 增一句：contrasts use patient-level paired bootstrap (B=10,000) on
  the seed-ensemble out-of-fold scores。

## F3（within-sepsis）——显著性标注待 Job Q

- 当前 PCA 对 4/5 FM 的 q=0.0037 为 per-seed Wilcoxon（作废）。
  待 `paired_inference_full/within_sepsis_15pairs.csv` 落地后按 bootstrap p
  重标；落地前**先移除星号**，只留数值与 0.5 参照线。

## F4（monocyte）——换推断口径 + UCE 例外句

- 森林图 p 值待 `paired_inference_full/monocyte_vs_pca_5pairs.csv`
  （替换原 per-seed 的 STATE 0.002 / scGPT 0.037）。
- 图注加一句：UCE favors monocytes across all four configurations
  (seed × disease-pool); other models revert to the platelet-favoring
  direction (or zero) under the 10-seed ensemble。数据 `monocyte_matched/`。

## F5（cross-disease）——必须换 Job J 终版数字

- 数据源改 `cross_disease/cross_disease_patient_jobJ_20260907_154957.csv`
  （Jiang 已剔除）。**逆向锚点全换**：STATE-LR 0.799 [0.718, 0.876]、
  scPRINT-2-LR 0.695 / RF 0.701、Geneformer-LR 0.668、scGPT-LR 0.661、
  UCE-LR 0.653、STATE-RF 0.653。
- 正向保持 91-patient GSE279452 测试集数字不变（0.644–0.684，
  `cross_disease_chongqing_20260904_233237.csv`，Jiang 不在该测试集）。
- 5-patient GSE167363 的 1.000 AUC 等值**只进 supplementary note**（标注小样本），
  不进此图。

## F6（adaptation curve）——换 Job J 曲线数据

- 数据源改 `adaptation_curve_chongqing/`（test_set=expanded_96pt）：
  每方法 CSV（STATE/UCE/Geneformer/scPRINT-2/scGPT/PCA-200+LR）+
  `summary.csv` 校验。终点（n=all）全体 ≤0.62，仍无救援。
- ⚠ 校正（2026-09-08）：检验集 Jiang 剔除后实为 **32 patients**（非本文件早前
  版本写的 34）——已在图生成时按 CSV 实际值验证通过，仅本文档旧文字未同步，
  现已更正。

## S2（pairwise 热图）——显著性待 Job Q

- 15 对的星标待 `paired_inference_full/within_covid_15pairs.csv` 替换
  （现值为 per-seed Wilcoxon，作废）。数值（ΔAUC）不变，只换显著性层。

## 新增两张 supplementary（可选）

- **S-QC**：QC 分层分布（`qc_baseline/patient_qc_distributions.csv`，
  nUMI/nGene/%mito by outcome by disease；正文只引 QC-only AUC 0.528）。
- **S-Age**：GSE279452 年龄分布直方图 + 分层 fatal（median 67，19–93；
  <40: 0/6、40–60: 7/30、>60: 10/54[已知结局]，n=56 含 2 例 mort28 未知——
  见下方详规格的 ⚠ 校正）。

---

# 新增补充图详细规格（S-QC / S-Age，2026-09-08）

## Figure S-QC — QC 分层分布与 QC-only 基线（评审 A1 的图证）

**目的**：直观展示"QC 指标不能解释预测信号"——fatal 病人的 QC 有偏移
（nUMI/mito 更高），但 QC-only 模型只有 0.528 ≈ 随机。

**数据**：
- `results_patient_v2/qc_baseline/patient_qc_distributions.csv`
  （252 行：patient_id, y, disease, n_cells, median_counts, median_genes,
  median_pct_mito）
- `results_patient_v2/qc_baseline/qc_baseline_auc.csv`
  （QC-only AUC 0.5282±0.0260, 10 seeds）

**⚠ 口径核对（先做）**：`qc_baseline_auc.csv` 为 252 patients / 27 fatal
（COVID 240 + legacy sepsis 12）。评审 A1 要求的是 **within-COVID**。
请按 disease 列分层重报 within-COVID（240/25）的 QC-only AUC 一并进图注；
若无法分层，图注写 "pooled analysis cohort (252 patients, 27 fatal)"。

**布局（2×2 面板，violin+box+jitter，仿 FigureYa227boxdensity）**：
- A：median_counts（nUMI，log y 轴）按结局分组（fatal vs survived），
  按 disease 分面（COVID-19 / Sepsis）
- B：median_genes 同布局
- C：median_pct_mito 同布局
- D：QC-only AUC（0.528，带 repeat-sd 误差条）与 within-COVID 各模型 AUC
  参照点（取 P0 §3a 值）的 lollipop 对比 + 0.5 随机线；
  图注一句："Fatal patients show higher nUMI and %mito, yet QC alone does
  not predict outcome (AUC 0.53)."

输出 `FigureS_QC_distributions.pdf/png`。

## Figure S-Age — GSE279452 年龄分布与分层结局

**目的**：记录 sepsis 队列为纯成人组成及年龄-结局结构（评审 A2 的图证）。

**数据**：`results_patient_v2/monocyte_matched/gse279452_age_final.csv`
（92 行：donor_id, gsm, age, gender, mort28；中位 67，范围 19–93）

**布局（双面板）**：
- A：年龄直方图（10 岁 bin），fatal（mort28=1，橙）与 survivor（灰）
  堆叠或并排；标 median 67 虚线
- B：分层 fatal 率点图（<40 / 40–60 / >60 三层），每层标 n 与 fatal 数
  （<40: 0/6、40–60: 7/30、>60: **10/54**——该层共 56 例，其中 2 例
  mort28 未知，排除在 rate 分母外，图中需注明），点+误差条（二项 SE 或
  Wilson CI，分母用已知结局数）
  ⚠ 校正（2026-09-08）：早前版本写 "10/56" 把 2 例未知结局计入分母但未计入
  分子，会低估真实 rate（10/56=17.9% vs 正确 10/54=18.5%）；图生成时已按
  正确口径（分母=已知结局数）计算并在图内注明 n=56 含 2 例未知，此文字同步更正。

输出 `FigureS_Age_gse279452.pdf/png`。

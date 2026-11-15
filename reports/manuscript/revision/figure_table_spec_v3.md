# 图表重生成规格 v3(2026-09-09,给作图 agent)

> 依据:`revision/P0_onepager.md`(round-17 终稿,唯一数字权威)。
> 项目根:`/Users/xinruqiu/Desktop/302006/302006_updategit/`。
> 数据根(除注明外):`scripts/31_donor_level_revision/results_patient_v2/`。
> 输出到:`reports/manuscript/figures/v2/`(覆盖同名旧文件;旧文件先备份到 `v2/backup_v2_20260908/`)。

## 〇、总则(每张图/表都必须遵守)

1. **术语**:禁止裸 "AUC"——必须带限定(patient-level AUC / ensemble-OOF AUC / pooled LOBO AUC / OOF nested ΔAUC)。禁用 "OOD/out-of-distribution" 字样,改用 **underrepresented**;monocyte 对照组可写 "in-distribution (monocytes)"。图面分面标签用 "COVID-19" / "sepsis"(GSE 编号只进 caption 与正文)。
2. **两套 AUC 口径**:repeat-mean±repeat-sd(10 seeds 的稳定性视图)与 seed-ensemble OOF(10 个 repeat 分数按 patient 平均后单次评分)。凡出现推断(p 值)的图,一律用 ensemble 口径并标注。
3. **显著性口径**:pairwise 一律 exact DeLong;max/range 统计量一律 permutation。图注里写清检验名。
4. **风格**(参照用户给的 INDEX.html 美学):dot-and-whisker、forest/dumbbell、heatmap 为主,**不要全是 bar**;主色不超过 4 个;基线(null/cohort/QC/PCA-200)与 FM 用形状或灰度区分;网格线浅、无边框。
5. **矛盾并录不解释**:图注只陈述数值与检验,不写"这说明了什么"的解释句。

---

## 一、主图(5 张需动,1 张保留)

### Figure 2 — within-COVID 主图(重画,两面板)

**Panel A:六方法 + 三基线的 patient-level AUC(repeat-mean ± repeat-sd)**
- 数据:
  - 主表:`within_covid/donor_level_auc_repeated.csv`(mean/LR 行)
  - 数值(校核用):Geneformer 0.786±0.026、STATE 0.778±0.033、scPRINT-2 0.740±0.027、PCA-200+LR 0.725±0.033、scGPT 0.720±0.047、UCE 0.661±0.046
  - 基线:null(count+disease)0.554(`within_covid/null_model_auc_repeated.csv`);cohort-only 0.486±0.038(`paired_inference_full/cohort_only_baseline.csv` + `_per_seed.csv`);QC+cell count 0.618(as-fit,无 sd,来自 P0 §3a 基线表)
- 形式:横向 dot-and-whisker(点=repeat-mean,须=repeat-sd),按值排序;基线用灰色方块,FM 用彩色圆点,PCA-200 用空心圆。
- 标注:
  - 主检验:Geneformer 与 PCA-200 之间画括号,"pre-specified winner-corrected permutation P = 0.181 (n.s.)"
  - 双雄 vs UCE 画括号,"confirmatory range-statistic permutation P = 0.011"
  - chance 线 0.5 虚线。

**Panel B:增量价值(incremental value)**
- 数据:`paired_inference_full/incremental_value_test.csv`
  - 基线(QC+计数+队列)AUC=0.627;augmented:STATE 0.6919、Geneformer 0.6802、PCA-200 0.6772、scPRINT-2 0.6733、scGPT 0.6534、UCE 0.643
  - Δ 与 DeLong p:STATE +0.065 (p=0.00015)、Geneformer +0.053 (0.00014)、PCA-200 +0.050 (0.0034)、scPRINT-2 +0.046 (0.0054)、scGPT +0.026 (0.0063)、UCE +0.016 (0.072)
  - nested-LR P(成对展示):五者 P<0.0001,UCE P=0.0022
- 形式:森林图——每行一个表示,点=ΔAUC,须=1.96×DeLong SE(CSV 有 delong_SE 列);右侧两列文字标注 DeLong p 与 nested-LR P;UCE 行用灰色(两个口径都最弱)。
- caption 必须含:"baseline = QC metrics + cell count + cohort (AUC 0.627); increments over it by exact DeLong; nested likelihood-ratio P against an empirical permutation null; the pre-specified primary for this question is the PCA-200 row (nested-LR P < 0.0001)"。

### Figure 3 — within-sepsis(数据不动,只改图注)

- 保留现有 `Figure3_within_sepsis.pdf/png`(2026-09-08 版)。
- 图注/标注改为:"no pairwise contrast survives Benjamini–Hochberg correction; at 17 fatal events the benchmark resolves only AUC ≳ 0.71 (80% power) — a bound, not a proof of no signal"。

### Figure 4 — in-distribution 对照(重标 + 核对)

- 数据不变:`monocyte/donor_level_auc_repeated.csv`(STATE 0.696、scGPT 0.677、Geneformer 0.674、scPRINT-2 0.672、UCE 0.648、PCA 0.654)与 within-COVID 对照;显著性用 `paired_inference_full/monocyte_vs_pca_5pairs.csv`。
- **所有 "OOD" 字样改 "underrepresented (platelets)"**;"in-distribution (monocytes)" 保留。
- 标注:"no foundation model significantly beats the from-scratch baseline in-distribution (all exact DeLong p > 0.23); every method, including the baseline, scores lower on monocytes than on platelets (task property)"。

### Figure 5 — cross-disease(逆向重画 + 正向保留 + adaptation 加不可分辨带)

**Panel A(双向 dumbbell)**:
- 正向(91 测试池,不变):UCE-LR 0.684 [0.543, 0.805]、STATE 0.653、scGPT 0.650、scPRINT-2 0.647、Geneformer 0.644;数据 `cross_disease_chongqing/cross_disease_chongqing_20260904_233237.csv`。
- 逆向(**96 训练池,更新**):STATE-LR 0.7209 [0.621, 0.809]、Geneformer-LR 0.7118 [0.612, 0.814]、UCE 0.6928、scGPT 0.6417、scPRINT-2 0.6186;RF:STATE 0.7113、Geneformer 0.6358、UCE 0.6074、scGPT 0.6221、scPRINT-2 0.5167。数据 `cross_disease/cross_disease_patient_96pool_20260908_142954.csv`。
- 标注:"forward = uninformative (mortality not learnable in the sepsis pool even in-domain); reverse = rank-indeterminate (top-2 reversed under severe-excluded; UCE shifts +0.086 with 5 added training patients)"。

**Panel B(adaptation curve)**:
- 数据不变:`adaptation_curve_chongqing/summary.csv` 与各模型 csv(n_sepsis_req = −1/0/10/50/100/500/1000,point_auc_mean 与 boot CI)。
- **加灰色横带 [0.15, 0.85]**:标注 "at ~3 fatal events in 32 test patients the 95% interval spans ≈0.7 — unresolvable, not null"。

### Figure 7 — calibration(保留)

- 保留 `Figure7_calibration.pdf/png`(2026-09-07)。数值核对:LR 中 STATE Brier 0.118/ECE 0.168;STATE-RF 0.095/0.079;RF 全体 ECE 0.054–0.080。

---

## 二、补充图(2 张重画,1 张新建,其余保留)

### FigureS2 — within-COVID pairwise 热图(**重画,补全 6 模型**)

- **上一版缺 STATE 与 PCA-200,这次 6×6 全方阵**(STATE、UCE、Geneformer、scPRINT-2、scGPT、PCA-200)。
- 数据:`paired_inference_full/delong_exact.csv`(15 对 ΔAUC、delong_exact_SE、delong_exact_p)。
- 形式:行−列 ΔAUC 热图(发散色,中心 0),格内标 Δ 与 p;p<0.05 加粗边框,过 BH(m=15)加星(**共 6 对**:UCE−Geneformer、STATE−UCE、UCE−scPRINT-2、Geneformer−scGPT、UCE−scGPT、STATE−scGPT——全部以 UCE 或 scGPT 为输家;leaders−UCE 两格为最大差距)。
- caption:"pairwise differences in ensemble-OOF AUC, exact DeLong; six contrasts survive Benjamini–Hochberg (q ≤ 0.037), all with UCE or scGPT on the losing side, the two leaders-versus-UCE gaps being the largest; the leaders-versus-UCE separation is the paper's confirmatory result (range-statistic permutation, P = 0.011)"。

### FigureS-LOBO — LOBO 面板(新建)

- 数据:`lobo_v2/lobo_per_cohort_auc.csv`(5 cohort × 5 FM 的 held_out_auc、n_patients、n_fatal)、`paired_inference_full/lobo_pooled.csv`(拼接 pooled:STATE 0.7525、scPRINT-2 0.7514、Geneformer 0.7481、UCE 0.6813、scGPT 0.6376)、`paired_inference_full/lobo_pooled_per_cohort_n.csv`。
- 形式:左面板=每 FM 一排,5 个 cohort 点(点大小∝fatal 数)+ 一条 pooled 横线;右侧面板=per-cohort n/fatal 小表。
- caption:"LOBO is outcome-selected: folds require ≥1 fatal and ≥1 survivor in the held-out cohort; survivor-only cohorts (34 patients) train only. Pooled AUC over concatenated held-out predictions (206 patients, 25 fatal)."

### FigureS-null — 置换零分布(新建,**数据已到位,必做**)

- 数据:`paired_inference_full/primary_and_range_permutation_null_vectors.csv`(两列:`max_delta_vs_pca` 与 `range_leaders_vs_uce`,各 10,000 个置换值;已验证再生 P 与正文一致:0.181 / 0.0112)。
- 形式:并排两个直方图(bin≈50,可叠加 KDE):左=primary max-statistic 零分布,竖线标观测 0.0644,注 "P = 0.181";右=range 零分布,竖线标观测 0.1302,注 "P = 0.011"。
- caption:"null distributions from 10,000 patient-level label permutations on the frozen seed-ensemble scores; vertical lines mark the observed statistics (max FM−baseline ΔAUC = 0.064, P = 0.181; max−min AUC across the five foundation models = 0.130, P = 0.011)"。

### 保留(不动):FigureS1_cellcount、FigureS3_cohort_timeline、FigureS_Age、FigureS_QC、FigureS3_dropout(若已存在)。

---

## 三、主文表(3 张新建/重建,其余已存在)

### Table 3 — within-COVID 主表(重建,双口径)

列:模型 | repeat-mean ± repeat-sd | ensemble-OOF AUC | vs PCA-200(exact DeLong p)| vs null(exact DeLong p)。
- repeat-mean 源:`within_covid/donor_level_auc_repeated.csv`
- ensemble 源:`paired_inference_full/winner_test_ensemble.csv`(STATE 0.7946、UCE 0.6765、Geneformer 0.8067、scPRINT-2 0.7561、scGPT 0.7392、PCA-200 0.7423)
- vs PCA 源:`paired_inference_full/delong_exact.csv`(Geneformer 0.036、STATE 0.048、scPRINT-2 0.65、scGPT 0.91、UCE 0.081)
- vs null 源:`paired_inference_full/each_model_vs_null.csv`(STATE 0.0014、Geneformer 0.0006、scPRINT-2 0.011、scGPT 0.023、PCA-200 0.013、UCE 0.157)
- 基线三行:null 0.554、cohort-only 0.486±0.038、QC+count 0.618。
- 表注:预声明主检验 permutation P=0.181(cross-check 0.073),双雄−UCE range P=0.011。

### Table 5 — transfer 表(重建)

列:方向 | 模型 | LR AUC [95% CI] | RF AUC | severe-excluded LR AUC。
- 逆向 96 池(主)+ 91 池(敏感性行):源 `cross_disease/cross_disease_patient_96pool_20260908_142954.csv` 与 `cross_disease_chongqing/cross_disease_chongqing_20260904_233237.csv`。
- severe-excluded:STATE 0.765、Geneformer 0.775(其余从 96pool CSV 的 sens_* 列取)。
- 正向 91 池同 CSV。

### Table 8 — 增量价值表(新建)

列:模型 | augmented AUC | ΔAUC | exact DeLong SE | DeLong p | nested-LR P。
- 全部来自 `paired_inference_full/incremental_value_test.csv`;nested-LR P 列:五者 P<0.0001、UCE P=0.0022(P0 §3a);基线行 AUC=0.627 单独列在表头。
- 表注必须写:"the pre-specified primary is the PCA-200 row (no model selection); the five FMs are exploratory with both readings paired; nested-LR establishes existence, ΔAUC quantifies magnitude"。

### Table 9 — 综述对照表(新建,进 Introduction)

按 P0 §5 的三行(EGSP / 9-模型 benchmark / TissueFormer),列:研究 | scFM | 单细胞输入方式 | 任务 | 队列 | 事件数 | from-scratch 基线 | 病人级聚合 | 校正 | 备注(embedding 迁移/预印本)。
- 数字照 P0 §5 表誊(EGSP:8,776 样本/3,095 death;9-模型:1,170 supervised 实验,事件数未按任务披露;TissueFormer:376 donors,severe 46/83/28)。
- 表注:搜索策略一句(PubMed+bioRxiv,查询串与日期 2026-09-09,4 召回 2 纳入 2 排除;TissueFormer 为定向复核纳入);并注明"三项均未报告 selection/multiplicity 校正,且无一按比较粒度披露事件数"。

---

## 四、补充表(数据就位,直接生成)

| 表 | 源 |
|---|---|
| tab:score_correlation | `paired_inference/score_correlation_pearson.csv` |
| tab:winner_ensemble | `paired_inference_full/winner_test_ensemble.csv` |
| tab:se_check(DeLong vs bootstrap)| `paired_inference_full/delong_exact.csv` + `paired_inference_full/within_covid_15pairs.csv` |
| tab:lobo_folds | `lobo_v2/lobo_per_cohort_auc.csv` + `paired_inference_full/lobo_pooled_per_cohort_n.csv` |
| tab:qc_strat | `qc_baseline/patient_qc_distributions_covid_only.csv` |
| tab:fold_balance | `paired_inference_full/fold_cohort_balance_seed42.csv` |
| tab:qc_by_cohort | `paired_inference_full/qc_by_cohort_association.csv` |
| tab:cohort_rates | `paired_inference_full/item4_per_cohort_table.csv` |
| tab:severe_sens | `cross_disease/cross_disease_patient_96pool_20260908_142954.csv` 的 sens_* 列 |
| tab:aggregation_sens | `within_covid/donor_level_auc_repeated.csv` 的 mean/median 行 |
| tab:cellcount | `within_covid/cellcount_association.csv` |
| tab:incremental 诊断(univariate/coef signs)| `paired_inference_full/item3_univariate_association.csv`、`item3_multivariate_coef_signs.csv` |

**tab:nested_oof(3×5)与 tab:nested_deviance:真实 CSV 已回传(2026-09-09 23:06),直接生成**:
- nested_oof 源:`paired_inference_full/round16_delong_scaled_comparator.csv`(三个 reduced_config × 5 FM,含 auc_reduced/auc_full/delta_auc/delong_se/delong_p/pred_corr;三个 reduced AUC=0.7098/0.7241/0.7131,BH 全部 0/5 通过)。BH rank-by-rank 阈值表按此文件 p 列现算。
- nested_deviance 源:`paired_inference_full/round12_headline_deviance_perm_task0..3.csv`(4 task 各 2,500 次,合并 10,000 个 per-permutation max)+ `round12_headline_deviance_obs.csv`;零分布百分位 median 0.188、95th 0.790、99th 1.230、max 3.819,观测 STATE 4.3284、Geneformer 3.2227、scPRINT-2 3.1843、scGPT 1.5365、UCE 0.8749;表注写 "in-sample cross-check, carries no headline weight"。

---

## 五、BioRender 两份(描述给 BioRender 或手绘)

### Figure 1(更新现有版本)

保持现有布局(队列 → 统一输入 → 5 模型 frozen embedding → 患者级协议 → 四个分析臂),只改文字:
- 标题/标签:"underrepresented platelets"(删掉一切 OOD)。
- 队列数字:COVID-19 240 patients/25 fatal;sepsis 91/17(GSE279452,28 天)+ 5(GSE167363);GSE151263 7 人无结局(灰显,标注 "excluded from all mortality analyses");healthy 70(reference only);monocyte 对照 148/19。
- 协议:10×5 repeated-CV,patient 分组;五个固定 AUC 术语。
- 四臂:within-COVID / within-sepsis / cross-disease(双向)/ monocyte 对照。

### Graphical Abstract(重做,全-null 框架)

三幕,左→右:
1. **设计**:5 个 scFM(图标)+ 同一 raw platelet 输入 + "identical preprocessing, frozen embeddings, patient-level repeated-CV";侧标 from-scratch PCA+LR 基线。
2. **结果**:大数字 "0"——no FM beats baseline(pre-specified permutation P = 0.181;LOFO 25/25);唯一存活分离:双雄 vs UCE(range P = 0.011,标注 "about the weakest model");sepsis/transfer/monocyte 三个 null 小图标。
3. **标尺**:横轴 = 事件数,标 "25 events → resolves ~0.10 AUC" 与 "0.03–0.05 AUC needs ~10³ patients";底句:"for patient-level outcome prediction at this scale, a from-scratch baseline suffices"。

---

## 六、明确不做

- 不做 S4 mapping audit 图(用户已拍板:映射只在正文与 Figure 1 文字说明)。
- 不重画 Figure3(sepsis)、Figure7(校准)、FigureS1/S3/S-Age/S-QC。
- 不做任何含 "OOD" 字样的素材;不做 3D、不做双 y 轴。

---

## 七、Figure: gene attribution + GSEA(masking 线程,2026-09-11 新增;编号待 P2 组装)

脚本:`code/figures_v2/fig_attribution_gsea_v1.py`(stdlib+numpy+matplotlib,本地 python3.12 可跑);
数据:`results_gene_attribution_gsea/<model>_masking/` ranked CSV + 3 GSEA CSV × 5 模型,
`cross_model/`(由 `code/compare_masking_attribution_5models.py` 生成:Spearman/top-N Jaccard/
consensus/sign-flip/NES 矩阵/显著项 Jaccard)。
输出:`figures/v2/Figure_gene_attribution_gsea.{pdf,png}`。

- **Panel A**:五模型各自 top-7 基因的并集(34 行),格 = 该基因在该模型 list 内的 rank
  percentile(viridis_r,0=top,0.5+ 截断),+/− 标 attribution 符号;缺失(不在该模型 list)灰格。
  **只比秩不比数值**(各模型 logit 尺度不同,across-gene SD 0.012–0.061)。
- **Panel B**:GSEA 通路×模型 NES 矩阵(RdBu_r,●=FDR q<0.25):共享行(hallmark 5 + kegg 3 +
  reactome 3,按 n_models_sig 与 max|NES| 排序)+ 分隔线 + 每模型最强 model-specific 一项
  (scPRINT-2 无显著项,缺)。scPRINT-2 整列无 ● = 全 null 的直接可视。
- **校核锚点**(脚本内打印):显著项计数 = GF 6/23/2、STATE 5/7/2、scGPT 28/63/257、
  UCE 23/70/82、scPRINT-2 0/0/0(H/K/R);共同基因 3,917;Spearman −0.028..0.401;
  ≥4/5 模型共享通路 = 0 条。
- **脚注**:in-sample 分类器(归因方向,非 OOF 性能);scPRINT-2 null 经 N=10/30/128 联合遮蔽
  剂量对照判定为架构遮蔽不敏感(冗余编码已排除);UCE 用 bool-mask 修复版。

# BioRender 作图说明 — Figure 1（Study design overview，终版 v4，2026-09-09）

> 用途：绘制并替换 `reports/manuscript/figures/Figure1_overview.pdf`。
> 图面三原则：**只讲现在的设计**（无任何历史/版本痕迹）、**队列人数与文字结论徽章
> 可以出现，但不出任何 AUC 数值**、**不出现数据集编号**（正文用 GSE accession，图面
> 只写 COVID-19 / sepsis / healthy / monocytes）、**禁用 "out-of-distribution / OOD"**
> （用 underrepresented）。
> 输出：PDF（矢量）+ PNG（≥300 dpi），命名 `Figure1_overview_v4.*`。

## 整体布局

横向 2×2 四段（A 左上 → B 右上 → C 左下 → D 右下），每段一个带底色圆角框，
粗箭头 A→B→C→D。小标题：
**A. Cohorts ｜ B. Unified input, models & controls ｜ C. Patient-level
evaluation ｜ D. Four analysis arms**。

## Panel A — Cohorts（左上，淡蓝底）

- 图标流（左→右）：采血管 → 无核血小板图标 → 液滴 scRNA-seq 图标；
  下方标注 "Blood draw / anucleate platelets / single-cell RNA sequencing"。
- 四个病人徽章（人形 icon + 两行字，横向排列）：
  - **COVID-19**（红 #C0392B）："240 patients" / "25 fatal"
  - **Sepsis**（橙 #E67E22）："103 patients" / "19 fatal"，徽章下挂一个灰色
    小脚注框："7 without outcome records — excluded from all mortality analyses"
  - **Healthy**（灰绿 #7F8C8D）："70 donors" / "reference only"
  - **Monocytes**（紫 #8E44AD，虚线侧支框）："148 patients / 19 fatal" /
    "in-distribution positive control"
- 徽章下一行小字："Multiple draws per patient merged to patient level"
- Panel 底部加粗一行："**61,759 platelet transcriptomes** across the COVID-19
  and sepsis analysis pools"

## Panel B — Unified input, models & controls（右上，淡橙底）

- 顶部数据文件图标 + 文字："**Identical raw input (one AnnData)**"
- 向下分 5 条平行箭头，各落一个模型框（每框：模型名加粗 + "native
  preprocessing, frozen embeddings" + embedding 矩阵小图标）：
  **STATE (SE-600M) ｜ UCE ｜ Geneformer ｜ scPRINT-2 ｜ scGPT**
  （注意最后一个拼写是 **scGPT**，s-c-G-P-T）
  STATE 框右上角加一个小注标签："embedding component of a two-component system"
- 底部两条虚线框（baseline / controls）：
  - "from-scratch baseline: **PCA-200 + logistic regression**"
  - "control families: **null** (cell count + disease) · **cohort-only** ·
    **QC-only** · **incremental-value** (QC + count + cohort ± representation)"

## Panel C — Patient-level evaluation（左下，淡青底）

- 上半：repeated-CV 示意——一排 ~10 个小圆（10 个 seeds），每个圆内 5 个扇区
  （5 folds）；关键视觉：**同一病人的人形（同色）始终落在同一折内、绝不跨折**。
  下方标注："**Repeated patient-level CV (10 seeds × 5-fold, patients grouped)**"
- 中部流程：一堆彩色细胞点（"cell-level predictions"）→ 聚合漏斗（mean）→
  单个分数条（"one score per patient"），下挂小字 "every prediction aggregated
  to the patient before scoring"
- 下半：推断层级三阶徽章（自上而下）：
  1. "**Pre-specified primary** — best FM vs baseline, winner-corrected
     permutation test"
  2. "**Confirmatory secondary** — range statistic"
  3. "**Exploratory** — exact DeLong, BH within each family"
- 底部醒目边框徽章："**Primary metric: patient-level out-of-fold ROC-AUC**"

## Panel D — Four analysis arms（右下，淡粉底，2×2 小格）

每格：迷你图标 + 臂名 + 队列行 + **一句文字结论徽章**（无数字）：
1. **Within-COVID-19 mortality (primary)** — "240 patients / 25 fatal"
   （下行趋势线 + 星标）——徽章："learnable by most methods; no FM beats the
   baseline (n.s.); one separation: two leaders > UCE"
2. **Within-sepsis mortality** — "91 patients / 17 fatal"
   （平线图标）——徽章："complete null at late windows; resolution bound stated"
3. **Cross-disease transfer** — "COVID-19 ⇄ sepsis"（双向箭头）——徽章：
   "uninformative forward; rank-indeterminate reverse"
4. **Monocyte positive control** — "in-distribution"（单核细胞图标）——徽章：
   "no FM beats the baseline there either"
- Panel 底部通栏小字："Adaptation curve, calibration, LOBO, dropout →
  supplementary"

## 风格

- 字体：无衬线（Arial/Helvetica）；标题加粗，数字与正文同字号。
- 疾病色固定：COVID 红 / Sepsis 橙 / Healthy 灰绿 / Monocyte 紫，与正文数据图一致。
- 禁止出现：AUC 等结果数字、GSE 编号、"legacy / original / corrected / new /
  old" 等版本词、"OOD / out-of-distribution"、"413 donors / 47,977 / 64,356"。

## caption（进 main.tex 的终版）

**Figure 1. Study design overview. (A)** Cohorts: platelet single-cell
transcriptomes from 240 COVID-19 patients (25 fatal), 103 sepsis patients
(19 fatal; 7 without outcome records are excluded from all mortality
analyses), and 70 healthy donors (reference only); multiple draws per patient
are merged to patient level. Monocytes from the same source cohorts (148
patients, 19 fatal) serve as an in-distribution positive control. **(B)**
Unified input: every model receives the identical raw AnnData object and
applies its own native preprocessing; five foundation models (STATE via its
state-embedding component, UCE, Geneformer, scPRINT-2, scGPT) are benchmarked
against a from-scratch PCA-200 + logistic-regression baseline and null,
cohort-only, QC-only, and incremental-value control families. **(C)**
Patient-level evaluation: repeated cross-validation (10 seeds × 5-fold) with
patients grouped so that a patient's cells and draws never span folds; all
predictions are aggregated to the patient before scoring; the primary metric
is patient-level out-of-fold ROC-AUC, with one pre-specified winner-corrected
primary test, one confirmatory secondary test, and exploratory contrasts
controlled within each family. **(D)** The four analysis arms: within-COVID-19
mortality (primary), within-sepsis mortality, bidirectional cross-disease
transfer, and the in-distribution monocyte control; the adaptation curve,
calibration, and robustness probes appear in the supplement. Figure created
with BioRender.com.

---

# BioRender 作图说明 — Figure 2（Cohort assembly and analysis sets，流向图）

> 用途：回应评审 B5——在**正文**放一张 CONSORT 式流向图（GSM → samples →
> patients → 各分析集）。与 Table 1（cohort 特征表）互为图形/表格互补。
> 图面规则同 Figure 1：无结果数字、无历史/版本词、不出现 AUC。

## 整体布局

纵向三段（上=Sources → 中=Samples/Patients → 下=Analysis sets），主流程粗箭头
自上而下；剔除分支走右侧虚线。

## 顶层 — Sources（4+1 框，横向）

- "**9 COVID-19 cohorts**"（红边）
- "**GSE167363** (sepsis)"（橙边）
- "**GSE151263** (sepsis)"（橙边）
- "**GSE279452** (sepsis)"（橙边）
- 独立灰框（虚线，不与主流程相连）："**Healthy donors — 70, reference only**"

## 中层 — Samples → Patients（两个主流程框 + 两个剔除分支）

- 左主框（蓝底）："**322 Sample_IDs → 252 patients**"，下方小字
  "multiple draws per patient merged to patient level"
- 右主框（蓝底）："**93 adult sepsis donors → 91 patients**"，下方小字
  "with 28-day outcome labels"
- 右侧两个带垃圾桶图标的虚线框：
  ① "**3 pooled libraries excluded** — 2,597 cells"
  ② "**2 donors excluded** — missing 28-day outcome"
- GSE151263 单独虚线框（橙虚线边）："**7 patients — dataset only, excluded from
  mortality analyses** (no per-patient outcome)"；箭头只指向"数据集"说明，
  不进入任何下方分析框。

## 底层 — Analysis sets（4 个彩色终点框，横向）

1. 红框："**Within-COVID-19 mortality (primary)** — 240 patients / 25 fatal"
2. 橙框："**Within-sepsis mortality** — 91 patients / 17 fatal"
3. 双向箭头框："**Cross-disease transfer** — COVID-19 (240) ⇄ sepsis (96)"
4. 紫框："**Monocyte positive control** — 148 patients / 19 fatal"

每框上角标一行极小字注明输入来源（如红框 "from COVID-19 cohorts"，
橙框 "from GSE279452"，双向框 "both directions, patient-level"）。

---

# BioRender 作图说明 — Graphical Abstract（投稿系统用，单幅横版）

## 画布与布局

单幅 16:9（建议 1200×675 px），三段式左→右：The test → The protocol →
The answer，浅色底。

## 左段 — The test（约 30% 宽）

- 中央一个大的无核血小板图标，旁标小字："**a cell type the models never saw**"
- 血小板上方 5 个模型名小徽章（STATE / UCE / Geneformer / scPRINT-2 / scGPT），
  下方两个虚线小框："PCA + LR" 和 "null"
- 底部一行小字："64,356 platelet transcriptomes · 343 patients · COVID-19 + sepsis"

## 中段 — The protocol（约 25% 宽）

- 一个圆形循环箭头内含 10 个小点（= 10 seeds），环绕一个 5 扇区饼
  （= 5 folds），下方标注："**repeated patient-level CV**"
- 流程箭头：彩色细胞点 → 漏斗（Σ）→ 单个人形，标注 "one score per patient"

## 右段 — The answer（约 45% 宽，2×2 结论格）

每格一个小图标 + 不超过两行字：
1. 绿色勾 — "**Within COVID-19**: top models beat simple baselines (AUC ~0.78)"
2. 灰色空圆 — "**Within sepsis (late windows)**: no model above chance"
3. 断裂双向箭头 — "**Cross-disease**: no transfer, no rescue"
4. 天平 — "**Conclusion**: pretraining value is *conditional & model-specific*"

## 底部横幅（整宽，深色条带白字）

"**Stress-testing foundation models where they were never trained — some win,
most don't, and it depends on the model, the task, and the protocol.**"

## 风格规则（三张图通用）

无衬线字体；除指定短句外不加文字；不出现表格、坐标轴、p 值、AUC 数字
（Graphical abstract 的 "AUC ~0.78" 除外）；疾病色仅用于结论格 1、2 的小图标
与 Figure 1 的队列徽章：COVID 红 #C0392B / Sepsis 橙 #E67E22 / Healthy 灰绿
#7F8C8D。

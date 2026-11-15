# Sepsis 臂扩展：GSE279452（Chongqing sepsis）整合评估

> 来源：用户新找的数据（`Downloads/sepsis-fatal/`，已归档到
> `literature_audit/sepsis_fatal_scout/`）+ 本人对 GEO 的一手核查（2026-09-03）。

## 数据集结论（一手核实）

**GSE279452**（Nat Immunol 2025, doi:10.1038/s41590-025-02345-x；重庆医科大学附属第一医院，
East Asian 队列）：

- **规模**：133 个 scRNA 样本（PBMC，10x MTX 矩阵，开放下载）+ 40 CITE-seq；
  成人 sepsis **93** 个 scRNA donors（呼吸 37 / 腹腔 36 / 泌尿 9 / 皮肤 8 / CNS 3；
  全队列口径下另有成人 ICU 对照 26、健康对照 43，但两者中分别只有 8 与 7 个
  做了 scRNA-seq——`cohort=adult`裸过滤会把这 8+7 个也误纳入"成人 108
  scRNA donors"，需按 `group=='Adult Sepsis'`排除）、儿童队列（**不纳入**）。
- **结局**：supplement 逐 donor 给出 14/28/90 天死亡率（`chongqing_sepsis_donor_mortality.csv`，
  295 行已抽好）。成人 sepsis 93 donors：**mort28 = 17 dead**（2 缺失 → 用于
  mortality 任务的队列 = **91 donors**），mort90 = 28 dead（3 缺失），
  mort14 = 9 dead（1 缺失）。
- **platelet 闸门 ✅**：实测一个成人 sepsis 样本（GSM8571191 Res-S22，7,814 细胞）：
  血小板标志物（PPBP/PF4/GP9/ITGA2B 等 12 个）≥2 counts 的细胞 **18.3%**（1,432），
  ≥10 counts 的 **10.5%**（822）。与 Qiu 2021 的 PBMC 类似存在 platelet–monocyte
  聚合体，注释流程需与现有整合一致（聚类+marker，不是裸计数）。
- **GEO 对接**：样本标题即 supplement 的 Donor ID（246/266 逐字匹配）→ 结局可直接 join。

## 对 benchmark 的意义

| | 现在 | 加入后（成人 scRNA） |
|---|---|---|
| sepsis patients | 12（Qiu 5 + Jiang 7） | **103**（+91） |
| sepsis fatal | 2 | **19**（28 天死亡口径，2+17） |
| within-sepsis CV | 不可做 | **可做**（17 fatal / 91） |
| cross-disease 测试集 | 描述性 | 有真实 n |

→ sepsis 结论可从 "descriptive only" 升级；COVID→sepsis transfer 变成有统计意义的检验。

## 必须先定的两个决策点

1. **结局口径**：该队列是 28/90 天死亡率，不是 in-hospital mortality。
   建议：主口径 = **28 天死亡**（sepsis 标准终点），Methods 一句话注明与
   in-hospital 的差异；mort90 作敏感性。mort28 缺失 2 例按删失处理（剔除或标 NA）。
2. **队列定位**：作为**第三个 sepsis 来源队列**并入（新 Data_NO，与 Qiu/Jiang 并列），
   不混合 donor；healthy/ICU 对照不纳入 mortality 任务（与现有一致）。

## 未关闭的闸门（整合前必查）

- **采样时点**：GEO characteristics 无采血时间字段——需查论文/supplement
  （mortality 表来源的同一个 supplement 大概率有）。Adam 的"sepsis 时点 critical"
  同样适用于新数据。→ 追加给 claude science（Task 5）或我直接查。
- **细胞注释方案**：与现有整合相同的 platelet 判定流程（聚类 + marker 模块分），
  并报告 platelet–monocyte 聚合体比例（双保险：marker 计数 + 聚类）。

## HPCC Job S（新 arm，排在 Job H2 之后、Job F 之前视资源而定）

1. 下载 91 个成人 sepsis scRNA MTX 样本（~4 GB）+ mortality CSV 对接；
   跳过 VDJ/CITE/儿童/对照。
2. 逐样本 QC → 合并 → 与 option-C 相同 schema（Sample_ID=Donor ID，
   Disease=Sepsis，covid_severity：fatal(mort28=1)/severe-survivor，
   Data_NO=chongqing2025）。
3. platelet 注释（同现有流程）→ platelet 细胞子集 h5ad。
4. 5 个 FM embedding（GPU；platelet 细胞量估计 20–80k，几小时级）。
5. 评估（全部 patient-level，沿用 v2 规则）：
   - within-sepsis mortality CV（91 patients / 17 fatal）——新主结果之一；
   - COVID→sepsis transfer（测试集 = 本队列）；
   - sepsis→COVID transfer（训练集含本队列）；
   - adaptation curve 重跑（sepsis 训练池从 12 → 103 patients）。
6. 报告：每 donor platelet 细胞数分布、fatal 患者清单与 mort28 对接表、
   上述四类数字。

## 风险与说明

- 28 天死亡 vs in-hospital 的口径差会让"fatal"标签在两病种间定义不完全对齐
  （COVID 侧是 in-hospital）——如实写，不辩护；必要时加 in-hospital-only 敏感性
  （若 supplement 有住院结局列）。
- 单中心、单族群（East Asian）队列——写进 cohort 表与 Limitations，不辩护。
- platelet–monocyte 聚合体：与现有数据相同的处理方式 + 在 cohort 表注明。

---

## 采样时点核查（Task 5，2026-09-03 已回）

**结论：这不是一个早期 sepsis 队列——多数采血在住院 1–4 周后。**

`gse279452_donor_sampling_dates.csv`（295 行，逐 donor 采血日历日期 + `hosp_days`）：
- 成人 sepsis scRNA 且带日期的 93 donors：`hosp_days` min 3 / **median 16** /
  max 121（Q1=9, Q3=26, P90≈56）
- ≤3 天仅 2 人；4–14 天 39 人；**>14 天 52 人（56%）**
- **fatal（m28=1, 17 人）采血更早**：median 10 天 vs survivor median 18.5 天
  ——结局相关的时点偏倚，within-sepsis 分析需当作混杂处理（同细胞数混杂的
  处理方式：报告 + 分层敏感性，不解释）
- 采血日历年跨 2021-11 → 2022-11；每 donor 单次采血

⚠ 两个待办：
1. **`hosp_days` 语义**（入院→采血的天数 vs 总住院时长）需从 supplement 一句话
   确认；字段位置（与采血日期同列）强烈支持前者。
2. 日期文件 93 donors vs mortality 文件 108 scRNA donors 有 15 人出入
   （m28 dead 17 vs 20）——以 mortality 表为标签权威，日期仅用于时点描述；
   对接时按 donor_id join 后报告缺失数。

**对论文口径的影响**：Chongqing 的时点异质性与 COVID 各队列（1–126 天窗口）
同级，可被同一句 Limitations 覆盖；但与 Qiu sepsis（T0/T+6h）时点锚完全不同
——cohort 表"采样窗口"列照实写 "median 16 d from admission (range 3–121)
[若 hosp_days 语义证实]"，不辩护。

---

## Job S 首批结果（2026-09-04，未含三项修复；scGPT 在途）

**within-sepsis CV（91 patients / 17 fatal，repeated-CV，mort28）：**

| 模型 | repeat-mean±sd | seed-42 |
|---|---|---|
| UCE | 0.548±0.038 | 0.595 |
| scPRINT-2 | 0.523±0.049 | 0.568 |
| Geneformer | 0.515±0.043 | 0.542 |
| STATE | 0.511±0.040 | 0.556 |
| PCA-200 | 0.428±0.082 | 0.540 |

- **所有 FM 基本在随机水平；PCA-200 显著低于 4/5 FM（STATE/UCE/Geneformer/
  scPRINT-2，q=0.0037），与 scGPT 不可分（q=0.40）；UCE 显著领先 STATE
  （q=0.010）——近随机水平下的组间差异如实报**
  ——FM > PCA 显著（与 within-COVID 故事方向相反），但绝对水平无用。
- COVID→Sepsis transfer（Chongqing 测试集）：STATE ≈0.58，vs 旧 12-patient 数字 0.80+
  ——两个候选解释（旧测试集 small-n artifact / Qiu-T0 与 Chongqing-d16 时点错配），
  按规则并录、不裁定。
- 混杂已记录：fatal 采血更早（median 10 vs 18.5 天）。

**对框架的影响（写作口径）：**
1. within-COVID 是唯一的真实信号来源；sepsis 28 天死亡从这些采血时点下不可预测。
2. "可转移的炎症-死亡生物学"假说被扩展数据拒绝——这**强化**了 scFM 批判主线
   （FM 不转移），同时削弱 sepsis 侧的临床 mortality 角度。
3. adaptation curve 按原计划跑（见下）：它是唯一还没回答的问题——
   "有目标疾病数据后，跨病迁移能否被救回来"。

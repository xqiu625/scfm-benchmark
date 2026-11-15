# 七模型 Architecture / Training Objective / Learning Paradigm 分析

> 用途：Discussion 的模型表 + 一段维度解读。事实层全部经原文核实（2026-09-06，
> 核实报告含逐字引用与页码，存档于本目录）：
> STATE=bioRxiv 2025.06.26.661135v2；scPRINT-2=bioRxiv 2025.12.11.693702v3；
> UCE=Nature 2026 (s41586-026-10689-z)；STACK=bioRxiv 2026.01.09.698608v2；
> TranscriptFormer 见 `docs/models/transcriptformer-best-practices.md`；
> Geneformer V1/scGPT 原文已核实；Geneformer V2-316M 规格暂依 HF model card（见文末清单）。
> 解读层纪律：每属性 n=1，只许 "consistent with"，不写 "explains"。

## 档案表（事实层，已核实）

| 模型 | Architecture | Training Objective | Learning Paradigm / 预训练数据 |
|---|---|---|---|
| **STATE** | **dense bidirectional transformer**（SE encoder：16 层、16 头、h=2048、600M 参数）。我们用 **SE 模块 frozen embedding（2,058 = 2048 CLS + 10 dataset-proj，Eq.31）**。⚠ `main.tex:268` 旧描述 "state-space architecture" 错误，已改 | SE：双轴重建损失（细胞内逐基因 L2 + 跨细胞群体一致性，Eq.33–35）+ 辅助 dataset-classification 损失（Eq.36–38）；ST：扰动状态转移（cell sets 上 bidirectional self-attention + MMD loss） | 自监督；SE=167M 观测人体细胞（scBaseCount 71.7M + CxG 59.2M + **Tahoe-100M 36.2M**，Table 2）；ST=>100M 扰动细胞（Tahoe-100M 50 癌细胞系、Replogle 4 细胞系、Parse-PBMC 9.7M **PBMC**） |
| **Geneformer** | BERT 式 encoder。**V1（2023 Nature）：6 层、4 头、256-dim、输入 2,048**（~10M 参数级；细胞嵌入 = 倒数第二层基因嵌入均值）；**V2-316M（我们 benchmark 用）：18 层、1,152-dim、~316M 参数、~20K 词表、input 4096**（HF 官方页权威："Geneformer-V2-104M and Geneformer-V2-316M: updated model trained Dec 2024 on ~104M human single cell transcriptomes… The current default model is Geneformer-V2-316M"；并经实际 checkpoint safetensors 配置核实；存储目录名 `geneformer-12L-30M` 为早年遗留误命名，三处 embedding 同一权重，无混用；V1-10M 仅用于 LoRA 臂） | **rank-value** 掩码基因建模（mask 15%，细胞内基因按全语料归一化表达排序成 token，丢弃绝对计数；原文明确此设计"more robust against technical artefacts"——为我们 batch 鲁棒性表述的直接依据） | 自监督 MLM；**Genecorpus-30M** = 29.9M cells / 561 datasets（仅 droplet 平台，**排除恶性/永生化细胞**——高突变负担） |
| **scPRINT-2** | multi-task bidirectional transformer + XPressor（10 latent cell tokens）+ 层级本体分类器（cell type/disease/assay 层级 + sex/ethnicity/organism/sequencer；**无 tissue 头**——v1 Methods 明确 omitted tissue/age）；small 版 **20M active params**。基因 token 用 **ESM3** 嵌入。v1 hidden dim 主文未报（Table S2 在补充材料，**1024-dim 待核实**）；我们的 424-dim = 256 emb + ontology heads 拼接（pipeline 记录） | 三大任务：denoising（ZIP 下采样重建，v1 r=0.6 / v2 r≈0.7）+ 本体层级分类 + bottleneck learning；聚合损失 Eq.20；imputation/counterfactual 为架构使然的推理时能力 | 自监督 + 监督多任务；v2：**350M cells / 16 organisms**（CxG LTS 2024-07-01 + Tahoe-100M + scBaseCount）；v1：**54,084,961 cells / 548 datasets**（CxG 1200→过滤：仅 primary cells、人+鼠、去 spatial，原文 Methods p.16 逐字确认） |
| **scGPT** | 生成式 transformer（**12 层、8 头、512-dim**；层块为 encoder 型）。注意力为**专门设计的生成掩码**：unknown 基因仅 attend 到 known 基因 + 自身 + `<cls>`，预测顺序由**预测置信度**迭代决定（原文明确：基因无序，不能直接套 GPT 因果掩码，仅 "reflects the causal masking design"；"随机排列"一说原文不存在） | 掩码基因表达预测（GEP，MSE）+ **GEPC**（以细胞表示预测基因表达）；ECS/DAR/DSBN 仅为微调目标 | 自监督；**33M 仅 non-disease** 人体细胞（51 organs / 441 studies / CELLxGENE Census 2023-05-15）。注：含 COVID/influenza 的 blood model 是**另一个 checkpoint**，不在 whole-human 语料内 |
| **UCE** | **33 层 transformer，674.7M 参数**（SI Table 2：d_emb 1,280、d_hid 5,120、20 头、dropout 0.05）；基因 token = **ESM2-15B** 蛋白嵌入；每细胞采样 1,024 基因（512 表达 + 512 未表达）；输出 = 末层 CLS + MLP 的 **1,280-dim** | **掩码基因建模**（mask 20% 表达基因，重建二值化表达，BCE 损失；20% 经消融最优）。**无显式跨物种对齐损失**——物种无关性完全来自 ESM2 蛋白 token 化（genome agnostic），对齐是 emergent 性质 | 全自监督（无细胞类型/数据集标签）；IMA：**36.2M cells / 313 datasets**（33.9M/285 来自 CxG Census v2023-07-10 + 2.3M/28 跨 8 物种；SI Fig.8b：Human 30.9M、Mouse 4.1M）；训练 43.5 天 × 24×A100 |
| **TranscriptFormer** | **自回归生成式 transformer**（"cell sentences" = 基因-表达对；表达量感知 multi-head attention；ESM-2 物种无关基因 token） | 基因身份与表达计数的**自回归联合建模**（next-token） | 跨物种生成式预训练；112M cells / 12 物种 / 1.53B 年进化跨度；本文用 tf_sapiens |
| **STACK** | **tabular attention**（intra-cellular MHA × inter-cellular MHA 交替；细胞=行、100 个 gene-module tokens=列；Large 9L/217M 参数）——每个细胞的表示由其所在 256 细胞 cell set 共同决定，**无 context-independent 固定向量**（Fig 1D + p.5 shuffle 实验证实） | 预训练：rectangular-mask 基因分布重建（NB 负对数似然 + Sliced-Wasserstein 正则）；后训练：**自蒸馏条件反事实重建**（frozen teacher + EMA student，mask-diffusion 推理） | **in-context learning**；148.8M 训练细胞（scBaseCount QC 后 189M 的 80%，19,978 SRX；摘要口径 ~149M） |

## 三个维度对照结果

### 1) Architecture：评价范式不同，骨架家族不预测成绩

- **逐细胞固定编码器**（Geneformer/scPRINT-2/scGPT/UCE/STATE-SE 的 transformer 系）：
  embedding 冻结后接任意下游头——unified-input 协议天然适配
- **自回归生成式**（TranscriptFormer）：为生成/提示设计，embedding 是副产物，
  判别任务 mid-pack——生成保真 ≠ 判别几何
- **跨细胞 tabular attention**（STACK）：表示依赖上下文集合，没有可冻结的固定
  表示——从主集移除的程序性原因（论文 Fig 1D 与 shuffle 实验直接支持）

家族内部证据：bidirectional transformer encoder 家族内，Geneformer 0.786 ↔ UCE
0.661 的跨度（0.125）是该家族与唯一非 encoder（scGPT 0.720）差距（0.066）的两倍
——**家族内部差 > 家族间差，骨架不预测成绩**。

### 2) Training Objective：与结果的关联（"consistent with" 纪律内）

| 目标 | 模型 | 在本 benchmark 的表现签名 |
|---|---|---|
| 连续程序几何（rank-value MLM / 双轴表达重建） | Geneformer、STATE-SE | within-COVID 并列第一（0.786/0.778）；Geneformer 兼队列最稳（LOBO σ=0.062） |
| 掩码重建（BCE/NB 分布）+ 本体多任务 | UCE、scPRINT-2 | UCE 细胞级判别垫底（0.661）但表示冗余（70% 置零不掉点）；scPRINT-2 弹性好但 seed 敏感（0.740） |
| 生成式重建（因果掩码 / 自回归） | scGPT、TranscriptFormer | 重建保真优先于判别几何：scGPT 0.720、TF mid-pack；scGPT 另有 non-disease-only 语料的结构性疾病短板 |
| 上下文范例预测（ICL） | STACK | zero-shot 崩溃（0.346）→ 随目标疾病数据陡升（0.832，legacy 测试集）——ICL 签名 |

注意：UCE 与 scGPT 的预训练目标同属"掩码表达重建"，结果一差一中——目标函数
的区分度同样有限，每格 n=1，以上仅作 "consistent with" 陈述。

### 3) Learning Paradigm：STACK 与 TranscriptFormer 补上的两个对照

- **STACK（ICL）**：zero-shot 0.346 → 有目标疾病上下文后陡升（旧曲线 0.832 未饱和）；
  它回答"给范例多快能适应"，不是"零样本表示好不好"。注意 2026-09-06 核实：
  其 adaptation 曲线数字来自 **legacy 4-patient 测试集**，与主 5 FM 的扩展池
  34-patient 检验集不可并列（见 `hpcc_data_verification_2026-09-06.md`）。
- **TranscriptFormer（跨物种生成式）**：覆盖最广（12 物种），单一人体疾病判别
  mid-pack——广度换深度。
- within-COVID 排序 ≈ 预训练分布与"真实人体患病免疫生物学"的距离；
  **跨病转移七个模型全部失败**—— paradigm 层面最一致的阴性。

## 排除逻辑与 OOD 防线的一致性

统一协议对所有模型一视同仁地测量——包括它们没被优化的任务，因为那正是
benchmark 的意义。模型只在两种情况下不进主集：协议在程序上无法可比地测量它
（STACK，表示不可冻结），或它对本文问题没有边际信息（TranscriptFormer，保留为
adaptation 曲线 secondary）。同一逻辑也回击"platelet 太远 OOD 所以结论不重要"的
审稿人质疑。

## 核实来源与待下载

**已核实（本文档数字均出自这些原文，含页码/引文的核实报告存档于本目录）**：
STATE (bioRxiv 661135v2)、scPRINT-2 (bioRxiv 693702v3)、UCE (Nature 2026 正文+SI)、
STACK (bioRxiv 698608v2)、Geneformer V1 (Nature 2023 正文)、scGPT (Nat Methods 2024)、
scPRINT v1 (Nat Commun 2025)、scBaseCount (Youngblut 2025, bioRxiv 640494v2)。

**P1-E platelet/megakaryocyte 逐模型状态（2026-09-06 全文核实后）：**

| 模型（语料） | platelet/MK 状态 | 依据 |
|---|---|---|
| STATE（scBaseCount 71.7M + CxG 59.2M + Tahoe 36.2M） | scBaseCount 部分 **ABSENT（注释层 + QC 双重层面）**；Tahoe 为癌细胞系（不可能含）；CxG 部分 ≤0.1% 上限 | scBaseCount Fig.3B 无 platelet/MK 类别（注释空间=Tabula Sapiens，本身无 platelet）；QC ≥300 genes/≥500 UMIs 必然滤除；Census 0.120–0.125% |
| STACK（scBaseCount 148.8M） | 同上（scBaseCount ABSENT） | 同上 |
| UCE（IMA：33.9M CxG + 2.3M 多物种） | **ABSENT**（SI Table 19 决定性：2,969,114 细胞代表样本的 51 个粗粒度类型完整清单中无 platelet/megakaryocyte；SI Table 16 训练集清单无 platelet 命名数据集） | Nature SI Excel（MOESM7/4），2026-09-06 核实 |
| Geneformer（Genecorpus-30M/104M，排除恶性） | **文档层已穷尽：任何出版物（正文/全部 15 张补充表/HF 页）均无细胞类型级构成**。器官级血液/免疫丰富：190/561（immune 111、bone_marrow 49、cord_blood 12 等）——巨核细胞（骨髓驻留）很可能在训，且其补充材料的网络分析**实际使用了小鼠 megakaryocyte 扰动数据集**（gata1 MK lof 等，属评估数据非预训练语料）；循环无核血小板不可判定 | Supp. Table 1（器官级）+ Supp. Tables 3–15（hits 均为 GO/HP 通路名与小鼠扰动集标签，非语料构成），2026-09-06 |
| scGPT（33M non-disease CxG 2023-05-15） | **文档层 INDETERMINATE**；Census 上限 0.125% | 主文零命中（megakaryocyte 仅出现在下游评估数据） |
| scPRINT v1（54M/548，CxG primary） | **ABSENT（正文+补充材料全文零命中；~400 类标签空间亦无 platelet 类）** | v1 主文 + scPRINT_supplementary.pdf 全文检索 0 命中，2026-09-06 |
| scPRINT-2（350M：CxG + Tahoe + scBaseCount） | Tahoe/scBaseCount 部分同 STATE；CxG 部分 ≤0.1% | 同上 |

→ OOD 论证新口径（比 "essentially absent" 更强且全部有据）：**没有任何模型报告其
语料含血小板；对 scBaseCount 系语料（STATE 的 ~43%、STACK 的 100%），血小板在
注释层结构上不存在且会被 QC 阈值（≥300 genes/≥500 UMIs）滤除；对 Census 系语料
（UCE/scGPT/scPRINT 系），血小板 ≤0.1% 且各模型子集过滤后的实际含量从未测量；
UCE 的 51 类粗粒度构成（SI Table 19）与 scPRINT v1 的全文+补充检索均证实无
platelet/megakaryocyte。唯一的灰色地带是 Geneformer：其语料富含血液/免疫
（190/561 数据集，含 49 个骨髓集），巨核细胞（有核前体）很可能在训，
但循环无核血小板与巨核细胞转录本并不等价（注释错位论点同样适用于此）。**

**仍待获取（均为可选增强，不阻塞写作）：**
1. scPRINT-2 `small-v2.ckpt` 内部 hidden state 维度（"1,024-dim" 一说需以
   checkpoint 配置为准；v1 补充 Table S2 文本层损坏无法确读）

**已闭环**：UCE SI Excel（Table 19 决定性 ABSENT）、Geneformer V2 规格
（HF 官方页 + checkpoint 配置双重核实）、Genecorpus-30M Supp. Table 1
（器官级构成：血液/免疫 190/561）、scPRINT v1 补充材料（platelet 零命中）、
scBaseCount（结构性 ABSENT）、scGPT（33M non-disease 确认）。

**核实发现的其他修正**：scPRINT-2 论文自曝 Open Problems 多数数据集在其
预训练语料中（其 zero-shot SOTA 声明需带此限定引用）；UCE 论文的模型对比
仅含 Geneformer/scGPT/tGPT/scVI/scArches/SATURN/SAMap/PCA（无
scFoundation/TranscriptFormer/STATE）；STATE 论文的 SE 对比仅含
scFoundation/scGPT/UCE/TF（无 Geneformer）；scGPT 的 "random permutation +
causal mask" 流行描述**与原文不符**（实为置信度排序的专用生成掩码）；
Geneformer V1 为 6 层/256-dim（勿与 V2-316M 混淆）；scBaseCount 实为
**502M cells**（非早前记录的 189M），blood+bone marrow 占比最大、癌症 30.7%。

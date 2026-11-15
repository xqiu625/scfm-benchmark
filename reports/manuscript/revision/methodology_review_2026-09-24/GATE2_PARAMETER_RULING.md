# Gate 2：十一参数一次性裁决

**日期：2026-09-24。角色：独立审阅者。裁决：11/11 方法学参数获批；尚未授权 pilot。**

本裁决是用户要求的单次 AI 方法学审阅；既往参与文件核验的情况见 [REVIEW.md](REVIEW.md)，不声称外部人类评审、盲审或预注册。适用目的为 **Gate 2 prior 敏感性的小型 scoring-only pilot**，不是最优检测器选择、真实 prevalence 校准或细胞删除许可。

## 资料边界与优先关系

开始时工作区只有制作资料包的任务说明；工作期间新增了 [参数资料包](../gate2_parameter_adjudication_2026-09-24/PACKET.md) 与 [资格草稿](../gate2_cohort_eligibility_2026-09-24/ELIGIBILITY_DRAFT.md)，本次已全部阅读并核对。裁决同时依据仓库 runner、现有方法学裁决与 **Scanpy 1.11.3 官方 tag 的三个源码文件**。源码地址、行号、资料包 SHA-256 见 [SOURCE_VERIFICATION.json](SOURCE_VERIFICATION.json)。下表补全资料包未提供的批准值、范围、理由和失效条件；不要求再制作一轮包才能作决定。

资料包实质修订：① Stage 31 待决的 PCA 不能作为 Gate 2 参数自动待决的依据；② prior-SD 进入误差/Z，不拓宽用于点分数计算的 prior 网格；③ Gate 2 面向完整 capture，不是 platelet-only counts。资格草稿中“上游移除了 doublets”证明残余模拟分布需要额外限制，不能单凭此断言所有 singlet 父状态已缺失；C3 的输入完整性/父状态支持与残余概率校准应分别记账。目前无 capture 完成 C1–C3 验收的结论不变。

本文仅替代 REVIEW.md 中“另外十一参数尚未裁决”的状态；其八点 prior 网格、cohort 分类、C1–C5、unknown/NA、无删除、无 outcome 调参等规定继续有效。**Gate 2 的 PC=30 不适用于 Stage 31**；Stage 31 新总体仍须先完成来源和 identity 闭合，再授权 dimension-only audit，之后另行选择固定维度。

## 统一适用范围

- Scanpy **1.11.3**，`gate2-preprocessing/1`：每个已验收的物理 capture 独立处理，输入为已证明来源的未归一化 RNA counts；完整保留所声明 cell-called/basic-QC 阶段的人群，不能先挑目标谱系或平衡抽样。基因/细胞最小过滤、HVG 顺序沿用核对过的 runner。
- Qiu、Chongqing：条件满足后的 standard scoring；Cillo、Liu、COMBAT：仅 residual-input diagnostic，必须保留 `residual_null_unvalidated`；Ren：单列 `residual_scoring_only`。这些分类并不使任何 capture 自动合格。Stephenson、Wilk 当前不适用。
- 本文的“批准范围”是允许执行的值或确定性规则，**不是可自行搜索的调参区间**。除 k 随合格 capture 的输入大小解析外，每个参数只批准一个值。改变任何值、算法版本、输入阶段或用途须重开相应裁决。
- 以下默认均指公开入口 `scanpy.pp.scrublet` 的真实签名，不混用内部函数或独立 Scrublet 包默认。默认相同只是选择结果；理由是固定可解释的参照实验，而非“库默认必然正确”。

## 逐项裁决（全部 APPROVED_SCOPED）

| 参数 | Scanpy 1.11.3 默认 | 本次批准值 | 批准范围 | 作用及理由 | 失效条件 |
|---|---|---|---|---|---|
| `sim_doublet_ratio` | `2.0` | `2.0` | 仅 `{2.0}`；`N_sim=2*N_obs` | 在原始 RNA counts 的 HVG 子集上模拟；固定模拟规模，避免把模拟量也变成 prior 网格的变化来源。2 倍为明确的参照设计，不声称已证明蒙特卡洛误差最小。 | 实际模拟数量不符、跨 prior 重抽父细胞、输入/父状态不满足 C1–C3、内存不足；不得自动降低比值。 |
| `synthetic_doublet_umi_subsampling` | `1.0` | `1.0` | 仅 `{1.0}` | 相加父细胞 UMI，不额外 thinning；当前没有可识别的 capture 效率依据支持另加采样率。 | 输入不是未归一化 RNA counts、进行了未声明的 thinning，或有新实验依据要求另一分子采样模型。 |
| `log_transform` | `False` | `False` | 仅 `{False}` | 固定进入 PCA 前的变换；不增加可能改变双细胞几何结构的另一个因素。HVG 选择所需临时 log1p 层仍保留，此值不代表整个流程没有 log1p。 | observed/simulated 处理不一致、上游已 log 后再次当 counts、HVG/归一化顺序改变，或实际运行路径忽略本值。 |
| `n_prin_comps` | `30` | `30` | 仅整数 `{30}`；要求 `30 < min(N_eff,G_HVG)`，PCA 输入有限且具有支持 30 个非退化分量的秩 | 为 Gate 2 完整 capture 的 prior 控制实验固定一个公开参照维度；不是用旧 Stage 31 的 301/161 覆盖率选值，也不宣称 30 是生物学最优。 | 维度不满足严格不等式、秩退化/非有限、PCA 不收敛、细胞总体或预处理改变；该 capture 记 unknown，不能降到 29 或按样本改值。Stage 31 不在批准范围。 |
| `n_neighbors` | `None` | `None`，运行前按规则解析整数 | 唯一规则 `k=round(0.5*sqrt(N_eff))`；`k_adj=round(k*(1+N_sim/N_eff))`；要求 `1<=k` 且 `k_adj<N_eff+N_sim` | 使用实际过滤后 observed 细胞数，保留原生规模适应规则；每个 capture/seed 的 k 在八个 prior 上固定。允许不同 capture 的 k 不同，但这不赋予跨 capture prevalence 可比性。 | 用 donor 数、原始 droplets 数、模拟加 observed 总数或抽样后目标谱系数替代 N_eff；库修改邻居数/舍入规则；跨 prior 改 k；解析/回显不符；不可静默 clip。 |
| `stdev_doublet_rate` | `0.02` | `0.0` | 仅 `{0.0}` | 网格每点是固定假设，不是一个已估计 prevalence 的不确定性分布；故将 prior-SD 项置零。源码中该参数进入 score-error/Z 计算，不进入点分数或 score-threshold calls；邻居采样误差项仍存在。 | 将输出误差/Z 当真实概率的置信度、引入实测 prior 不确定性或以误差/Z 做筛选；此时需重新定义误差模型。0 不表示真实 doublet rate 已知。 |
| `knn_dist_metric` | `'euclidean'` | `'euclidean'` | 仅 `{'euclidean'}` | 在同一居中、按基因方差标准化的 PCA 空间内固定距离，保持明确的参照几何结构。 | 实际使用 cosine/其他距离，换 embedding/尺度定义，或出现非有限距离；不得因结果不理想换 metric。 |
| `normalize_variance` | `True` | `True` | 仅 `{True}` | 使用 observed 基因方差并同样应用到 simulated，避免高方差尺度单独主导距离；与 mean_center=True 构成同一标准化。 | 保留的基因有零/非有限方差、observed 与 simulated 独立拟合尺度，或临时关闭以挽救数值错误。 |
| `mean_center` | `True` | `True` | 仅 `{True}` | 使用 observed 均值处理两组，拟合 observed PCA 后投影模拟组；锁定 PCA 路径，不将不同分解混入 prior 比较。 | 改为未居中的 TruncatedSVD、联合或仅模拟数据拟合，或内存不足导致实现切换。1.11.3 会转 dense，必须先满足资源预算。 |
| `use_approx_neighbors` | `None` | `False` | 仅 `{False}` | 指定 sklearn 精确邻居分支，避免 `None` 自动选择后端及近似搜索的额外变化；仍须固定顺序、依赖与检查 ties/图，不声称跨版本比特一致。 | 回落到近似/自动后端、依赖或 tie/self-neighbor 语义改变、预算不可承受、同一 capture 跨 prior 的图不同；不能静默改 True。 |
| `get_doublet_neighbor_parents` | `False` | `False` | 仅 `{False}` | 本 pilot 不解释逐细胞父状态，此开关只生成额外邻居父细胞诊断，不改变点分数。模拟父对本身仍需保存/校验，以检查跨 prior 复用。 | 要求父状态推断或诊断该额外映射；需另开明确保存/审阅方案。False 不能作为免除模拟父对 provenance 的理由。 |

`N_eff` 是经过该固定预处理后实际进入 PCA 的 observed 细胞数；`G_HVG` 是实际保留的基因数。规模条件只是技术可行性，不能替代 C3。所有单值范围均适用于上节允许的 cohort/use class，而不是所有来源资料。

机器可读逐项表：[GATE2_PARAMETER_RULING.csv](GATE2_PARAMETER_RULING.csv)。已批准值的代码记录：`scripts/32_rebuild_upfront_qc/code/gate2_parameter_ruling.py`，版本 `gate2-parameter-ruling/1`。

## 实际调用路径与源码核验

| 参数组 | runner 实际落点 | 提供外部 adata_sim 时的语义 |
|---|---|---|
| 两个 simulation 参数 | `prepare_capture → sc.pp.scrublet_simulate_doublets(..., layer='raw')` | 若仅交给 scoring wrapper 会被跳过；必须在模拟时实际应用，且每 capture/seed 只模拟一次。 |
| `log_transform` | `prepare_capture` 中 observed/simulated 同步条件分支 | 若仅交给外部 adata_sim 的 scoring wrapper 会被跳过；本次 False 冻结该分支不执行。 |
| 其余八个 | `scrublet_call_parameters → sc.pp.scrublet → _scrublet_call_doublets → core/pipeline` | 仍生效；`expected_doublet_rate` 与 `stdev_doublet_rate` 也生效，不能依文档的“where adata_sim not supplied”字样误判忽略。 |

源码核验的三个区别：公开入口 PC 默认为 30，内部 `pipeline.pca` 默认 50 **不适用**，因为 wrapper 显式传值；`use_approx_neighbors=True` 实际分支为 PyNNDescent，并非 docstring 中的 Annoy；`stdev_doublet_rate=0.0` 不改变 score-threshold 判定但改变误差/Z，且当前 wrapper 不返回这些误差数组，不能声称已审阅它们。

发现并修改了一个转发缺陷：原 runner 请求 `n_neighbors=None`，Scanpy 回显解析后的整数，原样比较会把正确运行标成 `parameter_echo_mismatch`。现在在调用前按已批准规则使用 post-filter N_eff 解析，传入与记录同一整数，再与库回显比较。没有弱化不一致检查。

## 配套固定项与 config 审计结论

继续使用已批准 prior 网格 `0.005,0.01,0.02,0.04,0.08,0.12,0.16,0.20`，标签 `assumed_algorithmic_stress_setting`；不把高 prior 当真实率的上下界。种子沿用 `DESIGN_SEED=20260923` 和现有 capture_seed 推导，所有 prior 相同。未批准第二组 seed；seed stability 必须保持未测标记，不能宣传已验证稳定。

阈值规则明确固定为 Scanpy 1.11.3 自动 `threshold=None`；失败保持 unknown/NA，不能手调至期望阳性比例。score 与 call 状态分别记录，不因阈值失败抹掉有效分数。只报告敏感性，不挑最佳网格点，不以模拟/observed AUC 当真值准确率。

**本轮 config 结论：NOT READY；没有可授权的运行配置。** 十一参数裁决已结束，剩余验收是具体实现/输入：

1. 至少一个允许 cohort 的具体 capture 提交自己的 C1–C3 来源、物理单元、raw barcode 和完整输入阶段证据；仅写 limitation 不能把未满足条件变成满足。
2. 配置绑定 capture/input hash、批准的十一值/规则、Scanpy 与数值依赖版本、阈值、种子、输出位置、资源上限和严格有效维度检查；参数类型/值必须与裁决一致，测试注入接口不是生产许可。
3. 同一 capture/seed 的 HVG、模拟父对、PCA 和邻居图在 prior 间固定并有可审计记录。**现 runner 只保存模拟矩阵 checksum，且每个 prior 重做 PCA/kNN；相同 seed 不足以证明图相同。** 实现应复用同一表示/图，或直接记录并核验重算结果一致；仅 rank/AUC 不变或模拟 hash 一致均不能替代图证据。无临床 outcome 的技术 fixture 可用于此核验，不得拿 cohort 结果调参。

在这些事项完成前，`METHODOLOGY_PARAMETER_RULING` 记录 11/11 已批准，`PRODUCTION_PARAMETER_RULING={}` 仍表示**尚未启用的部署槽位**，不再表示参数未裁决。draft 会显示 0 个未裁决参数，并维持全部 `runtime_eligible=false`。本次没有生成运行配置或提交任务。

## 两条启动规则（采用用户最新授权）

**Gate 2：** `11/11 方法学参数获批 AND 至少一个 capture 的 C1–C3 满足 AND 该配置审计通过 → 独立审阅者记录限定范围的小型 pilot 授权 → 才运行`。许可须指定合格 capture、输入/config hash、完整八点网格、资源限制、独立输出路径和停止条件；首个 pilot 限一个完整 capture，不能为缩小任务先做细胞抽样。当前第一项已满足，后两项和最终许可未满足，因此**不能跑**。用户已把条件满足后的 pilot 许可交由审阅者裁定，不另设重复用户批准步骤；全队列扩展及删除不在这项许可内。

**Stage 31：** 逐来源 raw-count provenance、物理 capture、raw barcode 与 canonical 细胞的一对一对应闭合，并确认完整 pre-subsample/cell-called 输入后，审阅者才可记录**新 dimension-only audit** 许可；新目录写测量/status，不运行 doublet detection，不读取 outcome、不预填 PCA 值。本次未验收这三项，**不能跑**。旧 301/161 结果、canonical token 唯一性、Gate 2 的 PC=30 和 DecontX 抽样背景均不满足该启动条件。

## 本次验证范围

十一行 CSV、方法学参数映射和官方签名已交叉核对；`git diff --check` 通过。与本次变更直接相关的 **35 项测试通过**，覆盖参数裁决完整但生产仍拦截、缺任一参数的测试输入被拒、全部 cohort 仍不可运行，以及不同 N_eff/八个 prior 的 k 解析和错误回显不被放过。日志见 [parameter_ruling_tests.txt](parameter_ruling_tests.txt)。

本地没有 pandas/anndata/Scanpy。最初常规测试收集因缺 pandas 中止；绕过与本次断言无关的共享 conftest 后，较大两文件测试尝试为 53 passed、52 failed、10 deselected，失败路径缺 anndata。这不记作全套通过。最终 35 项使用 `--noconftest` 运行不依赖这些库的定向断言；本次未跑真实 Scanpy、cohort scoring、HPCC 任务或新的 dimension audit。源码核验与本地控制流测试不能替代待完成的真实配置审计。

## 主要依据

- [Scanpy 1.11.3 wrapper / simulation](https://github.com/scverse/scanpy/blob/1.11.3/src/scanpy/preprocessing/_scrublet/__init__.py)：公开签名、外部模拟对象分支、转发与返回值。
- [Scanpy 1.11.3 core](https://github.com/scverse/scanpy/blob/1.11.3/src/scanpy/preprocessing/_scrublet/core.py)：k 规则、后端、score/error 公式与阈值。
- [Scanpy 1.11.3 pipeline](https://github.com/scverse/scanpy/blob/1.11.3/src/scanpy/preprocessing/_scrublet/pipeline.py)：observed 拟合的尺度/PCA、dense 转换及 ARPACK。
- [Scrublet 原论文](https://pmc.ncbi.nlm.nih.gov/articles/PMC6625319/)：模拟亲本支持和可检测双细胞的限制，沿用前次 REVIEW.md 已核对的依据。本文对固定实验设置的选择是审阅者判断，未宣称原论文验证了这些 cohort 的全部参数。

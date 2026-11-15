# Gate 2 单 capture pilot 候选：最终审查

**独立审阅者结论：三项通过，批准本候选 commit/push。**

范围仅为用户指定的最终字节/证据、permit、pilot fail-closed。本记录来自自动化审阅；本次也实施了一处必要修复，不声称外部人类盲审。本文没有创建生产 permit，也没有授权或提交评分任务。下一步保持既定的 scoped permit → 干净部署 → dry-run → 唯一 scoring pilot，不要求额外架构审查。

## 1. 最终字节与证据：PASS

- 两处旧图捕获描述已核对为当前 consumed-graph helper `_get_indices_distances_from_sparse_matrix`，历史 `NearestNeighbors.kneighbors` 只保留为失败方案的说明。
- 审计文档的旧 `HEAD == implementation_commit` 表述已改为祖先 commit + 七个文件 blob 一致；旧 preflight commit 不再标为 scoring deploy commit。
- **529 个冻结源码、脚本、测试文件**：测试前、每组测试后、每次 mutation 恢复后、最终 SHA 全部逐项一致；回传后也与本地同一清单逐项核对。

| 最终验收运行 | 结果 | exit |
|---|---|---|
| focused | 73 passed | 0 |
| 真实 Scanpy integration | 12 passed，714 deselected | 0 |
| 完整 tests/ | 726 passed | 0 |
| 既定 M4–M12 | 9/9 实际注入、被指定断言杀死、恢复后字节一致 | 每个 mutation 的预期失败 exit 1 |

环境：Python 3.10.18、Scanpy 1.11.3、anndata 0.11.4、NumPy 2.2.6、SciPy 1.15.3、sklearn 1.7.1、pytest 9.1.1。

首次完整测试为 724 passed / 2 failed，原因是测试调用的 PYTHONPATH 未包含仓库 preflight 目录，两个 import 报 ModuleNotFoundError。**未改任何冻结文件**；显式补齐仓库内模块路径后，完整测试 726 passed。首次失败完整保留，不作为通过证据。

精确命令、环境、前后 SHA 和退出码见 [测试 transcript](SCORING_TEST_TRANSCRIPT.txt)、[mutation transcript](SCORING_MUTATION_TRANSCRIPT.txt)；原始日志、每阶段 SHA、两次调用脚本及机器记录见 [原始证据目录](final_evidence_20260925T204838Z/) 和 [SUMMARY.json](final_evidence_20260925T204838Z/SUMMARY.json)。SUMMARY 中 `full_suite` 是保留的失败尝试，最终完整测试验收项是 `full_suite_retry`。

## 2. Permit：PASS

现有实现要求 implementation commit 存在并为 deployed HEAD 祖先；七个实现文件在两 commit 间 blob 相同，且工作树/index 与 HEAD 相同；permit 位于部署仓库内、tracked、已提交且无修改。capture、config/input/implementation hash、output root、scope、resource ceiling 必须匹配。

除了套件中的正反例，本次还使用完整的两 commit 临时仓库，调用**真实 permit repository verifier**，未传 `repo_verifier` 替身：commit A 是实现，commit B 提交完整 TEST_ONLY permit 并引用 A，验证通过；未跟踪、未提交修改、已暂存但未提交的 permit 都被拒绝。见 [核验记录](PERMIT_SATISFIABILITY_CHECK.txt)。这份临时 fixture 不属于部署仓库，也不是生产授权。

真实 wrapper 的行为测试证明：缺 permit 的普通模式 exit 3；dry-run exit 0 并在评分前返回。已移除的环境 bypass 即使值为 `0` 也被拒绝。部署时仍需用实际 permit commit 的干净 worktree 验证；不能把本地未提交草稿当部署状态。

批准的后续 permit scope 不变：**one job, one capture, eight approved priors, scoring only; no outcome access, no cell deletion, no retry, no second job**。资源 ceiling：short、2:00:00、48G、4 CPU、max_jobs=1。只适用于 Qiu `GSM5102900_HC1` 的完整 7,666-cell deposit。

## 3. Pilot fail-closed：PASS

- 图捕获次数、k_adj、shape、非空 checksum、跨 prior topology 都进入问题列表；有问题就写 stopped / graph_unacceptable 并非零退出。真实 Scanpy integration 验证实际捕获点。距离只记录，不作为一致性 gate。
- threshold failure 保留身份已核对的分数；calls 不进入 boolean sidecar，manifest 用 unknown/null，不能当 False/0。自动阈值失败不需要伪装成评分成功。
- sidecar 临时写入、长度/内容核对与 promotion 在 manifest complete 之前；未完成的 sidecar 集合不能依据部分文件自称完成。不扩展为目录级事务要求。
- 入口只读三个 deposit 文件，obs 只保留细胞身份；无 outcome 读取、无删除动作、无 production ruling 激活或写入。全局 `PRODUCTION_PARAMETER_RULING` 保持 `{}`，permit 通过后仅传 scoped ruling。

**本次唯一新增的实质修复：身份不一致的分数不得发布。** 原 runner 可能保留已发现细胞顺序不匹配的分数向量，入口也可能随后写 complete。本次限制 score sidecar 仅接收 identity-verified 向量，并让入口遇到 `identity_verified=false` 即 stopped / drift_detected。入口行为测试验证不会发布分数 sidecar 或 complete；真实 Scanpy reordered-result 测试验证不输出错配分数。阈值失败且身份正确的分数保留逻辑不变。

## 批准对象及下一步

| 绑定 | 值 |
|---|---|
| config SHA-256 | `0f02df7c7e2a81da1efa2f371d6890878ad8f747146986ae5bc5db7a94f1dd41` |
| implementation SHA-256（七文件） | `30831273d4d3030ef5f9b123ecefbb44af2ae37aa9b7a95aab18ffeb8de80d84` |
| 派生输出目录 basename | `GSM5102900_HC1__cfg0f02df7c7e2a__impl30831273d4d3` |
| 冻结清单 | `final_evidence_20260925T204838Z/FROZEN_BYTES.json` |

本文件随实现提交进入的 commit 即下一步 permit 应引用的 reviewed implementation commit；不要引用旧的 `9aab7be` preflight commit。随后将 scoped permit 提交为后续 commit，部署至全新干净 worktree，保存 dry-run 的 commit/hash/permit/output 验证。只核对该 permit 与 dry-run，即可按已授予的范围授权唯一 pilot。本文没有增加新的方法学或架构门槛。

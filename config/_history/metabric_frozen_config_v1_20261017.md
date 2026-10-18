# METABRIC 配置 v1（2026-10-17 18:40，**已停用**）

> ⚠️ **本文件仅供追溯，不再生效。** 当前生效的配置是 `config/metabric_frozen_config.md`（v2）。
> v1 有两处事实错误：把化生性癌当作"非浸润性癌"排除；未做基因符号规范化。
> 修正见 v2 的修订说明。

冻结时间：2026-10-17（在任何 METABRIC 模型被拟合之前）。
依据：第四轮评审第 51 行——"运行前冻结 OS 定义、入排标准、5 次划分种子、候选清单、规则及缺失处理"。

---

## 1. 数据

| 项 | 值 |
|---|---|
| 来源 | cBioPortal 公共 API，study `brca_metabric` |
| 表达谱 | `brca_metabric_mrna` — "mRNA expression (Illumina HT-12 v3 microarray)" |
| 下载日期 | 2026-10-17 |
| 表达矩阵 | `data/cohorts/METABRIC_expr.tsv`（209 基因 × 1,975 样本） |
| 生存表 | `data/cohorts/METABRIC_surv.tsv` |
| 溯源 | `data/cohorts/METABRIC_provenance.json`（含两个文件的 md5） |

## 2. 终点

| 项 | 值 |
|---|---|
| 时间 | `OS_MONTHS` — "Overall survival in months **since initial diagnosis**" |
| 状态 | `OS_STATUS`；**`1:DECEASED` = 事件**，`0:LIVING` = 删失 |
| 其他终点 | `RFS_MONTHS`/`RFS_STATUS` 存在但**本实验不使用** |

## 3. 入排标准（永久固定）

1. `SAMPLE_TYPE == "Primary"`（防止转移/复发灶混入）
2. `OS_MONTHS` 非缺失且 **> 0**
3. `OS_STATUS ∈ {1:DECEASED, 0:LIVING}`（其他取值按缺失处理）
4. **排除 4 例非浸润性癌**：`ONCOTREE_CODE` 为 `PBS`（乳腺血管肉瘤）或 `MBC`（化生性癌）

> 决策依据：与 TCGA-BRCA（Breast Invasive Carcinoma）的组织学定义保持一致。
> 该决策在可用性检查阶段提出、**在任何模型拟合之前**确定，不是事后选择。

**结果**：1,980 个表达样本中，排除 1 例（OS_MONTHS ≤ 0）与 4 例非浸润性癌，
**最终 1,975 例，1,141 事件（57.8%）**。

## 4. 划分

- 每个（签名 × 种子）独立做一次**按事件分层**的 50% / 25% / 25% 划分；
- 种子与主分析**完全相同**：`20261014, 20261016, 20261017, 20261018, 20261019`；
- 划分算法与 `19_main_corrected.R` **逐字相同**；
- 5 次划分共享患者，**不是** 5 个独立队列。

## 5. 候选池

与主分析**逐字一致**的 18 个配置（`scripts/_lib_fit.R`，镜像自 19）：

```
Lasso, Ridge, Enet0.9/0.7/0.5/0.3/0.1, Ridge_1se, UniCox+Lasso,
LassoTop100/50/25/10, Enet0.5Top25, RSF100, GBM100, GBM250, RSF250
```

- 全部 `cox.ties = "breslow"`，`nfolds = 5`，`set.seed(1)`（与主分析同）；
- 可用候选 = 在该任务上拟合成功的那些（与主分析同样的 90% 规则不适用——
  此处每个（签名×划分）只跑一次，因此"可用"= 拟合成功）。

## 6. 预处理

| 步骤 | 设定 |
|---|---|
| **log2** | **不做** —— 该谱已是 log2 强度（4.68–14.64，无负值）。主分析的 `log2(x+1)` **必须关掉** |
| 缺失过滤 | 训练集缺失率 > 20% 的基因剔除 |
| 缺失填补 | 训练集该基因的中位数 |
| 零方差剔除 | 训练集 sd 为 0 或非有限的基因剔除 |
| 标准化 | 训练集均值/标准差，原样应用到验证与测试 |
| 实测缺失 | 412,775 个取值中 **1 个缺失**（由上述规则处理并计数） |

## 7. 规则

`scripts/_lib_analysis.R::rules_one()`，15 条规则全部计算：

`none, argmax_val, argmax_nonreference, runner_up, linear_only, tree_only,
top2/3/5/10_ensemble, all_ensemble, majority_vote, worst_half, penalized, oracle_test`

- `B_VOTE = 12`、`B_HALF = 6`（与主分析相同）；
- 并列用 `make_prio(mods, seed)`（按候选名排序后生成，与输入列顺序无关）；
- **正文只报 5 条**：固定 Lasso、固定 Ridge、固定 Enet0.5、`argmax_val`、`top5_ensemble`；
  其余放补充（评审要求）。

## 8. 报告口径（预先规定）

- 报告 held-out 性能（`Ct`）、相对各固定基线的配对差（`Gv`/`Gt`/`Orel`）；
- 先按签名汇总 5 个划分，再对 6 个签名等权汇总；**保留每个签名的结果**；
- **不把 30 个任务当 30 个独立队列做显著性检验**；
- 以**透明的描述性结果**为主；
- **无论结果支持、削弱还是反转当前观察，都照实报告**；
- 不为新队列重新搜索 `k` 或阈值。

## 9. 不在本实验范围内

- 不把 TCGA 的冻结模型迁移过来做外部预测（那是另一类研究，本研究不做）；
- 不复现原商业检测的风险评分（只取基因列表重新训练）；
- 不扩展到 MSigDB 全池（只跑 6 个已发表签名）。

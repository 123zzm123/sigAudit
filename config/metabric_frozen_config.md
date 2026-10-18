# METABRIC 独立队列：**当前生效**的配置

版本：**v2**（2026-10-17 晚修订，第五轮评审后）。这是**唯一生效**的配置。
v1 已移入 `config/_history/metabric_frozen_config_v1_20261017.md`，仅供追溯，**不再生效**。

> **修订说明（2026-10-17，第五轮评审后）**
> v1 中有两处错误，属于**纠错**而非按结果调整设计：
> 1. **组织学筛选**：v1 把化生性癌（MBC）与血管肉瘤一起当作"非浸润性/非癌"排除。
>    化生性乳腺癌**本身就是浸润性乳腺癌**，这是事实错误。v2 只排除血管肉瘤。
> 2. **基因符号**：v1 让源列表原样取交集。旧符号→当前符号属于**标识规范化**
>    （HGNC 维护 approved / previous / alias），不是按结局挑基因。v2 先映射再取交集。
>
> **修订日期不得倒写成原始冻结日期。** 两处修改均在第五轮评审后作出。
>
> **修订导致的变化（如实记录）**：纳入 1,975 → **1,977**，事件 1,141 → **1,143**；
> 六个签名的可用基因变为 11/11、21/21、50/50、51/51、60/60、99/99；
> `argmax − ridge` 由 −0.0011 变为 **+0.0002**，`top5 − ridge` 由 +0.0031 变为 **+0.0039**。
> v1 的结果保留在 `results/tables/_superseded/`，正文只使用 v2。

---

## 1. 数据

| 项 | 值 |
|---|---|
| 来源 | cBioPortal 公共 API，study `brca_metabric` |
| 表达谱 | `brca_metabric_mrna` — "mRNA expression (Illumina HT-12 v3 microarray)" |
| 下载日期 | 2026-10-17（v2 用 HGNC 批准符号重新取数） |
| 表达矩阵 | `data/cohorts/METABRIC_expr.tsv`（**247 基因 × 1,977 样本**） |
| 生存表 | `data/cohorts/METABRIC_surv.tsv` |
| 溯源 | `data/cohorts/METABRIC_provenance.json`（含两文件 md5、修订标记） |
| 取数与重建脚本 | `scripts/35_rebuild_metabric.py` |

## 2. 终点

| 项 | 值 |
|---|---|
| 时间 | `OS_MONTHS` — "Overall survival in months **since initial diagnosis**" |
| 状态 | `OS_STATUS`；**`1:DECEASED` = 事件**，`0:LIVING` = 删失 |
| 其他终点 | `RFS_MONTHS`/`RFS_STATUS` 存在但**本实验不使用** |

## 3. 入排标准

1. `SAMPLE_TYPE == "Primary"`
2. `OS_MONTHS` 非缺失且 **> 0**
3. `OS_STATUS ∈ {1:DECEASED, 0:LIVING}`（其他取值按缺失处理）
4. **只排除血管肉瘤**（`ONCOTREE_CODE == "PBS"`，属肉瘤，不是癌）

> 化生性癌（`MBC`）**予以保留**：它属于浸润性乳腺癌，与本研究"浸润性乳腺癌"的
> 预设范围一致。v1 排除它是事实错误。

**结果**：1,980 个表达样本中，排除 1 例（`OS_MONTHS ≤ 0`）与 2 例血管肉瘤，
**最终 1,977 例，1,143 事件（57.8%）**。

## 4. 基因标识

- 六个发表签名的源符号先用 **HGNC** 的 `prev_symbol` / `alias_symbol` 映射到当前批准符号
  （`scripts/34_map_signature_symbols.py`；38/247 个符号发生映射）。
- **只对在两个表达矩阵中都找不到的符号查询 HGNC**；因此准确表述是
  "对未匹配的旧符号做了 HGNC 映射"，而不是"247 个都已验证为当前批准符号"。
- 映射逐条记录在 `results/manifest/signature_symbol_map.tsv`
  （源符号、批准符号、映射依据）；最终使用列表在
  `results/manifest/public_signatures_genes_mapped.tsv`。
- **多命中不得按字母序取值**：若一个旧符号命中多个批准符号，应保留原始稳定 ID 判断，
  必要时标注为不可唯一映射。当前 38 条映射的唯一性由脚本报告。

## 5. 划分

- 每个（签名 × 种子）独立做一次**按事件分层**的 50% / 25% / 25% 划分；
- 种子与主分析**完全相同**：`20261014, 20261016, 20261017, 20261018, 20261019`；
- 划分算法与 `19_main_corrected.R` **逐字相同**；
- 5 次划分共享患者，**不是** 5 个独立队列。

## 6. 候选池

与主分析**逐字一致**的 18 个配置（`scripts/_lib_fit.R`，镜像自 19）：

```
Lasso, Ridge, Enet0.9/0.7/0.5/0.3/0.1, Ridge_1se, UniCox+Lasso,
LassoTop100/50/25/10, Enet0.5Top25, RSF100, GBM100, GBM250, RSF250
```

- 全部 `cox.ties = "breslow"`，`nfolds = 5`，`set.seed(1)`（与主分析同）；
- "可用候选" = 在该任务上拟合成功的那些（每个（签名×划分）只跑一次）。

## 7. 预处理

| 步骤 | 设定 |
|---|---|
| **log2** | **不做** —— 按 cBioPortal 发布值原样使用（4.68–14.64，无负值）；范围本身不是"已取对数"的证明 |
| 缺失过滤 | 训练集缺失率 > 20% 的基因剔除 |
| 缺失填补 | 训练集该基因的中位数 |
| 零方差剔除 | 训练集 sd 为 0 或非有限的基因剔除 |
| 标准化 | 训练集均值/标准差，原样应用到验证与测试 |
| 实测缺失 | 极少数取值（逐任务记录在 `31_metabric_task_info.tsv` 的 `n_imputed`） |

## 8. 规则

`scripts/_lib_analysis.R::rules_one()`，15 条规则全部计算：

`none, argmax_val, argmax_nonreference, runner_up, linear_only, tree_only,
top2/3/5/10_ensemble, all_ensemble, majority_vote, worst_half, penalized, oracle_test`

- `B_VOTE = 12`、`B_HALF = 6`（与主分析相同）；
- 并列用 `make_prio(mods, seed)`（按候选名排序后生成，与输入列顺序无关）；
- **正文只报 5 条**：固定 LASSO、固定 ridge、固定 elastic net、`argmax_val`、`top5_ensemble`。

## 9. 报告口径

- 报告 held-out 性能（`Ct`）、相对各固定基线的配对差（`Gv`/`Gt`/`Orel`）；
- 先按签名汇总 5 个划分，再对 6 个签名等权汇总；**保留每个签名的结果**；
- **不把 30 个任务当 30 个独立队列做显著性检验**；
- 以**透明的描述性结果**为主；不为新队列重新搜索 `k` 或阈值。

## 10. 不在本实验范围内

- 不把 TCGA 的冻结模型迁移过来做外部预测；
- 不复现原商业检测的风险评分（只取基因列表重新训练）；
- 不扩展到 MSigDB 全池（只跑 6 个已发表签名）。

## 11. 依赖顺序（首次运行）

```
34_map_signature_symbols.py   # 需要 data/cohorts/{METABRIC,BRCA}_expr.tsv
        ↓
35_rebuild_metabric.py        # 用映射后的批准符号重新取数并重建 METABRIC
        ↓
31_metabric_audit.R           # 30 个任务；保存结果表 + 逐患者预测
        ↓
32_metabric_baseline_sensitivity.R  → 33_figure_metabric.R
```

**当前交付已包含冻结的映射表与重建后的数据**，因此首次复现不必联网重跑 34。

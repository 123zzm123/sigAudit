# README_current —— 当前状态（自动生成，请以本文件为准）

生成时间：**2026-10-19 14:03:41（时区 CST）**
本文件由 `scripts/27_status_readme.R` 从实际产物生成；**不要手工编辑**。

## 1. 唯一正式管线（只有这一条）

```
数据获取     scripts/02_fetch_cohorts.py, scripts/04_entrez2symbol.py
基因集构建   scripts/01b_build_genesets.py        （重建并逐集合核对一致）
正式主实验   run_v8_parallel.sh  ->  scripts/19_main_corrected.R（逐患者预测落盘）
清单         scripts/23_manifest.R（注册网格校验）, 24_manifests.R
主分析       scripts/20_analysis_corrected.R      （P1/P2 规则核算）
共用函数库   scripts/_lib_analysis.R              （20 与 29 共用，杜绝口径漂移）
辅助分析     scripts/21_entry_predictors.R        （进入效应 + 留一癌种敏感性）
基线敏感性   scripts/29_baseline_sensitivity.R    （Ridge / Enet0.5，第三轮审查后新增）
非恒定诊断   scripts/47_nonconstant_diagnostic.R  （非恒定 LASSO 分层；本次修订新增）
去重敏感性   scripts/48_dedup_sensitivity.R       （验证预测去重；本次修订新增）
预测级审计   scripts/49_audit_prediction_api.R   （R/audit_predictions.R 一致性核对）
现实任务     scripts/26a_extract_public_signatures.R -> 26b_public_signature_tasks.R
模拟         scripts/25_simulation.R pilot|final
图表         scripts/13_figures.R                 （v3：比例/裁切/COADREAD/oracle 已修）
状态文档     scripts/27_status_readme.R, scripts/28_freeze_current.sh
一键复现     reproduce_current.sh tables|signatures|metabric|train|all|check
验收         scripts/_verify_19.R, _verify_response_numbers.R, _verify_final_numbers.R
```

**不要**再把 08/09/16/17/18 当作正式管线：它们要么使用**旧的评价口径**
（在评价集上重新拟合方向、构成泄漏），要么只覆盖旧的三队列/旧池子。
`17_remedies.R` 已退役（功能并入 20）；根目录的 `reproduce.sh` 也已作废
（它跑的是 v6 + 09/10/11/12/13 旧流水线，**无法生成当前论文结果**），
统一改用 `reproduce_current.sh`。

**源码哈希（含共用库与冻结计划）由 `scripts/28_freeze_current.sh` 固化到**
`results/manifest/frozen_v1_<日期>.tsv`，那才是可核验的构建记录。

## 2. 源码哈希（MD5 前 10 位）

| 脚本 | 哈希 |
|---|---|
| `scripts/19_main_corrected.R` | `5f22b39406` |
| `scripts/20_analysis_corrected.R` | `6e4a739b99` |
| `scripts/21_entry_predictors.R` | `0da0cdec2a` |
| `scripts/23_manifest.R` | `f86c18d89c` |
| `scripts/29_baseline_sensitivity.R` | `da138d9ae0` |
| `scripts/_lib_analysis.R` | `4cee3e0a90` |
| `scripts/25_simulation.R` | `3d58a26cbc` |
| `scripts/13_figures.R` | `32c9dd6bb9` |
| `run_v8_parallel.sh` | `eae5c68f34` |
| `reproduce_current.sh` | `f29fa406cc` |
| `R/utils.R` | `26918b936f` |
| `R/audit_selection.R` | `bdcd9e6323` |
| `R/audit_modules.R` | `d847db14c8` |
| `R/audit_predictions.R` | `2d78a00f30` |
| `scripts/47_nonconstant_diagnostic.R` | `41b10b1c93` |
| `scripts/48_dedup_sensitivity.R` | `c5967b3618` |
| `scripts/49_audit_prediction_api.R` | `e8911d13b7` |
| `R/sigAudit.R` | `8ea470e960` |

## 3. 正式主实验：完成情况

注册网格：40 个（队列×种子）= 8 队列 × 5 划分种子
已完成：**40 / 40**
缺失：无
不在网格内的多余产物：无

| 队列 | 已完成划分种子 |
|---|---|
| MESO | 20261014, 20261016, 20261017, 20261018, 20261019 |
| PAAD | 20261014, 20261016, 20261017, 20261018, 20261019 |
| STAD | 20261014, 20261016, 20261017, 20261018, 20261019 |
| KIRC | 20261014, 20261016, 20261017, 20261018, 20261019 |
| LUAD | 20261014, 20261016, 20261017, 20261018, 20261019 |
| COADREAD | 20261014, 20261016, 20261017, 20261018, 20261019 |
| LGG | 20261014, 20261016, 20261017, 20261018, 20261019 |
| BRCA | 20261014, 20261016, 20261017, 20261018, 20261019 |

## 4. 失败与已修问题

链日志中当前没有未处理的 FAIL/ABORTED 记录。

已知并已修复的一例（保留记录）：
- **LUAD 首次失败**：`OS_STATUS` 含**空字符串**，被「显式取值映射」判为未定义取值而中止。
  修法：空串/纯空白**按缺失处理并计数**（实测 52 例）。链在失败时**按设计立即中止**，
  未把错误分母带进汇总；修复后 LUAD 已重跑成功。

## 5. 尚未完成 / 已知偏离

本节的每一条都由脚本按**文件是否真的存在**生成，不再手写进度。

- [x] 正文成稿（Methods / Results / Discussion）—— 短稿 `docs/paper_cmpb_short.md` 已成稿（摘要 248 词、正文 2,987 词）
- [ ] 投稿期刊、分区体系与 APC 需**官方来源核验**（见 `docs/journal_check.md`）
- [ ] LGG / PAAD 的 PFS 终点**时间起点**未能独立核实（门户元数据未给出）
- [ ] 下载日期与批次失败清单未在数据文件中记录
- [x] `S3`（家族对比）主张已**撤下**，偏离记录见 `config/frozen_analysis_plan.md`
- [x] `grade_credibility()` 已移除 A–D 字母评级；optimism / transport / concordance
      三个模块已标为 **experimental 并移出默认入口**，未参与本文任何结论
- [x] 公开仓库 `https://github.com/123zzm123/sigAudit` 已导出（代码、配置、结果表；
      **不含**全量逐患者预测与原始表达矩阵）
- [ ] 本次修订新增的三个模块（非恒定分层、去重、预测级入口）为**事后诊断**，
      已在正文标注；正式提交与对外发信仍需作者执行

## 6. 当前可提供的文件

- `results/manifest/sample_manifest.tsv` （存在）
- `results/manifest/geneset_manifest.tsv` （存在）
- `results/manifest/candidate_diversity.tsv` （存在）
- `results/manifest/public_signatures_genes.tsv` （存在）
- `results/manifest/public_signatures_meta.tsv` （存在）
- `results/manifest/public_signatures_excluded.tsv` （存在）
- `results/manifest/main_input_files.txt` （存在）
- `results/manifest/run_manifest.tsv` （存在）
- `results/manifest/frozen_v1_20261017.tsv` （存在）
- `config/data_dictionary.md` （存在）
- `config/design_grid.txt` （存在）
- `config/frozen_analysis_plan.md` （存在）
- `results/tables/20_rule_performance_long.tsv` （存在）
- `results/tables/20_rule_members_long.tsv` （存在）
- `results/tables/29_baseline_sensitivity.tsv` （存在）
- `results/tables/29_paired_vs_fixed.tsv` （存在）
- `results/tables/31_metabric_rule_performance.tsv` （存在）
- `results/tables/32_metabric_baseline_sensitivity.tsv` （存在）
- `results/tables/32_metabric_vs_tcga.tsv` （存在）
- `results/figures/fig5_metabric.png` （存在）
- `results/tables/36_metabric_prediction_check.tsv` （存在）
- `results/tables/25_simulation_final_replicates.tsv` （存在）
- `results/tables/47_nonconstant_overall.tsv` （存在）
- `results/tables/48_dedup_summary.tsv` （存在）
- `results/tables/49_audit_api_check.tsv` （存在）
- `docs/paper_cmpb_short.md` （存在）
- `R/audit_predictions.R` （存在）
- `results/figures/fig_captions.md` （存在）
- `reproduce_current.sh` （存在）

## 7. 复现入口

```bash
bash reproduce_current.sh tables      # 主分析表图（由已保存预测重建；不训练）
bash reproduce_current.sh signatures  # 26a+26b：会在签名基因集上【重新拟合】候选池
bash reproduce_current.sh metabric    # METABRIC（31 训练 → 32 → 33）——**会训练**
bash reproduce_current.sh train       # 从原始数据重训 40 组（约 2.5 小时）
bash reproduce_current.sh all         # train + signatures + tables
bash reproduce_current.sh check       # 单元测试 + R CMD check + 哈希固化
SIGAUDIT_ROOT=/path/to/sigAudit bash reproduce_current.sh tables   # 自定义根目录
```

关键词区分：`tables` **不训练任何模型**；`signatures` 与 `train` **会训练**。
根目录通过参数或 `SIGAUDIT_ROOT` 传给所有子进程（脚本内不再写死路径）；
运行记录与脚本哈希写入 `logs/reproduce_current_provenance.txt`。

**异地运行已验证**：在 `D:/bioinfo/sigAudit_repro_test` 独立根目录跑通 `tables`
（7 步全 OK），其产物与本目录**除路径字符串外字节相同**。

> 读取这些 TSV 时请用 `quote = ""`：数据里含 `van 't Veer` 的孤立单引号，
> `read.table` 默认 `quote = "\"'"` 会把它当作引号起点并**吞掉后续换行**
> （实测 6 行只读到 1 行）。

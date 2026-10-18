# 数据字典与来源取证（2026-10-16）

> 应第二轮评审要求建立。**每一条都给出实际验证方式**，不写未经核对的断言。
> 核对脚本：`scripts/_audit_provenance.py`（本地审计）、`scripts/_verify_ci2.R`（评价口径）。

---

## 1. 队列与数据来源

| 分析用名 | cBioPortal study ID | **实际癌种** | 表达谱 ID |
|---|---|---|---|
| BRCA | `brca_tcga_pan_can_atlas_2018` | Breast Invasive Carcinoma | `brca_tcga_pan_can_atlas_2018_rna_seq_v2_mrna` |
| **COADREAD**（原写作 COAD） | `coadread_tcga_pan_can_atlas_2018` | **Colorectal Adenocarcinoma**（结肠 + 直肠） | `coadread_tcga_pan_can_atlas_2018_rna_seq_v2_mrna` |
| KIRC | `kirc_tcga_pan_can_atlas_2018` | Kidney Renal Clear Cell Carcinoma | `..._rna_seq_v2_mrna` |
| LGG | 本地文件（见下） | Brain Lower Grade Glioma | 既有宽矩阵 |
| LUAD | `luad_tcga_pan_can_atlas_2018` | Lung Adenocarcinoma | `..._rna_seq_v2_mrna` |
| MESO | `meso_tcga_pan_can_atlas_2018` | Mesothelioma | `..._rna_seq_v2_mrna` |
| PAAD | `paad_tcga_pan_can_atlas_2018` | Pancreatic Adenocarcinoma | `..._rna_seq_v2_mrna` |
| STAD | `stad_tcga_pan_can_atlas_2018` | Stomach Adenocarcinoma | `..._rna_seq_v2_mrna` |

### ⚠️ COAD 必须改名为 COADREAD

`02_fetch_cohorts.py` 第 27 行把 "COAD" 指向 `coadread_tcga_pan_can_atlas_2018`，
**代码中没有按结肠/直肠诊断类型再筛选**。官方项目名是 Colorectal Adenocarcinoma。

**处理**：正文一律改称 **COADREAD**（或 CRC），保留完整来源。
不为了名字硬删病例。若审查要求纯 COAD，需另有明确筛选代码（本版未做）。
**待补**：表达矩阵内样本的结肠/直肠构成（本轮未能从 API 取到 `CANCER_TYPE_DETAILED`，
下一轮补查）。

### LGG 的数据来源链（已补全，2026-10-17）

LGG 由**另一个本地项目的脚本**产出，但用的是**同一个 cBioPortal API**，
现给出来源链（第三轮审查要求"补全这一链"，现已可完整追溯）：

| 环节 | 具体内容 |
|---|---|
| 生成脚本 | `D:/bioinfo/BRCA_TLS/scripts/15_lgg_broad_genes.py`（sha256 见 provenance 记录） |
| 研究 ID | `lgg_tcga_pan_can_atlas_2018` |
| 表达谱 ID | `lgg_tcga_pan_can_atlas_2018_rna_seq_v2_mrna` |
| 样本列表 | `lgg_tcga_pan_can_atlas_2018_all` |
| 取值接口 | `POST /molecular-profiles/{PROF}/molecular-data/fetch`，传 `entrezGeneIds`，按 symbol 透视成 样本 × 基因 |
| 基因范围 | Hallmark ∪ C2(KEGG_medicus) ∪ C2(Reactome) ∪ C5(GO-BP)（v2024.1 symbols）的并集 → 17,010 基因（与其余 7 队列一致） |
| 生存来源 | `GET /studies/{STUDY}/clinical-data`（`clinicalDataType=PATIENT`，`projection=SUMMARY`）→ 取 `PFS_MONTHS` / `PFS_STATUS` |
| 落盘位置 | `D:/bioinfo/BRCA_TLS/data/cancers/LGG_broad_expr.tsv`、`LGG_broad_surv.tsv`；本项目 `data/cohorts/LGG_expr.tsv` 是其**逐字节副本**（63,930,776 字节，md5 一致） |

**因此 LGG 不是"另一个来源的独立外部队列"**：它与其余 7 个队列同源（TCGA PanCanAtlas 2018，
同一 API、同一谱类型），差别只在于由哪个脚本抓取。**不得把它包装成独立外部验证。**

> 说明：旧脚本在生成矩阵时对每个基因取 `pivot_table(aggfunc="mean")`，
> 即同一 symbol 对应的多个 Entrez 取均值；该步骤**未被本项目修改**。

---

## 2. 表达单位（评审要求确认是否二次 log）

| 队列 | 行×列 | 最小值 | 中位 | 99 分位 | 最大值 | 判断 |
|---|---|---|---|---|---|---|
| BRCA | 1082×17010 | 0.000 | 275.1 | 11957 | 1.84e6 | 原始值 |
| COADREAD | 592×17010 | −0.955 | 347.0 | 12682 | 6.14e5 | 原始值 |
| KIRC | 510×17010 | 0.000 | 285.5 | 12614 | 1.26e6 | 原始值 |
| LGG | 514×17010 | 0.000 | 283.1 | 10975 | 2.55e6 | 原始值 |
| LUAD | 510×17010 | 0.000 | 290.6 | 12662 | 1.29e6 | 原始值 |
| MESO | 87×17010 | 0.000 | 267.5 | 15310 | 2.98e6 | 原始值 |
| PAAD | 177×17010 | 0.000 | 310.5 | 13469 | 6.70e6 | 原始值 |
| STAD | 412×15928 | −0.940 | 366.4 | 13006 | 1.89e6 | 原始值 |

### 表达谱的官方元数据（2026-10-17 由 API 直接取得）

`GET /molecular-profiles/{profileId}` 对全部 8 个队列返回**完全相同的描述**：

```
name        : "mRNA Expression, RSEM (Batch normalized from Illumina HiSeq_RNASeqV2)"
datatype    : CONTINUOUS
alteration  : MRNA_EXPRESSION
```

**这解释了 COADREAD / STAD 出现负值的原因**：该谱是 **RSEM 经批次归一化**
（PanCanAtlas 的 EB++ 批次校正）后的结果，批次校正会产生小幅负值。
**因此不能说它是"原始 counts"**（原始计数不可能为负）——
正确表述是 **"批次归一化后的 RSEM 表达值（线性尺度）"**。

**结论与变换依据**：
- 数值仍是**线性尺度**（中位约 280、最大约 10⁶；若已 log 通常落在 0–20），
  因此 `log2(expr + 1)` 是必要的单次变换，**不存在二次 log**；
- 负值**不截断为 0**（审查明确要求"不要将负值直接截成 0 以求好看"）。
  可直接验证其安全性：全 8 队列的**最小值为 −0.955 > −1**，
  因此 `log2(x + 1)` 对全部实际取值**都有定义**（结果落在约 −4.5 以上）。
- 该谱在不同 study 上的预处理可能不同；本表只证明**本次实际取得的矩阵**如上。

---

## 3. 样本清单审计（评审要求）

**方法**：TCGA 条形码第 14–15 位编码样本类型；逐队列统计。

| 队列 | 样本数 | 完整 ID | 唯一患者 | 重复样本 | 样本类型分布 |
|---|---|---|---|---|---|
| BRCA | 1082 | 1082 | 1082 | 0 | **01 原发实体瘤 100%** |
| COADREAD | 592 | 592 | 592 | 0 | **01 原发实体瘤 100%** |
| KIRC | 510 | 510 | 510 | 0 | **01 原发实体瘤 100%** |
| **LGG** | **514** | **514** | **514** | **0** | **01 原发实体瘤 100%** |
| LUAD | 510 | 510 | 510 | 0 | **01 原发实体瘤 100%** |
| MESO | 87 | 87 | 87 | 0 | **01 原发实体瘤 100%** |
| PAAD | 177 | 177 | 177 | 0 | **01 原发实体瘤 100%** |
| STAD | 412 | 412 | 412 | 0 | **01 原发实体瘤 100%** |

**结论**：**不存在正常组织、转移灶或同一患者多样本的情况**，
因此"截断到 12 字符后保留第一个重复样本"这一操作在实际数据上**没有丢弃任何样本**。
（LGG 于 2026-10-17 补齐核查：矩阵为 样本 × 基因 布局，514 行全部为 `-01`，
条形码形如 `TCGA-CS-4938-01`。）

**样本数与 API 一致性核验**（取覆盖良好的基因 A1BG / entrez=1）：

| 队列 | study `_all` | A1BG 有值样本 | 本地矩阵 |
|---|---|---|---|
| COADREAD | 594 | **592** | 592 ✓ |
| STAD | 440 | **412** | 412 ✓ |
| KIRC | 512 | **510** | 510 ✓ |

→ 本地矩阵样本数与该谱实际有值的样本数**完全一致**，无遗漏、无多余。

---

## 4. 表达值完整性（评审要求：不得静默置零）

| 队列 | 数值缺失 | 负值 | 缺失率 |
|---|---|---|---|
| BRCA / KIRC / LUAD / LGG / MESO / PAAD | **0** | **0** | 0% |
| **COADREAD** | **494,633** | 134,333 | **4.91%** |
| **STAD** | **65,436** | 214,476 | **1.00%** |

**结构**：缺失**按基因聚集**，不是按样本：
- COADREAD：17,010 个基因中约 2,179 个基因在**恰好 227 个样本**上缺值；
  另有 3,177 个基因存在负值（最小 −0.955），同样集中在那 227 个样本
- STAD：1,722 个基因在 38 个样本上缺值；2,789 个基因存在负值（最小 −0.940），
  分布在全部 412 个样本

**来源核实**：直接向 cBioPortal API 重新取这些高缺失基因，
**API 同样只返回 365 个样本的值**（592 − 227 = 365），与本地完全对应；
且取两个样本 × 20 个基因逐值比对，**与本地最大差为 0**。
→ **这是源数据本身的性质，不是下载或拼接引入的错误。**

**⚠️ 必须显式化的处理规则**（旧脚本 `Z[!is.finite(Z)] <- 0` 是**静默置零**，
会把缺失变成"该基因的平均表达"，必须改）：

1. **基因过滤**：在**训练集**上计算每个基因的缺失比例，
   缺失比例超过阈值（拟 20%）的基因**排除出该（队列×划分）的候选宇宙**；
2. **剩余缺失的填补**：用**训练集该基因的中位数**填补；
3. 验证集与测试集原样应用上述规则，**不得用各自集合的统计量**；
4. 每个（队列×划分）记录：过滤掉的基因数、填补的数值个数。

---

## 5. 终点字典（评审要求：证明事件语义，不能只看列名）

| 队列 | 时间列 | 状态列 | **原始取值及频数** | 终点 | 事件编码 |
|---|---|---|---|---|---|
| BRCA | OS_MONTHS | OS_STATUS | `0:LIVING` 933 / `1:DECEASED` 151 | OS | 1 = 死亡 |
| COADREAD | OS_MONTHS | OS_STATUS | `0:LIVING` 471 / `1:DECEASED` 120 / NA 3 | OS | 1 = 死亡 |
| KIRC | OS_MONTHS | OS_STATUS | `0:LIVING` 342 / `1:DECEASED` 170 | OS | 1 = 死亡 |
| LUAD | OS_MONTHS | OS_STATUS | `0:LIVING` 328 / `1:DECEASED` 186 / NA 52 | OS | 1 = 死亡 |
| MESO | OS_MONTHS | OS_STATUS | `1:DECEASED` 74 / `0:LIVING` 13 | OS | 1 = 死亡 |
| STAD | OS_MONTHS | OS_STATUS | `0:LIVING` 269 / `1:DECEASED` 171 / NA 5 | OS | 1 = 死亡 |
| **LGG** | PFS_MONTHS | PFS_STATUS | **`0:CENSORED` 321 / `1:PROGRESSION` 192** / NA 1 | **PFS** | 1 = 进展 |
| **PAAD** | PFS_MONTHS | PFS_STATUS | **`1:PROGRESSION` 109 / `0:CENSORED` 75** | **PFS** | 1 = 进展 |

**结论**：
- 6 个队列是**定义清楚的总生存（OS）**，事件 = 死亡；
- LGG 与 PAAD 的门户列是 `PFS_STATUS` / `PFS_MONTHS`，取值只有 `1:PROGRESSION` / `0:CENSORED`。
- ✅ **更正**：`02_对评审意见的逐条回应.md` 中把 LGG/PAAD 的事件编码写成
  `1:DECEASED / 0:LIVING` 是**我的笔误**，实际为 `1:PROGRESSION / 0:CENSORED`。

### ⚠️ 但"PFS_STATUS"不等于终点语义已被证明（第三轮审查要求补证）

**门户元数据本身说明不足**。API 返回的全部信息只有：

```
PFS_MONTHS : displayName "Progress Free Survival (Months)"   datatype NUMBER
PFS_STATUS : displayName "Progression Free Status"           datatype STRING
```

**没有**时间起点、**没有**事件构成、**没有**与原始资源的映射说明。
因此**不能**仅凭列名判定它是否涵盖"复发、新肿瘤事件、带瘤死亡"。

**可核对的部分（本轮已核，2026-10-17）**：

| 核对项 | LGG | PAAD | 说明 |
|---|---|---|---|
| 该 study 提供的生存终点的**全部**列 | `OS_*`、`DSS_*`、`PFS_*` | 同左 | **没有 PFI / DFI 列** → 不能改称 PFI |
| `PFS_MONTHS ≤ OS_MONTHS` 成立比例 | 512/512 = 100% | 184/184 = 100% | 与"以进展或死亡为事件"一致 |
| **死亡但 PFS 被删失**的患者数 | **9 / 125** | **18 / 99** | **死亡并非自动计为 PFS 事件** |
| 有 PFS 事件但未记录死亡 | 76 | 28 | 与"进展后仍存活"一致 |
| `PFS_MONTHS > OS_MONTHS` 的异常 | 1 例（差 3.91 月） | 0 | 记录为数据质量异常 |

**据实表述（正文必须照此写）**：
> 我们使用 cBioPortal PanCanAtlas 2018 发布的 `PFS_STATUS` / `PFS_MONTHS`，
> 事件标签为 `1:PROGRESSION`。该发布**同时提供 OS 与 DSS**，但**不提供 PFI/DFI 列**，
> 因此我们**不**把它改称 PFI。其门户元数据未说明时间起点与事件构成；我们实证发现
> 有 9/125（LGG）与 18/99（PAAD）的死亡患者仍被记为删失，即**死亡并非一律计为事件**。
> 故本文把它描述为 **"该发布所给的进展终点"**，不主张它等同于 TCGA-CDR 的
> PFS 或 PFI 定义，也不把事件简化为"只有进展"。

**时间起点**：门户元数据未给出；正文标注为"由该发布定义，未在本研究中独立核实"。

### ⚠️ 终点选择的实际依据必须如实写明

`02_fetch_cohorts.py` 第 31 行注释原文：
> `# 各自事件最多的终点（前面筛查过）`

即：每个队列的终点是**按事件数多少挑出来的**（LGG/PAAD 的 PFS 事件多于 OS），
**不是事先的临床预设**。这一条必须如实写进方法，**不能包装成独立的临床预设**。

### 混合终点的影响（评审指出原断言未经证明）

原回应中"所有比较在队列内配对，因此混合终点不影响结论方向"**已撤回** ——
配对只保证同一任务内规则之间终点一致，**不能**保证 OS 与进展终点有相同的事件率、
删失结构与预测信号。

**处理**（通常不需重新拟合）：
- 主表**标注每个队列的终点类型**；
- 提供 **6 个 OS 队列的单独汇总**；
- 2 个进展终点队列（LGG、PAAD）**分别展示**。

---

## 6. 生存表遗留问题

`02_fetch_cohorts.py` 只保留了时间列与状态列，其余临床属性**未落盘**，
因此本轮无法从本地核对 `CANCER_TYPE_DETAILED`（COADREAD 的结肠/直肠构成）。
**待补**：重新拉取并保存完整临床表（至少含癌种细分、分期、诊断年龄），
并保留原始取值以便复核。

---

## 7. 复现所需但尚未固化的内容（诚实列出）

| 项目 | 状态 |
|---|---|
| `msigdb_genesets.json` 的构建脚本 | ✅ `scripts/01b_build_genesets.py`（本版已随项目提供） |
| 表达文件哈希 | ✅ 已计算（见下表）；脚本 `scripts/28_freeze_current.sh` + `_audit_provenance.py` |
| 下载日期 | ⚠️ 可从文件时间戳近似，未在文件中记录（**如实标注为未完成**） |
| 每次下载的批次失败记录 | ⚠️ 原下载脚本失败后重试、最终跳过时不写遗漏清单。本版通过 A1BG 抽样核验确认无遗漏，但机制上仍未补齐 |
| LGG 来源文件的生成脚本 | ✅ 已补全：`D:/bioinfo/BRCA_TLS/scripts/15_lgg_broad_genes.py`（见第 1 节来源链） |
| 表达谱的官方元数据 | ✅ 已由 API 取得（RSEM, batch normalized；见第 2 节） |
| 生存终点的门户元数据 | ✅ 已由 API 取得，并附实证核对（见第 5 节）；**时间起点未能独立核实** |

**表达文件 SHA-256（前 16 位）**

```
LGG       2f457ff04a7291f3        BRCA     3879a01e1892e9fb
COADREAD  a762d9b08bfb8b08        KIRC     c68f28171bb4195a
LUAD      3f94551377587e28        MESO     ef3c4ffc8fc8f46e
PAAD      64995baef1451fb7        STAD     eb8dbfce77ec1a8d
```

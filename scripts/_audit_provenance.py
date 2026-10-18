# -*- coding: utf-8 -*-
"""本地样本清单审计：
 1. 从表达矩阵的样本轴取出完整 sample ID
 2. 按 TCGA 条形码第 14-15 位判定样本类型（01=原发瘤, 11=实体组织正常, 06=转移, 02=复发…）
 3. 统计同一患者是否有多个样本、以及这些样本的类型分布
 4. 报告生存表里的事件取值频数（用于核对终点语义）
不联网，只读本地文件。
"""
import io, os, sys, glob
import pandas as pd
sys.stdout.reconfigure(encoding="utf-8")

OUT = "D:/bioinfo/sigAudit/data/cohorts"
COH = ["BRCA", "COAD", "KIRC", "LUAD", "MESO", "PAAD", "STAD"]

STYPE = {"01": "原发实体瘤", "02": "复发瘤", "03": "原发血癌", "06": "转移瘤",
         "07": "附加转移", "11": "实体组织正常", "12": "血液来源正常", "14": "骨髓正常"}


def sample_axis(f):
    """返回该表达文件里的 sample ID 列表（自动判断样本在行还是列）"""
    hdr = pd.read_csv(f, sep="\t", nrows=3)
    cols = [str(c) for c in hdr.columns]
    # 基因名一般是字母开头的符号；样本 ID 形如 TCGA-xx-xxxx
    is_sample = [c.startswith("TCGA-") for c in cols]
    if sum(is_sample) > len(cols) / 2:
        return cols[1:] if cols[0].lower().startswith(("sample", "unnamed")) else cols
    # 否则样本在行索引
    idx = pd.read_csv(f, sep="\t", usecols=[0], index_col=0).index
    return [str(i) for i in idx]


print("=" * 92)
print("A. 样本类型审计（按 TCGA 条形码第 14-15 位）")
print("=" * 92)
for c in COH:
    f = f"{OUT}/{c}_expr.tsv"
    if not os.path.exists(f):
        print(f"  {c:5s} 缺文件"); continue
    ids = sample_axis(f)
    full = [i for i in ids if len(i) >= 15]
    code = [i[13:15] for i in full]
    dist = pd.Series(code).value_counts().to_dict()
    dist_cn = {f"{k}({STYPE.get(k,'?')})": v for k, v in dist.items()}
    pat = [i[:12] for i in full]
    npat = len(set(pat))
    dup = len(pat) - npat
    n_short = len(ids) - len(full)
    print(f"  {c:5s} 样本 {len(ids):5d}（完整ID {len(full)}，短ID {n_short}）  "
          f"唯一患者 {npat:5d}  重复样本 {dup:4d}")
    print(f"        类型分布: {dist_cn}")
    if dup > 0:
        s = pd.Series(pat); d = s[s.duplicated(keep=False)].value_counts()
        print(f"        有多样本的患者数: {len(d)}，最多一个患者有 {d.max()} 个样本")

print()
print("=" * 92)
print("B. 生存表事件取值频数（用于核对终点语义）")
print("=" * 92)
for c in COH:
    sf = f"{OUT}/{c}_surv.tsv"
    if not os.path.exists(sf):
        print(f"  {c:5s} 缺文件"); continue
    sv = pd.read_csv(sf, sep="\t", index_col=0)
    print(f"  {c:5s} 列: {list(sv.columns)}   n={len(sv)}")
    for col in sv.columns:
        if col.endswith("STATUS") or col.endswith("MONTHS"):
            vc = sv[col].value_counts(dropna=False).to_dict()
            if col.endswith("STATUS"):
                print(f"        {col:14s} 取值频数: {vc}")
            else:
                print(f"        {col:14s} 中位 {pd.to_numeric(sv[col], errors='coerce').median():.1f}  "
                      f"缺失 {int(pd.to_numeric(sv[col], errors='coerce').isna().sum())}")

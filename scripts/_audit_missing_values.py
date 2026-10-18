# -*- coding: utf-8 -*-
import sys, numpy as np, pandas as pd
sys.stdout.reconfigure(encoding="utf-8")
for c in ["COAD", "STAD"]:
    f = f"D:/bioinfo/sigAudit/data/cohorts/{c}_expr.tsv"
    d = pd.read_csv(f, sep="\t", index_col=0)
    v = d.to_numpy(dtype=float)
    nan = np.isnan(v)
    per_row = nan.sum(axis=1)      # 每个样本缺失多少基因
    per_col = nan.sum(axis=0)      # 每个基因缺失多少样本
    n = v.shape[0]
    print("="*80)
    print(f"{c}: {v.shape[0]} 样本 × {v.shape[1]} 基因；总缺失 {nan.sum():,} ({100*nan.mean():.2f}%)")
    print(f"  【按样本看】缺失基因数: 0 的样本 {(per_row==0).sum()}；"
          f"缺失>50%基因的样本 {(per_row>0.5*v.shape[1]).sum()}；"
          f"缺失>90%的样本 {(per_row>0.9*v.shape[1]).sum()}")
    print(f"             最大单样本缺失基因数 {per_row.max()}")
    print(f"  【按基因看】完全无缺失的基因 {(per_col==0).sum()}；"
          f"缺失>90%样本的基因 {(per_col>0.9*n).sum()}；最大 {per_col.max()}")
    # 缺失是否集中在少数基因？
    top = np.sort(per_col)[::-1][:10]
    print(f"             缺失最多的 10 个基因的缺失样本数: {list(top)}")
    # 负值结构
    neg = (v < 0)
    print(f"  【负值】总数 {neg.sum():,}；最小 {np.nanmin(v):.3f}；"
          f"负值占比 {100*neg.sum()/np.isfinite(v).sum():.2f}%")
    print(f"         负值所在基因数 {int((neg.sum(axis=0)>0).sum())}；负值所在样本数 {int((neg.sum(axis=1)>0).sum())}")

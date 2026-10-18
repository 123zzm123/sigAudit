# -*- coding: utf-8 -*-
import sys, numpy as np, pandas as pd
sys.stdout.reconfigure(encoding="utf-8")
OUT = "D:/bioinfo/sigAudit/data/cohorts"
COH = ["BRCA","COAD","KIRC","LUAD","MESO","PAAD","STAD"]
print("表达值分布（判断是否已经 log 变换）")
print(f"{'队列':6s} {'min':>10s} {'q50':>8s} {'q99':>10s} {'max':>12s}  {'若已是log2':>10s}")
for c in COH:
    f = f"{OUT}/{c}_expr.tsv"
    d = pd.read_csv(f, sep="\t", nrows=200)
    vals = d.iloc[:, 1:].to_numpy(dtype=float).ravel()
    vals = vals[np.isfinite(vals)]
    q = np.percentile(vals, [0, 50, 99, 100])
    flag = "是" if q[3] < 25 else ("否(需log)" if q[3] > 100 else "不确定")
    print(f"{c:6s} {q[0]:10.3f} {q[1]:8.3f} {q[2]:10.3f} {q[3]:12.3f}  {flag:>10s}")
# LGG（另一路径，来自 STAR counts）
f = "D:/bioinfo/BRCA_TLS/data/cancers/LGG_broad_expr.tsv"
d = pd.read_csv(f, sep="\t", nrows=200)
vals = d.iloc[:, 1:].to_numpy(dtype=float).ravel()
vals = vals[np.isfinite(vals)]
q = np.percentile(vals, [0, 50, 99, 100])
print(f"{'LGG':6s} {q[0]:10.3f} {q[1]:8.3f} {q[2]:10.3f} {q[3]:12.3f}  {'否(需log)' if q[3]>100 else '是':>10s}")
print()
print("判据：cBioPortal 的 rna_seq_v2_mrna 若已是 log2(RSEM+1)，最大值通常在 5~20；")
print("      原始 RSEM/counts 最大值会是几千到几十万。")

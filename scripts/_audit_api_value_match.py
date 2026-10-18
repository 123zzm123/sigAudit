# -*- coding: utf-8 -*-
import sys, numpy as np, pandas as pd, requests
sys.stdout.reconfigure(encoding="utf-8")
B = "https://www.cbioportal.org/api"
PROF = "coadread_tcga_pan_can_atlas_2018_rna_seq_v2_mrna"; STUDY = "coadread_tcga_pan_can_atlas_2018"
d = pd.read_csv("D:/bioinfo/sigAudit/data/cohorts/COAD_expr.tsv", sep="\t", index_col=0)
v = d.to_numpy(dtype=float)
per_col_nan = np.isnan(v).sum(axis=0)
worst = np.argsort(per_col_nan)[::-1][:8]
genes = [d.columns[i] for i in worst]
print("COAD 缺失最多的 8 个基因:", genes, " 缺失样本数:", [int(per_col_nan[i]) for i in worst])
r = requests.post(f"{B}/genes/fetch?geneIdType=HUGO_GENE_SYMBOL", json=genes, timeout=120)
e2s = {x["hugoGeneSymbol"]: x["entrezGeneId"] for x in r.json()}
ent = [e2s[g] for g in genes if g in e2s]
r = requests.post(f"{B}/molecular-profiles/{PROF}/molecular-data/fetch",
                  json={"sampleListId": f"{STUDY}_all", "entrezGeneIds": ent}, timeout=600)
r.raise_for_status()
api = pd.DataFrame(r.json()); api["gene"] = api["entrezGeneId"].map({v_: k for k, v_ in e2s.items()})
w = api.pivot_table(index="sampleId", columns="gene", values="value", aggfunc="mean")
print(f"  API 返回样本数: {len(w)}  （该 profile 在 _all 里有的样本）")
for g in genes[:4]:
    if g not in w.columns: print(f"  {g}: API 无此基因"); continue
    a = w[g]; loc = d[g]
    common = loc.index.intersection(a.index)
    loc_nan = int(loc.isna().sum()); api_nan = int(a.isna().sum())
    print(f"  {g:12s} 本地缺失 {loc_nan:4d}/{len(loc)} | API 缺失 {api_nan:4d}/{len(a)} | "
          f"共有样本 {len(common)}，其中本地缺失 {int(loc[common].isna().sum())}，API 非缺失 {int(a[common].notna().sum())}")
print()
print("判读：若 API 在该基因上也大量缺失/为负，则源数据如此，不是下载问题。")

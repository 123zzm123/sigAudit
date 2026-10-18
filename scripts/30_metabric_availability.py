# -*- coding: utf-8 -*-
"""
30 — METABRIC 数据可用性检查（只读，不训练）

对应 docs/metabric_availability_check.md。用于在任何建模之前回答三个问题：
  1. 数据能否取得、样本量/事件数是否足够；
  2. 六个发表签名的基因覆盖是否达到门槛（>=8）；
  3. 表达尺度与终点字段是否可用。

输出：results/manifest/metabric_availability.tsv + 控制台报告
用法： python scripts/30_metabric_availability.py
"""
import io, os, json, datetime
import requests
import pandas as pd
import numpy as np

ROOT = os.environ.get("SIGAUDIT_ROOT", "D:/bioinfo/sigAudit")
MAN = os.path.join(ROOT, "results", "manifest")
B = "https://www.cbioportal.org/api"
SID, PROF = "brca_metabric", "brca_metabric_mrna"

def g(path, **p):
    r = requests.get(B + path, params=p, timeout=180); r.raise_for_status(); return r.json()
def post(path, body, timeout=600):
    r = requests.post(B + path, json=body, timeout=timeout); r.raise_for_status(); return r.json()

print("=" * 68)
print("1) 研究与分子谱")
print("=" * 68)
studies = [s for s in g("/studies") if "metabric" in s["studyId"].lower()]
print(f"  study: {studies[0]['studyId']}  all samples: {studies[0]['allSampleCount']}")
profs = g(f"/studies/{SID}/molecular-profiles")
p = [x for x in profs if x["molecularProfileId"] == PROF][0]
print(f"  profile: {p['molecularProfileId']}")
print(f"  name   : {p['name']}")
print(f"  type   : {p['molecularAlterationType']} / {p['datatype']}")

print()
print("=" * 68)
print("2) 表达样本数与样本类型")
print("=" * 68)
sl = g(f"/sample-lists/{SID}_mrna")
print(f"  {SID}_mrna  sampleCount = {sl.get('sampleCount')}")
cs = g(f"/studies/{SID}/clinical-data", clinicalDataType="SAMPLE", projection="SUMMARY")
ws = pd.DataFrame(cs).pivot_table(index="sampleId", columns="clinicalAttributeId",
                                  values="value", aggfunc="first")
print("  SAMPLE_TYPE:")
print("    " + ws["SAMPLE_TYPE"].value_counts(dropna=False).to_string().replace("\n", "\n    "))

print()
print("=" * 68)
print("3) 终点字段")
print("=" * 68)
cl = g(f"/studies/{SID}/clinical-data", clinicalDataType="PATIENT", projection="SUMMARY")
wp = pd.DataFrame(cl).pivot_table(index="patientId", columns="clinicalAttributeId",
                                  values="value", aggfunc="first")
om = pd.to_numeric(wp["OS_MONTHS"], errors="coerce")
osst = wp["OS_STATUS"].astype(str)
print(f"  患者 {len(wp)}；OS_MONTHS 非缺失 {int(om.notna().sum())}")
print("  OS_STATUS:")
print("    " + wp["OS_STATUS"].value_counts(dropna=False).to_string().replace("\n", "\n    "))
print(f"  OS_MONTHS: min {om.min():.1f}  median {om.median():.1f}  max {om.max():.1f}")

print()
print("=" * 68)
print("4) 六个签名的基因覆盖（用 HGNC 规范化后的符号）")
print("=" * 68)
mapped = pd.read_csv(os.path.join(MAN, "public_signatures_genes_mapped.tsv"), sep="\t")
SIGS = {k: sorted(set(v["gene_used"])) for k, v in mapped.groupby("signature")}
allg = sorted({x for v in SIGS.values() for x in v})
m = post("/genes/fetch?geneIdType=HUGO_GENE_SYMBOL", allg, timeout=180)
sym2ent = {x["hugoGeneSymbol"]: x["entrezGeneId"] for x in m}
d = post(f"/molecular-profiles/{PROF}/molecular-data/fetch",
         {"sampleListId": f"{SID}_mrna", "entrezGeneIds": sorted(sym2ent.values())})
df = pd.DataFrame(d)
present = set(df["entrezGeneId"].unique())
rows = []
for k, gs in SIGS.items():
    hit = [x for x in gs if x in sym2ent and sym2ent[x] in present]
    rows.append(dict(signature=k, n_source=len(gs), n_available=len(hit),
                     pass_threshold=int(len(hit) >= 8)))
    print(f"  {k:20s} {len(hit):3d} / {len(gs):3d}  {'OK' if len(hit) >= 8 else '**不足 8**'}")
print(f"  表达矩阵: {df['entrezGeneId'].nunique()} 基因 × {df['sampleId'].nunique()} 样本")
v = pd.to_numeric(df["value"], errors="coerce")
print(f"  取值: min {v.min():.3f}  median {v.median():.3f}  max {v.max():.3f}  负值 {int((v<0).sum())}")
print("  尺度：按发布值使用，不再取对数（见 data/cohorts/METABRIC_provenance.json）")

print()
print("=" * 68)
print("5) 入排（与 35_rebuild_metabric.py 一致：只排除血管肉瘤 PBS）")
print("=" * 68)
surv = pd.DataFrame(index=sorted(df["sampleId"].unique())); surv.index.name = "sampleId"
s2p = dict(zip(df["sampleId"], df["patientId"]))
surv["patientId"] = [s2p[s] for s in surv.index]
surv["OS_MONTHS"] = surv["patientId"].map(lambda p_: wp.at[p_, "OS_MONTHS"] if p_ in wp.index else np.nan)
surv["OS_STATUS"] = surv["patientId"].map(lambda p_: wp.at[p_, "OS_STATUS"] if p_ in wp.index else np.nan)
surv["OT"] = [ws.at[s, "ONCOTREE_CODE"] if s in ws.index else None for s in surv.index]
surv["OS_MONTHS"] = pd.to_numeric(surv["OS_MONTHS"], errors="coerce")
surv["STYPE"] = [ws.at[s, "SAMPLE_TYPE"] if s in ws.index else None for s in surv.index]
ok = ((surv["STYPE"] == "Primary") & surv["OS_MONTHS"].notna() & (surv["OS_MONTHS"] > 0)
      & surv["OS_STATUS"].isin(["1:DECEASED", "0:LIVING"]) & (surv["OT"] != "PBS"))
ev = int((surv.loc[ok, "OS_STATUS"] == "1:DECEASED").sum())
print(f"  表达样本 {len(surv)}；排除血管肉瘤 {int((surv['OT']=='PBS').sum())}；"
      f"**纳入 {int(ok.sum())}，事件 {ev}（{100*ev/int(ok.sum()):.1f}%）**")

out = pd.DataFrame(rows)
out.to_csv(os.path.join(MAN, "metabric_availability.tsv"), sep="\t", index=False)
print(f"\n清单写入 {MAN}/metabric_availability.tsv")
print("检查日期:", datetime.date.today().isoformat())

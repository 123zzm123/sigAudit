# -*- coding: utf-8 -*-
"""
35 — 重建 METABRIC 表达与生存表（第五轮审查后的两处上游修正）

改动 1（基因标识规范化）：改用 public_signatures_genes_mapped.tsv 的【批准符号】，
   重新从 cBioPortal 取表达值。原下载只覆盖源列表符号，因此映射后的 NDC80/CDK1 等
   根本不在矩阵里——不重新取数就得不到映射的收益。

改动 2（组织学筛选修正）：化生性乳腺癌（MBC）**本身就是浸润性乳腺癌**，不属于
   "非浸润性/非癌"。原筛选把它与血管肉瘤一起排除，理由是错的。
   现在只排除血管肉瘤（PBS，属肉瘤），保留 MBC。
   原版解释见 config/metabric_frozen_config.md 的修订记录。

输出：
   data/cohorts/METABRIC_expr.tsv   基因 × 样本（批准符号）
   data/cohorts/METABRIC_surv.tsv
   data/cohorts/METABRIC_provenance.json

用法： python scripts/35_rebuild_metabric.py
"""
import io, os, json, hashlib, datetime
import requests
import pandas as pd
import numpy as np

ROOT = os.environ.get("SIGAUDIT_ROOT", "D:/bioinfo/sigAudit")
OUT = os.path.join(ROOT, "data", "cohorts")
MAN = os.path.join(ROOT, "results", "manifest")
B = "https://www.cbioportal.org/api"
SID, PROF = "brca_metabric", "brca_metabric_mrna"

def g(path, **p):
    r = requests.get(B + path, params=p, timeout=180); r.raise_for_status(); return r.json()
def post(path, body, timeout=600):
    r = requests.post(B + path, json=body, timeout=timeout); r.raise_for_status(); return r.json()

mapped = pd.read_csv(os.path.join(MAN, "public_signatures_genes_mapped.tsv"), sep="\t")
genes = sorted(set(mapped["gene_used"]))
print("需要下载的批准符号:", len(genes))

m = post("/genes/fetch?geneIdType=HUGO_GENE_SYMBOL", genes, timeout=180)
sym2ent = {x["hugoGeneSymbol"]: x["entrezGeneId"] for x in m}
ent2sym = {v: k for k, v in sym2ent.items()}
print("  映射到 Entrez:", len(sym2ent), "/", len(genes))
missing_hgnc = [x for x in genes if x not in sym2ent]
if missing_hgnc:
    print("  未找到 Entrez 的符号:", missing_hgnc)

d = post(f"/molecular-profiles/{PROF}/molecular-data/fetch",
         {"sampleListId": f"{SID}_mrna", "entrezGeneIds": sorted(sym2ent.values())})
df = pd.DataFrame(d)
df["symbol"] = df["entrezGeneId"].map(ent2sym)
print("  返回记录:", len(df), " 唯一基因:", df["symbol"].nunique(),
      " 唯一样本:", df["sampleId"].nunique())

mat = df.pivot_table(index="symbol", columns="sampleId", values="value", aggfunc="mean").sort_index()
print("表达矩阵:", mat.shape)

# ---- 临床与入排 ----
cl = g(f"/studies/{SID}/clinical-data", clinicalDataType="PATIENT", projection="SUMMARY")
wp = pd.DataFrame(cl).pivot_table(index="patientId", columns="clinicalAttributeId",
                                  values="value", aggfunc="first")
cs = g(f"/studies/{SID}/clinical-data", clinicalDataType="SAMPLE", projection="SUMMARY")
ws = pd.DataFrame(cs).pivot_table(index="sampleId", columns="clinicalAttributeId",
                                  values="value", aggfunc="first")

samples = list(mat.columns)
surv = pd.DataFrame(index=samples); surv.index.name = "sampleId"
surv["patientId"] = [df.loc[df["sampleId"] == s, "patientId"].iloc[0] for s in samples]
for c in ("OS_MONTHS", "OS_STATUS"):
    surv[c] = surv["patientId"].map(lambda p: wp.at[p, c] if p in wp.index and c in wp.columns else np.nan)
for c in ("SAMPLE_TYPE", "CANCER_TYPE_DETAILED", "ONCOTREE_CODE", "CANCER_TYPE"):
    surv[c] = [ws.at[s, c] if s in ws.index and c in ws.columns else np.nan for s in samples]
surv["OS_MONTHS"] = pd.to_numeric(surv["OS_MONTHS"], errors="coerce")

# 【修正】只排除血管肉瘤（PBS，属肉瘤）；保留化生性乳腺癌（MBC，属浸润性乳腺癌）
EXCL = {"PBS"}
surv["include"] = (
    (surv["SAMPLE_TYPE"] == "Primary")
    & surv["OS_MONTHS"].notna() & (surv["OS_MONTHS"] > 0)
    & surv["OS_STATUS"].isin(["1:DECEASED", "0:LIVING"])
    & ~surv["ONCOTREE_CODE"].isin(EXCL)
)
inc = surv[surv["include"]]
ev = int((inc["OS_STATUS"] == "1:DECEASED").sum())
print("\n=== 入排（修正后）===")
print("  表达样本:", len(surv))
print("  排除血管肉瘤(PBS):", int(surv["ONCOTREE_CODE"].isin(EXCL).sum()))
print("  保留化生性癌(MBC):", int((surv["ONCOTREE_CODE"] == "MBC").sum()), "（其中有 OS 的:",
      int(((surv["ONCOTREE_CODE"] == "MBC") & surv["include"]).sum()), "）")
print("  **最终纳入:", len(inc), " 事件:", ev, f"({100*ev/len(inc):.1f}%)**")

use = list(inc.index)
mat2 = mat[use]
mat2.to_csv(os.path.join(OUT, "METABRIC_expr.tsv"), sep="\t")
inc.drop(columns=["include"]).to_csv(os.path.join(OUT, "METABRIC_surv.tsv"), sep="\t")

def md5(p): return hashlib.md5(open(p, "rb").read()).hexdigest()
prov = {
    "study": SID, "profile": PROF,
    "profile_name": "mRNA expression (Illumina HT-12 v3 microarray)",
    "download_date": datetime.date.today().isoformat(),
    "revision": "round-5: HGNC-mapped symbols + metaplastic carcinoma re-included",
    "genes_in_matrix": int(mat2.shape[0]), "samples_in_matrix": int(mat2.shape[1]),
    "events": ev,
    "excluded_histology": sorted(EXCL),
    "retained_rare_histology": "MBC (metaplastic carcinoma, an invasive carcinoma)",
    "value_scale": "values used as released by cBioPortal; no additional log transform applied",
    "files": {
        "data/cohorts/METABRIC_expr.tsv": md5(os.path.join(OUT, "METABRIC_expr.tsv")),
        "data/cohorts/METABRIC_surv.tsv": md5(os.path.join(OUT, "METABRIC_surv.tsv")),
    },
}
io.open(os.path.join(OUT, "METABRIC_provenance.json"), "w", encoding="utf-8").write(
    json.dumps(prov, ensure_ascii=False, indent=2))
print("\n写入 provenance:", prov["files"])

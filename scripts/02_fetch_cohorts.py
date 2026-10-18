#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
02 — 下载多癌种队列（宽基因矩阵 + 生存），供大规模实验用

为什么需要多癌种：本研究的核心是「模型选择规则的偏倚」，
              需要在【信号强】和【信号弱】的队列上都验证，
              否则结论可能只适用于某些数据。
              这是"结论能不能推广"的关键，不是锦上添花。

基因范围：与 LGG 一致的 Hallmark∪C2∪C5 并集（约 17k），保证跨队列可比
"""
import os, time, json
import requests, pandas as pd, numpy as np

B = "https://www.cbioportal.org/api"
OUT = "D:/bioinfo/sigAudit/data/cohorts"
os.makedirs(OUT, exist_ok=True)

# 覆盖不同信号强度的癌种
CAND = {
    "LGG":  "lgg_tcga_pan_can_atlas_2018",     # 强信号（已有，跳过下载）
    "MESO": "meso_tcga_pan_can_atlas_2018",    # 最强信号 0.72
    "KIRC": "kirc_tcga_pan_can_atlas_2018",    # 中等偏强
    "STAD": "stad_tcga_pan_can_atlas_2018",    # 中等
    "LUAD": "luad_tcga_pan_can_atlas_2018",    # 中等
    "COAD": "coadread_tcga_pan_can_atlas_2018",# 中等
    "PAAD": "paad_tcga_pan_can_atlas_2018",    # 中等
    "BRCA": "brca_tcga_pan_can_atlas_2018",    # 弱信号
}
# 各自事件最多的终点（前面筛查过）
ENDPOINT = {
    "LGG":  ("PFS_MONTHS","PFS_STATUS"),
    "MESO": ("OS_MONTHS","OS_STATUS"),
    "KIRC": ("OS_MONTHS","OS_STATUS"),
    "STAD": ("OS_MONTHS","OS_STATUS"),
    "LUAD": ("OS_MONTHS","OS_STATUS"),
    "COAD": ("OS_MONTHS","OS_STATUS"),
    "PAAD": ("PFS_MONTHS","PFS_STATUS"),
    "BRCA": ("OS_MONTHS","OS_STATUS"),
}

def log(*a): print(time.strftime("%H:%M:%S"), *a, flush=True)

# 基因表：必须用与 LGG 一致的【四个 GMT 的并集】(~17k)，
# 而不是 hallmark-only (~4.4k)，否则跨队列基因覆盖不一致、不可比。
GS = "D:/bioinfo/BRCA_TLS/data/genesets"
_files = ["h.all.v2024.1.Hs.symbols.gmt",
          "c2.cp.kegg_medicus.v2024.1.Hs.symbols.gmt",
          "c2.cp.reactome.v2024.1.Hs.symbols.gmt",
          "c5.go.bp.v2024.1.Hs.symbols.gmt"]
_syms = set()
for _f in _files:
    _p = os.path.join(GS, _f)
    if not os.path.exists(_p): continue
    for line in open(_p, encoding="utf-8"):
        _syms.update(g for g in line.rstrip("\n").split("\t")[2:] if g)
log(f"GMT 并集基因: {len(_syms)}")
_r = requests.post(f"{B}/genes/fetch?geneIdType=HUGO_GENE_SYMBOL", json=sorted(_syms), timeout=180)
ent = sorted({x["entrezGeneId"] for x in _r.json()})
log(f"映射 Entrez: {len(ent)}")

def fetch(study, prof, ids, batch=60):
    frames = []
    for i in range(0, len(ids), batch):
        for a in range(4):
            try:
                r = requests.post(f"{B}/molecular-profiles/{prof}/molecular-data/fetch",
                    json={"sampleListId": f"{study}_all", "entrezGeneIds": ids[i:i+batch]}, timeout=600)
                r.raise_for_status(); d = r.json()
                if isinstance(d, list) and d:
                    frames.append(pd.DataFrame(d)[["sampleId","entrezGeneId","value"]])
                break
            except Exception:
                time.sleep(3)
    if not frames: return None
    raw = pd.concat(frames, ignore_index=True)
    wide = raw.pivot_table(index="sampleId", columns="entrezGeneId", values="value", aggfunc="mean")
    # 关键：cBioPortal 返回 Entrez ID，但 MSigDB 基因集用 HGNC 符号 ——
    # 必须在这里就转成符号，否则下游与基因集的交集恒为 0（本项目实际踩过）
    _e2s = {x["entrezGeneId"]: x["hugoGeneSymbol"] for x in _r.json()}
    wide.columns = [_e2s.get(int(c), str(c)) for c in wide.columns]
    wide = wide.loc[:, ~pd.Index(wide.columns).duplicated()]
    return wide

for nm, study in CAND.items():
    ef = f"{OUT}/{nm}_expr.tsv"; sf = f"{OUT}/{nm}_surv.tsv"
    if os.path.exists(ef) and os.path.exists(sf):
        log(f"[cache] {nm}"); continue
    log(f"=== {nm} ({study}) ===")
    try:
        profs = requests.get(f"{B}/studies/{study}/molecular-profiles", timeout=120).json()
        mrna = [p["molecularProfileId"] for p in profs
                if p.get("molecularAlterationType")=="MRNA_EXPRESSION" and "Zscore" not in p["molecularProfileId"]]
        if not mrna: log("  无 mRNA 谱"); continue
        ex = fetch(study, mrna[0], ent)
        if ex is None: log("  下载失败"); continue
        ex.to_csv(ef, sep="\t")
        log(f"  表达: {ex.shape}")
        # 生存
        r = requests.get(f"{B}/studies/{study}/clinical-data",
            params={"clinicalDataType":"PATIENT","projection":"SUMMARY"}, timeout=600)
        w = pd.DataFrame(r.json()).pivot_table(index="patientId", columns="clinicalAttributeId",
                                               values="value", aggfunc="first")
        tc, ec = ENDPOINT[nm]
        keep = [c for c in [tc, ec] if c in w.columns]
        w[keep].to_csv(sf, sep="\t")
        log(f"  生存: {w.shape}, 终点 {tc}")
    except Exception as e:
        log(f"  ERR {str(e)[:80]}")

log("全部完成")

# -*- coding: utf-8 -*-
import sys, pandas as pd, requests
sys.stdout.reconfigure(encoding="utf-8")
B = "https://www.cbioportal.org/api"
for study, prof, nm, local in [
    ("coadread_tcga_pan_can_atlas_2018","coadread_tcga_pan_can_atlas_2018_rna_seq_v2_mrna","COADREAD","COAD"),
    ("stad_tcga_pan_can_atlas_2018","stad_tcga_pan_can_atlas_2018_rna_seq_v2_mrna","STAD","STAD"),
    ("kirc_tcga_pan_can_atlas_2018","kirc_tcga_pan_can_atlas_2018_rna_seq_v2_mrna","KIRC","KIRC")]:
    d = pd.read_csv(f"D:/bioinfo/sigAudit/data/cohorts/{local}_expr.tsv", sep="\t", usecols=[0])
    loc = set(str(x)[:12] for x in d.iloc[:,0])
    n_all = len(set(requests.get(f"{B}/sample-lists/{study}_all/sample-ids", timeout=300).json()))
    n_all12 = len(set(str(x)[:12] for x in requests.get(f"{B}/sample-lists/{study}_all/sample-ids", timeout=300).json()))
    r = requests.post(f"{B}/molecular-profiles/{prof}/molecular-data/fetch",
                      json={"sampleListId": f"{study}_all", "entrezGeneIds": [1]}, timeout=600)
    r.raise_for_status()
    apiset = set(str(x["sampleId"])[:12] for x in r.json())
    print(f"{nm:8s}  _all 样本 {n_all}（患者级 {n_all12}）  A1BG(entrez=1) 有值样本 {len(apiset)}"
          f"   本地 {len(loc)}")
    print(f"          本地有、API(A1BG) 无: {len(loc-apiset)}    本地无、API 有: {len(apiset-loc)}")

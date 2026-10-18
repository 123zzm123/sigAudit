#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
04 — 把队列表达矩阵的基因标识符从 Entrez ID 换回 HGNC 基因符号

背景：cBioPortal 的分子数据接口只认 Entrez ID；而 MSigDB 基因集用 HGNC 符号。
      两者交集为 0 时，实验会把每一个基因集都"安全跳过"——表现为脚本几十秒
      跑完并打印"没有完成的基因集！"，极容易被误当成正常完成。

⚠️ 本脚本旧版有一处致命假设错误（2026-10-16 修正）：
   旧版只看 d.columns 是否是数字来判断"是否需要转换"，即假设基因在【列】上。
   但当前 *_expr.tsv 的布局是【行】=基因、列=样本，于是 BRCA 被误判成
   "已是符号，跳过"，一直没被修，导致 BRCA 全程空跑。
   现在两种布局都处理。

用法：python scripts/04_entrez2symbol.py
"""
import os
import glob
import shutil
import sys

import pandas as pd
import requests

B = "https://www.cbioportal.org/api"
OUT = "D:/bioinfo/sigAudit/data/cohorts"
GS = os.environ.get("MSIGDB_DIR", os.path.join(os.environ.get("SIGAUDIT_ROOT", "."), "data", "genesets"))
FILES = ["h.all.v2024.1.Hs.symbols.gmt",
         "c2.cp.kegg_medicus.v2024.1.Hs.symbols.gmt",
         "c2.cp.reactome.v2024.1.Hs.symbols.gmt",
         "c5.go.bp.v2024.1.Hs.symbols.gmt"]


def log(*a):
    print(*a)
    sys.stdout.flush()


def fetch_map():
    """用 MSigDB 里真正出现的符号，去 cBioPortal 反查 Entrez -> 符号映射表。"""
    syms = set()
    for f in FILES:
        p = os.path.join(GS, f)
        if not os.path.exists(p):
            continue
        for line in open(p, encoding="utf-8"):
            syms.update(g for g in line.rstrip("\n").split("\t")[2:] if g)
    syms = sorted(syms)
    log(f"MSigDB 基因符号: {len(syms)}")
    r = requests.post(f"{B}/genes/fetch?geneIdType=HUGO_GENE_SYMBOL",
                      json=syms, timeout=300)
    r.raise_for_status()
    e2s = {str(x["entrezGeneId"]): x["hugoGeneSymbol"] for x in r.json()}
    log(f"Entrez->Symbol 映射表: {len(e2s)}")
    return e2s


def is_entrez(labels):
    s = [str(x) for x in labels[:50]]
    return bool(s) and all(x.isdigit() for x in s)


def convert(labels, e2s):
    """返回 (保留位置, 新标签)。无法映射的位置被丢掉。"""
    keep, new = [], []
    for i, c in enumerate(labels):
        c = str(c)
        if c.isdigit():
            g = e2s.get(c)
            if g is None:
                continue
        else:
            g = c
        keep.append(i)
        new.append(g)
    return keep, new


def main():
    e2s = fetch_map()
    for f in sorted(glob.glob(os.path.join(OUT, "*_expr.tsv"))):
        name = os.path.basename(f)
        d = pd.read_csv(f, sep="\t", index_col=0)
        rows_are_genes = is_entrez(d.index)
        cols_are_genes = (not rows_are_genes) and is_entrez(d.columns)

        if rows_are_genes:
            keep, new = convert(d.index, e2s)
            if not keep:
                log(f"  {name:26s} ★ 行名是 Entrez 但一个都没映射上，保持原样")
                continue
            d = d.iloc[keep, :]
            d.index = new
            axis = "行"
        elif cols_are_genes:
            keep, new = convert(d.columns, e2s)
            if not keep:
                log(f"  {name:26s} ★ 列名是 Entrez 但一个都没映射上，保持原样")
                continue
            d = d.iloc[:, keep]
            d.columns = new
            axis = "列"
        else:
            log(f"  {name:26s} {d.shape[0]} 基因 x {d.shape[1]} 样本 —— 已是符号，跳过")
            continue

        before = d.shape[0] if axis == "行" else d.shape[1]
        # 注意：DataFrame 的 d[mask] 选的是【行】。列方向必须用 d.loc[:, mask]，
        # 否则布尔掩码长度（列数）和行数不匹配 -> ValueError。
        if axis == "行":
            d = d[~pd.Index(d.index).duplicated()]
        else:
            d = d.loc[:, ~pd.Index(d.columns).duplicated()]
        after = d.shape[0] if axis == "行" else d.shape[1]
        shutil.copy2(f, f + ".entrez.bak")       # 留一份原文件，别把数据弄丢
        d.to_csv(f, sep="\t")
        log(f"  {name:26s} {axis}名 Entrez -> {after} 个符号基因"
            f"（去重丢弃 {before - after}）")
    log("完成")


if __name__ == "__main__":
    main()

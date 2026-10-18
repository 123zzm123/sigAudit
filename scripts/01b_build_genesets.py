#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
01b — 【重建】MSigDB 基因集 JSON 的构建脚本（原构建脚本未随项目保留）

第二轮评审指出："MSigDB JSON 的构建脚本没有随包提供，四个 GMT 并集用于表达下载，
但报告还包含 ImmuneSigDB，需说明基因宇宙与基因集来源的关系及交集覆盖。"

本脚本把该构建过程显式化，并**逐字节核对**能否复现原有的 msigdb_genesets.json。

经逆向核对（见下方 SOURCES），原 JSON = 下列 6 个 GMT 的并集，并保留基因数 ≥ 10 的集合：
  h.all            Hallmark
  c2.cp            Canonical pathways（含 KEGG / WikiPathways / BioCarta / PID）
  c2.cp.kegg_medicus
  c2.cp.reactome
  c5.go.bp         GO Biological Process
  c7.immunesigdb   ImmuneSigDB

⚠️ 关键事实（评审要求的说明）：
  用于【定义表达矩阵基因宇宙】的是另外 4 个 GMT（h.all / kegg_medicus / reactome / c5.go.bp），
  其并集约 18,604 个符号。而基因集本身还来自 c2.cp 与 c7.immunesigdb，
  **因此部分基因集会包含不在表达矩阵中的基因** —— 覆盖率必须单独报告，
  不能假设基因集被完整评估。

用法：
  python scripts/01b_build_genesets.py            # 重建并核对
  python scripts/01b_build_genesets.py --write    # 核对通过后写回 JSON
"""
import glob
import io
import json
import os
import sys

GS_DIR = os.environ.get("MSIGDB_DIR", os.path.join(os.environ.get("SIGAUDIT_ROOT", "."), "data", "genesets"))
OUT = os.path.join(os.environ.get("SIGAUDIT_ROOT", "."), "data", "signatures", "msigdb_genesets.json")

SOURCES = [
    "h.all.v2024.1.Hs.symbols.gmt",
    "c2.cp.v2024.1.Hs.symbols.gmt",
    "c2.cp.kegg_medicus.v2024.1.Hs.symbols.gmt",
    "c2.cp.reactome.v2024.1.Hs.symbols.gmt",
    "c5.go.bp.v2024.1.Hs.symbols.gmt",
    "c7.immunesigdb.v2024.1.Hs.symbols.gmt",
]
MIN_GENES = 10


def build():
    sets = {}
    for name in SOURCES:
        p = os.path.join(GS_DIR, name)
        if not os.path.exists(p):
            raise SystemExit("缺少 GMT 文件: " + p)
        n = 0
        for line in io.open(p, encoding="utf-8"):
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            genes = sorted({g for g in parts[2:] if g})
            if len(genes) >= MIN_GENES:
                sets[parts[0]] = genes
            n += 1
        print(f"  {name:48s} 读入 {n:6d} 个集合")
    return sets


def main():
    print("按下列 6 个 GMT 重建基因集库（保留基因数 >= %d）：" % MIN_GENES)
    new = build()
    print(f"\n重建得到 {len(new)} 个基因集")

    if not os.path.exists(OUT):
        print("原 JSON 不存在，直接写出")
    else:
        old = json.load(io.open(OUT, encoding="utf-8"))
        old = {k: sorted(v) for k, v in old.items()}
        print(f"原 JSON 含 {len(old)} 个基因集")
        only_new = set(new) - set(old)
        only_old = set(old) - set(new)
        diff_genes = [k for k in (set(new) & set(old)) if new[k] != old[k]]
        print(f"  重建有、原无: {len(only_new)}")
        print(f"  原有、重建无: {len(only_old)}")
        print(f"  同名但基因集合不同: {len(diff_genes)}")
        ok = not (only_new or only_old or diff_genes)
        print("\n逐集合核对结果:", "完全一致（构建过程已复原）" if ok else "不一致")
        if not ok:
            if only_new: print("   例（重建有原无）:", sorted(only_new)[:3])
            if only_old: print("   例（原有重建无）:", sorted(only_old)[:3])
            if diff_genes: print("   例（基因不同）:", sorted(diff_genes)[:3])
        if not ok and "--write" not in sys.argv:
            print("\n未通过核对，不写回。（加 --write 可强制写回）")
            return

    if "--write" in sys.argv:
        os.makedirs(os.path.dirname(OUT), exist_ok=True)
        json.dump(new, io.open(OUT, "w", encoding="utf-8"), ensure_ascii=False)
        print("\n已写回:", OUT)


if __name__ == "__main__":
    main()

# -*- coding: utf-8 -*-
"""
34 — 把已发表签名里的旧基因符号映射到 HGNC 当前批准符号

为什么需要（第五轮审查）：
  agent 曾以"避免引入自由度"为由不做别名映射。这个理由不成立——
  把旧符号映射到当前符号是**数据标识规范化**，不是根据生存结局挑选特征。
  HGNC 维护 approved / previous / alias 三套符号并公开查询接口。

做法：
  对每个"在表达矩阵里找不到"的源符号，依次查 HGNC 的 prev_symbol 与 alias_symbol；
  命中即记录映射依据；仍找不到的记为 unmapped（不根据模型表现决定去留）。

输出：
  results/manifest/signature_symbol_map.tsv        源符号 -> 批准符号 + 依据
  results/manifest/public_signatures_genes_mapped.tsv  每个（签名 × 基因）的最终使用符号
  results/manifest/signature_symbol_map_report.md  人可读的核对表

用法： python scripts/34_map_signature_symbols.py
"""
import io, os, time
import requests
import pandas as pd

ROOT = os.environ.get("SIGAUDIT_ROOT", "D:/bioinfo/sigAudit")
MAN = os.path.join(ROOT, "results", "manifest")
H = {"Accept": "application/json"}
BASES = ["https://rest.genenames.org", "https://rest.genenames.org"]

src = pd.read_csv(os.path.join(MAN, "public_signatures_genes.tsv"), sep="\t")
# 用 METABRIC 与 TCGA-BRCA 两个矩阵的并集作为"可测量"的参考
def rownames(p):
    """取首表头与首列内容的并集：不同队列布局不同
    （BRCA 是 样本 × 基因，METABRIC 是 基因 × 样本），取并集可安全判成员。"""
    with io.open(p, encoding="utf-8") as fh:
        header = fh.readline().rstrip("\n").split("\t")[1:]
        first_col = [ln.split("\t", 1)[0] for ln in fh if ln.strip()]
    return set(header) | set(first_col)
met = rownames(os.path.join(ROOT, "data/cohorts/METABRIC_expr.tsv"))
brca = rownames(os.path.join(ROOT, "data/cohorts/BRCA_expr.tsv"))
universe = met | brca

genes = sorted(set(src["gene"]))
print("签名基因总数（唯一符号）:", len(genes))
missing = [g for g in genes if g not in universe]
print("两个矩阵里都没有的符号:", len(missing))


def hgnc_lookup(sym):
    """依次查 prev_symbol / alias_symbol。

    返回 (approved, basis, n_hits)。
    **多命中一律不映射**：不能按字母序决定基因身份。唯一命中才采用；
    多命中记为 ambiguous 并在表里保留全部候选，供人工用原始稳定 ID 判断。
    """
    for ep in ("prev_symbol", "alias_symbol"):
        for base in BASES:
            try:
                r = requests.get(f"{base}/fetch/{ep}/{sym}", headers=H, timeout=25)
                if r.status_code != 200:
                    continue
                d = r.json()
                docs = d.get("response", {}).get("docs", [])
                if docs:
                    hits = sorted({x["symbol"] for x in docs})
                    if len(hits) == 1:
                        return hits[0], ep, 1
                    return None, "ambiguous:" + ",".join(hits), len(hits)
            except Exception:
                time.sleep(1)
    return None, None, 0


rows = []
for i, g in enumerate(missing, 1):
    appr, basis, nh = hgnc_lookup(g)
    rows.append(dict(source_symbol=g, approved_symbol=appr,
                     basis=basis or "not_found", n_hits=nh))
    if i % 10 == 0:
        print(f"  {i}/{len(missing)}")
    time.sleep(0.15)

MAP = pd.DataFrame(rows)
MAP["in_metabric"] = [bool(a) and a in met for a in MAP["approved_symbol"]]
MAP["in_brca"] = [bool(a) and a in brca for a in MAP["approved_symbol"]]
MAP.to_csv(os.path.join(MAN, "signature_symbol_map.tsv"), sep="\t", index=False)

lookup = {r.source_symbol: r.approved_symbol for r in MAP.itertuples()
          if isinstance(r.approved_symbol, str) and r.approved_symbol}

out = []
for r in src.itertuples():
    used = lookup.get(r.gene, r.gene)
    out.append(dict(signature=r.signature, gene_source=r.gene, gene_used=used,
                    mapped=int(used != r.gene)))
MAPPED = pd.DataFrame(out).drop_duplicates()
MAPPED.to_csv(os.path.join(MAN, "public_signatures_genes_mapped.tsv"), sep="\t", index=False)

print("\n=== 映射结果 ===")
n_ok = int(((MAP.basis != "not_found") & MAP.approved_symbol.notna()).sum())
n_amb = int(MAP.basis.str.startswith("ambiguous").sum())
n_nf = int((MAP.basis == "not_found").sum())
print(f"  缺失符号 {len(missing)} 个：唯一命中 {n_ok}，多命中(不映射) {n_amb}，未找到 {n_nf}")
if n_amb:
    print("  多命中明细（需人工用原始稳定 ID 判断）：")
    for r in MAP[MAP.basis.str.startswith("ambiguous")].itertuples():
        print(f"    {r.source_symbol}: {r.basis}")
print(f"  映射后仍在 METABRIC 中的 {int(MAP.in_metabric.sum())} 个；在 BRCA 中 {int(MAP.in_brca.sum())} 个")
print("\n=== 各签名映射后的可用基因（相对源列表）===")
for sg, d in MAPPED.groupby("signature"):
    n_src = len(d)
    n_met = d["gene_used"].isin(met).sum()
    n_brca = d["gene_used"].isin(brca).sum()
    print(f"  {sg:20s} 源 {n_src:3d} | METABRIC {n_met:3d} | BRCA {n_brca:3d}")

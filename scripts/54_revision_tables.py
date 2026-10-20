#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
54 —— 把 52/53 的输出整理成两张紧凑补充表（不重新计算任何数值）。

输入：results/tables/{52_uno_*,53_*}.tsv
输出：results/tables/54_uno_compact.tsv
      results/tables/54_same_signature_compare.tsv
"""
import csv, io, os

ROOT = os.environ.get("SIGAUDIT_ROOT", r"D:\bioinfo\sigAudit")
TAB = os.path.join(ROOT, "results", "tables")


def rd(name):
    p = os.path.join(TAB, name)
    with io.open(p, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def f(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


# ---------- Uno 紧凑表 ----------
by = {r["cohort"]: r for r in rd("52_uno_by_cohort.tsv")}
sup = {(r["cohort"], r["seed"]): r for r in rd("52_uno_support.tsv")}
tasks = rd("52_uno_task.tsv")
summ = {r["scope"]: r for r in rd("52_uno_summary.tsv")}

OS = ["BRCA", "COAD", "KIRC", "LUAD", "MESO", "STAD"]
cols = ["d_uno_argmax_lasso", "d_uno_top5_lasso", "d_uno_argmax_ridge", "d_uno_top5_ridge",
        "d_har_argmax_lasso", "d_har_argmax_ridge", "d_har60_argmax_lasso"]

rows = []
for coh in OS + ["METABRIC"]:
    r = by[coh]
    s = [v for (c, seed), v in sup.items() if c == coh]
    n_test = f(s[0]["n_test"]) if s else None
    ev = f(s[0]["ev_test_before_tau"]) if s else None
    cens = f(s[0]["cens_test"]) if s else None
    g = f(s[0]["G_train_tau"]) if s else None
    out = {"scope": coh, "n_tasks": r["n_task"], "n_test": int(n_test) if n_test else "",
           "events_before_tau": int(ev) if ev else "", "censoring": "%.3f" % cens if cens else "",
           "G_train_60": "%.3f" % g if g else "", "tau_months": 60}
    for c in cols:
        out[c] = ("%+.4f" % f(r[c])) if f(r[c]) is not None else ""
    rows.append(out)

for scope, key in (("6 OS cohorts (macro)", "6 OS cohorts (macro)"),
                   ("METABRIC (macro)", "METABRIC (macro)")):
    r = summ[key]
    out = {"scope": scope, "n_tasks": r["n_task"], "n_test": "", "events_before_tau": "",
           "censoring": "", "G_train_60": "", "tau_months": 60}
    for c in cols:
        out[c] = ("%+.4f" % f(r[c])) if f(r[c]) is not None else ""
    rows.append(out)

fields = ["scope", "tau_months", "n_tasks", "n_test", "events_before_tau", "censoring", "G_train_60"] + cols
with io.open(os.path.join(TAB, "54_uno_compact.tsv"), "w", encoding="utf-8", newline="\n") as fh:
    w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n")
    w.writeheader()
    for r in rows:
        w.writerow(r)

# ---------- 同签名对照紧凑表 ----------
grp = rd("53_same_signature_groups.tsv")
bysig = rd("53_same_signature_by_signature.tsv")
ov = {r["signature"]: r for r in rd("53_signature_gene_overlap.tsv")}
taskall = rd("53_same_signature_tasks_all.tsv")

# 每组测试集规模
size = {}
for t in taskall:
    key = t["group"]
    size.setdefault(key, [0, 0, 0])
    size[key][0] += 1
    size[key][1] += f(t["n_test"]) or 0
    size[key][2] += f(t["ev_test"]) or 0

fields2 = ["scope", "group", "n_task", "mean_n_test", "mean_events", "n_gene_tcga", "n_gene_metabric",
           "n_shared", "Ct_lasso", "Ct_ridge", "Ct_argmax", "Ct_top5", "d_argmax_ridge",
           "d_top5_ridge", "d_top5_argmax", "Gv", "Gt", "Orel", "lasso_constant_pct"]

rows2 = []
for g in ["TCGA-BRCA / MSigDB gene sets", "TCGA-BRCA / six signatures", "METABRIC / six signatures"]:
    r = [x for x in grp if x["group"] == g][0]
    n, nt, ev = size[g]
    rows2.append({"scope": "group", "group": g, "n_task": r["n_task"],
                  "mean_n_test": "%.0f" % (nt / n), "mean_events": "%.1f" % (ev / n),
                  "n_gene_tcga": "", "n_gene_metabric": "", "n_shared": "",
                  "Ct_lasso": "%.4f" % f(r["Ct_lasso"]), "Ct_ridge": "%.4f" % f(r["Ct_ridge"]),
                  "Ct_argmax": "%.4f" % f(r["Ct_argmax"]), "Ct_top5": "%.4f" % f(r["Ct_top5"]),
                  "d_argmax_ridge": "%+.4f" % f(r["d_argmax_ridge"]),
                  "d_top5_ridge": "%+.4f" % f(r["d_top5_ridge"]),
                  "d_top5_argmax": "%+.4f" % f(r["d_top5_argmax"]),
                  "Gv": "%+.4f" % f(r["Gv"]), "Gt": "%+.4f" % f(r["Gt"]), "Orel": "%+.4f" % f(r["Orel"]),
                  "lasso_constant_pct": "%.1f" % f(r["lasso_constant_pct"])})

for r in bysig:
    o = ov.get(r["signature"], {})
    rows2.append({"scope": "signature", "group": r["group"], "n_task": r["n_task"],
                  "mean_n_test": "", "mean_events": "",
                  "n_gene_tcga": o.get("n_tcga_brca", ""), "n_gene_metabric": o.get("n_metabric", ""),
                  "n_shared": o.get("n_shared", ""),
                  "Ct_lasso": "%.4f" % f(r["Ct_lasso"]), "Ct_ridge": "%.4f" % f(r["Ct_ridge"]),
                  "Ct_argmax": "%.4f" % f(r["Ct_argmax"]), "Ct_top5": "%.4f" % f(r["Ct_top5"]),
                  "d_argmax_ridge": "%+.4f" % f(r["d_argmax_ridge"]),
                  "d_top5_ridge": "%+.4f" % f(r["d_top5_ridge"]),
                  "d_top5_argmax": "%+.4f" % f(r["d_top5_argmax"]),
                  "Gv": "%+.4f" % f(r["Gv"]), "Gt": "%+.4f" % f(r["Gt"]), "Orel": "%+.4f" % f(r["Orel"]),
                  "lasso_constant_pct": ""})

with io.open(os.path.join(TAB, "54_same_signature_compare.tsv"), "w", encoding="utf-8", newline="\n") as fh:
    w = csv.DictWriter(fh, fieldnames=fields2, delimiter="\t", lineterminator="\n")
    w.writeheader()
    for r in rows2:
        w.writerow(r)

print("54_uno_compact.tsv:", len(rows), "rows")
print("54_same_signature_compare.tsv:", len(rows2), "rows")

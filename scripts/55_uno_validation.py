#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
55（独立实现）—— 用 Python 重写同一套 Uno/IPCW 定义，与 R 版逐值比对。

刻意不复用 R 代码路径：删失 Kaplan–Meier 与加权可对比对都在这里独立实现，
只共享“书面定义”：训练折估计 G、tau = 60、严格 T_i < T_j、分数并列取 0.5、
事件时刻 G < 0.05 时不进入分母。
"""
import bisect
import csv
import io
import os

ROOT = os.environ.get("SIGAUDIT_ROOT", r"D:\bioinfo\sigAudit")
TAB = os.path.join(ROOT, "results", "tables")


def rd(name):
    with io.open(os.path.join(TAB, name), newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def km_censoring(times, events):
    """Kaplan–Meier of the censoring distribution (event = censored)."""
    ts = sorted(set(times))
    out_t, out_s = [], []
    s = 1.0
    for t in ts:
        at_risk = sum(1 for x in times if x >= t)
        d = sum(1 for x, e in zip(times, events) if x == t and e == 0)
        if d:
            s *= (1.0 - d / at_risk)
        out_t.append(t)
        out_s.append(s)
    return out_t, out_s


def g_fun(times, events):
    ts, ss = km_censoring(times, events)

    def G(x):
        i = bisect.bisect_right(ts, x) - 1
        return 1.0 if i < 0 else ss[i]
    return G


def uno(t, e, s, tau, G):
    n = len(t)
    ev = [i for i in range(n) if e[i] == 1 and t[i] < tau and G(t[i]) >= 0.05]
    num = 0.0
    den = 0.0
    for i in ev:
        w = 1.0 / G(t[i]) ** 2
        for j in range(n):
            if t[j] > t[i]:
                den += w
                if s[i] > s[j]:
                    num += w
                elif s[i] == s[j]:
                    num += 0.5 * w
    return (num / den) if den else None


rows = rd("55_uno_cases.tsv")
params = {r["dataset"]: r for r in rd("55_uno_params.tsv")}
rvals = {(r["dataset"], r["model"]): float(r["uno_r"]) for r in rd("55_uno_r_values.tsv")}

datasets = {}
for r in rows:
    d = datasets.setdefault(r["dataset"], {"train": [], "test": []})
    d[r["split"]].append(r)

models = sorted({m for m in rows[0].keys() if m not in
                 ("dataset", "split", "id", "time", "event")})
models = [m for m in models if any(r.get(m) not in (None, "", "NA") for r in rows)]

out = []
for ds, parts in datasets.items():
    tau = float(params[ds]["tau"])
    tr_t = [float(r["time"]) for r in parts["train"]]
    tr_e = [int(float(r["event"])) for r in parts["train"]]
    te_t = [float(r["time"]) for r in parts["test"]]
    te_e = [int(float(r["event"])) for r in parts["test"]]
    G = g_fun(tr_t, tr_e)
    for m in [m for m in params[ds]["models"].split(",") if m]:
        s = [float(r[m]) for r in parts["test"]]
        val = uno(te_t, te_e, s, tau, G)
        ref = rvals.get((ds, m))
        out.append({"dataset": ds, "model": m,
                    "uno_r": "" if ref is None else "%.12f" % ref,
                    "uno_python": "" if val is None else "%.12f" % val,
                    "abs_diff": "" if (val is None or ref is None) else "%.3e" % abs(val - ref)})

with io.open(os.path.join(TAB, "55_uno_validation.tsv"), "w", encoding="utf-8", newline="\n") as fh:
    w = csv.DictWriter(fh, fieldnames=["dataset", "model", "uno_r", "uno_python", "abs_diff"],
                       delimiter="\t", lineterminator="\n")
    w.writeheader()
    for r in out:
        w.writerow(r)

mx = max(float(r["abs_diff"]) for r in out if r["abs_diff"])
n_bad = sum(1 for r in out if r["abs_diff"] and float(r["abs_diff"]) > 1e-10)
print("compared:", len(out), " max abs diff:", mx, " mismatches>1e-10:", n_bad)
for r in out:
    print("  %-22s %-7s R=%s  Py=%s  d=%s" % (r["dataset"], r["model"], r["uno_r"], r["uno_python"], r["abs_diff"]))

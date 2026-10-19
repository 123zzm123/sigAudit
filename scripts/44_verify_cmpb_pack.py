#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
44 — 校验 CMPB 短稿：硬性要求、引用顺序、数字未失真、声明齐全
判据来自官方 Guide for Authors（2026-10-18 版），不自行发明规则。
"""
import io, os, re, sys, collections

ROOT = os.environ.get("SIGAUDIT_ROOT", "D:/bioinfo/sigAudit")
GEN   = os.path.join(ROOT, "docs", "paper_draft.md")          # 长稿：数字母本
SHORT = os.path.join(ROOT, "docs", "paper_cmpb_short.md")     # 短稿：唯一投稿文本
CMPB  = os.path.join(ROOT, "submission_cmpb")

fails, passes = [], []
def chk(name, cond, detail=""):
    (passes if cond else fails).append(name)
    print(f"  {'OK  ' if cond else 'FAIL'}  {name}" + (f"  {detail}" if detail and not cond else ""))

def read(p):
    with io.open(p, encoding="utf-8") as fh:
        return fh.read()

def wc(t):
    t = re.sub(r"\|[^\n]*\|", "", t)
    t = re.sub(r"[*`^]", "", t)
    return len(t.split())

gen   = read(GEN)
short = read(SHORT)
man   = read(os.path.join(CMPB, "02_manuscript_cmpb.md"))
hl    = read(os.path.join(CMPB, "03_highlights.md"))
tp    = read(os.path.join(CMPB, "01_title_page.md"))
cl    = read(os.path.join(CMPB, "04_cover_letter.md"))

i_abs  = short.index("## Abstract")
i_body = short.index("## 1. Introduction")
i_decl = short.index("## Declarations")
abstract = short[i_abs:i_body]
body     = short[i_body:i_decl]
decl     = short[i_decl:short.index("## Figure legends")]

print("=== 1. 篇幅（官方上限） ===")
chk("摘要 ≤250 词", wc(abstract) <= 250, str(wc(abstract)))
chk("正文 ≤3500 词", wc(body) <= 3500, str(wc(body)))

print("=== 2. 结构式摘要 ===")
heads = ["**Background and Objectives.**", "**Methods.**", "**Results.**", "**Conclusions.**"]
pos = [abstract.find(h) for h in heads]
chk("四段标题齐全且有序", all(p >= 0 for p in pos) and pos == sorted(pos))
chk("摘要不含引用标记", not re.search(r"\[\d", abstract))

print("=== 3. Highlights 与关键词 ===")
items = [l[2:].strip() for l in hl.splitlines() if l.startswith("- ")]
chk("Highlights 3–5 条", 3 <= len(items) <= 5, str(len(items)))
chk("每条 ≤85 字符", all(len(i) <= 85 for i in items), str([len(i) for i in items]))
kws = re.search(r"\*\*Keywords:\*\*\s*(.+)", short)
chk("关键词 3–6 个", bool(kws) and 3 <= len(kws.group(1).split(";")) <= 6,
    str(len(kws.group(1).split(";")) if kws else 0))

print("=== 4. 引用：编号顺序与覆盖 ===")
seen = []
for m in re.finditer(r"\[([0-9,\u2013\- ]+)\]", body):
    for part in m.group(1).split(","):
        part = part.strip()
        if "\u2013" in part or "-" in part:
            a, b = re.split(r"[\u2013\-]", part)
            rng = range(int(a), int(b) + 1)
        elif part.isdigit():
            rng = [int(part)]
        else:
            continue
        for n in rng:
            if n not in seen:
                seen.append(n)
chk("首次出现顺序连续", seen == list(range(1, len(seen) + 1)), str(seen))
n_ref = len(re.findall(r"^\d+\.\s", short[short.index("## References"):], re.M))
chk(f"文末 {n_ref} 条文献且与正文引用数一致",
    n_ref == max(seen) if seen else False, f"{n_ref} vs max_ref={max(seen) if seen else 0}")

print("=== 5. 数字与冻结结果表一致 ===")
def rdtab(name):
    p = os.path.join(ROOT, "results", "tables", name)
    out = {}
    if not os.path.exists(p):
        return out
    with io.open(p, encoding="utf-8") as fh:
        lines = [l.rstrip("\n").split("\t") for l in fh if l.strip()]
    hdr = lines[0]
    for row in lines[1:]:
        out[row[0]] = dict(zip(hdr, row))
    return out
def f4(x): return f"{float(x):.4f}"
ov = rdtab("47_nonconstant_overall.tsv")
new_vals = []
if "nonconstant_Lasso" in ov:
    r = ov["nonconstant_Lasso"]
    new_vals += [f4(r["Gv"]), f4(r["Gt"]), f4(r["Orel"])]
missing_new = [v for v in new_vals if v not in short]
chk("非恒定分层关键值已进入短稿（对照 47 表）", not missing_new, str(missing_new))
# 预测级入口的一致性核对必须存在且全部通过
p49 = os.path.join(ROOT, "results", "tables", "49_audit_api_check.tsv")
ok49 = False
if os.path.exists(p49):
    with io.open(p49, encoding="utf-8") as fh:
        rows = [l.rstrip("\n").split("\t") for l in fh if l.strip()][1:]
    ok49 = bool(rows) and all(r[-1].lower() == "true" for r in rows)
chk("预测级审计入口与原流水线一致性核对通过（对照 49 表）", ok49)
key = ["0.0610", "0.0237", "0.0374", "0.0377", "0.0378", "0.0064", "0.0016",
       "0.0048", "0.0084", "0.0002", "0.0039", "3,798", "3,919", "66,483",
       "1,977", "1,143", "61.2", "0.0494", "0.0109", "0.0385", "3,057",
       "0.0211", "0.0052", "72.9%"]
missing = [k for k in key if k not in short]
chk("关键结果值齐全", not missing, str(missing))
chk("Ridge 的 Orel 为 0.0377（非 0.0378）",
    "0.0374–0.0378" in short
    and "| fixed ridge | +0.0314 | −0.0064 | +0.0377 |" in short)
chk("去重百分比按任务计数（68 / 699）",
    "in 68 of 3,919 tasks (1.7%)" in short
    and "changed in 699 of 3,919 tasks (17.8%)" in short)
chk("表 2 含舍入说明",
    "Values were calculated from the full-precision results before rounding" in short)

print("=== 6. 声明齐全 ===")
for key in ["**Ethics approval", "**Competing interests", "**Funding",
            "**CRediT", "**Declaration of generative AI"]:
    chk(f"含 {key.strip('*')}", key in decl)
chk("伦理说明为不适用并给出原因", "Not applicable" in decl and "de-identified" in decl)

print("=== 7. 结构与图 ===")
chk("图 5 张、表 5 张",
    len(re.findall(r"^\*\*Fig\.\s*\d+\.", short, re.M)) == 5
    and len(re.findall(r"^\*\*Table\s*\d+\.", short, re.M)) == 5)
chk("Limitations 为 Discussion 子节", "### 4.1 Limitations" in short)
chk("含六队列 p = 0.084", "p = 0.084" in short)
chk("含搜索规模曲线一节与 Fig. 5",
    "### 3.5 Gains keep rising on validation" in short and "**Fig. 5**" in short
    and "**Fig. 5. Reference-anchored search-budget curves.**" in short)
chk("含实时软件示例与 BBC-CV 相关工作",
    "A packaged real-data example audits a stored METABRIC task" in short
    and "Bootstrapping the out-of-sample predictions" in short)

print("=== 8. 文字卫生 ===")
for name, txt in [("short", short), ("manuscript", man), ("title page", tp),
                  ("cover letter", cl), ("highlights", hl)]:
    chk(f"{name} 无中文", not re.search(r"[\u4e00-\u9fff]", txt))
chk("投稿信含四要素", all(k in cl for k in
    ["What is known", "What the study asks", "What we did", "has not been published elsewhere"]))
chk("投稿信面向 CMPB", "Computer Methods and Programs in Biomedicine" in cl)
chk("无内部文件引用", "refs_editor_note" not in man)

print("\n" + "=" * 52)
print(f"通过 {len(passes)} 项，失败 {len(fails)} 项")
for f in fails:
    print("  失败:", f)
sys.exit(1 if fails else 0)

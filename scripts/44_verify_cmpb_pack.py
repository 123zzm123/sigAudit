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
chk("首次出现顺序为 1–17", seen == list(range(1, 18)), str(seen))
n_ref = len(re.findall(r"^\d+\.\s", short[short.index("## References"):], re.M))
chk("文末 17 条文献", n_ref == 17, str(n_ref))

print("=== 5. 数字未失真（对照长稿） ===")
def dec(t):
    # 去掉补充材料编号（S1.5、S2.2 等）与标题行，避免把编号当作数值
    t = re.sub(r"\bS\d+(?:\.\d+)?\b", " ", t)
    t = re.sub(r"^#.*$", "", t, flags=re.M)
    return collections.Counter(re.findall(r"\d+\.\d+", t))
extra = set(dec(short)) - set(dec(gen))
chk("短稿未引入长稿中不存在的数值", not extra, str(sorted(extra)))
key = ["0.0610", "0.0237", "0.0374", "0.0377", "0.0378", "0.0064", "0.0016",
       "0.0048", "0.0084", "0.0002", "0.0039", "3,798", "3,919", "66,483",
       "1,977", "1,143", "61.2"]
missing = [k for k in key if k not in short]
chk("关键结果值齐全", not missing, str(missing))
chk("Ridge 的 Orel 为 0.0377（非 0.0378）",
    "0.0374, 0.0377 and 0.0378" in short
    and "| fixed ridge | +0.0314 | −0.0064 | +0.0377 |" in short)

print("=== 6. 声明齐全 ===")
for key in ["**Ethics approval", "**Competing interests", "**Funding",
            "**CRediT", "**Declaration of generative AI"]:
    chk(f"含 {key.strip('*')}", key in decl)
chk("伦理说明为不适用并给出原因", "Not applicable" in decl and "de-identified" in decl)

print("=== 7. 结构与图 ===")
chk("图 4 张、表 4 张",
    len(re.findall(r"^\*\*Fig\.\s*\d+\.", short, re.M)) == 4
    and len(re.findall(r"^\*\*Table\s*\d+\.", short, re.M)) == 4)
chk("Limitations 为 Discussion 子节", "### 4.1 Limitations" in short)
chk("含六队列 p = 0.084", "p = 0.084" in short)

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

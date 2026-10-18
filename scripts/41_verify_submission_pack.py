#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
41 — 独立校验 submission/：逐块比对生成稿与 docs/paper_draft.md
（与 40 分开实现，不复用 40 的中间变量，避免"自己验自己"）

校验项：
  1  摘要块与正文块逐字一致（正文仅允许一处补充材料引用插入）
  2  图注块逐字一致
  3  参考文献逐字一致（原稿内部工作说明已移除）
  4  Supplementary Note 1 与原 §9 正文一致
  5  图片与 results/figures 字节一致
  6  面向期刊的文件不含中文；占位符只有约定那一处
  7  结构：章节 1–8 齐全且有序、5 图注、5 表、17 条文献
  8  作者与单位在三处一致
  9  不含内部说明、不含虚构的 DOI/仓库链接
"""
import io, os, re, sys, hashlib

ROOT = os.environ.get("SIGAUDIT_ROOT", "D:/bioinfo/sigAudit")
SRC  = os.path.join(ROOT, "docs", "paper_draft.md")
SUB  = os.path.join(ROOT, "submission")

AUTHOR = "Ziming Zhao"
AFFIL  = "College of Chemistry, Fuzhou University, Fuzhou, Fujian, China"
PTR_FROM = "`Ct` is an estimate on a held-out set, not a known truth. **Fig. 1** summarises"
PTR_TO   = ("`Ct` is an estimate on a held-out set, not a known truth; the reporting conventions\n"
            "are summarised in Supplementary Note 1. **Fig. 1** summarises")

fails, passes = [], []
def chk(name, cond, detail=""):
    (passes if cond else fails).append(name)
    print(f"  {'OK  ' if cond else 'FAIL'}  {name}" + (f"  {detail}" if detail and not cond else ""))

def read(p):
    with io.open(p, encoding="utf-8") as fh:
        return fh.read()

def block(text, start, end=None):
    i = text.index(start)
    j = text.index(end, i) if end else len(text)
    return text[i:j].strip()

def strip_rule(t):
    return re.sub(r"\n+\s*---\s*$", "", t.strip()).strip()

def wc(t):
    # 与 40 同口径：只去 Markdown 标记，保留全部内容，按空白分词
    t = re.sub(r"\|[^\n]*\|", "", t)
    t = re.sub(r"[*`^]", "", t)
    return len(t.split())

for f in ["docs/paper_draft.md", "submission/02_manuscript_generic.md"]:
    if not os.path.exists(os.path.join(ROOT, f)):
        sys.exit(f"缺少 {f}")

src = read(SRC)
man = read(os.path.join(SUB, "02_manuscript_generic.md"))

# ---------- 原稿分块 ----------
s_abs = block(src, "## Abstract", "## 1. Introduction")
s_main = block(src, "## 1. Introduction", "## 7. Figure legends")
s_leg = block(src, "## 7. Figure legends", "## 8. References")
s_ref = block(src, "## 8. References", "## 9. Supplementary notes on reporting")
s_sup = block(src, "## 9. Supplementary notes on reporting")

s_ref_clean = re.sub(
    r"\*\*Note for the authors \(remove before submission\)\.\*\*.*?\n\n(?=### 8\.1)",
    "", s_ref, flags=re.S)
s_ref_clean = re.sub(r"\n> \*\*Note\.\*\*.*\Z", "", s_ref_clean, flags=re.S)
s_ref_clean = re.sub(
    r"\s*See\s+`refs_editor_note\.md` for the verified method details\.",
    "", s_ref_clean).strip()

# ---------- 生成稿分块 ----------
m_abs = strip_rule(block(man, "## Abstract", "## 1. Introduction"))
m_main = strip_rule(block(man, "## 1. Introduction", "## Declarations"))
m_leg = strip_rule(block(man, "## 7. Figure legends", "## 8. References"))
m_ref = block(man, "## 8. References")
m_decl = block(man, "## Declarations", "## 7. Figure legends")

print("=== 1. 正文块逐字比对 ===")
expected_main = s_main.replace(PTR_FROM, PTR_TO)
chk("摘要块一致", m_abs == strip_rule(s_abs))
chk("正文块一致（含唯一一处 S1 引用插入）",
    m_main == expected_main.strip(),
    f"长度 生成 {len(m_main)} / 期望 {len(expected_main.strip())}")
chk("补充材料引用只插一次", man.count("Supplementary Note 1") == 1)
chk("无连续分隔线", not re.search(r"\n---\s*\n\s*---", man))

print("=== 2-3. 图注与参考文献 ===")
chk("图注块一致", m_leg == strip_rule(s_leg))
chk("参考文献一致（内部说明已移除）", m_ref == s_ref_clean)
chk("参考文献无内部说明", "Note for the authors" not in man)
chk("无写给作者的批注", "> **Note.**" not in man and "the authors should" not in man)
chk("无内部文件引用", "refs_editor_note" not in man)

print("=== 4. 补充材料 ===")
sup_file = read(os.path.join(SUB, "05_supplementary", "S1_reporting_notes.md"))
sup_body = sup_file.split("\n", 1)[1].split("\n---")[0].strip()
s_sup_body = s_sup.split("\n", 1)[1].strip()
chk("S1 正文与原 §9 一致", sup_body == s_sup_body)
chk("S1 标题已改名", sup_file.startswith("# Supplementary Note 1"))

print("=== 5. 图片字节一致 ===")
FIGS = ["fig1_design", "fig2_gain_decomposition", "fig3_cross_split",
        "fig4_rule_pairing", "fig5_metabric"]
def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as fh:
        for c in iter(lambda: fh.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()
img_ok = True
for f in FIGS:
    for ext in ("pdf", "png"):
        a = os.path.join(ROOT, "results", "figures", f"{f}.{ext}")
        b = os.path.join(SUB, "figures", f"{f}.{ext}")
        if not (os.path.exists(b) and md5(a) == md5(b)):
            img_ok = False
            print(f"        不一致: {f}.{ext}")
chk("10 个图片文件全部字节一致", img_ok)
chk("图片数量为 10", len(os.listdir(os.path.join(SUB, "figures"))) == 10)

print("=== 6. 中文字符与占位符 ===")
for name in ["01_title_page.md", "02_manuscript_generic.md", "03_figure_legends.md",
             "04_tables_index.md", "06_cover_letter_draft.md",
             "05_supplementary/S1_reporting_notes.md"]:
    t = read(os.path.join(SUB, name))
    chk(f"{name} 无中文", not re.search(r"[一-鿿]", t))
chk("正文无 [待作者确认]", "[待作者确认]" not in man)
chk("正文占位符仅为约定项", man.count("[to be completed]") == 2,
    str(man.count("[to be completed]")))

print("=== 7. 结构 ===")
heads = [l.split()[1].rstrip(".") for l in man.splitlines() if re.match(r"^## \d+\. ", l)]
chk("章节 1–8 有序", heads == [str(i) for i in range(1, 9)], str(heads))
chk("5 条图注", m_leg.count("**Fig.") == 5)
chk("5 张表", man.count("**Table ") == 5)
nums = [int(x) for x in re.findall(r"^(\d+)\. ", m_ref, re.M)]
chk("17 条文献且编号 1–17 连续", nums[:17] == list(range(1, 18)) and len(set(nums)) == 17, str(nums))
chk("Declarations 位于 §6 与 §7 之间",
    man.index("## 6. Data and code availability") < man.index("## Declarations")
    < man.index("## 7. Figure legends"))

print("=== 8. 作者与单位 ===")
for name in ["01_title_page.md", "02_manuscript_generic.md", "06_cover_letter_draft.md"]:
    t = read(os.path.join(SUB, name))
    chk(f"{name} 含作者与单位", AUTHOR in t and AFFIL in t)

print("=== 9. 补充表与字数 ===")
sup_dir = os.path.join(SUB, "05_supplementary", "tables")
chk("补充表目录存在", os.path.isdir(sup_dir))
files = sorted(os.listdir(sup_dir)) if os.path.isdir(sup_dir) else []
chk("补充表 15 个", len(files) == 15, str(len(files)))
toc = read(os.path.join(SUB, "05_supplementary", "00_补充材料目录.md"))
rows = re.findall(r"\| `tables/([^`]+)` \|[^|]*\| ([\d,]+) \|", toc)
chk("目录列出 15 张", len(rows) == 15, str(len(rows)))
toc_ok = len(rows) == 15
for fn, n in rows:
    fp = os.path.join(sup_dir, fn)
    if not os.path.exists(fp):
        toc_ok = False; continue
    actual = sum(1 for _ in io.open(fp, encoding="utf-8")) - 1
    if actual != int(n.replace(",", "")):
        toc_ok = False
chk("补充表行数与目录一致", toc_ok)
tp = read(os.path.join(SUB, "01_title_page.md"))
m_abs_n = int(re.search(r"Abstract: (\d+) words", tp).group(1))
m_main_n = int(re.search(r"Main text[^:]*: (\d+) words", tp).group(1))
chk("标题页摘要词数与正文一致",
    m_abs_n == wc(re.sub(r"^##\s*Abstract\s*", "", s_abs)), f"{m_abs_n}")
chk("标题页正文词数与正文一致",
    m_main_n == wc(expected_main), f"{m_main_n}")
cover = read(os.path.join(SUB, "06_cover_letter_draft.md")).lower()
chk("投稿信无 DOI/仓库链接",
    "doi:" not in cover and "zenodo" not in cover and "github.com" not in cover)

# ---------- 汇总 ----------
print("\n" + "=" * 52)
print(f"通过 {len(passes)} 项，失败 {len(fails)} 项")
for f in fails:
    print("  失败:", f)
sys.exit(1 if fails else 0)

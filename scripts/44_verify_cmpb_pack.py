#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
44 — 校验 CMPB 短稿：硬性要求、引用顺序、数字未失真、声明齐全
判据来自官方 Guide for Authors（2026-10-18 版），不自行发明规则。
"""
import io, os, re, sys, glob, collections

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
_flat5 = re.sub(r"\s+", " ", short)
chk("去重百分比按任务计数（1.7% / 17.8%）",
    "changed the validation-selected winner in 1.7% of tasks" in _flat5
    and "top-five membership changed in 17.8% of tasks" in _flat5)
chk("方法写明 λ.min 与训练内调参",
    "cv.glmnet" in short and "λ.min" in short and "No penalty was tuned on validation or test data" in _flat5)
chk("主表不再列 p，改指向补充表 S31",
    "listed in Supplementary Table S31" in _flat5
    and bool(re.search(r"\| Comparison \| ΔCt vs reference \| 95% CI \|", short)))
chk("新增两段补强分析与对应表号",
    "### 3.7 Two checks of the comparison layer" in short
    and "Supplementary Table S34" in _flat5 and "Supplementary Tables S31 and S32" in _flat5)
chk("Uno 结果与同签名结果写入正文",
    "+0.0253" in short and "+0.0290" in short
    and ("-0.0093" in short or "−0.0093" in short) and "+0.0002" in short)
chk("近零不等于等效的表述",
    "are near zero and do not establish equivalence" in _flat5)
chk("常数预测表述改为全零系数",
    "coefficients were all zero" in _flat5 and "intercept-only" not in short)
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
    "A packaged example audits a stored" in short
    and "METABRIC task (18 candidates, 989 patients)" in short
    and "Bootstrapping the out-of-sample predictions" in short)
print("=== 8. 文字卫生 ===")
_flat = re.sub(r"\s+", " ", short)
_hlflat = re.sub(r"\s+", " ", hl)
chk("摘要写出 TCGA 搜索规模模式与 METABRIC 近零",
    "expanding the candidate pool increased validation gains" in _flat
    and "ridge-anchored gains remained near zero" in _flat)
chk("引言把搜索规模问题写入研究目的",
    "how it changes as the pool expands" in _flat
    and "reference-anchored search budgets" in _flat and "search budgets" in _flat)
chk("§3.5 措辞为 evaluated K greater than one",
    "at every evaluated K greater than one" in _flat)
chk("§3.5 METABRIC 段写明小而非零、签名间不均一",
    "ridge-anchored gains" in _flat and "not uniform across signatures" in _flat)
chk("Discussion 的 Ridge 结论限定在 TCGA 与 K>1",
    "in TCGA the expected held-out performance stayed below the fixed ridge start at every evaluated" in _flat)
chk("Highlights 含 TCGA 搜索规模一条",
    "In TCGA, expanding ridge-anchored search raises validation gains" in _hlflat)
chk("正文引用 Supplementary Figs. S2–S5", "Supplementary Figs. S2–S5" in _flat)
_supfigs = os.path.join(CMPB, "06_supplementary")
for _n in range(1, 6):
    chk(f"补充材料有 Fig. S{_n} 的 PDF 与 PNG",
        bool(glob.glob(os.path.join(_supfigs, f"figS{_n}_*.pdf")))
        and bool(glob.glob(os.path.join(_supfigs, f"figS{_n}_*.png"))))
_figfiles = os.path.join(CMPB, "figures")
for _n in range(1, 6):
    chk(f"投稿包 figures/ 有 Fig. {_n} 的 PDF 与 PNG",
        bool(glob.glob(os.path.join(_figfiles, f"fig{_n}_*.pdf")))
        and bool(glob.glob(os.path.join(_figfiles, f"fig{_n}_*.png"))))
chk("投稿包 figures/ 只有 5 组图文件",
    len([f for f in glob.glob(os.path.join(_figfiles, "*")) if f.endswith(".png")]) == 5)

print("=== 8. 文字卫生 ===")
for name, txt in [("short", short), ("manuscript", man), ("title page", tp),
                  ("cover letter", cl), ("highlights", hl)]:
    chk(f"{name} 无中文", not re.search(r"[\u4e00-\u9fff]", txt))
chk("投稿信含四要素", all(k in cl for k in
    ["What is known", "What the study asks", "What we did", "has not been published elsewhere"]))
chk("投稿信面向 CMPB", "Computer Methods and Programs in Biomedicine" in cl)
chk("无内部文件引用", "refs_editor_note" not in man)

print("=== 9. 补强分析的表与补充章节 ===")
_supdir = os.path.join(CMPB, "06_supplementary")
_suptab = os.path.join(_supdir, "tables")
for _f, _n in [("fixed_reference_paired.tsv", 45), ("uno_sensitivity.tsv", 9),
               ("uno_task.tsv", 2977), ("same_signature_compare.tsv", 15),
               ("same_signature_tasks.tsv", 555), ("uno_support.tsv", 60),
               ("uno_attrition.tsv", 31), ("uno_validation.tsv", 15),
               ("uno_validation_cases.tsv", 2611),
               ("same_signature_selection_change.tsv", 30),
               ("same_signature_split_identity.tsv", 5),
               ("same_signature_gene_overlap.tsv", 6)]:
    _p = os.path.join(_suptab, _f)
    _rows = 0
    if os.path.exists(_p):
        with io.open(_p, encoding="utf-8") as _fh:
            _rows = sum(1 for _l in _fh if _l.strip()) - 1
    chk(f"补充表 {_f} 存在且行数为 {_n}", _rows == _n, str(_rows))
chk("补充目录写四十二张表",
    "Forty-two tables" in io.open(os.path.join(_supdir, "00_contents.md"), encoding="utf-8").read())
_s1 = io.open(os.path.join(_supdir, "S1_methods.md"), encoding="utf-8").read()
_s2 = io.open(os.path.join(_supdir, "S2_results.md"), encoding="utf-8").read()
chk("补充方法含 S1.8 与 S1.9",
    "## S1.8 Censoring-weighted concordance" in _s1 and "## S1.9 Signature-matched comparison" in _s1)
chk("补充结果含 S2.11 与 S2.12",
    "## S2.11 Censoring-weighted concordance" in _s2 and "## S2.12 Signature-matched comparison" in _s2)
_s1f = re.sub(r"\s+", " ", _s1)
chk("补充方法写明 tau 与训练侧删失分布",
    "tau = 60 months" in _s1f and "training partition" in _s1f and "G_train(tau) >= 0.05" in _s1f)
chk("S1.8 写明按队列×划分的支持规则与层级汇总",
    "cohort-by-split level" in _s1f and "tasks, then splits within a cohort, then cohorts" in _s1f)
chk("S1.8 区分同定义核对与删失来源差异",
    "same definition" in _s1f and "not the same definition" in _s1f and "1.5e-14" in _s1f)
chk("S1.9 写明划分一致性、并列规则与验证集常数定义",
    "Split identity" in _s1f and "Tie-breaking" in _s1f and "validation predictions are" in _s1f)
_s2f = re.sub(r"\s+", " ", _s2)
chk("S2.11 写明 2,947 与事件规则未额外排除",
    "2,947 entered the analysis" in _s2f and "excluded no additional task" in _s2f)
chk("S2.11 弹性网方向为 argmax - Enet = +0.0206",
    "+0.0206" in _s2f and "elastic net itself was +0.0046 above" in _s2f)
chk("S2.11 说明宏观序保持、单队列方向改变",
    "preserved the ordering of the macro averages" in _s2f
    and "LUAD +0.0008" in _s2f and "LUAD was -0.0023" in _s2f)
chk("S2.12 基因数范围为 11 到 99 且不称癌种差异",
    "from 11 to 99 genes" in _s2f and "cancer type" not in _s2f)

print("\n" + "=" * 52)
print(f"通过 {len(passes)} 项，失败 {len(fails)} 项")
for f in fails:
    print("  失败:", f)
sys.exit(1 if fails else 0)

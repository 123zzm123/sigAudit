#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
43 — 组装 CMPB 投稿包

**短稿是唯一投稿文本**：`docs/paper_cmpb_short.md`（人工撰写并确认，脚本不再改写正文）。
本脚本只做三件事：
  1. 校验短稿是否满足 CMPB 硬性要求（摘要 ≤250、正文 ≤3500、引用按序、图注数量）
  2. 由短稿派生标题页、Highlights、投稿信、声明文件
  3. 检查声明齐全（伦理/利益/基金/CRediT/生成式 AI）

此前的版本会从长稿重新生成摘要，导致已确认的短稿被覆盖（v15 实际发生过）。
现在改为：脚本只读不写正文，正文由作者/agent 直接维护。
"""
import io, os, re, sys

ROOT = os.environ.get("SIGAUDIT_ROOT", "D:/bioinfo/sigAudit")
DOCS = os.path.join(ROOT, "docs")
OUT  = os.path.join(ROOT, "submission_cmpb")
SRC  = os.path.join(DOCS, "paper_cmpb_short.md")

AUTHOR = "Ziming Zhao"
AFFIL  = "College of Chemistry, Fuzhou University, Fuzhou, Fujian, China"

s = io.open(SRC, encoding="utf-8").read()

i_abs  = s.index("## Abstract")
i_body = s.index("## 1. Introduction")
i_decl = s.index("## Declarations")
i_leg  = s.index("## Figure legends")
i_ref  = s.index("## References")

title    = s[:i_abs].strip().splitlines()[0].lstrip("# ").strip()
abstract = s[i_abs:i_body].strip()
body     = s[i_body:i_decl].strip()
decl     = s[i_decl:i_leg].strip()
legends  = s[i_leg:i_ref].strip()
refs     = s[i_ref:].strip()

def wc(t):
    t = re.sub(r"\|[^\n]*\|", "", t)
    t = re.sub(r"[*`^]", "", t)
    return len(t.split())

W_ABS  = wc(abstract)
W_BODY = wc(body)

# ---------- 硬性要求校验（不通过则中止，不生成包） ----------
errs = []
if W_ABS > 250:
    errs.append(f"摘要 {W_ABS} 词，超过 250")
if W_BODY > 3500:
    errs.append(f"正文 {W_BODY} 词，超过 3500")

heads = ["**Background and Objectives.**", "**Methods.**", "**Results.**", "**Conclusions.**"]
pos = [abstract.find(h) for h in heads]
if not all(p >= 0 for p in pos) or pos != sorted(pos):
    errs.append("结构式摘要四段标题缺失或顺序错误")
if re.search(r"\[\d", abstract):
    errs.append("摘要中含引用标记（CMPB 要求摘要不引文献）")

# 引用按首次出现顺序编号
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
n_refs = len(re.findall(r"^\d+\.\s", refs, re.M))
if seen != list(range(1, n_refs + 1)):
    errs.append(f"引用首次出现顺序不是 1–{n_refs}：{seen}")

# 声明齐全
for key in ["**Ethics approval", "**Competing interests", "**Funding",
            "**CRediT", "**Declaration of generative AI"]:
    if key not in decl:
        errs.append(f"声明缺少 {key.strip('*')}")

if errs:
    print("短稿不满足 CMPB 要求，未生成投稿包：")
    for e in errs:
        print("  -", e)
    sys.exit(1)

n_fig = len(re.findall(r"^\*\*Fig\.\s*\d+\.", legends, re.M))
n_tab = len(re.findall(r"^\*\*Table\s*\d+\.", s, re.M))
n_ref = len(re.findall(r"^\d+\.\s", refs, re.M))

# ---------- 派生材料 ----------
HIGHLIGHTS = [
    "Paired audit of algorithm search across 8 TCGA cohorts and 3,919 tasks",
    "Validation gains over a fixed LASSO only partly replicate on held-out data",
    "The unmet-gain fraction moves from 61% to above 100% with the reference",
    "Post hoc: fixed ridge beats validation selection in all eight cohorts",
    "Search-budget curves: held-out gains flatten while validation gains keep rising",
]
for h in HIGHLIGHTS:
    assert len(h) <= 85, f"Highlight 超过 85 字符（{len(h)}）：{h}"

title_page = f"""# {title}

**Article type:** Full length article

**Authors**

{AUTHOR}^1^

**Affiliations**

^1^ {AFFIL}

**Corresponding author** (the sole author)

{AUTHOR}
College of Chemistry, Fuzhou University
Qishan Campus, Minhou County, Fuzhou, Fujian 350108, China
Email: 13317834871@163.com
Telephone: +86 13317834871

**Keywords** (6)

model selection; selection bias; transcriptomic signatures; survival analysis; benchmarking;
reproducibility

**Counts**

- Structured abstract: {W_ABS} words (limit 250)
- Main text (Introduction to Data and code availability): {W_BODY} words (limit 3,500)
- Figures: {n_fig} · Tables: {n_tab} · References: {n_ref} (limit 50)
- Highlights: {len(HIGHLIGHTS)} bullets, each ≤ 85 characters
- Supplementary: Supplementary Methods, Supplementary Results, Supplementary Figs. S1–S3, 29 S-tables

**Declarations**

Ethics: not applicable (public de-identified secondary data) · Competing interests: none ·
Funding: none · CRediT statement included · Generative-AI use declared ·
Data and code availability statement included
"""

highlights = "# Highlights\n\n" + "\n".join(f"- {h}" for h in HIGHLIGHTS) + "\n"

cover = f"""# Cover letter — Computer Methods and Programs in Biomedicine

[Date]

The Editor-in-Chief
Computer Methods and Programs in Biomedicine

Dear Editor,

We submit our manuscript **“{title}”** for consideration as a Full Length Article.

**What is known.** Validation-selected models are frequently reported as improved transcriptomic
prognostic signatures. Selection bias in performance evaluation is well documented, and nested
resampling is the usual remedy; less clear is how much of the reported gain survives in
practice, and how much the answer depends on the comparator that is used.

**What the study asks.** We used a fully paired benchmark to quantify how much of those gains
survives held-out evaluation and comparison with strong fixed models, and how much of the answer
is a property of the comparator rather than the data.

**What we did.** We evaluated 3,919 jointly evaluable tasks (66,483 candidate-model evaluations)
across eight TCGA cohorts, five stratified splits per cohort, 100 MSigDB gene sets and 18
learner configurations, and repeated the workflow in an independent microarray cohort (METABRIC,
1,977 patients) with six published signatures. The gain accounting and rule comparison are
implemented as a documented R package and accompanying scripts with assertion-based tests, so
the procedure can be reused on other omics endpoints.

**What we found.** The validation gain replicated only in part (reference-adjusted optimism
0.0374). The relative unmet fraction was reference-dependent (61% against a fixed LASSO, above
100% against a fixed ridge), whereas the absolute magnitude was similar across references. In a
post-hoc comparator analysis, a single fixed ridge regression matched or exceeded validation-set
selection in all eight cohorts, while the pre-specified top-five ensemble remained below it. The
pattern persisted in diagnostic analyses restricted to tasks where the LASSO reference was not
degenerate and after exactly duplicated candidates were removed. Post-hoc search-budget curves
showed that held-out gains flatten as the pool grows while validation gains keep rising. The
independent cohort again
showed that the apparent value of search depends on the comparator.

**Why we believe it fits the journal.** The contribution is a computational workflow and its
systematic evaluation: a paired accounting procedure, a cross-cancer benchmark, and a tested
open implementation. It gives users a concrete basis for deciding whether algorithm search adds
predictive value, and it requires no new data to be adopted.

**Declarations.** The manuscript is original, has not been published elsewhere, and is not
under consideration by another journal. The author declares no competing interests. This
research received no specific funding. The analysis uses only public, de-identified data; code
and frozen artefacts are described in the Data and code availability section.

Thank you for your consideration.

Sincerely,

{AUTHOR}
{AFFIL}
Email: 13317834871@163.com
"""

declarations_file = "# Declarations\n\n" + decl.split("## Declarations", 1)[1].strip() + "\n"

readme = f"""# CMPB 投稿包说明

**短稿是唯一投稿文本**：`docs/paper_cmpb_short.md`。
本目录由 `43_build_cmpb_pack.py` 读取该文件生成，脚本**不改写正文**。

## 当前字数

- 结构式摘要：**{W_ABS} 词**（CMPB 上限 250）
- 正文：**{W_BODY} 词**（CMPB 上限 3,500）
- 图 {n_fig} 张、表 {n_tab} 张、文献 {n_ref} 条

## 相对长稿（`docs/paper_draft.md`）的结构变化

- 13 个 Methods 小节合并为 5 个，10 个 Results 小节合并为 6 个
- 候选多样性、进入效应、六 OS 子集、模拟细节移入 Supplementary
- 跨划分稳定性图（原 Fig. 3）移为 Supplementary Fig. S1
- Limitations 作为 Discussion 的编号小节保留
- 参考文献按首次出现顺序重编为 1–{n_ref}

## 作者信息（已填）

通讯作者为作者本人：Ziming Zhao，College of Chemistry, Fuzhou University，
Qishan Campus, Minhou County, Fuzhou, Fujian 350108, China；
邮箱 13317834871@163.com，电话 +86 13317834871。

- ORCID 未填：属可选，注册后可在投稿系统补充
- 致谢：按作者要求不写
- 伦理：学校不出具免除审查证明，采用 `Not applicable` 并说明原因（官方要求如此）
- 代码与数据：正文指向已公开仓库 `https://github.com/123zzm123/sigAudit`

## 交付状态

- 双倍行距、宽页边距的投稿 Word 由 `45_make_cmpb_word.py` 生成：
  `02_manuscript.docx`（正文）、`01_title_page.docx`、`03_highlights.docx`、
  `06_supplementary/07_supplementary.docx`（S1+S2 合并）。
- 公开仓库结构与 ZIP 由 `46_export_github_repo.py` 导出；打包/导出/推送三步见根目录 `publish.sh`。
"""

os.makedirs(OUT, exist_ok=True)
files = {
    "00_README_CMPB.md": readme,
    "01_title_page.md": title_page,
    "02_manuscript_cmpb.md": s.strip() + "\n",
    "03_highlights.md": highlights,
    "04_cover_letter.md": cover,
    "05_declarations.md": declarations_file,
}
for name, text in files.items():
    with io.open(os.path.join(OUT, name), "w", encoding="utf-8", newline="\n") as fh:
        fh.write(text)

print(f"输出: {OUT}")
print(f"摘要 {W_ABS} 词；正文 {W_BODY} 词；图 {n_fig}；表 {n_tab}；文献 {n_ref}")

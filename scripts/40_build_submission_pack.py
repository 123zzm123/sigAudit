#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
40 — 由论文正文生成【通用投稿包】（目标期刊未定，只做与期刊无关的部分）

输入：docs/paper_draft.md（唯一科学内容来源，本脚本不修改它）
输出：submission/
  00_README_投稿包说明.md
  01_title_page.md
  02_manuscript_generic.md
  03_figure_legends.md
  04_tables_index.md
  05_supplementary/S1_reporting_notes.md
  05_supplementary/00_补充材料目录.md
  06_cover_letter_draft.md
  07_declarations.md

原则：
  - 不编造作者、单位、邮箱、基金号、DOI 或仓库链接；未知项一律写 [待作者确认]
  - 不改变任何数字、结论与图表内容
  - 内部工作说明（"Note for the authors…"）不进入投稿稿，移到 README
"""
import io, os, re, sys

ROOT = os.environ.get("SIGAUDIT_ROOT", "D:/bioinfo/sigAudit")
SRC  = os.path.join(ROOT, "docs", "paper_draft.md")
OUT  = os.path.join(ROOT, "submission")

AUTHOR     = "Ziming Zhao"
AFFIL      = "College of Chemistry, Fuzhou University, Fuzhou, Fujian, China"
CORR       = "[to be completed]"     # 通讯作者姓名 / 邮箱
KEYWORDS   = ("model selection; selection bias; transcriptomic signatures; "
              "survival analysis; benchmarking; reproducibility")
RUNNING    = "Baseline dependence of algorithm-search gains"

def w(text):
    path = os.path.join(OUT, text[0])
    os.makedirs(os.path.dirname(path) or OUT, exist_ok=True)
    with io.open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(text[1])
    return path

s = io.open(SRC, encoding="utf-8").read()

# ---------- 切分 ----------
i_abs  = s.index("## Abstract")
i_intro = s.index("## 1. Introduction")
i_leg   = s.index("## 7. Figure legends")
i_ref   = s.index("## 8. References")
i_sup   = s.index("## 9. Supplementary notes on reporting")

title   = s[:i_abs].strip().splitlines()[0].lstrip("# ").strip()
abstract = s[i_abs:i_intro].strip()
main     = s[i_intro:i_leg].strip()
legends  = s[i_leg:i_ref].strip()
refs     = s[i_ref:i_sup].strip()
supnote  = s[i_sup:].strip()

# ---------- 在正文中给补充说明加一处引用 ----------
PTR_FROM = "`Ct` is an estimate on a held-out set, not a known truth. **Fig. 1** summarises"
PTR_TO   = ("`Ct` is an estimate on a held-out set, not a known truth; the reporting conventions\n"
            "are summarised in Supplementary Note 1. **Fig. 1** summarises")
assert main.count(PTR_FROM) == 1, "补充说明引用锚点未唯一命中"
main = main.replace(PTR_FROM, PTR_TO)

# ---------- 删除参考文献前的内部工作说明 ----------
NOTE_RE = re.compile(
    r"\*\*Note for the authors \(remove before submission\)\.\*\*.*?\n\n(?=### 8\.1)", re.S)
refs_clean, n_note = NOTE_RE.subn("", refs)
assert n_note == 1, "未找到内部工作说明段落，请检查正文格式"

# 参考文献末尾还有一段写给作者的批注（> **Note.** …），同样不得随稿提交
NOTE2_RE = re.compile(r"\n> \*\*Note\.\*\*.*\Z", re.S)
refs_clean, n_note2 = NOTE2_RE.subn("", refs_clean)
assert n_note2 == 1, "未找到文献末尾的作者批注，请检查正文格式"

# 文献注释里指向内部工作文件的句子，也必须在投稿版中删除
refs_clean, n_refnote = re.subn(
    r"\s*See\s+`refs_editor_note\.md` for the verified method details\.", "", refs_clean)
assert n_refnote == 1, "未找到指向内部工作文件的句子"
refs_clean = refs_clean.strip()

# ---------- 字数 ----------
# 只去掉 Markdown 标记，保留全部实际内容；按空白分词，便于与投稿系统口径对照
def wc(t):
    t = re.sub(r"\|[^\n]*\|", "", t)     # 去表格行
    t = re.sub(r"[*`^]", "", t)          # 只去标记，不去内容
    return len(t.split())

abs_body = re.sub(r"^##\s*Abstract\s*", "", abstract)
W_ABS = wc(abs_body)
W_MAIN = wc(main)
W_ALL = W_ABS + W_MAIN
assert 300 < W_ABS < 500, W_ABS
assert 6000 < W_MAIN < 11000, W_MAIN

# ---------- 01 标题页 ----------
title_page = f"""# {title}

**Authors**

{AUTHOR}^1^

**Affiliations**

^1^ {AFFIL}

**Corresponding author**

{CORR}

**Running title**

{RUNNING}

**Keywords**

{KEYWORDS}

**Article type**

Original research article

**Counts** (whitespace-delimited, after removing Markdown markers only; recount in the
submission system before final upload)

- Abstract: {W_ABS} words
- Main text (Introduction to Data availability): {W_MAIN} words
- Figures: 5
- Tables: 5
- References: 17
- Supplementary: 1 note + 15 tables
"""

# ---------- 07 声明（正文用干净版；工作文件另附中文待办） ----------
declarations = """## Declarations

**Funding.** This research received no specific grant from any funding agency in the public,
commercial, or not-for-profit sectors.

**Competing interests.** The authors declare that they have no competing interests.

**Ethics approval and consent to participate.** This study is a secondary analysis of
publicly available, de-identified data (TCGA Pan-Cancer Atlas 2018 through cBioPortal, and
METABRIC through cBioPortal). No new human participants or animals were involved, and no
institutional ethics approval or informed consent was required.

**Consent for publication.** Not applicable.

**Author contributions.** Z.Z. conceived and designed the study, performed the analyses and
wrote the manuscript. [Further contributions to be completed once the author list is
finalised.]

**Acknowledgements.** [to be completed]
"""

declarations_work = declarations.replace("## Declarations", "# Declarations", 1) + """
---

## 内部待确认（本段不要随稿提交）

1. **伦理声明**：现写法是“公共去标识数据的二次分析无需伦理审批”。请按福州大学的实际规定
   确认；若学校要求，应改为“经 XX 委员会审查，本研究豁免伦理审批”。
2. **其他作者**：若最终只有作者本人一位，Author contributions 保持单数即可；
   若有合作者，按实际分工逐项补充。
3. **致谢**：是否需要感谢导师或课题组，由作者决定；不需要可整段删除。
4. **基金**：现为“无特定资助”。若后续获得资助，替换为基金名称与编号。
5. **单位**：目前只写到 `Fuzhou University`，多数期刊要求写到学院/系。
"""

# ---------- 02 投稿正文 ----------
def cut_rule(t):
    """去掉块尾的分隔线，避免拼接后出现连续两条 ---"""
    return re.sub(r"\n+\s*---\s*$", "", t.strip()).strip()

manuscript = f"""# {title}

**Authors:** {AUTHOR}^1^

**Affiliations:** ^1^ {AFFIL}

**Corresponding author:** {CORR}

**Keywords:** {KEYWORDS}

---

{cut_rule(abstract)}

---

{cut_rule(main)}

---

{declarations}

---

{cut_rule(legends)}

---

{cut_rule(refs_clean)}
"""

# ---------- 03 图注 ----------
leg = legends.split("\n", 1)[1].strip()
assert leg.count("**Fig.") == 5, leg.count("**Fig.")
legends_file = ("# Figure legends\n\n"
                "Figures are provided as separate files (`fig1_design`, `fig2_gain_decomposition`,"
                " `fig3_cross_split`, `fig4_rule_pairing`, `fig5_metabric`; PDF and 300-dpi PNG).\n\n"
                + leg + "\n")

# ---------- 04 表格索引 ----------
tables_file = """# Tables

All five tables are typeset inline in the manuscript, immediately after the paragraph that
first refers to them. They are listed here so that a journal-specific layout (tables at the
end of the file, or separate table files) can be produced without searching the text.

| Table | Location in manuscript | Content |
|---|---|---|
| Table 1 | Methods §2.1 | Cohorts and endpoints (eight TCGA cohorts; n, events, endpoint) |
| Table 2 | Results §3.4 | Held-out test C-index relative to each fixed reference: `argmax_val − ridge`, `top5 − ridge`, `LASSO − ridge` |
| Table 3 | Results §3.4 | The same gain accounting under three fixed references (LASSO, ridge, elastic net α = 0.5) |
| Table 4 | Results §3.9 | Simulated `Orel` by condition (Monte Carlo standard error in brackets; n = 100 replicates) |
| Table 5 | Results §3.10 | Effect sizes in the TCGA benchmark and in the independent METABRIC check |

The comparison of the fifteen decision rules against `argmax_val` is presented in **Fig. 4**
and in the rule-performance supplement (S3); Table 2 reports differences against the fixed
references, not that rule comparison.

*When the target journal requires tables at the end of the manuscript, move the five blocks
verbatim; no table content needs to change.*
"""

# ---------- 05 补充材料 ----------
import shutil

sup_note = supnote.replace(
    "## 9. Supplementary notes on reporting",
    "# Supplementary Note 1. Terminology and reporting conventions")

# 补充表：使用当前分析口径的文件，构建时实际拷入并统计行数
SUP_TABLES = [
    ("S2",  "tcga_rule_summary.tsv",            "results/tables/20_summary.tsv",
     "Rule-level summary over all evaluable tasks", "§3.1–3.4"),
    ("S3",  "tcga_rule_performance_long.tsv",   "results/tables/20_rule_performance_long.tsv",
     "Task-level record for every decision rule", "§3.3, Fig. 4"),
    ("S4",  "model_metrics_long.tsv",           "results/tables/20_model_metrics_long.tsv",
     "Candidate-model evaluations (validation and test C-index)", "§3.2–3.4"),
    ("S5",  "baseline_sensitivity.tsv",         "results/tables/29_baseline_sensitivity.tsv",
     "Gain accounting under three fixed references", "§3.4, Table 3"),
    ("S6",  "candidate_diversity.tsv",          "results/manifest/candidate_diversity.tsv",
     "Per-task candidate-pool diversity", "§3.5"),
    ("S7",  "candidate_top5_dup.tsv",           "results/manifest/candidate_top5_dup.tsv",
     "Duplicate candidate pairs entering the top-5", "§3.5"),
    ("S8",  "entry_effect.tsv",                 "results/tables/21_entry_effect.tsv",
     "Entry-effect estimates", "§3.6"),
    ("S9",  "entry_summary.tsv",                "results/tables/21_entry_summary.tsv",
     "Entry-effect summary", "§3.6"),
    ("S10", "signature_tasks.tsv",              "results/tables/26_signature_tasks.tsv",
     "Published-signature task metrics", "§3.8"),
    ("S11", "signature_summary.tsv",            "results/tables/26_signature_summary.tsv",
     "Published-signature summary", "§3.8"),
    ("S12", "simulation_replicates.tsv",        "results/tables/25_simulation_final_replicates.tsv",
     "Simulation replicates behind Table 4", "§3.9"),
    ("S13", "metabric_model_metrics.tsv",       "results/tables/31_metabric_model_metrics.tsv",
     "METABRIC candidate-model metrics", "§3.10"),
    ("S14", "metabric_rule_performance.tsv",    "results/tables/31_metabric_rule_performance.tsv",
     "METABRIC decision-rule metrics", "§3.10, Table 5"),
    ("S15", "metabric_baseline_sensitivity.tsv", "results/tables/32_metabric_baseline_sensitivity.tsv",
     "METABRIC gain accounting under three references", "§3.10"),
    ("S16", "metabric_prediction_check.tsv",    "results/tables/36_metabric_prediction_check.tsv",
     "Independent reconstruction of METABRIC C-indices", "§6"),
]

SUP_DIR = os.path.join(OUT, "05_supplementary", "tables")
os.makedirs(SUP_DIR, exist_ok=True)
sup_rows = []
for sid, fname, rel, desc, where in SUP_TABLES:
    src = os.path.join(ROOT, rel.replace("/", os.sep))
    if not os.path.exists(src):
        sys.exit(f"缺少补充表源文件: {src}")
    shutil.copy2(src, os.path.join(SUP_DIR, fname))
    with io.open(src, encoding="utf-8") as fh:
        n = sum(1 for _ in fh) - 1
    if n <= 0:
        sys.exit(f"补充表 {rel} 行数为 {n}，请检查")
    sup_rows.append(f"| Supplementary Table {sid[1:]} | `tables/{fname}` | {desc} | {n:,} | {where} |")

sup_toc = ("# Supplementary material — contents\n\n"
           "All supplementary tables are the current-analysis files, copied unchanged from the\n"
           "frozen results directory. Temporary numbering will be remapped if the target journal\n"
           "requires a different scheme.\n\n"
           "| Item | File | Content | Rows | Supports |\n"
           "|---|---|---|---|---|\n"
           "| Supplementary Note 1 | `S1_reporting_notes.md` | Terminology (`Orel` vs `Oraw`), "
           "status of `Ct`, and the exploratory nature of cohort-level p-values | — | §2.6 |\n"
           + "\n".join(sup_rows) + "\n")

# ---------- 06 cover letter ----------
cover = f"""# Cover letter (draft)

[Date]

The Editor-in-Chief
[Target journal]

Dear Editor,

We submit our manuscript **“{title}”** for consideration as an Original Research article.

**What the study asks.** Validation-selected models are often presented as improved
transcriptomic prognostic signatures, yet apparent gains can reflect both selection on a noisy
criterion and comparison with an undemanding baseline. We used a fully paired benchmark to
quantify how much of those gains survives held-out evaluation and comparison with strong fixed
models, and how much of the answer is a property of the comparator rather than the data.

**What we did.** Using eight TCGA cohorts (n = 3,798), five stratified splits per cohort, 100
sampled MSigDB gene sets and 18 learner configurations, we evaluated 3,919 jointly evaluable
tasks (66,483 candidate-model evaluations) with a fully paired design, so that every selection rule is compared
on exactly the same tasks. We separated the validation gain (`Gv`), the held-out gain (`Gt`)
and their difference (`Orel`), and repeated the workflow in an independent microarray cohort
(METABRIC, 1,977 patients, 1,143 deaths) with six published signatures.

**What we found.** The validation gain replicated only in part (`Orel` = 0.0374). The widely
quoted *relative* unmet fraction was entirely reference-dependent (61% against a fixed LASSO,
above 100% against a fixed ridge), whereas the *absolute* magnitude was similar across
references. A single fixed ridge regression matched or exceeded validation-set selection in all
eight cohorts, and the top-five ensemble — the best implementable rule among those evaluated —
remained below the fixed ridge in this benchmark. The independent cohort reproduced the sign
and the reference dependence of the effect, while showing that the margin over a dense
comparator is cohort-dependent.

**Why we believe it fits.** The contribution is methodological and empirical at once: a paired
accounting framework, a cross-cancer benchmark, and a concrete reporting recommendation. It
speaks directly to authors, reviewers and readers who use multi-algorithm pipelines, and it
does not require new data to be adopted.

**Declarations.** The manuscript is original and is not under consideration elsewhere. The
authors declare no competing interests. This research received no specific funding. The
analysis uses only public, de-identified data; the analysis code and frozen artefacts are
described in the Data and code availability section of the manuscript.

Thank you for your consideration.

Sincerely,

{AUTHOR}
{AFFIL}
{CORR}
"""

# ---------- 00 README ----------
readme = f"""# 通用投稿包说明

生成时间：{__import__('datetime').datetime.now().strftime('%Y-%m-%d %H:%M')}
来源：`docs/paper_draft.md`（科学内容唯一来源，本包不修改它）

## 包里有什么

| 文件 | 用途 | 状态 |
|---|---|---|
| `01_title_page.md` | 标题、作者、单位、通讯、关键词、字数 | 作者与单位已填；通讯待补 |
| `02_manuscript_generic.md` | 完整投稿正文（含声明、参考文献、图注） | 通用稿；待通讯、伦理与目标刊格式确认 |
| `03_figure_legends.md` | 五张图的统一图注 | 可用 |
| `04_tables_index.md` | 五张表的位置索引 | 可用 |
| `05_supplementary/` | 补充说明 S1 **与 15 张已装入的补充表** | 已装入当前分析文件 |
| `06_cover_letter_draft.md` | 投稿信草稿 | 需填目标刊与日期 |
| `07_declarations.md` | 基金、利益冲突、伦理、贡献声明 | 英文部分可用；中文待办另存 |

**这不是"可直接投"的终稿。** 通讯作者、伦理表述、目标期刊格式三项未定，
标题页与声明中相应位置标为 `[to be completed]`。

## 基金与利益冲突的处理

**基金**：按作者本人说明"无资助"，采用标准事实句——*This research received no specific
grant from any funding agency in the public, commercial, or not-for-profit sectors.*
留空会触发编辑部质询；这句话陈述的是事实，不是免责。

**利益冲突**：指可能影响研究判断的经济或个人关系（任职、咨询、专利、持股、资助、
亲属关系等）。作者说明没有此类关系，采用标准句——*The authors declare that they have no
competing interests.*

## 还需要作者确认（共 3 项）

1. **通讯作者**——姓名与邮箱（`01_title_page.md`、`06_cover_letter_draft.md` 内为
   `[to be completed]`）
2. **伦理声明**——现写法是"公共去标识数据的二次分析，无需伦理审批"，
   请按福州大学的实际规定确认后再定稿
3. **其他作者**——若最终只有一位，Author contributions 保持单数即可

## 字数怎么算的

摘要 **{W_ABS} 词**、正文 **{W_MAIN} 词**，口径为：去掉 Markdown 标记后按空白分词，
**保留全部实际内容**（含行内指标名与数值）。

不同工具和不同计数方式会差 10–25 词，因此这个数是**编辑工作目标**，不是投稿系统读数。
目标刊确定后，把摘要放进投稿系统按其规则重新计数；压缩时优先删方法细节与重复数列，
结论与数字不动。

## 目标期刊定了之后要改什么（都不涉及内容）

- 摘要是否需压到该刊上限（以该刊当前指南为准，不要按传言数字执行）
- 参考文献格式（现为姓名-年份式）
- 章节编号（正文有 `§2.1` 交叉引用，改编号需同步）
- 表格位置（现内嵌于正文，多数期刊要求移到文末）
- 补充材料编号与文件格式（现用临时编号 S1–S16）

## 投稿信的使用注意

- 发送前替换 `[Target journal]` 与日期
- 部分期刊要求 100 词以内新颖性声明，可用"**What we found**"段
- 存档未公开前，不得在信中写"已公开可获取"；现文只写"described in the Data and code
  availability section"

## 本包不包含

- 作者签名、投稿系统账号、APC 支付信息
- 任何仓库链接或 DOI——**在存档真正公开前不得填写**
- 目标期刊的排版模板
- 完整 40 份主分析逐患者预测（只附 MESO 一份样例）
"""

for name, text in [
    ("00_README_投稿包说明.md", readme),
    ("01_title_page.md", title_page),
    ("02_manuscript_generic.md", manuscript),
    ("03_figure_legends.md", legends_file),
    ("04_tables_index.md", tables_file),
    ("05_supplementary/S1_reporting_notes.md", sup_note),
    ("05_supplementary/00_补充材料目录.md", sup_toc),
    ("06_cover_letter_draft.md", cover),
    ("07_declarations.md", declarations_work),
]:
    w((name, text))

# ---------- 图片随包 ----------
FIG_SRC = os.path.join(ROOT, "results", "figures")
FIG_OUT = os.path.join(OUT, "figures")
os.makedirs(FIG_OUT, exist_ok=True)
FIGS = ["fig1_design", "fig2_gain_decomposition", "fig3_cross_split",
        "fig4_rule_pairing", "fig5_metabric"]
for f in FIGS:
    for ext in ("pdf", "png"):
        src = os.path.join(FIG_SRC, f + "." + ext)
        if not os.path.exists(src):
            sys.exit(f"缺少图片文件: {src}")
        shutil.copy2(src, os.path.join(FIG_OUT, f + "." + ext))

# ---------- 自检 ----------
checks = {
    "标题页含作者": AUTHOR in title_page,
    "标题页含单位": AFFIL in title_page,
    "正文无内部说明": "Note for the authors" not in manuscript,
    "正文无作者批注": "> **Note.**" not in manuscript and "the authors should" not in manuscript,
    "正文无内部文件引用": "refs_editor_note" not in manuscript,
    "正文含 5 条图注": all(f"**Fig. {i}." in manuscript for i in range(1, 6))
                      and manuscript.count("**Fig. 5.") == 1,
    "正文含 5 张表": manuscript.count("**Table ") == 5,
    "正文含声明": "## Declarations" in manuscript,
    "正文含 S1 引用外移": "Supplementary notes on reporting" not in manuscript,
    "图注 5 条": leg.count("**Fig.") == 5,
    "补充说明已改名": sup_note.startswith("# Supplementary Note 1"),
    "cover 无 DOI": "doi:" not in cover.lower() and "zenodo" not in cover.lower(),
    "图片 10 个": len(os.listdir(FIG_OUT)) == 10,
    "正文无中文": not re.search(r"[一-鿿]", manuscript),
    "正文引用 S1": "Supplementary Note 1" in manuscript,
    "章节顺序 1-8": [l.split()[1].rstrip(".") for l in manuscript.splitlines()
                      if re.match(r"^## \d+\. ", l)] == [str(i) for i in range(1, 9)],
    "S1 无编辑说明": "This note was Section 9" not in sup_note,
    "标题页无内部脚注": "journal-neutral" not in title_page,
    "投稿信无作者备注": "Notes for the author" not in cover,
    "补充表 15 张": len(sup_rows) == 15,
}
bad = [k for k, v in checks.items() if not v]
print(f"输出目录: {OUT}")
print(f"摘要 {W_ABS} 词；正文 {W_MAIN} 词；合计 {W_ALL} 词")
for k, v in checks.items():
    print(f"  {'OK ' if v else 'FAIL'}  {k}")
sys.exit(1 if bad else 0)

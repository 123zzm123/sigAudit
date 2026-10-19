#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
45 — 生成 CMPB 投稿用 Word 文件
官方要求：double spacing, wide margins, 可编辑；标题页单独一页。

输出：
  submission_cmpb/01_title_page.docx
  submission_cmpb/02_manuscript.docx
"""
import io, os, re, sys
from docx import Document
from docx.shared import Pt, Cm
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml.ns import qn
from docx.oxml import OxmlElement

ROOT = os.environ.get("SIGAUDIT_ROOT", "D:/bioinfo/sigAudit")
DOCS = os.path.join(ROOT, "docs")
OUT  = os.path.join(ROOT, "submission_cmpb")
SRC  = os.path.join(DOCS, "paper_cmpb_short.md")
TITLE_PAGE = os.path.join(OUT, "01_title_page.md")

TOKEN = re.compile(r"(\*\*.+?\*\*|\*[^*]+?\*|`[^`]+?`|\^[^^]+?\^)")

def add_inline(par, text, size=12, italic=False, bold=False):
    for tok in TOKEN.split(text):
        if not tok:
            continue
        b, i, sup = bold, italic, False
        val = tok
        if tok.startswith("**") and tok.endswith("**") and len(tok) > 4:
            b, val = True, tok[2:-2]
        elif tok.startswith("*") and tok.endswith("*") and len(tok) > 2:
            i, val = True, tok[1:-1]
        elif tok.startswith("`") and tok.endswith("`") and len(tok) > 2:
            i, val = True, tok[1:-1]
        elif tok.startswith("^") and tok.endswith("^") and len(tok) > 2:
            sup, val = True, tok[1:-1]
        r = par.add_run(val)
        r.bold, r.italic, r.font.size = b, i, Pt(size)
        r.font.superscript = sup

def new_doc(double=True):
    doc = Document()
    sec = doc.sections[0]
    sec.page_width, sec.page_height = Cm(21.0), Cm(29.7)
    sec.left_margin = sec.right_margin = Cm(2.5)
    sec.top_margin = sec.bottom_margin = Cm(2.5)
    st = doc.styles["Normal"]
    st.font.name = "Times New Roman"
    st.font.size = Pt(12)
    st.element.rPr.rFonts.set(qn("w:eastAsia"), "Times New Roman")
    st.paragraph_format.space_after = Pt(0)
    st.paragraph_format.line_spacing = 2.0 if double else 1.15
    p = sec.footer.paragraphs[0]
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    fld = OxmlElement("w:fldSimple")
    fld.set(qn("w:instr"), "PAGE")
    p._p.append(fld)
    return doc

def para(doc, text="", align="justify", size=12, bold=False, italic=False,
         space_after=0, hang=False, keep=False):
    p = doc.add_paragraph()
    p.alignment = {"justify": WD_ALIGN_PARAGRAPH.JUSTIFY,
                   "left": WD_ALIGN_PARAGRAPH.LEFT,
                   "center": WD_ALIGN_PARAGRAPH.CENTER}[align]
    p.paragraph_format.space_after = Pt(space_after)
    if keep:
        p.paragraph_format.keep_with_next = True
    if hang:
        p.paragraph_format.left_indent = Cm(0.75)
        p.paragraph_format.first_line_indent = Cm(-0.75)
        p.paragraph_format.line_spacing = 1.0
    if text:
        add_inline(p, text, size=size, italic=italic, bold=bold)
    return p

def split_row(line):
    line = line.strip().strip("|")
    return [c.strip() for c in line.split("|")]

def is_sep(line):
    return bool(re.fullmatch(r"\|?[\s:\-|]+\|?", line.strip())) and "-" in line

def _trPr(row):
    return row._tr.get_or_add_trPr()

def set_cant_split(row):
    # 禁止同一行跨页断开
    el = OxmlElement("w:cantSplit"); el.set(qn("w:val"), "true"); _trPr(row).append(el)

def set_repeat_header(row):
    # 表头跨页时重复
    el = OxmlElement("w:tblHeader"); el.set(qn("w:val"), "true"); _trPr(row).append(el)

def add_table(doc, rows, tbl_no=0):
    ncol = max(len(r) for r in rows)
    t = doc.add_table(rows=0, cols=ncol)
    t.style = "Table Grid"
    for i, row in enumerate(rows):
        cells = t.add_row().cells
        for j in range(ncol):
            p = cells[j].paragraphs[0]
            p.paragraph_format.line_spacing = 1.0
            p.paragraph_format.space_after = Pt(2)
            add_inline(p, row[j] if j < len(row) else "", size=9, bold=(i == 0))
    # Table 1 的 study ID 列需要更宽，避免长标识符折行
    if tbl_no == 1 and ncol == 5:
        widths = [Cm(3.2), Cm(6.2), Cm(2.2), Cm(1.9), Cm(1.6)]
        for row in t.rows:
            for j, w in enumerate(widths):
                row.cells[j].width = w
    for i, row in enumerate(t.rows):
        set_cant_split(row)
        if i == 0:
            set_repeat_header(row)
    # 短表整体不拆页：除最后一行外，所有单元格段落设 keep_with_next
    if len(t.rows) <= 6:
        for row in t.rows[:-1]:
            for cell in row.cells:
                for p in cell.paragraphs:
                    p.paragraph_format.keep_with_next = True
    return t

# ---------- 标题页 ----------
tp_md = io.open(TITLE_PAGE, encoding="utf-8").read().splitlines()
doc = new_doc(double=False)
first = True
for line in tp_md:
    s = line.strip()
    if not s or s == "---":
        continue
    if s.startswith("# "):
        para(doc, s[2:], align="center", size=16, bold=True, space_after=18)
    elif s.startswith("**") and s.endswith("**") and ":" not in s:
        para(doc, s.strip("*"), align="left", size=12, bold=True, space_after=4)
    else:
        para(doc, s, align="left", size=12, space_after=4)
doc.save(os.path.join(OUT, "01_title_page.docx"))

# ---------- 通用 Markdown -> Word 渲染 ----------
def render_md(doc, md):
    i = 0
    n_tab = 0
    while i < len(md):
        s = md[i].strip()
        if not s or s == "---":
            i += 1
            continue
        if s.startswith("# ") and not s.startswith("## "):
            para(doc, s[2:], align="center", size=14, bold=True, space_after=12)
            i += 1
            continue
        if s.startswith("## "):
            para(doc, s[3:].strip(), align="left", size=12, bold=True, space_after=6, keep=True)
            i += 1
            continue
        if s.startswith("### "):
            para(doc, s[4:].strip(), align="left", size=12, bold=True, space_after=4, keep=True)
            i += 1
            continue
        if s.startswith("```"):
            code = []
            i += 1
            while i < len(md) and not md[i].strip().startswith("```"):
                code.append(md[i]); i += 1
            i += 1
            for line in code:
                p = para(doc, line, align="left", size=10, space_after=0)
                for r in p.runs:
                    r.font.name = "Consolas"
            continue
        if s.startswith("|"):
            rows = []
            while i < len(md) and md[i].strip().startswith("|"):
                if not is_sep(md[i]):
                    rows.append(split_row(md[i]))
                i += 1
            add_table(doc, rows, tbl_no=n_tab + 1)
            n_tab += 1
            para(doc, "", space_after=0)
            continue
        if re.match(r"^\d+\.\s", s):
            buf = [s]
            i += 1
            while i < len(md):
                nxt = md[i].strip()
                if not nxt or nxt.startswith("#") or re.match(r"^\d+\.\s", nxt):
                    break
                buf.append(nxt)
                i += 1
            para(doc, " ".join(buf), align="left", size=10, space_after=4, hang=True)
            doc.paragraphs[-1].paragraph_format.line_spacing = 2.0
            continue
        buf = [s]
        i += 1
        while i < len(md):
            nxt = md[i].strip()
            if (not nxt or nxt == "---" or nxt.startswith("#") or nxt.startswith("|")
                    or nxt.startswith("```") or re.match(r"^\d+\.\s", nxt)):
                break
            buf.append(nxt)
            i += 1
        keepit = bool(re.match(r"^\*\*(Table|Fig\.)\s", s))
        para(doc, " ".join(buf), align="justify", size=12, space_after=6, keep=keepit)
    return n_tab

# ---------- 正文 ----------
md = io.open(SRC, encoding="utf-8").read().splitlines()
doc = new_doc(double=True)
n_tab = render_md(doc, md)
doc.save(os.path.join(OUT, "02_manuscript.docx"))

# ---------- Highlights ----------
doc = new_doc(double=False)
render_md(doc, io.open(os.path.join(OUT, "03_highlights.md"), encoding="utf-8").read().splitlines())
doc.save(os.path.join(OUT, "03_highlights.docx"))

# ---------- 补充材料（S1 + S2 合成一个可上传 Word） ----------
sup_lines = []
for p in [os.path.join(OUT, "06_supplementary", "S1_methods.md"),
          os.path.join(OUT, "06_supplementary", "S2_results.md")]:
    if os.path.exists(p):
        sup_lines += io.open(p, encoding="utf-8").read().splitlines() + [""]
doc = new_doc(double=False)
render_md(doc, sup_lines)
doc.save(os.path.join(OUT, "06_supplementary", "07_supplementary.docx"))

print("生成:")
for f in ["01_title_page.docx", "02_manuscript.docx", "03_highlights.docx",
          os.path.join("06_supplementary", "07_supplementary.docx")]:
    p = os.path.join(OUT, f)
    print(f"  {f}  {os.path.getsize(p)/1024:.1f} KB")
print(f"  正文表格 {n_tab} 张")

# 自检：双倍行距、页边距
d = Document(os.path.join(OUT, "02_manuscript.docx"))
ls = d.styles["Normal"].paragraph_format.line_spacing
sec = d.sections[0]
print(f"  行距 {ls}（应为 2.0）；页边距 "
      f"{sec.left_margin.cm:.1f}/{sec.right_margin.cm:.1f} cm")
if abs(ls - 2.0) > 1e-6:
    sys.exit("行距不是双倍")

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
42 — 把投稿正文排成可阅读的 docx/PDF（仅供阅读与检查，不是期刊排版稿）

输入：submission/02_manuscript_generic.md + submission/figures/*.png
输出：submission/paper_reading_copy.docx
      <桌面>/sigAudit_论文_阅读版.pdf   （由 Word 导出）

排版：Times New Roman 11pt，A4，页边距 2.5cm，正文两端对齐，页脚页码；
      变量用斜体（不用等宽），图插在图注上方，参考文献左对齐 + 悬挂缩进。
"""
import io, os, re, sys, subprocess
from docx import Document
from docx.shared import Pt, Cm, Inches
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml.ns import qn
from docx.oxml import OxmlElement

ROOT = os.environ.get("SIGAUDIT_ROOT", "D:/bioinfo/sigAudit")
SUB  = os.path.join(ROOT, "submission")
SRC  = sys.argv[1] if len(sys.argv) > 1 else os.path.join(SUB, "02_manuscript_generic.md")
DOCX = os.path.join(SUB, "paper_reading_copy.docx")
PDF  = os.path.join(os.path.expanduser("~"), "Desktop",
                    "sigAudit_论文_CMPB版.pdf" if "cmpb" in SRC.lower()
                    else "sigAudit_论文_阅读版.pdf")

# 长稿：Fig.3 是跨划分稳定性；短稿：Fig.3 是规则比较（跨划分图移为 Fig. S1）
if "cmpb_short" in SRC.lower():
    FIGMAP = {1: "fig1_design.png", 2: "fig2_gain_decomposition.png",
              3: "fig4_rule_pairing.png", 4: "fig5_metabric.png"}
else:
    FIGMAP = {1: "fig1_design.png", 2: "fig2_gain_decomposition.png",
              3: "fig3_cross_split.png", 4: "fig4_rule_pairing.png",
              5: "fig5_metabric.png"}

# ---------------- 页面与样式 ----------------
doc = Document()
sec = doc.sections[0]
sec.page_width, sec.page_height = Cm(21.0), Cm(29.7)
sec.left_margin = sec.right_margin = Cm(2.5)
sec.top_margin = sec.bottom_margin = Cm(2.4)

st = doc.styles["Normal"]
st.font.name = "Times New Roman"
st.font.size = Pt(11)
st.element.rPr.rFonts.set(qn("w:eastAsia"), "Times New Roman")
st.paragraph_format.space_after = Pt(6)
st.paragraph_format.line_spacing = 1.15

def add_page_number(paragraph):
    fld = OxmlElement("w:fldSimple")
    fld.set(qn("w:instr"), "PAGE")
    paragraph._p.append(fld)

footer_p = sec.footer.paragraphs[0]
footer_p.alignment = WD_ALIGN_PARAGRAPH.CENTER
add_page_number(footer_p)

# ---------------- 行内格式 ----------------
TOKEN = re.compile(r"(\*\*.+?\*\*|\*[^*]+?\*|`[^`]+?`|\^[^^]+?\^)")

def soften_long_ids(text):
    """超长的下划线标识符无法断行，会在两端对齐时把整行拉稀；插入零宽断行点"""
    if re.search(r"[A-Za-z0-9]{25,}", text):
        return re.sub(r"(?<=_)(?=[A-Za-z])", "​", text)
    return text

def add_inline(par, text, size=11, base_italic=False, base_bold=False):
    """**粗体**、*斜体*、`变量`（渲染为斜体）、^上标^；支持嵌套"""
    text = soften_long_ids(text)
    for tok in TOKEN.split(text):
        if not tok:
            continue
        if tok.startswith("**") and tok.endswith("**") and len(tok) > 4:
            add_inline(par, tok[2:-2], size, base_italic, True)
        elif tok.startswith("*") and tok.endswith("*") and len(tok) > 2:
            add_inline(par, tok[1:-1], size, True, base_bold)
        elif tok.startswith("`") and tok.endswith("`") and len(tok) > 2:
            r = par.add_run(tok[1:-1])
            r.italic, r.bold, r.font.size = True, base_bold, Pt(size)
        elif tok.startswith("^") and tok.endswith("^") and len(tok) > 2:
            r = par.add_run(tok[1:-1])
            r.font.superscript, r.bold, r.font.size = True, base_bold, Pt(size)
        else:
            r = par.add_run(tok)
            r.italic, r.bold, r.font.size = base_italic, base_bold, Pt(size)

def para(text="", align="justify", size=11, space_after=6, italic=False,
         bold=False, keep=False):
    p = doc.add_paragraph()
    p.alignment = {"justify": WD_ALIGN_PARAGRAPH.JUSTIFY,
                   "center": WD_ALIGN_PARAGRAPH.CENTER,
                   "left": WD_ALIGN_PARAGRAPH.LEFT}[align]
    p.paragraph_format.space_after = Pt(space_after)
    if keep:
        p.paragraph_format.keep_with_next = True
    if text:
        add_inline(p, text, size=size, base_italic=italic, base_bold=bold)
    return p

# ---------------- 表格 ----------------
def add_table(rows):
    ncol = max(len(r) for r in rows)
    t = doc.add_table(rows=0, cols=ncol)
    t.style = "Table Grid"
    t.autofit = True
    for i, row in enumerate(rows):
        cells = t.add_row().cells
        for j in range(ncol):
            txt = row[j] if j < len(row) else ""
            p = cells[j].paragraphs[0]
            p.paragraph_format.space_after = Pt(2)
            p.alignment = WD_ALIGN_PARAGRAPH.LEFT
            add_inline(p, txt, size=8.5, base_bold=(i == 0))
    doc.add_paragraph().paragraph_format.space_after = Pt(4)

def split_row(line):
    line = line.strip()
    if line.startswith("|"):
        line = line[1:]
    if line.endswith("|"):
        line = line[:-1]
    return [c.strip() for c in line.split("|")]

def is_sep_row(line):
    return bool(re.fullmatch(r"\|?[\s:\-|]+\|?", line.strip())) and "-" in line

# ---------------- 解析 ----------------
md = io.open(SRC, encoding="utf-8").read().splitlines()
i = 0
n_tab = n_fig = 0
in_legends = in_refs = False
first_h1 = True

def collect(i, stop):
    """从第 i 行起收集连续行，直到 stop(nxt) 为真或遇到空行/标题"""
    buf, i = [md[i].strip()], i + 1
    while i < len(md):
        nxt = md[i].strip()
        if not nxt or nxt.startswith("#") or stop(nxt):
            break
        buf.append(nxt)
        i += 1
    return re.sub(r"\s+", " ", " ".join(buf)), i

while i < len(md):
    s = md[i].strip()

    if not s or s == "---":
        i += 1
        continue

    if s.startswith("# ") and not s.startswith("## "):
        if first_h1:
            para(s[2:].strip(), align="center", size=18, space_after=18, bold=True)
            first_h1 = False
        else:
            para(s[2:].strip(), align="left", size=15, space_after=10, bold=True, keep=True)
        i += 1
        continue

    if s.startswith("## "):
        h = s[3:].strip()
        if h.startswith("7. Figure legends") or h.startswith("Figure legends"):
            in_legends = True
            doc.add_page_break()
        if h.startswith("8. References") or h.startswith("References"):
            in_refs = True
        para(h, align="left", size=13, space_after=8, bold=True, keep=True)
        i += 1
        continue

    if s.startswith("### "):
        para(s[4:].strip(), align="left", size=11.5, space_after=6, bold=True, keep=True)
        i += 1
        continue

    if s.startswith("|"):
        rows = []
        while i < len(md) and md[i].strip().startswith("|"):
            if not is_sep_row(md[i]):
                rows.append(split_row(md[i]))
            i += 1
        add_table(rows)
        n_tab += 1
        continue

    # 图注：先插图，再整段图注（合并续行）
    if in_legends and s.startswith("**Fig."):
        m = re.match(r"\*\*Fig\.\s*(\d+)\.", s)
        if m:
            img = os.path.join(SUB, "figures", FIGMAP[int(m.group(1))])
            if not os.path.exists(img):
                sys.exit(f"缺少图片 {img}")
            pp = doc.add_paragraph()
            pp.alignment = WD_ALIGN_PARAGRAPH.CENTER
            pp.paragraph_format.space_before = Pt(10)
            pp.paragraph_format.keep_with_next = True
            pp.add_run().add_picture(img, width=Inches(6.1))
            n_fig += 1
        text, i = collect(i, lambda n: n.startswith("**Fig."))
        align = "left" if re.search(r"[A-Za-z0-9_]{30,}", text) else "justify"
        para(text, align=align, size=9.5, space_after=14)
        continue

    if s.startswith("- "):
        p = doc.add_paragraph(style="List Bullet")
        p.paragraph_format.space_after = Pt(3)
        add_inline(p, s[2:].strip(), size=11)
        i += 1
        continue

    if s.startswith("> "):
        text, i = collect(i, lambda n: False)
        para(text[2:] if text.startswith("> ") else text,
             align="justify", size=10, italic=True)
        continue

    # 标题页字段
    if s.startswith("**Authors:**"):
        para(s.split(":", 1)[1].strip().lstrip("*").strip(),
             align="center", size=12, space_after=4)
        i += 1
        continue
    if s.startswith("**Affiliations:**"):
        para(s.split(":", 1)[1].strip().lstrip("*").strip(), align="center", size=10.5,
             space_after=4, italic=True)
        i += 1
        continue
    if s.startswith("**Corresponding author:**"):
        val = s.split(":", 1)[1].strip().lstrip("*").strip()
        para("Corresponding author: " + val, align="center", size=10.5, space_after=12)
        i += 1
        continue
    if s.startswith("**Keywords:**"):
        text, i = collect(i, lambda n: False)
        para(text, align="justify", size=10.5, space_after=10)
        continue

    # 编号条目：左对齐 + 悬挂缩进，条目之间不合并（参考文献与正文列举通用）
    if re.match(r"^\d+\.\s", s):
        text, i = collect(i, lambda n: bool(re.match(r"^\d+\.\s", n)))
        p = para(text, align="left", size=(9.5 if in_refs else 11),
                 space_after=(4 if in_refs else 5))
        p.paragraph_format.left_indent = Cm(0.75)
        p.paragraph_format.first_line_indent = Cm(-0.75)
        if in_refs:
            p.paragraph_format.line_spacing = 1.0
        continue

    # 普通段落
    text, i = collect(i, lambda n: (n.startswith("|") or n.startswith("- ")
                                    or n.startswith("> ")
                                    or bool(re.match(r"^\d+\.\s", n))
                                    or n.startswith("**Fig.") or n.startswith("**Authors:")
                                    or n.startswith("**Affiliations:")
                                    or n.startswith("**Corresponding")
                                    or n.startswith("**Keywords:")))
    align = "left" if re.search(r"[A-Za-z0-9_]{30,}", text) else "justify"
    para(text, align=align, size=11, space_after=7)

doc.save(DOCX)
print(f"docx: {DOCX}")
print(f"  段落 {len(doc.paragraphs)}，表格 {n_tab}，插图 {n_fig}")
src_text = io.open(SRC, encoding="utf-8").read()
exp_tab = len(re.findall(r"^\*\*Table\s*\d+\.", src_text, re.M))
exp_fig = len(re.findall(r"^\*\*Fig\.\s*\d+\.", src_text, re.M))
if n_tab != exp_tab or n_fig != exp_fig:
    sys.exit(f"表格/插图数量异常：表格 {n_tab}/{exp_tab}，插图 {n_fig}/{exp_fig}")

# ---------------- Word 导出 PDF ----------------
DOCX_W, PDF_W = DOCX.replace("/", "\\"), PDF.replace("/", "\\")
ps = os.path.join(SUB, "_export_pdf.ps1")
with io.open(ps, "w", encoding="utf-8-sig", newline="\r\n") as fh:
    fh.write(f'''$ErrorActionPreference = "Stop"
$w = New-Object -ComObject Word.Application
$w.Visible = $false
$w.DisplayAlerts = 0
$doc = $w.Documents.Open("{DOCX_W}")
$doc.SaveAs2("{PDF_W}", 17)
$doc.Close(0)
$w.Quit()
''')
r = subprocess.run(["powershell.exe", "-NoProfile", "-ExecutionPolicy", "Bypass",
                    "-File", ps], capture_output=True, text=True)
os.remove(ps)
if r.returncode != 0 or not os.path.exists(PDF):
    print(r.stdout, r.stderr)
    sys.exit("Word 导出失败")

import fitz
d = fitz.open(PDF)
print(f"pdf : {PDF}")
print(f"  页数 {d.page_count}，大小 {os.path.getsize(PDF)/1024/1024:.2f} MB")
d.close()

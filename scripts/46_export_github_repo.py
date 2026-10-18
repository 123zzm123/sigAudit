#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
46 — 导出可直接上传 GitHub 的仓库结构

原则：只导出代码、配置、图与论文用到的结果表；
原始表达矩阵（615 MB）与日志不进入仓库，改为数据说明。
"""
import io, os, shutil, sys, zipfile, datetime

ROOT = os.environ.get("SIGAUDIT_ROOT", "D:/bioinfo/sigAudit")
OUT  = os.path.join(os.path.dirname(ROOT), "sigAudit_github")
ZIP  = os.path.join(os.path.expanduser("~"), "Desktop", "sigAudit_github_repo_v0.1.0.zip")

AUTHOR = "Ziming Zhao"
EMAIL  = "13317834871@163.com"
AFFIL  = "College of Chemistry, Fuzhou University, Fuzhou, Fujian, China"

if os.path.exists(OUT):
    shutil.rmtree(OUT)
os.makedirs(OUT)

def cp(src, dst=None):
    src = os.path.join(ROOT, src)
    dst = os.path.join(OUT, dst or src.replace(ROOT + os.sep, ""))
    if os.path.isdir(src):
        shutil.copytree(src, dst, dirs_exist_ok=True)
    else:
        os.makedirs(os.path.dirname(dst), exist_ok=True)
        shutil.copy2(src, dst)

# ---------- R 包 ----------
for item in ["DESCRIPTION", "NAMESPACE", "R", "man", "tests", "inst"]:
    cp(item)
cp("third_party/genefu/README.md", "third_party/genefu/README.md")

# ---------- 分析脚本与配置 ----------
cp("scripts")
cp("config")
cp("reproduce_current.sh")
cp("run_v8_parallel.sh")
cp("README_current.md")

# ---------- 图 ----------
os.makedirs(os.path.join(OUT, "results", "figures"), exist_ok=True)
for f in sorted(os.listdir(os.path.join(ROOT, "results", "figures"))):
    if f.startswith("_old"):
        continue
    src = os.path.join(ROOT, "results", "figures", f)
    if os.path.isfile(src) and (f.endswith(".pdf") or f.endswith(".png") or f.endswith(".tsv")):
        shutil.copy2(src, os.path.join(OUT, "results", "figures", f))

# ---------- 论文用到的结果表 ----------
TABLES = [
    "20_summary.tsv", "20_rule_performance_long.tsv", "20_model_metrics_long.tsv",
    "21_entry_effect.tsv", "21_entry_summary.tsv",
    "25_simulation_final_replicates.tsv",
    "26_signature_tasks.tsv", "26_signature_summary.tsv",
    "29_baseline_sensitivity.tsv",
    "31_metabric_model_metrics.tsv", "31_metabric_rule_performance.tsv",
    "32_metabric_baseline_sensitivity.tsv", "36_metabric_prediction_check.tsv",
]
MANIFEST = ["candidate_diversity.tsv", "candidate_top5_dup.tsv"]

os.makedirs(os.path.join(OUT, "results", "tables"), exist_ok=True)
os.makedirs(os.path.join(OUT, "results", "manifest"), exist_ok=True)
for t in TABLES:
    src = os.path.join(ROOT, "results", "tables", t)
    if not os.path.exists(src):
        sys.exit(f"缺少结果表 {t}")
    shutil.copy2(src, os.path.join(OUT, "results", "tables", t))
for t in MANIFEST:
    src = os.path.join(ROOT, "results", "manifest", t)
    if not os.path.exists(src):
        sys.exit(f"缺少清单表 {t}")
    shutil.copy2(src, os.path.join(OUT, "results", "manifest", t))

# ---------- DESCRIPTION 更新 ----------
desc = io.open(os.path.join(OUT, "DESCRIPTION"), encoding="utf-8").read()
desc = desc.replace(
    'Title: Credibility Audit for Omics-Based Prognostic Signatures',
    'Title: Gain Accounting for Multi-Algorithm Search in Omics Prognostics')
desc = desc.replace(
    'Authors@R: person("Anonymous", "Author", email = "author@example.org", role = c("aut", "cre"))',
    f'Authors@R: person("Ziming", "Zhao", email = "{EMAIL}",\n'
    f'    role = c("aut", "cre"))')
new_desc = ('Description: Implements reference-adjusted gain accounting for omics prognostic\n'
            '    modelling. For a candidate pool evaluated on validation and held-out test\n'
            '    partitions, it computes `Gv` (validation gain over a fixed reference), `Gt`\n'
            '    (held-out gain) and `Orel = Gv - Gt` (reference-adjusted selection optimism),\n'
            '    and compares selection rules (argmax, top-k ensembles, resampling remedies)\n'
            '    on a common task set. It separates selection on a noisy criterion from\n'
            '    comparison against a weak baseline, and keeps `Oraw` distinct from `Orel`.')
desc = re.sub(r'Description:.*?(?=\nLicense:)', new_desc, desc, flags=re.S) if False else desc
# 用显式行替换 Description 块
lines = desc.splitlines()
out, skip = [], False
for ln in lines:
    if ln.startswith("Description:"):
        out.append(new_desc)
        skip = True
        continue
    if skip:
        if ln.startswith("    ") or ln.strip() == "":
            continue
        skip = False
    out.append(ln)
desc = "\n".join(out) + "\n"
io.open(os.path.join(OUT, "DESCRIPTION"), "w", encoding="utf-8", newline="\n").write(desc)

# ---------- LICENSE ----------
mit = f"""MIT License

Copyright (c) 2026 {AUTHOR}

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""
io.open(os.path.join(OUT, "LICENSE"), "w", encoding="utf-8", newline="\n").write(mit)
if os.path.exists(os.path.join(OUT, "LICENSE.md")):
    os.remove(os.path.join(OUT, "LICENSE.md"))

# ---------- .gitignore ----------
io.open(os.path.join(OUT, ".gitignore"), "w", encoding="utf-8", newline="\n").write(
""".Rproj.user
.Rhistory
.RData
.Ruserdata
*.Rcheck/
*.tar.gz
data/cohorts/
data/signatures/
logs/
results/main/
results/main_v1_firstattempt/
results/metabric/*.rds
results/remedies/
submission*/
submission_cmpb*/
.DS_Store
Thumbs.db
""")

# ---------- data 说明 ----------
os.makedirs(os.path.join(OUT, "data"), exist_ok=True)
io.open(os.path.join(OUT, "data", "README.md"), "w", encoding="utf-8", newline="\n").write(
f"""# Data

The raw expression matrices and survival tables are **not redistributed** here.
TCGA Pan-Cancer Atlas 2018 and METABRIC data are available from cBioPortal:

- TCGA: `https://www.cbioportal.org/api` — eight study identifiers are listed in Table 1 of
  the manuscript (`brca_tcga_pan_can_atlas_2018`, `coadread_...`, `kirc_...`, `lgg_...`,
  `luad_...`, `meso_...`, `paad_...`, `stad_...`), profile `*_rna_seq_v2_mrna`.
- METABRIC: study `brca_metabric`, profile `brca_metabric_mrna`.

Gene sets: MSigDB v2024.1 (`https://www.gsea-msigdb.org/`), six collections listed in
Supplementary Methods S1.

Signature gene lists were read from `genefu` 2.44.0 (Bioconductor) and are not redistributed;
`third_party/genefu/README.md` records the objects used and their MD5 checksums.

After download, place files under `data/cohorts/` as described in `reproduce_current.sh`.
""")

# ---------- results 说明 ----------
io.open(os.path.join(OUT, "results", "README.md"), "w", encoding="utf-8", newline="\n").write(
"""# Results

- `figures/` — final figures (PDF and 300-dpi PNG) and the figure data behind them.
- `tables/` — result tables used in the manuscript:
  - `20_summary.tsv` — rule-level summary (15 rows)
  - `20_rule_performance_long.tsv` — task-level record for every rule (58,785 rows)
  - `20_model_metrics_long.tsv` — candidate-model evaluations (66,483 rows)
  - `21_entry_effect.tsv`, `21_entry_summary.tsv` — entry effects
  - `25_simulation_final_replicates.tsv` — simulation replicates
  - `26_signature_tasks.tsv`, `26_signature_summary.tsv` — published-signature module
  - `29_baseline_sensitivity.tsv` — three fixed references
  - `31_*`, `32_*`, `36_*` — METABRIC metrics, sensitivity and re-verification
- `manifest/` — candidate-pool diversity tables.

Numbers in the manuscript are generated from these files; see `scripts/` for the pipeline.
""")

# ---------- README ----------
io.open(os.path.join(OUT, "README.md"), "w", encoding="utf-8", newline="\n").write(
f"""# sigAudit

Reference-adjusted gain accounting for multi-algorithm search in omics prognostic modelling.

Multi-algorithm pipelines fit many learners, keep the one with the best validation
performance, and report that number. `sigAudit` separates the two things that number mixes:
the gain that survives on held-out data, and the gain that depends on the comparator.

## What it computes

For a candidate pool evaluated on a validation and a held-out test partition:

| Quantity | Definition | Meaning |
|---|---|---|
| `Gv` | `Cv(rule) − Cv(reference)` | validation gain over a fixed reference |
| `Gt` | `Ct(rule) − Ct(reference)` | held-out gain over the same reference |
| `Orel` | `Gv − Gt` | reference-adjusted selection optimism |
| `Oraw` | `Cv(rule) − Ct(rule)` | rule's own validation–test difference; never conflated with `Orel` |

`audit_rules()` compares selection rules (argmax, top-k ensembles, resampling remedies) on a
common task set; `audit_selection()` performs the paired accounting for a single rule.

## Install

```r
# install.packages("remotes")
remotes::install_github("<your-account>/sigAudit")
```

## Quick start

```r
library(sigAudit)

# val_ci, test_ci: matrices of C-index, rows = tasks, columns = candidates
res <- audit_rules(val_ci, test_ci, reference = "Lasso")
res$rules        # per-rule Gv, Gt, Orel, Oraw
res$members      # ensemble membership
```

The package ships assertion-based tests covering column-order invariance of tie-breaking,
row-identity checking, rejection of out-of-range C-indices, and the identity between the
comparable-pair count and `concordant + discordant + tied.x`.

## Reproducing the study

```bash
bash reproduce_current.sh tables      # rebuild TCGA tables from saved predictions
bash reproduce_current.sh metabric    # METABRIC: refit, sensitivity, figure
bash reproduce_current.sh signatures  # published-signature module
bash reproduce_current.sh train       # main training run + simulation
bash reproduce_current.sh check       # unit tests + R CMD check + hash freeze
```

Data are not redistributed; see `data/README.md` for access instructions and
`scripts/36_verify_metabric_predictions.R` for an independent reconstruction check.

## Repository layout

```
R/                     package source
man/, NAMESPACE         package documentation and exports
tests/                  assertion-based tests
scripts/                analysis pipeline (numbered by stage)
config/                 frozen analysis plans (fixed before the runs)
results/figures         final figures and figure data
results/tables          result tables used in the manuscript
data/README.md          how to obtain the public data
```

## Citation

If you use this software, please cite the accompanying article (under review). A
`CITATION.cff` file is included; it will be updated with the DOI on publication.

## License

MIT — see `LICENSE`.

## Contact

{AUTHOR}, {AFFIL} — {EMAIL}
""")

# ---------- CITATION ----------
io.open(os.path.join(OUT, "CITATION.cff"), "w", encoding="utf-8", newline="\n").write(
f"""cff-version: 1.2.0
message: "If you use this software, please cite the accompanying article."
title: "sigAudit"
version: 0.1.0
date-released: 2026-10-18
license: MIT
authors:
  - family-names: "Zhao"
    given-names: "Ziming"
    affiliation: "College of Chemistry, Fuzhou University, Fuzhou, Fujian, China"
    email: "{EMAIL}"
abstract: >-
  Reference-adjusted gain accounting for multi-algorithm search in omics prognostic
  modelling. Computes validation gain, held-out gain and reference-adjusted selection
  optimism, and compares selection rules on a common task set.
""")

# ---------- 统计与打包 ----------
total = 0
nfiles = 0
for root, dirs, fs in os.walk(OUT):
    for f in fs:
        total += os.path.getsize(os.path.join(root, f))
        nfiles += 1

with zipfile.ZipFile(ZIP, "w", zipfile.ZIP_DEFLATED, compresslevel=6) as z:
    for root, dirs, fs in os.walk(OUT):
        for f in fs:
            p = os.path.join(root, f)
            z.write(p, os.path.relpath(p, OUT))

print(f"仓库目录: {OUT}")
print(f"  {nfiles} 文件, {total/1024/1024:.1f} MB")
print(f"压缩包: {ZIP}  ({os.path.getsize(ZIP)/1024/1024:.1f} MB)")
print("顶层内容:", sorted(os.listdir(OUT)))

#!/bin/bash
# ============================================================
# 28_freeze_current.sh —— 冻结当前版（v6/19-27 流水线）的全部产物哈希。
#
# 为什么需要：第三轮审查指出 RDS 元数据里只有"输入哈希前缀 + 脚本名 + 完成标志"，
# 没有生成该产物时的【脚本哈希 / 基因集 JSON 哈希 / 配置哈希 / 环境版本】，
# 因此 README 里当前的脚本 MD5 不能证明历史产物由同一份脚本生成。
#
# 本脚本不重训任何东西，只【读取并记录】已有文件的哈希，
# 生成一份不可歧义的构建记录：results/manifest/frozen_v1_<日期>.tsv
# 以及环境记录 results/manifest/frozen_v1_<日期>_env.txt
#
# 注意：这是【事后构建记录】，不是运行时记录。审查明确要求不得冒充运行时记录，
# 因此文件里显式标注 frozen_at 与 recorded_retrospectively=TRUE。
#
# 用法: bash scripts/28_freeze_current.sh
# ============================================================
set -u
cd "$(dirname "$0")/.."
ROOT="$(pwd -W 2>/dev/null || pwd)"
# 标签可由参数指定：v1 = 修订前的基线冻结（用于证明"没有重训"），v2 = 本轮修订后的状态
TAG="${1:-frozen_v1_20261017}"
OUT="results/manifest/${TAG}.tsv"
ENV="results/manifest/${TAG}_env.txt"

echo "项目根目录: $ROOT"
echo "输出: $OUT"

{
  printf "path\tmd5\tbytes\tmtime\tcategory\n"
} > "$OUT"

emit() {  # $1=文件, $2=类别
  local f="$1" c="$2"
  [ -f "$f" ] || { echo "  [缺] $f" >&2; return 0; }
  local m s t
  m=$(md5sum "$f" | awk '{print $1}')
  s=$(stat -c %s "$f")
  t=$(stat -c %y "$f" | cut -d. -f1)
  printf "%s\t%s\t%s\t%s\t%s\n" "$f" "$m" "$s" "$t" "$c" >> "$OUT"
}

# --- 1. 主实验产物：40 个 RDS ---
for f in results/main/19_*_100.rds; do emit "$f" "main_rds"; done

# --- 2. 当前版分析表 ---
for f in results/tables/20_*.tsv results/tables/21_*.tsv \
         results/tables/25_*.tsv results/tables/26_*.tsv; do
  emit "$f" "table_current"
done

# --- 3. 清单 ---
for f in results/manifest/*.tsv results/manifest/main_input_files.txt; do
  case "$f" in
    *"$TAG"*|*frozen*) continue ;;
  esac
  emit "$f" "manifest"
done

# --- 4. 生成当前版的脚本（不含被取代的 _ 前缀与旧版）---
for f in scripts/19_main_corrected.R scripts/20_analysis_corrected.R \
         scripts/21_entry_predictors.R scripts/23_manifest.R \
         scripts/24_manifests.R scripts/25_simulation.R \
         scripts/26a_extract_public_signatures.R scripts/26b_public_signature_tasks.R \
         scripts/27_status_readme.R scripts/13_figures.R \
         scripts/01b_build_genesets.py scripts/04_entrez2symbol.py \
         scripts/02_fetch_cohorts.py; do
  emit "$f" "script"
done

# --- 5. R 包源码（论文软件产物）---
for f in R/*.R DESCRIPTION NAMESPACE; do emit "$f" "package_src"; done
for f in tests/testthat/*.R; do [ -f "$f" ] && emit "$f" "package_test"; done

# --- 6. 配置与设计文件 ---
for f in config/*.md config/*.txt config/*.tsv config/*.json; do
  [ -f "$f" ] && emit "$f" "config"
done

# --- 7. 数据输入（哈希前缀级，完整文件太大只记大小与哈希）---
for f in data/cohorts/*_expr.tsv data/cohorts/*_surv.tsv data/signatures/*.json; do
  [ -f "$f" ] && emit "$f" "data_input"
done

# --- 8. 运行脚本 ---
for f in run_v8_parallel.sh run_v8_addworker.sh run_v8_main.sh run_v8_resume.sh; do
  [ -f "$f" ] && emit "$f" "runner"
done

N=$(($(wc -l < "$OUT") - 1))
echo "已记录 $N 个文件的哈希"

# --- 环境记录 ---
R="D:/R-4.6.1/bin/Rscript.exe"
{
  echo "# 当前版构建记录（事后补录，非运行时记录）"
  echo "frozen_at: $(date '+%Y-%m-%d %H:%M:%S %z')"
  echo "recorded_retrospectively: TRUE"
  echo "note: 本文件由 scripts/28_freeze_current.sh 生成，只做哈希固化，不重训任何模型。"
  echo "note: 这些哈希是【冻结基线】；后续任何重算都必须与这里的值对照，说明哪些表变了、为什么。"
  echo
  echo "## 环境"
  date
  uname -a 2>/dev/null
  echo
  echo "## R sessionInfo"
} > "$ENV"
if [ -x "$R" ]; then
  "$R" -e 'cat(R.version.string, "\n"); cat("平台:", R.version$platform, "\n\n"); print(sessionInfo())' >> "$ENV" 2>&1
else
  echo "R 未在 $R 找到" >> "$ENV"
fi

echo "环境记录: $ENV"
echo "完成。"

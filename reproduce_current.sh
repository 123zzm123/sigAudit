#!/bin/bash
# ============================================================
# reproduce_current.sh —— 当前版（19–29）的端到端复现入口
#
# 第四轮审查指出的两个问题的修正：
#  ① "shell 换目录不等于脚本能在新根目录工作" —— 旧版虽然允许指定 ROOT，
#     但下游 R 脚本内部把 PROJ 写死成 D:/bioinfo/sigAudit。
#     现在：本脚本 `export SIGAUDIT_ROOT="$ROOT"`，**所有当前管线脚本**
#     都以 `PROJ <- Sys.getenv("SIGAUDIT_ROOT", "<默认>")` 读取根目录，
#     项目外部的绝对路径（LGG 表达/生存、genefu 数据对象）也已改为项目内路径。
#  ② "tables 路径调用 26a/26b，其中涉及公开签名模型训练，
#     不能整体标为'仅从已有预测生成表、无需重训'" ——
#     现在把路径**明确拆开**：
#       tables     只重建【不需要重新拟合任何模型】的表与图
#       signatures 26a+26b，**会在签名基因集上重新拟合候选池**（属训练）
#       train      从原始数据重训 40 组主实验
#       all        train + signatures + tables
#       check      单元测试 + R CMD check + 哈希固化（不重算分析）
#       metabric   METABRIC 独立队列的流程复现（31→32→33）——**会训练**（30 个任务，约 1 分钟，12 worker）
#
# 用法：
#   bash reproduce_current.sh [tables|signatures|metabric|train|all|check] [ROOT]
#   或  SIGAUDIT_ROOT=/path/to/sigAudit bash reproduce_current.sh tables
#
# 不传 ROOT 时使用本脚本所在目录。
# ============================================================
set -u
SELF_DIR="$(cd "$(dirname "$0")" && { pwd -W 2>/dev/null || pwd; })"
MODE="${1:-tables}"
ROOT="${2:-${SIGAUDIT_ROOT:-$SELF_DIR}}"
cd "$ROOT" || { echo "无法进入根目录: $ROOT"; exit 1; }

# 关键：把根目录导出给所有子进程（R 脚本据此定位 data/ results/ logs/ config/）
export SIGAUDIT_ROOT="$ROOT"

R="${SIGAUDIT_R:-D:/R-4.6.1/bin/Rscript.exe}"
RCMD="${SIGAUDIT_RCMD:-D:/R-4.6.1/bin/R.exe}"
PY="${SIGAUDIT_PY:-D:/Python/python.exe}"

LOGDIR="logs"; mkdir -p "$LOGDIR"
RUNLOG="$LOGDIR/reproduce_current_$(date +%Y%m%d_%H%M%S).log"
PROV="$LOGDIR/reproduce_current_provenance.txt"

step() { echo; echo "=========== $* ===========" | tee -a "$RUNLOG"; }
run()  {  # run <可读名> <命令...>
  local name="$1"; shift
  step "$name"
  echo "  \$ $*" | tee -a "$RUNLOG"
  local t0; t0=$(date +%s)
  if "$@" >> "$RUNLOG" 2>&1; then
    echo "  [OK]   $name  ($(( $(date +%s) - t0 ))s)" | tee -a "$RUNLOG"
    return 0
  else
    local rc=$?
    echo "  [FAIL] $name  rc=$rc  ($(( $(date +%s) - t0 ))s)" | tee -a "$RUNLOG"
    FAILED="${FAILED:-} $name(rc=$rc)"
    return $rc
  fi
}
need() { for f in "$@"; do [ -e "$f" ] || { echo "  [缺] $f"; return 1; }; done; }

echo "根目录: $ROOT"
echo "模式  : $MODE"
echo "R     : $R"
echo "日志  : $RUNLOG"

# ---------- 0. 环境与来源记录 ----------
step "0. 环境与来源记录"
{
  echo "# reproduce_current.sh 运行记录"
  echo "mode: $MODE"; echo "root: $ROOT"
  echo "started_at: $(date '+%Y-%m-%d %H:%M:%S %z')"
  echo "host: $(uname -a 2>/dev/null)"
  echo; echo "## 关键脚本哈希（本次运行所用的版本）"
} > "$PROV"
for f in scripts/19_main_corrected.R scripts/20_analysis_corrected.R scripts/21_entry_predictors.R \
         scripts/23_manifest.R scripts/24_manifests.R scripts/25_simulation.R \
         scripts/26a_extract_public_signatures.R scripts/26b_public_signature_tasks.R \
         scripts/29_baseline_sensitivity.R scripts/13_figures.R scripts/27_status_readme.R \
         scripts/_lib_analysis.R config/frozen_analysis_plan.md; do
  [ -f "$f" ] && printf "%-50s %s\n" "$f" "$(md5sum "$f" | awk '{print $1}')" >> "$PROV"
done
[ -x "$R" ] && "$R" -e 'cat("\n## R sessionInfo\n"); cat(R.version.string,"\n"); print(sessionInfo())' >> "$PROV" 2>&1
echo "  来源记录: $PROV"

FAILED=""

case "$MODE" in

tables)
  # 只重建【不需要重新拟合任何模型】的表与图。
  # 唯一前置条件 = 已保存的逐患者预测（results/main/19_*_100.rds）；
  # 23 会据此重建 main_input_files.txt，20 会重建 model_metrics_long 等长表。
  if ! ls results/main/19_*_100.rds >/dev/null 2>&1; then
    echo "  [缺] results/main/19_*_100.rds（已保存的逐患者预测）"
    echo "  缺少主实验产物；请先运行： bash reproduce_current.sh train"; exit 3
  fi
  echo "  找到已保存的预测文件: $(ls results/main/19_*_100.rds | wc -l) 个"
  run "23 清单与最终输入列表"        "$R" scripts/23_manifest.R
  run "20 主分析（重建全部主表）"     "$R" scripts/20_analysis_corrected.R
  run "21 进入效应与候选特征"         "$R" scripts/21_entry_predictors.R
  run "29 固定基线敏感性分析"         "$R" scripts/29_baseline_sensitivity.R
  # METABRIC 的 32 只重算表、不训练（31 才训练），因此放在 tables 路径
  if [ -f results/tables/31_metabric_model_metrics.tsv ]; then
    run "32 METABRIC 参照敏感性（仅重算）" "$R" scripts/32_metabric_baseline_sensitivity.R
    run "33 METABRIC 补充图"               "$R" scripts/33_figure_metabric.R
    if [ -f results/metabric/31_metabric_predictions.rds ]; then
      run "36 METABRIC 预测核验（只读，不训练）" "$R" scripts/36_verify_metabric_predictions.R
    fi
  fi
  run "24 样本/基因集/候选多样性清单"  "$R" scripts/24_manifests.R
  run "13 论文主图与图注"             "$R" scripts/13_figures.R
  run "27 状态 README"                "$R" scripts/27_status_readme.R
  echo
  echo "说明：本模式【不训练任何模型】。"
  echo "  · 公开签名模块（26a/26b）**会在签名基因集上重新拟合候选池**，"
  echo "    因此单独列为签名路径： bash reproduce_current.sh signatures"
  echo "  · 模拟（25）与主表无关且耗时： Rscript scripts/25_simulation.R final 100"
  ;;

signatures)
  # 公开签名衍生任务：会在签名的基因集上**重新训练**候选池。
  run "26a 提取公开签名基因列表"      "$R" scripts/26a_extract_public_signatures.R
  run "26b 公开签名衍生任务（含训练）" "$R" scripts/26b_public_signature_tasks.R
  ;;

train)
  if [ ! -f data/cohorts/LGG_expr.tsv ]; then
    run "02 获取 TCGA 队列" "$PY" scripts/02_fetch_cohorts.py
    run "04 Entrez -> symbol" "$PY" scripts/04_entrez2symbol.py
  fi
  [ -f data/signatures/msigdb_genesets.json ] || run "01b 构建基因集" "$PY" scripts/01b_build_genesets.py
  run "19 主实验 8 队列 x 5 划分（并行）" bash run_v8_parallel.sh
  run "25 模拟（final 100，并行）"        "$R" scripts/25_simulation.R final 100
  run "23 清单与最终输入列表"             "$R" scripts/23_manifest.R
  run "20 主分析"                         "$R" scripts/20_analysis_corrected.R
  run "21 进入效应与候选特征"              "$R" scripts/21_entry_predictors.R
  run "29 固定基线敏感性分析"              "$R" scripts/29_baseline_sensitivity.R
  run "24 清单"                           "$R" scripts/24_manifests.R
  run "13 论文主图"                       "$R" scripts/13_figures.R
  run "27 状态 README"                    "$R" scripts/27_status_readme.R
  ;;

all)
  run "训练"     bash "$0" train "$ROOT"
  run "签名模块" bash "$0" signatures "$ROOT"
  run "表与图"   bash "$0" tables "$ROOT"
  ;;

metabric)
  # 独立队列的**流程复现**：在 METABRIC 内部重新划分、重新拟合。
  # 注意：本模式**会训练**（与 tables 不同）。配置见 config/metabric_frozen_config.md。
  if ! need data/cohorts/METABRIC_expr.tsv data/cohorts/METABRIC_surv.tsv; then
    echo "  缺少 METABRIC 数据；请先按 docs/metabric_availability_check.md 取得数据"; exit 3
  fi
  run "31 METABRIC 审计（30 任务，含训练）" "$R" scripts/31_metabric_audit.R
  run "32 METABRIC 参照敏感性"             "$R" scripts/32_metabric_baseline_sensitivity.R
  run "33 METABRIC 补充图"                 "$R" scripts/33_figure_metabric.R
  ;;

check)
  run "R 包单元测试" "$R" -e 'pkgload::load_all(".", quiet=TRUE); library(testthat); test_dir("tests/testthat")'
  run "R CMD build"  "$RCMD" CMD build . --no-build-vignettes --no-manual
  run "R CMD check"  "$RCMD" CMD check "$(ls -t sigAudit_*.tar.gz | head -1)" --no-manual --no-build-vignettes
  run "哈希固化（冻结当前产物）" bash scripts/28_freeze_current.sh "frozen_$(date +%Y%m%d)"
  ;;

*)
  echo "未知模式: $MODE （可用: tables | signatures | metabric | train | all | check）"; exit 2 ;;
esac

echo
if [ -n "$FAILED" ]; then
  echo "有步骤失败:$FAILED" | tee -a "$RUNLOG"
  echo "详见 $RUNLOG"; exit 1
fi
echo "全部完成。表: results/tables/   图: results/figures/   日志: $RUNLOG"
{ echo "finished_at: $(date '+%Y-%m-%d %H:%M:%S %z')"; echo "status: success"; } >> "$PROV"

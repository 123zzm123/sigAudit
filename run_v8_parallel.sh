#!/bin/bash
# ============================================================
# v8 并行版：同一份注册网格，最多 NWORK 个 (队列×种子) 同时跑
#
# 为什么并行：本机 16 物理核 / 32 线程，原链串行只用一个核；
#   单个 19 进程实测约 676 MB 内存。串行预计 13.4 小时，4 路并行约 4 小时。
#
# 科学语义与串行版完全一致：
#   · 注册网格 8 队列 × 5 种子 = 40 个组合，一个不多一个不少；
#   · "完成"的判据仍是【走完全部计划基因集才写最终文件名】；
#   · 任一任务失败 → 记录 FAIL，最后不做汇总并以非零码退出；
#   · 结束前核验 40 个组合齐全、无重复。
#
# 调度：按预计耗时从大到小入队（LPT），减少尾部长任务造成的空等。
# 并发上限靠 NWORK 控制：单进程实测约 676 MB，4 路约 2.7 GB，本机可用内存充足。
#
# 用法: bash run_v8_parallel.sh [并发数=4]
# ============================================================
cd /d/bioinfo/sigAudit
R="D:/R-4.6.1/bin/Rscript.exe"

N=100
NWORK=${1:-4}

LOG=logs/chain_v8_par.log
LOCK=logs/v8_queue.lock
QUEUE=logs/v8_queue.txt
FAILLOG=logs/v8_failed.txt
rm -rf "$LOCK"
: > "$LOG"
: > "$FAILLOG"
: > "$QUEUE"

# 队列按预计耗时降序（实测分钟数）；只把尚未完成的组合入队
for pair in BRCA:52 COAD:32 KIRC:25 LUAD:23 LGG:22 STAD:13 PAAD:5 MESO:2; do
  C=${pair%%:*}
  for S in 20261014 20261016 20261017 20261018 20261019; do
    if [ ! -f "results/main/19_${C}_${S}_${N}.rds" ]; then echo "$C $S" >> "$QUEUE"; fi
  done
done

n_left() { wc -l < "$QUEUE" 2>/dev/null || echo 0; }

echo "[par] 启动 $(date '+%H:%M:%S')  并发=$NWORK  待办=$(n_left)" >> "$LOG"

worker() {
  local wid="$1"
  local job C S OUT rc guard
  while :; do
    # 用 mkdir 做原子锁来取任务
    if ! mkdir "$LOCK" 2>/dev/null; then sleep 2; continue; fi
    job=$(head -1 "$QUEUE" 2>/dev/null)
    if [ -z "$job" ]; then rmdir "$LOCK"; break; fi
    sed -i '1d' "$QUEUE"
    rmdir "$LOCK"

    C=$(echo "$job" | awk '{print $1}')
    S=$(echo "$job" | awk '{print $2}')
    OUT="results/main/19_${C}_${S}_${N}.rds"
    if [ -f "$OUT" ]; then continue; fi

    echo "[run ] W$wid $C seed=$S $(date '+%H:%M:%S')  剩余$(n_left)" >> "$LOG"
    "$R" scripts/19_main_corrected.R "$N" "$C" "$S" >> "logs/19_${C}_${S}_console.log" 2>&1
    rc=$?
    if [ "$rc" -ne 0 ] || [ ! -f "$OUT" ]; then
      echo "[FAIL] W$wid $C seed=$S rc=$rc 有文件=$([ -f "$OUT" ] && echo yes || echo no) $(date '+%H:%M:%S')" >> "$LOG"
      echo "$C $S" >> "$FAILLOG"
    else
      echo "[done] W$wid $C seed=$S $(date '+%H:%M:%S')  剩余$(n_left)" >> "$LOG"
    fi
  done
}

for w in $(seq 1 "$NWORK"); do
  worker "$w" &
done
wait

echo "[par] 本脚本的 worker 已结束，等待可能存在的追加 worker $(date '+%H:%M:%S')" >> "$LOG"
# 追加 worker（run_v8_addworker.sh）共享同一队列；必须等它们也退出，
# 否则会在产物尚未齐全时做完整性核验而误判失败。
for i in $(seq 1 240); do
  still=$(powershell -NoProfile -Command "(Get-CimInstance Win32_Process | Where-Object { $_.CommandLine -match '19_main_corrected' } | Measure-Object).Count" 2>/dev/null | tr -d ' ')
  [ "$still" = "0" ] && sleep 20 && still=$(powershell -NoProfile -Command "(Get-CimInstance Win32_Process | Where-Object { $_.CommandLine -match '19_main_corrected' } | Measure-Object).Count" 2>/dev/null | tr -d ' ')
  [ "$still" = "0" ] && break
  sleep 30
done
echo "[par] 所有 worker（含追加）结束 $(date '+%H:%M:%S')" >> "$LOG"

if [ -s "$FAILLOG" ]; then
  echo "[ABORTED] 有任务失败，不做汇总：" >> "$LOG"
  cat "$FAILLOG" >> "$LOG"
  exit 1
fi

NEXP=40
NGOT=$(ls results/main/19_*_[0-9]*.rds 2>/dev/null | grep -v partial | wc -l)
DUPS=$(ls results/main/19_*_[0-9]*.rds 2>/dev/null | grep -v partial | sed -E 's#.*/19_([A-Z]+)_[0-9]+_[0-9]+\.rds#\1#' | sort | uniq -c | awk '$1 > 1' | wc -l)
echo "[check] 期望 $NEXP 实际 $NGOT 重复组合 $DUPS" >> "$LOG"
if [ "$NGOT" -ne "$NEXP" ] || [ "$DUPS" -ne 0 ]; then
  echo "[ABORTED] 产物数量或唯一性不符" >> "$LOG"
  exit 1
fi

"$R" scripts/23_manifest.R           >> logs/23_console.log 2>&1 || { echo "[FAIL] manifest" >> "$LOG"; exit 1; }
"$R" scripts/20_analysis_corrected.R >> logs/20_console.log 2>&1 || { echo "[FAIL] 20" >> "$LOG"; exit 1; }
"$R" scripts/21_entry_predictors.R   >> logs/21_console.log 2>&1 || { echo "[FAIL] 21" >> "$LOG"; exit 1; }
echo "[ALL DONE] $(date '+%Y-%m-%d %H:%M:%S')" >> "$LOG"

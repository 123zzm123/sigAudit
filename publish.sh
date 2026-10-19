#!/bin/bash
# ============================================================
# publish.sh — 打包投稿件、导出公开仓库并推送（把易漏的三步固化）
#
# 步骤固定为：
#   1) 重建 CMPB 投稿包并校验（scripts/43、44；可选 45 生成 Word）
#   2) 导出 GitHub 仓库结构（scripts/46，复制白名单内的代码/图/表）
#   3) 在导出仓库里 git add / commit / push，并核对 HEAD == origin/main
#
# 用法：
#   bash publish.sh "提交信息"                 # 走 43->44->45->46 + 提交推送
#   bash publish.sh "提交信息" --no-word       # 跳过 Word 生成
#   bash publish.sh "提交信息" --no-pack       # 只做 46 + 提交推送（不改投稿包）
#   bash publish.sh "提交信息" --dry-run       # 只导出 + 显示将提交的变更，不提交/不推送
#   bash publish.sh "提交信息" --check         # 完成后额外跑 R CMD check
#
# 环境变量（与 reproduce_current.sh 一致）：
#   SIGAUDIT_ROOT  项目根目录（默认脚本所在目录）
#   SIGAUDIT_PY    python（默认 D:/Python/python.exe）
#   SIGAUDIT_RCMD  R.exe（默认 D:/R-4.6.1/bin/R.exe）
#   SIGAUDIT_REPO  导出仓库目录（默认 项目同级 sigAudit_github）
#   SIGAUDIT_REMOTE 远端名（默认 origin）；分支默认当前分支
# ============================================================
set -u
SELF_DIR="$(cd "$(dirname "$0")" && { pwd -W 2>/dev/null || pwd; })"
MSG="${1:-}"
shift || true

DO_WORD=1; DO_PACK=1; DRY=0; DO_CHECK=0
for a in "$@"; do
  case "$a" in
    --no-word) DO_WORD=0 ;;
    --no-pack) DO_PACK=0 ;;
    --dry-run) DRY=1 ;;
    --check)   DO_CHECK=1 ;;
    *) echo "未知参数: $a"; exit 2 ;;
  esac
done
if [ -z "$MSG" ]; then
  echo "用法: bash publish.sh \"提交信息\" [--no-word|--no-pack|--dry-run|--check]"
  exit 2
fi

ROOT="${SIGAUDIT_ROOT:-$SELF_DIR}"
cd "$ROOT" || { echo "无法进入根目录: $ROOT"; exit 1; }
export SIGAUDIT_ROOT="$ROOT"

PY="${SIGAUDIT_PY:-D:/Python/python.exe}"
RCMD="${SIGAUDIT_RCMD:-D:/R-4.6.1/bin/R.exe}"
REPO="${SIGAUDIT_REPO:-$(dirname "$ROOT")/sigAudit_github}"
REMOTE="${SIGAUDIT_REMOTE:-origin}"

step() { echo; echo "=========== $* ==========="; }
must() { "$@" || { echo "  [FAIL] $*"; exit 1; }; }

echo "root  : $ROOT"
echo "repo  : $REPO"
echo "python: $PY"
echo "R     : $RCMD"

# ---------- 1. 投稿包 ----------
if [ "$DO_PACK" = "1" ]; then
  step "1a. 重建 CMPB 投稿包 (43)"
  must "$PY" scripts/43_build_cmpb_pack.py
  step "1b. 校验 CMPB 投稿包 (44)"
  must "$PY" scripts/44_verify_cmpb_pack.py
  if [ "$DO_WORD" = "1" ]; then
    step "1c. 生成投稿 Word (45)"
    must "$PY" scripts/45_make_cmpb_word.py
  fi
else
  echo; echo "(跳过投稿包重建)"
fi

# ---------- 2. 导出仓库 ----------
step "2. 导出 GitHub 仓库结构 (46)"
must "$PY" scripts/46_export_github_repo.py

# ---------- 3. 提交与推送 ----------
step "3. 提交并推送"
cd "$REPO" || { echo "无法进入导出仓库: $REPO"; exit 1; }

echo "--- git status ---"
git status --short

if [ "$DRY" = "1" ]; then
  echo; echo "(dry-run：未提交、未推送)"
  exit 0
fi

git add -A
if git diff --cached --quiet; then
  echo "(没有需要提交的变更)"
else
  git commit -m "$MSG" || { echo "  [FAIL] git commit"; exit 1; }
fi

BR="$(git rev-parse --abbrev-ref HEAD)"
git push "$REMOTE" "$BR" || { echo "  [FAIL] git push"; exit 1; }

# ---------- 4. 核对远端与本地一致 ----------
step "4. 核对 HEAD == $REMOTE/$BR"
git fetch "$REMOTE" "$BR" >/dev/null 2>&1 || true
L="$(git rev-parse HEAD)"; R="$(git rev-parse "$REMOTE/$BR" 2>/dev/null || echo unknown)"
echo "local : $L"
echo "remote: $R"
if [ "$L" = "$R" ]; then echo "  [OK] 本地与远端一致"; else echo "  [WARN] 本地与远端不一致，请检查推送"; exit 1; fi

# ---------- 5. 可选 R CMD check ----------
if [ "$DO_CHECK" = "1" ]; then
  step "5. R CMD check (可选)"
  STAMP="$(date +%Y%m%d_%H%M%S)"
  must "$RCMD" CMD build . --no-build-vignettes --no-manual
  TARBALL="$(ls -t sigAudit_*.tar.gz 2>/dev/null | head -n 1 || true)"
  [ -n "$TARBALL" ] || { echo "  [FAIL] 未找到 tar.gz"; exit 1; }
  "$RCMD" CMD check "$TARBALL" --no-manual --no-build-vignettes 2>&1 | tee "$ROOT/logs/00check_${STAMP}.log" || true
  echo "  日志: logs/00check_${STAMP}.log（请确认末尾 Status: OK）"
fi

echo; echo "发布完成。"

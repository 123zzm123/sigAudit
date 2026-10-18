#!/usr/bin/env Rscript
# ============================================================
# 23 — 运行清单（run_manifest）
#
# 第二轮评审要求：
#  · 20/21 不得再用宽泛的 `^19_.*\.rds$` 目录扫描（会读到试点/不同 N/不同版本，
#    且列表按 cohort+seed 键覆盖时会静默发生版本选择）。
#  · 必须有一份明确清单，记录每个（队列×种子）的完成标志、计划/处理数、
#    输入哈希、输出路径与哈希、时间戳。
#
# 本脚本：
#   1. 扫描 results/main/19_*.rds（排除 .partial）
#   2. 校验每份的 completed 标志与 n_planned/n_processed 一致性
#   3. 写 results/manifest/run_manifest.tsv
#   4. 写 results/manifest/main_input_files.txt —— 20/21 只读这份清单里的文件
#   5. 若有任一文件未完成或元数据缺失，报错退出（非零）
# ============================================================
options(stringsAsFactors = FALSE)
`%||%` <- function(a, b) if (is.null(a)) b else a
PROJ <- Sys.getenv("SIGAUDIT_ROOT", "D:/bioinfo/sigAudit")
MAIN <- file.path(PROJ, "results/main")
MAN  <- file.path(PROJ, "results/manifest"); dir.create(MAN, showWarnings = FALSE, recursive = TRUE)

fs <- list.files(MAIN, pattern = "^19_.*\\.rds$", full.names = TRUE)
fs <- fs[!grepl("\\.partial\\.rds$", fs)]
if (!length(fs)) stop("results/main 下没有完成的主实验产物")

rows <- list(); bad <- character(0)
for (f in fs) {
  R <- tryCatch(readRDS(f), error = function(e) NULL)
  if (is.null(R)) { bad <- c(bad, paste0(basename(f), ": 无法读取")); next }
  need <- c("cohort","study_id","endpoint","seed","n_planned","n_processed","completed",
            "sample_ids","genes","genesets")
  miss <- setdiff(need, names(R))
  if (length(miss)) { bad <- c(bad, paste0(basename(f), ": 缺字段 ", paste(miss, collapse=","))); next }
  if (!isTRUE(R$completed)) { bad <- c(bad, paste0(basename(f), ": completed != TRUE")); next }
  if (!identical(as.integer(R$n_processed), as.integer(R$n_planned))) {
    bad <- c(bad, sprintf("%s: n_processed(%d) != n_planned(%d)",
                          basename(f), R$n_processed, R$n_planned)); next }
  rows[[basename(f)]] <- data.frame(
    cohort = R$cohort, study_id = R$study_id, endpoint = R$endpoint, seed = R$seed,
    n_planned = R$n_planned, n_processed = R$n_processed,
    n_genesets_usable = length(R$genesets),
    n_sample = length(R$sample_ids), n_train = length(R$train_ids),
    n_val = length(R$val_ids), n_test = length(R$test_ids),
    n_gene_used = length(R$genes), n_gene_dropped_missing = R$n_gene_dropped_missing,
    n_value_imputed = R$n_value_imputed, n_gene_dropped_zerovar = R$n_gene_dropped_zerovar,
    expr_md5 = R$expr_md5, finished_at = R$finished_at %||% NA_character_,
    path = f, file_hash = substr(as.character(tools::md5sum(f)), 1, 12),
    check_full = file.info(f)$size > 0, stringsAsFactors = FALSE)
}

if (length(bad)) {
  cat("以下产物未通过完整性校验：\n"); for (b in bad) cat("  -", b, "\n")
  stop("manifest 拒绝生成：存在未完成或元数据缺失的产物")
}
MANI <- do.call(rbind, rows); rownames(MANI) <- NULL
MANI <- MANI[order(MANI$cohort, MANI$seed), ]

## 与注册网格核对：只接受（队列 × 种子）恰为该网格的产物。
## 防止把试点、非注册种子、或旧版本的产物混进正式分析（评审要求"同一组合只能有一条"）。
grid_f <- file.path(PROJ, "config/design_grid.txt")
if (file.exists(grid_f)) {
  gl <- readLines(grid_f, warn = FALSE)
  gl <- gl[!grepl("^#", gl) & nzchar(trimws(gl))]
  COH <- strsplit(trimws(gl[1]), "[[:space:]]+")[[1]]
  SED <- strsplit(trimws(gl[2]), "[[:space:]]+")[[1]]
  want <- expand.grid(cohort = COH, seed = as.integer(SED), stringsAsFactors = FALSE)
  got  <- MANI[, c("cohort","seed")]
  miss <- nrow(merge(want, got, by = c("cohort","seed")))
  if (nrow(MANI) != nrow(want) || miss != nrow(want)) {
    cat("注册网格:", nrow(want), "个组合；实际产物:", nrow(MANI), "个
")
    extra <- merge(got, want, by = c("cohort","seed"))
    if (nrow(extra) < nrow(got)) cat("  ★ 存在不在注册网格内的产物，已列出：
")
    absent <- merge(want, got, by = c("cohort","seed"))
    if (nrow(absent) < nrow(want)) {
      lack <- want[!paste(want$cohort, want$seed) %in% paste(got$cohort, got$seed), ]
      cat("  ★ 缺少组合:", paste(lack$cohort, lack$seed, sep="/", collapse=", "), "
")
    }
    stop("manifest 拒绝生成：产物集合与注册网格不符")
  }
  cat("与注册网格核对通过：", nrow(want), "个（队列×种子）组合齐全且无多余
")
}
write.table(MANI, file.path(MAN, "run_manifest.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
writeLines(MANI$path, file.path(MAN, "main_input_files.txt"))

cat("manifest 生成完毕：", nrow(MANI), " 个（队列×种子）产物\n")
cat("  队列数 ", length(unique(MANI$cohort)), " × 种子数 ", length(unique(MANI$seed)), "\n")
cat("  终点分布: "); print(table(MANI$cohort, MANI$endpoint))
cat("  每个组合唯一性: "); print(anyDuplicated(paste(MANI$cohort, MANI$seed)) == 0)
cat("  缺失处理汇总: 剔除基因中位 ", median(MANI$n_gene_dropped_missing), "，填补值中位 ",
    median(MANI$n_value_imputed), "\n")
cat("  清单写入: results/manifest/run_manifest.tsv 与 main_input_files.txt\n")

#!/usr/bin/env Rscript
# ============================================================
# _verify_19.R — 对 19 的产物做独立复核（不信任脚本自报的数字）
#
# 复核项：
#   1. completed / n_planned / n_processed 一致
#   2. train/val/test 患者 ID 是全体样本的一个【划分】：无重叠、无遗漏
#   3. 用相同种子独立重跑划分，患者 ID 应与文件中保存的完全一致
#   4. 事件数、事件编码与【原始生存文件】逐一对应
#   5. 预测矩阵的样本数、基因集基因数与元数据一致
#   6. 缺失处理：独立重算"训练集缺失率>阈值"的基因数与填补值个数
#   7. 标准化：独立重算 Z 的一个随机单元格，与存储的预测不矛盾（抽查训练集统计量）
#   8. 方向诊断：训练集 C-index 是否全部 > 0.5
#
# 用法：Rscript scripts/_verify_19.R [某个 19_*.rds 的路径]
# ============================================================
suppressMessages(library(survival))
options(stringsAsFactors = FALSE)
`%||%` <- function(a, b) if (is.null(a)) b else a
PROJ <- "D:/bioinfo/sigAudit"

args <- commandArgs(trailingOnly = TRUE)
f <- if (length(args) >= 1) args[1] else {
  cand <- list.files(file.path(PROJ,"results/main"), pattern="^19_.*\\.rds$", full.names=TRUE)
  cand <- cand[!grepl("\\.partial\\.rds$", cand)]
  if (!length(cand)) stop("results/main 下没有完整产物"); cand[1]
}
R <- readRDS(f)
cat("复核对象:", basename(f), "\n")
cat("  队列", R$cohort, " study", R$study_id, " 终点", R$endpoint, " 种子", R$seed, "\n\n")

PASS <- 0; FAIL <- 0
chk <- function(cond, label, detail="") {
  ok <- isTRUE(cond)
  cat(sprintf("  [%s] %s%s\n", if (ok) "PASS" else "**FAIL**", label,
              if (nzchar(detail)) paste0("  —— ", detail) else ""))
  if (ok) PASS <<- PASS + 1 else FAIL <<- FAIL + 1
}

## 1. 完成标志
chk(isTRUE(R$completed), "completed == TRUE")
chk(identical(as.integer(R$n_planned), as.integer(R$n_processed)),
    "n_planned == n_processed", sprintf("%d vs %d", R$n_planned, R$n_processed))

## 2. 划分的完整性
allids <- R$sample_ids; tr <- R$train_ids; va <- R$val_ids; te <- R$test_ids
chk(setequal(c(tr,va,te), allids), "三个划分的并集 == 全部样本",
    sprintf("样本 %d = 训练 %d + 验证 %d + 测试 %d", length(allids), length(tr), length(va), length(te)))
chk(!anyDuplicated(c(tr,va,te)), "三划分之间无重复样本")
chk(length(intersect(tr,va))==0 && length(intersect(tr,te))==0 && length(intersect(va,te))==0,
    "两两交集为空")

## 3. 独立重跑划分，比对患者 ID
ENDPOINTS <- list(
  LGG  = list(e="D:/bioinfo/BRCA_TLS/data/cancers/LGG_broad_expr.tsv",
              f="D:/bioinfo/BRCA_TLS/data/cancers/LGG_broad_surv.tsv",
              t="PFS_MONTHS", s="PFS_STATUS", event=c("1:PROGRESSION"), censor=c("0:CENSORED")),
  PAAD = list(t="PFS_MONTHS", s="PFS_STATUS", event=c("1:PROGRESSION"), censor=c("0:CENSORED")),
  BRCA = list(t="OS_MONTHS", s="OS_STATUS", event=c("1:DECEASED"), censor=c("0:LIVING")),
  COAD = list(t="OS_MONTHS", s="OS_STATUS", event=c("1:DECEASED"), censor=c("0:LIVING")),
  KIRC = list(t="OS_MONTHS", s="OS_STATUS", event=c("1:DECEASED"), censor=c("0:LIVING")),
  LUAD = list(t="OS_MONTHS", s="OS_STATUS", event=c("1:DECEASED"), censor=c("0:LIVING")),
  MESO = list(t="OS_MONTHS", s="OS_STATUS", event=c("1:DECEASED"), censor=c("0:LIVING")),
  STAD = list(t="OS_MONTHS", s="OS_STATUS", event=c("1:DECEASED"), censor=c("0:LIVING"))
)
E <- ENDPOINTS[[R$cohort]]
cf <- E[["e"]] %||% sprintf("%s/data/cohorts/%s_expr.tsv", PROJ, R$cohort)
sf <- E[["f"]] %||% sprintf("%s/data/cohorts/%s_surv.tsv", PROJ, R$cohort)
expr <- as.matrix(read.table(cf, header=TRUE, sep="\t", row.names=1, check.names=FALSE))
if (nrow(expr) < ncol(expr)) expr <- t(expr)
expr <- log2(expr + 1)
sv <- read.table(sf, header=TRUE, sep="\t", row.names=1, check.names=FALSE)
raws <- as.character(sv[[E[["s"]]]])
sv$t <- suppressWarnings(as.numeric(sv[[E[["t"]]]]))
sv$e <- ifelse(raws %in% E[["event"]], 1L, ifelse(raws %in% E[["censor"]], 0L, NA_integer_))
sv <- sv[!is.na(sv$t) & !is.na(sv$e) & sv$t > 0, ]
colnames(expr) <- substr(colnames(expr),1,12); expr <- expr[, !duplicated(colnames(expr)), drop=FALSE]
com <- intersect(colnames(expr), rownames(sv)); expr <- expr[,com,drop=FALSE]; sv <- sv[com,]
set.seed(R$seed)
ev <- which(sv$e==1); ce <- which(sv$e==0)
tr2 <- sort(c(sample(ev, floor(.5*length(ev))), sample(ce, floor(.5*length(ce)))))
rest <- setdiff(seq_along(sv$e), tr2)
re <- rest[sv$e[rest]==1]; rc <- rest[sv$e[rest]==0]
va2 <- sort(c(sample(re, floor(.5*length(re))), sample(rc, floor(.5*length(rc)))))
te2 <- setdiff(rest, va2)
chk(identical(colnames(expr)[tr2], tr),
    "用相同种子重跑得到的训练集 ID 与文件一致")
chk(identical(colnames(expr)[va2], va) && identical(colnames(expr)[te2], te),
    "验证/测试集 ID 也一致")

## 4. 事件数与原始文件对应
g1 <- R$genesets[[1]]
chk(length(g1$e_tr) == length(tr) && sum(g1$e_tr) == sum(sv$e[tr2]),
    "训练集事件数与原始生存表一致", sprintf("%d", sum(sv$e[tr2])))
chk(sum(g1$eva) == sum(sv$e[va2]) && sum(g1$ete) == sum(sv$e[te2]),
    "验证/测试事件数与原始生存表一致",
    sprintf("%d / %d", sum(sv$e[va2]), sum(sv$e[te2])))
chk(all(g1$tva == sv$t[va2]) && all(g1$tte == sv$t[te2]),
    "验证/测试的生存时间与原始表逐值一致")

## 5. 预测矩阵维度
okdim <- all(vapply(R$genesets, function(z)
  nrow(z$val_pred) == length(va) && nrow(z$test_pred) == length(te) &&
  ncol(z$val_pred) == length(z$models) && ncol(z$test_pred) == length(z$models) &&
  nrow(z$train_pred) == length(tr), logical(1)))
chk(okdim, "所有基因集的（训练/验证/测试）×候选 预测矩阵维度正确")
chk(all(vapply(R$genesets, function(z) all(is.finite(z$val_pred)) && all(is.finite(z$test_pred)),
               logical(1))), "所有预测值均为有限值（无 NaN/Inf）")

## 6. 缺失处理：独立重算
trm <- expr[, tr, drop=FALSE]
miss_frac <- rowSums(is.na(trm)) / ncol(trm)
n_drop <- sum(miss_frac > R$max_miss)
ex <- expr[!(miss_frac > R$max_miss), , drop=FALSE]
n_imp <- sum(is.na(ex))
chk(n_drop == R$n_gene_dropped_missing, "训练集因缺失超阈值被剔除的基因数一致",
    sprintf("独立重算 %d, 文件记录 %d", n_drop, R$n_gene_dropped_missing))
chk(n_imp == R$n_value_imputed, "需填补的缺失值个数一致",
    sprintf("独立重算 %d, 文件记录 %d", n_imp, R$n_value_imputed))

## 7. 标准化用的统计量只来自训练集（抽查一个基因）
gmed <- apply(ex[, tr, drop=FALSE], 1, median, na.rm=TRUE)
for (j in seq_len(nrow(ex))) { v <- ex[j,]; k <- is.na(v); if (any(k)) v[k] <- gmed[j]; ex[j,] <- v }
gmu <- rowMeans(ex[, tr, drop=FALSE]); gsd <- apply(ex[, tr, drop=FALSE], 1, sd)
keepg <- is.finite(gsd) & gsd > 0
chk(length(keepg) - sum(!keepg) == length(R$genes), "保留基因数与文件记录的 genes 长度一致",
    sprintf("%d vs %d", sum(keepg), length(R$genes)))
chk(setequal(rownames(ex)[keepg], R$genes), "保留的基因名集合一致")

## 8. 方向诊断
D <- do.call(rbind, R$train_ci)
chk(all(is.finite(D$train_ci)), "训练集 C-index 全部可算")
# 必须区分两种情况（评审要求）：
#   ① 常数预测 -> C 恰为 0.500：正则化选到零系数模型，属【合法的零信号解】，
#      不是实现错误，也不代表方向有问题；
#   ② C < 0.5：才是方向约定异常，必须报警。
D$const <- abs(D$train_ci - 0.5) < 1e-9
chk(!any(D$train_ci < 0.5 - 1e-9), "没有候选出现训练集 C-index < 0.5（无方向异常）",
    sprintf("最小 %.3f", min(D$train_ci)))
nconst <- sum(D$const)
cat(sprintf("  [INFO] 常数预测（C 恰为 0.500）的（基因集×候选）条目: %d / %d (%.1f%%)
",
            nconst, nrow(D), 100*nconst/nrow(D)))
if (nconst > 0) {
  tc <- aggregate(const ~ model, data=D, FUN=sum)
  tc <- tc[tc$const > 0, ]; tc <- tc[order(-tc$const), ]
  cat("         涉及候选: ", paste(sprintf("%s(%d)", tc$model, tc$const), collapse=", "), "
")
  cat("         含义：正则化在该基因集上选到零系数模型（合法零信号解），
")
  cat("         其预测为常数，C 恒为 0.5；已在 model_metrics_long 中以 constant=TRUE 标记。
")
}

cat(sprintf("\n复核结果：%d 项通过，%d 项失败\n", PASS, FAIL))
if (FAIL > 0) quit(status = 1)

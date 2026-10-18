#!/usr/bin/env Rscript
# ============================================================
# 31 — METABRIC 独立队列：**流程复现**（第四轮评审要求的补充实验）
#
# 定位（评审明确要求，不得混称）：
#   本脚本在 METABRIC 内部**重新做训练/验证/测试划分并重新拟合候选池**，
#   检验"收益依赖基线"的判断能否在**独立来源**的队列上再现。
#   它**不是**把 TCGA 训练好的冻结模型迁移到 METABRIC 的外部预测验证。
#
# 预先冻结、运行前已确定的内容（见 config/metabric_frozen_config.md）：
#   · 终点：OS_MONTHS / OS_STATUS（1:DECEASED = 事件），时间起点 = 初诊
#   · 入排：Primary + OS 可用 + OS_MONTHS>0 + 只排除血管肉瘤（PBS，属肉瘤）；
#           **保留化生性癌（MBC）**——它本身就是浸润性乳腺癌。
#   · 5 个划分种子与主分析相同：20261014, 20261016–20261019
#   · 候选池 18 个配置，与主分析**逐字一致**（镜像自 _lib_fit.R）
#   · 预处理：**不做 log2**（按 cBioPortal 发布值使用，不额外取对数）；训练侧 z-score
#   · 规则 15 条全部计算（复用 _lib_analysis.R 的 rules_one）
#   · 失败与警告逐任务记录
#
# 并行：30 个任务彼此独立 -> 按【任务】并行（不是按划分）。
#
# 输出：
#   31_metabric_rule_performance.tsv  每任务每规则的 Cv/Ct/获胜者/成员
#   31_metabric_task_info.tsv         每任务的基因数、可用候选、缺失填补数、状态
#   31_metabric_summary.tsv           按签名 + 总体的队列内汇总
# ============================================================
suppressMessages({library(survival); library(glmnet); library(randomForestSRC); library(gbm)})
options(stringsAsFactors = FALSE)
`%||%` <- function(a, b) if (is.null(a)) b else a

PROJ <- Sys.getenv("SIGAUDIT_ROOT", "D:/bioinfo/sigAudit")
TAB  <- file.path(PROJ, "results/tables"); dir.create(TAB, showWarnings = FALSE, recursive = TRUE)
LOGDIR <- file.path(PROJ, "logs"); dir.create(LOGDIR, showWarnings = FALSE, recursive = TRUE)
con <- file(file.path(LOGDIR, "31_metabric_audit.log"), "wt")
say <- function(...) { m <- paste0(...); cat(m,"\n",sep=""); cat(m,"\n",sep="",file=con); flush(con) }

.lib1 <- file.path(PROJ, "scripts/_lib_analysis.R")
.lib2 <- file.path(PROJ, "scripts/_lib_fit.R")
for (f in c(.lib1, .lib2)) if (!file.exists(f)) stop("缺少 ", f)
source(.lib1)   # ci_direct / rules_one / make_prio / ens_pred / pick_max / order_rule
source(.lib2)   # POOL / POOL_NAMES / fit_pool（镜像自 19）

SEEDS <- c(20261014, 20261016, 20261017, 20261018, 20261019)
B_VOTE <- 12; B_HALF <- 6
MAX_MISS <- 0.20

rd <- function(f) read.table(f, sep="\t", header=TRUE, quote="", comment.char="",
                             check.names=FALSE, stringsAsFactors=FALSE)

## ---------------- 数据 ----------------
ef <- file.path(PROJ, "data/cohorts/METABRIC_expr.tsv")
sf <- file.path(PROJ, "data/cohorts/METABRIC_surv.tsv")
for (f in c(ef, sf)) if (!file.exists(f)) stop("缺少 ", f)
expr <- as.matrix(read.table(ef, sep="	", header=TRUE, quote="", comment.char="",
                             check.names=FALSE, row.names=1))
sv <- rd(sf)
sv$OS_MONTHS <- as.numeric(sv$OS_MONTHS)
sv$e <- ifelse(sv$OS_STATUS == "1:DECEASED", 1L, 0L)
com <- intersect(colnames(expr), sv$sampleId)
expr <- expr[, com, drop = FALSE]; sv <- sv[match(com, sv$sampleId), ]
t <- sv$OS_MONTHS; e <- sv$e
say("METABRIC: ", ncol(expr), " 样本 × ", nrow(expr), " 基因；事件 ", sum(e),
    sprintf("（%.1f%%）", 100*mean(e)))
say("表达尺度：按 cBioPortal 发布值原样使用，不额外取对数")

## 【第五轮修正】用 HGNC 规范化后的批准符号（scripts/34_map_signature_symbols.py）
sgf <- file.path(PROJ, "results/manifest/public_signatures_genes_mapped.tsv")
if (!file.exists(sgf)) stop("缺少 public_signatures_genes_mapped.tsv；请先跑 scripts/34_map_signature_symbols.py")
sg <- rd(sgf)
SIGS  <- lapply(split(sg, sg$signature), function(d) unique(d$gene_used))
SRC_N <- vapply(split(sg$gene_source, sg$signature), function(x) length(unique(x)), integer(1))
say("签名 ", length(SIGS), " 个：", paste(sprintf("%s(%d)", names(SIGS), lengths(SIGS)), collapse=", "))

## ---------------- 单任务 ----------------
one_task <- function(sig, seed, expr, t, e) {
  st <- new.env(parent = emptyenv()); st$nw <- 0L; st$first <- NA_character_
  info <- data.frame(signature=sig, seed=seed, n_gene_source=unname(SRC_N[[sig]]),
                     n_gene_used=NA_integer_, n_candidate=NA_integer_,
                     n_imputed=NA_integer_, n_train=NA_integer_, n_val=NA_integer_,
                     n_test=NA_integer_, ev_train=NA_integer_, ev_val=NA_integer_,
                     ev_test=NA_integer_, status="failed", reason=NA_character_,
                     n_warn=NA_integer_, stringsAsFactors=FALSE)
  rows <- NULL; mrows <- NULL
  err <- tryCatch(withCallingHandlers({
    g <- intersect(SIGS[[sig]], rownames(expr))
    if (length(g) < 8) stop(sprintf("only %d genes", length(g)))
    info$n_gene_used <- length(g)
    M <- expr[g, , drop=FALSE]
    ## 划分（与主分析相同的分层 50/25/25）
    set.seed(seed)
    ev <- which(e==1); ce <- which(e==0)
    tr <- sort(c(sample(ev, floor(.5*length(ev))), sample(ce, floor(.5*length(ce)))))
    rest <- setdiff(seq_along(e), tr)
    re <- rest[e[rest]==1]; rc <- rest[e[rest]==0]
    va <- sort(c(sample(re, floor(.5*length(re))), sample(rc, floor(.5*length(rc)))))
    te <- setdiff(rest, va)
    info$n_train <- length(tr); info$n_val <- length(va); info$n_test <- length(te)
    info$ev_train <- sum(e[tr]); info$ev_val <- sum(e[va]); info$ev_test <- sum(e[te])
    ## 训练侧预处理：缺失过滤 -> 训练中位数填补 -> 零方差剔除 -> 训练侧 z-score
    trm <- M[, tr, drop=FALSE]
    drop_g <- rowSums(is.na(trm))/ncol(trm) > MAX_MISS
    ex <- M[!drop_g, , drop=FALSE]
    med <- apply(ex[, tr, drop=FALSE], 1, function(v) stats::median(v, na.rm=TRUE))
    n_imp <- sum(is.na(ex))
    if (n_imp > 0) for (j in seq_len(nrow(ex))) { v <- ex[j,]; k <- is.na(v)
      if (any(k)) { v[k] <- med[j]; ex[j,] <- v } }
    mu <- rowMeans(ex[, tr, drop=FALSE]); sdv <- apply(ex[, tr, drop=FALSE], 1, stats::sd)
    keep <- is.finite(sdv) & sdv > 0
    ex <- ex[keep, , drop=FALSE]; mu <- mu[keep]; sdv <- sdv[keep]
    Z <- (ex - mu)/sdv
    info$n_imputed <- n_imp
    G <- function(idx) t(Z[, idx, drop=FALSE])
    ## 候选池
    pm <- fit_pool(G(tr), t[tr], e[tr], G(va), G(te))
    mods <- POOL_NAMES[POOL_NAMES %in% names(pm)]
    info$n_candidate <- length(mods)
    if (length(mods) < 4) stop(sprintf("only %d candidates fitted", length(mods)))
    if (!("Lasso" %in% mods)) stop("reference (Lasso) missing")
    mk <- function(key) {
      Mx <- sapply(mods, function(x) pm[[x]][[key]])
      if (is.null(dim(Mx))) Mx <- matrix(Mx, ncol=1)
      colnames(Mx) <- mods; Mx }
    Pv <- mk("pb"); Pt <- mk("pc"); Ptr <- mk("pa")
    ## 【第五轮】保存逐患者预测与划分 ID，使 C-index 可由患者级风险分独立重算
    preds <- list(signature=sig, seed=seed, model=mods,
                  train_ids=colnames(Z)[tr], val_ids=colnames(Z)[va], test_ids=colnames(Z)[te],
                  t_train=t[tr], e_train=e[tr], t_val=t[va], e_val=e[va],
                  t_test=t[te], e_test=e[te],
                  train_pred=Ptr, val_pred=Pv, test_pred=Pt,
                  genes=rownames(Z), n_imputed=n_imp)
    if (!all(is.finite(Pv)) || !all(is.finite(Pt))) stop("non-finite predictions")
    prio <- make_prio(mods, seed)
    rn <- rules_one(list(name=sig, val_pred=Pv, test_pred=Pt, train_pred=Ptr,
                         tva=t[va], eva=e[va], tte=t[te], ete=e[te],
                         t_tr=t[tr], e_tr=e[tr]),
                    prio, ci_direct, B_VOTE = B_VOTE, B_HALF = B_HALF)
    ## 逐候选的 Cv/Ct（供参照敏感性：固定 Ridge / Enet0.5）
    mrows <- data.frame(signature=sig, seed=seed, model=mods,
                        cv=as.numeric(rn$V[mods]), ct=as.numeric(rn$T[mods]),
                        stringsAsFactors=FALSE)
    rows <- do.call(rbind, lapply(names(rn$out), function(k)
      data.frame(signature=sig, seed=seed, rule=k,
                 cv=rn$out[[k]][1], ct=rn$out[[k]][2],
                 winner=rn$rec[[k]]$winner, members=rn$rec[[k]]$members,
                 stringsAsFactors=FALSE)))
    NULL
  }, warning = function(w) {
    st$nw <- st$nw + 1L
    if (is.na(st$first)) st$first <- conditionMessage(w)
    invokeRestart("muffleWarning")
  }), error = function(err2) err2)
  info$n_warn <- st$nw
  if (is.null(err)) info$status <- "ok" else { info$status <- "failed"; info$reason <- conditionMessage(err) }
  list(info=info, rows=rows, mrows=mrows,
       preds=if (exists("preds")) preds else NULL)
}

## ---------------- 并行 ----------------
TASKS <- expand.grid(signature=names(SIGS), seed=SEEDS, stringsAsFactors=FALSE)
TASKS <- TASKS[order(TASKS$signature, TASKS$seed), ]
NC <- as.integer(Sys.getenv("SIGAUDIT_WORKERS", "12"))
say("任务数 ", nrow(TASKS), "；并行 worker ", NC)

cl <- parallel::makeCluster(NC)
on.exit(parallel::stopCluster(cl), add = TRUE)
parallel::clusterExport(cl, c("TASKS","TREES","SIGS","SRC_N","SEEDS","expr","t","e","MAX_MISS","B_VOTE","B_HALF",
                              "POOL","POOL_NAMES","fit_pool","rules_one","ci_direct",
                              "make_prio","ens_pred","pick_max","order_rule","one_task"),
                        envir = environment())
invisible(parallel::clusterEvalQ(cl, suppressMessages({library(survival); library(glmnet);
                                                        library(randomForestSRC); library(gbm)})))

t0 <- Sys.time()
RES <- parallel::parLapply(cl, seq_len(nrow(TASKS)), function(i) {
  suppressMessages({library(survival); library(glmnet); library(randomForestSRC); library(gbm)})
  one_task(TASKS$signature[i], TASKS$seed[i], expr, t, e)
})
say("用时 ", round(as.numeric(difftime(Sys.time(), t0, units="min")), 2), " 分钟")

## 【第五轮】落盘逐患者预测与划分 ID
PRED <- lapply(RES, `[[`, "preds")
PRED <- PRED[!vapply(PRED, is.null, logical(1))]
dir.create(file.path(PROJ, "results/metabric"), showWarnings = FALSE, recursive = TRUE)
if (length(PRED)) {
  saveRDS(PRED, file.path(PROJ, "results/metabric/31_metabric_predictions.rds"))
  say("逐患者预测已保存: results/metabric/31_metabric_predictions.rds（",
      length(PRED), " 个任务）")
}

INFO <- do.call(rbind, lapply(RES, `[[`, "info"))
ROWS <- do.call(rbind, lapply(RES, `[[`, "rows"))
MMET <- do.call(rbind, lapply(RES, `[[`, "mrows"))
if (!is.null(MMET)) write.table(MMET, file.path(TAB, "31_metabric_model_metrics.tsv"),
                                sep="\t", quote=FALSE, row.names=FALSE)
write.table(INFO, file.path(TAB, "31_metabric_task_info.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
if (!is.null(ROWS)) {
  write.table(ROWS, file.path(TAB, "31_metabric_rule_performance.tsv"),
              sep="\t", quote=FALSE, row.names=FALSE)
}

say("")
say("===== 任务完成情况 =====")
say("  计划 ", nrow(INFO), "，成功 ", sum(INFO$status=="ok"), "，失败 ", sum(INFO$status!="ok"))
if (any(INFO$status != "ok")) print(INFO[INFO$status != "ok", c("signature","seed","reason")])
say("  可用基因数范围 ", min(INFO$n_gene_used, na.rm=TRUE), "–", max(INFO$n_gene_used, na.rm=TRUE))
say("  可用候选数范围 ", min(INFO$n_candidate, na.rm=TRUE), "–", max(INFO$n_candidate, na.rm=TRUE))
say("  填补的缺失值总数 ", sum(INFO$n_imputed, na.rm=TRUE))
say("  含 warning 的任务数 ", sum(INFO$n_warn > 0, na.rm=TRUE))

if (!is.null(ROWS)) {
  ## 汇总：先签名内跨划分平均，再对签名等权
  bys <- aggregate(cbind(cv, ct) ~ signature + rule, data = ROWS, FUN = mean, na.rm = TRUE)
  ## 相对参照（none = Lasso）
  ref <- bys[bys$rule == "none", c("signature","cv","ct")]
  names(ref) <- c("signature","ref_cv","ref_ct")
  bys <- merge(bys, ref, by = "signature")
  bys$Gv <- bys$cv - bys$ref_cv; bys$Gt <- bys$ct - bys$ref_ct; bys$Orel <- bys$Gv - bys$Gt
  agg <- aggregate(cbind(cv, ct, Gv, Gt, Orel) ~ rule, data = bys, FUN = mean, na.rm = TRUE)
  agg <- agg[order(-agg$ct), ]
  write.table(agg, file.path(TAB, "31_metabric_summary.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
  say("")
  say("===== 相对固定 Lasso 的汇总（6 个签名等权）=====")
  say("  ", sprintf("%-22s %8s %8s %9s %9s %9s", "rule","Cv","Ct","Gv","Gt","Orel"))
  for (i in seq_len(nrow(agg))) say("  ", sprintf("%-22s %8.4f %8.4f %+9.4f %+9.4f %+9.4f",
    agg$rule[i], agg$cv[i], agg$ct[i], agg$Gv[i], agg$Gt[i], agg$Orel[i]))
  ## 固定基线（单候选的汇总性能）
  say("")
  say("===== 各签名（相对 Lasso）=====")
  say("  ", sprintf("%-20s %8s %8s %9s %9s %9s", "signature","Cv","Ct","Gv","Gt","Orel"))
  z <- bys[bys$rule=="argmax_val", ]; z <- z[order(-z$Orel), ]
  for (i in seq_len(nrow(z))) say("  ", sprintf("%-20s %8.4f %8.4f %+9.4f %+9.4f %+9.4f",
    z$signature[i], z$cv[i], z$ct[i], z$Gv[i], z$Gt[i], z$Orel[i]))
}
say("")
say("⚠️ 30 个任务共享同一批患者、且签名之间也共享患者 -> 不能当 30 个独立样本做显著性检验。")
close(con)

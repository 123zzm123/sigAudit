#!/usr/bin/env Rscript
# ============================================================
# 14 — 划分稳健性：主结果全部来自同一固定划分（seed 20261014），
#      审稿人必然要问"换一个划分还成立吗"。
#
# run_v6c.sh 用 seed 20261016 重跑同一批基因集（只换划分），输出 09b_*.tsv。
# 本脚本比较两套划分下：
#   1) 每个候选的进入效应排名（Spearman / Pearson 相关）
#   2) "危害最大的候选"是否还是同一个
#   3) 家族对比（树集 − 线性，同 k）的符号是否翻转
#   4) 剂量曲线形状
# 若排名高度一致 -> 进入效应是候选的稳定属性；
# 若不一致        -> 结论必须限定为"依划分而定"，这本身也是重要发现。
# ============================================================
suppressMessages(source("D:/bioinfo/sigAudit/R/audit_modules.R"))
options(stringsAsFactors = FALSE)
PROJ <- "D:/bioinfo/sigAudit"

POOL_NAMES <- c("Lasso","Ridge","Enet0.9","Enet0.7","Enet0.5","Enet0.3","Enet0.1",
                "Ridge_1se","UniCox+Lasso","LassoTop100","LassoTop50","LassoTop25",
                "LassoTop10","Enet0.5Top25","RSF100","GBM100","GBM250","RSF250")
TREES <- c("RSF100","GBM100","GBM250","RSF250")
FAM <- setNames(ifelse(POOL_NAMES %in% TREES, "tree", "linear"), POOL_NAMES)

LOGF <- file.path(PROJ,"logs","14_robustness.log"); con <- file(LOGF,"wt")
say <- function(...) { m <- paste0(...); cat(m,"\n",sep=""); cat(m,"\n",sep="",file=con); flush(con) }

load_mat <- function(f) {
  D <- tryCatch(read.table(f, header=TRUE, sep="\t", check.names=FALSE), error=function(e) NULL)
  if (is.null(D) || !nrow(D)) return(NULL)
  D <- D[D$model %in% POOL_NAMES, ]; if (!nrow(D)) return(NULL)
  all_gs <- unique(D$geneset); n_all <- length(all_gs)
  cnt <- tapply(D$geneset, D$model, function(g) length(unique(g)))
  pnames <- POOL_NAMES[POOL_NAMES %in% names(cnt)[cnt/n_all >= 0.90]]
  if (!("Lasso" %in% pnames) || length(pnames) < 5) return(NULL)
  sub <- D[D$model %in% pnames, ]
  tab <- table(sub$geneset); keep <- names(tab)[tab == length(pnames)]
  if (length(keep) < 30) return(NULL)
  mk <- function(v) { M <- tapply(sub[[v]], list(sub$geneset, sub$model), function(x) x[1])
    M <- M[keep, pnames, drop=FALSE]
    matrix(as.numeric(M), nrow=nrow(M), dimnames=dimnames(M)) }
  V <- mk("val_ci"); Tt <- mk("test_ci")
  ok <- complete.cases(V) & complete.cases(Tt)
  list(V = V[ok,,drop=FALSE], Tt = Tt[ok,,drop=FALSE])
}

fam_contrast <- function(M, kk) {   # 同 k、均不含参照
  V <- M$V; Tt <- M$Tt
  ref_v <- V[,"Lasso"]; ref_t <- Tt[,"Lasso"]
  cand <- setdiff(colnames(V), "Lasso")
  lin <- setdiff(cand, TREES); tre <- intersect(cand, TREES)
  k <- min(kk, length(lin), length(tre)); if (k < 2) return(NA_real_)
  b <- function(cols) { sub <- cols[seq_len(k)]
    idx <- max.col(V[, sub, drop=FALSE], ties.method="first")
    vapply(seq_len(nrow(V)), function(i) { m <- sub[idx[i]]
      (V[i,m]-ref_v[i]) - (Tt[i,m]-ref_t[i]) }, numeric(1)) }
  mean(b(tre) - b(lin))
}

cohs <- unique(sub("^09b?_(.+)_[0-9]+\\.tsv$","\\1",
                   list.files(file.path(PROJ,"results/tables"), pattern="^09b?_")))
RES <- list()
for (cn in cohs) {
  f1 <- file.path(PROJ,"results/tables", sprintf("09_%s_100.tsv", cn))
  f2 <- file.path(PROJ,"results/tables", sprintf("09b_%s_100.tsv", cn))
  if (!file.exists(f1) || !file.exists(f2)) { say("[skip] ", cn, " 缺一套划分"); next }
  M1 <- load_mat(f1); M2 <- load_mat(f2)
  if (is.null(M1) || is.null(M2)) { say("[skip] ", cn, " 载入失败"); next }
  r1 <- audit_selection(M1$V, M1$Tt, reference="Lasso")
  r2 <- audit_selection(M2$V, M2$Tt, reference="Lasso")
  a <- merge(r1$entry_summary[, c("target","mean","optimism")],
             r2$entry_summary[, c("target","mean","optimism")],
             by="target", suffixes=c("_s1","_s2"))
  sp <- suppressWarnings(cor(a$mean_s1, a$mean_s2, method="spearman"))
  pe <- cor(a$mean_s1, a$mean_s2)
  top1 <- a$target[which.max(a$mean_s1)]; top2 <- a$target[which.max(a$mean_s2)]
  fc1 <- fam_contrast(M1, 4); fc2 <- fam_contrast(M2, 4)
  RES[[cn]] <- data.frame(cohort=cn, n_gs1=nrow(M1$V), n_gs2=nrow(M2$V),
                          spearman=sp, pearson=pe, top_s1=top1, top_s2=top2,
                          same_top=(top1==top2), fam1=fc1, fam2=fc2,
                          same_fam_sign=(sign(fc1)==sign(fc2)))
  say(""); say("===== ", cn, " =====")
  say(sprintf("  基因集数: 主划分 %d / 第二划分 %d", nrow(M1$V), nrow(M2$V)))
  say(sprintf("  进入效应排名一致性: Spearman %.3f  Pearson %.3f", sp, pe))
  say(sprintf("  危害最大候选: 主 %s / 第二 %s  -> %s", top1, top2,
              if (top1==top2) "相同" else "★不同★"))
  say(sprintf("  同样 k=4 家族对比: 主 %+.4f / 第二 %+.4f  -> %s",
              fc1, fc2, if (sign(fc1)==sign(fc2)) "符号一致" else "★符号翻转★"))
  say("  进入效应对照（按主划分排序）:")
  a <- a[order(-a$mean_s1), ]
  for (i in seq_len(nrow(a))) say(sprintf("    %-14s 主 %+.4f   第二 %+.4f",
                                          a$target[i], a$mean_s1[i], a$mean_s2[i]))
}

if (length(RES)) {
  R <- do.call(rbind, RES)
  write.table(R, file.path(PROJ,"results/tables","14_robustness.tsv"),
              sep="\t", quote=FALSE, row.names=FALSE)
  say(""); say("========== 汇总 ==========")
  say(sprintf("  队列数 %d；进入效应排名 Spearman 中位 %.3f（范围 %.3f~%.3f）",
              nrow(R), median(R$spearman), min(R$spearman), max(R$spearman)))
  say(sprintf("  危害最大候选相同的队列: %d / %d", sum(R$same_top), nrow(R)))
  say(sprintf("  家族对比符号一致的队列: %d / %d", sum(R$same_fam_sign), nrow(R)))
  # 队列数太少时 cor.test 会直接报错（"not enough finite observations"）。
  # 这在实际运行时可能发生——例如某个队列静默失败、只剩一两个队列。
  if (nrow(R) >= 4) {
    ct <- suppressWarnings(cor.test(R$spearman, R$n_gs1, method="spearman"))
    say(sprintf("  排名一致性与队列大小（基因集数）无关? Spearman %.3f p=%.2f",
                ct$estimate, ct$p.value))
  } else {
    say(sprintf("  队列数只有 %d 个，跳过「排名一致性 vs 队列大小」的相关检验", nrow(R)))
  }
}
close(con)

#!/usr/bin/env Rscript
# ============================================================
# 11 — 哪个候选层面的指标能预测"该候选对偏倚的贡献"？
#
# 背景：10_analysis 显示"进入效应 vs 候选自身乐观度"的相关在各队列间不稳定
#       （LGG +0.68 / MESO +0.63 / STAD +0.54 / KIRC −0.53）。
#       说明"边际平均乐观度"不是正确的解释变量。这里系统地比较一组候选，
#       并严格区分【验证集可观测】与【不可观测】两类：
#
#   可观测（实践中能从训练+验证算出的）：
#     mean_dval   验证集上相对参照的平均增益
#     sd_dval     验证集增益在基因集间的离散度  ← 决定"上尾"能被抬多高
#     winrate     在全池中赢得验证集 argmax 的频率
#   不可观测（需要测试集才有的）：
#     mean_dtest  测试集平均增益
#     cor_vt      验证增益与测试增益的跨基因集相关  ← 决定 val 优势能否传导到 test
#     optimism    mean_dval − mean_dtest
#     tail_opt    top-quartile 基因集上的 (dval − dtest)
#
# 关键问题：进入效应能不能由【可观测】量预测？
#   能 -> 可以在验证集阶段就给出可操作的预警
#   不能 -> 只能靠估计量事后审计（本身就是重要结论）
#
# 复用 R/audit_modules.R 的 audit_selection()，避免与 10 的独立实现重复。
# ============================================================
suppressMessages(source("D:/bioinfo/sigAudit/R/audit_modules.R"))
options(stringsAsFactors = FALSE)
PROJ <- "D:/bioinfo/sigAudit"

POOL_NAMES <- c("Lasso","Ridge","Enet0.9","Enet0.7","Enet0.5","Enet0.3","Enet0.1",
                "Ridge_1se","UniCox+Lasso","LassoTop100","LassoTop50","LassoTop25",
                "LassoTop10","Enet0.5Top25","RSF100","GBM100","GBM250","RSF250")
FAM <- setNames(ifelse(POOL_NAMES %in% c("RSF100","GBM100","GBM250","RSF250"),"tree","linear"), POOL_NAMES)

LOGF <- file.path(PROJ,"logs","11_predictors.log"); con <- file(LOGF,"wt")
say <- function(...) { m <- paste0(...); cat(m,"\n",sep=""); cat(m,"\n",sep="",file=con); flush(con) }

files <- list.files(file.path(PROJ,"results/tables"), pattern="^09_.*_[0-9]+\\.tsv$", full.names=TRUE)
FEAT <- list()
for (f in files) {
  coh <- sub("^09_(.+)_[0-9]+\\.tsv$","\\1", basename(f))
  D <- tryCatch(read.table(f, header=TRUE, sep="\t", check.names=FALSE), error=function(e) NULL)
  if (is.null(D) || !nrow(D)) { say("[skip] ", coh, " 读取失败"); next }
  D <- D[D$model %in% POOL_NAMES, ]
  if (!nrow(D)) next
  all_gs <- unique(D$geneset); n_all <- length(all_gs)
  cnt <- tapply(D$geneset, D$model, function(g) length(unique(g)))
  pnames <- POOL_NAMES[POOL_NAMES %in% names(cnt)[cnt/n_all >= 0.90]]
  if (!("Lasso" %in% pnames) || length(pnames) < 5) { say("[skip] ", coh, " 候选不足"); next }
  sub <- D[D$model %in% pnames, ]
  tab <- table(sub$geneset); keep <- names(tab)[tab == length(pnames)]
  if (length(keep) < 30) { say("[skip] ", coh, " 完整基因集 <30（可能是增量写入中）"); next }
  mk <- function(v) { M <- tapply(sub[[v]], list(sub$geneset, sub$model), function(x) x[1])
    M <- M[keep, pnames, drop=FALSE]
    matrix(as.numeric(M), nrow=nrow(M), dimnames=dimnames(M)) }
  V <- mk("val_ci"); Tt <- mk("test_ci")
  ok <- complete.cases(V) & complete.cases(Tt)
  V <- V[ok,,drop=FALSE]; Tt <- Tt[ok,,drop=FALSE]
  if (nrow(V) < 30) { say("[skip] ", coh, " 无缺失基因集 <30"); next }
  say(""); say("===== ", coh, " : ", nrow(V), " 基因集 × ", ncol(V), " 候选 =====")

  r <- audit_selection(V, Tt, reference = "Lasso")
  ag <- r$entry_summary                       # target / mean / sd_over_pos / optimism
  dv <- sweep(V, 1, V[,"Lasso"], "-")         # 相对参照的验证增益
  dt <- sweep(Tt, 1, Tt[,"Lasso"], "-")
  win <- apply(V, 1, function(z) colnames(V)[which.max(z)])

  res <- do.call(rbind, lapply(setdiff(pnames,"Lasso"), function(m) {
    a <- dv[,m]; b <- dt[,m]
    q <- a >= quantile(a, .75, na.rm=TRUE)
    data.frame(cohort=coh, target=m, family=FAM[m], n_gs=nrow(V),
               entry_mean = ag$mean[match(m, ag$target)],
               sd_over_pos= ag$sd_over_pos[match(m, ag$target)],
               mean_dval=mean(a, na.rm=TRUE), sd_dval=sd(a, na.rm=TRUE),
               winrate  =mean(win==m, na.rm=TRUE),
               mean_dtest=mean(b, na.rm=TRUE),
               cor_vt   =suppressWarnings(cor(a,b,use="complete.obs")),
               optimism =mean(a-b, na.rm=TRUE),
               tail_opt =mean((a-b)[q], na.rm=TRUE))
  }))
  FEAT[[coh]] <- res
  # 注意：res 是数据框，直接喂给 sprintf 会被向量化，17 行全挤在一行打印。
  # 必须逐行循环。
  for (i in seq_len(nrow(res))) say(sprintf(
    "  %-14s %-7s 进入效应 %+.4f | 可观测: sd_dval %.4f winrate %.2f | 不可观测: cor_vt %+.2f 乐观度 %+.4f",
    res$target[i], res$family[i], res$entry_mean[i], res$sd_dval[i], res$winrate[i],
    res$cor_vt[i], res$optimism[i]))
}

if (!length(FEAT)) stop("没有可用队列")
FT <- do.call(rbind, FEAT)
write.table(FT, file.path(PROJ,"results/tables","11_candidate_features.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)

say(""); say("========== 哪个候选层面指标能预测『进入效应』==========")
say(sprintf("  合并 %d 条（队列×候选），涉及 %d 个队列", nrow(FT), length(FEAT)))
obsv <- c("mean_dval","sd_dval","winrate")
unob <- c("mean_dtest","cor_vt","optimism","tail_opt")
say("  --- 验证集【可观测】的指标 ---")
for (v in obsv) { ct <- cor.test(FT[[v]], FT$entry_mean)
  say(sprintf("    %-12s vs 进入效应: r = %+.3f  p = %.2e", v, ct$estimate, ct$p.value)) }
say("  --- 【不可观测】的指标（需测试集）---")
for (v in unob) { ct <- cor.test(FT[[v]], FT$entry_mean)
  say(sprintf("    %-12s vs 进入效应: r = %+.3f  p = %.2e", v, ct$estimate, ct$p.value)) }

say(""); say("  --- 分队列看：最优的单一预测因子 ---")
for (cn in unique(FT$cohort)) {
  s <- FT[FT$cohort==cn, ]; best <- NA; bv <- 0
  for (v in c(obsv, unob)) { rr <- suppressWarnings(cor(s[[v]], s$entry_mean)); if (is.finite(rr) && abs(rr) > abs(bv)) { bv <- rr; best <- v } }
  say(sprintf("    %-6s 最强单一相关: %-12s r = %+.3f  (n=%d)", cn, best, bv, nrow(s)))
}

# ---- 关键稳健性：合并相关有多少来自【队列间】差异？----
# 把每个变量在每个队列内中心化后重算相关 = 纯队列内（within-cohort）相关。
# 若队列内相关远弱于合并相关，则"进入效应 ∝ 乐观度"只是队列间现象，
# 不能用来预测【同一队列内】哪个候选危险。
say(""); say("  --- 队列内中心化后的相关（剔除队列间差异）---")
cen <- function(v) v - ave(v, FT$cohort, FUN = function(z) mean(z, na.rm=TRUE))
ec <- cen(FT$entry_mean)
say("  可观测：")
for (v in obsv) { ct <- suppressWarnings(cor.test(cen(FT[[v]]), ec))
  say(sprintf("    %-12s r_within = %+.3f  p = %.2e", v, ct$estimate, ct$p.value)) }
say("  不可观测：")
for (v in unob) { ct <- suppressWarnings(cor.test(cen(FT[[v]]), ec))
  say(sprintf("    %-12s r_within = %+.3f  p = %.2e", v, ct$estimate, ct$p.value)) }
say("  队列内中心化后的多元回归（仅可观测指标）：")
fitw <- lm(as.formula(paste("entry_mean ~ 0 +", paste(sprintf("I(%s - ave(%s, cohort, FUN=function(z) mean(z, na.rm=TRUE)))", obsv, obsv), collapse="+"))), data=FT)
say(sprintf("    R^2 = %.3f   n = %d", summary(fitw)$r.squared, nrow(FT)))

say(""); say("  --- 多元：用【全部可观测】指标预测进入效应 ---")
fml <- as.formula(paste("entry_mean ~", paste(obsv, collapse="+")))
fit <- lm(fml, data = FT)
say(sprintf("    R^2 = %.3f   校正 R^2 = %.3f   n = %d",
            summary(fit)$r.squared, summary(fit)$adj.r.squared, nrow(FT)))
for (i in seq_len(nrow(summary(fit)$coefficients)))
  say(sprintf("      %-12s 系数 %+.4f  p = %.3f", rownames(summary(fit)$coefficients)[i],
              summary(fit)$coefficients[i,1], summary(fit)$coefficients[i,4]))
say("  --- 多元：加上不可观测指标 ---")
fml2 <- as.formula(paste("entry_mean ~", paste(c(obsv,unob), collapse="+")))
fit2 <- lm(fml2, data = FT)
say(sprintf("    R^2 = %.3f   校正 R^2 = %.3f", summary(fit2)$r.squared, summary(fit2)$adj.r.squared))

say(""); say("========== 各队列的『最危险候选』与『最乐观候选』是否一致 ==========")
for (cn in unique(FT$cohort)) {
  s <- FT[FT$cohort==cn, ]
  w <- s$target[which.max(s$entry_mean)]
  o <- s$target[which.max(s$optimism)]
  say(sprintf("  %-6s 危害最大: %-14s (+%.4f)   最乐观: %-14s (%+.4f)   %s",
              cn, w, max(s$entry_mean), o, max(s$optimism),
              if (w == o) "一致" else "★不一致★"))
}
close(con)

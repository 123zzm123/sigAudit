#!/usr/bin/env Rscript
# ============================================================
# 18 — 双划分完整对照：把论文的每一条主张都在第二个划分上重算一遍
#
# 为什么需要：14_robustness.R 只对比了"进入效应排名"和"家族对比"。
#   但从 LGG/MESO 的早期结果看，这两项在换划分后都不稳定，而它们并不是
#   论文的核心主张。论文真正依赖的是下面三条，必须逐条在第二划分上验证：
#     A. 代价核算：报告增益 vs 真实增益（61% 虚高）
#     B. 家族效应方向（是否只是噪声）
#     C. 候选层面预测因子（mean_dval / sd_dval / winrate / tail_opt 的相关）
#   另加 D：进入效应排名的一致性（已知不稳，作为诚实报告的一部分）。
#
# 输入：09_<C>_100.tsv（主划分）与 09b_<C>_100.tsv（第二划分）
# 输出：results/tables/18_robustness_full.tsv
# ============================================================
suppressMessages(source("D:/bioinfo/sigAudit/R/audit_modules.R"))
options(stringsAsFactors = FALSE)
PROJ <- "D:/bioinfo/sigAudit"; TAB <- file.path(PROJ, "results/tables")

POOL_NAMES <- c("Lasso","Ridge","Enet0.9","Enet0.7","Enet0.5","Enet0.3","Enet0.1",
                "Ridge_1se","UniCox+Lasso","LassoTop100","LassoTop50","LassoTop25",
                "LassoTop10","Enet0.5Top25","RSF100","GBM100","GBM250","RSF250")
TREES <- c("RSF100","GBM100","GBM250","RSF250")

LOGF <- file.path(PROJ,"logs","18_robustness_full.log"); con <- file(LOGF,"wt")
say <- function(...) { m <- paste0(...); cat(m,"\n",sep=""); cat(m,"\n",sep="",file=con); flush(con) }

load_mat <- function(f) {
  D <- tryCatch(read.table(f, header=TRUE, sep="\t", check.names=FALSE), error=function(e) NULL)
  if (is.null(D) || !nrow(D) || !all(c("geneset","model","val_ci","test_ci") %in% names(D))) return(NULL)
  D <- D[D$model %in% POOL_NAMES, ]
  if (!nrow(D)) return(NULL)
  n_all <- length(unique(D$geneset))
  cnt <- tapply(D$geneset, D$model, function(g) length(unique(g)))
  pn <- POOL_NAMES[POOL_NAMES %in% names(cnt)[cnt/n_all >= 0.90]]
  if (!("Lasso" %in% pn) || length(pn) < 5) return(NULL)
  sub <- D[D$model %in% pn, ]
  tab <- table(sub$geneset); keep <- names(tab)[tab == length(pn)]
  if (length(keep) < 30) return(NULL)
  mk <- function(v) { M <- tapply(sub[[v]], list(sub$geneset, sub$model), function(x) x[1])
    M <- M[keep, pn, drop=FALSE]
    matrix(as.numeric(M), nrow=nrow(M), dimnames=dimnames(M)) }
  V <- mk("val_ci"); Tt <- mk("test_ci")
  ok <- complete.cases(V) & complete.cases(Tt)
  list(V = V[ok,,drop=FALSE], Tt = Tt[ok,,drop=FALSE], models = pn)
}

# ---- A. 代价核算（复用包内 audit_rules，已与脚本实现交叉验证过）----
stat_rules <- function(M) {
  fam <- setNames(ifelse(M$models %in% TREES, "tree", "linear"), M$models)
  r <- audit_rules(M$V, M$Tt, reference = "Lasso", family = fam)
  t <- r$table
  rep_gain  <- t$reported[t$rule=="argmax_val"] - t$reported[t$rule=="none"]
  true_gain <- t$gain_vs_ref[t$rule=="argmax_val"]
  # 全池总偏倚 = 剂量曲线末端的值。这是"估计量本身"的量（不涉及候选之间的
  # 相对排名），必须与"归因排名"分开看：前者可以很稳，后者可能纯噪声。
  sel <- audit_selection(M$V, M$Tt, reference = "Lasso")
  pool_bias <- sel$dose$bias[length(sel$dose$bias)]
  list(infl = 100*(rep_gain - true_gain)/rep_gain, rep_gain = rep_gain,
       true_gain = true_gain, pool_bias = pool_bias,
       topk = t$true_test[t$rule=="argmax_val"])
}

# ---- B. 家族对比（同样 k、均不含参照）----
stat_family <- function(M) {
  cand <- setdiff(M$models, "Lasso")
  lin <- setdiff(cand, TREES); tre <- intersect(cand, TREES)
  k <- min(4, length(lin), length(tre)); if (k < 2) return(NA_real_)
  b <- function(cols) { sub <- cols[seq_len(k)]
    idx <- max.col(M$V[, sub, drop=FALSE], ties.method="first")
    vapply(seq_len(nrow(M$V)), function(i) { m <- sub[idx[i]]
      (M$V[i,m]-M$V[i,"Lasso"]) - (M$Tt[i,m]-M$Tt[i,"Lasso"]) }, numeric(1)) }
  mean(b(tre) - b(lin))
}

# ---- C. 候选层面预测因子 ----
stat_feats <- function(M) {
  sel <- audit_selection(M$V, M$Tt, reference = "Lasso")
  ag <- sel$entry_summary
  dv <- sweep(M$V, 1, M$V[,"Lasso"], "-"); dt <- sweep(M$Tt, 1, M$Tt[,"Lasso"], "-")
  win <- apply(M$V, 1, function(z) colnames(M$V)[which.max(z)])
  do.call(rbind, lapply(setdiff(M$models, "Lasso"), function(m) {
    a <- dv[,m]; b <- dt[,m]; q <- a >= quantile(a, .75, na.rm=TRUE)
    data.frame(target=m, entry_mean=ag$mean[match(m, ag$target)], optimism=ag$optimism[match(m, ag$target)],
               mean_dval=mean(a,na.rm=TRUE), sd_dval=sd(a,na.rm=TRUE),
               winrate=mean(win==m,na.rm=TRUE), tail_opt=mean((a-b)[q],na.rm=TRUE))
  }))
}

cohs <- unique(sub("^09b?_(.+)_[0-9]+\\.tsv$", "\\1",
                   list.files(TAB, pattern="^09b?_")))
RES <- list()
for (cn in cohs) {
  f1 <- file.path(TAB, sprintf("09_%s_100.tsv", cn)); f2 <- file.path(TAB, sprintf("09b_%s_100.tsv", cn))
  M1 <- if (file.exists(f1)) load_mat(f1) else NULL
  M2 <- if (file.exists(f2)) load_mat(f2) else NULL
  if (is.null(M1) || is.null(M2)) { say("[skip] ", cn, " 缺一套划分或数据不足"); next }
  A1 <- stat_rules(M1); A2 <- stat_rules(M2)
  B1 <- stat_family(M1); B2 <- stat_family(M2)
  F1 <- stat_feats(M1);   F2 <- stat_feats(M2)
  mg <- merge(F1, F2, by="target", suffixes=c("_s1","_s2"))
  sp <- suppressWarnings(cor(mg$entry_mean_s1, mg$entry_mean_s2, method="spearman"))
  cr <- function(v) { a <- cor(mg[[paste0(v,"_s1")]], mg$entry_mean_s1)
                      b <- cor(mg[[paste0(v,"_s2")]], mg$entry_mean_s2); c(a,b) }
  cd <- sapply(c("mean_dval","sd_dval","winrate","tail_opt"), cr)

  RES[[cn]] <- data.frame(cohort=cn, n1=nrow(M1$V), n2=nrow(M2$V),
    infl1=A1$infl, infl2=A2$infl,
    truegain1=A1$true_gain, truegain2=A2$true_gain,
    poolbias1=A1$pool_bias, poolbias2=A2$pool_bias,
    fam1=B1, fam2=B2, fam_same_sign=(sign(B1)==sign(B2)),
    spearman=sp,
    r_mean_dval_1=cd[1,1], r_mean_dval_2=cd[2,1],
    r_sd_dval_1=cd[1,2],   r_sd_dval_2=cd[2,2],
    r_winrate_1=cd[1,3],   r_winrate_2=cd[2,3],
    r_tail_opt_1=cd[1,4],  r_tail_opt_2=cd[2,4])
  say(""); say("===== ", cn, " (主 ", nrow(M1$V), " / 第二 ", nrow(M2$V), " 个基因集) =====")
  say(sprintf("  A. 虚高占报告增益:  主 %.0f%%   第二 %.0f%%     (真实增益 %.4f / %.4f)",
              A1$infl, A2$infl, A1$true_gain, A2$true_gain))
  say(sprintf("  A'. 全池总偏倚:      主 %+.4f   第二 %+.4f   （估计量本身的量，非相对排名）",
              A1$pool_bias, A2$pool_bias))
  say(sprintf("  B. 家族对比(树−线性): 主 %+.4f   第二 %+.4f   -> %s",
              B1, B2, if (sign(B1)==sign(B2)) "同号" else "★翻号★"))
  say(sprintf("  D. 进入效应排名一致性 Spearman: %+.3f", sp))
  say(sprintf("  C. 预测因子相关（主 / 第二）: mean_dval %+.2f/%+.2f  sd_dval %+.2f/%+.2f  winrate %+.2f/%+.2f  tail_opt %+.2f/%+.2f",
              cd[1,1], cd[2,1], cd[1,2], cd[2,2], cd[1,3], cd[2,3], cd[1,4], cd[2,4]))
}

if (length(RES)) {
  R <- do.call(rbind, RES)
  write.table(R, file.path(TAB, "18_robustness_full.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
  say(""); say("========== 跨队列汇总 ==========")
  say(sprintf("  队列数 %d", nrow(R)))
  say(sprintf("  A. 虚高占报告增益:  主划分均值 %.0f%% (范围 %.0f~%.0f)   第二划分均值 %.0f%% (范围 %.0f~%.0f)",
              mean(R$infl1), min(R$infl1), max(R$infl1), mean(R$infl2), min(R$infl2), max(R$infl2)))
  say(sprintf("     两套划分的逐队列相关: r = %+.3f", cor(R$infl1, R$infl2)))
  say(sprintf("     真实增益: 主 %.4f / 第二 %.4f", mean(R$truegain1), mean(R$truegain2)))
  say(sprintf("  A'. 全池总偏倚: 主 %.4f (范围 %+.4f~%+.4f)   第二 %.4f (范围 %+.4f~%+.4f)   (n=%d 队列)",
              mean(R$poolbias1), min(R$poolbias1), max(R$poolbias1),
              mean(R$poolbias2), min(R$poolbias2), max(R$poolbias2), nrow(R)))
  say(sprintf("  B. 家族对比同号的队列: %d / %d", sum(R$fam_same_sign), nrow(R)))
  say(sprintf("  D. 进入效应排名一致性: 中位 %+.3f (范围 %+.3f~%+.3f)；正相关 %d/%d",
              median(R$spearman), min(R$spearman), max(R$spearman),
              sum(R$spearman > 0), nrow(R)))
  say("  C. 预测因子相关（跨队列均值）:")
  for (v in c("mean_dval","sd_dval","winrate","tail_opt"))
    say(sprintf("     %-10s 主 %+.3f   第二 %+.3f   （差 %+.3f）", v,
                mean(R[[paste0("r_",v,"_1")]]), mean(R[[paste0("r_",v,"_2")]]),
                mean(R[[paste0("r_",v,"_2")]]) - mean(R[[paste0("r_",v,"_1")]])))
  say("")
  say("  注：A 与 C 是论文核心主张，必须跨划分稳定；B 与 D 已知不稳，如实报告。")
}
close(con)

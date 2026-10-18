# 核对 04_终版结果补充_40组完成.md 里引用的每一个数字
ok <- 0; bad <- 0
ck <- function(cond, lab) { v <- isTRUE(cond)
  cat(sprintf("  [%s] %s\n", if (v) "OK" else "**错**", lab)); if (v) ok <<- ok+1 else bad <<- bad+1 }
rd <- function(f) read.table(f, sep="\t", header=TRUE, quote="", comment.char="")
P <- "D:/bioinfo/sigAudit/results/tables/"

## 1 完成情况
cat("=== 1 完成情况 ===\n")
M <- rd("D:/bioinfo/sigAudit/results/manifest/run_manifest.tsv")
ck(nrow(M) == 40, "manifest 记录 40 个组合")
ck(all(M$n_planned == M$n_processed), "所有组合 n_planned == n_processed")
ck(length(unique(paste(M$cohort, M$seed))) == 40, "40 个组合唯一无重复")
ck(all(M$n_genesets_usable >= 94 & M$n_genesets_usable <= 99),
   sprintf("manifest 里拟合成功签名数范围 %d-%d", min(M$n_genesets_usable), max(M$n_genesets_usable)))
# "可评估任务"与"可用候选"以 20_coverage.tsv 为准（manifest 记录的是拟合成功数）
CV <- rd(paste0(P, "20_coverage.tsv"))
ck(min(CV$n_tasks) == 91 && max(CV$n_tasks) == 99,
   sprintf("可评估任务数范围 %d-%d（文中 91-99）", min(CV$n_tasks), max(CV$n_tasks)))
ck(min(CV$n_candidates) == 15 && max(CV$n_candidates) == 18,
   sprintf("可用候选范围 %d-%d（文中 15-18）", min(CV$n_candidates), max(CV$n_candidates)))

## 2 P1
cat("\n=== 2 P1 ===\n")
S <- rd(paste0(P, "20_summary.tsv"))
g <- function(r, col) S[[col]][S$rule == r]
ck(abs(g("argmax_val","Cv") - 0.6554) < 5e-5, sprintf("argmax_val Cv = %.4f（文中 0.6554）", g("argmax_val","Cv")))
ck(abs(g("argmax_val","Ct") - 0.6079) < 5e-5, sprintf("argmax_val Ct = %.4f（文中 0.6079）", g("argmax_val","Ct")))
ck(abs(g("argmax_val","Gv") - 0.0610) < 5e-5, sprintf("Gv = %+.4f（文中 +0.0610）", g("argmax_val","Gv")))
ck(abs(g("argmax_val","Gt") - 0.0237) < 5e-5, sprintf("Gt = %+.4f（文中 +0.0237）", g("argmax_val","Gt")))
ck(abs(g("argmax_val","Orel") - 0.0374) < 5e-5, sprintf("Orel = %+.4f（文中 +0.0374）", g("argmax_val","Orel")))
pct <- 100*g("argmax_val","Orel")/g("argmax_val","Gv")
ck(abs(pct - 61.2) < 0.2, sprintf("比例 = %.1f%%（文中 61.2%%）", pct))
ck(abs(g("none","Cv") - 0.5944) < 5e-5 && abs(g("none","Ct") - 0.5842) < 5e-5,
   sprintf("none Cv/Ct = %.4f/%.4f（文中 0.5944/0.5842）", g("none","Cv"), g("none","Ct")))

## 3 P2
cat("\n=== 3 P2 ===\n")
Q <- rd(paste0(P, "20_paired_vs_argmax.tsv"))
cl <- aggregate(dCt ~ cohort + rule, data=Q, FUN=function(x) mean(x, na.rm=TRUE))
ag <- do.call(rbind, lapply(split(cl, cl$rule), function(d){ x <- d$dCt[is.finite(d$dCt)]
  tt <- t.test(x); data.frame(rule=d$rule[1], n=length(x), mean=mean(x), lo=tt$conf.int[1], hi=tt$conf.int[2], p=tt$p.value) }))
chk <- function(r, m, lo, hi, pv, tol=5e-5) {
  z <- ag[ag$rule == r, ]
  if (!nrow(z)) { ck(FALSE, paste("缺规则", r)); return(invisible()) }
  ck(abs(z$mean-m)<tol && abs(z$lo-lo)<tol && abs(z$hi-hi)<tol && abs(z$p-pv)/pv < 0.05,
     sprintf("%-20s ΔCt %+.4f [%+.4f,%+.4f] p=%.2e", r, z$mean, z$lo, z$hi, z$p)) }
chk("top5_ensemble",   0.0048, 0.0012, 0.0085, 0.0170)
chk("top3_ensemble",   0.0037, 0.0014, 0.0059, 0.0067)
chk("top2_ensemble",   0.0023, 0.0012, 0.0034, 0.0017)
chk("top10_ensemble",  0.0046, -0.0007, 0.0099, 0.0805)
chk("all_ensemble",    0.0045, -0.0022, 0.0112, 0.155)
chk("linear_only",     0.0014, -0.0004, 0.0033, 0.107)
chk("runner_up",      -0.0010, -0.0020, 0.0001, 0.0610)
chk("majority_vote",  -0.0012, -0.0021, -0.0002, 0.0229)
chk("worst_half",     -0.0020, -0.0037, -0.0004, 0.0218)
chk("penalized",      -0.0025, -0.0042, -0.0008, 0.0108)
chk("tree_only",      -0.0043, -0.0102, 0.0015, 0.120)
chk("none",           -0.0237, -0.0352, -0.0121, 0.00185)
ck(all(ag$n == 8), "所有规则 n = 8（队列数）")

## 4 跨划分稳定性
cat("\n=== 4 跨划分稳定性 ===\n")
R <- rd(paste0(P, "20_rules.tsv")); R <- R[R$rule == "argmax_val", ]
ST <- do.call(rbind, lapply(split(R, R$cohort), function(d)
  data.frame(cohort=d$cohort[1], mean=mean(d$Orel), rng=diff(range(d$Orel)))))
exp_mean <- c(BRCA=0.0550, COAD=0.0396, KIRC=0.0170, LGG=0.0219, LUAD=0.0370, MESO=0.0484, PAAD=0.0430, STAD=0.0370)
exp_rng  <- c(BRCA=0.0546, COAD=0.0851, KIRC=0.0169, LGG=0.0187, LUAD=0.0311, MESO=0.0700, PAAD=0.0094, STAD=0.0456)
for (cn in names(exp_mean)) {
  z <- ST[ST$cohort == cn, ]
  ck(abs(z$mean-exp_mean[cn])<5e-5 && abs(z$rng-exp_rng[cn])<5e-5,
     sprintf("%-5s 均值 %+.4f 极差 %.4f（文中 %+.4f / %.4f）", cn, z$mean, z$rng, exp_mean[cn], exp_rng[cn]))
}
ck(all(ST$mean > 0), "8 个队列的 Orel 均值全部为正")

## 5 辅助分析
cat("\n=== 5 辅助 21 ===\n")
E <- rd(paste0(P, "21_entry_effect.tsv"))
a2 <- aggregate(entry ~ target + k, data=E, FUN=function(x) mean(x, na.rm=TRUE))
a3 <- aggregate(entry ~ k, data=a2, FUN=function(x) x[which.max(x)])
topk <- sapply(sort(unique(a2$k)), function(k){ s <- a2[a2$k==k,]; s$target[which.max(s$entry)] })
ck(all(topk == "RSF100"), sprintf("各规模下进入效应最高的都是 RSF100: %s", paste(unique(topk), collapse=",")))
mc <- E$mc_se[is.finite(E$mc_se)]
med_k <- aggregate(mc_se ~ k, data=E[is.finite(E$mc_se),], FUN=function(x) median(x))
ck(abs(med_k$mc_se[med_k$k==2] - 0.000811) < 1e-6 && abs(med_k$mc_se[med_k$k==12] - 0.000054) < 1e-6,
   sprintf("各规模中位 MCSE: k=2 %.6f ... k=12 %.6f（文中 0.000811/0.000054）",
           med_k$mc_se[med_k$k==2], med_k$mc_se[med_k$k==12]))
ck(abs(median(mc) - 0.000262) < 1e-6 && abs(max(mc) - 0.005045) < 1e-6,
   sprintf("MCSE 全体: 中位 %.6f 最大 %.6f（文中 0.000262/0.005045）", median(mc), max(mc)))
F <- rd(paste0(P, "21_predictors.tsv"))
fa <- aggregate(cbind(entry_full, mean_dval, sd_dval, winrate, optimism, tail_opt) ~ cohort + target,
                data=F, FUN=function(x) mean(x, na.rm=TRUE))
ck(nrow(fa) == 128, sprintf("候选层面 n = %d（文中 128）", nrow(fa)))
exp_r <- c(tail_opt=0.635, winrate=0.623, mean_dval=0.565, sd_dval=0.547, optimism=0.405)
for (v in names(exp_r)) { rr <- cor(fa[[v]], fa$entry_full)
  ck(abs(rr - exp_r[v]) < 0.002, sprintf("%-10s r = %+.3f（文中 %+.3f）", v, rr, exp_r[v])) }
LO <- rd(paste0(P, "21_loo_sensitivity.tsv"))
vals <- unlist(LO[, c("r_mean_dval","r_sd_dval","r_winrate","r_tail_opt")])
ck(min(vals) >= 0.44 && max(vals) <= 0.69, sprintf("留一后相关范围 %.2f–%.2f（文中 0.45–0.68）", min(vals), max(vals)))

## 6 图
cat("\n=== 6 图 ===\n")
for (f in c("fig1_design","fig2_gain_decomposition","fig3_cross_split","fig4_rule_pairing"))
  ck(file.exists(sprintf("D:/bioinfo/sigAudit/results/figures/%s.pdf", f)) &&
     file.exists(sprintf("D:/bioinfo/sigAudit/results/figures/%s.png", f)), paste(f, "PDF+PNG"))

cat(sprintf("\n核对结果: %d 项一致, %d 项不一致\n", ok, bad))
if (bad > 0) quit(status = 1)

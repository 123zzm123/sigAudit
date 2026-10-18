# 在【同一个 19 产物】上，比较"argmax 允许 vs 不允许 Lasso 参与"的差异。
# 这正是评审 B3 指出的口径问题，必须用同一数据严格对照。
suppressMessages(library(survival))
R <- readRDS("D:/bioinfo/sigAudit/results/main/19_MESO_20261014_100.rds")
G <- R$genesets
ci <- function(t,e,s){ s<-as.numeric(s); ok<-is.finite(s)&!is.na(t)&!is.na(e)
  if(sum(ok)<10) return(NA_real_)
  as.numeric(concordance(Surv(t[ok],e[ok])~s[ok], reverse=TRUE, timewt="n")$concordance) }
mods_all <- names(sort(table(unlist(lapply(G,function(z) colnames(z$val_pred)))), decreasing=TRUE))
cnt <- table(unlist(lapply(G,function(z) colnames(z$val_pred))))
keep <- names(cnt)[cnt/length(G) >= 0.90]
use <- Filter(function(z) all(keep %in% colnames(z$val_pred)), G)
cat("可评估任务:", length(use), " 保留候选:", length(keep), "\n\n")
res <- list()
for (z in use) {
  P <- z$val_pred[, keep, drop=FALSE]; Q <- z$test_pred[, keep, drop=FALSE]
  V <- sapply(keep, function(m) ci(z$tva, z$eva, P[,m])); Tt <- sapply(keep, function(m) ci(z$tte, z$ete, Q[,m]))
  cv <- V[["Lasso"]]; ct <- Tt[["Lasso"]]
  jA <- keep[which.max(V)];                     cand_all <- keep
  cand_no <- setdiff(keep, "Lasso")
  jB <- cand_no[which.max(V[cand_no])]
  res[[z$name]] <- c(
    GvA = V[jA]-cv, GtA = Tt[jA]-ct,   # 全池（含 Lasso）
    GvB = V[jB]-cv, GtB = Tt[jB]-ct,   # 排除 Lasso（原版口径）
    lasso_wins = as.integer(jA == "Lasso"))
}
D <- as.data.frame(do.call(rbind, res))
cat(sprintf("【全池（含 Lasso）】   Gv %+.4f  Gt %+.4f  Orel %+.4f\n",
            mean(D$GvA), mean(D$GtA), mean(D$GvA-D$GtA)))
cat(sprintf("【排除 Lasso（原版）】 Gv %+.4f  Gt %+.4f  Orel %+.4f\n",
            mean(D$GvB), mean(D$GtB), mean(D$GvB-D$GtB)))
cat(sprintf("\nargmax 选中 Lasso 的任务数: %d / %d (%.1f%%)\n",
            sum(D$lasso_wins), nrow(D), 100*mean(D$lasso_wins)))
cat("→ 若 Lasso 从不被选中，两种口径应当完全相同；实际差异即来自该点。\n")

cat(sprintf("\n精确差值: ΔGv %+.6f  ΔGt %+.6f  ΔOrel %+.6f\n",
            mean(D$GvA-D$GvB), mean(D$GtA-D$GtB), mean(D$Orel_A <- (D$GvA-D$GtA), na.rm=TRUE)*0 +
            mean((D$GvA-D$GtA)-(D$GvB-D$GtB))))
cat(sprintf("最大单任务差: ΔGt %+.4f\n", max(abs(D$GtA-D$GtB))))

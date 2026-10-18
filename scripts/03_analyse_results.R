#!/usr/bin/env Rscript
# ============================================================
# 03 — 汇总分析实验结果的跨队列脚本
#
# 输入：results/tables/01_<COHORT>_<N>.tsv （每个队列一份）
# 产出：
#   results/tables/03_summary.tsv     各队列的关键指标
#   results/figures/03_*.pdf          核心图
#   results/tables/03_stats.tsv       配对检验
# ============================================================
suppressMessages({library(ggplot2)})
options(stringsAsFactors=FALSE)
PROJ <- "D:/bioinfo/sigAudit"
TAB <- file.path(PROJ,"results/tables"); FIG <- file.path(PROJ,"results/figures")
dir.create(FIG, showWarnings=FALSE, recursive=TRUE)

fs <- list.files(TAB, pattern="^01_.*\\.tsv$", full.names=TRUE)
if (!length(fs)) stop("没有找到实验结果文件")
cat("发现结果文件:", length(fs), "\n")

all <- do.call(rbind, lapply(fs, function(f) {
  d <- read.table(f, header=TRUE, sep="\t", check.names=FALSE)
  if (!nrow(d)) return(NULL)
  # 从文件名解析队列与 N：01_<COHORT>_<N>.tsv
  b <- sub("^01_","",sub("\\.tsv$","",basename(f)))
  parts <- strsplit(b,"_")[[1]]
  d$cohort <- parts[1]; d$N_planned <- as.integer(tail(parts,1))
  d
}))
cat("总记录:", nrow(all), " 队列:", paste(unique(all$cohort), collapse=", "), "\n\n")

# ---- 各队列汇总 ----
summ <- do.call(rbind, lapply(split(all, all$cohort), function(d) {
  data.frame(cohort = d$cohort[1], n_geneset = nrow(d),
             mean_train = mean(d$max_train, na.rm=TRUE),
             mean_test  = mean(d$max_test,  na.rm=TRUE),
             mean_selA  = mean(d$c_mean,    na.rm=TRUE),
             mean_selB  = mean(d$c_innercv, na.rm=TRUE),
             mean_oracle= mean(d$c_oracle,  na.rm=TRUE),
             inflation_A_vs_B = mean(d$c_mean - d$c_innercv, na.rm=TRUE),
             inflation_best = mean(d$max_train - d$max_test, na.rm=TRUE))
}))
summ <- summ[order(-summ$inflation_A_vs_B), ]
write.table(summ, file.path(TAB,"03_summary.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
cat("=== 各队列汇总 ===\n"); print(summ, row.names=FALSE, digits=3)

# ---- 配对检验：选法A vs 诚实B ----
cat("\n=== 配对检验 (选法A 与 诚实B 的差) ===\n")
st <- do.call(rbind, lapply(split(all, all$cohort), function(d) {
  ok <- !is.na(d$c_mean) & !is.na(d$c_innercv)
  if (sum(ok) < 10) return(NULL)
  tt <- t.test(d$c_mean[ok], d$c_innercv[ok], paired=TRUE)
  data.frame(cohort=d$cohort[1], n=sum(ok),
             mean_diff=mean(d$c_mean[ok]-d$c_innercv[ok]),
             CI_low=tt$conf.int[1], CI_high=tt$conf.int[2], p=tt$p.value)
}))
write.table(st, file.path(TAB,"03_stats.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
print(st, row.names=FALSE, digits=3)

# ---- 图 1：三种选法的测试集表现对比 ----
lf <- reshape(all[, c("cohort","c_mean","c_innercv","c_oracle")],
              varying=c("c_mean","c_innercv","c_oracle"),
              v.names="cindex", times=c("A 训练+测试均值","B 内部CV(诚实)","C oracle"),
              timevar="rule", direction="long")
p1 <- ggplot(lf, aes(x=cohort, y=cindex, fill=rule)) +
  geom_boxplot(outlier.size=0.3, linewidth=0.3) + theme_classic() +
  labs(x=NULL, y="C-index on held-out test set", fill="模型选择规则",
       title="Model-selection rule inflates reported performance") +
  theme(legend.position="top")
ggsave(file.path(FIG,"03_rules_by_cohort.pdf"), p1, width=8, height=4.6)

# ---- 图 2：训练 vs 测试的过拟合散点 ----
p2 <- ggplot(all, aes(x=max_train, y=max_test)) +
  geom_point(alpha=0.35, size=0.8) + geom_abline(slope=1, linetype=2, colour="red") +
  facet_wrap(~cohort, scales="free") + theme_classic() +
  labs(x="Best C-index on training set", y="Best C-index on test set",
       title="Train-test gap by cohort (red line = no optimism)")
ggsave(file.path(FIG,"03_train_test_gap.pdf"), p2, width=9, height=4)

cat("\n图已输出到 results/figures/\n")

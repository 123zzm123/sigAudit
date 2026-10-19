#!/usr/bin/env Rscript
# ============================================================
# 从本地冻结预测导出一个【真实数据、可直接调用安装后包】的示例
#
# 动机（第 8 包评审）：仓库里的随机分数示例只能证明"接口能调用"，
# 不能展示软件对真实生物医学预测的用途。这里导出一个真实的、派生自
# METABRIC 冻结预测的示例（只有风险分数与生存结局，不含表达矩阵），
# 供安装后的 sigAudit 直接审计。
#
# 选择规则（固定、与效果无关）：
#   signature = 六签名中【字典序最小】者；seed = 五个划分种子中【最小】者。
#   -> EndoPredict_11 / 20261014。不按 C-index 好坏挑选。
#
# 去标识化：验证/测试患者 ID 重写为 V0001.. 与 T0001..；不含表达值、
#   不含原始样本条码。原数据为公开去标识数据，本示例只发布模型输出。
#
# 用法（在项目根目录，且已有 results/metabric/31_metabric_predictions.rds）：
#   Rscript inst/examples/export_audit_example.R
# 输出：inst/extdata/audit_example_val.tsv
#       inst/extdata/audit_example_test.tsv
#       inst/extdata/audit_example_survival.tsv
# ============================================================
options(stringsAsFactors = FALSE)
PROJ <- Sys.getenv("SIGAUDIT_ROOT", "D:/bioinfo/sigAudit")
PRF  <- file.path(PROJ, "results/metabric/31_metabric_predictions.rds")
OUT  <- file.path(PROJ, "inst/extdata")
if (!file.exists(PRF)) stop("缺少冻结预测：", PRF, "（本脚本只能导出本地已有的派生预测）")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
P <- readRDS(PRF)

sig  <- sort(unique(vapply(P, function(z) z$signature, character(1))))[1]
seeds <- sort(unique(vapply(P, function(z) as.character(z$seed), character(1))))
seed <- seeds[1]
z <- P[[which(vapply(P, function(z2) identical(z2$signature, sig) &&
                             identical(as.character(z2$seed), seed), logical(1)))[1]]]
stopifnot(!is.null(z), identical(z$signature, sig))

Pv <- as.matrix(z$val_pred); Pt <- as.matrix(z$test_pred)
nv <- nrow(Pv); nt <- nrow(Pt)
vid <- sprintf("V%04d", seq_len(nv)); tid <- sprintf("T%04d", seq_len(nt))

val <- data.frame(id = vid, Pv, check.names = FALSE, stringsAsFactors = FALSE)
tst <- data.frame(id = tid, Pt, check.names = FALSE, stringsAsFactors = FALSE)
sur <- rbind(
  data.frame(id = vid, set = "validation", time = z$t_val, event = z$e_val,
             stringsAsFactors = FALSE),
  data.frame(id = tid, set = "test", time = z$t_test, event = z$e_test,
             stringsAsFactors = FALSE))

wr <- function(x, f) write.table(x, file.path(OUT, f), sep = "\t", quote = FALSE,
                                 row.names = FALSE, na = "NA")
wr(val, "audit_example_val.tsv")
wr(tst, "audit_example_test.tsv")
wr(sur, "audit_example_survival.tsv")

cat("示例已导出：", OUT, "\n")
cat(sprintf("  signature=%s  seed=%s  candidates=%d  val=%d  test=%d  events=%d/%d\n",
            sig, seed, ncol(Pv), nv, nt, sum(z$e_val), sum(z$e_test)))
cat("  （ID 已重写为 V####/T####；仅含风险分数与生存结局）\n")

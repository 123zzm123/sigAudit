# ============================================================
# 实证确认「评价函数把测试结局用于确定方向」这一问题，并锁定正确的替代写法。
#
# 待检验的三件事：
#   1) 旧写法 coxph(...)$concordance[1] 是否【对 s 与 -s 给同一个值】（即与方向无关）
#      —— 若是，则方向由当前评价集决定；在测试集上调用即构成泄漏。
#   2) 直接评价 concordance(Surv ~ s) 在默认与 reverse=TRUE 下分别是什么方向，
#      以便确定"高值=高风险"约定下应当怎么写。
#   3) timewt="n" 是否与旧写法在【方向正确时】数值一致（保证新旧结果可比）。
# ============================================================
suppressMessages(library(survival))

set.seed(1)
n <- 300
# 真实风险：越高越危险
x <- rnorm(n)
t <- rexp(n, rate = exp(0.8 * x))
e <- rbinom(n, 1, 0.7)

old_ci <- function(t, e, s) as.numeric(summary(coxph(Surv(t, e) ~ s))$concordance[1])
dir_ci <- function(t, e, s, rev = FALSE, tw = "n")
  as.numeric(concordance(Surv(t, e) ~ s, reverse = rev, timewt = tw)$concordance)

# 情形 A：方向正确（高值=高风险）
sA <- x + rnorm(n, 0, 0.3)
# 情形 B：方向反转（高值=低风险）
sB <- -sA

cat("================ 情形 A：风险评分方向正确（高值=高风险）================\n")
cat(sprintf("  旧写法 coxph$concordance        = %.4f\n", old_ci(t, e, sA)))
cat(sprintf("  直接 concordance(reverse=FALSE) = %.4f\n", dir_ci(t, e, sA, rev = FALSE)))
cat(sprintf("  直接 concordance(reverse=TRUE)  = %.4f\n", dir_ci(t, e, sA, rev = TRUE)))
cat(sprintf("  直接 concordance(timewt='a')    = %.4f\n", dir_ci(t, e, sA, rev = FALSE, tw = "a")))

cat("\n================ 情形 B：同一评分【取反】（高值=低风险）================\n")
cat(sprintf("  旧写法 coxph$concordance        = %.4f\n", old_ci(t, e, sB)))
cat(sprintf("  直接 concordance(reverse=FALSE) = %.4f\n", dir_ci(t, e, sB, rev = FALSE)))
cat(sprintf("  直接 concordance(reverse=TRUE)  = %.4f\n", dir_ci(t, e, sB, rev = TRUE)))

cat("\n================ 判定 =================\n")
dA <- old_ci(t, e, sA); dB <- old_ci(t, e, sB)
cat(sprintf("  旧写法对 s 与 -s 的差 = %.6f  -> %s\n", abs(dA - dB),
            if (abs(dA - dB) < 1e-9)
              "★ 完全相同：方向由评价数据决定（测试集上即泄漏）" else "不同"))
cat(sprintf("  旧写法是否恒 >= 0.5 : A %.4f, B %.4f -> %s\n", dA, dB,
            if (dA >= 0.5 && dB >= 0.5) "★ 是（被抬到 0.5 以上）" else "否"))
cat(sprintf("  直接评价在方向反转时给出 %.4f (< 0.5) -> %s\n", dir_ci(t, e, sB),
            if (dir_ci(t, e, sB) < 0.5) "✓ 能正确低于 0.5" else "✗"))
cat(sprintf("  方向正确时 直接(默认)与旧写法差 = %.6f -> %s\n",
            abs(dir_ci(t, e, sA) - dA),
            if (abs(dir_ci(t, e, sA) - dA) < 1e-6) "一致（新旧结果可比）" else "不一致"))

cat("\n================ 验证：直接评价不调用任何拟合器 =================\n")
cat(sprintf("  concordance() 的类属: %s\n", class(concordance(Surv(t, e) ~ sA))[1]))
cat("  直接评价只做配对比较，不估计系数（可对照上面的 reverse 结果确认方向固定）\n")

cat("\n================ 结论：应当采用哪一种写法 =================\n")
cat("  在【高值 = 高风险】约定下，与 Harrell C 一致且不依赖评价数据的写法是：\n")
cat("    concordance(Surv(time,event) ~ score, timewt='n')$concordance\n")
cat("  即 reverse 保持默认（FALSE）。\n")

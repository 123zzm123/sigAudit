suppressMessages(library(survival))
set.seed(1); n <- 300
x <- rnorm(n); t <- rexp(n, rate = exp(0.8 * x)); e <- rbinom(n, 1, 0.7)
s <- x + rnorm(n, 0, 0.3)          # 约定：高值 = 高风险

old_ci <- function(t, e, s) as.numeric(summary(coxph(Surv(t, e) ~ s))$concordance[1])
new_ci <- function(t, e, s) as.numeric(concordance(Surv(t, e) ~ s, reverse = TRUE, timewt = "n")$concordance)

cat(sprintf("旧写法  coxph$concordance(s)     = %.6f\n", old_ci(t, e, s)))
cat(sprintf("旧写法  coxph$concordance(-s)    = %.6f\n", old_ci(t, e, -s)))
cat(sprintf("  -> 两者之差 = %.2e   %s\n", abs(old_ci(t, e, s) - old_ci(t, e, -s)),
    if (abs(old_ci(t, e, s) - old_ci(t, e, -s)) < 1e-12) "★ 完全相同：方向由评价数据决定" else "不同"))
cat("\n")
cat(sprintf("新写法  concordance(reverse=TRUE)(s)  = %.6f\n", new_ci(t, e, s)))
cat(sprintf("新写法  concordance(reverse=TRUE)(-s) = %.6f\n", new_ci(t, e, -s)))
cat(sprintf("  -> 与旧写法(方向正确时)之差 = %.2e   %s\n", abs(new_ci(t, e, s) - old_ci(t, e, s)),
    if (abs(new_ci(t, e, s) - old_ci(t, e, s)) < 1e-9) "✓ 方向正确时数值一致（新旧可比）" else "✗"))
cat(sprintf("  -> 取反后 = %.6f，是否等于 1-原值 %.6f : %s\n",
    new_ci(t, e, -s), 1 - new_ci(t, e, s),
    if (abs(new_ci(t, e, -s) - (1 - new_ci(t, e, s))) < 1e-9) "✓ 能把方向反过来的评分正确打到 0.5 以下" else "✗"))
cat("\n")
cat("=== 再看一个【与风险无关的随机评分】（C 真值应为 0.5 附近）===\n")
r <- rnorm(n)
cat(sprintf("  旧写法 = %.4f   新写法 = %.4f\n", old_ci(t, e, r), new_ci(t, e, r)))
set.seed(7); r2 <- rnorm(n)
cat(sprintf("  换一个随机评分: 旧写法 = %.4f   新写法 = %.4f\n", old_ci(t, e, r2), new_ci(t, e, r2)))
cat("  -> 旧写法应稳定地 >= 0.5（被抬到 0.5 以上）；新写法可自由低于 0.5\n")

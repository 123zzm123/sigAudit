## 内部工具函数的断言式测试。
## 重点覆盖第三轮审查确认的两个错误与一个约定：
##   ① `.n_pairs()` 必须只计 concordant + discordant + tied.x
##      （旧实现把 tied.y / tied.xy 也加进分母，把 MESO 测试集报成 224 而非 223）
##   ② `.c_index()` 方向固定，不做 Cox 重拟合（否则反向评分会同样拿到高分）
##   ③ `.fixed_prio()` 只依赖候选名字集合与 seed，与输入顺序无关

test_that(".c_index 方向固定：高分为高风险时 C > 0.5，反向评分恰好互补", {
  set.seed(6)
  t <- rexp(80, 1); e <- rbinom(80, 1, 0.7)
  s <- -t + stats::rnorm(80, 0, .2)          # 生存时间越短 -> 分数越高 -> 高风险
  c1 <- .c_index(t, e, s); c0 <- .c_index(t, e, -s)
  expect_gt(c1, 0.5)
  expect_lt(c0, 0.5)
  # 方向固定意味着不会出现"反向评分也拿到高分"（旧 Cox 重拟合写法的泄漏表现）
  expect_equal(c1 + c0, 1, tolerance = 1e-8)
})

test_that(".n_pairs 只计 concordant + discordant + tied.x，不含 tied.y / tied.xy", {
  set.seed(4)
  n <- 90
  # 刻意造并列：三例同时点、两例同时点、以及大量重复时间
  time  <- c(5, 5, 5, 8, 8, 10, sort(sample(1:20, n - 6, replace = TRUE)))
  event <- c(1, 1, 0, 1, 0, 1, stats::rbinom(n - 6, 1, 0.6))
  cd  <- survival::concordance(survival::Surv(time, event) ~ seq_along(time), timewt = "n")
  cts <- cd$count
  # 与 survival 自己的 count 字段严格一致
  expect_equal(.n_pairs(time, event),
               unname(cts["concordant"] + cts["discordant"] + cts["tied.x"]))
  # 旧口径（把所有 tied.* 相加）在这里必然更大
  expect_lt(.n_pairs(time, event), unname(sum(cts)))
  # 存在并列时间时，两者必须真的不同（否则这个测试没有检出能力）
  expect_gt(unname(cts["tied.y"]), 0)
})

test_that(".fixed_prio 只依赖候选名字集合与 seed，与输入顺序无关", {
  a <- .fixed_prio(c("B", "A", "C"), 3)
  b <- .fixed_prio(c("C", "B", "A"), 3)
  expect_identical(a, b)
  expect_identical(names(a), c("A", "B", "C"))
  expect_false(isTRUE(all.equal(.fixed_prio(c("A","B","C"), 3),
                                .fixed_prio(c("A","B","C"), 4))))
})

test_that(".validate_ci 把非有限值记为 NA 并计数，越界值报错", {
  m <- matrix(c(0.6, Inf, NaN, 0.7), 2, 2)
  r <- .validate_ci(m)
  expect_equal(r$n_nonfinite, 2)
  expect_true(all(is.na(r$m[!is.finite(m)])))
  expect_error(.validate_ci(matrix(c(0.6, 1.2), 1, 2)), "outside")
})

test_that(".zscore_genes 的标准化参数只在训练侧估计", {
  set.seed(9)
  X <- matrix(stats::rnorm(6 * 40), 6, 40, dimnames = list(paste0("G", 1:6), paste0("S", 1:40)))
  X[1, 21:40] <- X[1, 21:40] * 100          # 测试侧整体平移放大
  tr <- 1:20
  Z <- .zscore_genes(X, genes = rownames(X), train = tr)
  # 训练侧标准化后均值≈0、sd≈1；测试侧被训练参数缩放，绝不可能恰好也≈0/1
  expect_lt(abs(mean(Z[1, tr])), 1e-8)
  expect_gt(abs(mean(Z[1, 21:40])), 1)
})

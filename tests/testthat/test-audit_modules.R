## 覆盖 audit_optimism / audit_transport / audit_concordance /
## grade_credibility / print 的出错路径。这些模块此前完全没有测试。

make_cohort <- function(n = 150, p = 60, seed = 1, shift = 0) {
  set.seed(seed)
  X <- matrix(rnorm(n * p), n, p,
              dimnames = list(paste0("S", seq_len(n)), paste0("G", seq_len(p))))
  lp <- X[, 1] * 1.2 + X[, 2] * -0.9 + X[, 3] * 0.7
  tt <- rexp(n, rate = exp(lp * 0.5 + shift))
  ev <- rbinom(n, 1, 0.7)
  tt[ev == 0] <- tt[ev == 0] * 1.6          # 造一些删失
  surv <- data.frame(time = tt, status = ev, row.names = rownames(X))
  list(expr = t(X), surv = surv, genes = paste0("G", 1:6))
}

test_that("audit_optimism 在合成数据上给出有限值", {
  d <- make_cohort(seed = 11)
  r <- audit_optimism(d$expr, d$surv, "time", "status", d$genes, B = 25)
  expect_s3_class(r, "sigAudit_module")
  expect_true(is.finite(r$apparent))
  expect_true(is.finite(r$corrected))
  expect_gt(r$apparent, 0.5)                 # 植入的 3 个基因确实有信号
  expect_equal(r$corrected, r$apparent - r$optimism)
  expect_output(print(r), "corrected")
})

test_that("audit_transport 训练/外部队列能跑通并报出落差", {
  a <- make_cohort(seed = 21); b <- make_cohort(seed = 22, shift = 0.3)
  r <- audit_transport(a$expr, a$surv, "time", "status",
                       b$expr, b$surv, "time", "status", a$genes)
  expect_s3_class(r, "sigAudit_module")
  expect_true(is.finite(r$apparent) && is.finite(r$external))
  expect_equal(r$gap, r$apparent - r$external)
  expect_identical(r$missing_in_external, character(0))   # 两队列基因相同
  expect_output(print(r), "external")
})

test_that("audit_transport 报出外部队列缺失的签名基因", {
  a <- make_cohort(seed = 31); b <- make_cohort(seed = 32)
  b$expr <- b$expr[rownames(b$expr) != "G2", , drop = FALSE]   # 故意删掉 G2
  r <- audit_transport(a$expr, a$surv, "time", "status",
                       b$expr, b$surv, "time", "status", a$genes)
  expect_true("G2" %in% r$missing_in_external)
  expect_output(print(r), "G2")
})

test_that("audit_concordance 能算出 mRNA-蛋白相关", {
  d <- make_cohort(seed = 41)
  set.seed(5)
  prot <- t(d$expr[1:6, ]) + matrix(rnorm(6 * 150, 0, 0.3), 150, 6)  # 蛋白≈mRNA+噪声
  dimnames(prot) <- list(colnames(d$expr), rownames(d$expr)[1:6])
  r <- audit_concordance(d$expr, prot, d$genes)
  expect_s3_class(r, "sigAudit_module")
  expect_true(!is.null(r$table))
  expect_true(all(c("gene", "rho", "p") %in% names(r$table)))
  expect_gt(mean(r$table$rho), 0.5)          # 人工造的高相关应当被检出
})

test_that("audit_concordance 在样本/基因不足时优雅返回而不是报错", {
  d <- make_cohort(seed = 51)
  prot <- t(d$expr[1:6, 1:10, drop = FALSE])
  dimnames(prot) <- list(colnames(d$expr)[1:10], rownames(d$expr)[1:6])
  r <- audit_concordance(d$expr, prot, d$genes, min_n = 30)
  expect_s3_class(r, "sigAudit_module")
  expect_null(r$table)
  expect_true(!is.null(r$message))
})

test_that("grade_credibility 遇到 NA 不崩（小队列的 C-index 就是 NA）", {
  rep <- structure(list(
    genes = "g", n_gene = 1, method = "lasso",
    optimism = structure(list(module = "optimism", optimism = NA_real_,
                              corrected = NA_real_, n = 10, events = 3,
                              n_gene = 3, B = 0, apparent = NA_real_),
                         class = "sigAudit_module"),
    transport = list(cohort1 = structure(list(module = "transportability",
                              gap = NA_real_, external = NA_real_,
                              missing_in_external = character(0),
                              n_gene = 3, apparent = NA_real_,
                              n_train = 10, n_external = 10, events_external = 3),
                         class = "sigAudit_module")),
    concordance = NULL), class = "sigAudit_report")
  g <- grade_credibility(rep)
  expect_s3_class(g, "sigAudit_summary")
  expect_true(is.list(g$estimates) && is.list(g$coverage))
  expect_false("grade" %in% names(g))     # 已移除 A–D 总评级
  expect_false("points" %in% names(g))
  expect_output(print(g), "no letter grade")
})

test_that("出错模块打印出错误原因，而不是一片空白", {
  m <- structure(list(module = "optimism", error = "boom"), class = "sigAudit_module")
  expect_output(print(m), "boom")
  r <- structure(list(genes = "g", n_gene = 1, method = "lasso",
                      optimism = m, transport = list(), concordance = NULL),
                 class = "sigAudit_report")
  expect_output(print(r), "boom")
})

test_that("grade_credibility 会把选择偏倚模块计入分级", {
  set.seed(7); n_gs <- 40
  nm <- c("REF", "A", "B", "C", "BAD")
  V  <- sapply(nm, function(z) switch(z,
    REF = rnorm(n_gs, .60, .03), A = rnorm(n_gs, .61, .03), B = rnorm(n_gs, .61, .03),
    C = rnorm(n_gs, .60, .03), BAD = rnorm(n_gs, .66, .07)))
  Tt <- sapply(nm, function(z) switch(z,
    REF = rnorm(n_gs, .58, .03), A = rnorm(n_gs, .59, .03), B = rnorm(n_gs, .59, .03),
    C = rnorm(n_gs, .58, .03), BAD = rnorm(n_gs, .55, .07)))
  dimnames(V) <- dimnames(Tt) <- list(paste0("gs", 1:n_gs), nm)
  sel <- audit_selection(V, Tt, reference = "REF")
  r0 <- structure(list(genes = "g", n_gene = 1, method = "lasso",
                       optimism = NULL, transport = list(), concordance = NULL),
                  class = "sigAudit_report")
  r1 <- r0; r1$selection <- sel
  g1 <- grade_credibility(r1)
  # 不再有分数/等级；改为检查池级偏倚估计与覆盖信息被正确带出
  expect_true("pool_bias" %in% names(g1$estimates))
  expect_true(is.finite(g1$estimates$pool_bias))
  expect_true("selection" %in% names(g1$coverage))
  expect_true(any(grepl("marginal", g1$notes)))          # 归因的限定说明必须存在
  expect_false("grade" %in% names(g1))
})

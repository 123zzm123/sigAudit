## 预测级审计入口 audit_predictions() 的断言式测试。
## 覆盖会改变结论/行为的点（对照大修后评审提出的真实风险）：
##   ① Orel = Gv − Gt 逐任务成立；none 规则恒为 0
##   ② top_k 真正生效（top_k=4 会生成并计算 top4_ensemble）
##   ③ 家族规则可请求且会返回；未知规则名报错，不静默丢弃
##   ④ 默认规则不含 oracle_test
##   ⑤ by_cohort/summary 用"先队列×划分、再队列"，不等于任务层平均
##   ⑥ 输入契约：事件长度、事件 0/1、候选名唯一、ID 唯一且两侧不重叠
##   ⑦ 任务级共同有限掩码：丢弃被记录，不按候选各自删缺失
##   ⑧ 验证/测试样本数可不同；更换参照改变参照水平

mk_task <- function(seed = 1, nv = 60, nt = 55, cohort = "C1", split = "s1") {
  set.seed(seed)
  val  <- cbind(Lasso = rnorm(nv), A = rnorm(nv), B = rnorm(nv),
                C = rnorm(nv), D = rnorm(nv))
  test <- cbind(Lasso = rnorm(nt), A = rnorm(nt), B = rnorm(nt),
                C = rnorm(nt), D = rnorm(nt))
  val[, "C"] <- 0                      # 常数验证预测（无信息）
  val[, "D"] <- val[, "A"]             # D 与 A 精确重复
  list(val_pred = val, test_pred = test,
       time_val = rexp(nv), event_val = rbinom(nv, 1, .6),
       time_test = rexp(nt), event_test = rbinom(nt, 1, .6),
       cohort = cohort, split = split, val_ids = paste0(cohort, split, "p", seq_len(nv)),
       test_ids = paste0(cohort, split, "q", seq_len(nt)),
       name = paste0(cohort, "|", split))
}

test_that("Orel 恒等于 Gv − Gt，且 none 规则为 0", {
  r <- audit_predictions(mk_task(1), reference = "Lasso")
  expect_s3_class(r, "sigAudit_prediction")
  d <- r$rules[r$rules$evaluable, ]
  expect_equal(d$Gv, d$Gt + d$Orel, tolerance = 1e-12)
  z <- d[d$rule == "none", ]
  expect_equal(z$Gv, 0, tolerance = 1e-12)
  expect_equal(z$Orel, 0, tolerance = 1e-12)
})

test_that("top_k 真正生效：top_k=4 生成并计算 top4_ensemble", {
  r <- audit_predictions(mk_task(2), reference = "Lasso", top_k = 4L)
  expect_true("top4_ensemble" %in% r$rules_run)
  expect_true("top4_ensemble" %in% unique(r$rules$rule))
  expect_equal(unique(r$diagnostics$topk), 4L)
  expect_true("top5_ensemble" %in% r$rules_run)     # 2/3/5/10 仍在
  expect_error(audit_predictions(mk_task(2), top_k = 0L), "positive integer")
  expect_error(audit_predictions(mk_task(2), top_k = 2.5), "positive integer")
})

test_that("家族规则可请求；未知规则报错", {
  tk <- mk_task(3)
  fam <- c(Lasso = "linear", A = "linear", B = "tree", C = "tree", D = "linear")
  r <- audit_predictions(tk, reference = "Lasso",
                         rules = c("linear_only", "tree_only"), family = fam)
  expect_true(all(c("linear_only", "tree_only") %in% unique(r$rules$rule)))
  expect_error(audit_predictions(tk, reference = "Lasso", rules = "not_a_rule"),
               "unknown rule")
  expect_error(audit_predictions(tk, reference = "Lasso", rules = "linear_only"),
               "family")
})

test_that("家族内只有一个候选时返回该候选（与原 rules_one 一致）", {
  tk <- mk_task(31)
  fam <- c(Lasso = "linear", A = "linear", B = "linear", C = "linear", D = "tree")
  r <- audit_predictions(tk, reference = "Lasso",
                         rules = c("linear_only", "tree_only"), family = fam)
  d <- r$rules[r$rules$evaluable, ]
  expect_equal(d$winner[d$rule == "tree_only"], "D")
  expect_true(d$evaluable[d$rule == "tree_only"])
  expect_false(is.na(d$winner[d$rule == "linear_only"]))
})

test_that("没有任何可评估记录时给出明确错误，而不是空表或 aggregate 报错", {
  tk <- mk_task(32)
  tk$time_val[] <- NA_real_          # 任务级掩码后验证侧为空 -> 所有 Cv 为 NA
  expect_error(audit_predictions(tk, reference = "Lasso"),
               "no evaluable task-rule records")
})

test_that("默认规则不含 oracle_test；显式请求时才有", {
  r <- audit_predictions(mk_task(4), reference = "Lasso")
  expect_false("oracle_test" %in% r$rules_run)
  r2 <- audit_predictions(mk_task(4), reference = "Lasso", rules = "oracle_test")
  expect_true("oracle_test" %in% r2$rules_run)
})

test_that("by_cohort/summary 用先队列×划分、再队列，而非任务层平均", {
  ## 同一队列两个划分：s1 有 1 个任务，s2 有 3 个任务
  t1 <- mk_task(11, cohort = "C1", split = "s1")
  t2 <- mk_task(12, cohort = "C1", split = "s2")
  t3 <- mk_task(13, cohort = "C1", split = "s2")
  t4 <- mk_task(14, cohort = "C1", split = "s2")
  r <- audit_predictions(list(t1, t2, t3, t4), reference = "Lasso",
                         rules = c("none", "argmax_val"))
  d <- r$rules[r$rules$rule == "argmax_val", ]
  gv_task_weighted <- mean(d$Gv)
  gv_split_then_cohort <- mean(tapply(d$Gv, d$split, mean))
  expect_false(isTRUE(all.equal(gv_task_weighted, gv_split_then_cohort)))   # 两者确实不同
  bc <- r$by_cohort[r$by_cohort$cohort == "C1" & r$by_cohort$rule == "argmax_val", ]
  expect_equal(bc$Gv, gv_split_then_cohort, tolerance = 1e-12)              # 队列层是划分等权
  s <- r$summary[r$summary$rule == "argmax_val", ]
  expect_equal(s$Gv, gv_split_then_cohort, tolerance = 1e-12)
})

test_that("输入契约：事件长度/取值、候选名唯一、ID 规则", {
  tk <- mk_task(5)
  bad <- tk; bad$event_val <- bad$event_val[-1]
  expect_error(audit_predictions(bad, reference = "Lasso"), "event_val")
  bad2 <- tk; bad2$event_val[1] <- 2
  expect_error(audit_predictions(bad2, reference = "Lasso"), "0/1")
  bad3 <- tk; colnames(bad3$val_pred) <- c("Lasso","Lasso","B","C","D")
  expect_error(audit_predictions(bad3, reference = "Lasso"), "duplicated candidate")
  bad4 <- tk; bad4$val_ids <- c(bad4$test_ids[1], bad4$val_ids[-1])
  bad5 <- tk; bad5$val_ids[1] <- bad5$test_ids[1]   # 与 test 重叠
  expect_error(audit_predictions(bad5, reference = "Lasso"), "overlap")
  expect_silent(audit_predictions(tk, reference = "Lasso"))
})

test_that("任务级共同掩码：丢弃被记录，不按候选各自删缺失", {
  tk <- mk_task(6, nv = 60, nt = 50)
  tk$val_pred[3, "A"] <- NA          # A 有 1 个缺失
  tk$time_val[8] <- NA               # 时间缺失
  r <- audit_predictions(tk, reference = "Lasso")
  dg <- r$diagnostics
  expect_equal(dg$n_drop_val, 2)     # 第 3（A 缺失）与第 8（时间缺失）位被丢弃
  expect_true(all(is.finite(r$rules$Cv[r$rules$rule == "argmax_val"])))
  expect_true(r$coverage[["n_drop_val"]] >= 2)
})

test_that("验证与测试样本数可以不同；缺少参照会报错", {
  tk <- mk_task(7, nv = 70, nt = 41)
  expect_silent(audit_predictions(tk, reference = "Lasso"))
  expect_error(audit_predictions(tk, reference = "NOPE"), "reference")
})

test_that("去重：唯一化、代表按名称字典序（A 保留、D 丢弃）", {
  tk <- mk_task(8)
  r <- audit_predictions(tk, reference = "Lasso", deduplicate = TRUE,
                         rules = c("none", "argmax_val", "top5_ensemble"))
  dg <- r$diagnostics
  expect_equal(dg$n_unique, 4)       # Lasso, A, B, C
  mem <- r$rules$members[r$rules$rule == "top5_ensemble"]
  expect_true(grepl("A", mem))
  expect_false(grepl("D", mem))
  expect_true("none" %in% r$rules$rule)
})

test_that("更换固定参照改变参照水平与 gain 口径", {
  tk <- mk_task(9)
  rL <- audit_predictions(tk, reference = "Lasso")
  rA <- audit_predictions(tk, reference = "A")
  gL <- rL$rules[rL$rules$rule == "argmax_val", ]
  gA <- rA$rules[rA$rules$rule == "argmax_val", ]
  expect_false(isTRUE(all.equal(gA$Gt, gL$Gt, tolerance = 1e-12)))
  expect_equal(rA$rules$Orel[rA$rules$rule == "none"], 0, tolerance = 1e-12)
})

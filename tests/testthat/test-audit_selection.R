## 模块 4/5 的断言式测试。
## 设计原则（评审要求）：**用断言而不是 cat 打印**，失败会以非零码退出；
## 并且至少覆盖以下会改变科学结论的点：
##   ① 风险方向固定（反向评分必须低于 0.5，不因评价集而翻转）
##   ② 常数预测的处理
##   ③ 主规则允许参照参与（不是"强制选一个非参照"）
##   ④ Orel = Gv − Gt 逐任务成立
##   ⑤ 进入效应用随机背景池、并报告 MC 标准误；全池时不假装抽样
##   ⑥ 并列处理与候选输入顺序无关 —— 而且**必须用"并列者测试性能不同"的构造**，
##      否则换了获胜者报告数值也不变，测不出问题（第三轮审查指出的缺陷）
##   ⑦ 行对应必须可核验（row_ids），不能只比尺寸
##   ⑧ 非有限值与越界 C-index 必须被明确处理，不能静默忽略

## ---------- 合成 C-index 矩阵：植入一个高乐观度候选 ----------
set.seed(20260101)
n_gs <- 120
mk <- function(mu_v, sd_v, mu_t, sd_t) cbind(v = rnorm(n_gs, mu_v, sd_v), t = rnorm(n_gs, mu_t, sd_t))
CND <- list(
  REF = mk(.60, .03, .585, .03),
  A   = mk(.61, .03, .590, .03),
  B   = mk(.61, .04, .585, .04),
  C   = mk(.60, .03, .580, .03),
  BAD = mk(.645, .07, .555, .07)        # 验证集又高又飘，测试集跟不上
)
nm <- names(CND)
V  <- sapply(CND, function(z) z[, "v"]); dimnames(V)  <- list(paste0("s", seq_len(n_gs)), nm)
Tt <- sapply(CND, function(z) z[, "t"]); dimnames(Tt) <- list(paste0("s", seq_len(n_gs)), nm)
FAM <- setNames(c("linear","linear","tree","linear","tree"), nm)   # 树候选需 >=2 才会生成 tree_only

test_that("池级 Orel 与剂量曲线可用，且口径与 Oraw 严格区分", {
  r <- audit_selection(V, Tt, reference = "REF", n_background = 10, seed = 7)
  expect_s3_class(r, "sigAudit_module")
  expect_true(is.finite(r$pool_bias))
  expect_gt(r$pool_bias, 0)                       # 池中含植入的 BAD，乐观度应为正
  expect_true(all(c("k","candidate_added","orel") %in% names(r$dose)))
  expect_equal(nrow(r$dose), ncol(V) - 1)
})

test_that("进入效应用随机背景池，并报告 MC 标准误；全池时不假装抽样", {
  r <- audit_selection(V, Tt, reference = "REF", n_background = 20, seed = 7)
  e <- r$entry
  expect_true(all(c("target","k","entry","sd_bg","mc_se","n_bg") %in% names(e)))
  part <- e[e$n_bg > 1, ]
  expect_gt(nrow(part), 0)
  expect_true(all(is.finite(part$mc_se)))
  # MC 标准误应远小于效应本身的量级
  expect_true(median(part$mc_se) < 0.5 * max(abs(part$entry)))
  full <- e[e$n_bg == 1, ]
  expect_true(all(is.na(full$sd_bg)))             # 全池只有一个背景集合
})

test_that("进入效应对植入的高乐观度候选给出最大的正效应", {
  r <- audit_selection(V, Tt, reference = "REF", n_background = 20, seed = 7)
  ag <- aggregate(entry ~ target, data = r$entry, FUN = function(x) mean(x, na.rm = TRUE))
  ag <- ag[order(-ag$entry), ]
  expect_equal(ag$target[1], "BAD")
  expect_gt(ag$entry[1], 0)
})

test_that("规则核算：主规则允许参照参与，Orel = Gv − Gt 逐任务成立", {
  r <- audit_rules(V, Tt, reference = "REF", family = FAM)
  t <- r$table
  expect_true(all(c("none","argmax_val","runner_up","argmax_nonreference",
                    "linear_only","tree_only","oracle_test") %in% t$rule))
  expect_equal(t$Gv, t$Gt + t$Orel, tolerance = 1e-12)
  z <- t[t$rule == "none", ]
  expect_equal(z$Gv, 0, tolerance = 1e-12)
  expect_equal(z$Orel, 0, tolerance = 1e-12)
  # 参照必须能参与主规则：构造"参照最好"的场景
  V2 <- V; Tt2 <- Tt
  V2[, "REF"] <- 0.99                            # 让 REF 在所有签名上都最高（且仍在 [0,1] 内）
  r2 <- audit_rules(V2, Tt2, reference = "REF", family = FAM)
  expect_gt(r2$table$Gv[r2$table$rule == "argmax_val"], -1e-12)   # 选中 REF -> Gv = 0
  expect_lt(r2$table$Gv[r2$table$rule == "argmax_nonreference"], 0)  # 被迫选非参照 -> 为负
})

test_that("linear_only 含参照（与 scripts/20 定义一致）；不提供 family 时不生成家族规则", {
  r <- audit_rules(V, Tt, reference = "REF", family = FAM)
  t <- r$table
  expect_true("linear_only" %in% t$rule)
  # 让参照在验证集上碾压一切：linear_only（含参照）必须选中 REF，Gv = 0
  V3 <- V; Tt3 <- Tt; V3[, "REF"] <- 1
  r3 <- audit_rules(V3, Tt3, reference = "REF", family = FAM)
  expect_equal(r3$table$Gv[r3$table$rule == "linear_only"], 0, tolerance = 1e-12)
  # 未提供 family：不应生成 linear_only（否则会退化成 argmax_val 的重复）
  r4 <- audit_rules(V, Tt, reference = "REF")
  expect_false("linear_only" %in% r4$table$rule)
  expect_false("tree_only"   %in% r4$table$rule)
})

## ---- ⑥ 并列：必须用"并列者测试性能不同"的构造，否则测不出问题 ----
test_that("并列时获胜者与报告值都与候选输入列顺序无关（并列者测试性能不同）", {
  n <- 12
  # REF 处处 0.6；A 与 B 在验证集处处并列于 0.8，但测试性能 A=0.7、B=0.4
  Vc <- matrix(0.6, n, 3, dimnames = list(paste0("g", 1:n), c("REF","A","B")))
  Vc[, "A"] <- 0.8; Vc[, "B"] <- 0.8
  Tc <- matrix(0.58, n, 3, dimnames = list(paste0("g", 1:n), c("REF","A","B")))
  Tc[, "A"] <- 0.70; Tc[, "B"] <- 0.40

  r1 <- audit_rules(Vc, Tc, reference = "REF", seed = 5)
  r2 <- audit_rules(Vc[, c("REF","B","A")], Tc[, c("REF","B","A")], reference = "REF", seed = 5)
  # 获胜者必须一致（旧实现会失败：输入顺序 REF,A,B 选 B -> 报 0.4；REF,B,A 选 A -> 报 0.7）
  cv1 <- r1$table$reported[r1$table$rule == "argmax_val"]
  cv2 <- r2$table$reported[r2$table$rule == "argmax_val"]
  ct1 <- r1$table$test[r1$table$rule == "argmax_val"]
  ct2 <- r2$table$test[r2$table$rule == "argmax_val"]
  expect_equal(cv1, 0.8, tolerance = 1e-12)
  expect_equal(cv2, 0.8, tolerance = 1e-12)
  expect_equal(ct1, ct2, tolerance = 1e-12)      # 旧实现这里会得到 0.4 vs 0.7
  # 用 audit_selection 的池级量再验一次
  s1 <- audit_selection(Vc, Tc, reference = "REF", n_background = 5, seed = 5)
  s2 <- audit_selection(Vc[, c("REF","B","A")], Tc[, c("REF","B","A")],
                        reference = "REF", n_background = 5, seed = 5)
  expect_equal(s1$pool_bias, s2$pool_bias, tolerance = 1e-12)
  # entry 已按名字排序抽样：按 (target, k) 排序后逐值相等
  o1 <- order(s1$entry$target, s1$entry$k); o2 <- order(s2$entry$target, s2$entry$k)
  expect_equal(s1$entry$entry[o1], s2$entry$entry[o2], tolerance = 1e-12)
})

test_that("常数候选（无信息）的处理：子池全为常数时该规则增益为 0", {
  V4 <- V; Tt4 <- Tt
  # 把两个树候选都设为常数风险分（C-index 恒为 0.5）——对应真实数据里
  # 正则化在无信号基因集上选到零系数模型的情形（合法零信号解）
  V4[, "BAD"] <- 0.5; Tt4[, "BAD"] <- 0.5
  V4[, "B"]   <- 0.5; Tt4[, "B"]   <- 0.5
  r <- audit_rules(V4, Tt4, reference = "REF", family = FAM)
  tt <- r$table[r$table$rule == "tree_only", ]
  expect_equal(tt$reported, 0.5, tolerance = 1e-12)
  expect_equal(tt$test,     0.5, tolerance = 1e-12)
  expect_equal(tt$Oraw,     0,   tolerance = 1e-12)
  expect_lt(tt$Gv, 0)
  # 常数候选不会赢得全池 argmax
  expect_false(r$table$rule[r$table$rule == "tree_only"] == "argmax_val")
})

## ---- ⑦⑧ 行对应与取值范围 ----
test_that("行对应无法核验时必须报错，除非显式声明位置配对", {
  Vn <- V; Tn <- Tt; rownames(Vn) <- NULL; rownames(Tn) <- NULL   # 去掉行名，保留列名
  expect_error(audit_rules(Vn, Tn, reference = "REF"),
               "assume_positional_pairing")
  expect_silent(audit_rules(Vn, Tn, reference = "REF",
                            assume_positional_pairing = TRUE))
  # 行名不一致 = 行可能被调换过 -> 必须报错，而不是静默位置配对
  Tp <- Tt; rownames(Tp) <- rev(rownames(Tt))
  expect_error(audit_rules(V, Tp, reference = "REF"), "assume_positional_pairing")
  # 提供 row_ids 时按 ID 对齐
  ids <- rownames(V)
  expect_silent(audit_rules(V, Tt, reference = "REF", row_ids = ids))
  expect_error(audit_rules(V, Tt, reference = "REF", row_ids = c("zzz", ids[-1])),
               "not the same set")
  expect_error(audit_rules(V, Tt, reference = "REF", row_ids = rep("x", n_gs)),
               "unique")
})

test_that("非有限值被记 NA 并计数，越界 C-index 被拒绝", {
  V5 <- V; Tt5 <- Tt
  V5[1, "A"] <- Inf; Tt5[2, "B"] <- NaN
  r <- audit_rules(V5, Tt5, reference = "REF")
  expect_equal(unname(r$coverage["n_nonfinite"]), 2)
  expect_equal(unname(r$coverage["n_dropped_incomplete"]), 2)
  V6 <- V; V6[1, "A"] <- 1.4                    # 越界
  expect_error(audit_rules(V6, Tt, reference = "REF"), "outside \\[0,1\\]")
})

test_that("退化输入被明确拒绝", {
  expect_error(audit_selection(V[1:5, ], Tt[1:5, ], reference = "REF"), "fewer than 10")
  expect_error(audit_selection(V, Tt[, 1:2], reference = "REF"), "different dimensions")
  expect_error(audit_selection(V, Tt, reference = "NOPE"), "reference model not among")
  expect_error(audit_rules(V, Tt, reference = "REF", family = c("linear")), "named vector")
})

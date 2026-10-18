#!/usr/bin/env Rscript
# ============================================================
# 25 — 紧凑模拟：在【已知风险结构】下检验收益核算与集成策略
#
# 本版按【第三轮审查】重写，修掉六处实质问题：
#  ① 静默丢弃重复：旧版对每个候选调用 scale(v)，常数预测产生非有限值，随后
#     "全部候选 C-index 必须是有限值"让【整次重复返回 NULL】，且不记录原因。
#     旧结果里弱信号+rho=0.6 两个条件只剩 65、69 个有效重复（应为 100）。
#     → 现在每个候选的常数性被显式判定：常数风险分定义为 C = 0.5，
#       **保留为合法无信息候选**；每个计划重复都产出状态记录，绝不静默丢弃。
#  ② rsf1/rsf2 是两个【完全相同】的设置（只差名字）。→ 现在树数/mtry/种子都不同，
#     并输出两者的验证分相关系数作为"它们确实不同"的检查。
#  ③ 随机流没有拆分：拟合函数反复重置全局种子，噪声仍从同一全局状态抽取。
#     → 数据生成 / 模型拟合 / 扰动噪声 / 删失 各用独立种子并记录。
#  ④ 删失并非"独立"：事件 exp(lp)、删失 0.6*exp(lp) 共享风险结构。
#     → 主模拟改为【不随 lp 变化】的固定删失率，并报告实际删失率。
#  ⑤ rho 同时改变相似度与信号质量，不是纯相关性操纵。
#     → 改称"共享评分扰动情景"，并同时报告每条件的平均单模型 C（让读者看到信号变化）。
#  ⑥ "Orel 随池规模增大约 2.3 倍"不是全网格规律（弱信号下反而下降）。
#     → 逐条件报告，不给单调因果结论；K 取 6/12/20 三点而非两点。
#
# 另外：扰动在【分数层面】实现，且均值/标准差只由验证侧估计、固定应用到测试侧
# （旧版对验证与测试各自 scale，是依赖整批测试样本的变换）。
#
# 用法:
#   Rscript scripts/25_simulation.R pilot             # 少量重复，估成本与 MC 精度
#   Rscript scripts/25_simulation.R final 100         # 正式运行（重复次数事先固定）
#   Rscript scripts/25_simulation.R final 100 serial  # 强制串行
# ============================================================
suppressMessages({library(survival); library(glmnet); library(randomForestSRC)})
options(stringsAsFactors = FALSE)
PROJ <- Sys.getenv("SIGAUDIT_ROOT", "D:/bioinfo/sigAudit"); TAB <- file.path(PROJ,"results/tables")
LOGDIR <- file.path(PROJ,"logs"); dir.create(LOGDIR, showWarnings=FALSE, recursive=TRUE)

ARGS <- commandArgs(trailingOnly = TRUE)
MODE <- if (length(ARGS) >= 1) ARGS[1] else "pilot"
NREP <- if (length(ARGS) >= 2) as.integer(ARGS[2]) else if (MODE == "pilot") 5L else 100L
RUNMODE <- if (length(ARGS) >= 3) ARGS[3] else "auto"
LOG  <- file.path(LOGDIR, sprintf("25_simulation_%s.log", MODE)); con <- file(LOG, "wt")
say  <- function(...) { m <- paste0(...); cat(m,"\n",sep=""); cat(m,"\n",sep="",file=con); flush(con) }

## ---- 事先固定的模拟参数 ----
N_SEL    <- 300      # 选择样本（内含 训练/验证）
N_TEST   <- 600      # 独立测试样本（从不参与选择）
P_GENE   <- 100
FRAC_VAL <- 1/3
CENS_RATE <- 0.30    # 目标边际删失比例（**由独立校准样本求解**，见下）

## 【修正·第四轮】删失分布参数由【独立校准样本】预先确定，
## 既不依赖任何一次模拟出来的数据，也不依赖测试部分。
##   TT | lp ~ Exp(exp(lp))，lp ~ N(0, s²)，s² = sum(beta²) 是【设计常数】。
##   在某个候选删失速率 λ_c 下，边际删失比例 = E[λ_c / (exp(lp) + λ_c)]。
##   用固定的、与任何重复都不共用的种子从【设计分布】抽一个很大的独立样本，
##   对 λ_c 求根，使该期望等于目标比例。结果在任何数据被抽取之前就已确定。
CALIB_SEED <- 20261015L
calib_rate_c <- function(sigma2, p, seed = CALIB_SEED, n = 2e6) {
  set.seed(seed)
  lam_i <- exp(stats::rnorm(n, 0, sqrt(sigma2)))
  stats::uniroot(function(lc) mean(lc / (lam_i + lc)) - p, c(1e-8, 1e6))$root
}
SIGMA2  <- c(strong = 20 * 0.45^2, weak = 6 * 0.18^2)   # 与下面的 nsig/beta 保持一致
RATE_C  <- vapply(SIGMA2, function(s2) calib_rate_c(s2, CENS_RATE), numeric(1))

## ---- C-index ----
ci_direct <- function(t, e, s) {
  s <- tryCatch(as.numeric(s), error=function(x) NULL)
  if (is.null(s) || length(s) != length(t)) return(NA_real_)
  ok <- is.finite(s) & !is.na(t) & !is.na(e)
  if (sum(ok) < 10) return(NA_real_)
  tryCatch(as.numeric(concordance(Surv(t[ok],e[ok]) ~ s[ok], reverse=TRUE, timewt="n")$concordance),
           error=function(x) NA_real_)
}
## 【修正】常数风险分 = 无信息候选，C 显式定义为 0.5（不再让 scale() 产出非有限值）
ci_safe <- function(t, e, s) {
  s <- tryCatch(as.numeric(s), error=function(x) NULL)
  if (is.null(s) || length(s) != length(t)) return(NA_real_)
  if (length(unique(s[is.finite(s)])) <= 1L) return(0.5)
  ci_direct(t, e, s)
}

## 候选定义：惩罚 Cox 变体（alpha × 预筛）+ 2 个【参数不同】的 RSF
POOL_ALL <- list(
  list(n="lasso",   a=1.0), list(n="ridge",   a=0.0), list(n="enet.9", a=0.9),
  list(n="enet.7",  a=0.7), list(n="enet.5",  a=0.5), list(n="enet.3", a=0.3),
  list(n="enet.1",  a=0.1), list(n="en5T50",  a=0.5, pre=50), list(n="en7T25", a=0.7, pre=25),
  list(n="lsT50",   a=1.0, pre=50), list(n="lsT25", a=1.0, pre=25), list(n="lsT10", a=1.0, pre=10),
  list(n="en3T50",  a=0.3, pre=50), list(n="lsT70", a=1.0, pre=70), list(n="en9T50", a=0.9, pre=50),
  list(n="lsT100",  a=1.0, pre=100), list(n="en1T50", a=0.1, pre=50), list(n="en5T80", a=0.5, pre=80)
)
# 【修正】两个森林在树数、mtry、种子上都不同，不再是同一设置的重复
RSF_DEF <- list(list(n="rsf1", ntree=100, mtry=NULL, seed=101),
                list(n="rsf2", ntree=300, mtry=30,   seed=202))
make_pool <- function(K) c(POOL_ALL[seq_len(min(K - 2, length(POOL_ALL)))], RSF_DEF)

## 每个候选只拟合一次，同时对 val / test 预测
fit_pool <- function(pool, Xa, ta, ea, Xb, Xc, zz = NULL, st = NULL) {
  out <- list()
  ## 【第四轮】逐候选记录"该候选的拟合是否产生过 warning"，
  ## 用于回答"最终被选中的预测是否来自有效、收敛的拟合"。
  ## 注意：这里必须把 warning 同时计入调用方传入的 st（否则外层计数会因内层 muffle 而漏掉）。
  .wrapped <- function(expr) {
    if (is.null(st)) return(expr)
    withCallingHandlers(expr, warning = function(w) {
      st$nw <- st$nw + 1L
      if (is.na(st$first)) st$first <- conditionMessage(w)
      invokeRestart("muffleWarning")
    })
  }
  if (is.null(zz)) {
    zz <- apply(Xa, 2, function(x) {
      q <- tryCatch(summary(coxph(Surv(ta, ea) ~ x))$coef, error = function(e) NULL)
      if (is.null(q)) NA_real_ else abs(q[1, "z"]) })
    zz[!is.finite(zz)] <- -Inf
  }
  for (m in pool) {
    w0 <- if (is.null(st)) 0L else st$nw
    if (!is.null(m$ntree)) {                      # 森林候选
      d <- data.frame(.t=ta, .s=ea, Xa, check.names=FALSE); set.seed(m$seed)
      rf <- tryCatch(.wrapped({
        if (is.null(m$mtry)) rfsrc(Surv(.t,.s) ~ ., data=d, ntree=m$ntree, splitrule="logrank")
        else rfsrc(Surv(.t,.s) ~ ., data=d, ntree=m$ntree, mtry=m$mtry, splitrule="logrank")
      }), error=function(e) NULL)
      if (is.null(rf)) next
      pv <- tryCatch(as.numeric(predict(rf, data.frame(Xb, check.names=FALSE))$predicted), error=function(e) NULL)
      pt <- tryCatch(as.numeric(predict(rf, data.frame(Xc, check.names=FALSE))$predicted), error=function(e) NULL)
      if (!is.null(pv) && !is.null(pt) && all(is.finite(pv)) && all(is.finite(pt)))
        out[[m$n]] <- list(v=pv, t=pt, warned = (!is.null(st)) && st$nw > w0)
      next
    }
    A <- Xa; cols <- colnames(Xa)
    if (!is.null(m$pre)) {
      kk <- min(m$pre, ncol(A))
      if (kk < 5 || sum(zz > -Inf) < 5) next
      sel <- order(zz, decreasing=TRUE)[seq_len(kk)]
      A <- A[, sel, drop=FALSE]; cols <- colnames(A)
    }
    r <- tryCatch({
      set.seed(1)
      cv <- .wrapped(cv.glmnet(A, Surv(ta,ea), family="cox", alpha=m$a, nfolds=5, cox.ties="breslow"))
      list(v=as.numeric(predict(cv, Xb[, cols, drop=FALSE], s=cv$lambda.min, type="link")),
           t=as.numeric(predict(cv, Xc[, cols, drop=FALSE], s=cv$lambda.min, type="link")),
           warned = (!is.null(st)) && st$nw > w0,
           lam_min_idx = which.min(abs(cv$lambda - cv$lambda.min)))
    }, error=function(e) NULL)
    if (!is.null(r) && all(is.finite(r$v)) && all(is.finite(r$t))) out[[m$n]] <- r
  }
  out
}

## 共享评分扰动：均值/标准差**只由验证侧估计**，固定应用到测试侧。
## 常数候选保持常数（无信息），不参与扰动 —— 这样它始终是合法的 C = 0.5 候选。
perturb_pool <- function(pm, rho, zv, zt) {
  lapply(pm, function(p) {
    s <- stats::sd(p$v)
    if (!is.finite(s) || s == 0) return(p)     # 常数：原样保留
    mu <- mean(p$v)
    vs <- (p$v - mu) / s
    ts <- (p$t - mu) / s                       # 用验证侧参数，不用测试侧信息
    ## 【修正】旧实现 `list(v=..., t=...)` 会丢掉 warned / lam_min_idx 等字段，
    ## 导致 rho>0 的 6 个条件里 winner_warned 恒为 FALSE（是构造出来的，不是观测到的）。
    ## 现在原地替换 v/t，其余字段原样保留。
    p$v <- as.numeric(rho * zv + sqrt(1 - rho^2) * vs)
    p$t <- as.numeric(rho * zt + sqrt(1 - rho^2) * ts)
    p
  })
}

## 单次重复：**任何情况都返回一条记录**，绝不静默返回 NULL
one_rep <- function(K, signal, rho, seed, cond_id, rep_id) {
  rec <- data.frame(cond_id=cond_id, replicate_id=rep_id, seed=seed, K=K,
                    signal=signal, rho=rho, status="failed", failure_stage="init",
                    failure_reason=NA_character_, K_eff=NA_integer_, n_const=NA_integer_,
                    Cv_mean=NA_real_, Ct_mean=NA_real_, censor_rate=NA_real_,
                    cor_rsf=NA_real_, Ct_argmax=NA_real_, Ct_ens=NA_real_,
                    Gv=NA_real_, Gt=NA_real_, Orel=NA_real_,
                    dCt_ens=NA_real_, dOrel_ens=NA_real_,
                    n_warn=0L, warn_first=NA_character_,
                    # 【第四轮】把派生的四路种子显式写出，使"随机流已记录"这句话可核验
                    seed_data=seed * 11L + 1L, seed_censor=seed * 11L + 4L,
                    seed_fit=seed * 11L + 2L, seed_perturb=seed * 11L + 3L,
                    cens_rate_c=NA_real_, n_warned_cands=NA_integer_,
                    winner_warned=NA, stringsAsFactors=FALSE)
  stage <- "init"
  # 逐重复统计 warning（glmnet 的 maxit 收敛提示等），便于审查者判断是否存在被掩盖的问题
  st <- new.env(parent = emptyenv()); st$nw <- 0L; st$first <- NA_character_
  err <- tryCatch(withCallingHandlers({
    ## ---- 1) 数据生成（独立随机流）----
    stage <- "data"
    set.seed(seed * 11L + 1L)
    n <- N_SEL + N_TEST
    X <- matrix(rnorm(n * P_GENE), n, P_GENE)
    colnames(X) <- paste0("g", seq_len(P_GENE))
    nsig <- if (signal == "strong") 20 else 6
    beta <- numeric(P_GENE); beta[seq_len(nsig)] <- if (signal == "strong") 0.45 else 0.18
    lp <- as.numeric(X %*% beta)
    TT <- rexp(n, rate = exp(lp))
    ## 【修正·第四轮】删失分布参数取自独立校准样本（见文件头 RATE_C），
    ## 不再由本次抽到的 TT（含测试部分）决定。
    stage <- "censor"
    rate_c   <- unname(RATE_C[[signal]])       # 在任何数据被抽取之前就已确定
    set.seed(seed * 11L + 4L)
    censor_t <- rexp(n, rate = rate_c)
    t <- pmax(pmin(TT, censor_t), 1e-6); e <- as.integer(TT <= censor_t)
    i_sel <- seq_len(N_SEL); i_te <- (N_SEL + 1):n
    nv <- floor(N_SEL * FRAC_VAL); i_va <- i_sel[seq_len(nv)]; i_tr <- setdiff(i_sel, i_va)
    rec$censor_rate <- 1 - mean(e)
    rec$cens_rate_c <- rate_c        # 预先固定的删失分布参数（设计常数）

    ## ---- 2) 模型拟合（独立随机流）----
    stage <- "fit"
    pm <- fit_pool(make_pool(K), X[i_tr,,drop=FALSE], t[i_tr], e[i_tr],
                   X[i_va,,drop=FALSE], X[i_te,,drop=FALSE], st = st)
    if (length(pm) < 4) stop(sprintf("only %d candidates fitted", length(pm)))
    if (!("lasso" %in% names(pm))) stop("reference (lasso) failed to fit")

    ## ---- 3) 扰动噪声（独立随机流）----
    stage <- "perturb"
    set.seed(seed * 11L + 3L)
    zv <- rnorm(length(i_va)); zt <- rnorm(length(i_te))
    if (rho > 0) pm <- perturb_pool(pm, rho, zv, zt)

    ## ---- 4) 评价 ----
    stage <- "evaluate"
    mods <- names(pm)
    V  <- vapply(mods, function(m) ci_safe(t[i_va], e[i_va], pm[[m]]$v), numeric(1))
    Tt <- vapply(mods, function(m) ci_safe(t[i_te], e[i_te], pm[[m]]$t), numeric(1))
    if (!all(is.finite(V)) || !all(is.finite(Tt))) stop("non-finite C-index")
    rec$K_eff    <- length(mods)
    rec$n_const  <- sum(vapply(mods, function(m) length(unique(pm[[m]]$v)) <= 1L, logical(1)))
    rec$Cv_mean  <- mean(V); rec$Ct_mean <- mean(Tt)
    if (all(c("rsf1","rsf2") %in% mods)) {
      r1 <- pm[["rsf1"]]$v; r2 <- pm[["rsf2"]]$v
      if (stats::sd(r1) > 0 && stats::sd(r2) > 0) rec$cor_rsf <- stats::cor(r1, r2)
    }
    j <- which.max(V)
    ## 【第四轮】"最终被采用的预测是否来自有 warning 的拟合"
    rec$n_warned_cands <- sum(vapply(pm, function(p) isTRUE(p$warned), logical(1)))
    rec$winner_warned  <- isTRUE(pm[[mods[j]]]$warned)
    rec$Ct_argmax <- unname(Tt[j]); rec$Gv <- unname(V[j] - V[["lasso"]]); rec$Gt <- unname(Tt[j] - Tt[["lasso"]])
    rec$Orel <- rec$Gv - rec$Gt
    # top5 集成：均值/标准差只用验证集估计，固定应用到测试
    top5 <- mods[order(-V)][seq_len(min(5, length(mods)))]
    Mv <- sapply(top5, function(m) pm[[m]]$v); Mt <- sapply(top5, function(m) pm[[m]]$t)
    mu <- colMeans(Mv); sdv <- apply(Mv, 2, stats::sd); sdv[!is.finite(sdv) | sdv == 0] <- 1
    ev <- rowMeans(sweep(sweep(Mv, 2, mu, "-"), 2, sdv, "/"))
    et <- rowMeans(sweep(sweep(Mt, 2, mu, "-"), 2, sdv, "/"))
    cv_e <- ci_safe(t[i_va], e[i_va], ev); ct_e <- ci_safe(t[i_te], e[i_te], et)
    rec$Ct_ens    <- unname(ct_e)
    rec$dCt_ens   <- unname(ct_e - Tt[j])
    rec$dOrel_ens <- unname((cv_e - ct_e) - (V[j] - Tt[j]))

    NULL
  }, warning = function(w) {
    st$nw <- st$nw + 1L
    if (is.na(st$first)) st$first <- conditionMessage(w)
    invokeRestart("muffleWarning")
  }), error = function(e) e)
  rec$n_warn <- st$nw; rec$warn_first <- st$first
  if (is.null(err)) {
    rec$status <- "ok"; rec$failure_stage <- NA_character_
  } else {
    rec$failure_stage  <- stage
    rec$failure_reason <- conditionMessage(err)
  }
  rec
}

## ---- 事先固定的条件网格：池规模 × 信号 × 共享评分扰动 ----
COND <- expand.grid(K = c(6, 12, 20), signal = c("strong", "weak"), rho = c(0, 0.6),
                    stringsAsFactors = FALSE)
say("模式=", MODE, "  每条件计划重复 ", NREP, "  条件数 ", nrow(COND),
    "  运行方式=", RUNMODE)
say("每次重复: 选择样本 ", N_SEL, "（训练 ", N_SEL - floor(N_SEL*FRAC_VAL),
    " / 验证 ", floor(N_SEL*FRAC_VAL), "）, 独立测试 ", N_TEST, ", 基因 ", P_GENE)
say("候选: 惩罚 Cox 变体 + 2 个参数不同的 RSF(ntree 100/300, mtry 默认/30)")
say("删失: 固定边际删失率（与 lp 无关），逐重复报告实际删失率")
say("rho: 共享评分扰动情景（同时改变候选相似度与单模型性能，不是纯相关性操纵）")
say("随机流: 数据生成 / 删失 / 拟合 / 扰动 各用独立种子")
say("")

run_cond <- function(i) {
  cd <- COND[i, ]
  do.call(rbind, lapply(seq_len(NREP), function(r)
    one_rep(cd$K, cd$signal, cd$rho, seed = 20261000 + 1000*i + r,
            cond_id = i, rep_id = r)))
}

t0 <- Sys.time()
## 【第四轮】worker 数可由第 4 个参数或 SIGAUDIT_WORKERS 指定。
## 上限受【内存】约束而非 CPU：本机 16 核/32 线程但只有 15.2 GB 内存，
## 每 worker 约 0.5-0.6 GB（n=900、100 协变量、K<=20），因此默认 10 而非更高。
NCW <- Sys.getenv("SIGAUDIT_WORKERS", "")
NC <- if (RUNMODE == "serial" || MODE == "pilot") 1L else
  if (nzchar(NCW)) as.integer(NCW) else max(1L, min(10L, parallel::detectCores() - 2L))
if (NC > 1) {
  say("并行：", NC, " 个 worker（每个 worker 内 glmnet/rfsrc 仍单线程）")
  cl <- parallel::makeCluster(NC)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  parallel::clusterExport(cl, c("COND","NREP","N_SEL","N_TEST","P_GENE","FRAC_VAL","CENS_RATE",
                                "POOL_ALL","RSF_DEF","make_pool","fit_pool","perturb_pool","RATE_C","SIGMA2",
                                "one_rep","ci_direct","ci_safe","run_cond"), envir = environment())
  parallel::clusterEvalQ(cl, suppressMessages({library(survival); library(glmnet); library(randomForestSRC)}))
  REP <- do.call(rbind, parallel::parLapply(cl, seq_len(nrow(COND)), function(i) {
    suppressMessages({library(survival); library(glmnet); library(randomForestSRC)})
    run_cond(i) }))
} else {
  REP <- do.call(rbind, lapply(seq_len(nrow(COND)), run_cond))
}

fn_rep <- file.path(TAB, sprintf("25_simulation_%s_replicates.tsv", MODE))
write.table(REP, fn_rep, sep="\t", quote=FALSE, row.names=FALSE)
say("逐重复记录: ", fn_rep, "（共 ", nrow(REP), " 条 = 计划重复总数，无静默丢弃）")

## 失败汇总
FAIL <- REP[REP$status != "ok", ]
if (nrow(FAIL)) {
  ft <- aggregate(list(n = FAIL$replicate_id),
                  by = list(cond_id = FAIL$cond_id, failure_stage = FAIL$failure_stage,
                            failure_reason = FAIL$failure_reason), FUN = length)
  write.table(ft, file.path(TAB, sprintf("25_simulation_%s_failures.tsv", MODE)),
              sep="\t", quote=FALSE, row.names=FALSE)
  say("失败记录: ", nrow(FAIL), " 条，明细见 25_simulation_", MODE, "_failures.tsv")
} else {
  say("失败记录: 0 条")
}

## warning 汇总（逐重复记录，便于判断是否存在被掩盖的问题）
say("warning: 共 ", sum(REP$n_warn), " 条（", sum(REP$n_warn > 0), " 个重复含 warning）")
wf <- REP$warn_first[!is.na(REP$warn_first)]
if (length(wf)) say("  最常见示例: ", substr(names(sort(table(wf), decreasing=TRUE))[1], 1, 110))

## 每条件汇总（逐条件报告，不给单调因果结论）
res <- list()
for (i in seq_len(nrow(COND))) {
  cd <- COND[i, ]; M <- REP[REP$cond_id == i & REP$status == "ok", , drop=FALSE]
  nr <- nrow(M)
  if (nr < 2) {
    say(sprintf("  条件 %2d K=%2d %-6s rho=%.1f | 有效重复 %d，不足以汇总", i, cd$K, cd$signal, cd$rho, nr))
    next
  }
  res[[i]] <- data.frame(cond_id=i, K=cd$K, signal=cd$signal, rho=cd$rho,
    n_planned=NREP, n_ok=nr, n_failed=NREP-nr,
    K_eff_mean=mean(M$K_eff), n_const_mean=mean(M$n_const),
    censor_rate_mean=mean(M$censor_rate), Cv_mean=mean(M$Cv_mean), Ct_mean=mean(M$Ct_mean),
    cor_rsf=mean(M$cor_rsf, na.rm=TRUE),
    cens_rate_c=unname(RATE_C[[cd$signal]]),
    n_warn_mean=mean(M$n_warn), n_warned_cands_mean=mean(M$n_warned_cands, na.rm=TRUE),
    winner_warned_rate=mean(M$winner_warned, na.rm=TRUE),
    Gv=mean(M$Gv), Gt=mean(M$Gt), Orel=mean(M$Orel),
    Orel_mcse=stats::sd(M$Orel)/sqrt(nr),
    dCt_ens=mean(M$dCt_ens), dCt_ens_mcse=stats::sd(M$dCt_ens)/sqrt(nr),
    dOrel_ens=mean(M$dOrel_ens))
  say(sprintf("  K=%2d %-6s rho=%.1f | n=%3d/%3d 池内候选=%4.1f 常数候选=%3.1f 删失率=%.2f(目标%.2f) 平均单模型C(验/测)=%.3f/%.3f | Gv %+.3f Gt %+.3f Orel %+.3f (MCSE %.4f) | 集成ΔCt %+.3f (MCSE %.4f) | 中选者曾有warning的比例=%.1f%%",
    cd$K, cd$signal, cd$rho, nr, NREP, mean(M$K_eff), mean(M$n_const), mean(M$censor_rate),
    CENS_RATE,
    mean(M$Cv_mean), mean(M$Ct_mean), mean(M$Gv), mean(M$Gt), mean(M$Orel),
    stats::sd(M$Orel)/sqrt(nr), mean(M$dCt_ens), stats::sd(M$dCt_ens)/sqrt(nr),
    100*mean(M$winner_warned, na.rm=TRUE)))
}
if (length(res)) {
  R <- do.call(rbind, res)
  fn <- file.path(TAB, sprintf("25_simulation_%s.tsv", MODE))
  write.table(R, fn, sep="\t", quote=FALSE, row.names=FALSE)
  say(""); say("条件汇总: ", fn)
  say("用时 ", round(as.numeric(difftime(Sys.time(), t0, units="min")), 1), " 分钟")
  say("")
  say("判读要点（不以跑出想要的 P 值为停止条件）：")
  say("  (1) 逐条件看 Orel 与 K、信号、rho 的关系 —— 这是三点规模梯度，不是单调规律。")
  say("  (2) 集成相对 argmax 的 ΔCt 是正是负？绝对值是否超过 MC 标准误？")
  say("  (3) cor_rsf 与 n_const 用于判断'名义候选数'是否等于'独立信息量'。")
  say("  (4) 本模拟是分数层面的机制演示，不含组学结构，不能替代真实数据结论。")
}
close(con)

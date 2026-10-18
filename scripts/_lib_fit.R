# ============================================================
# _lib_fit.R —— 候选池与候选拟合（供 31_metabric_audit.R 使用）
#
# ⚠️ 本文件中的 POOL 与 fit_pool【必须与 scripts/19_main_corrected.R 中的定义逐字一致】。
#    为什么复制而不直接引用 19：19 是已冻结的主实验脚本，其产物已落盘；
#    把它改成 source() 会让"随包提供的 19"与"实际产出结果的 19"不一致。
#    因此这里保留一份**镜像**，并在 31 的日志里打印两者的 md5 供核对。
#
#    19 中的相关行：POOL 定义在 "候选池" 注释之后；fit_pool 紧接其后。
#    与 19 的唯一差别：本文件不在此处定义 `%||%`（由调用方提供或 base 版本）。
# ============================================================

POOL <- list(
  list(n="Lasso",        k="lin", a=1.0, lam="min"),
  list(n="Ridge",        k="lin", a=0.0, lam="min"),
  list(n="Enet0.9",      k="lin", a=0.9, lam="min"),
  list(n="Enet0.7",      k="lin", a=0.7, lam="min"),
  list(n="Enet0.5",      k="lin", a=0.5, lam="min"),
  list(n="Enet0.3",      k="lin", a=0.3, lam="min"),
  list(n="Enet0.1",      k="lin", a=0.1, lam="min"),
  list(n="Ridge_1se",    k="lin", a=0.0, lam="1se"),
  list(n="UniCox+Lasso", k="uni", lam="min"),
  list(n="LassoTop100",  k="lin", a=1.0, lam="min", pre=100),
  list(n="LassoTop50",   k="lin", a=1.0, lam="min", pre=50),
  list(n="LassoTop25",   k="lin", a=1.0, lam="min", pre=25),
  list(n="LassoTop10",   k="lin", a=1.0, lam="min", pre=10),
  list(n="Enet0.5Top25", k="lin", a=0.5, lam="min", pre=25),
  list(n="RSF100",       k="rsf", ntree=100),
  list(n="GBM100",       k="gbm", nt=100),
  list(n="GBM250",       k="gbm", nt=250),
  list(n="RSF250",       k="rsf", ntree=250)
)
POOL_NAMES <- vapply(POOL, function(m) m$n, character(1))

fit_pool <- function(Xa, ta, ea, Xb, Xc) {
  y <- Surv(ta, ea); out <- list()
  for (m in POOL) {
    nm <- m$n; s <- NULL
    Aa <- Xa; Ab <- Xb; Ac <- Xc
    if (!is.null(m$pre)) {
      zz <- apply(Aa, 2, function(x){
        q <- tryCatch(summary(coxph(Surv(ta,ea) ~ x))$coef, error=function(e) NULL)
        if (is.null(q)) NA_real_ else abs(q[1, "z"]) })
      zz[!is.finite(zz)] <- -Inf
      kk <- min(m$pre, ncol(Aa))
      if (kk < 5 || sum(zz > -Inf) < 5) next
      ix <- order(zz, decreasing = TRUE)[seq_len(kk)]
      Aa <- Aa[, ix, drop=FALSE]; Ab <- Ab[, ix, drop=FALSE]; Ac <- Ac[, ix, drop=FALSE]
    }
    if (m$k == "lin") {
      set.seed(1)
      cv <- tryCatch(cv.glmnet(Aa, y, family="cox", alpha=m$a, nfolds=5, cox.ties="breslow"),
                     error=function(e) NULL)
      if (!is.null(cv)) { l <- if (m$lam == "min") cv$lambda.min else cv$lambda.1se
        s <- list(pa=as.numeric(predict(cv, Aa, s=l, type="link")),
                  pb=as.numeric(predict(cv, Ab, s=l, type="link")),
                  pc=as.numeric(predict(cv, Ac, s=l, type="link"))) }
    } else if (m$k == "uni") {
      up <- apply(Aa, 2, function(x){
        z <- tryCatch(summary(coxph(Surv(ta,ea) ~ x))$coef, error=function(e) NULL)
        if (is.null(z)) NA_real_ else z[1, "Pr(>|z|)"] })
      k <- which(!is.na(up) & up < 0.05)
      if (length(k) >= 5) { set.seed(1)
        cv6 <- tryCatch(cv.glmnet(Aa[,k,drop=FALSE], y, family="cox", alpha=1, nfolds=5, cox.ties="breslow"),
                        error=function(e) NULL)
        if (!is.null(cv6)) s <- list(
          pa=as.numeric(predict(cv6, Aa[,k,drop=FALSE], s=cv6$lambda.min, type="link")),
          pb=as.numeric(predict(cv6, Ab[,k,drop=FALSE], s=cv6$lambda.min, type="link")),
          pc=as.numeric(predict(cv6, Ac[,k,drop=FALSE], s=cv6$lambda.min, type="link"))) }
    } else if (m$k == "rsf") {
      d <- data.frame(.time = ta, .status = ea, Aa, check.names = FALSE); set.seed(1)
      rf <- tryCatch(rfsrc(Surv(.time, .status) ~ ., data=d, ntree=m$ntree,
                           nodesize=m$ns %||% 5, splitrule="logrank"), error=function(e) NULL)
      if (!is.null(rf)) s <- list(
        pa=as.numeric(predict(rf, data.frame(Aa, check.names=FALSE))$predicted),
        pb=as.numeric(predict(rf, data.frame(Ab, check.names=FALSE))$predicted),
        pc=as.numeric(predict(rf, data.frame(Ac, check.names=FALSE))$predicted))
    } else if (m$k == "gbm") {
      set.seed(1)
      gb <- tryCatch(gbm.fit(x=Aa, y=Surv(ta,ea), distribution="coxph", n.trees=m$nt,
                             shrinkage=m$sh %||% 0.01, interaction.depth=m$d %||% 3,
                             n.minobsinnode=10, verbose=FALSE), error=function(e) NULL)
      if (!is.null(gb)) s <- list(
        pa=as.numeric(predict(gb, Aa, n.trees=m$nt, type="link")),
        pb=as.numeric(predict(gb, Ab, n.trees=m$nt, type="link")),
        pc=as.numeric(predict(gb, Ac, n.trees=m$nt, type="link")))
    }
    if (!is.null(s)) out[[nm]] <- s
  }
  out
}

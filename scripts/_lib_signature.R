# ============================================================
# _lib_signature.R —— 发表签名派生任务的候选池与拟合函数
#
# 出处：原样取自 scripts/26b_public_signature_tasks.R（2026-10-20 重构）。
# 目的：让 26b（全队列签名任务）与 53（TCGA-BRCA 同签名对照）共用同一份
#       候选池定义与 fit_pool 实现，避免两处拟合代码漂移。
# 只搬运，未改动任何逻辑或数值。
# ============================================================
suppressMessages({library(survival); library(glmnet); library(randomForestSRC); library(gbm)})
options(stringsAsFactors = FALSE)
`%||%` <- function(a, b) if (is.null(a)) b else a

TREES <- c("RSF100","GBM100","GBM250","RSF250")
MAX_MISS <- 0.20
ci_direct <- function(t, e, s) {
  s <- tryCatch(as.numeric(s), error=function(x) NULL)
  if (is.null(s) || length(s) != length(t)) return(NA_real_)
  ok <- is.finite(s) & !is.na(t) & !is.na(e)
  if (sum(ok) < 10) return(NA_real_)
  tryCatch(as.numeric(concordance(Surv(t[ok],e[ok]) ~ s[ok], reverse=TRUE, timewt="n")$concordance),
           error=function(x) NA_real_)
}
pick_max <- function(v, prio) { m <- suppressWarnings(max(v, na.rm=TRUE))
  if (!is.finite(m)) return(NA_integer_); cc <- which(is.finite(v) & v == m); cc[which.max(prio[cc])] }

POOL <- list(
  list(n="Lasso", k="lin", a=1.0), list(n="Ridge", k="lin", a=0.0),
  list(n="Enet0.9", k="lin", a=0.9), list(n="Enet0.7", k="lin", a=0.7),
  list(n="Enet0.5", k="lin", a=0.5), list(n="Enet0.3", k="lin", a=0.3),
  list(n="Enet0.1", k="lin", a=0.1), list(n="Ridge_1se", k="lin", a=0.0, l1se=TRUE),
  list(n="UniCox+Lasso", k="uni"),
  list(n="LassoTop100", k="lin", a=1.0, pre=100), list(n="LassoTop50", k="lin", a=1.0, pre=50),
  list(n="LassoTop25", k="lin", a=1.0, pre=25),   list(n="LassoTop10", k="lin", a=1.0, pre=10),
  list(n="Enet0.5Top25", k="lin", a=0.5, pre=25),
  list(n="RSF100", k="rsf", ntree=100), list(n="GBM100", k="gbm", nt=100),
  list(n="GBM250", k="gbm", nt=250),    list(n="RSF250", k="rsf", ntree=250))
POOL_NAMES <- vapply(POOL, function(m) m$n, character(1))

fit_pool <- function(Xa, ta, ea, Xb, Xc) {
  y <- Surv(ta, ea); out <- list()
  for (m in POOL) {
    nm <- m$n; s <- NULL; Aa <- Xa; Ab <- Xb; Ac <- Xc
    if (!is.null(m$pre)) {
      zz <- apply(Aa, 2, function(x){ q <- tryCatch(summary(coxph(Surv(ta,ea) ~ x))$coef, error=function(e) NULL)
        if (is.null(q)) NA_real_ else abs(q[1,"z"]) })
      zz[!is.finite(zz)] <- -Inf; kk <- min(m$pre, ncol(Aa))
      if (kk < 5 || sum(zz > -Inf) < 5) next
      ix <- order(zz, decreasing=TRUE)[seq_len(kk)]
      Aa <- Aa[,ix,drop=FALSE]; Ab <- Ab[,ix,drop=FALSE]; Ac <- Ac[,ix,drop=FALSE]
    }
    if (m$k == "lin") {
      set.seed(1)
      cv <- tryCatch(cv.glmnet(Aa, y, family="cox", alpha=m$a, nfolds=5, cox.ties="breslow"),
                     error=function(e) NULL)
      if (!is.null(cv)) { l <- if (isTRUE(m$l1se)) cv$lambda.1se else cv$lambda.min
        s <- list(pa=as.numeric(predict(cv,Aa,s=l,type="link")),
                  pb=as.numeric(predict(cv,Ab,s=l,type="link")),
                  pc=as.numeric(predict(cv,Ac,s=l,type="link"))) }
    } else if (m$k == "uni") {
      up <- apply(Aa, 2, function(x){ z <- tryCatch(summary(coxph(Surv(ta,ea)~x))$coef, error=function(e) NULL)
        if (is.null(z)) NA_real_ else z[1,"Pr(>|z|)"] })
      k <- which(!is.na(up) & up < 0.05)
      if (length(k) >= 5) { set.seed(1)
        cv6 <- tryCatch(cv.glmnet(Aa[,k,drop=FALSE], y, family="cox", alpha=1, nfolds=5, cox.ties="breslow"),
                        error=function(e) NULL)
        if (!is.null(cv6)) s <- list(pa=as.numeric(predict(cv6,Aa[,k,drop=FALSE],s=cv6$lambda.min,type="link")),
                                     pb=as.numeric(predict(cv6,Ab[,k,drop=FALSE],s=cv6$lambda.min,type="link")),
                                     pc=as.numeric(predict(cv6,Ac[,k,drop=FALSE],s=cv6$lambda.min,type="link"))) }
    } else if (m$k == "rsf") {
      d <- data.frame(.time=ta, .status=ea, Aa, check.names=FALSE); set.seed(1)
      rf <- tryCatch(rfsrc(Surv(.time,.status) ~ ., data=d, ntree=m$ntree, splitrule="logrank"),
                     error=function(e) NULL)
      if (!is.null(rf)) s <- list(pa=as.numeric(predict(rf, data.frame(Aa,check.names=FALSE))$predicted),
                                  pb=as.numeric(predict(rf, data.frame(Ab,check.names=FALSE))$predicted),
                                  pc=as.numeric(predict(rf, data.frame(Ac,check.names=FALSE))$predicted))
    } else if (m$k == "gbm") {
      set.seed(1)
      gb <- tryCatch(gbm.fit(x=Aa, y=Surv(ta,ea), distribution="coxph", n.trees=m$nt,
                             shrinkage=0.01, interaction.depth=3, n.minobsinnode=10, verbose=FALSE),
                     error=function(e) NULL)
      if (!is.null(gb)) s <- list(pa=as.numeric(predict(gb,Aa,n.trees=m$nt,type="link")),
                                  pb=as.numeric(predict(gb,Ab,n.trees=m$nt,type="link")),
                                  pc=as.numeric(predict(gb,Ac,n.trees=m$nt,type="link")))
    }
    if (!is.null(s)) out[[nm]] <- s
  }
  out
}

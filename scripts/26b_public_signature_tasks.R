#!/usr/bin/env Rscript
# ============================================================
# 26b — 现实任务模块：以【已发表乳腺癌预后签名】的基因集为特征集
#
# 命名与口径（评审明确要求）：
#   本模块做的是"**发表签名衍生任务**"——取已发表签名的【基因列表】，
#   在这个基因集上**重新训练**候选池，而不是复现原模型。
#   要复现原模型必须遵守其系数、标准化与使用场景；那不在本模块范围内。
#
# 做法：对每个（签名 × 队列 × 划分种子）
#   · 取该签名的基因与队列表达谱的交集；要求 >= 8 个基因
#   · 使用与主分析【完全相同】的划分种子与训练集标准化规则
#   · 在**该签名的基因集**上拟合候选池，按同一套规则核算 Gv / Gt / Orel
# 从而回答：把特征集从"MSigDB 基因集合"换成"真实发表签名的基因"，
# 收益核算的结论是否改变。
#
# 输入：results/manifest/public_signatures_genes.tsv（由 26a 生成）
# 输出：results/tables/26_signature_tasks.tsv, 26_signature_summary.tsv
# ============================================================
suppressMessages({library(survival); library(glmnet); library(randomForestSRC); library(gbm)})
options(stringsAsFactors = FALSE)
`%||%` <- function(a, b) if (is.null(a)) b else a
PROJ <- Sys.getenv("SIGAUDIT_ROOT", "D:/bioinfo/sigAudit"); TAB <- file.path(PROJ,"results/tables")
LOGDIR <- file.path(PROJ,"logs"); con <- file(file.path(LOGDIR,"26_signature_tasks.log"), "wt")
say <- function(...) { m <- paste0(...); cat(m,"\n",sep=""); cat(m,"\n",sep="",file=con); flush(con) }

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

ENDPOINTS <- list(
  LGG  = list(e=file.path(PROJ, "data/cohorts/LGG_expr.tsv"),
              f=file.path(PROJ, "data/cohorts/LGG_surv.tsv"),
              t="PFS_MONTHS", s="PFS_STATUS", endpoint="PFS",
              event=c("1:PROGRESSION"), censor=c("0:CENSORED")),
  PAAD = list(t="PFS_MONTHS", s="PFS_STATUS", endpoint="PFS",
              event=c("1:PROGRESSION"), censor=c("0:CENSORED")),
  BRCA = list(t="OS_MONTHS", s="OS_STATUS", endpoint="OS", event=c("1:DECEASED"), censor=c("0:LIVING")),
  COAD = list(t="OS_MONTHS", s="OS_STATUS", endpoint="OS", event=c("1:DECEASED"), censor=c("0:LIVING")),
  KIRC = list(t="OS_MONTHS", s="OS_STATUS", endpoint="OS", event=c("1:DECEASED"), censor=c("0:LIVING")),
  LUAD = list(t="OS_MONTHS", s="OS_STATUS", endpoint="OS", event=c("1:DECEASED"), censor=c("0:LIVING")),
  MESO = list(t="OS_MONTHS", s="OS_STATUS", endpoint="OS", event=c("1:DECEASED"), censor=c("0:LIVING")),
  STAD = list(t="OS_MONTHS", s="OS_STATUS", endpoint="OS", event=c("1:DECEASED"), censor=c("0:LIVING"))
)
SEEDS <- c(20261014, 20261016, 20261017, 20261018, 20261019)

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

## ---------- 签名清单 ----------
## 【第五轮修正】改用 HGNC 规范化后的基因符号（scripts/34_map_signature_symbols.py）。
## 旧符号（如 KNTC2）映射为批准符号（NDC80）属于**标识规范化**，不是按结局挑基因。
LF <- file.path(PROJ,"results/manifest/public_signatures_genes_mapped.tsv")
if (!file.exists(LF)) {
  LF0 <- file.path(PROJ,"results/manifest/public_signatures_genes.tsv")
  if (!file.exists(LF0)) stop("缺少签名清单，请先跑 scripts/26a_extract_public_signatures.R")
  warning("未找到 mapped 版本，回退到未映射清单；请先运行 scripts/34_map_signature_symbols.py")
  L <- read.table(LF0, sep="\t", header=TRUE, quote="", check.names=FALSE)
  L$gene_source <- L$gene; L$gene_used <- L$gene
} else L <- read.table(LF, sep="\t", header=TRUE, quote="", check.names=FALSE)
SIGS  <- lapply(split(L, L$signature), function(d) unique(d$gene_used))
SRC_N <- vapply(split(L$gene_source, L$signature), function(x) length(unique(x)), integer(1))
say("纳入的已发表签名: ", paste(sprintf("%s(%d)", names(SIGS), sapply(SIGS, length)), collapse=", "))

## ---------- 主循环 ----------
res <- list()
for (cn in names(ENDPOINTS)) {
  E <- ENDPOINTS[[cn]]
  cf <- E[["e"]] %||% sprintf("%s/data/cohorts/%s_expr.tsv", PROJ, cn)
  sf <- E[["f"]] %||% sprintf("%s/data/cohorts/%s_surv.tsv", PROJ, cn)
  expr <- as.matrix(read.table(cf, header=TRUE, sep="\t", row.names=1, check.names=FALSE))
  if (nrow(expr) < ncol(expr)) expr <- t(expr)
  expr <- log2(expr + 1)
  sv <- read.table(sf, header=TRUE, sep="\t", row.names=1, check.names=FALSE)
  raw <- as.character(sv[[E[["s"]]]])
  sv$t <- suppressWarnings(as.numeric(sv[[E[["t"]]]]))
  sv$e <- ifelse(raw %in% E[["event"]], 1L, ifelse(raw %in% E[["censor"]], 0L, NA_integer_))
  sv <- sv[!is.na(sv$t) & !is.na(sv$e) & sv$t > 0, ]
  colnames(expr) <- substr(colnames(expr),1,12); expr <- expr[,!duplicated(colnames(expr)),drop=FALSE]
  com <- intersect(colnames(expr), rownames(sv)); expr <- expr[,com,drop=FALSE]; sv <- sv[com,]
  say(""); say("===== ", cn, " (", E[["endpoint"]], ", ", ncol(expr), " 样本) =====")
  for (sd in SEEDS) {
    set.seed(sd)
    ev <- which(sv$e==1); ce <- which(sv$e==0)
    tr <- sort(c(sample(ev, floor(.5*length(ev))), sample(ce, floor(.5*length(ce)))))
    rest <- setdiff(seq_along(sv$e), tr)
    re <- rest[sv$e[rest]==1]; rc <- rest[sv$e[rest]==0]
    va <- sort(c(sample(re, floor(.5*length(re))), sample(rc, floor(.5*length(rc)))))
    te <- setdiff(rest, va)
    for (sg in names(SIGS)) {
      ## 基因数的四个阶段必须分别记录（第三轮审查：旧版只记初始交集，
      ## 掩盖了之后缺失过滤与零方差剔除的损失）
      n_gene_source  <- unname(SRC_N[[sg]])
      g <- intersect(SIGS[[sg]], rownames(expr))
      n_gene_initial <- length(g)
      if (length(g) < 8) { say(sprintf("  跳过 %s/%s: 基因交集仅 %d", cn, sg, length(g))); next }
      sub <- expr[g, , drop=FALSE]
      trm <- sub[, tr, drop=FALSE]; miss <- rowSums(is.na(trm))/ncol(trm)
      ex <- sub[miss <= MAX_MISS, , drop=FALSE]
      n_gene_after_miss <- nrow(ex)
      if (nrow(ex) < 8) { say(sprintf("  跳过 %s/%s: 通过缺失过滤后 <8 基因", cn, sg)); next }
      med <- apply(ex[, tr, drop=FALSE], 1, function(v) median(v, na.rm=TRUE))
      for (j in seq_len(nrow(ex))) { v <- ex[j,]; k <- is.na(v); if (any(k)) { v[k] <- med[j]; ex[j,] <- v } }
      mu <- rowMeans(ex[, tr, drop=FALSE]); sdd <- apply(ex[, tr, drop=FALSE], 1, sd)
      kp <- is.finite(sdd) & sdd > 0; ex <- ex[kp,,drop=FALSE]; mu <- mu[kp]; sdd <- sdd[kp]
      n_gene_final <- nrow(ex)
      genes_final  <- paste(rownames(ex), collapse = "|")
      Z <- (ex - mu)/sdd; Z[!is.finite(Z)] <- 0
      G <- function(idx) t(Z[, colnames(Z)[idx], drop=FALSE])
      pm <- tryCatch(fit_pool(G(tr), sv$t[tr], sv$e[tr], G(va), G(te)), error=function(e) NULL)
      if (is.null(pm) || length(pm) < 4) { say(sprintf("  跳过 %s/%s: 可用候选 <4", cn, sg)); next }
      mods <- POOL_NAMES[POOL_NAMES %in% names(pm)]
      if (!("Lasso" %in% mods)) next
      mk <- function(key) { M <- sapply(mods, function(x) pm[[x]][[key]]); if (is.null(dim(M))) M <- matrix(M,ncol=1)
        colnames(M) <- mods; M }
      Pv <- mk("pb"); Pt <- mk("pc")
      V <- vapply(mods, function(m) ci_direct(sv$t[va], sv$e[va], Pv[,m]), numeric(1))
      Tt<- vapply(mods, function(m) ci_direct(sv$t[te], sv$e[te], Pt[,m]), numeric(1))
      if (!all(is.finite(V)) || !all(is.finite(Tt))) next
      set.seed(sd); prio <- runif(length(mods)); names(prio) <- mods
      j <- pick_max(V, prio); w <- mods[j]
      gv <- unname(V[w] - V["Lasso"]); gt <- unname(Tt[w] - Tt["Lasso"])
      top5 <- mods[order(-V, -prio[mods])][seq_len(min(5,length(mods)))]
      mu5 <- colMeans(Pv[, top5, drop=FALSE]); sd5 <- apply(Pv[, top5, drop=FALSE],2,sd); sd5[!is.finite(sd5)|sd5==0] <- 1
      ev5 <- rowMeans(sweep(sweep(Pv[,top5,drop=FALSE],2,mu5,"-"),2,sd5,"/"))
      et5 <- rowMeans(sweep(sweep(Pt[,top5,drop=FALSE],2,mu5,"-"),2,sd5,"/"))
      ct5 <- ci_direct(sv$t[te], sv$e[te], et5)
      res[[length(res)+1]] <- data.frame(cohort=cn, endpoint=E[["endpoint"]], seed=sd,
        signature=sg, n_gene=length(g), n_candidate=length(mods),
        n_gene_source=n_gene_source, n_gene_initial=n_gene_initial,
        n_gene_after_miss=n_gene_after_miss, n_gene_final=n_gene_final,
        genes_final=genes_final,
        Cv=V[w], Ct=Tt[w], Gv=gv, Gt=gt, Orel=gv-gt,
        Ct_ens=ct5, dCt_ens=unname(ct5 - Tt[w]),
        argmax_picks_ref=as.integer(w == "Lasso"))
    }
  }
}
if (!length(res)) stop("没有任何可用结果")
R <- do.call(rbind, res)
write.table(R, file.path(TAB,"26_signature_tasks.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

## ---------- 汇总（先队列内跨划分平均 -> 再队列等权）----------
bc <- aggregate(cbind(Cv,Ct,Gv,Gt,Orel,dCt_ens,argmax_picks_ref) ~ signature + cohort,
                data=R, FUN=function(x) mean(x, na.rm=TRUE))
ag <- aggregate(cbind(Ct,Gv,Gt,Orel,dCt_ens,argmax_picks_ref) ~ signature, data=bc,
                FUN=function(x) mean(x, na.rm=TRUE))
ag$n_cohort <- as.numeric(table(bc$signature)[ag$signature])
write.table(ag, file.path(TAB,"26_signature_summary.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
say(""); say("================ 汇总（发表签名衍生任务；先队列内跨划分平均，再队列等权）================")
say(sprintf("  %-22s %4s %9s %9s %9s %9s %10s", "签名","队列","Ct","Gv","Gt","Orel","集成ΔCt"))
for (i in seq_len(nrow(ag))) say(sprintf("  %-22s %4.0f %9.4f %+9.4f %+9.4f %+9.4f %+10.4f",
  ag$signature[i], ag$n_cohort[i], ag$Ct[i], ag$Gv[i], ag$Gt[i], ag$Orel[i], ag$dCt_ens[i]))
say("")
say("  口径说明：取已发表签名的【基因列表】，在该基因集上【重新训练】候选池，")
say("  属「发表签名衍生任务」，**不是**原模型复现（复现需遵守原系数与标准化）。")
say(sprintf("  共 %d 个（签名×队列×划分）单元，%d 个签名 × %d 个队列 × %d 个划分",
            nrow(R), length(unique(R$signature)), length(unique(R$cohort)), length(unique(R$seed))))
close(con)

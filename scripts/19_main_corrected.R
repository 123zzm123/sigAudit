#!/usr/bin/env Rscript
# ============================================================
# 19 — 修正版主实验数据生成（第二轮评审意见落地）
#
# 相对上一版的改动：
#
#  【B1 阻断｜部分产物被当成完成结果】
#     上一版每处理 10 个基因集就写【最终文件名】，运行链见到文件即跳过，
#     因此中途中断会留下一个"看起来合法"的部分结果（20 只要求 >=30 个基因集）。
#     现在：进度写 <tag>.partial.rds；只有遍历完全部计划基因集后才写最终文件，
#     且文件中记录 completed / n_planned / n_processed / 时间戳 / 输入哈希。
#
#  【B5 阻断｜缺失值被静默置零】
#     上一版 `Z[!is.finite(Z)] <- 0` 把缺失变成"该基因的平均表达"。
#     现在：缺失比例只在【训练集】上统计；超过阈值的基因整列剔除；
#     剩余缺失用【训练集中位数】填补；全过程计数并落盘。
#     （实测：COADREAD 4.91%、STAD 1.00% 缺失，其余 6 队列为 0）
#
#  【B5 阻断｜终点事件语义】
#     不再用宽松正则，改为数据字典里的【显式取值映射】；
#     未知取值直接报错，不让"0:"与文字标签冲突时由正则裁决。
#
#  【元数据】
#     保存样本 ID、train/val/test 归属、基因名、队列/种子/终点/输入哈希，
#     以便事后核对预测与结局的对齐（评审要求）。
#
# 输出：results/main/19_<COHORT>_<SEED>_<N>.rds      仅当全部计划基因集处理完
#       results/main/19_<COHORT>_<SEED>_<N>.partial.rds  中途进度（可删）
# ============================================================
suppressMessages({library(survival); library(glmnet); library(randomForestSRC); library(gbm); library(jsonlite)})
options(stringsAsFactors = FALSE)
`%||%` <- function(a, b) if (is.null(a)) b else a

PROJ    <- Sys.getenv("SIGAUDIT_ROOT", "D:/bioinfo/sigAudit")
LOGDIR  <- file.path(PROJ, "logs");            dir.create(LOGDIR, showWarnings = FALSE, recursive = TRUE)
OUTDIR  <- file.path(PROJ, "results/main");    dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
MAX_MISS <- 0.20        # 训练集缺失比例超过此值的基因整列剔除

ARGS   <- commandArgs(trailingOnly = TRUE)
N_GS   <- if (length(ARGS) >= 1) as.integer(ARGS[1]) else 100
COHORT <- if (length(ARGS) >= 2) ARGS[2] else "MESO"
SEED   <- if (length(ARGS) >= 3) as.integer(ARGS[3]) else 20261014L
TAG    <- sprintf("19_%s_%d_%d", COHORT, SEED, N_GS)
LOG    <- file.path(LOGDIR, paste0(TAG, ".log")); con <- file(LOG, "wt")
say    <- function(...) { m <- paste0(...); cat(m,"\n",sep=""); cat(m,"\n",sep="",file=con); flush(con) }

# ---- 终点字典：显式取值映射（不再用正则裁决事件语义）----
ENDPOINTS <- list(
  LGG  = list(e = file.path(PROJ, "data/cohorts/LGG_expr.tsv"),
              f = file.path(PROJ, "data/cohorts/LGG_surv.tsv"),
              t = "PFS_MONTHS", s = "PFS_STATUS", endpoint = "PFS",
              event = c("1:PROGRESSION"), censor = c("0:CENSORED")),
  PAAD = list(t = "PFS_MONTHS", s = "PFS_STATUS", endpoint = "PFS",
              event = c("1:PROGRESSION"), censor = c("0:CENSORED")),
  BRCA = list(t = "OS_MONTHS", s = "OS_STATUS", endpoint = "OS",
              event = c("1:DECEASED"), censor = c("0:LIVING")),
  COAD = list(t = "OS_MONTHS", s = "OS_STATUS", endpoint = "OS",
              event = c("1:DECEASED"), censor = c("0:LIVING")),
  KIRC = list(t = "OS_MONTHS", s = "OS_STATUS", endpoint = "OS",
              event = c("1:DECEASED"), censor = c("0:LIVING")),
  LUAD = list(t = "OS_MONTHS", s = "OS_STATUS", endpoint = "OS",
              event = c("1:DECEASED"), censor = c("0:LIVING")),
  MESO = list(t = "OS_MONTHS", s = "OS_STATUS", endpoint = "OS",
              event = c("1:DECEASED"), censor = c("0:LIVING")),
  STAD = list(t = "OS_MONTHS", s = "OS_STATUS", endpoint = "OS",
              event = c("1:DECEASED"), censor = c("0:LIVING"))
)
if (!COHORT %in% names(ENDPOINTS)) stop("终点字典里没有该队列: ", COHORT)
E  <- ENDPOINTS[[COHORT]]
STUDY_ID <- if (COHORT == "LGG") "local:LGG_broad" else sprintf("%s_tcga_pan_can_atlas_2018",
                                                                 tolower(if (COHORT=="COAD") "coadread" else COHORT))
cf <- E[["e"]] %||% sprintf("%s/data/cohorts/%s_expr.tsv", PROJ, COHORT)
sf <- E[["f"]] %||% sprintf("%s/data/cohorts/%s_surv.tsv", PROJ, COHORT)

say("载入 ", COHORT, "  study=", STUDY_ID, "  终点=", E[["endpoint"]], "  划分种子=", SEED,
    "  计划基因集=", N_GS)
expr <- as.matrix(read.table(cf, header=TRUE, sep="\t", row.names=1, check.names=FALSE))
if (nrow(expr) < ncol(expr)) expr <- t(expr)          # -> 基因 × 样本
expr <- log2(expr + 1)                                # 实测原始 RSEM 量级，单次 log，非二次
sv <- read.table(sf, header=TRUE, sep="\t", row.names=1, check.names=FALSE)
for (nm in c(E[["t"]], E[["s"]])) if (!nm %in% names(sv)) stop("生存表缺少列: ", nm)
raw_s <- as.character(sv[[E[["s"]]]])
# 空串/纯空白 视为缺失（实测 LUAD 的 OS_STATUS 含空字符串），
# 不能当成"字典未定义的类别"而报错；同时也避免被当成 0 事件。
blank <- is.na(raw_s) | !nzchar(trimws(raw_s))
miss_status <- sum(blank)
raw_s[blank] <- NA_character_
bad <- setdiff(unique(raw_s[!is.na(raw_s)]), c(E[["event"]], E[["censor"]]))
if (length(bad)) stop("状态列出现字典未定义的取值: ", paste(bad, collapse=" / "))
say("  状态列空值/缺失: ", miss_status, " 例（按缺失剔除）")
sv$t <- suppressWarnings(as.numeric(sv[[E[["t"]]]]))
sv$e <- ifelse(raw_s %in% E[["event"]], 1L, ifelse(raw_s %in% E[["censor"]], 0L, NA_integer_))
sv <- sv[!is.na(sv$t) & !is.na(sv$e) & sv$t > 0, ]
colnames(expr) <- substr(colnames(expr),1,12)
expr <- expr[, !duplicated(colnames(expr)), drop=FALSE]
com <- intersect(colnames(expr), rownames(sv)); expr <- expr[,com,drop=FALSE]; sv <- sv[com,]
say("  ", ncol(expr), " 样本 × ", nrow(expr), " 基因, 事件 ", sum(sv$e),
    " (", E[["endpoint"]], ", 事件取值=", paste(E[["event"]], collapse=","), ")")

# ---------------- 三分法 ----------------
set.seed(SEED)
ev <- which(sv$e==1); ce <- which(sv$e==0)
tr <- sort(c(sample(ev, floor(.5*length(ev))), sample(ce, floor(.5*length(ce)))))
rest <- setdiff(seq_along(sv$e), tr)
re <- rest[sv$e[rest]==1]; rc <- rest[sv$e[rest]==0]
va <- sort(c(sample(re, floor(.5*length(re))), sample(rc, floor(.5*length(rc)))))
te <- setdiff(rest, va)
say("  训练 ", length(tr), " / 验证 ", length(va), " / 测试 ", length(te),
    "  事件 ", sum(sv$e[tr]), "/", sum(sv$e[va]), "/", sum(sv$e[te]))

# ---------------- B5：缺失处理规则（只在训练集上确定）----------------
say("  缺失处理：阈值 ", MAX_MISS, "（仅用训练集统计）")
trm     <- expr[, tr, drop=FALSE]
n_na_tr <- rowSums(is.na(trm))
miss_frac <- n_na_tr / ncol(trm)
drop_g  <- miss_frac > MAX_MISS
ex <- expr[!drop_g, , drop=FALSE]
gmed <- apply(ex[, tr, drop=FALSE], 1, median, na.rm = TRUE)
n_imp <- sum(is.na(ex))
if (n_imp > 0) {
  for (j in seq_len(nrow(ex))) { v <- ex[j, ]; k <- is.na(v)
    if (any(k)) v[k] <- gmed[j]; ex[j, ] <- v }
}
gmu <- rowMeans(ex[, tr, drop=FALSE])
gsd <- apply(ex[, tr, drop=FALSE], 1, sd)
keep_g <- is.finite(gsd) & gsd > 0
ex <- ex[keep_g, , drop=FALSE]; gmu <- gmu[keep_g]; gsd <- gsd[keep_g]; gsd[gsd==0] <- 1
Z <- (ex - gmu) / gsd
Z[!is.finite(Z)] <- 0        # 此时只可能是 0/0 之类的残渣；已在下方计数
say("   训练集缺失>", MAX_MISS, " 而剔除的基因: ", sum(drop_g),
    " ；填补的缺失值个数: ", n_imp,
    " ；零方差剔除: ", sum(!keep_g), " ；保留 ", nrow(Z), " 个基因")
say("   标准化参数（均值/标准差/中位数）全部由训练集估计，原样应用于验证与测试")
G <- function(idx) t(Z[, colnames(Z)[idx], drop=FALSE])

# ---------------- 候选池 ----------------
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
say("候选池 ", length(POOL), " 个（算法 × 超参 × 筛选策略的组合，非独立算法家族）")

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
ci_direct <- function(t, e, s) {
  s <- tryCatch(as.numeric(s), error=function(x) NULL)
  if (is.null(s) || length(s) != length(t)) return(NA_real_)
  ok <- is.finite(s) & !is.na(t) & !is.na(e)
  if (sum(ok) < 10) return(NA_real_)
  tryCatch(as.numeric(concordance(Surv(t[ok],e[ok]) ~ s[ok], reverse=TRUE, timewt="n")$concordance),
           error=function(x) NA_real_)
}

# ---------------- 基因集 ----------------
gs <- fromJSON(file.path(PROJ,"data/signatures/msigdb_genesets.json"))
gs <- gs[sapply(gs,length)>=10 & sapply(gs,length)<=500]
set.seed(20261015); pick <- sample(names(gs), min(N_GS, length(gs)))
say("基因集 ", length(pick), " 个（抽样种子固定 20261015，跨划分一致）")

OUT     <- file.path(OUTDIR, paste0(TAG, ".rds"))
OUTPART <- file.path(OUTDIR, paste0(TAG, ".partial.rds"))
in_hash <- paste0(substr(as.character(tools::md5sum(cf)), 1, 8), "-",
                   substr(as.character(tools::md5sum(sf)), 1, 8))
if (file.exists(OUTPART)) file.remove(OUTPART)

out <- list(); t0 <- Sys.time(); n_skip <- 0; diag_tr <- list()
Xtr <- G(tr); Xva <- G(va); Xte <- G(te)
ttr<-sv$t[tr]; etr<-sv$e[tr]; tva<-sv$t[va]; eva<-sv$e[va]; tte<-sv$t[te]; ete<-sv$e[te]
META <- list(cohort=COHORT, study_id=STUDY_ID, endpoint=E[["endpoint"]],
             seed=SEED, n_planned=length(pick), completed=FALSE,
             sample_ids=colnames(Z), train_ids=colnames(Z)[tr], val_ids=colnames(Z)[va],
             test_ids=colnames(Z)[te], genes=rownames(Z),
             n_gene_dropped_missing=sum(drop_g), n_value_imputed=n_imp,
             n_gene_dropped_zerovar=sum(!keep_g), max_miss=MAX_MISS,
             expr_md5=in_hash, script="19_main_corrected.R")

for (i in seq_along(pick)) {
  nm <- pick[i]; g <- intersect(gs[[nm]], rownames(Z))
  if (length(g) < 10) { n_skip <- n_skip + 1; next }
  pm <- tryCatch(fit_pool(Xtr[,g,drop=FALSE], ttr, etr, Xva[,g,drop=FALSE], Xte[,g,drop=FALSE]),
                 error=function(e) NULL)
  if (is.null(pm) || length(pm) < 3) { n_skip <- n_skip + 1; next }
  mods <- POOL_NAMES[POOL_NAMES %in% names(pm)]
  if (!("Lasso" %in% mods)) { n_skip <- n_skip + 1; next }
  mk <- function(key) { M <- sapply(mods, function(x) pm[[x]][[key]])
    if (is.null(dim(M))) M <- matrix(M, ncol=1); colnames(M) <- mods; M }
  PA <- mk("pa"); PB <- mk("pb"); PC <- mk("pc")
  diag_tr[[nm]] <- data.frame(geneset=nm, model=mods,
    train_ci = vapply(seq_along(mods), function(j) ci_direct(ttr, etr, PA[,j]), numeric(1)))
  out[[nm]] <- list(name=nm, n_gene=length(g), models=mods,
                    train_pred=PA, val_pred=PB, test_pred=PC,
                    t_tr=ttr, e_tr=etr, tva=tva, eva=eva, tte=tte, ete=ete)
  if (i %% 10 == 0) {
    META$n_processed <- i
    saveRDS(c(list(genesets=out, train_ci=diag_tr), META), OUTPART)   # 进度只写 .partial
    say("  进度 ", i, "/", length(pick), "  已完成 ", length(out), "  跳过 ", n_skip,
        "  用时 ", round(as.numeric(difftime(Sys.time(), t0, units="min")),1), " 分钟")
    gc(verbose = FALSE)
  }
}

if (length(out)) {
  META$n_processed <- length(pick)
  META$completed   <- TRUE                       # 只有走完全部计划基因集才置 TRUE
  META$finished_at <- format(Sys.time())
  saveRDS(c(list(genesets=out, train_ci=diag_tr), META), OUT)
  if (file.exists(OUTPART)) file.remove(OUTPART)
  say(""); say("==== 19 完成 (", COHORT, " seed=", SEED, ") ====")
  say("基因集 ", length(out), "，跳过 ", n_skip, "，n_planned=", length(pick),
      "，completed=TRUE")
  say("验证 ", length(tva), " / 测试 ", length(tte), " 例；终点 ", E[["endpoint"]])
  D <- do.call(rbind, diag_tr)
  a <- aggregate(train_ci ~ model, data=D, FUN=function(x) mean(x, na.rm=TRUE))
  a <- a[order(a$train_ci), ]
  say("方向诊断：各候选【训练集】C-index 均值（应全部 > 0.5）")
  for (i in seq_len(nrow(a))) say(sprintf("    %-14s %.3f", a$model[i], a$train_ci[i]))
  bad <- a$model[a$train_ci < 0.5]
  if (length(bad)) say("  ★ 警告：以下候选训练集 C<0.5，方向约定可能不一致: ", paste(bad, collapse=", "))
  else say("  ✓ 全部候选在训练集上 C>0.5，方向约定一致")
} else { say("没有完成的基因集！"); if (file.exists(OUTPART)) file.remove(OUTPART) }
close(con)

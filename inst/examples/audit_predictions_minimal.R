## Minimal, self-contained example for the prediction-level audit.
## It runs after installation and needs no study data:
##   source(system.file("examples", "audit_predictions_minimal.R", package = "sigAudit"))
library(sigAudit)

set.seed(1)
n <- 120
mk <- function() {
  m <- cbind(Lasso = rnorm(n), Ridge = rnorm(n), RSF100 = rnorm(n), GBM100 = rnorm(n))
  rownames(m) <- paste0("p", seq_len(n))
  m
}
val  <- mk()
test <- mk()

tk <- list(list(
  val_pred  = val,
  test_pred = test,
  time_val  = rexp(n), event_val  = rbinom(n, 1, 0.6),
  time_test = rexp(n), event_test = rbinom(n, 1, 0.6),
  cohort = "C1", split = "s1"
))

fam <- c(Lasso = "linear", Ridge = "linear", RSF100 = "tree", GBM100 = "tree")

res <- audit_predictions(tk, reference = "Lasso", top_k = 4L, family = fam)
print(res)
print(res$summary)

# Examples

Two entry points are illustrated here. All of them use the **installed** package
(`library(sigAudit)`); none of them sources the study scripts under `scripts/`.

| File | Kind | What it does |
|---|---|---|
| `audit_predictions_minimal.R` | smoke test | builds three random candidates in memory and calls `audit_predictions()`. Fast, but synthetic. |
| `audit_predictions_metabric.R` | real data | audits one **saved METABRIC task** (18 candidate risk scores, 494 validation and 495 test patients) and then re-audits the same prediction against a different fixed reference. |
| `export_audit_example.R` | provenance | shows how the real-data files in `inst/extdata/` were derived from the frozen local predictions. Only needed to regenerate them. |

## Running the real-data example

```r
# install.packages("remotes")
remotes::install_local(".")          # or: R CMD INSTALL .
source(system.file("examples", "audit_predictions_metabric.R", package = "sigAudit"))
```

Measured on the reference machine (R 4.6.1, x86_64-w64-mingw32, Windows): one
`audit_predictions()` call for the requested rules (`none`, `argmax_val`, `top5_ensemble`) on
this **already-loaded** 18-candidate, 989-record task takes about **0.55–0.62 s on the first
call in a fresh session**, and about **0.1 s** for further calls in the same session. The timer
is started after the score matrices and outcomes have been read from disk and assembled into a
task, so the figure excludes file reading and excludes model fitting; it is also not a benchmark
of every rule (the resampling remedies refit nothing but call the concordance routine
repeatedly and are not requested here).

## What the example shows

1. A single call returns candidate diagnostics, the fixed-reference comparison and the rule
   results (`res$summary`, `res$diagnostics`, `res$rules`).
2. Switching the fixed reference from `Lasso` to `Ridge` leaves the selected candidate's `Cv`
   and `Ct` unchanged for `argmax_val` and `top5_ensemble` (the candidate pool is the same),
   while `none` *is* the reference and so its `Cv`/`Ct` change with it. For `argmax_val` and
   `top5_ensemble`, the gains and `Orel` change with the reference; for `none`, they remain zero
   by definition.

## Data

`inst/extdata/audit_example_val.tsv`, `_test.tsv` and `_survival.tsv` hold only derived model
outputs (risk scores) and survival outcomes for one signature/split; patient identifiers were
rewritten to `V####` / `T####`. The signature and split were chosen by a fixed, result-independent
rule (lexicographically smallest signature name and smallest split seed), not by performance.
See `export_audit_example.R` for the exact derivation.

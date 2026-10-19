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
`audit_predictions()` call on this task takes about **0.5 s** including reading the two score
matrices from disk. No model is fitted, and the wall time is not a training time.

## What the example shows

1. A single call returns candidate diagnostics, the fixed-reference comparison and the rule
   results (`res$summary`, `res$diagnostics`, `res$rules`).
2. Switching the fixed reference from `Lasso` to `Ridge` leaves each selected model's `Cv` and
   `Ct` unchanged, because those are properties of the selected model, while `Gv`, `Gt` and
   `Orel` move with the reference.

## Data

`inst/extdata/audit_example_val.tsv`, `_test.tsv` and `_survival.tsv` hold only derived model
outputs (risk scores) and survival outcomes for one signature/split; patient identifiers were
rewritten to `V####` / `T####`. The signature and split were chosen by a fixed, result-independent
rule (lexicographically smallest signature name and smallest split seed), not by performance.
See `export_audit_example.R` for the exact derivation.

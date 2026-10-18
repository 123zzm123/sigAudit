# sigAudit

Reference-adjusted gain accounting for multi-algorithm search in omics prognostic modelling.

Repository: <https://github.com/123zzm123/sigAudit>

Multi-algorithm pipelines fit many learners, keep the one with the best validation
performance, and report that number. `sigAudit` separates the two things that number mixes:
the gain that survives on held-out data, and the gain that depends on the comparator.

## What it computes

For a candidate pool evaluated on a validation and a held-out test partition:

| Quantity | Definition | Meaning |
|---|---|---|
| `Gv` | `Cv(rule) − Cv(reference)` | validation gain over a fixed reference |
| `Gt` | `Ct(rule) − Ct(reference)` | held-out gain over the same reference |
| `Orel` | `Gv − Gt` | reference-adjusted selection optimism |
| `Oraw` | `Cv(rule) − Ct(rule)` | rule's own validation–test difference; never conflated with `Orel` |

`audit_rules()` compares selection rules (argmax, top-k ensembles, resampling remedies) on a
common task set; `audit_selection()` performs the paired accounting for a single rule.

## Install

```r
# install.packages("remotes")
remotes::install_github("123zzm123/sigAudit")
```

## Quick start

```r
library(sigAudit)

d  <- read.csv(system.file("extdata", "example_ci.csv", package = "sigAudit"))
V  <- xtabs(val_ci  ~ geneset + model, data = d)   # rows = tasks, columns = candidates
Tt <- xtabs(test_ci ~ geneset + model, data = d)

res <- audit_rules(V, Tt, reference = "Lasso")
res$table   # rule, n, reported, test, Gv, Gt, Orel, Oraw
```

`audit_rules()` covers the selection-type rules (no selection, validation argmax, runner-up,
non-reference argmax). Top-k ensembles and the resampling remedies need per-candidate risk
scores rather than C-index matrices, and are implemented in the accompanying analysis scripts
under `scripts/`.

The package ships assertion-based tests covering column-order invariance of tie-breaking,
row-identity checking, rejection of out-of-range C-indices, and the identity between the
comparable-pair count and `concordant + discordant + tied.x`.

## Reproducing the study

**Prerequisites.** `tables` and `metabric` need inputs that are not in this repository:
`results/main/19_*_100.rds` (saved per-patient predictions, produced by `train`) and
`data/cohorts/` (expression and survival tables, see `data/README.md`). Without them the code
and the released result tables are available, but the tables cannot be rebuilt from scratch.

```bash
bash reproduce_current.sh tables      # rebuild TCGA tables from saved predictions
bash reproduce_current.sh metabric    # METABRIC: refit, sensitivity, figure
bash reproduce_current.sh signatures  # published-signature module
bash reproduce_current.sh train       # main training run + simulation
bash reproduce_current.sh check       # unit tests + R CMD check + hash freeze
```

Data are not redistributed; see `data/README.md` for access instructions and
`scripts/36_verify_metabric_predictions.R` for an independent reconstruction check.

## Repository layout

```
R/                     package source
man/, NAMESPACE         package documentation and exports
tests/                  assertion-based tests
scripts/                analysis pipeline (numbered by stage)
config/                 frozen analysis plans (fixed before the runs)
results/figures         final figures and figure data
results/tables          result tables used in the manuscript
data/README.md          how to obtain the public data
```

## Citation

If you use this software, please cite the accompanying article (under review). A
`CITATION.cff` file is included; it will be updated with the DOI on publication.

## License

MIT — see `LICENSE`.

## Contact

Ziming Zhao, College of Chemistry, Fuzhou University, Fuzhou, Fujian, China — 13317834871@163.com

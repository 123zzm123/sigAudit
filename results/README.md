# Results

- `figures/` — final figures (PDF and 300-dpi PNG) and the figure data behind them.
  File names carry a historical offset relative to the manuscript numbering:
  `fig1_design` (Fig. 1), `fig2_gain_decomposition` (Fig. 2), `fig4_rule_pairing` (Fig. 3),
  `fig5_metabric` (Fig. 4), `fig5_search_budget` (Fig. 5, revision analysis), while
  `fig3_cross_split` is Supplementary Fig. S1 and `figs2_*`, `figs3_*`, `figs4_*` and `figs5_*`
  are Supplementary Figs. S2–S5 (search-budget curves per TCGA cohort and per METABRIC
  signature, for the LASSO and ridge starts respectively). The matching data files are
  `figdata_*.tsv` (`figdata_search_budget_*.tsv` for the search-budget figures).
- `tables/` — result tables used in the manuscript:
  - `20_summary.tsv` — rule-level summary (15 rows)
  - `20_rule_performance_long.tsv` — task-level record for every rule (58,785 rows)
  - `20_model_metrics_long.tsv` — candidate-model evaluations (66,483 rows)
  - `21_entry_effect.tsv`, `21_entry_summary.tsv` — entry effects
  - `25_simulation_final_replicates.tsv` — simulation replicates
  - `26_signature_tasks.tsv`, `26_signature_summary.tsv` — published-signature module
  - `29_baseline_sensitivity.tsv` — three fixed references
  - `47_*` — non-degenerate-reference stratification (all vs non-constant LASSO)
  - `48_*` — validation-prediction duplicate-collapse sensitivity
  - `49_*` — prediction-level interface consistency check and reference switch
  - `50_*` — reference-anchored search-budget curves (expected Cv/Ct by pool size K)
  - `52_*` — censoring-weighted (Uno/IPCW) concordance at a fixed five-year horizon
  - `53_*` — signature-matched comparison (TCGA-BRCA and METABRIC, same six signatures)
  - `54_*` — compact summary tables assembled from the above without recomputation
  - `31_*`, `32_*`, `36_*` — METABRIC metrics, sensitivity and re-verification
- `manifest/` — candidate-pool diversity tables.

Numbers in the manuscript are generated from these files; see `scripts/` for the pipeline.

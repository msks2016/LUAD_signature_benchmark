# LUAD lncRNA signature benchmark - v2 (clean rebuild)

Starter kit for rebuilding the analysis from raw public data. It replaces every analysis
script of v0.3, which stays archived (read-only) as the audit trail of the corrections.

**Status of this kit:** the functions and tests were written and desk-checked but could not
be executed where they were prepared (no R there). Run the tests first (step 2 below); if
any test fails, fix that before touching data.

## Folder layout

```
config/analysis.yaml        frozen decisions (copy them into the OSF amendment)
config/cohorts.yaml         GEO clinical field names (every CHECK line must be confirmed)
config/gene_overrides.csv   manual Ensembl IDs, each with its evidence
data/signatures/            signature_master.csv and independence.csv (you fill these)
data/raw/                   downloads (not in git)
data/processed/             one .rds per cohort (not in git; deposit on Zenodo at the end)
data/reference/             gene tables and the GPL570 probe map
R/                          functions only - nothing runs when sourced
scripts/                    numbered steps; each reads only earlier outputs
tests/testthat/             unit tests; need no data
results/                    tables, figures, numbers.csv (every manuscript number)
run_all.R                   runs the tests, then every step in order
```

## Getting started

1. Open the project in RStudio (R 4.6.x) and run `scripts/00_setup.R` once.
2. Run the tests: `testthat::test_dir("tests/testthat")`. All must pass.
3. Copy `data/signatures/signature_master_TEMPLATE.csv` to `signature_master.csv`.
   Column `coef_extractor1` holds the coefficients from v0.3. A second person fills
   `coef_extractor2` blind from the PDFs; then fill `coefficient` (agreed value),
   `annotation` (`v22` or `v36`: the gene model the paper used), `source_location`.
4. Copy `independence_TEMPLATE.csv` to `independence.csv` and replace every `to_check`
   with `yes` or `no` from each paper's methods (development, feature selection, validation).
5. Run the scripts in the order listed in `run_all.R`, checking each output before moving on.

## Working rules

- One scoring function (`score_signature()`); no script computes a score any other way.
- No number is typed into the manuscript: every number comes from `results/numbers.csv`.
- No decision in `config/analysis.yaml` changes after the OSF amendment is filed; a later
  change is a new dated deviation.
- Commit after every working step with a message that says what changed and why.
- Never edit files in `data/raw/`; never hand-edit anything in `results/`.

## What each step produces

| Script | Output |
|---|---|
| 01_signatures.R | results/tables/table1_signatures.csv |
| 10-12 (13 optional) | data/processed/<cohort>.rds |
| 19_cohort_qc.R | results/qc/cohort_qc.csv, sample_flow.csv, expression_boxplots.pdf |
| 21-22 | data/reference/gene_tables.rds, data/signatures/gene_map.csv |
| 23_coverage.R | results/tables/coverage_matrix.csv, data/processed/measurable_ids.rds |
| 30_h1_null.R | results/tables/h1_null_test.csv |
| 31_h2_stage.R | results/tables/h2_added_value.csv |
| 32_h3_cross_cohort.R | results/tables/h3_per_pair.csv, h3_pooled_exploratory.csv |
| 33_transportability.R | results/tables/transportability_per_pair.csv, transportability_pooled.csv |
| 34_sensitivity.R | results/tables/sensitivity_*.csv |

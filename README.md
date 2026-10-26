# Benchmark of published lncRNA prognostic signatures for lung adenocarcinoma

Code for the study *Published lncRNA prognostic signatures for lung adenocarcinoma share a modest
transportable signal, but none is individually validated: a pre-registered benchmark*
(Karim MS, Ma Z, Li Y; Shanghai University).

Fifteen published lncRNA prognostic signatures for lung adenocarcinoma (LUAD) are applied with their
published genes and coefficients, without refitting, to TCGA-LUAD (apparent performance) and six
independent public cohorts (CPTAC-3, GSE31210, GSE50081, GSE30219, GSE81089, GSE72094). The code
rebuilds every cohort from raw public data, maps every signature gene to an Ensembl identifier,
and writes every number reported in the manuscript to `results/numbers.csv`.

Registration: Open Science Framework, [osf.io/qckvr](https://osf.io/qckvr) (6 October 2026; Updates 1-4).

## Versions and tags

This repository holds two versions of the analysis code. Only the v2 pipeline produces the results
in the manuscript.

| Tag or commit | Version | Meaning |
|---|---|---|
| `v2.0-submission` | v2 | Code and outputs submitted with the manuscript |
| `v2-preanalysis` (commit 871d199) | v2 | Inputs frozen (signatures, gene mapping, cohorts, coverage) before any v2 outcome analysis, 24 Oct 2026 |
| commit 229fc8f | v0.3 | Last state of the earlier code, kept unchanged as the record of the corrections (OSF Update 3 and 4) |
| `v0.3-corrected-pipeline`, `v1.0`, `v1.1-submission`, `v2.0` to `v2.6` | v0.3 | Tags of the earlier code. Despite their names, `v2.0` to `v2.6` belong to the earlier code, not to the v2 pipeline |

The files of the earlier code are in `archive_v1/`. No v2 script reads them, and renv ignores them
(`.renvignore`). Why the code was rebuilt, and every change since the registration, are recorded in
`docs/deviations_v2.md` and in the OSF updates.

## Repository layout

```
config/analysis.yaml            analysis settings (coverage rule, seeds, resamples, decision rules)
config/cohorts.yaml             cohort sources and GEO clinical field names
config/gene_overrides.csv       six manual Ensembl mappings, each with its evidence
config/gene_id_alternatives.csv successor Ensembl IDs from the official ID history
data/signatures/                signature genes and coefficients (double extraction), independence table
data/raw/, data/processed/      downloads and one .rds per cohort (not in git)
data/reference/                 gene tables and the GPL570 probe map
R/                              functions only; nothing runs when sourced
scripts/                        numbered steps; each reads only earlier outputs
tests/testthat/                 unit tests; need no data
results/tables/                 analysis tables
results/paper/                  manuscript tables, figures and the supplementary workbook
results/numbers.csv             every number reported in the manuscript
docs/deviations_v2.md           deviations log of the v2 pipeline
docs/coding/                    reporting-completeness coding of the 76 papers (see its README)
docs/osf/                       text of OSF Update 4
docs/reporting_checklist.md     reporting checklist (Supplementary File 2)
archive_v1/                     files of the earlier code (v0.3); not used
run_all.R                       runs the tests, then every step in order
```

## Reproducing the analysis

Requirements: R 4.6.1 (RStudio recommended) and an internet connection for the downloads. Package
versions are recorded in `renv.lock`.

1. Open the project in RStudio and restore the packages:
   ```r
   renv::restore()
   ```
2. Run the unit tests (they need no data):
   ```r
   testthat::test_dir("tests/testthat")
   ```
3. Run the pipeline, either all at once with `source("run_all.R")` or step by step in the order
   listed in `run_all.R`. Downloads are cached in `data/raw/`; the bootstrap steps (scripts 33 and
   34) take the longest.

Raw data are not redistributed. The scripts download them from the NCI Genomic Data Commons (TCGA-LUAD
and CPTAC-3, GDC Data Release 46.0) and the Gene Expression Omnibus (GSE31210, GSE50081, GSE30219 as
raw CEL files; GSE81089 FPKM file; GSE72094 series matrix). All random draws use fixed seeds (base
seed 2026).

## Steps and outputs

| Script | What it does | Main output |
|---|---|---|
| 00_setup.R | Package setup with renv | renv.lock |
| 01_signatures.R | Checks the double-extracted signatures | results/tables/table1_signatures.csv |
| 10-14 | Builds each cohort from raw data | data/processed/\<cohort\>.rds |
| 19_cohort_qc.R | Cohort quality checks and patient flow | results/qc/ |
| 21-22 (22b, 22c) | Gene tables, Ensembl mapping, gene-name history and audit | data/signatures/gene_map.csv |
| 23 (23b, 23c) | Measurable genes and coverage; reasons for missing genes; predecessor check | results/tables/coverage_*.csv |
| 29_positive_control.R | Stage positive control | results/tables/positive_control.csv |
| 30_h1_null.R | H1: comparison with expression-matched random gene sets | results/tables/h1_null_test.csv |
| 31_h2_stage.R | H2: added value beyond stage, sex and age | results/tables/h2_added_value.csv |
| 32_h3_cross_cohort.R | H3: per-pair metrics and exploratory pooling | results/tables/h3_*.csv |
| 37_pooling_sensitivity.R | Pooling-interval sensitivity (post hoc) | results/tables/h3_pooled_sensitivity.csv |
| 35_direction_genes.R | Direction consistency, signature overlap, gene-level concordance | results/tables/direction_*.csv |
| 36_power.R | Power simulation at each cohort's size | results/tables/power_simulation.csv |
| 33_transportability.R | Aggregate transportability (added after registration) | results/tables/transportability_*.csv |
| 34_sensitivity.R | Sensitivity analyses | results/tables/sensitivity_*.csv |
| 40_paper_outputs.R | Manuscript tables and figures | results/paper/tables, results/paper/figures |
| 41_reporting_completeness.R | Reporting completeness of 76 papers | results/paper/tables/supp_reporting_completeness_76.csv |
| 42_supplementary_workbook.R | Supplementary Tables S1-S19 in one workbook | results/paper/Supplementary_Tables_S1-S19.xlsx |
| 43_graphical_abstract.R | Graphical abstract, drawn from the result tables | results/paper/figures/graphical_abstract.tif |

## Manuscript items

| Item | File |
|---|---|
| Tables 1-5 | results/paper/tables/table_*.csv |
| Figures 2-6 | results/paper/figures/fig_*.pdf |
| Supplementary Tables S1-S19 | results/paper/Supplementary_Tables_S1-S19.xlsx |
| Supplementary Figures S1-S2 | results/paper/figures/figS_*.pdf |
| Every number in the text | results/numbers.csv |

Figure 1 (selection of signatures) is not built by the code.

## Conventions

- One scoring function (`score_signature()`): risk score = sum of published coefficient x
  within-cohort z-scored expression.
- No number is typed into the manuscript; every number comes from `results/numbers.csv`.
- Files in `data/raw/` and `results/` are never edited by hand.
- Any change after the freeze is recorded, with its date and effect, in `docs/deviations_v2.md`.

#!/usr/bin/env bash
# End-to-end reproduction of all manuscript numbers.
# Requires: R 4.6.1, renv, raw data downloaded (see README).
# Usage: bash scripts/bash/run_all.sh

set -euo pipefail
cd "$(dirname "$0")/../.."

echo "=== 1. Restoring R environment ==="
Rscript --vanilla -e "renv::restore(prompt = FALSE)"

echo "=== 2. Preparing cohorts ==="
Rscript --vanilla scripts/R/06_prepare_TCGA.R
test -f data/processed/CPTAC3_LUAD_expr_matched.rds || (echo "ERROR: data/processed/CPTAC3_LUAD_expr_matched.rds not found."; echo "Reconstruct via scripts/build_CPTAC3_matrix.R or restore from the release archive."; exit 1)
Rscript --vanilla scripts/14_load_GSE31210.R

echo "=== 3. Scoring signatures ==="
Rscript --vanilla scripts/09_score_CPTAC3.R
Rscript --vanilla scripts/15_score_GSE31210.R

echo "=== 4. Primary metrics ==="
Rscript --vanilla scripts/10_metrics_CPTAC3.R
Rscript --vanilla scripts/16_H2_added_value.R
Rscript --vanilla scripts/17_H3_cross_cohort.R
Rscript --vanilla scripts/18_both_Cindexes.R
Rscript --vanilla scripts/23_transportability.R

echo "=== 5. Null tests and robustness ==="
Rscript --vanilla scripts/11_null_matched_CPTAC3.R
Rscript --vanilla scripts/11b_sensitivity_80pct_CPTAC3.R
Rscript --vanilla scripts/19_robustness.R
Rscript --vanilla scripts/13_power_simulation.R

echo "=== 6. Tables and figures ==="
Rscript --vanilla scripts/20_master_tables.R
Rscript --vanilla scripts/21_prisma_diagram.R
Rscript --vanilla scripts/22_extra_tables.R
Rscript --vanilla scripts/25_restructure_Table3.R


echo "=== Done. Compare results/tables to the committed versions. ==="

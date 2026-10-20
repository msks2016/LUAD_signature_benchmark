# run_all.R - the whole pipeline, in order. Run from the project root:
#   Rscript run_all.R
# Each script reads only files written by earlier scripts or by you (data/signatures, config).

testthat::test_dir("tests/testthat", stop_on_failure = TRUE)

steps <- c(
  "scripts/01_signatures.R",        # Step 1: double-extracted signatures
  "scripts/10_cohort_tcga.R",       # Step 2: cohorts (downloads; slow the first time)
  "scripts/11_cohort_cptac3.R",
  "scripts/12_cohort_gpl570.R",
  "scripts/13_cohort_gse81089.R",   # restored in the OSF amendment
  "scripts/19_cohort_qc.R",         # Step 2 gate
  "scripts/21_gene_tables.R",       # Step 3: gene identifiers
  "scripts/22_gene_map.R",
  "scripts/23_coverage.R",
  "scripts/30_h1_null.R",           # Step 5: analyses
  "scripts/31_h2_stage.R",
  "scripts/32_h3_cross_cohort.R",
  "scripts/33_transportability.R",
  "scripts/34_sensitivity.R"
)
for (s in steps) {
  message("\n==== ", s, " ====")
  source(s, local = new.env())
}
dir.create("results", showWarnings = FALSE)
writeLines(utils::capture.output(utils::sessionInfo()), "results/session_info.txt")
message("Done. Every manuscript number is in results/numbers.csv")

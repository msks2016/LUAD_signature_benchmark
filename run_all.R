# run_all.R - the whole pipeline, in order. Run from the project root:
#   Rscript run_all.R
# Each script reads only files written by earlier scripts or by you (data/signatures, config).

testthat::test_dir("tests/testthat", stop_on_failure = TRUE)

steps <- c(
  "scripts/01_signatures.R",        # Step 1: double-extracted signatures
  "scripts/10_cohort_tcga.R",       # Step 2: cohorts (downloads; slow the first time)
  "scripts/11_cohort_cptac3.R",
  "scripts/12_cohort_gpl570.R",
  "scripts/13_cohort_gse81089.R",
  if (!is.null(yaml::read_yaml("config/cohorts.yaml")$GSE72094)) "scripts/14_cohort_gse72094.R",
  "scripts/19_cohort_qc.R",         # Step 2 gate
  "scripts/21_gene_tables.R",       # Step 3: gene identifiers and coverage
  "scripts/22_gene_map.R",
  "scripts/23_coverage.R",
  "scripts/29_positive_control.R",  # Step 4: stage positive control (gates scripts 30-32)
  "scripts/30_h1_null.R",           # Step 5: registered analyses
  "scripts/31_h2_stage.R",
  "scripts/32_h3_cross_cohort.R",
  "scripts/37_pooling_sensitivity.R", # interval-method sensitivity for the pooled estimates
  "scripts/35_direction_genes.R",
  "scripts/36_power.R",
  "scripts/33_transportability.R",  # added after registration (not pre-registered)
  "scripts/34_sensitivity.R",
  "scripts/40_paper_outputs.R",     # manuscript tables and figures (results/paper)
  "scripts/41_reporting_completeness.R", # literature reporting completeness (76 papers)
  "scripts/42_supplementary_workbook.R",  # Supplementary Tables S1-S19 (Excel; needs openxlsx)
  "scripts/43_graphical_abstract.R"       # graphical abstract for the journal
)
for (s in steps) {
  message("\n==== ", s, " ====")
  source(s, local = new.env())
}
dir.create("results", showWarnings = FALSE)
writeLines(utils::capture.output(utils::sessionInfo()), "results/session_info.txt")
message("Done. Every manuscript number is in results/numbers.csv")

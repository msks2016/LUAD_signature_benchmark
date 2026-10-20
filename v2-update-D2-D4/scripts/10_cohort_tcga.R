# 10_cohort_tcga.R - TCGA-LUAD from GDC STAR-Counts (reference cohort: every signature was
# developed here, so it only shows apparent performance).
# One primary tumour per patient: the largest total read count (registered rule "highest total
# expression"; OSF Update 3, decision D4).
# The GDC data saved by the first run (data/raw/GDC/TCGA-LUAD_se.rds and _clinical_raw.rds) are
# reused, so a re-run downloads nothing and keeps the GDC data release of the first download.
# Output to check: results/qc/TCGA-LUAD_build_summary.txt

for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) source(f)
suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(SummarizedExperiment)
})
dir.create("data/raw/GDC", recursive = TRUE, showWarnings = FALSE)
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
dir.create("results/qc", recursive = TRUE, showWarnings = FALSE)

se_file <- "data/raw/GDC/TCGA-LUAD_se.rds"
clin_file <- "data/raw/GDC/TCGA-LUAD_clinical_raw.rds"
release_file <- "data/raw/GDC/TCGA-LUAD_gdc_release.txt"   # line 1 release, line 2 download date
out_file <- "data/processed/TCGA-LUAD.rds"

summary_file <- "results/qc/TCGA-LUAD_build_summary.txt"
report <- c("TCGA-LUAD build summary", paste("Run:", format(Sys.time(), "%Y-%m-%d %H:%M")))
add <- function(title, x) {
  report <<- c(report, "", paste("==", title), utils::capture.output(print(x)))
  writeLines(report, summary_file)
}

old <- if (file.exists(out_file)) readRDS(out_file) else NULL   # previous build, for comparison

if (file.exists(se_file) && file.exists(clin_file)) {
  message("Reusing the GDC data saved on the first run (nothing is downloaded)")
  se <- readRDS(se_file)
  clin_raw <- readRDS(clin_file)
  if (!file.exists(release_file)) {
    rel <- if (!is.null(old)) c(as.character(old$meta$gdc_data_release), as.character(old$meta$downloaded))
           else c("not recorded", "not recorded")
    writeLines(rel, release_file)
  }
} else {
  q <- GDCquery(project = "TCGA-LUAD", data.category = "Transcriptome Profiling",
                data.type = "Gene Expression Quantification", workflow.type = "STAR - Counts")
  GDCdownload(q, directory = "data/raw/GDC", files.per.chunk = 50)
  se <- GDCprepare(q, directory = "data/raw/GDC")
  saveRDS(se, se_file)
  clin_raw <- GDCquery_clinic("TCGA-LUAD", type = "clinical")
  saveRDS(clin_raw, clin_file)
  writeLines(c(gdc_release(), format(Sys.Date())), release_file)
}
rel <- readLines(release_file)

co <- build_gdc_cohort(se, clin_raw, name = "TCGA-LUAD")
co$meta$gdc_data_release <- rel[1]
co$meta$downloaded <- rel[2]
check_cohort(co)
saveRDS(co, out_file)
choice <- co$meta$sample_choice
utils::write.csv(choice, "data/processed/TCGA-LUAD_sample_choice.csv", row.names = FALSE)

add("GDC data release / download date", rel)
add("Sample flow (patients)", co$meta$sample_flow)
add("RNA-seq files / primary tumour samples / patients with a primary tumour sample",
    c(files = co$meta$files, primary_tumour_samples = co$meta$primary_tumour_samples,
      patients = co$meta$patients_with_rna))
add("One sample per patient (largest total read count)",
    c(patients_with_several_primary_samples = sum(choice$n_primary_samples > 1),
      chosen_sample_differs_from_first_barcode = sum(choice$changed_by_rule)))
add("Patients whose chosen sample changed (listed in data/processed/TCGA-LUAD_sample_choice.csv)",
    utils::head(choice[choice$changed_by_rule, c("patient", "barcode", "first_barcode", "total_reads")], 20))
add("Cohort", c(patients = ncol(co$expr), deaths = sum(co$event),
                median_follow_up_months = round(reverse_km_median(co$time, co$event), 1),
                genes = nrow(co$expr)))
add("Stage group", table(co$clin$stage, useNA = "ifany"))
add("Expression range, log2(TPM + 1)", round(range(co$expr), 2))
if (!is.null(old)) {
  same_patients <- identical(sort(colnames(old$expr)), sort(colnames(co$expr)))
  add("Compared with the previous build (first barcode per patient)",
      c(previous_patients = ncol(old$expr), previous_deaths = sum(old$event),
        same_patients = same_patients))
}
cat(report, sep = "\n")
message("\nTCGA-LUAD: ", ncol(co$expr), " patients, ", sum(co$event), " deaths. Summary saved to ",
        summary_file)

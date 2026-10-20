# 10_cohort_tcga.R - TCGA-LUAD from GDC STAR-Counts (reference cohort: every signature was
# developed here, so it only shows apparent performance).

for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) source(f)
suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(SummarizedExperiment)
})
dir.create("data/raw/GDC", recursive = TRUE, showWarnings = FALSE)
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)

q <- GDCquery(project = "TCGA-LUAD", data.category = "Transcriptome Profiling",
              data.type = "Gene Expression Quantification", workflow.type = "STAR - Counts")
GDCdownload(q, directory = "data/raw/GDC", files.per.chunk = 50)
se <- GDCprepare(q, directory = "data/raw/GDC")
saveRDS(se, "data/raw/GDC/TCGA-LUAD_se.rds")

clin_raw <- GDCquery_clinic("TCGA-LUAD", type = "clinical")
saveRDS(clin_raw, "data/raw/GDC/TCGA-LUAD_clinical_raw.rds")

co <- build_gdc_cohort(se, clin_raw, name = "TCGA-LUAD")
co$meta$downloaded <- format(Sys.Date())
co$meta$gdc_data_release <- gdc_release()   # GDC data release on the download day
saveRDS(co, "data/processed/TCGA-LUAD.rds")

print(co$meta$sample_flow)
message("TCGA-LUAD: ", ncol(co$expr), " patients, ", sum(co$event), " deaths")

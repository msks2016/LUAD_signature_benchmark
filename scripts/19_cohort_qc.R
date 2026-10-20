# 19_cohort_qc.R - Step 2 gate: every cohort passes check_cohort() and gets a QC row.

for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) source(f)
cfg <- load_config()
cohort_names <- c(cfg$cohorts$reference, cfg$cohorts$independent, cfg$cohorts$candidate)
dir.create("results/qc", recursive = TRUE, showWarnings = FALSE)

# Gene types from any GDC STAR-Counts file (GENCODE v36), for the TPM-composition check.
star_files <- list.files("data/raw/GDC", pattern = "star_gene_counts\\.tsv$", recursive = TRUE,
                         full.names = TRUE)
keep_ids <- character(0)
if (length(star_files)) {
  gt <- read_star_gene_table(star_files[1])
  keep_ids <- gt$gene_id[gt$gene_type %in% cfg$sensitivity$tpm_rescale_gene_types]
}

qc <- list()
flow <- list()
grDevices::pdf("results/qc/expression_boxplots.pdf", width = 11, height = 5)
for (n in cohort_names) {
  co <- readRDS(file.path("data/processed", paste0(n, ".rds")))
  check_cohort(co)
  st <- co$clin$stage
  share <- if (identical(co$meta$platform, "RNA-seq") && length(keep_ids)) {
    tpm_share_outside(co$expr, keep_ids)
  } else NA_real_
  qc[[n]] <- data.frame(
    cohort = n, platform = co$meta$platform, patients = ncol(co$expr), deaths = sum(co$event),
    median_followup_months = round(reverse_km_median(co$time, co$event), 1),
    stage_I_II_pct = round(100 * sum(st %in% c("I", "II")) / sum(!is.na(st)), 1),
    stage_missing = sum(is.na(st)),
    age_recorded_pct = round(100 * mean(!is.na(co$clin$age)), 1),
    expr_min = round(min(co$expr), 2), expr_max = round(max(co$expr), 2),
    tpm_share_other_rna_median = round(stats::median(share), 3),
    tpm_share_other_rna_max = round(max(share), 3),
    stringsAsFactors = FALSE)
  flow[[n]] <- data.frame(cohort = n, step = names(co$meta$sample_flow),
                          n = as.integer(co$meta$sample_flow), stringsAsFactors = FALSE)
  show <- seq_len(min(60, ncol(co$expr)))
  graphics::boxplot(co$expr[, show], outline = FALSE, las = 2, cex.axis = 0.5,
                    main = paste(n, "- first", length(show), "patients"), ylab = "log2 expression")
  record_number(paste0("patients_", n), ncol(co$expr))
  record_number(paste0("deaths_", n), sum(co$event))
}
grDevices::dev.off()

qc <- do.call(rbind, qc)
utils::write.csv(qc, "results/qc/cohort_qc.csv", row.names = FALSE)
utils::write.csv(do.call(rbind, flow), "results/qc/sample_flow.csv", row.names = FALSE)
print(qc)

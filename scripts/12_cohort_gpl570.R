# 12_cohort_gpl570.R - GSE31210, GSE50081 and GSE30219 from raw CEL files.
# Same steps for all three: select samples -> RMA per cohort -> one probe map and one probe
# choice shared by the three cohorts -> Ensembl gene-level matrix.
# Clinical fields come from config/cohorts.yaml (confirmed against logs/geo_fields.txt).
# The CEL archives are large (about 1-3 GB each); if a download stops, run the script again.
# Output to check: results/qc/GEO_build_summary.txt

for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) source(f)
suppressPackageStartupMessages({
  library(GEOquery)
  library(affy)
})
fields <- yaml::read_yaml("config/cohorts.yaml")
gses <- c("GSE31210", "GSE50081", "GSE30219")
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
dir.create("results/qc", recursive = TRUE, showWarnings = FALSE)

summary_file <- "results/qc/GEO_build_summary.txt"
report <- c("GEO (GPL570) build summary", paste("Run:", format(Sys.time(), "%Y-%m-%d %H:%M")))
add <- function(title, x) {       # add one block to the summary file (rewritten each time)
  report <<- c(report, "", paste("==", title), utils::capture.output(print(x)))
  writeLines(report, summary_file)
}
count_table <- function(...) {     # cross-table as a list of the combinations that occur
  d <- as.data.frame(table(..., useNA = "ifany"), responseName = "n")
  d[d$n > 0, , drop = FALSE]
}

probe_level <- list()
clin <- list()
for (gse in gses) {
  f <- fields[[gse]]
  pd <- geo_pheno(gse)
  if (is.null(pd[[f$include$field]])) stop(gse, ": include field not found: ", f$include$field)
  keep <- tolower(trimws(pd[[f$include$field]])) %in% tolower(f$include$values)
  add(paste(gse, "- samples selected / all"), c(selected = sum(keep), all = nrow(pd)))
  pd <- pd[keep, , drop = FALSE]
  clin[[gse]] <- make_clin_geo(pd, f)
  cl <- clin[[gse]]
  add(paste(gse, "- event (1 = dead) by raw value"),
      count_table(raw = pd[[f$event$field]], event = cl$os_event))
  add(paste(gse, "- stage: raw value -> group"), count_table(raw = cl$stage_raw, group = cl$stage))
  add(paste(gse, "- follow-up in months (summary)"), summary(cl$os_months))
  add(paste(gse, "- age recorded / sex"), list(age_recorded = sum(!is.na(cl$age)),
                                                sex = table(cl$sex, useNA = "ifany")))
  if (!is.null(cl$sensitivity_exclude)) {
    add(paste(gse, "- flagged for the sensitivity analysis"), sum(cl$sensitivity_exclude))
  }
  message(gse, ": ", sum(keep), " samples selected; reading CEL files ...")
  probe_level[[gse]] <- rma_from_cel(gse, pd$geo_accession)
}

probes <- Reduce(intersect, lapply(probe_level, rownames))
pm <- gpl570_probe_map(probes)
shared_means <- rowMeans(sapply(probe_level, function(m) rowMeans(m[probes, , drop = FALSE])))

for (gse in gses) {
  col <- collapse_probes(probe_level[[gse]][probes, , drop = FALSE], pm, shared_means)
  meta <- list(source = "GEO CEL files; affy::justRMA on the selected samples of this cohort",
               probe_map = "data/reference/gpl570_probe_map_ensembl_primary.csv (Ensembl BioMart, release 114, primary-assembly genes)",
               probe_choice = col$probe_choice, downloaded = format(Sys.Date()))
  co <- finalise_cohort(col$expr, clin[[gse]], gse, platform = "GPL570", meta = meta)
  check_cohort(co)
  co$probe_expr <- probe_level[[gse]][, co$clin$patient_id, drop = FALSE]   # for diagnostics
  saveRDS(co, file.path("data/processed", paste0(gse, ".rds")))
  add(paste(gse, "- sample flow (patients)"), co$meta$sample_flow)
  add(paste(gse, "- cohort"), c(patients = ncol(co$expr), deaths = sum(co$event),
                                median_follow_up_months = round(reverse_km_median(co$time, co$event), 1),
                                max_follow_up_months = round(max(co$time), 1), genes = nrow(co$expr)))
  add(paste(gse, "- stage group"), table(co$clin$stage, useNA = "ifany"))
  add(paste(gse, "- expression range"), round(range(co$expr), 2))
}
add("GPL570 probes: on the array / mapped to one Ensembl gene / genes kept",
    c(probes = length(probes), probe_gene_pairs = nrow(pm), genes = nrow(col$probe_choice)))
utils::write.csv(col$probe_choice, "data/reference/gpl570_probe_choice.csv", row.names = FALSE)
cat(report, sep = "\n")
message("\nSummary saved to ", summary_file)

# 23_coverage.R - Step 3c: which signature genes are measurable in each cohort.

for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) source(f)
cfg <- load_config()
gm <- utils::read.csv("data/signatures/gene_map.csv", stringsAsFactors = FALSE)
cohort_names <- c(cfg$cohorts$reference, cfg$cohorts$independent, cfg$cohorts$candidate)
cohorts <- load_cohorts(cohort_names)

measurable <- lapply(cohorts, function(co) {
  frac <- if (identical(co$meta$platform, "RNA-seq")) cfg$presence$rnaseq_min_expressed_fraction else 0
  measurable_rows(co$expr, min_sd = cfg$presence$min_sd, min_expressed_frac = frac)
})
saveRDS(measurable, "data/processed/measurable_ids.rds")

cov <- do.call(rbind, lapply(names(cohorts), function(n) coverage_table(gm, measurable[[n]], n)))
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
utils::write.csv(cov, "results/tables/coverage_long.csv", row.names = FALSE)
wide <- stats::reshape(cov[, c("signature", "cohort", "frac_present")], idvar = "signature",
                       timevar = "cohort", direction = "wide")
utils::write.csv(wide, "results/tables/coverage_matrix.csv", row.names = FALSE)

# Same platform, same probe map: coverage must agree (differences can only come from
# genes that are constant in one cohort - inspect any warning).
gpl <- intersect(c("GSE31210", "GSE50081", "GSE30219"), names(cohorts))
if (length(gpl) > 1) {
  ref <- cov$frac_present[cov$cohort == gpl[1]]
  for (g in gpl[-1]) {
    if (!isTRUE(all.equal(ref, cov$frac_present[cov$cohort == g]))) {
      warning("Coverage differs between ", gpl[1], " and ", g, ": inspect coverage_long.csv")
    }
  }
}

for (n in names(cohorts)) {
  fr <- cov$frac_present[cov$cohort == n]
  record_number(paste0("evaluable_primary_", n), sum(fr >= cfg$coverage$primary_min_fraction - 1e-9))
  record_number(paste0("evaluable_100pct_", n), sum(fr >= 1 - 1e-9))
}
print(wide)

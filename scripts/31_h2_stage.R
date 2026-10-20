# 31_h2_stage.R - H2: added value beyond stage (and age when >= 80% recorded), in every
# independent cohort with enough deaths. Benjamini-Hochberg within cohort.
# Classification below uses the BH-adjusted likelihood-ratio test; replace it with the exact
# H2 criterion in your registration if that differs.

for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) source(f)
cfg <- load_config()
gm <- utils::read.csv("data/signatures/gene_map.csv", stringsAsFactors = FALSE)
ind <- utils::read.csv("data/signatures/independence.csv", stringsAsFactors = FALSE)
meas <- readRDS("data/processed/measurable_ids.rds")
cohorts <- load_cohorts(cfg$cohorts$independent)
P <- build_pairs(gm, meas, ind, names(cohorts), cfg$coverage$primary_min_fraction)

out <- list()
for (k in seq_len(nrow(P$pairs))) {
  p <- P$pairs[k, ]
  co <- cohorts[[p$cohort]]
  if (sum(co$event) < cfg$cohorts$min_events_h2) next
  g <- P$genes[[p$pair_id]]
  s <- score_signature(co$expr, g$ensembl_id, g$coefficient, cfg$scoring$scale)
  use_age <- mean(!is.na(co$clin$age)) >= 0.8
  r <- h2_added_value(co$time, co$event, s, co$clin$stage, if (use_age) co$clin$age else NULL)
  out[[p$pair_id]] <- data.frame(p, r, stringsAsFactors = FALSE)
}
res <- do.call(rbind, out)
res$p_bh <- stats::ave(res$p_lrt, res$cohort, FUN = function(x) stats::p.adjust(x, "BH"))
res$h2_result <- ifelse(res$p_bh < 0.05, "Adds value", "No added value")

utils::write.csv(res, "results/tables/h2_added_value.csv", row.names = FALSE)
for (cn in unique(res$cohort)) {
  record_number(paste0("h2_tested_", cn), sum(res$cohort == cn))
  record_number(paste0("h2_adds_value_", cn), sum(res$cohort == cn & res$h2_result == "Adds value"))
}
print(res[, c("signature", "cohort", "covariates", "c_base", "c_full", "delta_c", "p_lrt", "p_bh", "hr_adj")])

# 31_h2_stage.R - H2: added value beyond stage, sex and age, in every independent cohort.
# Model A: stage (I, II, III, IV) + sex + age (age only where recorded for >= 80% of patients);
# Model B: Model A + scale(risk score); likelihood-ratio test on 1 df; Holm correction within
# cohort (decision D5). Delta C-index with a 1,000-resample bootstrap 95% CI.
# Output: results/tables/h2_added_value.csv

for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) source(f)
cfg <- load_config()
gm <- utils::read.csv("data/signatures/gene_map.csv", stringsAsFactors = FALSE)
ind <- utils::read.csv("data/signatures/independence.csv", stringsAsFactors = FALSE)
meas <- readRDS("data/processed/measurable_ids.rds")
gate <- cohort_gate()
cohorts <- load_cohorts(cfg$cohorts$independent)
P <- build_pairs(gm, meas, ind, names(cohorts), cfg$coverage$primary_min_fraction)
if (is.null(P$pairs)) stop("No independent signature x cohort pair passes the coverage rule")

out <- lapply(seq_len(nrow(P$pairs)), function(k) {
  p <- P$pairs[k, ]
  co <- cohorts[[p$cohort]]
  g <- P$genes[[p$pair_id]]
  s <- score_signature(co$expr, g$ensembl_id, g$coefficient, cfg$scoring$scale)
  age <- if (use_age_in_h2(co$clin$age, cfg$h2$age_min_recorded)) co$clin$age else NULL
  r <- h2_added_value(co$time, co$event, s, co$clin$stage, sex = co$clin$sex, age = age,
                      B = cfg$bootstrap$cindex_resamples, seed = cfg$bootstrap$seed + k)
  data.frame(p, r, stringsAsFactors = FALSE)
})
res <- do.call(rbind, out)
res$p_holm <- stats::ave(res$p_lrt, res$cohort,
                         FUN = function(x) stats::p.adjust(x, cfg$h2$multiple_testing))
res$h2_result <- ifelse(res$p_holm < cfg$alpha, "Adds value", "No added value")
res$cohort_status <- unname(gate[res$cohort])

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
utils::write.csv(res, "results/tables/h2_added_value.csv", row.names = FALSE)
for (cn in unique(res$cohort)) {
  record_number(paste0("h2_tested_", cn), sum(res$cohort == cn))
  record_number(paste0("h2_adds_value_", cn), sum(res$cohort == cn & res$h2_result == "Adds value"))
  record_number(paste0("h2_covariates_", cn), res$covariates[res$cohort == cn][1])
}
print(res[, c("signature", "cohort", "covariates", "n", "c_base", "c_full", "delta_c",
              "p_lrt", "p_holm", "hr_adj", "h2_result")])

# 37_pooling_sensitivity.R - sensitivity of the exploratory pooled estimates (script 32) to the
# interval method. The registered REML + Hartung-Knapp interval becomes narrower than the
# common-effect interval when the cohorts agree more than expected by chance (q < 1), as for
# Yao 2021 (q = 0.19). Reported beside it (not pre-registered): common-effect, REML Wald and
# modified Hartung-Knapp (q* = max(1, q)) intervals, for the full set and leave-one-cohort-out.
# Run after script 32. Output: results/tables/h3_pooled_sensitivity.csv

for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) source(f)
cfg <- load_config()
pp <- utils::read.csv("results/tables/h3_per_pair.csv", stringsAsFactors = FALSE)
pp <- pp[pp$role == "independent", ]
n_coh <- table(pp$signature)
sigs <- names(n_coh)[n_coh >= cfg$pooling$min_independent_cohorts]
out <- do.call(rbind, lapply(sigs, function(s) {
  d <- pp[pp$signature == s, ]
  sets <- c(list(all = seq_len(nrow(d))),
            if (nrow(d) > cfg$pooling$min_independent_cohorts)
              stats::setNames(lapply(seq_len(nrow(d)), function(i) seq_len(nrow(d))[-i]),
                              paste("without", d$cohort)))
  do.call(rbind, lapply(names(sets), function(nm) {
    i <- sets[[nm]]
    data.frame(signature = s, cohorts = nm, pool_ci_variants(d$log_hr[i], d$se[i]),
               stringsAsFactors = FALSE)
  }))
}))
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
utils::write.csv(out, "results/tables/h3_pooled_sensitivity.csv", row.names = FALSE)
for (s in sigs) {
  m <- out[out$signature == s & out$cohorts == "all" & out$interval == "reml_hk_modified", ]
  record_number(paste0("h3_pooled_hk_modified_", s),
                sprintf("%.2f (%.2f-%.2f)", m$hr, m$hr_low, m$hr_high))
}
print(out[out$cohorts == "all", c("signature", "interval", "k", "q", "hr", "hr_low", "hr_high")], row.names = FALSE)

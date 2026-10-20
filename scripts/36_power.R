# 36_power.R - power simulation (registered; Update 3: at each independent cohort's own numbers
# of patients and deaths). For a score with true C-index 0.60, 0.65 or 0.70 (proportional
# hazards, normal score), the share of simulated cohorts in which C = 0.5 is rejected:
# at alpha = 0.05 and at the strictest Holm step (0.05 / evaluable signatures in the cohort).
# The H1 test compares with random genes, whose null spread is wider than chance alone, so
# these figures are upper bounds for H1. Run after script 30 (it reads the evaluable counts).
# Output: results/tables/power_simulation.csv

for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) source(f)
cfg <- load_config()
cohorts <- load_cohorts(cfg$cohorts$independent)
h1 <- if (file.exists("results/tables/h1_null_test.csv")) {
  utils::read.csv("results/tables/h1_null_test.csv", stringsAsFactors = FALSE)
} else NULL

res <- do.call(rbind, lapply(names(cohorts), function(cn) {
  co <- cohorts[[cn]]
  m <- if (is.null(h1)) 1 else max(1, sum(h1$cohort == cn))
  alphas <- unique(c(cfg$alpha, cfg$alpha / m))
  do.call(rbind, lapply(seq_along(cfg$power$target_c), function(i) {
    r <- power_simulation(n = length(co$time), deaths = sum(co$event), target_c = cfg$power$target_c[i],
                          alphas = alphas, B = cfg$power$replicates, seed = cfg$power$seed + i)
    data.frame(cohort = cn, evaluable_signatures = m, r, stringsAsFactors = FALSE)
  }))
}))
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
utils::write.csv(res, "results/tables/power_simulation.csv", row.names = FALSE)
for (i in which(res$alpha == cfg$alpha)) {
  record_number(sprintf("power_%s_C%.2f", res$cohort[i], res$target_c[i]), fmt_num(res$power[i], 2))
}
print(res)

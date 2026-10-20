# 33_transportability.R - aggregate transportability: mean log-HR per SD over independent pairs,
# patient-level cluster bootstrap, decision rule with delta = log(1.25). Added after registration
# (rules written in the v0.3 deviations log on 11 Oct 2026, declared in OSF Update 3): report it
# as NOT pre-registered, separately from H1-H3.

for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) source(f)
cfg <- load_config()
gm <- utils::read.csv("data/signatures/gene_map.csv", stringsAsFactors = FALSE)
ind <- utils::read.csv("data/signatures/independence.csv", stringsAsFactors = FALSE)
meas <- readRDS("data/processed/measurable_ids.rds")
P <- build_pairs(gm, meas, ind, cfg$cohorts$independent, cfg$coverage$primary_min_fraction)
cohorts <- load_cohorts(unique(P$pairs$cohort))

tr <- transport_bootstrap(cohorts, P$pairs, P$genes, B = cfg$bootstrap$resamples,
                          scale = cfg$scoring$scale, delta = cfg$transportability$delta,
                          seed = cfg$bootstrap$seed)

per_pair <- tr$per_pair
per_pair$se <- vapply(seq_len(nrow(P$pairs)), function(k) {
  p <- P$pairs[k, ]
  g <- P$genes[[p$pair_id]]
  co <- cohorts[[p$cohort]]
  hr_per_sd(co$time, co$event,
            score_signature(co$expr, g$ensembl_id, g$coefficient, cfg$scoring$scale))[["se"]]
}, numeric(1))

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
utils::write.csv(per_pair, "results/tables/transportability_per_pair.csv", row.names = FALSE)
pooled <- data.frame(n_pairs = tr$n_pairs, n_cohorts = tr$n_cohorts, pooled_log_hr = tr$pooled,
                     ci_low = tr$ci_low, ci_high = tr$ci_high, verdict = tr$verdict,
                     failed_replicates = tr$failed_replicates, resamples = cfg$bootstrap$resamples)
utils::write.csv(pooled, "results/tables/transportability_pooled.csv", row.names = FALSE)
saveRDS(tr$boot, "results/transportability_bootstrap_matrix.rds")

rc <- robust_crosscheck(per_pair)
if (!is.null(rc)) utils::write.csv(rc, "results/tables/transportability_robust_crosscheck.csv", row.names = FALSE)

record_number("transport_n_pairs", tr$n_pairs)
record_number("transport_pooled_log_hr", fmt_num(tr$pooled, 3))
record_number("transport_ci", paste(fmt_num(tr$ci_low, 3), fmt_num(tr$ci_high, 3), sep = " to "))
record_number("transport_verdict", tr$verdict)
print(pooled)

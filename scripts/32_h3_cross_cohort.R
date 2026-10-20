# 32_h3_cross_cohort.R - every registered per-pair metric, and H3 (descriptive, as registered).
# Pairs: TCGA-LUAD (apparent performance; every signature was developed there) and every
# independent evaluable pair. Per pair: C-index of the published score (primary) with the
# Cox-model C-index beside it, Uno's C-index at 60 months, HR per SD, time-dependent AUC at 12
# and 36 months, log-rank test at the median score, proportional-hazards test and, where it
# fails, the RMST difference between the median-split groups.
# Pooling (exploratory): only signatures evaluable in >= 3 independent cohorts (decision D6),
# REML + Hartung-Knapp with leave-one-cohort-out; TCGA-LUAD is never pooled.
# Output: results/tables/h3_per_pair.csv, h3_pooled_exploratory.csv, h3_leave_one_out.csv,
#         results/figures/h3_cindex_by_cohort.pdf

for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) source(f)
cfg <- load_config()
gm <- utils::read.csv("data/signatures/gene_map.csv", stringsAsFactors = FALSE)
ind <- utils::read.csv("data/signatures/independence.csv", stringsAsFactors = FALSE)
meas <- readRDS("data/processed/measurable_ids.rds")
gate <- cohort_gate()
primary <- cfg$coverage$primary_min_fraction

P_ref <- build_pairs(gm, meas, ind, cfg$cohorts$reference, primary, include_used = TRUE)
P_ind <- build_pairs(gm, meas, ind, cfg$cohorts$independent, primary)
add_role <- function(x, role) if (is.null(x)) NULL else cbind(x, role = role, stringsAsFactors = FALSE)
pairs_all <- rbind(add_role(P_ref$pairs, "reference"), add_role(P_ind$pairs, "independent"))
genes_all <- c(P_ref$genes, P_ind$genes)
cohorts <- load_cohorts(unique(pairs_all$cohort))
tt <- cfg$secondary$tdauc_months

per_pair <- do.call(rbind, lapply(seq_len(nrow(pairs_all)), function(k) {
  p <- pairs_all[k, ]
  co <- cohorts[[p$cohort]]
  g <- genes_all[[p$pair_id]]
  s <- score_signature(co$expr, g$ensembl_id, g$coefficient, cfg$scoring$scale)
  e <- hr_per_sd(co$time, co$event, s)
  ci <- boot_cindex_ci(co$time, co$event, s, B = cfg$bootstrap$cindex_resamples,
                       seed = cfg$bootstrap$seed + k)
  lr <- logrank_median(co$time, co$event, s)
  ph_p <- ph_test_p(co$time, co$event, s)
  rm <- if (!is.na(ph_p) && ph_p < cfg$secondary$ph_alpha) rmst_median_split(co$time, co$event, s) else NULL
  auc <- do.call(cbind, lapply(tt, function(t) {
    a <- td_auc(co$time, co$event, s, t)
    b <- boot_td_auc_ci(co$time, co$event, s, t, B = cfg$bootstrap$cindex_resamples,
                        seed = cfg$bootstrap$seed + k)
    stats::setNames(data.frame(a, b[1], b[2]),
                    paste0(c("td_auc_", "td_auc_low_", "td_auc_high_"), t, "m"))
  }))
  data.frame(p, patients = length(s), deaths = sum(co$event),
             log_hr = e[["log_hr"]], se = e[["se"]], hr = e[["hr"]],
             hr_low = e[["lo"]], hr_high = e[["hi"]],
             c_index = cindex(co$time, co$event, s), c_low = ci[1], c_high = ci[2],
             c_index_cox = cindex_cox(co$time, co$event, s),
             c_uno = cindex_uno(co$time, co$event, s, cfg$survival$uno_tau_months),
             auc,
             logrank_chisq = lr[["chisq"]], logrank_p = lr[["p"]],
             ph_p = ph_p, ph_violated = !is.na(ph_p) & ph_p < cfg$secondary$ph_alpha,
             rmst_tau = if (is.null(rm)) NA_real_ else rm[["tau"]],
             rmst_diff = if (is.null(rm)) NA_real_ else rm[["diff"]],
             rmst_low = if (is.null(rm)) NA_real_ else rm[["lo"]],
             rmst_high = if (is.null(rm)) NA_real_ else rm[["hi"]],
             rmst_p = if (is.null(rm)) NA_real_ else rm[["p"]],
             cohort_status = unname(gate[p$cohort]),
             stringsAsFactors = FALSE)
}))
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
utils::write.csv(per_pair, "results/tables/h3_per_pair.csv", row.names = FALSE)

# Exploratory pooling: >= 3 independent cohorts only; with 2 cohorts both estimates are shown
# (h3_per_pair.csv) without pooling.
ind_pp <- per_pair[per_pair$role == "independent", ]
k_min <- cfg$pooling$min_independent_cohorts
n_coh <- table(ind_pp$signature)
pooled <- list()
loo <- list()
for (s in names(n_coh)[n_coh >= k_min]) {
  d <- ind_pp[ind_pp$signature == s, ]
  pooled[[s]] <- data.frame(signature = s, pool_log_hr(d$log_hr, d$se, d$cohort),
                            deaths = sum(d$deaths), cohorts = paste(d$cohort, collapse = ";"),
                            label = "exploratory", stringsAsFactors = FALSE)
  loo[[s]] <- data.frame(signature = s, pool_leave_one_out(d$log_hr, d$se, d$cohort),
                         stringsAsFactors = FALSE)
}
utils::write.csv(if (length(pooled)) do.call(rbind, pooled) else
                   data.frame(note = paste("no signature evaluable in >=", k_min, "independent cohorts")),
                 "results/tables/h3_pooled_exploratory.csv", row.names = FALSE)
if (length(loo)) utils::write.csv(do.call(rbind, loo), "results/tables/h3_leave_one_out.csv", row.names = FALSE)
record_number("h3_signatures_pooled_exploratory", length(pooled))
record_number("h3_signatures_two_cohorts_not_pooled", sum(n_coh == 2))
for (cn in unique(ind_pp$cohort)) {
  record_number(paste0("ph_violated_", cn), sum(ind_pp$ph_violated[ind_pp$cohort == cn]))
}

# Figure: C-index of the published score with 95% CI, by signature and cohort.
if (requireNamespace("ggplot2", quietly = TRUE)) {
  dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
  pp <- per_pair
  pp$cohort <- factor(pp$cohort, levels = c(cfg$cohorts$reference, cfg$cohorts$independent))
  fig <- ggplot2::ggplot(pp, ggplot2::aes(x = c_index, y = signature, colour = cohort)) +
    ggplot2::geom_vline(xintercept = 0.5, linetype = 2, colour = "grey50") +
    ggplot2::geom_pointrange(ggplot2::aes(xmin = c_low, xmax = c_high),
                             position = ggplot2::position_dodge(width = 0.6), size = 0.3) +
    ggplot2::labs(x = "C-index of the published score (95% CI)", y = NULL, colour = "Cohort") +
    ggplot2::theme_bw()
  ggplot2::ggsave("results/figures/h3_cindex_by_cohort.pdf", fig, width = 7, height = 6)
}
print(per_pair[, c("signature", "cohort", "role", "c_index", "c_index_cox", "hr", "hr_low",
                   "hr_high", "logrank_p", "ph_p")])

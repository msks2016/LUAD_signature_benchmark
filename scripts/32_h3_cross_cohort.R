# 32_h3_cross_cohort.R - H3: per-cohort performance of every evaluable pair (descriptive, as
# registered), with TCGA-LUAD as the apparent-performance reference. Any pooled estimate is
# EXPLORATORY: the registered pooling rule excludes signatures developed in TCGA-LUAD.

for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) source(f)
cfg <- load_config()
gm <- utils::read.csv("data/signatures/gene_map.csv", stringsAsFactors = FALSE)
ind <- utils::read.csv("data/signatures/independence.csv", stringsAsFactors = FALSE)
meas <- readRDS("data/processed/measurable_ids.rds")
primary <- cfg$coverage$primary_min_fraction

P_ref <- build_pairs(gm, meas, ind, cfg$cohorts$reference, primary, include_used = TRUE)
P_ind <- build_pairs(gm, meas, ind, cfg$cohorts$independent, primary)
add_role <- function(x, role) if (is.null(x)) NULL else cbind(x, role = role, stringsAsFactors = FALSE)
pairs_all <- rbind(add_role(P_ref$pairs, "reference"), add_role(P_ind$pairs, "independent"))
genes_all <- c(P_ref$genes, P_ind$genes)
cohorts <- load_cohorts(unique(pairs_all$cohort))

per_pair <- do.call(rbind, lapply(seq_len(nrow(pairs_all)), function(k) {
  p <- pairs_all[k, ]
  co <- cohorts[[p$cohort]]
  g <- genes_all[[p$pair_id]]
  s <- score_signature(co$expr, g$ensembl_id, g$coefficient, cfg$scoring$scale)
  e <- hr_per_sd(co$time, co$event, s)
  ci <- boot_cindex_ci(co$time, co$event, s, B = cfg$bootstrap$cindex_resamples,
                       seed = cfg$bootstrap$seed + k)
  data.frame(p, patients = length(s), deaths = sum(co$event),
             log_hr = e[["log_hr"]], se = e[["se"]], hr = e[["hr"]],
             hr_low = e[["lo"]], hr_high = e[["hi"]],
             c_index = cindex(co$time, co$event, s), c_low = ci[1], c_high = ci[2],
             c_uno = cindex_uno(co$time, co$event, s, cfg$survival$uno_tau_months),
             stringsAsFactors = FALSE)
}))
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
utils::write.csv(per_pair, "results/tables/h3_per_pair.csv", row.names = FALSE)

# Exploratory pooling for signatures evaluable in >= 2 independent cohorts.
ind_pp <- per_pair[per_pair$role == "independent", ]
pooled <- list()
loo <- list()
for (s in unique(ind_pp$signature)) {
  d <- ind_pp[ind_pp$signature == s, ]
  if (nrow(d) < 2) next
  pooled[[s]] <- data.frame(signature = s, pool_log_hr(d$log_hr, d$se, d$cohort),
                            deaths = sum(d$deaths), label = "exploratory", stringsAsFactors = FALSE)
  if (nrow(d) >= 3) {
    loo[[s]] <- data.frame(signature = s, pool_leave_one_out(d$log_hr, d$se, d$cohort),
                           stringsAsFactors = FALSE)
  }
}
if (length(pooled)) utils::write.csv(do.call(rbind, pooled), "results/tables/h3_pooled_exploratory.csv", row.names = FALSE)
if (length(loo)) utils::write.csv(do.call(rbind, loo), "results/tables/h3_leave_one_out.csv", row.names = FALSE)
record_number("h3_signatures_pooled_exploratory", length(pooled))
print(per_pair[, c("signature", "cohort", "role", "hr", "hr_low", "hr_high", "c_index")])

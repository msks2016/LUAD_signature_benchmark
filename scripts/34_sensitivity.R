# 34_sensitivity.R - sensitivity analyses declared in OSF Update 3 (reported as not
# pre-registered, except the cohort-composition analyses registered for CPTAC-3, now applied
# to every cohort). Every setting changes one thing and reuses the same functions.
# Part A: aggregate transportability under each setting (B = 1000 resamples).
# Part B: H1 (C-index, expression-matched null test, Holm, classification) under the null-pool
#         settings and the cohort-composition settings (stage I-II only; early deaths
#         excluded; GSE31210 without the 22 flagged cases).
# Measurable genes are recomputed in every setting, because a gene can become constant in a
# subset of patients. Slow: allow 30-60 minutes; progress is printed.
# Output: results/tables/sensitivity_transportability.csv, sensitivity_h1.csv

for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) source(f)
cfg <- load_config()
gm <- utils::read.csv("data/signatures/gene_map.csv", stringsAsFactors = FALSE)
ind <- utils::read.csv("data/signatures/independence.csv", stringsAsFactors = FALSE)
base_cohorts <- load_cohorts(cfg$cohorts$independent)
B <- 1000   # fewer resamples than the primary transportability analysis; state this in the supplement
v36 <- readRDS("data/reference/gene_tables.rds")$v36
rescale_ids <- v36$gene_id[v36$gene_type %in% cfg$sensitivity$tpm_rescale_gene_types]

subset_all <- function(cohorts, rule) lapply(cohorts, function(co) subset_cohort(co, rule(co)))
stage_i_ii <- subset_all(base_cohorts, function(co) co$clin$stage %in% c("I", "II"))
no_early_deaths <- subset_all(base_cohorts, function(co) {
  !(co$event == 1 & co$time < cfg$survival$early_death_months)
})
gse31210_22 <- NULL
if ("GSE31210" %in% names(base_cohorts) && !is.null(base_cohorts$GSE31210$clin$sensitivity_exclude)) {
  gse31210_22 <- base_cohorts
  gse31210_22$GSE31210 <- subset_cohort(gse31210_22$GSE31210, !gse31210_22$GSE31210$clin$sensitivity_exclude)
}

# ---- Part A: aggregate transportability ---------------------------------------------------
run_transport <- function(label, cohorts = base_cohorts, min_fraction = cfg$coverage$primary_min_fraction,
                          scale = cfg$scoring$scale) {
  message("Transportability, setting: ", label)
  meas <- measurable_by_cohort(cohorts, cfg)
  P <- build_pairs(gm, meas, ind, names(cohorts), min_fraction)
  if (is.null(P$pairs)) return(data.frame(setting = label, n_pairs = 0, pooled_log_hr = NA, ci_low = NA,
                                          ci_high = NA, verdict = "no evaluable pair", failed_replicates = NA))
  tr <- transport_bootstrap(cohorts, P$pairs, P$genes, B = B, scale = scale,
                            delta = cfg$transportability$delta, seed = cfg$bootstrap$seed)
  data.frame(setting = label, n_pairs = tr$n_pairs, pooled_log_hr = tr$pooled,
             ci_low = tr$ci_low, ci_high = tr$ci_high, verdict = tr$verdict,
             failed_replicates = tr$failed_replicates, stringsAsFactors = FALSE)
}
settings <- list(
  run_transport("primary (repeated with B = 1000)"),
  run_transport("raw log2 scale", scale = cfg$scoring$sensitivity_scale),
  run_transport("100% coverage", min_fraction = 1),
  run_transport("70% coverage", min_fraction = 0.7),
  run_transport("RNA-seq TPM rescaled to protein-coding + lncRNA genes",
                cohorts = lapply(base_cohorts, rescale_tpm_cohort, keep_ids = rescale_ids)),
  run_transport("stage I-II only", cohorts = stage_i_ii),
  run_transport("deaths within 3 months excluded", cohorts = no_early_deaths)
)
if (!is.null(gse31210_22)) {
  settings[[length(settings) + 1]] <- run_transport("GSE31210 without the 22 flagged cases", cohorts = gse31210_22)
}
res <- do.call(rbind, settings)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
utils::write.csv(res, "results/tables/sensitivity_transportability.csv", row.names = FALSE)
print(res)

# ---- Part B: H1 ---------------------------------------------------------------------------
run_h1 <- function(label, cohorts = base_cohorts, pool_types = cfg$null_test$pool_gene_types) {
  meas <- measurable_by_cohort(cohorts, cfg)
  P <- build_pairs(gm, meas, ind, names(cohorts), cfg$coverage$primary_min_fraction)
  if (is.null(P$pairs)) return(NULL)
  message("H1, setting: ", label, " (", nrow(P$pairs), " pairs)")
  out <- do.call(rbind, lapply(seq_len(nrow(P$pairs)), function(k) {
    message("  pair ", k, " of ", nrow(P$pairs))
    p <- P$pairs[k, ]
    co <- cohorts[[p$cohort]]
    g <- P$genes[[p$pair_id]]
    pool <- null_pool(meas[[p$cohort]], v36, pool_types, g$ensembl_id)
    nt <- expression_matched_null(co$expr, g$ensembl_id, g$coefficient, co$time, co$event,
                                  pool = pool, scale = cfg$scoring$scale,
                                  n_draws = cfg$null_test$draws, n_bins = cfg$null_test$bins,
                                  seed = cfg$null_test$seed + k)
    data.frame(setting = label, p, patients = ncol(co$expr), deaths = sum(co$event),
               c_index = nt$observed, p_two_sided = nt$p_two_sided, stringsAsFactors = FALSE)
  }))
  out$p_holm <- stats::ave(out$p_two_sided, out$cohort, FUN = function(x) stats::p.adjust(x, "holm"))
  out$classification <- classify_h1(out$c_index, out$p_holm, alpha = cfg$alpha)
  out
}
h1 <- list(
  run_h1("null pool: lncRNA genes only", pool_types = "lncRNA"),
  run_h1("null pool: all measurable genes", pool_types = NULL),
  run_h1("stage I-II only", cohorts = stage_i_ii),
  run_h1("deaths within 3 months excluded", cohorts = no_early_deaths)
)
if (!is.null(gse31210_22)) {
  h1[[length(h1) + 1]] <- run_h1("GSE31210 without the 22 flagged cases",
                                 cohorts = gse31210_22["GSE31210"])
}
h1 <- do.call(rbind, h1)
utils::write.csv(h1, "results/tables/sensitivity_h1.csv", row.names = FALSE)
print(stats::aggregate(cbind(Replicated = classification == "Replicated",
                             Contradicted = classification == "Contradicted",
                             Inconclusive = classification == "Inconclusive") ~ setting + cohort,
                       data = h1, FUN = sum))

# Not automated here: GPL570 series matrices instead of CEL-RMA; alternative mappings of
# individual genes (edit config/gene_overrides.csv in a copy of the project).

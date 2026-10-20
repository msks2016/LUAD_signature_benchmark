# 34_sensitivity.R - pre-specified sensitivity analyses for the aggregate transportability
# result and H1. Every setting changes one thing and reuses the same functions.

for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) source(f)
cfg <- load_config()
gm <- utils::read.csv("data/signatures/gene_map.csv", stringsAsFactors = FALSE)
ind <- utils::read.csv("data/signatures/independence.csv", stringsAsFactors = FALSE)
meas <- readRDS("data/processed/measurable_ids.rds")
base_cohorts <- load_cohorts(cfg$cohorts$independent)
B <- 1000   # fewer resamples than the primary analysis; state this in the supplement
tables <- readRDS("data/reference/gene_tables.rds")
rescale_ids <- tables$v36$gene_id[tables$v36$gene_type %in% cfg$sensitivity$tpm_rescale_gene_types]

run_setting <- function(label, cohorts = base_cohorts, min_fraction = cfg$coverage$primary_min_fraction,
                        scale = cfg$scoring$scale) {
  P <- build_pairs(gm, meas, ind, names(cohorts), min_fraction)
  tr <- transport_bootstrap(cohorts, P$pairs, P$genes, B = B, scale = scale,
                            delta = cfg$transportability$delta, seed = cfg$bootstrap$seed)
  data.frame(setting = label, n_pairs = tr$n_pairs, pooled_log_hr = tr$pooled,
             ci_low = tr$ci_low, ci_high = tr$ci_high, verdict = tr$verdict,
             failed_replicates = tr$failed_replicates, stringsAsFactors = FALSE)
}

subset_all <- function(cohorts, rule) {
  lapply(cohorts, function(co) subset_cohort(co, rule(co)))
}

settings <- list(
  run_setting("primary (repeated with B = 1000)"),
  run_setting("raw log2 scale", scale = cfg$scoring$sensitivity_scale),
  run_setting("100% coverage", min_fraction = 1),
  run_setting("70% coverage", min_fraction = 0.7),
  run_setting("RNA-seq TPM rescaled to protein-coding + lncRNA genes",
              cohorts = lapply(base_cohorts, rescale_tpm_cohort, keep_ids = rescale_ids)),
  run_setting("stage I-II only", cohorts = subset_all(base_cohorts, function(co) co$clin$stage %in% c("I", "II"))),
  run_setting("deaths within 3 months excluded",
              cohorts = subset_all(base_cohorts, function(co) !(co$event == 1 & co$time < cfg$survival$early_death_months)))
)
if ("GSE31210" %in% names(base_cohorts) && !is.null(base_cohorts$GSE31210$clin$sensitivity_exclude)) {
  alt <- base_cohorts
  alt$GSE31210 <- subset_cohort(alt$GSE31210, !alt$GSE31210$clin$sensitivity_exclude)
  settings[[length(settings) + 1]] <- run_setting("GSE31210 without the 22 flagged cases", cohorts = alt)
}
res <- do.call(rbind, settings)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
utils::write.csv(res, "results/tables/sensitivity_transportability.csv", row.names = FALSE)
print(res)

# H1 with an lncRNA-only random-gene pool (gene types from GENCODE v36).
lnc_ids <- tables$v36$gene_id[tables$v36$gene_type == "lncRNA"]
P <- build_pairs(gm, meas, ind, names(base_cohorts), cfg$coverage$primary_min_fraction)
h1_lnc <- do.call(rbind, lapply(seq_len(nrow(P$pairs)), function(k) {
  p <- P$pairs[k, ]
  co <- base_cohorts[[p$cohort]]
  g <- P$genes[[p$pair_id]]
  pool <- union(intersect(meas[[p$cohort]], lnc_ids), g$ensembl_id)
  nt <- expression_matched_null(co$expr, g$ensembl_id, g$coefficient, co$time, co$event, pool = pool,
                                scale = cfg$scoring$scale, n_draws = cfg$null_test$draws,
                                n_bins = cfg$null_test$bins, seed = cfg$null_test$seed + k)
  data.frame(p, c_index = nt$observed, p_two_sided = nt$p_two_sided, stringsAsFactors = FALSE)
}))
h1_lnc$p_holm <- stats::ave(h1_lnc$p_two_sided, h1_lnc$cohort, FUN = function(x) stats::p.adjust(x, "holm"))
h1_lnc$classification <- classify_h1(h1_lnc$c_index, h1_lnc$p_holm)
utils::write.csv(h1_lnc, "results/tables/sensitivity_h1_lncRNA_pool.csv", row.names = FALSE)

# Not automated here (add if registered): GPL570 series matrices instead of CEL-RMA;
# alternative C20orf197 / MGC32805 mappings (edit config/gene_overrides.csv in a copy).

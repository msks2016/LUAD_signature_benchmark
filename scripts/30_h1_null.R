# 30_h1_null.R - H1: are published signatures better than expression-matched random genes?
# Independent pairs only. Primary C-index = that of the published score (Update 3); the
# Cox-model C-index is reported beside it. Random genes come from measurable protein-coding and
# lncRNA genes (GENCODE v36 types; decision D3) in the same mean-expression decile.
# Holm correction within each cohort over its evaluable signatures; three-way classification.
# Cohorts that failed the stage positive control (script 29) are flagged exploratory.
# Output: results/tables/h1_null_test.csv

for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) source(f)
cfg <- load_config()
gm <- utils::read.csv("data/signatures/gene_map.csv", stringsAsFactors = FALSE)
ind <- utils::read.csv("data/signatures/independence.csv", stringsAsFactors = FALSE)
meas <- readRDS("data/processed/measurable_ids.rds")
v36 <- readRDS("data/reference/gene_tables.rds")$v36
gate <- cohort_gate()
cohorts <- load_cohorts(cfg$cohorts$independent)
P <- build_pairs(gm, meas, ind, names(cohorts), cfg$coverage$primary_min_fraction)
if (is.null(P$pairs)) stop("No independent signature x cohort pair passes the coverage rule")

res <- do.call(rbind, lapply(seq_len(nrow(P$pairs)), function(k) {
  p <- P$pairs[k, ]
  co <- cohorts[[p$cohort]]
  g <- P$genes[[p$pair_id]]
  pool <- null_pool(meas[[p$cohort]], v36, cfg$null_test$pool_gene_types, g$ensembl_id)
  nt <- expression_matched_null(co$expr, g$ensembl_id, g$coefficient, co$time, co$event,
                                pool = pool, scale = cfg$scoring$scale,
                                n_draws = cfg$null_test$draws, n_bins = cfg$null_test$bins,
                                seed = cfg$null_test$seed + k)
  s <- score_signature(co$expr, g$ensembl_id, g$coefficient, cfg$scoring$scale)
  ci <- boot_cindex_ci(co$time, co$event, s, B = cfg$bootstrap$cindex_resamples,
                       seed = cfg$bootstrap$seed + k)
  data.frame(p, patients = ncol(co$expr), deaths = sum(co$event),
             c_index = nt$observed, c_low = ci[1], c_high = ci[2],
             c_index_cox = cindex_cox(co$time, co$event, s),
             null_pool_genes = length(pool) - nrow(g),
             null_q025 = nt$null_q025, null_q975 = nt$null_q975,
             p_two_sided = nt$p_two_sided, stringsAsFactors = FALSE)
}))
res$p_holm <- stats::ave(res$p_two_sided, res$cohort, FUN = function(x) stats::p.adjust(x, "holm"))
res$holm_family_size <- stats::ave(res$p_two_sided, res$cohort, FUN = length)
res$classification <- classify_h1(res$c_index, res$p_holm, alpha = cfg$alpha)
res$cohort_status <- unname(gate[res$cohort])

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
utils::write.csv(res, "results/tables/h1_null_test.csv", row.names = FALSE)

for (cn in unique(res$cohort)) {
  d <- res[res$cohort == cn, ]
  record_number(paste0("h1_evaluable_", cn), nrow(d))
  record_number(paste0("h1_replicated_", cn), sum(d$classification == "Replicated"))
  record_number(paste0("h1_contradicted_", cn), sum(d$classification == "Contradicted"))
  record_number(paste0("h1_inconclusive_", cn), sum(d$classification == "Inconclusive"))
}
print(res[, c("signature", "cohort", "c_index", "c_index_cox", "p_two_sided", "p_holm",
              "classification", "cohort_status")])

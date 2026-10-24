# 35_direction_genes.R - registered secondary analyses on signature genes:
# (1) signature overlap: Jaccard index of the mapped gene sets; clusters of signatures with
#     Jaccard > 0.3 (the clusters of the direction-consistency test);
# (2) direction consistency (secondary exploratory outcome, Update 1): sign of the HR per SD in
#     every independent evaluable pair (from script 32), binomial sign test against 50% with one
#     observation per cluster per cohort (mean log-HR of the cluster's signatures);
# (3) gene-level concordance: for every signature gene measurable in a cohort, the sign of its
#     univariable Cox coefficient versus the sign of its published coefficient.
# Run after scripts 23 and 32.
# Output: results/tables/jaccard_overlap.csv, signature_clusters.csv,
#         direction_consistency.csv, direction_consistency_observations.csv,
#         gene_concordance_long.csv, gene_concordance_by_pair.csv, gene_concordance_by_gene.csv

for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) source(f)
cfg <- load_config()
gm <- utils::read.csv("data/signatures/gene_map.csv", stringsAsFactors = FALSE)
ind <- utils::read.csv("data/signatures/independence.csv", stringsAsFactors = FALSE)
meas <- readRDS("data/processed/measurable_ids.rds")
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

# (1) Overlap and clusters ---------------------------------------------------------------------
sets <- lapply(split(gm$ensembl_id, gm$signature), function(x) unique(stats::na.omit(x)))
J <- jaccard_matrix(sets)
clusters <- jaccard_clusters(J, cfg$secondary$jaccard_cluster)
utils::write.csv(data.frame(signature = rownames(J), round(J, 3), check.names = FALSE),
                 "results/tables/jaccard_overlap.csv", row.names = FALSE)
utils::write.csv(data.frame(signature = names(clusters), cluster = unname(clusters)),
                 "results/tables/signature_clusters.csv", row.names = FALSE)
off <- J[upper.tri(J)]
record_number("jaccard_max", fmt_num(max(off), 3))
record_number("jaccard_pairs_above_threshold", sum(off > cfg$secondary$jaccard_cluster))
record_number("signature_clusters", length(unique(clusters)))

# (2) Direction consistency ----------------------------------------------------------------------
if (!file.exists("results/tables/h3_per_pair.csv")) stop("Run scripts/32_h3_cross_cohort.R first")
pp <- utils::read.csv("results/tables/h3_per_pair.csv", stringsAsFactors = FALSE)
pp <- pp[pp$role == "independent", c("signature", "cohort", "log_hr")]
dc <- direction_consistency(pp, clusters)
utils::write.csv(dc$summary, "results/tables/direction_consistency.csv", row.names = FALSE)
utils::write.csv(dc$observations, "results/tables/direction_consistency_observations.csv", row.names = FALSE)
record_number("direction_effective_n", dc$summary$n[1])
record_number("direction_effective_consistent", dc$summary$consistent[1])
record_number("direction_effective_p", fmt_num(dc$summary$p_two_sided[1], 3))
record_number("direction_nominal_n", dc$summary$n[2])
print(dc$summary)

# (3) Gene-level concordance ------------------------------------------------------------------
cohort_names <- c(cfg$cohorts$reference, cfg$cohorts$independent)
cohorts <- load_cohorts(cohort_names)
P_ind <- build_pairs(gm, meas, ind, cfg$cohorts$independent, cfg$coverage$primary_min_fraction)
evaluable <- if (is.null(P_ind$pairs)) character(0) else P_ind$pairs$pair_id
long <- do.call(rbind, lapply(cohort_names, function(cn) {
  co <- cohorts[[cn]]
  sm <- cohort_gene_map(gm, meas[[cn]])   # alternative IDs where the original is not measurable
  d <- sm[!is.na(sm$ensembl_id) & sm$ensembl_id %in% meas[[cn]], , drop = FALSE]
  if (!nrow(d)) return(NULL)
  gc <- gene_cox(co$expr, unique(d$ensembl_id), co$time, co$event)
  d <- merge(d[, c("signature", "published_symbol", "ensembl_id", "coefficient")], gc, by = "ensembl_id")
  used <- ind$used_by_authors[match(paste(d$signature, cn), paste(ind$signature, ind$cohort))]
  data.frame(cohort = cn, d, sign_agrees = sign(d$cox_beta) == sign(d$coefficient),
             used_by_authors = used,
             signature_evaluable_independent = paste(d$signature, cn, sep = "__") %in% evaluable,
             stringsAsFactors = FALSE)
}))
utils::write.csv(long, "results/tables/gene_concordance_long.csv", row.names = FALSE)
by_pair <- stats::aggregate(cbind(genes = 1, agree = sign_agrees) ~ signature + cohort +
                              used_by_authors + signature_evaluable_independent, data = long, FUN = sum)
by_pair$agreement_rate <- by_pair$agree / by_pair$genes
utils::write.csv(by_pair, "results/tables/gene_concordance_by_pair.csv", row.names = FALSE)
indep <- long[long$cohort %in% cfg$cohorts$independent & long$used_by_authors %in% "no", ]
if (nrow(indep)) {
  by_gene <- stats::aggregate(cbind(cohorts = 1, agree = sign_agrees) ~ signature + published_symbol + ensembl_id,
                              data = indep, FUN = sum)
  by_gene$agreement_rate <- by_gene$agree / by_gene$cohorts
  utils::write.csv(by_gene, "results/tables/gene_concordance_by_gene.csv", row.names = FALSE)
  record_number("gene_sign_agreement_independent", fmt_num(mean(indep$sign_agrees), 3),
                note = "share of (gene, independent cohort) pairs with matching sign")
}
print(by_pair)

# coverage.R - which signature genes count as present in a cohort.

# Rows that are measurable: no missing values, not constant, and (for RNA-seq) above
# `floor_value` in at least `min_expressed_frac` of patients. For arrays use the default 0.
measurable_rows <- function(expr, min_sd = 1e-6, min_expressed_frac = 0, floor_value = 0) {
  sdv <- apply(expr, 1, stats::sd)
  no_na <- rowSums(is.na(expr)) == 0
  expressed <- rowMeans(expr > floor_value)
  ok <- no_na & !is.na(sdv) & sdv >= min_sd & expressed >= min_expressed_frac
  rownames(expr)[ok]
}

# Coverage of every signature in one cohort.
# sig_map: one row per signature gene with columns signature, published_symbol,
#          coefficient, ensembl_id (NA if unmapped).
coverage_table <- function(sig_map, measurable_ids, cohort) {
  sig_map$present <- !is.na(sig_map$ensembl_id) & sig_map$ensembl_id %in% measurable_ids
  rows <- lapply(split(sig_map, sig_map$signature), function(d) {
    data.frame(cohort = cohort, signature = d$signature[1], n_genes = nrow(d),
               n_present = sum(d$present), frac_present = mean(d$present),
               missing = paste(d$published_symbol[!d$present], collapse = ";"),
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

# Genes and published coefficients to score for one signature-cohort pair, or NULL when the
# signature fails the coverage rule. Coefficients of absent genes are dropped, never rescaled.
pair_genes <- function(sig_map, signature, measurable_ids, min_fraction) {
  d <- sig_map[sig_map$signature == signature, , drop = FALSE]
  present <- !is.na(d$ensembl_id) & d$ensembl_id %in% measurable_ids
  if (mean(present) < min_fraction - 1e-9) return(NULL)
  d[present, c("published_symbol", "ensembl_id", "coefficient"), drop = FALSE]
}

# All signature-cohort pairs that pass the coverage rule.
# independence: data.frame(signature, cohort, used_by_authors = "yes"/"no").
# Pairs the authors used (development or validation) are skipped unless include_used = TRUE,
# which is how TCGA-LUAD enters as the apparent-performance reference.
build_pairs <- function(sig_map, measurable, independence, cohort_names, min_fraction,
                        include_used = FALSE) {
  rows <- list()
  genes <- list()
  for (cn in cohort_names) {
    for (s in sort(unique(sig_map$signature))) {
      used <- independence$used_by_authors[independence$signature == s & independence$cohort == cn]
      if (length(used) != 1 || !used %in% c("yes", "no")) {
        stop("independence.csv needs exactly one yes/no entry for ", s, " x ", cn)
      }
      if (used == "yes" && !include_used) next
      g <- pair_genes(sig_map, s, measurable[[cn]], min_fraction)
      if (is.null(g)) next
      id <- paste(s, cn, sep = "__")
      rows[[id]] <- data.frame(pair_id = id, signature = s, cohort = cn,
                               n_genes = sum(sig_map$signature == s), n_present = nrow(g),
                               used_by_authors = used, stringsAsFactors = FALSE)
      genes[[id]] <- g
    }
  }
  pairs <- if (length(rows)) do.call(rbind, rows) else NULL
  if (!is.null(pairs)) rownames(pairs) <- NULL
  list(pairs = pairs, genes = genes)
}

# Measurable genes of every cohort in a list, with the presence rule of config/analysis.yaml
# (used again on subsets in the sensitivity analyses, where a gene can become constant).
measurable_by_cohort <- function(cohorts, cfg) {
  lapply(cohorts, function(co) {
    frac <- if (identical(co$meta$platform, "RNA-seq")) cfg$presence$rnaseq_min_expressed_fraction else 0
    measurable_rows(co$expr, min_sd = cfg$presence$min_sd, min_expressed_frac = frac)
  })
}

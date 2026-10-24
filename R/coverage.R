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

# Alternative Ensembl IDs (config/gene_id_alternatives.csv, one row per ID with its evidence):
# the official Ensembl successor of a signature gene's ID (Ensembl ID history). It is used in a
# cohort only where the original ID is missing from that cohort's gene model or not measurable
# there, and only if the successor is measurable.
read_gene_alternatives <- function(file = "config/gene_id_alternatives.csv") {
  if (!file.exists(file)) return(data.frame(ensembl_id = character(0), alternative_id = character(0)))
  utils::read.csv(file, stringsAsFactors = FALSE, colClasses = "character")
}

# The signature gene map as used in one cohort: original IDs, or alternatives where needed.
cohort_gene_map <- function(sig_map, measurable_ids, alternatives = read_gene_alternatives()) {
  sig_map$original_id <- sig_map$ensembl_id
  sig_map$id_used <- "original"
  if (!nrow(alternatives)) return(sig_map)
  need <- which(!is.na(sig_map$ensembl_id) & !sig_map$ensembl_id %in% measurable_ids)
  for (i in need) {
    alt <- alternatives$alternative_id[alternatives$ensembl_id == sig_map$ensembl_id[i]]
    alt <- alt[alt %in% measurable_ids]
    if (length(alt)) {
      sig_map$ensembl_id[i] <- alt[1]
      sig_map$id_used[i] <- "alternative"
    }
  }
  dup <- !is.na(sig_map$ensembl_id) & duplicated(paste(sig_map$signature, sig_map$ensembl_id))
  if (any(dup)) {
    stop("An alternative ID duplicates another gene of the same signature: ",
         paste(sig_map$signature[dup], sig_map$published_symbol[dup], collapse = ", "))
  }
  sig_map
}

# Coverage of every signature in one cohort.
# sig_map: one row per signature gene with columns signature, published_symbol,
#          coefficient, ensembl_id (NA if unmapped).
coverage_table <- function(sig_map, measurable_ids, cohort, alternatives = read_gene_alternatives()) {
  sig_map <- cohort_gene_map(sig_map, measurable_ids, alternatives)
  sig_map$present <- !is.na(sig_map$ensembl_id) & sig_map$ensembl_id %in% measurable_ids
  rows <- lapply(split(sig_map, sig_map$signature), function(d) {
    data.frame(cohort = cohort, signature = d$signature[1], n_genes = nrow(d),
               n_present = sum(d$present), frac_present = mean(d$present),
               missing = paste(d$published_symbol[!d$present], collapse = ";"),
               via_alternative_id = paste(d$published_symbol[d$present & d$id_used == "alternative"], collapse = ";"),
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
                        include_used = FALSE, alternatives = read_gene_alternatives()) {
  rows <- list()
  genes <- list()
  base_map <- sig_map
  for (cn in cohort_names) {
    sig_map <- cohort_gene_map(base_map, measurable[[cn]], alternatives)
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

# Official Ensembl ID history (rest.ensembl.org/archive/id). `arch`: data.frame with columns id,
# current (logical) and replaced_by (";"-separated possible replacements of a retired ID).
# history_walk() returns every ID on the history chains that start at `id` (`id` included) and
# the current IDs those chains end in.
history_walk <- function(id, current, replaced_by) {
  seen <- character(0)
  ends <- character(0)
  queue <- id
  while (length(queue)) {
    x <- queue[1]
    queue <- queue[-1]
    if (x %in% seen) next
    seen <- c(seen, x)
    if (isTRUE(current[x])) {
      ends <- c(ends, x)
      next
    }
    nx <- replaced_by[x]
    nx <- if (is.na(nx)) character(0) else strsplit(nx, ";", fixed = TRUE)[[1]]
    queue <- c(queue, nx[nzchar(nx)])
  }
  list(path = seen, ends = unique(ends))
}

# Older (retired) cohort IDs whose official history leads to a signature gene that is missing
# from the cohort - the successor rule of gene_id_alternatives.csv applied in the other direction,
# for a cohort quantified with an older Ensembl release. One row per (signature gene, cohort ID).
# Decision rule (`accept`): the cohort ID is measurable, its history passes through the signature
# gene's own ID, ends in exactly one current gene, and it is the only such measurable cohort ID.
predecessor_hits <- function(sig_map, cohort_ids, measurable_ids, arch,
                             alternatives = read_gene_alternatives()) {
  empty <- data.frame(signature = character(0), published_symbol = character(0),
                      ensembl_id = character(0), cohort_id = character(0), history = character(0),
                      ends_in = character(0), n_current = integer(0), through_original = logical(0),
                      measurable = logical(0), accept = logical(0), note = character(0))
  sm <- cohort_gene_map(sig_map, measurable_ids, alternatives)
  miss <- sm[!is.na(sm$original_id) & !sm$ensembl_id %in% measurable_ids, , drop = FALSE]
  current <- stats::setNames(arch$current, arch$id)
  replaced_by <- stats::setNames(arch$replaced_by, arch$id)
  old <- intersect(cohort_ids, arch$id[arch$current %in% FALSE])
  if (!nrow(miss) || !length(old)) return(empty)
  targets <- lapply(miss$original_id, function(s) {
    c(s, alternatives$alternative_id[alternatives$ensembl_id == s])
  })
  rows <- lapply(old, function(y) {
    w <- history_walk(y, current, replaced_by)
    later <- setdiff(w$path, y)
    if (!length(later)) return(NULL)
    i <- which(vapply(targets, function(tg) any(later %in% tg), logical(1)))
    if (!length(i)) return(NULL)
    data.frame(signature = miss$signature[i], published_symbol = miss$published_symbol[i],
               ensembl_id = miss$original_id[i], cohort_id = y,
               history = paste(w$path, collapse = " > "), ends_in = paste(w$ends, collapse = ";"),
               n_current = length(w$ends), through_original = miss$original_id[i] %in% later,
               measurable = y %in% measurable_ids, stringsAsFactors = FALSE)
  })
  hits <- do.call(rbind, rows)
  if (is.null(hits) || !nrow(hits)) return(empty)
  usable <- hits$measurable & hits$through_original
  n_usable <- stats::ave(as.integer(usable), hits$signature, hits$ensembl_id, FUN = sum)
  in_sig <- paste(hits$signature, hits$cohort_id) %in% paste(sig_map$signature, sig_map$ensembl_id)
  hits$accept <- usable & hits$n_current == 1 & n_usable == 1 & !in_sig
  hits$note <- ifelse(hits$accept, "accepted: official predecessor",
               ifelse(!hits$measurable, "not measurable in the cohort",
               ifelse(!hits$through_original, "history reaches only the listed alternative ID",
               ifelse(hits$n_current > 1, "split: history ends in several current genes",
               ifelse(in_sig, "ID is already another gene of this signature",
                      "merge: several cohort IDs lead to this gene")))))
  hits <- hits[order(hits$signature, hits$published_symbol, hits$cohort_id), ]
  rownames(hits) <- NULL
  hits
}

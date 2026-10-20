# ids.R - gene identifiers. Everything is matched by Ensembl gene ID without version,
# never by gene name.

strip_version <- function(ids) sub("\\.[0-9]+$", "", ids)

# Gene records of a GENCODE GTF (any release): gene_id / gene_name / gene_type.
read_gencode_genes <- function(gtf_path) {
  gr <- rtracklayer::import(gtf_path, feature.type = "gene")
  m <- as.data.frame(S4Vectors::mcols(gr))
  keep <- !grepl("_PAR_Y$", m$gene_id)
  data.frame(gene_id = strip_version(m$gene_id[keep]), gene_name = m$gene_name[keep],
             gene_type = m$gene_type[keep], stringsAsFactors = FALSE)
}

# Gene table from one GDC STAR-Counts file (GENCODE v36): no GTF needed.
read_star_gene_table <- function(star_tsv) {
  d <- utils::read.delim(star_tsv, comment.char = "#", stringsAsFactors = FALSE)
  d <- d[grepl("^ENSG", d$gene_id) & !grepl("_PAR_Y$", d$gene_id),
         c("gene_id", "gene_name", "gene_type")]
  d$gene_id <- strip_version(d$gene_id)
  d
}

# Look up each published gene name in several gene tables (for example GENCODE v22 and v36).
# Returns one row per name, the Ensembl ID(s) found in each table, and a status:
#   consistent = one ID across all tables that know the name
#   conflict   = different IDs (or two IDs in one table): decide by hand
#   not_found  = no table knows the name: decide by hand
resolve_symbols <- function(symbols, gene_tables) {
  stopifnot(is.list(gene_tables), !is.null(names(gene_tables)))
  symbols <- unique(symbols)
  out <- data.frame(published_symbol = symbols, stringsAsFactors = FALSE)
  for (v in names(gene_tables)) {
    g <- gene_tables[[v]]
    out[[paste0("ids_", v)]] <- unname(vapply(symbols, function(s) {
      ids <- sort(unique(g$gene_id[g$gene_name == s]))
      if (length(ids)) paste(ids, collapse = ";") else NA_character_
    }, character(1)))
  }
  id_cols <- grep("^ids_", names(out), value = TRUE)
  out$status <- apply(out[, id_cols, drop = FALSE], 1, function(r) {
    r <- r[!is.na(r)]
    if (!length(r)) return("not_found")
    all_ids <- unique(unlist(strsplit(r, ";")))
    if (length(all_ids) == 1) "consistent" else "conflict"
  })
  out
}

# Final Ensembl ID per signature gene.
#   1. a row in config/gene_overrides.csv wins (record the evidence there);
#   2. otherwise the single ID found in the annotation the paper used (column `annotation`,
#      which must name one of the gene tables, e.g. "v22" or "v36");
#   3. otherwise, if all tables agree, that ID.
# Genes left without an ID are flagged needs_review = TRUE.
choose_ids <- function(sig_genes, resolved, overrides = NULL) {
  m <- merge(sig_genes, resolved, by = "published_symbol", all.x = TRUE, sort = FALSE)
  m$ensembl_id <- NA_character_
  m$id_source <- NA_character_
  id_cols <- grep("^ids_", names(m), value = TRUE)
  for (i in seq_len(nrow(m))) {
    col <- paste0("ids_", m$annotation[i])
    ids <- if (col %in% names(m)) m[[col]][i] else NA_character_
    if (!is.na(ids) && !grepl(";", ids)) {
      m$ensembl_id[i] <- ids
      m$id_source[i] <- m$annotation[i]
    } else if (identical(m$status[i], "consistent")) {
      any_id <- unique(stats::na.omit(unlist(m[i, id_cols])))
      m$ensembl_id[i] <- any_id[1]
      m$id_source[i] <- "other_version"
    }
  }
  if (!is.null(overrides) && nrow(overrides) > 0) {
    hit <- match(paste(m$signature, m$published_symbol),
                 paste(overrides$signature, overrides$published_symbol))
    m$ensembl_id[!is.na(hit)] <- overrides$ensembl_id[hit[!is.na(hit)]]
    m$id_source[!is.na(hit)] <- "manual_override"
  }
  m$needs_review <- is.na(m$ensembl_id)
  m[order(m$signature, m$published_symbol), , drop = FALSE]
}

# Collapse a probe-level array matrix to Ensembl genes with rules fixed in advance:
#   - drop probes that map to more than one Ensembl gene;
#   - per gene keep one probe: "_at" before "_s_at" before "_x_at", then the highest value
#     in `probe_means` (pass the same probe_means to every cohort of a platform).
collapse_probes <- function(expr, probe_map, probe_means = rowMeans(expr)) {
  pm <- unique(probe_map[probe_map$probe_id %in% rownames(expr) &
                           !is.na(probe_map$ensembl_id) & probe_map$ensembl_id != "",
                         c("probe_id", "ensembl_id")])
  multi <- names(which(table(pm$probe_id) > 1))
  pm <- pm[!pm$probe_id %in% multi, , drop = FALSE]
  pm$suffix_rank <- ifelse(grepl("_x_at$", pm$probe_id), 3L,
                           ifelse(grepl("_s_at$", pm$probe_id), 2L, 1L))
  pm$mean_expr <- unname(probe_means[pm$probe_id])
  pm <- pm[order(pm$ensembl_id, pm$suffix_rank, -pm$mean_expr), , drop = FALSE]
  chosen <- pm[!duplicated(pm$ensembl_id), , drop = FALSE]
  out <- expr[chosen$probe_id, , drop = FALSE]
  rownames(out) <- chosen$ensembl_id
  list(expr = out, probe_choice = chosen)
}

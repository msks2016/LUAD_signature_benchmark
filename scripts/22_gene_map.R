# 22_gene_map.R - Step 3b: one Ensembl ID for every published signature gene.
# Stops until every gene has an ID; unresolved genes go into config/gene_overrides.csv
# with the evidence (HGNC previous symbols, Ensembl ID history, coordinates, the paper).

for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) source(f)

sig <- utils::read.csv("data/signatures/signature_master.csv", stringsAsFactors = FALSE)
tables <- readRDS("data/reference/gene_tables.rds")
bad_annot <- setdiff(unique(sig$annotation), names(tables))
if (length(bad_annot)) stop("`annotation` must name a gene table (", paste(names(tables), collapse = ", "),
                            "); found: ", paste(bad_annot, collapse = ", "))

resolved <- resolve_symbols(sig$published_symbol, tables)
overrides <- utils::read.csv("config/gene_overrides.csv", stringsAsFactors = FALSE)
gene_map <- choose_ids(sig[, c("signature", "published_symbol", "coefficient", "annotation")],
                       resolved, overrides)
utils::write.csv(gene_map, "data/signatures/gene_map.csv", row.names = FALSE)

# Review tables. (1) For each signature, how many of its gene names each gene table knows: the
# `annotation` column should name the table whose naming the paper follows (the one that knows
# the most names). (2) Names that point to different Ensembl IDs in different tables
# ("conflict"): the paper's annotation decides; every such gene is listed for review.
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
id_cols <- grep("^ids_", names(gene_map), value = TRUE)
chk <- do.call(rbind, lapply(split(gene_map, gene_map$signature), function(d) {
  found <- vapply(id_cols, function(cc) sum(!is.na(d[[cc]])), 1L)
  names(found) <- sub("^ids_", "names_in_", id_cols)
  data.frame(signature = d$signature[1], annotation = d$annotation[1], n_genes = nrow(d), t(found),
             best_table = paste(sub("^ids_", "", id_cols)[found == max(found)], collapse = "/"),
             conflicts = sum(d$status %in% "conflict"), not_found = sum(d$status %in% "not_found"),
             stringsAsFactors = FALSE, check.names = FALSE)
}))
chk$annotation_knows_fewer_names <- !mapply(function(a, b) a %in% strsplit(b, "/")[[1]], chk$annotation, chk$best_table)
utils::write.csv(chk, "results/tables/gene_map_annotation_check.csv", row.names = FALSE)
conf <- gene_map[gene_map$status %in% "conflict",
                 c("signature", "published_symbol", "annotation", id_cols, "ensembl_id", "id_source")]
for (v in sub("^ids_", "", id_cols)) {
  conf[[paste0("type_", v)]] <- tables[[v]]$gene_type[match(conf$ensembl_id, tables[[v]]$gene_id)]
}
utils::write.csv(conf, "results/tables/gene_map_conflicts.csv", row.names = FALSE)
print(chk)
if (nrow(conf)) {
  message("\n", nrow(conf), " gene name(s) point to different Ensembl IDs in different tables ",
          "(listed in results/tables/gene_map_conflicts.csv): send that file for review.")
}

if (any(gene_map$needs_review)) {
  print(gene_map[gene_map$needs_review, c("signature", "published_symbol", "status")])
  stop("Add the genes above to config/gene_overrides.csv with evidence, then re-run.")
}
dup <- duplicated(paste(gene_map$signature, gene_map$ensembl_id))
if (any(dup)) {
  print(gene_map[dup, c("signature", "published_symbol", "ensembl_id")])
  stop("Two published names map to the same Ensembl ID within one signature.")
}

record_number("genes_resolved_directly", sum(gene_map$id_source %in% names(tables)))
record_number("genes_resolved_other_version", sum(gene_map$id_source == "other_version"))
record_number("genes_manual_override", sum(gene_map$id_source == "manual_override"))
record_number("genes_version_conflict", sum(gene_map$status == "conflict", na.rm = TRUE))
message("Gene map written: ", nrow(gene_map), " genes")

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

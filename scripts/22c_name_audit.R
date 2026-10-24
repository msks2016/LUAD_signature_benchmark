# 22c_name_audit.R - audit of every signature gene name across the GENCODE releases available
# locally: v22 and v36 (gene tables) plus any other data/reference/gencode.vNN.annotation.gtf.gz
# (script 22b downloaded v29, v31 and v32). For each name: the Ensembl ID(s) and gene type in each
# release, the ID chosen by script 22, and two flags:
#   type_flag   - the chosen gene is not a long non-coding RNA type in GENCODE v36;
#   id_unstable - the name points to different Ensembl IDs in different releases.
# Flagged genes are reviewed by hand (decisions go into config/gene_overrides.csv with evidence).
# Output: results/tables/gene_name_audit.csv

for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) source(f)
gm <- utils::read.csv("data/signatures/gene_map.csv", stringsAsFactors = FALSE)
tables <- readRDS("data/reference/gene_tables.rds")
gtfs <- list.files("data/reference", pattern = "^gencode\\.v[0-9]+\\.annotation\\.gtf\\.gz$", full.names = TRUE)
versions <- as.integer(sub("^gencode\\.v([0-9]+)\\..*$", "\\1", basename(gtfs)))
for (k in order(versions)) {
  v <- paste0("v", versions[k])
  if (v %in% names(tables)) next
  message("Reading GENCODE ", v, " ...")
  gr <- rtracklayer::import(gtfs[k], feature.type = "gene")
  m <- as.data.frame(S4Vectors::mcols(gr))
  tables[[v]] <- data.frame(gene_id = strip_version(m$gene_id), gene_name = m$gene_name,
                            gene_type = m$gene_type, stringsAsFactors = FALSE)
}
tables <- tables[order(as.integer(sub("^v", "", names(tables))))]
lnc_types <- c("lncRNA", "lincRNA", "antisense", "processed_transcript", "sense_intronic",
               "sense_overlapping", "3prime_overlapping_ncRNA", "bidirectional_promoter_lncRNA",
               "macro_lncRNA", "non_coding")

names_q <- unique(gm$published_symbol)
audit <- data.frame(published_symbol = names_q, stringsAsFactors = FALSE)
for (v in names(tables)) {
  g <- tables[[v]]
  audit[[paste0("ids_", v)]] <- vapply(names_q, function(s) {
    x <- g[g$gene_name == s, , drop = FALSE]
    if (!nrow(x)) NA_character_ else paste(paste0(x$gene_id, " (", x$gene_type, ")"), collapse = "; ")
  }, character(1))
}
chosen <- unique(gm[, c("published_symbol", "ensembl_id", "id_source")])
chosen <- stats::aggregate(cbind(ensembl_id, id_source) ~ published_symbol, data = chosen,
                           FUN = function(x) paste(unique(x), collapse = " | "))
audit <- merge(audit, chosen, by = "published_symbol", all.x = TRUE)
audit$signatures <- vapply(audit$published_symbol, function(s) paste(gm$signature[gm$published_symbol == s], collapse = ";"), "")
first_id <- sub(" \\|.*$", "", audit$ensembl_id)
audit$type_chosen_v36 <- tables$v36$gene_type[match(first_id, tables$v36$gene_id)]
audit$type_chosen_v22 <- tables$v22$gene_type[match(first_id, tables$v22$gene_id)]
id_cols <- grep("^ids_", names(audit), value = TRUE)
audit$id_unstable <- apply(audit[, id_cols, drop = FALSE], 1, function(r) {
  ids <- unique(unlist(regmatches(r[!is.na(r)], gregexpr("ENSG[0-9]+", r[!is.na(r)]))))
  length(ids) > 1
})
audit$type_flag <- !is.na(audit$type_chosen_v36) & !audit$type_chosen_v36 %in% lnc_types
audit <- audit[order(!audit$type_flag, !audit$id_unstable, audit$published_symbol), ]
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
utils::write.csv(audit, "results/tables/gene_name_audit.csv", row.names = FALSE)
message("Releases compared: ", paste(names(tables), collapse = ", "))
print(audit[audit$type_flag | audit$id_unstable, c("published_symbol", "signatures", "ensembl_id",
                                                     "type_chosen_v36", "type_flag", "id_unstable")],
      row.names = FALSE)

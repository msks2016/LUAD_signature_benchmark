# 22b_name_history.R - evidence for gene names that script 22 could not place (in neither the
# GENCODE v22 nor the v36 gene table). Two independent lookups:
#  1. NCBI Gene via org.Hs.eg.db (installed with the GPL570 packages): the name as a current
#     symbol or a previous symbol/alias -> Entrez ID -> Ensembl ID;
#  2. intermediate GENCODE releases (v33, v32, v31, v29, v27, v35, v38): the GTF is downloaded
#     once into data/reference/ (about 45 MB each; kept out of git) until every name is found.
# For every Ensembl ID found it shows the gene's name and type in v22 and v36. The result is the
# evidence for config/gene_overrides.csv.
# Output: results/tables/gene_name_history.csv

for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) source(f)
gm <- utils::read.csv("data/signatures/gene_map.csv", stringsAsFactors = FALSE)
tables <- readRDS("data/reference/gene_tables.rds")
names_q <- unique(gm$published_symbol[gm$needs_review])
if (!length(names_q)) stop("No gene needs review in data/signatures/gene_map.csv")
message("Looking up: ", paste(names_q, collapse = ", "))
hits <- list()

# 1. NCBI Gene (org.Hs.eg.db)
for (key in c("SYMBOL", "ALIAS")) {
  for (nm in names_q) {
    r <- tryCatch(suppressMessages(AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, keys = nm,
                                                         keytype = key, columns = c("ENTREZID", "SYMBOL", "GENETYPE", "ENSEMBL"))),
                  error = function(e) NULL)
    if (!is.null(r) && nrow(r)) {
      hits[[paste("ncbi", key, nm)]] <- data.frame(query = nm, source = paste0("NCBI Gene (", tolower(key), ")"),
                                                   ensembl_id = r$ENSEMBL, name_in_source = r$SYMBOL,
                                                   type_in_source = r$GENETYPE, entrez_id = r$ENTREZID,
                                                   location = NA_character_, stringsAsFactors = FALSE)
    }
  }
}

# 2. Intermediate GENCODE releases
gtf_genes <- function(v) {
  f <- sprintf("data/reference/gencode.v%d.annotation.gtf.gz", v)
  if (!file.exists(f)) {
    url <- sprintf("https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_%d/gencode.v%d.annotation.gtf.gz", v, v)
    message("  downloading GENCODE v", v, " ...")
    old <- options(timeout = 3600)
    ok <- tryCatch(utils::download.file(url, paste0(f, ".part"), mode = "wb", quiet = TRUE) == 0,
                   error = function(e) FALSE)
    options(old)
    if (!ok) { unlink(paste0(f, ".part")); message("  GENCODE v", v, " download failed"); return(NULL) }
    file.rename(paste0(f, ".part"), f)
  }
  gr <- rtracklayer::import(f, feature.type = "gene")
  d <- as.data.frame(gr)
  data.frame(gene_id = strip_version(d$gene_id), gene_name = d$gene_name, gene_type = d$gene_type,
             location = paste0(d$seqnames, ":", d$start, "-", d$end, "(", d$strand, ")"),
             stringsAsFactors = FALSE)
}
left <- names_q
for (v in c(33, 32, 31, 29, 27, 35, 38)) {
  if (!length(left)) break
  g <- gtf_genes(v)
  if (is.null(g)) next
  m <- g[g$gene_name %in% left, , drop = FALSE]
  if (nrow(m)) {
    hits[[paste("gencode", v)]] <- data.frame(query = m$gene_name, source = paste0("GENCODE v", v),
                                              ensembl_id = m$gene_id, name_in_source = m$gene_name,
                                              type_in_source = m$gene_type, entrez_id = NA_character_,
                                              location = m$location, stringsAsFactors = FALSE)
    left <- setdiff(left, m$gene_name)
  }
}

res <- if (length(hits)) do.call(rbind, hits) else
  data.frame(query = names_q, source = "not found", ensembl_id = NA, name_in_source = NA,
             type_in_source = NA, entrez_id = NA, location = NA)
res <- unique(res)
res$name_in_v22 <- tables$v22$gene_name[match(res$ensembl_id, tables$v22$gene_id)]
res$name_in_v36 <- tables$v36$gene_name[match(res$ensembl_id, tables$v36$gene_id)]
res$type_in_v36 <- tables$v36$gene_type[match(res$ensembl_id, tables$v36$gene_id)]
res <- res[order(res$query, res$source), ]
rownames(res) <- NULL
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
utils::write.csv(res, "results/tables/gene_name_history.csv", row.names = FALSE)
print(res)
if (length(left)) message("Not found in any GENCODE release tried: ", paste(left, collapse = ", "))

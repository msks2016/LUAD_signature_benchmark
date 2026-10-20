# 23b_presence_check.R - why is a signature gene missing in a cohort? For every signature gene x
# cohort: is the Ensembl ID in the cohort's matrix, is it measurable, how often is it expressed,
# does the ID still exist in Ensembl 114 (the release of the microarray probe maps), and, if
# not, which ID replaced it (Ensembl REST ID history). Needs internet for the last two.
# Run after script 23. Output: results/tables/presence_check.csv (one row per gene x cohort)

for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) source(f)
cfg <- load_config()
gm <- utils::read.csv("data/signatures/gene_map.csv", stringsAsFactors = FALSE)
meas <- readRDS("data/processed/measurable_ids.rds")
tables <- readRDS("data/reference/gene_tables.rds")
cohort_names <- c(cfg$cohorts$reference, cfg$cohorts$independent)
cohorts <- load_cohorts(cohort_names)
ids <- unique(gm$ensembl_id)

# 1. Which IDs exist in Ensembl 114
e114 <- tryCatch({
  mart <- biomaRt::useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", version = 114)
  biomaRt::getBM(attributes = c("ensembl_gene_id", "external_gene_name", "gene_biotype"),
                 filters = "ensembl_gene_id", values = ids, mart = mart)
}, error = function(e) { message("Ensembl 114 BioMart not reachable: ", conditionMessage(e)); NULL })

# 2. Ensembl REST ID history for IDs missing from Ensembl 114
history <- function(id) {
  r <- tryCatch(httr::GET(paste0("https://rest.ensembl.org/archive/id/", id), httr::accept_json(),
                          httr::timeout(60)), error = function(e) NULL)
  if (is.null(r) || httr::status_code(r) != 200) return(c(id_current = NA, id_last_release = NA, replaced_by = NA))
  j <- jsonlite::fromJSON(httr::content(r, as = "text", encoding = "UTF-8"))
  rep <- if (length(j$possible_replacement) && NROW(j$possible_replacement)) {
    paste(unique(unlist(j$possible_replacement$stable_id)), collapse = ";")
  } else ""
  c(id_current = as.character(j$is_current), id_last_release = as.character(j$release), replaced_by = rep)
}
in_e114 <- if (is.null(e114)) rep(NA, length(ids)) else ids %in% e114$ensembl_gene_id
names(in_e114) <- ids
hist_tab <- do.call(rbind, lapply(ids[in_e114 %in% FALSE], function(id) data.frame(ensembl_id = id, t(history(id)))))

# 3. One row per signature gene x cohort
rows <- lapply(cohort_names, function(cn) {
  co <- cohorts[[cn]]
  rnaseq <- identical(co$meta$platform, "RNA-seq")
  data.frame(cohort = cn, platform = co$meta$platform,
             gm[, c("signature", "published_symbol", "ensembl_id", "id_source")],
             in_matrix = gm$ensembl_id %in% rownames(co$expr),
             measurable = gm$ensembl_id %in% meas[[cn]],
             frac_expressed = vapply(gm$ensembl_id, function(id) {
               if (!id %in% rownames(co$expr) || !rnaseq) NA_real_ else round(mean(co$expr[id, ] > 0), 3)
             }, numeric(1)),
             mean_log2 = vapply(gm$ensembl_id, function(id) {
               if (!id %in% rownames(co$expr)) NA_real_ else round(mean(co$expr[id, ]), 3)
             }, numeric(1)),
             stringsAsFactors = FALSE)
})
pc <- do.call(rbind, rows)
pc$name_v36 <- tables$v36$gene_name[match(pc$ensembl_id, tables$v36$gene_id)]
pc$in_v36 <- pc$ensembl_id %in% tables$v36$gene_id
pc$in_ensembl114 <- unname(in_e114[pc$ensembl_id])
if (!is.null(hist_tab) && nrow(hist_tab)) pc <- merge(pc, hist_tab, by = "ensembl_id", all.x = TRUE)
pc$reason <- with(pc, ifelse(measurable, "present",
  ifelse(in_matrix, "in the matrix but below the presence rule",
  ifelse(platform == "RNA-seq", "ID not in this cohort's gene model",
  ifelse(in_ensembl114 %in% FALSE, "ID not in Ensembl 114 (retired or replaced): no probe can map to it",
         "no probe on this array maps to this gene")))))
pc <- pc[order(pc$cohort, pc$signature, pc$published_symbol), ]
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
utils::write.csv(pc, "results/tables/presence_check.csv", row.names = FALSE)
print(table(pc$reason, pc$cohort))
miss <- unique(pc[pc$reason != "present" & !pc$cohort %in% "GSE72094",
                  c("cohort", "published_symbol", "ensembl_id", "reason", "frac_expressed")])
print(miss, row.names = FALSE)

# 23c_predecessor_check.R - older Ensembl IDs in GSE81089.
# GSE81089 was quantified with Ensembl 73 (2013); the signature genes carry GENCODE v22-v36 IDs.
# A signature gene missing from GSE81089 can be there under an older ID that Ensembl later
# retired and replaced (merge or re-annotation). This script follows the official Ensembl ID
# history (rest.ensembl.org/archive/id) of every retired ID in the GSE81089 matrix and lists the
# ones whose history leads to a missing signature gene: the successor rule of
# config/gene_id_alternatives.csv applied in the other direction.
# Decision rule (fixed 20 Oct 2026, before any outcome analysis; see R/coverage.R,
# predecessor_hits): an older ID is accepted when it is measurable in GSE81089, its history passes
# through the signature gene's own ID and ends in exactly one current gene, and it is the only
# such GSE81089 ID. Accepted IDs are added to config/gene_id_alternatives.csv with the evidence.
# Needs internet (rest.ensembl.org); about 70 requests, 5-15 minutes; resumes if interrupted.
# Run after script 23, then run script 23 again.
# Output: results/tables/gse81089_predecessor_check.csv
#         data/reference/gse81089_ensembl_archive.rds (cache of the ID history)

for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) source(f)
cohort <- "GSE81089"
gm <- utils::read.csv("data/signatures/gene_map.csv", stringsAsFactors = FALSE)
alts <- read_gene_alternatives()
meas <- readRDS("data/processed/measurable_ids.rds")[[cohort]]
co <- load_cohorts(cohort)[[cohort]]
ids <- rownames(co$expr)

cache_file <- "data/reference/gse81089_ensembl_archive.rds"
arch <- if (file.exists(cache_file)) readRDS(cache_file) else
  data.frame(id = character(0), current = logical(0), release = character(0),
             replaced_by = character(0), stringsAsFactors = FALSE)

archive_post <- function(x) {
  for (attempt in 1:5) {
    r <- tryCatch(httr::POST("https://rest.ensembl.org/archive/id", httr::accept_json(),
                             httr::content_type_json(), httr::timeout(600),
                             body = as.character(jsonlite::toJSON(list(id = x)))),
                  error = function(e) NULL)
    if (!is.null(r) && httr::status_code(r) == 200) {
      j <- jsonlite::fromJSON(httr::content(r, as = "text", encoding = "UTF-8"), simplifyVector = FALSE)
      out <- lapply(j, function(e) {
        rep <- vapply(e$possible_replacement, function(p) as.character(p$stable_id), "")
        data.frame(id = as.character(e$id),
                   current = as.character(if (is.null(e$is_current)) "" else e$is_current) %in% c("1", "TRUE", "true"),
                   release = as.character(if (is.null(e$release)) NA else e$release),
                   replaced_by = paste(unique(rep), collapse = ";"), stringsAsFactors = FALSE)
      })
      return(do.call(rbind, out))
    }
    message("  retry ", attempt, " (", if (is.null(r)) "no connection" else httr::status_code(r), ")")
    Sys.sleep(15 * attempt)
  }
  stop("rest.ensembl.org/archive/id not reachable: run this script again later (finished batches are kept)")
}

lookup <- function(x) {
  x <- setdiff(unique(x), arch$id)
  if (!length(x)) return(invisible(NULL))
  batches <- split(x, ceiling(seq_along(x) / 1000))
  for (i in seq_along(batches)) {
    message("Ensembl ID history: batch ", i, " of ", length(batches))
    res <- archive_post(batches[[i]])
    unknown <- setdiff(batches[[i]], res$id)
    if (length(unknown)) {
      res <- rbind(res, data.frame(id = unknown, current = NA, release = NA_character_,
                                   replaced_by = "", stringsAsFactors = FALSE))
    }
    arch <<- rbind(arch, res)
    saveRDS(arch, cache_file)
  }
}

lookup(ids)
for (round in 1:10) {                        # follow replacements that are themselves retired
  nxt <- unlist(strsplit(arch$replaced_by[arch$current %in% FALSE], ";", fixed = TRUE))
  nxt <- setdiff(nxt[nzchar(nxt)], arch$id)
  if (!length(nxt)) break
  lookup(nxt)
}
message(sum(arch$id %in% ids & arch$current %in% FALSE), " of ", length(ids),
        " GSE81089 IDs are retired in the current Ensembl release")

hits <- predecessor_hits(gm, ids, meas, arch, alts)
hits$cohort_id_last_release <- arch$release[match(hits$cohort_id, arch$id)]
hits$frac_expressed <- vapply(hits$cohort_id, function(y) round(mean(co$expr[y, ] > 0), 3), numeric(1))
hits$cohort_id_used_by <- vapply(hits$cohort_id, function(y) {
  paste(unique(gm$signature[gm$ensembl_id %in% y]), collapse = ";")
}, "")
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
utils::write.csv(hits, "results/tables/gse81089_predecessor_check.csv", row.names = FALSE)

acc <- hits[hits$accept, , drop = FALSE]
if (nrow(acc)) {
  acc <- do.call(rbind, lapply(split(acc, paste(acc$ensembl_id, acc$cohort_id)), function(d) {
    data.frame(ensembl_id = d$ensembl_id[1], alternative_id = d$cohort_id[1],
               published_symbol = paste(unique(d$published_symbol), collapse = "/"),
               evidence = sprintf(paste0(
                 "Ensembl REST ID history (rest.ensembl.org/archive/id, queried %s by ",
                 "scripts/23c_predecessor_check.R): %s, quantified in GSE81089 (Ensembl 73), was ",
                 "last in release %s; its history (%s) leads only to this gene. Official predecessor; ",
                 "used only where the original ID is not measurable."),
                 format(Sys.Date(), "%d %b %Y"), d$cohort_id[1], d$cohort_id_last_release[1], d$history[1]),
               stringsAsFactors = FALSE)
  }))
  old_alts <- read_gene_alternatives()
  new <- acc[!paste(acc$ensembl_id, acc$alternative_id) %in%
               paste(old_alts$ensembl_id, old_alts$alternative_id), , drop = FALSE]
  if (nrow(new)) {
    utils::write.csv(rbind(old_alts, new[, names(old_alts)]), "config/gene_id_alternatives.csv",
                     row.names = FALSE)
  }
  message(nrow(new), " predecessor ID(s) added to config/gene_id_alternatives.csv: run scripts/23_coverage.R again")
} else {
  message("No predecessor IDs accepted: config/gene_id_alternatives.csv unchanged")
}
record_number("gse81089_predecessor_ids_accepted", sum(hits$accept))
print(hits[, c("signature", "published_symbol", "cohort_id", "measurable", "accept", "note")], row.names = FALSE)

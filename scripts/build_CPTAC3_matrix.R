# scripts/build_CPTAC3_matrix.R
# Reads nested STAR-Counts TSVs from gdc-client download.
# Filters to primary tumors via manifest. Normalization: log2(TPM + 1).

manifest <- read.csv("data/raw/cptac3_luad_manifest.csv", stringsAsFactors = FALSE)
manifest <- manifest[manifest$sample_type == "Primary Tumor", ]

all_tsv <- list.files("data/raw/CPTAC3_LUAD",
                      pattern = "star_gene_counts\\.tsv$",
                      recursive = TRUE, full.names = TRUE)

manifest$path <- sapply(manifest$file_id, function(id) {
  hit <- all_tsv[grepl(id, all_tsv, fixed = TRUE)]
  if (length(hit) > 0) hit[1] else NA_character_
})
manifest <- manifest[!is.na(manifest$path), ]

cat("Primary tumor files on disk:", nrow(manifest), "\n")
stopifnot(nrow(manifest) >= 150)

read_star <- function(path) {
  d <- read.delim(path, comment.char = "#", stringsAsFactors = FALSE)
  d[, c("gene_id", "gene_name", "gene_type", "tpm_unstranded")]
}

first  <- read_star(manifest$path[1])
counts <- matrix(0, nrow = nrow(first), ncol = nrow(manifest))
rownames(counts) <- first$gene_name
colnames(counts) <- substr(manifest$cases_submitter_id, 1, 12)

for (i in seq_len(nrow(manifest))) {
  d <- read_star(manifest$path[i])
  counts[, i] <- as.numeric(d$tpm_unstranded)
  if (i %% 25 == 0) cat(i, "/", nrow(manifest), "\n")
}

keep   <- first$gene_type %in% c("protein_coding", "lncRNA")
counts <- counts[keep, ]
cat("Matrix:", nrow(counts), "genes x", ncol(counts), "samples\n")

expr <- log2(counts + 1)
saveRDS(expr, "data/processed/CPTAC3_LUAD_expr.rds")
cat("Saved:", nrow(expr), "genes x", ncol(expr), "patients\n")
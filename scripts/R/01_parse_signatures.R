# ============================================================
# 01_parse_signatures.R - Parse all signature files
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
})

sig_files <- list.files("signatures", pattern = "\\.tsv$", full.names = TRUE)

cat("========================================\n")
cat("Found", length(sig_files), "signature files\n")
cat("========================================\n\n")

parse_signature <- function(file) {
  fname <- basename(file)
  paper_id <- stringr::str_extract(fname, "Paper\\d+")
  author_year <- fname %>%
    stringr::str_remove("\\.tsv$") %>%
    stringr::str_remove("Paper\\d+_")
  df <- fread(file, sep = "\t", header = TRUE)
  df$paper_id <- paper_id
  df$author_year <- author_year
  return(df)
}

all_sigs <- purrr::map_dfr(sig_files, parse_signature)

sig_summary <- all_sigs %>%
  group_by(paper_id, author_year) %>%
  summarise(
    n_genes = n(),
    n_risk = sum(coefficient > 0),
    n_protective = sum(coefficient < 0),
    .groups = "drop"
  ) %>%
  arrange(paper_id)

cat("=== Signature Summary ===\n")
print(sig_summary, n = 20)

dir.create("data/processed", showWarnings = FALSE, recursive = TRUE)
dir.create("results/tables", showWarnings = FALSE, recursive = TRUE)

write_csv(all_sigs, "data/processed/all_signatures_combined.csv")
write_csv(sig_summary, "results/tables/signature_summary.csv")

gene_freq <- all_sigs %>%
  count(gene, sort = TRUE) %>%
  filter(n > 1)

cat("\n=== Genes in Multiple Signatures ===\n")
print(as_tibble(gene_freq), n = 30)

write_csv(gene_freq, "results/tables/duplicate_genes.csv")

cat("\n========================================\n")
cat("PARSING COMPLETE\n")
cat("========================================\n")
cat("Total signatures:", n_distinct(all_sigs$paper_id), "\n")
cat("Total unique genes:", n_distinct(all_sigs$gene), "\n")
cat("Total gene-signature pairs:", nrow(all_sigs), "\n")

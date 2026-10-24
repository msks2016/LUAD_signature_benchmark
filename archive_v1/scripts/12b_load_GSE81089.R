# ============================================================
# 12b_load_GSE81089.R (FIXED)
# ============================================================
suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
})

# --- Load ---
counts <- fread("data/raw/geo/GSE81089_readcounts.tsv", sep = "\t")
cat("=== GSE81089 ===\n")
cat("Dimensions:", nrow(counts), "genes x", ncol(counts), "columns\n\n")

# --- Gene IDs ---
gene_ids <- counts[[1]]
cat("Gene IDs are:", ifelse(grepl("^ENSG", gene_ids[1]), "Ensembl", "Symbols"), "\n\n")

# --- Signature coverage ---
sigs <- read.csv("data/processed/signature_coefficients.csv")
sig_genes <- unique(sigs$gene)

gmap <- read.csv("docs/gene_mapping.csv")
# Build named vector: paper_name -> paper_ensembl
paper_to_ens <- setNames(as.character(gmap$paper_ensembl), as.character(gmap$paper_name))
cat("Gene mapping table loaded:", length(paper_to_ens), "entries\n\n")

# --- FIXED check_gene function ---
check_gene <- function(g) {
  # Direct symbol match
  if (g %in% gene_ids) return(TRUE)
  
  # Look up paper Ensembl ID safely
  if (g %in% names(paper_to_ens)) {
    eid <- paper_to_ens[g]  # single bracket (safe)
    if (!is.na(eid) && nchar(eid) > 0) {
      return(sub("\\..*$", "", eid) %in% gene_ids)
    }
  }
  return(FALSE)
}

# --- Test on each signature gene ---
cat("=== Checking each signature gene ===\n")
match_list <- sapply(sig_genes, check_gene)
names(match_list) <- sig_genes

present <- sig_genes[match_list]
missing <- sig_genes[!match_list]

cat("\nPresent:", length(present), "/", length(sig_genes), 
    "(", round(100 * length(present) / length(sig_genes), 1), "%)\n")
cat("\nMissing genes:\n")
print(missing)

# --- Per-signature coverage ---
cov <- sigs %>%
  filter(!duplicated(paste(record_no, gene))) %>%
  group_by(record_no) %>%
  summarise(
    n_genes = n(),
    n_present = sum(sapply(gene, check_gene)),
    pct = round(100 * n_present / n_genes, 1),
    .groups = "drop"
  ) %>%
  arrange(desc(pct))

cat("\n=== Per-signature coverage ===\n")
print(cov, n = 20)

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
write.csv(cov, "results/tables/GSE81089_coverage.csv", row.names = FALSE)

# --- Save ---
saveRDS(counts, "data/raw/geo/GSE81089_counts.rds")
cat("\nSaved: data/raw/geo/GSE81089_counts.rds\n")

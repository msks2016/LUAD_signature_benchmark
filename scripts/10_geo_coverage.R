# ============================================================
# 10_geo_coverage.R
# GEO probe coverage check — how many signature genes are
# represented on Affymetrix U133 Plus 2.0 (GPL570)?
# ============================================================
suppressPackageStartupMessages({
  library(tidyverse)
  library(GEOquery)
})

dir.create("data/raw/geo", recursive = TRUE, showWarnings = FALSE)

message("Loading GPL570 annotation (1-2 min)...")
gpl <- getGEO("GPL570", destdir = "data/raw/geo")
gpl_df <- Table(gpl)

cat("GPL570 has", nrow(gpl_df), "probes\n")
cat("Columns:", paste(head(colnames(gpl_df), 8), collapse = ", "), "\n\n")

# --- Map probes to gene symbols ---
probe2gene <- gpl_df %>%
  select(ID, `Gene Symbol`) %>%
  rename(probe_id = ID, gene_symbol = `Gene Symbol`) %>%
  separate_rows(gene_symbol, sep = " /// ") %>%
  filter(gene_symbol != "", !is.na(gene_symbol)) %>%
  distinct(gene_symbol) %>%
  pull(gene_symbol)

cat("Unique gene symbols on GPL570:", length(probe2gene), "\n\n")

# --- Coverage per signature ---
sigs <- read.csv("data/processed/signature_coefficients.csv")

coverage <- sigs %>%
  filter(!duplicated(paste(record_no, gene))) %>%
  group_by(record_no) %>%
  summarise(
    n_genes = n(),
    present = sum(gene %in% probe2gene),
    coverage_pct = round(100 * sum(gene %in% probe2gene) / n(), 1),
    missing = paste(gene[!gene %in% probe2gene], collapse = ";"),
    .groups = "drop"
  ) %>%
  arrange(desc(coverage_pct))

cat("=========================================\n")
cat("GEO COVERAGE (GPL570 — U133 Plus 2.0)\n")
cat("=========================================\n")
print(coverage, n = 20)

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
write.csv(coverage, "results/tables/geo_coverage.csv", row.names = FALSE)
cat("\nSaved: results/tables/geo_coverage.csv\n")
# 1. Log Deviation 3
deviation_text <- c(
  "",
  "## Deviation 3: GEO microarray coverage fallback (2026-10-09)",
  "",
  "Trigger: Phase 2 GEO coverage check on GPL570.",
  "Result: Only 2 of 15 signatures (Paper305, Paper62) had 100% coverage.",
  "Average coverage: ~28%. Five signatures had 0%.",
  "Decision: pivot to RNA-seq external cohorts (GSE81089, GSE72094, CPTAC-3).",
  "Table: results/tables/geo_coverage.csv"
)
existing <- readLines("docs/PROTOCOL_DEVIATIONS.md")
writeLines(c(existing, deviation_text), "docs/PROTOCOL_DEVIATIONS.md")
cat("Deviation 3 logged\n")

# 2. Overlap check
overlap_table <- data.frame(
  paper = c("Paper62", "Paper305"),
  dev_cohort = c("TCGA-LUAD", "TCGA-LUAD"),
  author_validation = c("GSE3141 + GSE37745", "GSE30219, GSE31210, GSE37745+GSE50081"),
  GSE31210 = c("No", "YES"),
  GSE50081 = c("No", "YES"),
  GSE30219 = c("No", "YES"),
  action = c("Keep", "DROP from GPL570 cohorts")
)
write.csv(overlap_table, "docs/cohort_overlap_verified.csv", row.names = FALSE)
print(overlap_table)

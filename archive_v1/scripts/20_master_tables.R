# scripts/20_master_tables.R
# Assembles the four master tables for the manuscript.

suppressPackageStartupMessages({
  library(dplyr)
})

dir.create("results/tables", showWarnings = FALSE, recursive = TRUE)

# ============================================================
# TABLE 1 — Signature characteristics
# ============================================================
sig_files <- list.files("signatures", pattern = "\\.tsv$", full.names = TRUE)

sig_info <- do.call(rbind, lapply(sig_files, function(f) {
  s <- read.delim(f, stringsAsFactors = FALSE)
  fname <- tools::file_path_sans_ext(basename(f))
  parts <- strsplit(fname, "_")[[1]]
  author_year <- if (length(parts) >= 2) parts[2] else NA
  data.frame(
    signature_id  = fname,
    short_id      = parts[1],
    author_year   = author_year,
    n_genes       = nrow(s),
    n_protective  = sum(s$coefficient < 0),
    n_risk        = sum(s$coefficient > 0),
    stringsAsFactors = FALSE
  )
}))

# Fill in publication metadata (edit if any are wrong — these are from your screening)
pub_meta <- data.frame(
  short_id  = c("Paper9",  "Paper51", "Paper53", "Paper62", "Paper91",
                "Paper141","Paper151","Paper159","Paper165","Paper192",
                "Paper193","Paper298","Paper305","Paper352","Paper369"),
  first_author = c("Zhang","Song","Yang","Yao","Lin",
                   "Yang","Diao","He","Jin","Geng",
                   "Zheng","Li","Yu","Zhang","Wu"),
  year      = c(2021, 2021, 2022, 2021, 2024,
                2023, 2022, 2021, 2020, 2021,
                2021, 2023, 2024, 2025, 2022),
  stringsAsFactors = FALSE
)

table1 <- merge(sig_info, pub_meta, by = "short_id", all.x = TRUE)
table1 <- table1[order(table1$short_id), ]
table1 <- table1[, c("signature_id","first_author","year","n_genes",
                     "n_protective","n_risk")]

write.csv(table1, "results/tables/Table1_signatures.csv", row.names = FALSE)
cat("Table 1 saved.\n")
print(table1)

# ============================================================
# TABLE 2 — Cohort characteristics
# ============================================================
table2 <- data.frame(
  cohort          = c("TCGA-LUAD", "CPTAC-3", "GSE31210", "GSE72094"),
  platform        = c("RNA-seq (STAR-Counts)", "RNA-seq (STAR-Counts)",
                      "GPL570 (microarray)", "GPL15048 (targeted panel)"),
  n_samples       = c(NA, 150, 226, 442),
  n_events        = c(NA, 38, 35, NA),
  median_fu_mo    = c(NA, 19.7, 57.3, NA),
  stage_I_II_pct  = c(NA, 76.7, 100.0, NA),
  evaluable_sigs  = c(13, 15, 2, 0),
  role            = c("Development/apparent", "Independent RNA-seq validation",
                      "Independent microarray validation",
                      "Coverage analysis only"),
  stringsAsFactors = FALSE
)

write.csv(table2, "results/tables/Table2_cohorts.csv", row.names = FALSE)
cat("\nTable 2 saved.\n")
print(table2)

# ============================================================
# TABLE 3 — Both C-indexes per signature-cohort pair
# ============================================================
both <- read.csv("results/tables/H3_both_Cindexes.csv", stringsAsFactors = FALSE)

table3 <- both[order(both$signature, both$cohort),
               c("signature","cohort","platform","n","events",
                 "C_raw","C_cox","HR_SD","HR_lo","HR_hi",
                 "dir_ok","divergence","interpretation")]

write.csv(table3, "results/tables/Table3_both_Cindexes.csv", row.names = FALSE)
cat("\nTable 3 saved. Rows:", nrow(table3), "\n")
print(head(table3, 10))

# ============================================================
# TABLE 4 — H2 added value beyond stage (CPTAC-3)
# ============================================================
h2 <- read.csv("results/tables/CPTAC3_H2_added_value.csv", stringsAsFactors = FALSE)

table4 <- h2[order(h2$signature),
             c("signature","n","events",
               "C_stage_only","C_full","delta_C",
               "LRT_p","LRT_p_BH","p_sig_BH",
               "HR_sig_adj","HR_sig_adj_lo","HR_sig_adj_hi",
               "H2_result")]

write.csv(table4, "results/tables/Table4_H2_added_value.csv", row.names = FALSE)
cat("\nTable 4 saved. Rows:", nrow(table4), "\n")
print(table4)

# ============================================================
# Done
# ============================================================
cat("\n========== All four master tables saved ==========\n")
cat("  results/tables/Table1_signatures.csv\n")
cat("  results/tables/Table2_cohorts.csv\n")
cat("  results/tables/Table3_both_Cindexes.csv\n")
cat("  results/tables/Table4_H2_added_value.csv\n")
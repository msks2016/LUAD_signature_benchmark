suppressPackageStartupMessages(library(dplyr))

# ── 1. Load the two tables ──────────────────────────────────────────
transport <- read.csv("results/tables/transportability_by_pair.csv",
                      stringsAsFactors = FALSE)
paired_25 <- read.csv("results/tables/H3_both_Cindexes.csv",
                      stringsAsFactors = FALSE)

# ── 2. Build new main-text Table 3 (transportability) ──────────────
# Rename columns for the paper
table3 <- transport[, c("record_no", "cohort", "n_genes",
                        "raw_C", "slope", "slope_se",
                        "HR_per_SD", "refit_C", "refit_C_corrected",
                        "obs_div", "null_div_prob")]

# Rename for the paper
colnames(table3) <- c(
  "signature", "cohort", "n_genes",
  "raw_C_index", "calibration_slope", "slope_SE",
  "HR_per_SD", "refit_C_apparent", "refit_C_corrected",
  "observed_divergence", "null_divergence_prob"
)

# Add a direction-of-effect column based on slope sign
table3$slope_direction <- ifelse(table3$calibration_slope > 0, "Positive",
                                 ifelse(table3$calibration_slope < 0, "Negative", "Null"))

# Add a classification column based on the raw-C / slope consistency
table3$pair_interpretation <- ifelse(
  table3$calibration_slope > 0.3 & table3$raw_C_index > 0.5, "Consistent positive",
  ifelse(table3$calibration_slope < -0.1, "Consistent negative",
         ifelse(abs(table3$calibration_slope) < 0.15, "Near null",
                "Weak positive")))

write.csv(table3, "results/tables/Table3_transportability.csv",
          row.names = FALSE)

cat("=== Table 3 (main text) — transportability by pair ===\n")
print(table3)

# ── 3. Demote the 25-row paired table to Supplementary Table S3 ────
paired_supp <- paired_25
# Rename columns to match supplementary style
colnames(paired_supp)[colnames(paired_supp) == "signature"]    <- "signature"
colnames(paired_supp)[colnames(paired_supp) == "cohort"]       <- "cohort"
colnames(paired_supp)[colnames(paired_supp) == "C_raw"]        <- "raw_C_index"
colnames(paired_supp)[colnames(paired_supp) == "C_cox"]        <- "cox_model_C_index"
colnames(paired_supp)[colnames(paired_supp) == "divergence"]   <- "raw_minus_cox"

write.csv(paired_supp,
          "results/tables/Supplementary_Table_S3_paired_C_indexes.csv",
          row.names = FALSE)

cat("\n=== Supplementary Table S3 — paired C-indexes (25 rows) ===\n")
print(head(paired_supp, 5))
cat("...\n")
cat("Rows:", nrow(paired_supp), "\n")

# ── 4. Print the summary of changes ─────────────────────────────────
cat("\n=== Summary ===\n")
cat("Table 3 (main text, NEW):  ", nrow(table3), "rows — transportability by pair\n")
cat("Supplementary Table S3:    ", nrow(paired_supp), "rows — paired C-indexes (demoted)\n")
cat("\nFiles written:\n")
cat("  results/tables/Table3_transportability.csv\n")
cat("  results/tables/Supplementary_Table_S3_paired_C_indexes.csv\n")
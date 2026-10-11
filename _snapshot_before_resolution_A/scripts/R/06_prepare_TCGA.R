# ============================================================
# 06_prepare_TCGA.R - Clean expression + survival table
# ============================================================
suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(tidyverse)
})
message("Loading se_rna and clin...")
se_rna <- readRDS("data/raw/TCGA-LUAD/se_rna.rds")
clin   <- readRDS("data/raw/TCGA-LUAD/clinical.rds")
tpm <- assay(se_rna, "tpm_unstrand")
message("TPM matrix: ", nrow(tpm), " x ", ncol(tpm))
gene_info <- as.data.frame(rowData(se_rna)) %>%
  select(gene_id, gene_name, gene_type)
expr_df <- as.data.frame(tpm)
expr_df$gene_id <- rownames(expr_df)
expr_df <- left_join(expr_df, gene_info, by = "gene_id")
keep_types <- c("lncRNA", "antisense", "lincRNA",
                "sense_intronic", "sense_overlapping",
                "3prime_overlapping_ncRNA", "processed_transcript",
                "protein_coding")
expr_filtered <- expr_df %>% filter(gene_type %in% keep_types)
message("Filtered genes: ", nrow(expr_filtered))
sample_cols <- setdiff(colnames(expr_filtered),
                       c("gene_id", "gene_name", "gene_type"))
expr_collapsed <- expr_filtered %>%
  select(gene_name, all_of(sample_cols)) %>%
  group_by(gene_name) %>%
  summarise(across(everything(), ~ max(.x, na.rm = TRUE)),
            .groups = "drop") %>%
  column_to_rownames("gene_name")
message("Final matrix: ", nrow(expr_collapsed), " x ", ncol(expr_collapsed))
rename_map <- c(
  "AL034397.1"  = "AL034397.3",
  "AC010980.1"  = "AC010980.2",
  "AL590226.1"  = "HSPC324",
  "AC079766.1"  = "KIAA0125",
  "AC110792.3"  = "RP11-21L23.2"
)
for (old in names(rename_map)) {
  new <- rename_map[old]
  if (old %in% rownames(expr_collapsed)) {
    if (!(new %in% rownames(expr_collapsed))) {
      rownames(expr_collapsed)[rownames(expr_collapsed) == old] <- new
      message("Renamed ", old, " → ", new)
    } else {
      message("Duplicate: ", old, " and ", new, " both present, keeping ", new)
      expr_collapsed <- expr_collapsed[rownames(expr_collapsed) != old, ]
    }
  }
}
saveRDS(expr_collapsed, "data/processed/tcga_luad_expr_clean.rds")
write.csv(expr_collapsed, "data/processed/tcga_luad_expr_clean.csv")
id_col <- if ("submitter_id" %in% colnames(clin)) "submitter_id" else "bcr_patient_barcode"
message("Using patient ID column: ", id_col)
surv <- clin %>%
  mutate(
    patient = .data[[id_col]],
    OS      = ifelse(vital_status == "Dead", 1, 0),
    OS_time = ifelse(vital_status == "Dead",
                     days_to_death,
                     days_to_last_follow_up)
  ) %>%
  filter(!is.na(OS), !is.na(OS_time), OS_time > 0)
optional_cols <- c("age_at_index", "sex_at_birth",
                   "ajcc_pathologic_stage", "ajcc_pathologic_t",
                   "ajcc_pathologic_n", "ajcc_pathologic_m")
present_optional <- intersect(optional_cols, colnames(surv))
message("Optional clinical columns kept: ", paste(present_optional, collapse = ", "))
surv <- surv %>% select(patient, OS, OS_time, all_of(present_optional))
saveRDS(surv, "data/processed/tcga_luad_survival.rds")
write.csv(surv, "data/processed/tcga_luad_survival.csv", row.names = FALSE)
message("Survival table: ", nrow(surv), " patients")
message("  Alive: ", sum(surv$OS == 0))
message("  Dead:  ", sum(surv$OS == 1))
message("\n=== Script 06 COMPLETE ===")

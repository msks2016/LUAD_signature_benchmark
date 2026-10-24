# scripts/22_extra_tables.R
# Builds Table 5 (robustness summary) + Supplementary Table 2 (coverage matrix).

suppressPackageStartupMessages(library(dplyr))

# ============================================================
# PART A — Table 5: Robustness summary (unchanged)
# ============================================================
cat("===== Table 5: Robustness summary =====\n\n")

r_missing <- read.csv("results/tables/robustness_missing_gene_thresholds.csv",
                      stringsAsFactors = FALSE)
r_stage   <- read.csv("results/tables/robustness_stage_restricted.csv",
                      stringsAsFactors = FALSE)
r_loo     <- read.csv("results/tables/robustness_leave_one_cohort_out.csv",
                      stringsAsFactors = FALSE)
r_early   <- read.csv("results/tables/robustness_early_death.csv",
                      stringsAsFactors = FALSE)

check1 <- data.frame(
  check   = "Missing-gene threshold",
  detail  = paste0(
    "10% rule: ", sum(r_missing$eval_10pct), " evaluable; ",
    "20% rule: ", sum(r_missing$eval_20pct), " evaluable; ",
    "30% rule: ", sum(r_missing$eval_30pct), " evaluable"),
  finding = "Pattern stable; only Paper352 flips at 10/20 boundary",
  stringsAsFactors = FALSE
)

stage_wide <- reshape(r_stage[, c("signature","subset","C_raw","C_cox")],
                      idvar = "signature", timevar = "subset",
                      direction = "wide")
stage_delta <- abs(stage_wide$C_raw.all - stage_wide$C_raw.stage_I_II_only)
check2 <- data.frame(
  check   = "Stage restriction",
  detail  = paste0("Max change in raw C-index across all / stage I-II: ",
                   round(max(stage_delta, na.rm = TRUE), 3)),
  finding = "Directional pattern preserved in stage I-II subset",
  stringsAsFactors = FALSE
)

check3 <- data.frame(
  check   = "Leave-one-cohort-out (Yao 2021)",
  detail  = paste0("Pooled HR range across subsets: ",
                   round(min(r_loo$HR_pooled), 2), "-",
                   round(max(r_loo$HR_pooled), 2),
                   "; I2 = 0% in all subsets"),
  finding = "Pooled HR stable; no single cohort drives the effect",
  stringsAsFactors = FALSE
)

early_wide <- reshape(r_early[, c("signature","subset","C_raw","C_cox")],
                      idvar = "signature", timevar = "subset",
                      direction = "wide")
early_delta <- abs(early_wide$C_raw.all - early_wide$C_raw.no_death_within_3mo)
check4 <- data.frame(
  check   = "Early-death exclusion (no death <3 mo)",
  detail  = paste0("Max change in raw C-index: ",
                   round(max(early_delta, na.rm = TRUE), 3)),
  finding = "Pattern unchanged after removing early deaths",
  stringsAsFactors = FALSE
)

table5 <- rbind(check1, check2, check3, check4)
print(table5)
write.csv(table5, "results/tables/Table5_robustness_summary.csv", row.names = FALSE)
cat("\nSaved Table5_robustness_summary.csv\n\n")

# ============================================================
# PART B — Supplementary Table 2: Coverage matrix
# ============================================================
cat("===== Supplementary Table 2: Coverage matrix =====\n\n")

sig_files <- list.files("signatures", pattern = "\\.tsv$", full.names = TRUE)

# Load expression matrices
tcga_expr     <- readRDS("data/processed/tcga_luad_expr_clean_v3.rds")
cptac_expr    <- readRDS("data/processed/CPTAC3_LUAD_expr_matched.rds")
gse31210_expr <- readRDS("data/processed/GSE31210_expr_matched.rds")

# Coverage: percentage of signature genes present in each matrix
compute_cov <- function(expr_matrix, sig_file) {
  s <- read.delim(sig_file, stringsAsFactors = FALSE)
  pct <- 100 * sum(s$gene %in% rownames(expr_matrix)) / nrow(s)
  round(pct, 1)
}

coverage_matrix <- do.call(rbind, lapply(sig_files, function(f) {
  sig_name <- tools::file_path_sans_ext(basename(f))
  s <- read.delim(f, stringsAsFactors = FALSE)
  data.frame(
    signature     = sig_name,
    n_genes       = nrow(s),
    pct_TCGA      = compute_cov(tcga_expr, f),
    pct_CPTAC3    = compute_cov(cptac_expr, f),
    pct_GSE31210  = compute_cov(gse31210_expr, f),
    stringsAsFactors = FALSE
  )
}))

coverage_matrix$eval_TCGA     <- coverage_matrix$pct_TCGA     >= 80
coverage_matrix$eval_CPTAC3   <- coverage_matrix$pct_CPTAC3   >= 80
coverage_matrix$eval_GSE31210 <- coverage_matrix$pct_GSE31210 >= 80

coverage_matrix <- coverage_matrix[order(coverage_matrix$signature), ]

print(coverage_matrix)
write.csv(coverage_matrix,
          "results/tables/Supplementary_Table_2_coverage_matrix.csv",
          row.names = FALSE)
cat("\nSaved Supplementary_Table_2_coverage_matrix.csv\n")
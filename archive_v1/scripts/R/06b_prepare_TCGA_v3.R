# 06b_prepare_TCGA_v3.R - Log-transform TCGA-LUAD v3 matrix
# TCGA-LUAD v3 was TPM-normalized but not log-transformed; the
# independent cohorts are on log scale. Applies log2(x + 1) so all
# cohorts are on the same scale for per-gene z-scoring.

v3_path <- "data/processed/tcga_luad_expr_clean_v3.rds"
v3 <- readRDS(v3_path)
cat("Before:", range(as.matrix(v3), na.rm = TRUE), "\n")
v3_log <- as.data.frame(log2(as.matrix(v3) + 1))
cat("After: ", range(as.matrix(v3_log), na.rm = TRUE), "\n")
saveRDS(v3_log, v3_path)
message("Done.")

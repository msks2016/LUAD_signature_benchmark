# ============================================================
# 08_random_baseline.R
# Random gene set null distribution (n=1000 per signature)
# Tests whether each signature's C-index beats chance
# ============================================================
suppressPackageStartupMessages({
  library(tidyverse)
  library(survival)
})

set.seed(42)  # reproducibility

message("Loading data...")
expr <- readRDS("data/processed/tcga_luad_expr_clean_v3.rds")
surv <- readRDS("data/processed/tcga_luad_survival.rds")
sigs <- read.csv("data/processed/signature_coefficients.csv")

# --- Filter + align (same as Script 07b) ---
expr <- expr[, substr(colnames(expr), 14, 15) == "01", drop = FALSE]
colnames(expr) <- substr(colnames(expr), 1, 12)
expr <- expr[, !duplicated(colnames(expr)), drop = FALSE]

common <- intersect(colnames(expr), surv$patient)
expr <- expr[, common, drop = FALSE]
surv <- surv[match(common, surv$patient), ]

cat("Patients:", ncol(expr), "| Events:", sum(surv$event), "\n")

# --- Log2 + z-score ---
lexpr <- log2(expr + 1)
sd_vec <- apply(lexpr, 1, sd, na.rm = TRUE)
lexpr  <- lexpr[sd_vec > 0, , drop = FALSE]
zexpr  <- t(scale(t(lexpr)))

# --- Function: compute C-index for a gene set ---
compute_cindex <- function(genes, coefs, mat, surv_df) {
  present <- genes %in% rownames(mat)
  if (sum(present) < 1) return(NA)
  s <- coefs[present]
  sc <- as.numeric(crossprod(mat[genes[present], , drop = FALSE], s))
  d <- data.frame(time = surv_df$time_months, event = surv_df$event, z = sc)
  tryCatch({
    concordance(Surv(time, event) ~ z, data = d, reverse = TRUE)$concordance
  }, error = function(e) NA_real_)
}

# --- Test each signature ---
n_perm <- 1000
results <- list()

for (rid in unique(sigs$record_no)) {
  sig <- sigs %>% filter(record_no == !!rid)
  sig <- sig[!duplicated(sig$gene), ]
  
  # Observed C-index
  obs_c <- compute_cindex(sig$gene, sig$coefficient, zexpr, surv)
  
  # Null distribution: 1000 random gene sets of the SAME size
  n_genes <- nrow(sig)
  all_genes <- rownames(zexpr)
  null_c <- numeric(n_perm)
  
  for (i in seq_len(n_perm)) {
    random_genes <- sample(all_genes, n_genes, replace = FALSE)
    random_coefs <- rnorm(n_genes, 0, 1)  # random directions
    null_c[i] <- compute_cindex(random_genes, random_coefs, zexpr, surv)
  }
  
  null_c <- null_c[!is.na(null_c)]
  p_value <- mean(null_c >= obs_c)
  
  cat(sprintf("%-10s | n_genes=%2d | obs_C=%.3f | null_median=%.3f | null_95pct=%.3f | p=%.4f\n",
              rid, n_genes, obs_c,
              median(null_c, na.rm = TRUE),
              quantile(null_c, 0.95, na.rm = TRUE),
              p_value))
  
  results[[rid]] <- data.frame(
    record_no = rid,
    n_genes = n_genes,
    observed_C = obs_c,
    null_median = median(null_c, na.rm = TRUE),
    null_95pct = quantile(null_c, 0.95, na.rm = TRUE),
    null_max = max(null_c, na.rm = TRUE),
    empirical_p = p_value,
    n_valid_perm = length(null_c)
  )
}

res_df <- do.call(rbind, results)
res_df <- res_df[order(-res_df$observed_C), ]

write.csv(res_df, "results/tables/random_baseline_results.csv", row.names = FALSE)

cat("\n=========================================\n")
cat("RANDOM GENE SET BASELINE\n")
cat("=========================================\n")
print(res_df, row.names = FALSE)

cat("\nSaved: results/tables/random_baseline_results.csv\n")

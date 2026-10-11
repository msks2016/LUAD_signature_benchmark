# scripts/11_null_matched_CPTAC3.R
# Expression-matched, coefficient-weighted null test for CPTAC-3 signatures.

suppressPackageStartupMessages(library(survival))
set.seed(2026)

# --- 1. Load data ---------------------------------------------------
expr   <- readRDS("data/processed/CPTAC3_LUAD_expr_matched.rds")
clin   <- readRDS("data/processed/CPTAC3_LUAD_clin_matched.rds")
scores <- read.csv("results/tables/CPTAC3_LUAD_signature_scores.csv",
                   stringsAsFactors = FALSE)
sig_files <- list.files("signatures", pattern = "\\.tsv$", full.names = TRUE)

stopifnot(all(colnames(expr) == clin$submitter_id))

# (z-scoring disabled; raw log2 expression used for scoring)
# zexpr <- t(scale(t(expr)))  # disabled: raw expression used
expr <- expr_z[complete.cases(expr), ]

# Mean expression per gene for the expression-matched null
mu   <- rowMeans(expr, na.rm = TRUE)
pool <- intersect(rownames(expr_z), names(mu))

# Deciles of mean expression
q   <- quantile(mu[pool], probs = 0:10 / 10, include.lowest = TRUE)
bin <- cut(mu[pool], breaks = q, labels = FALSE, include.lowest = TRUE)
names(bin) <- pool
pool_by_bin <- split(pool, bin)

# --- 2. Helper functions --------------------------------------------
cstat <- function(sc) {
  concordance(Surv(clin$os_time, clin$os_event) ~ sc)$concordance
}

boot_c <- function(sc, B = 1000) {
  n <- length(sc)
  replicate(B, {
    i <- sample(n, replace = TRUE)
    if (sum(clin$os_event[i]) < 5) return(NA_real_)
    concordance(Surv(clin$os_time[i], clin$os_event[i]) ~ sc[i])$concordance
  })
}

# --- 3. Per-signature analysis --------------------------------------
B <- 2000

run_sig <- function(sig_path) {
  sig <- read.delim(sig_path, stringsAsFactors = FALSE)
  n_present <- sum(sig$gene %in% rownames(expr_z)); n_total <- nrow(sig); if (n_present / n_total < 0.80) return(NULL)  # require 100% coverage
  sig <- sig[sig$gene %in% rownames(expr_z), ]
  k   <- nrow(sig)
  if (k == 0) return(NULL)
  
  coefs <- sig$coefficient
  
  # Observed score and C-index
  sc  <- as.numeric(crossprod(coefs, expr_z[sig$gene, , drop = FALSE]))
  obs <- cstat(sc)
  
  # Bootstrap 95% CI for observed
  boot <- boot_c(sc, B = 1000)
  lo   <- quantile(boot, 0.025, na.rm = TRUE)
  hi   <- quantile(boot, 0.975, na.rm = TRUE)
  
  # Expression-matched null
  sig_bins <- bin[sig$gene]
  null <- replicate(B, {
    g <- vapply(sig_bins,
                function(b) sample(pool_by_bin[[as.character(b)]], 1),
                character(1))
    cstat(as.numeric(crossprod(coefs, expr_z[g, , drop = FALSE])))
  })
  
  p_two <- (1 + sum(abs(null - 0.5) >= abs(obs - 0.5))) / (B + 1)
  
  data.frame(
    signature   = tools::file_path_sans_ext(basename(sig_path)),
    n_genes     = k,
    C_observed  = round(obs, 3),
    C_lo        = round(lo, 3),
    C_hi        = round(hi, 3),
    null_2.5    = round(quantile(null, 0.025), 3),
    null_97.5   = round(quantile(null, 0.975), 3),
    p_two_sided = p_two,
    stringsAsFactors = FALSE
  )
}

res <- do.call(rbind, lapply(sig_files, run_sig))

# Holm correction + 3-way classification
res$p_holm <- p.adjust(res$p_two_sided, method = "holm")
res$class  <- ifelse(res$p_holm < 0.05 & res$C_observed > 0.5, "Replicated",
                     ifelse(res$p_holm < 0.05 & res$C_observed < 0.5, "Contradicted",
                            "Inconclusive"))

print(res)
write.csv(res, "results/tables/CPTAC3_signature_null_matched.csv",
          row.names = FALSE)
cat("\n=== Classification summary ===\n")
print(table(res$class))

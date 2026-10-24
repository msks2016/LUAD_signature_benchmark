# scripts/15_score_GSE31210.R
# Applies the two evaluable signatures to GSE31210 and runs the
# expression-matched coefficient null test (2000 draws per signature).
# zexpr <- t(scale(t(expr_31210)))  # disabled: raw expression used
expr_31210 <- expr_31210[complete.cases(expr_31210), ]
cat("Expression:", nrow(expr_31210), "genes x", ncol(expr_31210), "patients\n")

expr_31210 <- readRDS("data/processed/GSE31210_expr_matched.rds")
pool <- intersect(rownames(expr_31210), names(mu))
  sig <- sig[sig$gene %in% rownames(expr_31210), ]
  sc  <- as.numeric(crossprod(coefs, expr_31210[sig$gene, , drop = FALSE]))
    cstat(as.numeric(crossprod(coefs, expr_31210[g, , drop = FALSE])))
    n_patients  = ncol(expr_31210),
zexpr <- t(scale(t(expr_31210)))
zexpr <- zexpr[complete.cases(zexpr), ]
cat("Expression:", nrow(zexpr), "genes x", ncol(zexpr), "patients\n")
cat("Events:", sum(clin_31210$os_event), "\n\n")

# mean-expression bins for the null
mu   <- rowMeans(expr_31210, na.rm = TRUE)
pool <- intersect(rownames(zexpr), names(mu))
q    <- quantile(mu[pool], probs = 0:10 / 10, include.lowest = TRUE)
bin  <- cut(mu[pool], breaks = q, labels = FALSE, include.lowest = TRUE)
names(bin) <- pool
pool_by_bin <- split(pool, bin)

# Only the two evaluable signatures in GSE31210
sig_files <- c("signatures/Paper62_Yao2021.tsv",
               "signatures/Paper305_Yu2024.tsv")

cstat <- function(sc) concordance(Surv(clin_31210$os_time, clin_31210$os_event) ~ sc, reverse = TRUE)$concordance

B <- 2000

run_sig <- function(sig_path) {
  sig <- read.delim(sig_path, stringsAsFactors = FALSE)
  sig <- sig[sig$gene %in% rownames(zexpr), ]
  k   <- nrow(sig)
  coefs <- sig$coefficient
  
  sc  <- as.numeric(crossprod(coefs, zexpr[sig$gene, , drop = FALSE]))
  obs <- cstat(sc)
  
  boot <- replicate(1000, {
    i <- sample(length(sc), replace = TRUE)
    if (sum(clin_31210$os_event[i]) < 5) return(NA_real_)
    concordance(Surv(clin_31210$os_time[i], clin_31210$os_event[i]) ~ sc[i], reverse = TRUE)$concordance
  })
  lo <- quantile(boot, 0.025, na.rm = TRUE)
  hi <- quantile(boot, 0.975, na.rm = TRUE)
  
  sig_bins <- bin[sig$gene]
  null <- replicate(B, {
    g <- vapply(sig_bins,
                function(b) sample(pool_by_bin[[as.character(b)]], 1),
                character(1))
    cstat(as.numeric(crossprod(coefs, zexpr[g, , drop = FALSE])))
  })
  
  p_two <- (1 + sum(abs(null - 0.5) >= abs(obs - 0.5))) / (B + 1)
  
  data.frame(
    signature   = tools::file_path_sans_ext(basename(sig_path)),
    n_genes     = k,
    n_patients  = ncol(zexpr),
    n_events    = sum(clin_31210$os_event),
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
res$p_holm <- p.adjust(res$p_two_sided, method = "holm")
res$class  <- ifelse(res$p_holm < 0.05 & res$C_observed > 0.5, "Replicated",
                     ifelse(res$p_holm < 0.05 & res$C_observed < 0.5, "Contradicted",
                            "Inconclusive"))

print(res)
write.csv(res, "results/tables/GSE31210_signature_null_matched.csv",
          row.names = FALSE)
cat("\n=== Classification summary ===\n")
print(table(res$class))

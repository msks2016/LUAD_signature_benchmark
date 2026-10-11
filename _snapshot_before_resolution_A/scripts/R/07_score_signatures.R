# ============================================================
# 07_score_signatures.R
# Score all signatures in TCGA-LUAD, evaluate apparent performance
# ============================================================
suppressPackageStartupMessages({
  library(tidyverse)
  library(survival)
})
message("Loading data...")
expr <- readRDS("data/processed/tcga_luad_expr_clean_v2.rds")
surv <- readRDS("data/processed/tcga_luad_survival.rds")
sigs <- read.csv("data/processed/signature_coefficients.csv")
message("Expression: ", nrow(expr), " x ", ncol(expr))
message("Survival: ", nrow(surv), " patients")
message("Signatures: ", n_distinct(sigs$record_no))
## ---- 1. Filter to primary tumor (code 01) ----
message("\nSample type breakdown:")
print(table(substr(colnames(expr), 14, 15)))
expr <- expr[, substr(colnames(expr), 14, 15) == "01", drop = FALSE]
colnames(expr) <- substr(colnames(expr), 1, 12)
expr <- expr[, !duplicated(colnames(expr)), drop = FALSE]
message("After filter to 01 + dedup: ", ncol(expr), " samples")
## ---- 2. Align with survival ----
common <- intersect(colnames(expr), surv$patient)
expr   <- expr[, common, drop = FALSE]
surv   <- surv[match(common, surv$patient), ]
stopifnot(identical(colnames(expr), surv$patient))
cat("\n=========================================\n")
cat("Patients:", ncol(expr), " | Events:", sum(surv$event), "\n")
cat("=========================================\n")
## ---- 3. log2(TPM + 1), then z-score per gene ----
message("\nApplying log2(TPM + 1)...")
lexpr <- log2(expr + 1)
message("Removing zero-variance genes...")
sd_vec <- apply(lexpr, 1, sd, na.rm = TRUE)
lexpr  <- lexpr[sd_vec > 0, , drop = FALSE]
message("Genes after filtering: ", nrow(lexpr))
message("Z-scoring per gene...")
zexpr <- t(scale(t(lexpr)))
## ---- 4. Signature scoring function ----
score_sig <- function(sig, mat, max_missing = 0.10) {
  sig <- sig[!duplicated(sig$gene), ]
  present <- sig$gene %in% rownames(mat)
  frac_missing <- mean(!present)
  if (frac_missing > max_missing) {
    return(list(score = NULL, n_present = sum(present),
                n_total = nrow(sig),
                missing = paste(sig$gene[!present], collapse = ";")))
  }
  s <- sig[present, ]
  sc <- as.numeric(crossprod(mat[s$gene, , drop = FALSE], s$coefficient))
  list(score = sc, n_present = sum(present),
       n_total = nrow(sig),
       missing = paste(sig$gene[!present], collapse = ";"))
}
## ---- 5. Evaluate every signature ----
message("\nScoring ", n_distinct(sigs$record_no), " signatures...")
res_list <- list()
for (rid in unique(sigs$record_no)) {
  sig <- sigs %>% filter(record_no == !!rid)
  out <- score_sig(sig, zexpr)
base <- data.frame(
    record_no = rid,
    n_genes   = out$n_total,
    n_present = out$n_present,
    missing   = out$missing,
    stringsAsFactors = FALSE
  )
if (is.null(out$score)) {
    res_list[[rid]] <- cbind(base,
                             status = "not evaluable (missing > 10%)",
                             HR_per_SD = NA, HR_low = NA, HR_high = NA,
                             C_index = NA, logrank_p = NA)
    next
  }
sc <- out$score
  d <- data.frame(
    time  = surv$time_months,
    event = surv$event,
    z     = as.numeric(scale(sc)),
    high  = sc > median(sc)
  )
fit  <- coxph(Surv(time, event) ~ z, data = d)
  ci   <- exp(confint(fit))
  cidx <- concordance(Surv(time, event) ~ z, data = d, reverse = TRUE)$concordance
  lr   <- survdiff(Surv(time, event) ~ high, data = d)
  lr_p <- 1 - pchisq(lr$chisq, 1)
status <- ifelse(out$n_present < out$n_total,
                   "evaluable, genes missing",
                   "evaluable")
res_list[[rid]] <- cbind(
    base,
    status    = status,
    HR_per_SD = exp(coef(fit)),
    HR_low    = ci[1],
    HR_high   = ci[2],
    C_index   = cidx,
    logrank_p = lr_p
  )
}
res <- do.call(rbind, res_list)
res <- res[order(-res$C_index, na.last = TRUE), ]
## ---- 6. Save ----
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
write.csv(res, "results/tables/tcga_apparent_performance.csv", row.names = FALSE)
message("\n=========================================")
message("TCGA APPARENT PERFORMANCE")
message("=========================================")
print(res, row.names = FALSE)
message("\nSaved: results/tables/tcga_apparent_performance.csv")

# ============================================================
# 07b_score_signatures_relaxed.R
# ============================================================
suppressPackageStartupMessages({
  library(tidyverse)
  library(survival)
})

message("Loading data...")
expr <- readRDS("data/processed/tcga_luad_expr_clean_v3.rds")
surv <- readRDS("data/processed/tcga_luad_survival.rds")
sigs <- read.csv("data/processed/signature_coefficients.csv")

message("Expression: ", nrow(expr), " x ", ncol(expr))
message("Survival: ", nrow(surv), " patients")

expr <- expr[, substr(colnames(expr), 14, 15) == "01", drop = FALSE]
patient_id  <- substr(colnames(expr), 1, 12)
aliquot_let <- substr(colnames(expr), 14, 14)
ord <- order(patient_id, aliquot_let, colnames(expr))
expr <- expr[, ord, drop = FALSE]
patient_id <- patient_id[ord]
keep <- !duplicated(patient_id)
expr <- expr[, keep, drop = FALSE]
colnames(expr) <- patient_id[keep]
message("After filter: ", ncol(expr), " samples")

common <- intersect(colnames(expr), surv$patient)
expr <- expr[, common, drop = FALSE]
surv <- surv[match(common, surv$patient), ]

cat("Patients:", ncol(expr), "| Events:", sum(surv$OS), "\n")

lexpr <- expr   # v3 is already log2(TPM + 1)
sd_vec <- apply(lexpr, 1, sd, na.rm = TRUE)
lexpr  <- lexpr[sd_vec > 0, , drop = FALSE]
zexpr  <- t(scale(t(lexpr)))

score_sig <- function(sig, mat, max_coef_missing = 0.25) {
  sig <- sig[!duplicated(sig$gene), ]
  present <- sig$gene %in% rownames(mat)
  coef_missing <- sum(abs(sig$coefficient[!present])) / sum(abs(sig$coefficient))
  if (coef_missing > max_coef_missing) {
    return(list(score = NULL, n_present = sum(present),
                n_total = nrow(sig),
                missing = paste(sig$gene[!present], collapse = ";"),
                coef_missing_frac = coef_missing))
  }
  s <- sig[present, ]
  sc <- as.numeric(crossprod(mat[s$gene, , drop = FALSE], s$coefficient))
  list(score = sc, n_present = sum(present),
       n_total = nrow(sig),
       missing = paste(sig$gene[!present], collapse = ";"),
       coef_missing_frac = coef_missing)
}

res_list <- list()
for (rid in unique(sigs$record_no)) {
  sig <- sigs %>% filter(record_no == !!rid)
  out <- score_sig(sig, zexpr)
  
  base <- data.frame(
    record_no = rid, n_genes = out$n_total,
    n_present = out$n_present, missing = out$missing,
    coef_missing_frac = round(out$coef_missing_frac, 3),
    stringsAsFactors = FALSE)
  
  if (is.null(out$score)) {
    res_list[[rid]] <- cbind(base,
                             status = "not evaluable",
                             HR_per_SD = NA, HR_low = NA, HR_high = NA,
                             C_index = NA, logrank_p = NA)
    next
  }
  
  sc <- out$score
  d <- data.frame(time = surv$OS_time, event = surv$OS,
                  z = as.numeric(scale(sc)),
                  high = sc > median(sc))
  fit <- coxph(Surv(time, event) ~ z, data = d)
  ci <- exp(confint(fit))
  cidx <- concordance(Surv(time, event) ~ z, data = d, reverse = TRUE)$concordance
  lr <- survdiff(Surv(time, event) ~ high, data = d)
  lr_p <- 1 - pchisq(lr$chisq, 1)
  
  status <- ifelse(out$n_present < out$n_total,
                   "evaluable, genes missing", "evaluable")
  
  res_list[[rid]] <- cbind(base, status = status,
                           HR_per_SD = exp(coef(fit)),
                           HR_low = ci[1], HR_high = ci[2],
                           C_index = cidx, logrank_p = lr_p)
}

res <- do.call(rbind, res_list)
res <- res[order(-res$C_index, na.last = TRUE), ]

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
write.csv(res, "results/tables/tcga_apparent_performance_relaxed.csv", row.names = FALSE)

cat("\n=========================================\n")
cat("TCGA APPARENT PERFORMANCE (relaxed rule)\n")
cat("=========================================\n")
print(res, row.names = FALSE)

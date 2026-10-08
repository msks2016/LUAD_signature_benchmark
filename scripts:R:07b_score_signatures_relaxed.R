# ============================================================
# 07b_score_signatures_relaxed.R
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

expr <- expr[, substr(colnames(expr), 14, 15) == "01", drop = FALSE]
colnames(expr) <- substr(colnames(expr), 1, 12)
expr <- expr[, !duplicated(colnames(expr)), drop = FALSE]
message("After filter: ", ncol(expr), " samples")

common <- intersect(colnames(expr), surv$patient)
expr <- expr[, common, drop = FALSE]
surv <- surv[match(common, surv$patient), ]

cat("Patients:", ncol(expr), "| Events:", sum(surv$event), "\n")

lexpr <- log2(expr + 1)
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
  d <- data.frame(time = surv$time_months, event = surv$event,
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
se_rna <- readRDS("data/raw/TCGA-LUAD/se_rna.rds")
gi <- as.data.frame(rowData(se_rna))
genes <- c("AL034397.3", "AC010980.2", "HSPC324",
           "KIAA0125", "RP11-21L23.2", "AC018629.1")
for (g in genes) {
  idx <- which(gi$gene_name == g)
  if (length(idx) > 0) {
    cat(sprintf("%-15s → %s | %s\n", g, gi$gene_id[idx[1]], gi$gene_type[idx[1]]))
  } else {
    cat(sprintf("%-15s → not in current GENCODE\n", g))
  }
}
se_rna <- readRDS("data/raw/TCGA-LUAD/se_rna.rds")
gi <- as.data.frame(rowData(se_rna))
ids <- gi$gene_id
names <- gi$gene_name
types <- gi$gene_type
paper_ids <- c(
  "AL034397.3"    = "ENSG00000272630",
  "AC010980.2"    = "ENSG00000272631",
  "HSPC324"       = "ENSG00000228401",
  "KIAA0125"      = "ENSG00000272646",
  "RP11-21L23.2"  = "ENSG00000272650",
  "AC018629.1"    = "ENSG00000279806"
)
for (pname in names(paper_ids)) {
  base_id <- paper_ids[[pname]]
  idx <- grep(paste0("^", base_id), ids)
  if (length(idx) > 0) {
    cat(sprintf("%-15s (%s) → %s | %s\n",
                pname, base_id, names[idx[1]], types[idx[1]]))
  } else {
    cat(sprintf("%-15s (%s) → NOT in current GENCODE\n", pname, base_id))
  }
}
expr2 <- readRDS("data/processed/tcga_luad_expr_clean_v2.rds")
gn <- rownames(expr2)
for (n in c("AL731563.3", "AL731563.1",
            "AC067750.1",
            "AL034397.1", "AL034397.3",
            "AC010980.1", "AC010980.2",
            "AL590226.1", "HSPC324",
            "AC079766.1", "KIAA0125",
            "AC110792.3", "RP11-21L23.2",
            "AC018629.1")) {
  cat(sprintf("  %-15s present: %s\n", n, n %in% gn))
}
se_rna <- readRDS("data/raw/TCGA-LUAD/se_rna.rds")
gi <- as.data.frame(rowData(se_rna))
for (n in c("AL731563.3", "AL034397.1", "AL034397.3",
            "AC067750.1", "AC010980.1", "AC010980.2")) {
  idx <- which(gi$gene_name == n)
  if (length(idx) > 0) {
    cat(sprintf("%-15s → %s | %s\n",
                n, gi$gene_id[idx[1]], gi$gene_type[idx[1]]))
  } else {
    cat(sprintf("%-15s → NOT in se_rna\n", n))
  }
}
tpm <- assay(se_rna, "tpm_unstrand")
gi <- as.data.frame(rowData(se_rna))
for (n in c("AL731563.3", "AL034397.1",
            "AC067750.1", "AC010980.1")) {
  idx <- which(gi$gene_name == n)
  if (length(idx) > 0) {
    x <- tpm[idx[1], ]
    cat(sprintf("%-15s: max=%g, mean=%g, non-zero=%d/%d\n",
                n, max(x, na.rm=TRUE), mean(x, na.rm=TRUE),
                sum(x > 0, na.rm=TRUE), length(x)))
  }
}

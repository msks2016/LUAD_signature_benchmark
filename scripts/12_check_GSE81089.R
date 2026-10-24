# ============================================================
# 12_check_GSE81089.R
# Inspect GSE81089 before committing as validation cohort
# ============================================================
suppressPackageStartupMessages({
  library(GEOquery)
  library(tidyverse)
})

dir.create("data/raw/geo", recursive = TRUE, showWarnings = FALSE)

message("Downloading GSE81089 (may take 3-10 min from China)...")
gse <- getGEO("GSE81089", destdir = "data/raw/geo", GSEMatrix = TRUE)
if (is.list(gse)) gse <- gse[[1]]

pheno <- pData(gse)
cat("=== GSE81089 ===\n")
cat("Samples:", nrow(pheno), "\n")
cat("Platform:", annotation(gse), "\n\n")

# (a) Survival columns
cat("=== Survival-related columns ===\n")
surv_cols <- grep("surv|death|vital|follow|status|time|relapse|os\\.|dfs",
                  colnames(pheno), value = TRUE, ignore.case = TRUE)
print(surv_cols)

if (length(surv_cols) > 0) {
  cat("\n--- Sample values ---\n")
  for (c in head(surv_cols, 3)) {
    cat("\n", c, ":\n", sep = "")
    print(head(pheno[[c]], 10))
  }
}

# (b) Histology
cat("\n=== Histology-related columns ===\n")
hist_cols <- grep("hist|subtype|diagnos|type|characteristics|tumor",
                  colnames(pheno), value = TRUE, ignore.case = TRUE)
for (c in head(hist_cols, 5)) {
  cat("\n---", c, "---\n")
  print(table(pheno[[c]], useNA = "ifany"))
}

# (c) Expression matrix
cat("\n=== Expression matrix ===\n")
expr <- exprs(gse)
cat("Dimensions:", nrow(expr), "x", ncol(expr), "\n")
cat("Range:", range(expr, na.rm = TRUE), "\n")
cat("Row name format (first 5):\n")
print(head(rownames(expr), 5))

# (d) Signature coverage
sigs <- read.csv("data/processed/signature_coefficients.csv")
sig_genes <- unique(sigs$gene)

cat("\n=== Signature gene coverage ===\n")
present_direct <- sig_genes[sig_genes %in% rownames(expr)]
cat("Direct symbol match:", length(present_direct), "/", length(sig_genes), "\n")

if (grepl("^ENSG", rownames(expr)[1])) {
  expr_base <- sub("\\..*$", "", rownames(expr))
  present_ens <- sig_genes[sig_genes %in% expr_base]
  cat("Ensembl-base match:", length(present_ens), "/", length(sig_genes), "\n")
}

saveRDS(gse, "data/raw/geo/GSE81089.rds")
cat("\nSaved: data/raw/geo/GSE81089.rds\n")

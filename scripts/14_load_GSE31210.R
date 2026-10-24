# scripts/14_load_GSE31210.R
# Harmonize GSE31210 expression and clinical data to HGNC symbols,
# matching the CPTAC-3 pipeline format.

suppressPackageStartupMessages({
  library(GEOquery)
  library(Biobase)
})

dir.create("data/raw/GSE31210", showWarnings = FALSE, recursive = TRUE)
dir.create("data/processed", showWarnings = FALSE, recursive = TRUE)

# --- 1. Load from local cache ---------------------------------------
cat("Loading GSE31210...\n")
gset <- getGEO("GSE31210", GSEMatrix = TRUE,
               destdir = "data/raw/GSE31210")
if (length(gset) > 1) gset <- gset[[1]]
gset <- gset[[1]]

expr  <- exprs(gset)
pheno <- pData(gset)
fdata <- fData(gset)

cat("Expression:", nrow(expr), "probes x", ncol(expr), "samples\n")

# --- 2. Map probes -> HGNC symbols ----------------------------------
sym_col <- grep("symbol", colnames(fdata), value = TRUE, ignore.case = TRUE)[1]
cat("Using symbol column:", sym_col, "\n")

# Split multi-gene mappings and take the first symbol only
genes_first <- vapply(strsplit(fdata[[sym_col]], " ?/// ?"),
                      function(x) if (length(x) > 0) x[1] else NA_character_,
                      character(1))
genes_first[genes_first == ""] <- NA_character_

keep <- !is.na(genes_first)
expr <- expr[keep, , drop = FALSE]
genes_first <- genes_first[keep]

# Collapse probes -> genes: keep the probe with highest mean expression
ord <- order(-rowMeans(expr, na.rm = TRUE))
expr <- expr[ord, , drop = FALSE]
genes_first <- genes_first[ord]
keep_first <- !duplicated(genes_first)
expr <- expr[keep_first, , drop = FALSE]
rownames(expr) <- genes_first[keep_first]

cat("After gene collapse:", nrow(expr), "genes\n")

# --- 3. Subset to signature genes -----------------------------------
sig_files <- list.files("signatures", pattern = "\\.tsv$", full.names = TRUE)
sig_genes <- unique(unlist(lapply(sig_files, function(f) read.delim(f)$gene)))
keep_sig  <- intersect(sig_genes, rownames(expr))
cat("Signature genes present:", length(keep_sig), "/", length(sig_genes), "\n")

# Save the FULL gene matrix (not just signature genes)
expr_log <- log2(expr + 1)   # ensure log scale

saveRDS(expr_log, "data/processed/GSE31210_expr.rds")
cat("Saved expr:", nrow(expr_log), "x", ncol(expr_log), "\n")

# --- 4. Clinical -----------------------------------------------------
clin <- data.frame(
  sample_id    = pheno$geo_accession,
  os_days      = suppressWarnings(as.numeric(pheno$`days before death/censor:ch1`)),
  vital_status = pheno$`death:ch1`,
  stage        = pheno$`pathological stage:ch1`,
  sex          = pheno$gender,
  age          = suppressWarnings(as.numeric(pheno$`age:ch1`)),
  stringsAsFactors = FALSE
)

clin$os_time  <- clin$os_days / 30.44
clin$os_event <- as.integer(tolower(clin$vital_status) == "dead")

# Remove samples with missing survival or invalid times
clin <- clin[!is.na(clin$os_time) & clin$os_time > 0, ]
clin <- clin[clin$os_event %in% c(0L, 1L), ]

# Align expr to remaining samples
expr_final <- expr_log[, colnames(expr_log) %in% clin$sample_id, drop = FALSE]
clin_final <- clin[match(colnames(expr_final), clin$sample_id), ]

saveRDS(expr_final, "data/processed/GSE31210_expr_matched.rds")
saveRDS(clin_final, "data/processed/GSE31210_clin_matched.rds")

cat("\nFinal matched cohort:\n")
cat("  Expression:", nrow(expr_final), "genes x", ncol(expr_final), "samples\n")
cat("  Clinical:  ", nrow(clin_final), "patients\n")
cat("  Events:    ", sum(clin_final$os_event), "\n")
cat("  Median OS: ", round(median(clin_final$os_time), 1), "months\n")
cat("  OS range:  ", round(range(clin_final$os_time), 1), "months\n")
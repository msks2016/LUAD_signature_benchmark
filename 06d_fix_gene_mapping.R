suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(tidyverse)
})
message("Loading se_rna + v2 matrix...")
se_rna <- readRDS("data/raw/TCGA-LUAD/se_rna.rds")
expr2 <- readRDS("data/processed/tcga_luad_expr_clean_v2.rds")
tpm <- assay(se_rna, "tpm_unstrand")
gi <- as.data.frame(rowData(se_rna))
cat("v2 has:", nrow(expr2), "genes x", ncol(expr2), "samples\n")
get_gene_row <- function(gene_name, target_cols) {
  idx <- which(gi$gene_name == gene_name)
  if (length(idx) == 0) {
    warning("Gene not in se_rna: ", gene_name)
    return(NULL)
  }
  x <- tpm[idx[1], ]
  names(x) <- substr(names(x), 1, 12)
  x <- x[!duplicated(names(x))]
  lookup <- setNames(as.numeric(x), names(x))
  result <- lookup[substr(target_cols, 1, 12)]
  names(result) <- target_cols
  return(result)
}
wrong_names <- c("AL034397.1", "AL034397.3",
                 "AC010980.1", "AC010980.2")
present_wrong <- intersect(wrong_names, rownames(expr2))
wrong_names <- c("AL034397.1", "AL034397.3",
                 "AC010980.1", "AC010980.2")
present_wrong <- intersect(wrong_names, rownames(expr2))
if (length(present_wrong) > 0) {
  expr2 <- expr2[!rownames(expr2) %in% present_wrong, ]
  cat("Removed wrong entries:",
      paste(present_wrong, collapse = ", "), "\n")
}
mapping_to_add <- c(
  "AL034397.3"  = "AL731563.3",
  "AC010980.2"  = "AC067750.1"
)
for (paper_name in names(mapping_to_add)) {
  current_name <- mapping_to_add[paper_name]
  row <- get_gene_row(current_name, colnames(expr2))
  if (!is.null(row)) {
    expr2 <- rbind(expr2, row)
    rownames(expr2)[nrow(expr2)] <- paper_name
    cat(sprintf("Added %s (values from %s)\n", paper_name, current_name))
  }
}
for (current_name in c("AL731563.3", "AC067750.1")) {
  if (!(current_name %in% rownames(expr2))) {
    row <- get_gene_row(current_name, colnames(expr2))
    if (!is.null(row)) {
      expr2 <- rbind(expr2, row)
      rownames(expr2)[nrow(expr2)] <- current_name
      cat("Added current-name entry:", current_name, "\n")
    }
  }
}
cat("\nFinal dims:", nrow(expr2), "x", ncol(expr2), "\n")
saveRDS(expr2, "data/processed/tcga_luad_expr_clean_v3.rds")
write.csv(expr2, "data/processed/tcga_luad_expr_clean_v3.csv")
cat("\nVerification (all should be TRUE):\n")
for (n in c("AL034397.3", "AL731563.3",
            "AC010980.2", "AC067750.1",
            "HSPC324", "AL590226.1",
            "KIAA0125", "AC079766.1",
            "RP11-21L23.2", "AC110792.3",
            "AC018629.1")) {
  cat(sprintf("  %-15s present: %s\n", n, n %in% rownames(expr2)))
}
message("=== v3 COMPLETE ===")
expr3 <- readRDS("data/processed/tcga_luad_expr_clean_v3.rds")
cat("Dims:", dim(expr3), "\n\n")
for (n in c("AL034397.3", "AC010980.2",
            "HSPC324", "AL590226.1",
            "KIAA0125", "AC079766.1",
            "RP11-21L23.2", "AC110792.3",
            "AC018629.1")) {
  cat(sprintf("  %-15s present: %s\n", n, n %in% rownames(expr3)))
}

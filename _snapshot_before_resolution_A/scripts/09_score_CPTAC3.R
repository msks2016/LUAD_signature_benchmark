# scripts/09_score_CPTAC3.R
# Matches signature genes to matrix by exact gene name (no suffix stripping).

suppressPackageStartupMessages(library(dplyr))

expr <- readRDS("data/processed/CPTAC3_LUAD_expr_matched.rds")
clin <- readRDS("data/processed/CPTAC3_LUAD_clin_matched.rds")
cat("Expression:", nrow(expr), "genes x", ncol(expr), "patients\n")

# expr_z <- t(scale(t(expr)))  # disabled: raw expression used
expr <- expr[complete.cases(expr), ]
cat("After NA filter:", nrow(expr), "genes\n")

sig_files <- list.files("signatures", pattern = "\\.tsv$", full.names = TRUE)
cat("Signatures found:", length(sig_files), "\n")

score_one <- function(sig_path) {
  sig <- read.delim(sig_path, stringsAsFactors = FALSE)
  sig_name <- tools::file_path_sans_ext(basename(sig_path))
  
  present     <- intersect(sig$gene, rownames(expr))
  missing_pct <- 1 - length(present) / nrow(sig)
  cat(sprintf("  %-25s  %d/%d present (%.1f%% missing)\n",
              sig_name, length(present), nrow(sig), 100 * missing_pct))
  
  if (missing_pct > 0.20) {
    return(data.frame(submitter_id = colnames(expr),
                      signature    = sig_name,
                      risk_score   = NA_real_,
                      evaluable    = FALSE))
  }
  
  sig_sub <- sig[sig$gene %in% present, ]
  mat <- expr[sig_sub$gene, , drop = FALSE]
  scores <- as.numeric(crossprod(sig_sub$coefficient, mat))
  
  data.frame(submitter_id = colnames(expr),
             signature    = sig_name,
             risk_score   = scores,
             evaluable    = TRUE)
}

all_scores <- bind_rows(lapply(sig_files, score_one))

clin_join <- clin[, c("submitter_id", "os_time", "os_event",
                      "sex_at_birth", "ajcc_pathologic_stage")]
all_scores <- merge(all_scores, clin_join, by = "submitter_id", all.x = TRUE)
all_scores <- all_scores[all_scores$evaluable, ]

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
write.csv(all_scores,
          "results/tables/CPTAC3_LUAD_signature_scores.csv",
          row.names = FALSE)

cat("\nSaved:", nrow(all_scores), "rows\n")
cat("Unique patients:", length(unique(all_scores$submitter_id)), "\n")
cat("Unique signatures:", length(unique(all_scores$signature)), "\n")
